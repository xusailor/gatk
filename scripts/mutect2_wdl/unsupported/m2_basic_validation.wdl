import "mutect2.wdl" as m2

#
# This is the validation step of a workflow that requires a set of bam quadruplets.  Each quadruplet comprises a
# discovery tumor-normal pair to validate and a validation tumor-normal pair.  Typically, the validation pair would be
# separate replicates from a different sequencing process e.g. matched WGS vs discovery exome bams.

# In the first step of the workflow we run M2 on both pairs, generating a reassembly bamout for the validation pair and a
# vcf callset for the discovery pair.  In this step we validate this callset against the pileup information in the validation bamout.
#
# We use pileups from the bamout and not from the original bam in order to obtain a standardized and consistent representation
# for indels.
#
#
# This workflow is not recommended for use with RNA as validation, due to biases in RNA pileups.
#
# Attempts are made to reconcile the pileup validation with the M2 haplotype creation process.
#
# The output is a tar file of validation reports (tsv) for the
#
workflow m2_validation {
    #### M2 parameters
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai

    File validation_bamout
    File validation_bamout_bai
    String tumor_sample_name
    File? normal_bam
    File? normal_bai
    String? normal_sample_name


    Int? preemptible_attempts

    File? gatk_override

    String gatk_docker

    ### Validation parameters
    # A name to identify this group of samples.  This can be arbitrary, but should not contain special characters, nor whitespace.
    String group_id

    # Only use reads with a minimum base quality for the base at the variant.
    Int? base_quality_cutoff

        # Delete the reads from the normal and HC sample from the bamout.
        call rewrite_bam_by_sample as m2_rewrite_bam_by_sample {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                bam = m2_validation_bamout.bamout,
                new_sample_name = m2_validation_bamout.tumor_bam_sample_name,
                output_bam_basename = m2_validation_bamout.tumor_bam_sample_name
        }

        call basic_validator as m2_basic_validator {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                call_intervals = m2_tn.filtered_vcf,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                discovery_tumor_sample_name = m2_tn.tumor_bam_sample_name,
                discovery_normal_sample_name = m2_tn.normal_bam_sample_name,
                validation_tumor_bam = m2_rewrite_bam_by_sample.sample_bam,
                validation_tumor_bai = m2_rewrite_bam_by_sample.sample_bai,
                validation_normal_bam = validation_normal_bam_files[i],
                validation_normal_bai = validation_normal_bam_indices[i],
                vcf_calls = m2_tn.filtered_vcf,
                vcf_calls_idx = m2_tn.filtered_vcf_index,
                entity_id = "m2_" + m2_tn.tumor_bam_sample_name,
                base_quality_cutoff = base_quality_cutoff
        }

        call m2.CollectSequencingArtifactMetrics as validation_normal_CollectSequencingArtifactMetrics {
            input:
                gatk_docker = gatk_docker,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                preemptible_attempts = preemptible_attempts,
                tumor_bam = validation_normal_bam_files[i],
                tumor_bai = validation_normal_bam_indices[i]
        }

    call tar_results as tar_results_m2 {
        input:
            to_be_tarred = m2_basic_validator.validation_output,
            group_id = group_id
    }

    output {
        File m2_tar_file = tar_results_m2.tar_file
        Array[File] validation_normal_pre_adapter_metrics = validation_normal_CollectSequencingArtifactMetrics.pre_adapter_metrics
    }
}

# Validation bams should *not* be RNA.
task basic_validator {
    File? gatk_override
    String gatk_docker
    # Same calls as what is in the VCF
    File call_intervals

    File ref_fasta
    File ref_fai
    File ref_dict

    String discovery_tumor_sample_name
    String discovery_normal_sample_name

    # For validating M2 these should generally be the bamout, not the bam used for input.
    File validation_tumor_bam
    File validation_tumor_bai
    File validation_normal_bam
    File validation_normal_bai

    File vcf_calls
    File vcf_calls_idx

    # Unique name for the entity.  Only used for naming output files.
    String entity_id

    Int? base_quality_cutoff

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb

    Int final_mem = select_first([mem, 7])

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk_override}

        echo "Getting sample names...."
        java -Xmx${final_mem-1}g -jar $GATK_JAR GetSampleName -I ${validation_normal_bam} -O validation_normal_name.txt
        java -Xmx${final_mem-1}g -jar $GATK_JAR GetSampleName -I ${validation_tumor_bam} -O validation_tumor_name.txt
        echo ${discovery_tumor_sample_name}
        echo ${discovery_normal_sample_name}
        echo "Sample names (for validation files):"
        echo `cat validation_tumor_name.txt`
        echo `cat validation_normal_name.txt`
        echo "Running BasicValidatorWalker (incl. filtering info).... "
        java -Xmx${final_mem-1}g -jar $GATK_JAR ValidateBasicSomaticShortMutations \
            -discv ${discovery_tumor_sample_name} \
            -V ${vcf_calls} \
            -I ${validation_tumor_bam} \
            -I ${validation_normal_bam} \
            -valcase  `cat validation_tumor_name.txt` \
            -valcontrol `cat validation_normal_name.txt` \
            -O output_validation_${entity_id}.tsv \
            -R ${ref_fasta} \
            -bqcutoff ${default=20 base_quality_cutoff}

    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: "${final_mem}" + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File validation_output = "output_validation_${entity_id}.tsv"
    }
}

task tar_results {
    Array[File] to_be_tarred
    String group_id

    # Needs to be an image with python installed, but defaults to a bash shell.
    String basic_python_docker="broadinstitute/oncotator:1.9.3.0"
    Int preemptible_attempts=2

    command <<<
    set -e
    python <<CODE
	import shutil
	import os
	os.makedirs("${group_id}_results")
	files = "${sep="," to_be_tarred}".split(",")

	for f in files:
		print(f)
		if f.strip() != "":
			shutil.copy(f, "${group_id}_results/")
	CODE
    tar zcvfh "${group_id}_results.tar.gz" "${group_id}_results"
    >>>

    runtime {
        docker: "${basic_python_docker}"
        preemptible: "${preemptible_attempts}"
        memory: "3 GB"
        disks: "local-disk " + 50 + " HDD"
    }

    output {
        File tar_file = "${group_id}_results.tar.gz"
    }
}

task rewrite_bam_by_sample {
    File? gatk_override
    String gatk_docker

    # Also, removes samples not in the list from the header
    String new_sample_name
    File bam
    String output_bam_basename

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int final_mem = select_first([mem, 3])

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk_override}

        java -Xmx${final_mem-1}g -jar $GATK_JAR PrintReads -I ${bam} -O ${output_bam_basename}.tmp.bam -RF SampleReadFilter -sample ${sep=" -sample " new_sample_name}

        samtools view -H ${output_bam_basename}.tmp.bam > tmpheader.txt

        egrep -v "^\@RG" tmpheader.txt > new_header.txt
        egrep "^\@RG" tmpheader.txt | egrep "${new_sample_name}" >> new_header.txt

        java -Xmx${final_mem-1}g -jar $GATK_JAR ReplaceSamHeader --HEADER new_header.txt -I ${output_bam_basename}.tmp.bam -O ${output_bam_basename}.bam
        samtools index ${output_bam_basename}.bam ${output_bam_basename}.bai
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 3]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File sample_bam = "${output_bam_basename}.bam"
        File sample_bai = "${output_bam_basename}.bai"
    }
}
