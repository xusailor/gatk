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

    # Only use reads with a minimum base quality for the base at the variant.
    Int? base_quality_cutoff

        # Delete the reads from the normal and HC sample from the bamout.
        call SelectSingleSample {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                bam = m2_validation_bamout.bamout,
                sample_to_keep = m2_validation_bamout.tumor_bam_sample_name,
                output_basename = m2_validation_bamout.tumor_bam_sample_name
        }

        call Validate {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                intervals = m2_tn.filtered_vcf,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                discovery_tumor_sample_name = m2_tn.tumor_bam_sample_name,
                validation_tumor_bam = SelectSingleSample.single_sample_bam,
                validation_tumor_bai = SelectSingleSample.single_sample_bai,
                validation_normal_bam = validation_normal_bam_files[i],
                validation_normal_bai = validation_normal_bam_indices[i],
                vcf_calls = m2_tn.filtered_vcf,
                vcf_calls_idx = m2_tn.filtered_vcf_index,
                base_quality_cutoff = base_quality_cutoff
        }

    output {
        File m2_tar_file = tar_results_m2.tar_file
    }
}

# Validation bams should *not* be RNA.
task Validate {
    File? gatk_override
    String gatk_docker
    File intervals

    File ref_fasta
    File ref_fai
    File ref_dict

    String discovery_tumor_sample_name

    # For validating M2 these should generally be the bamout, not the bam used for input.
    File validation_tumor_bam
    File validation_tumor_bai
    File validation_normal_bam
    File validation_normal_bai

    File vcf_calls
    File vcf_calls_idx

    Int? base_quality_cutoff

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb

    Int final_mem = select_first([mem, 7])

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}



        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR GetSampleName -I ${validation_normal_bam} -O validation_normal_name.txt
        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR GetSampleName -I ${validation_tumor_bam} -O validation_tumor_name.txt

        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR ValidateBasicSomaticShortMutations \
            -discv ${discovery_tumor_sample_name} \
            -V ${vcf_calls} \
            -I ${validation_tumor_bam} \
            -I ${validation_normal_bam} \
            -valcase  `cat validation_tumor_name.txt` \
            -valcontrol `cat validation_normal_name.txt` \
            -O validation.tsv \
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
        File validation_output = "validation.tsv"
    }
}

task SelectSingleSample {
    File? gatk_override
    String gatk_docker

    # Also, removes samples not in the list from the header
    String sample_to_keep
    File bam
    String output_basename

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int final_mem = select_first([mem, 3])

    command <<<
        set -e

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${final_mem-1}g" PrintReads -I ${bam} -O ${output_basename}.tmp.bam -RF SampleReadFilter -sample ${sep=" -sample " sample_to_keep}

        samtools view -H ${output_basename}.tmp.bam > tmpheader.txt
        egrep -v "^\@RG" tmpheader.txt > new_header.txt
        egrep "^\@RG" tmpheader.txt | egrep "${sample_to_keep}" >> new_header.txt
        gatk --java-options "-Xmx${final_mem-1}g" ReplaceSamHeader --HEADER new_header.txt -I ${output_basename}.tmp.bam -O ${output_basename}.bam

        samtools index ${output_basename}.bam ${output_basename}.bai
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 3]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File single_sample_bam = "${output_basename}.bam"
        File single_sample_bai = "${output_basename}.bai"
    }
}
