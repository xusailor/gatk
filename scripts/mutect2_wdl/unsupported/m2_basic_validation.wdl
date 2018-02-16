#
# This is the validation step of a workflow that requires four bams, a discovery tumor-normal pair and a validation
# tumor-normal pair.  Typically, the validation pair would be separate replicates from a different sequencing process
# e.g. matched WGS vs discovery exome bams.

# In the first step of the workflow we run M2 on both pairs, generating a reassembly bamout for the validation pair and a
# vcf callset for the discovery pair.  The bamout is important for indel validation because it yields pileups that are
# consistent with M2's local reassembly.
#
# In this step we validate the discovery callset against the validation bamout pileup
#
# This workflow is not recommended for use with RNA as validation, due to biases in RNA pileups.
#
workflow m2_validation {
    File ref_fasta
    File ref_fai
    File ref_dict

    File calls_vcf
    File calls_vcf_index

    File validation_bamout
    File validation_bamout_bai
    File validation_normal_bam
    File validation_normal_bai

    String discovery_tumor_sample
    String validation_tumor_sample

    Int? preemptible_attempts
    File? gatk_override
    String gatk_docker

    # Only use reads with a minimum base quality for the base at the variant.
    Int? base_quality_cutoff

    # Delete the reads from the normal and HC sample from the bamout.
    call SelectSingleSample as TumorOnlyValidationBamout {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            bam = validation_bamout,
            sample_to_keep = validation_tumor_sample
    }

    call Validate {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            discovery_tumor_sample = discovery_tumor_sample,
            validation_tumor_bam = TumorOnlyValidationBamout.single_sample_bam,
            validation_tumor_bai = TumorOnlyValidationBamout.single_sample_bai,
            validation_normal_bam = validation_normal_bam,
            validation_normal_bai = validation_normal_bai,
            calls_vcf = calls_vcf,
            calls_vcf_idx = calls_vcf_index,
            base_quality_cutoff = base_quality_cutoff
    }

    output {
        File validation = Validate.validation_output
    }
}

# Validation bams should *not* be RNA.
task Validate {
    File? gatk_override
    String gatk_docker

    File ref_fasta
    File ref_fai
    File ref_dict

    String discovery_tumor_sample

    # For validating M2 these should generally be the bamout, not the bam used for input.
    File validation_tumor_bam
    File validation_tumor_bai
    File validation_normal_bam
    File validation_normal_bai

    File calls_vcf
    File calls_vcf_idx

    Int? base_quality_cutoff

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb

    Int final_mem = select_first([mem, 7])

    String output_file_name = "validation-${discovery_tumor_sample}.tsv"

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR GetSampleName -I ${validation_normal_bam} -O validation_normal_name.txt
        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR GetSampleName -I ${validation_tumor_bam} -O validation_tumor_name.txt

        gatk --java-options "-Xmx${final_mem-1}g" $GATK_JAR ValidateBasicSomaticShortMutations \
            -discv ${discovery_tumor_sample} \
            -V ${calls_vcf} \
            -I ${validation_tumor_bam} \
            -I ${validation_normal_bam} \
            -valcase  `cat validation_tumor_name.txt` \
            -valcontrol `cat validation_normal_name.txt` \
            -O ${output_file_name} \
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
        File validation_output = "${output_file_name}"
    }
}

task SelectSingleSample {
    File? gatk_override
    String gatk_docker

    # Also, removes samples not in the list from the header
    String sample_to_keep
    File bam

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int final_mem = select_first([mem, 3])

    String output_basename = "single_sample"

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
