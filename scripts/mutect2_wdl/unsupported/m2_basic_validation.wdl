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

    File discovery_calls_vcf
    File discovery_calls_vcf_index

    #Note we don't use the validation vcf to validate, but it contains the bam sample names, whihc we do need
    # We could also get these from the bams, but a vcf is a *much* smaller file to localize
    File validation_calls_vcf
    File validation_calls_vcf_index

    # Note that this contains reads from both the tumor and the normal
    File validation_bamout
    File validation_bamout_bai

    Int? preemptible_attempts
    File? gatk_override
    String gatk_docker

    # Only use reads with a minimum base quality for the base at the variant.
    Int? base_quality_cutoff

    call GetSampleNames { input: discovery_calls_vcf = discovery_calls_vcf, validation_calls_vcf = validation_calls_vcf}

    # Delete the reads from the normal and HC sample from the bamout.
    call SingleSampleBam as TumorValidationBamout {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            input_bam = validation_bamout,
            sample_to_keep = GetSampleNames.validation_tumor_sample
    }

    call SingleSampleBam as NormalValidationBamout {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            input_bam = validation_bamout,
            sample_to_keep = GetSampleNames.validation_normal_sample
    }

    call Validate {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            validation_tumor_bamout = TumorValidationBamout.single_sample_bam,
            validation_tumor_bamout_bai = TumorValidationBamout.single_sample_bai,
            validation_normal_bamout = NormalValidationBamout.single_sample_bam,
            validation_normal_bamout_bai = NormalValidationBamout.single_sample_bai,
            discovery_tumor_sample = GetSampleNames.discovery_tumor_sample,
            validation_tumor_sample = GetSampleNames.validation_tumor_sample,
            validation_normal_sample = GetSampleNames.validation_normal_sample,
            discovery_calls_vcf = discovery_calls_vcf,
            discovery_calls_vcf_idx = discovery_calls_vcf_index,
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

    # For validating M2 these should generally be the bamout, not the bam used for input.
    File validation_tumor_bamout
    File validation_tumor_bamout_bai
    File validation_normal_bamout
    File validation_normal_bamout_bai

    File discovery_calls_vcf
    File discovery_calls_vcf_idx

    String discovery_tumor_sample
    String validation_tumor_sample
    String validation_normal_sample

    Int? base_quality_cutoff

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb

    Int final_mem = select_first([mem, 7])

    String output_file_name = "validation.tsv"

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${final_mem-1}g" ValidateBasicSomaticShortMutations \
            -discv ${discovery_tumor_sample} \
            -V ${discovery_calls_vcf} \
            -I ${validation_tumor_bamout} \
            -I ${validation_normal_bamout} \
            -valcase  ${validation_tumor_sample} \
            -valcontrol ${validation_normal_sample} \
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

task GetSampleNames {
    File discovery_calls_vcf
    File validation_calls_vcf

    command {
        grep '##tumor_sample=' ${discovery_calls_vcf} | sed 's/=/ /g' | while read junk sample; do echo $sample; done > discovery_tumor_sample.txt
        grep '##normal_sample=' ${discovery_calls_vcf} | sed 's/=/ /g' | while read junk sample; do echo $sample; done > discovery_normal_sample.txt
        grep '##tumor_sample=' ${validation_calls_vcf} | sed 's/=/ /g' | while read junk sample; do echo $sample; done > validation_tumor_sample.txt
        grep '##normal_sample=' ${validation_calls_vcf} | sed 's/=/ /g' | while read junk sample; do echo $sample; done > validation_normal_sample.txt
    }

    output {
        String discovery_tumor_sample = read_string("discovery_tumor_sample.txt")
        String discovery_normal_sample = read_string("discovery_normal_sample.txt")
        String validation_tumor_sample = read_string("validation_tumor_sample.txt")
        String validation_normal_sample = read_string("validation_normal_sample.txt")

    }
}

task SingleSampleBam {
    File? gatk_override
    String gatk_docker

    # Also, removes samples not in the list from the header
    File input_bam
    String sample_to_keep

    # Runtime parameters
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int final_mem = select_first([mem, 3])

    command <<<
        set -e

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${final_mem-1}g" PrintReads -I ${input_bam} -O TMP.bam -RF SampleReadFilter -sample ${sample_to_keep}

        samtools view -H TMP.bam > tmpheader.txt
        egrep -v "^\@RG" tmpheader.txt > new_header.txt
        egrep "^\@RG" tmpheader.txt | egrep "${sample_to_keep}" >> new_header.txt
        gatk --java-options "-Xmx${final_mem-1}g" ReplaceSamHeader --HEADER new_header.txt -I TMP.bam -O single_sample.bam

        samtools index single_sample.bam single_sample.bai
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 3]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(input_bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File single_sample_bam = "single_sample.bam"
        File single_sample_bai = "single_sample.bai"
    }
}
