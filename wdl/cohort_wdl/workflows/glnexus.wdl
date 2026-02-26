version 1.0

workflow run_glnexus{
    input {
        Array[File] vcfFiles = []
        File reference
        File? regional_bed
        String out_name
        String configuration = "DeepVariantWGS"
        String? extraArgs = ""
        String dockerImage = "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"

    }

    # Scatter the filterVCF task over the input VCF files to remove SV lines
    scatter (input_vcf in vcfFiles) {
        call filterVCF as filter_vcf {
            input:
                input_vcf = input_vcf,
                reference = reference
        }
    }

    call glnexus as glnexux_merge{
        input:
					vcfFiles = filter_vcf.filtered_gvcf,
					out_name = out_name,
					regional_bed = regional_bed,
                    dockerImage = dockerImage
    }

    output{
        File? merged_gvcf = glnexux_merge.merged_gvcf
        Array[File] gvcfFiles = filter_vcf.filtered_gvcf
        Array[File] svVcfFiles = filter_vcf.filtered_sv_vcf
    }
}

task glnexus {
    input {
        Array[File] vcfFiles = []
        File? regional_bed
        String out_name
        String configuration = "DeepVariantWGS"
        String? extraArgs = ""
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 5 * round(size(vcfFiles, 'G')) + 300
        String dockerImage = "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"

    }

    Boolean regionaly =  defined(regional_bed)

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # note output format: https://nanoporetech.github.io/modkit/intro_bedmethyl.html
        # filtering threshold default 10-th percentile of calls https://github.com/nanoporetech/modkit/blob/master/filtering.md
        if [ ~{regionaly} == true ]
        then
          glnexus_cli \
          --config ~{configuration} \
          --bed ~{regional_bed} \
          ~{sep=" " select_all(vcfFiles)} | bcftools view \
         | bgzip -@ ~{threadCount} > ~{out_name}.deepvariant.cohort.vcf.gz

        else
          glnexus_cli \
          --config ~{configuration} \
          ~{sep=" " vcfFiles} | bcftools view \
         | bgzip -@ ~{threadCount} > ~{out_name}.deepvariant.cohort.vcf.gz
        fi


    >>>

        output {
            File? merged_gvcf = "~{out_name}.deepvariant.cohort.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task filterVCF {
    input {
        File input_vcf
        File reference
        Int memSizeGB = 128
        Int threadCount = 4
        Int diskSizeGB = 3 * round(size(input_vcf, 'G')) + 30
        String dockerImage = "quay.io/mlin/glnexus:v1.2.7"
    }
    
    String filtFile = basename(input_vcf)
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        
        # Remove sv lines and write to the gvcf file
        bgzip -dc ~{input_vcf} | awk -F "\t" '(/^#/ || $3 !~ /svim_asm/)' | bcftools norm -f ~{reference} -m -any --threads 8 -Oz ~{filtFile}.filt_gvcf.norm.vcf.gz
        tabix ~{filtFile}.filt_gvcf.norm.vcf.gz
        

        # select sv lines 
        bgzip -dc ~{input_vcf} | awk -F "\t" '(/^#/ || $3 ~ /svim_asm/)' | bgzip > ~{filtFile}.svim_asm.vcf.gz
        tabix ~{filtFile}.svim_asm.vcf.gz

    >>>

    output {
        File filtered_gvcf = "~{filtFile}.filt_gvcf.norm.vcf.gz"
        File filtered_gvcf_idx = "~{filtFile}.filt_gvcf.norm.vcf.gz.tbi"
        File filtered_sv_vcf = "~{filtFile}.svim_asm.vcf.gz"
        File filtered_sv_vcf_idx = "~{filtFile}.svim_asm.vcf.gz.tbi"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
