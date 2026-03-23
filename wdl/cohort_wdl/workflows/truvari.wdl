version 1.0

workflow run_truvari_collapse{
    input {
        Array[File] vcfFiles 
        File reference
        File reference_index
        File? regional_bed
        String out_name
        String out_prefix

        # Matching parameters
        Float refdist = 1000
        Float pctsize = 0.75
        Float pctseq = 0.75
        Float pctovl = 0.0
        Int typeignore = 0 
        String passonly = "--passonly"

        # Collapse parameters
        String keep = "first"
        String? extraArgs = ""
        String dockerImage = "meredith705/truvari"

    }

    # Scatter the filterVCF task over the input VCF files to select SVs >=50 bps
    scatter (input_vcf in vcfFiles) {
        call filterVCF as filter_vcf {
            input:
                input_vcf = input_vcf,
                reference = reference
        }
    }

    call mergeVCFs {
        input:
            vcfFiles = filter_vcf.fiftybp_sv_vcf,
            vcfFilesIdxs = filter_vcf.ffiftybp_sv_vcf_idx,
            out_prefix = out_prefix,
            dockerImage = dockerImage

    }

    call truvari as truvari_merge{
        input:
                    vcfFile = mergeVCFs.combined_sv_vcf,
                    vcfFileIdxs = mergeVCFs.combined_sv_vcf_idx,
                    out_name = out_name,
                    refdist = refdist, 
                    pctsize = pctsize, 
                    pctseq = pctseq,
                    pctovl = pctovl,
                    typeignore = typeignore,
                    passonly = passonly,
                    keep = keep,
                    extraArgs = extraArgs,
                    reference = reference,
                    reference_index = reference_index,
                    dockerImage = dockerImage
    }

    output{
        File? merged_truvari_vcf = truvari_merge.merged_vcf
        File? collapsed_vcf = truvari_merge.collapsed_vcf
        File combined_sv_vcf = mergeVCFs.combined_sv_vcf
        #File combined_sv_vcf_idx = mergeVCFs.combined_sv_vcf_idx
        Array[File] vcf50Files = filter_vcf.fiftybp_sv_vcf
        Array[File] vcf50FilesIdxs = filter_vcf.ffiftybp_sv_vcf_idx
    }
}

task truvari {
    input {
        File vcfFile
        File vcfFileIdxs 

        # Matching parameters
        Float refdist 
        Float pctsize 
        Float pctseq 
        Float pctovl 
        Int typeignore  
        String passonly 

        # Collapse parameters
        String keep 
        String? extraArgs  

        File reference
        File reference_index
        String out_name
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 3 * round(size(vcfFile, 'G')) + 300
        String dockerImage 

    }


    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # add the options from input:

        truvari collapse \
            --pctsize ~{pctsize} \
            --pctseq ~{pctseq} \
            ~{passonly} \
            --keep ~{keep} \
            -i ~{vcfFile} \
            -o ~{out_name}.truvari_merged.cohort.vcf \
            -c ~{out_name}.truvari_collapsed.vcf \
            -f ~{reference} ~{extraArgs}

        bgzip ~{out_name}.truvari_merged.cohort.vcf
        bgzip ~{out_name}.truvari_collapsed.vcf

        tabix ~{out_name}.truvari_merged.cohort.vcf.gz
        #tabix ~{out_name}.truvari_collapsed.vcf.gz

    >>>

        output {
            File? merged_vcf = "~{out_name}.truvari_merged.cohort.vcf.gz"
            File? merged_vcf_idx = "~{out_name}.truvari_merged.cohort.vcf.gz.tbi"
            File? collapsed_vcf = "~{out_name}.truvari_collapsed.vcf.gz"
            #File? collapsed_vcf_idx = "~{out_name}.truvari_collapsed.vcf.gz.tbi"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task mergeVCFs {
    input {
        Array[File] vcfFiles = []
        Array[File] vcfFilesIdxs = []
        String out_prefix
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 2 * round(size(vcfFiles, 'G')) + 50
        String dockerImage 

    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # make a space separated list of file names for bcftools merge
        bcftools merge \
            --force-samples \
            --threads ~{threadCount} \
            -m none \
            ~{sep=" " vcfFiles} | bgzip -@ ~{threadCount} > ~{out_prefix}_harmonized_merge.vcf.gz

        bcftools index --tbi ~{out_prefix}_harmonized_merge.vcf.gz

    >>>
    output {
        File combined_sv_vcf = "~{out_prefix}_harmonized_merge.vcf.gz"
        File combined_sv_vcf_idx = "~{out_prefix}_harmonized_merge.vcf.gz.tbi"
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
        String dockerImage = "meredith705/truvari"
    }
    
    String filtFile = basename(input_vcf)
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # select for >=50 bp SVs
        bcftools view -i 'INFO/SVLEN >= 50 | INFO/SVLEN <= -50' ~{input_vcf}| bgzip  > ~{input_vcf}.50bps.vcf.gz

        # index
        tabix ~{input_vcf}.50bps.vcf.gz

    >>>

    output {
        File fiftybp_sv_vcf = "~{input_vcf}.50bps.vcf.gz"
        File ffiftybp_sv_vcf_idx = "~{input_vcf}.50bps.vcf.gz.tbi"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
