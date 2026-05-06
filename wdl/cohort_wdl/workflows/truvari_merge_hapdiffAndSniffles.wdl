version 1.0

workflow run_truvari_collapse{
    input {
        File snifflesSvVcf
        File snifflesSvVcfidx
        File assemblySvVcf
        File assemblySvVcfidx 
        File reference
        File reference_index
        File? regional_bed
        String out_name_truvari
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

    call concatVCFs {
        input:
            snifflesSvVcf = snifflesSvVcf,
            snifflesSvVcfidx = snifflesSvVcfidx,
            assemblySvVcf = assemblySvVcf,
            assemblySvVcfidx = assemblySvVcfidx,
            out_prefix = out_prefix,
            dockerImage = dockerImage

    }

    call truvari as truvari_merge{
        input:
                    vcfFile = concatVCFs.combined_sv_vcf,
                    vcfFileIdxs = concatVCFs.combined_sv_vcf_idx,
                    out_name = out_name_truvari,
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
        File combined_sv_vcf = concatVCFs.combined_sv_vcf
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

task concatVCFs {
    input {
        File snifflesSvVcf
        File snifflesSvVcfidx
        File assemblySvVcf
        File assemblySvVcfidx
        String out_prefix
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 150
        String dockerImage 

    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # concat bcftools merge
        bcftools concat --allow-overlaps -o ~{out_prefix}_snf_hapdiff_concat.vcf.gz ~{assemblySvVcf} ~{snifflesSvVcf}

        bcftools index --tbi ~{out_prefix}_snf_hapdiff_concat.vcf.gz

    >>>
    output {
        File combined_sv_vcf = "~{out_prefix}_snf_hapdiff_concat.vcf.gz"
        File combined_sv_vcf_idx = "~{out_prefix}_snf_hapdiff_concat.vcf.gz.tbi"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }

}

