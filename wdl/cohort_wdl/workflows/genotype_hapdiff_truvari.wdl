version 1.0

workflow truvari_fill_genotypes{
    input {
        Array[File] bedFiles 
        File inVCF
        String out_vcf_name

        String dockerImage = "meredith705/truvari:2.0"

    }

    call fillInGenotypes {
        input:
            bedFiles = bedFiles,
            inVcfFile = inVCF,
            outVcfFilename = out_vcf_name,
            dockerImage = dockerImage

    }


    output{
        File vcf_genotyped = fillInGenotypes.genotyped_sv_vcf
        File vcf_genotyped_idx = fillInGenotypes.genotyped_sv_vcf_idx
    }
}


task fillInGenotypes {
    input {
        Array[File] bedFiles 
        File inVcfFile
        String outVcfFilename
        Int memSizeGB = 128
        Int threadCount = 4
        Int diskSizeGB = 2 * round(size(bedFiles, 'G')) + 50
        String dockerImage 

    }

    String inVcfBasename = basename(inVcfFile, ".vcf.gz") + ".vcf"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # decompress the vcf for use in the genotype script
        bgzip -dc ~{inVcfFile} > ~{inVcfBasename}


        # store variables
        IN_VCF=~{inVcfBasename}
        OUT_VCF=~{outVcfFilename}
        IT=0

        # loop through each bed file and fill in genotypes based on coverage
        while IFS= read -r b
            do
            echo $b
            python3 /opt/genotype_truvari.py $b $IN_VCF $OUT_VCF
            IN_VCF=$OUT_VCF
            OUT_VCF=$IN_VCF
            echo 
            echo "Iteration: $IT | IN: $IN_VCF | OUT: $OUT_VCF"

            (( IT++ ))
        done < ~{write_lines(bedFiles)}

        bgzip ~{outVcfFilename}

        tabix "~{outVcfFilename}.gz"


    >>>
    output {
        File genotyped_sv_vcf = "~{outVcfFilename}.gz"
        File genotyped_sv_vcf_idx = "~{outVcfFilename}.gz.tbi"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }

}
