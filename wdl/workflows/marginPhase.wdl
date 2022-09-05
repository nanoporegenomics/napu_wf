version 1.0

workflow runMarginPhase {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage = "miramastoras/marginphase_sv:latest"
    }

    call combineVcfs {
        input:
            smallVariantsFile = smallVariantsFile,
            structuralVariantsFile = structuralVariantsFile,
            sampleName = sampleName,
            dockerImage = dockerImage
    }

    call marginPhase {
        input:
            combinedVcfFile = combineVcfs.outVcf,
            refFile = refFile,
            bamFile = bamFile,
            sampleName = sampleName,
            dockerImage = dockerImage
    }

    output {
        File out_margin_phase_svs = marginPhase.phasedVcf
    }
}

task combineVcfs {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        String sampleName
        String dockerImage
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # check if input is bgzipped or not
        SV_FILENAME=$(basename -- "~{structuralVariantsFile}")
        SMALLV_FILENAME=$(basename -- "~{smallVariantsFile}")

        if [[ ! $SV_FILENAME =~ \.gz$ ]]; then
            cp ~{structuralVariantsFile} .
            bgzip $SV_FILENAME
            SV_FILENAME="${SV_FILENAME}.gz"
        else
            ln -s ~{structuralVariantsFile}
        fi

        if [[ ! $SMALLV_FILENAME =~ \.gz$ ]]; then
            cp ~{smallVariantsFile} .
            bgzip $SMALLV_FILENAME
            SMALLV_FILENAME="${SMALLV_FILENAME}.gz"
        else
            ln -s ~{smallVariantsFile}
        fi

        tabix -p vcf $SV_FILENAME
        tabix -p vcf $SMALLV_FILENAME
        bcftools concat -a $SMALLV_FILENAME $SV_FILENAME -o ~{sampleName}.merged_small_svs.vcf
    >>>
    output {
        File outVcf = "~{sampleName}.merged_small_svs.vcf"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task marginPhase {
    input {
        File combinedVcfFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} ~{combinedVcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} -o output/~{sampleName} -M
    >>>
    output {
    File phasedVcf = "output/~{sampleName}.phased.vcf"
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
