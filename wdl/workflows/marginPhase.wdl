version 1.0

workflow runMarginPhase {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        File? harmonizedVariantFile
        File gvcfFile
        File refFile
        File bamFile
        String sampleName
        Int preemptible_count = 0
        Int threads = 64
        String dockerImage = "mkolmogo/card_harmonize_vcf:0.1"
        File? resourceLogScript
    }

    if(!defined(harmonizedVariantFile)){
        call combineVcfs {
            input:
                smallVariantsFile = gvcfFile,
                structuralVariantsFile = structuralVariantsFile,
                sampleName = sampleName,
                preemptible_count = preemptible_count,
                threads = threads,
                dockerImage = dockerImage
        }
    }

    File combinedVariantVCF = select_first([harmonizedVariantFile, combineVcfs.outVcf])

    call marginPhase {
        input:
        combinedVcfFile = combinedVariantVCF,
        gVcfFile = gvcfFile,
        refFile = refFile,
        bamFile = bamFile,
        sampleName = sampleName,
        dockerImage = dockerImage,
        preemptible_count = preemptible_count,
        threads = threads,
        resourceLogScript = resourceLogScript
    }

    output {
        File out_margin_phase_svs = marginPhase.phasedVcf
        File out_phasedVcfIdx = marginPhase.phasedVcfIdx
        File out_phasedVCFPhaseSetBED = marginPhase.phasedVCFPhaseSetBED
        #File out_margin_phasedgVcf = marginPhase.phasedgVcf
        #File out_margin_phasedgVCFPhaseSetBED = marginPhase.phasedgVCFPhaseSetBED
        File out_margin_phase_bam = marginPhase.haplotaggedBam
        File out_margin_phase_bam_bai = marginPhase.haplotaggedBamIdx
    }
}

task combineVcfs {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        String sampleName
        String dockerImage
        Int preemptible_count
        Int svLengthCutoff = 25
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
        SV_FILTERED=~{structuralVariantsFile}_size_filtered.vcf
        SMALL_FILTERED=~{smallVariantsFile}_size_filtered.vcf

        echo ~{sampleName} > samplename.txt
        
        #-f option supports unzgipped input
        zcat -f ~{structuralVariantsFile} | python3 /opt/vcf_filter_size.py greater ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SV_FILTERED
        tabix -p vcf $SV_FILTERED
        zcat -f ~{smallVariantsFile} | python3 /opt/vcf_filter_size.py less ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SMALL_FILTERED
        tabix -p vcf $SMALL_FILTERED

        bcftools concat -a $SMALL_FILTERED $SV_FILTERED -o ~{sampleName}.merged_small_svs.vcf
    >>>
    output {
        File outVcf = "~{sampleName}.merged_small_svs.vcf"
    }
    runtime {
        preemptible: preemptible_count
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
        File gVcfFile
        String sampleName
        String dockerImage
        String marginOtherArgs = ""
        Int preemptible_count
        Int threads = 64
        Int memSizeGb = 2 * round(size(bamFile, 'G')) + 200
        Int diskSizeGb = 2 * round(size(bamFile, 'G')) + round(size(refFile, 'G')) + 100
        File? resourceLogScript
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## run a recurrent "top" in the background to monitor resource usage
        if [ ~{resourceLogScript} != "" ]
        then
            bash ~{resourceLogScript} 20 top.log &
        fi
        
        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} ~{combinedVcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_hvcf 

        # gzip vcf and index bam
        bgzip output/~{sampleName}_hvcf.phased.vcf
        tabix output/~{sampleName}_hvcf.phased.vcf.gz
        samtools index -@ ~{threads} output/~{sampleName}_hvcf.haplotagged.bam

        # Don't output a bam (-M) for gVCF phasing
        #margin phase output/~{sampleName}_hvcf.haplotagged.bam ~{refFile} ~{gVcfFile} /opt/margin/params/phase/allParams.haplotag.ont-r104q20.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}.g -M

        #bgzip output/~{sampleName}.g.phased.vcf

    >>>
    output {
        File phasedVcf = "output/~{sampleName}_hvcf.phased.vcf.gz"
        File phasedVcfIdx = "output/~{sampleName}_hvcf.phased.vcf.gz.tbi"
        File phasedVCFPhaseSetBED = "output/~{sampleName}_hvcf.phaseset.bed"
        #File phasedgVcf = "output/~{sampleName}.g.phased.vcf.gz"
        #File phasedgVCFPhaseSetBED = "output/~{sampleName}.g.phaseset.bed"
        File haplotaggedBam = "output/~{sampleName}_hvcf.haplotagged.bam"
        File haplotaggedBamIdx = "output/~{sampleName}_hvcf.haplotagged.bam.bai"
        File? toplog = "top.log"
    }

    runtime {
        preemptible: preemptible_count
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
