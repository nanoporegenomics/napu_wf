version 1.0

workflow runMarginPhase {
    input {
        File smallVariantsgVCFFile
        File structuralVariantsFile
        File? harmonizedVariantFile
        #File gvcfFile
        File refFile
        File bamFile
        String sampleName
        Int preemptible_count = 0
        Int threads = 64
        String dockerImage = "meredith705/card_harmonize_vcf@sha256:5a0ef5a7a4b58a502b9ac2510eac16ef02db9aae8aabc850d12e183c317729f4"
        File? resourceLogScript
    }

    if(!defined(harmonizedVariantFile)){
        call combineVcfs {
            input:
                smallVariantsFile = smallVariantsgVCFFile,
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
        Int memSizeGb = 2 * round(size(smallVariantsFile, 'G')) + round(size(structuralVariantsFile, 'G')) + 200
        Int diskSizeGb = 2 * round(size(smallVariantsFile, 'G')) + 2 * round(size(structuralVariantsFile, 'G')) + 500
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # check if input is bgzipped or not
        SV_FILTERED=~{structuralVariantsFile}_size_filtered.vcf.gz
        SMALL_FILTERED=~{smallVariantsFile}_size_filtered.vcf.gz

        echo ~{sampleName} > samplename.txt
        
        #-f option supports unzgipped input
        zcat -f ~{structuralVariantsFile} | python3 /opt/vcf_filter_size.py greater ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SV_FILTERED
        tabix -p vcf $SV_FILTERED
        zcat -f ~{smallVariantsFile} | python3 /opt/vcf_filter_size.py less ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SMALL_FILTERED
        tabix -p vcf $SMALL_FILTERED

        bcftools concat -a $SMALL_FILTERED $SV_FILTERED -Oz -o ~{sampleName}.merged_small_svs.vcf.gz
        tabix -p vcf ~{sampleName}.merged_small_svs.vcf.gz


    >>>
    output {
        File outVcf = "~{sampleName}.merged_small_svs.vcf.gz"
        File outVcfIdx = "~{sampleName}.merged_small_svs.vcf.gz.tbi"
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
        String sampleName
        String dockerImage
        String marginOtherArgs = ""
        Int filter_window_size = 100000
        Int filter_min_cluster_size = 10
        Int filter_threshold_SD = 3
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

        #filter the VCF by depth or do this in combineVcf task?
        ./opt/filter_vcf.sh ~{combinedVcfFile} ~{sampleName} ~{filter_window_size} ~{filter_min_cluster_size} ~{filter_threshold_SD}
        # Make the name of the filterd VCF
        filtVcf="${sampleName}.merged_small_svs.${filter_threshold_SD}_sd_depthFilt.vcf.gz"
        
        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} $filtVcf /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_hvcf 

        # gzip vcf and index bam
        bgzip output/~{sampleName}_hvcf.phased.vcf
        tabix output/~{sampleName}_hvcf.phased.vcf.gz
        samtools index -@ ~{threads} output/~{sampleName}_hvcf.haplotagged.bam


    >>>
    output {
        File phasedVcf = "output/~{sampleName}_hvcf.phased.vcf.gz"
        File phasedVcfIdx = "output/~{sampleName}_hvcf.phased.vcf.gz.tbi"
        File phasedVCFPhaseSetBED = "output/~{sampleName}_hvcf.phaseset.bed"
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
