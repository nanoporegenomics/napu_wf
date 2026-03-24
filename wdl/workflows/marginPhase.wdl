version 1.0

workflow runMarginPhase {
    input {
        File smallVariantsgVCFFile
        File structuralVariantsFile
        File? harmonizedVariantFile
        File? harmonizedVariantFileIdx
        File refFile
        File bamFile
        String sampleName
        Int preemptible_count = 0
        Int threads = 96
        String? machineType
        String dockerImage = "meredith705/card_harmonize_vcf:0.2"
        File? resourceLogScript
    }

    if(!defined(harmonizedVariantFile)){
        call combineVcfs {
            input:
                smallVariantsFile = smallVariantsgVCFFile,
                structuralVariantsFile = structuralVariantsFile,
                sampleName = sampleName,
                preemptible_count = preemptible_count,
                #threads = threads,
                dockerImage = dockerImage
        }
    }

    File combinedVariantVCF = select_first([harmonizedVariantFile, combineVcfs.outVcf])
    File combinedVariantVCFIdx = select_first([harmonizedVariantFileIdx, combineVcfs.outVcfIdx])

    call marginPhase {
        input:
        combinedVcfFile = combinedVariantVCF,
        combinedVcfFileIdx = combinedVariantVCFIdx,
        refFile = refFile,
        bamFile = bamFile,
        sampleName = sampleName,
        dockerImage = dockerImage,
        preemptible_count = preemptible_count,
        threads = threads,
        machineType = machineType,
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
        File? out_exclusionBed = marginPhase.exclusionBed
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
        File combinedVcfFileIdx
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
        String? machineType
        # reducing 4 * to 3 * temp for large samples
        Int memSizeGb = 3 * round(size(bamFile, 'G')) + 200
        Int diskSizeGb = 2 * round(size(bamFile, 'G')) + round(size(refFile, 'G')) + 100
        Int mem_mb = memSizeGb + 1024
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

        #filter the VCF by depth 
        bash /opt/filter_vcf.sh ~{combinedVcfFile} ~{sampleName} ~{filter_window_size} ~{filter_min_cluster_size} ~{filter_threshold_SD}


        # Make the name of the filterd VCF
        filtVcf="~{sampleName}.merged_small_svs.~{filter_threshold_SD}_sd_depthFilt.vcf.gz"
        mergedFilteredBed="~{sampleName}.merged_small_svs.filt~{filter_window_size}bp_~{filter_threshold_SD}_sds.100kbmerged.bed"
        
        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/

        # if any dense clusters are filtered out use the filtered vcf
        if [ -f "$mergedFilteredBed"] && [ -f "$filtVcf"]; then
            margin phase ~{bamFile} ~{refFile} $filtVcf /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_harm_gvcf 
        else
            margin phase ~{bamFile} ~{refFile} ~{combinedVcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_harm_gvcf
        fi

        # gzip vcf and index bam
        bgzip -@ ~{threads} output/~{sampleName}_harm_gvcf.phased.vcf
        tabix output/~{sampleName}_harm_gvcf.phased.vcf.gz
        samtools index -@ ~{threads} output/~{sampleName}_harm_gvcf.haplotagged.bam

        # consider separating gVCF and SV vcf here


    >>>
    output {
        File phasedVcf = "output/~{sampleName}_harm_gvcf.phased.vcf.gz"
        File phasedVcfIdx = "output/~{sampleName}_harm_gvcf.phased.vcf.gz.tbi"
        File phasedVCFPhaseSetBED = "output/~{sampleName}_harm_gvcf.phaseset.bed"
        File haplotaggedBam = "output/~{sampleName}_harm_gvcf.haplotagged.bam"
        File haplotaggedBamIdx = "output/~{sampleName}_harm_gvcf.haplotagged.bam.bai"
        File? exclusionBed = "~{sampleName}.merged_small_svs.filt~{filter_window_size}bp_~{filter_threshold_SD}_sds.100kbmerged.bed"
        File? toplog = "top.log"
    }

    runtime {
        preemptible: preemptible_count
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        machineType : select_first([machineType, "custom-" + threads + "-" + mem_mb]) 
        docker: dockerImage
    }
}
