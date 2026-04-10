version 1.0

workflow runMarginPhase {
    input {
        File VCFFile
        File refFile
        File bamFile
        String sampleName
        Int preemptible_count = 0
        Int threads = 96
        String dockerImage = "meredith705/card_harmonize_vcf:0.2"
        File? monitoring_script
    }

    call marginPhase {
        input:
        VcfFile = VCFFile,
        refFile = refFile,
        bamFile = bamFile,
        sampleName = sampleName,
        dockerImage = dockerImage,
        preemptible_count = preemptible_count,
        threads = threads,
        monitoring_script = monitoring_script
    }

    output {
        File out_margin_phase_vcf = marginPhase.phasedVcf
        File out_phasedVcfIdx = marginPhase.phasedVcfIdx
        File out_phasedVCFPhaseSetBED = marginPhase.phasedVCFPhaseSetBED
        File out_margin_phase_bam = marginPhase.haplotaggedBam
        File out_margin_phase_bam_bai = marginPhase.haplotaggedBamIdx
        File margin_out_monitor = marginPhase.monitoring_log
    }
}


task marginPhase {
    input {
        File VcfFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage
        String marginOtherArgs = ""
        Int preemptible_count
        Int threads = 96
        Int memSizeGb = 2 * round(size(bamFile, 'G')) + 200
        Int diskSizeGb = 2 * round(size(bamFile, 'G')) + round(size(refFile, 'G')) + 100
        File? monitoring_script
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # create this empty log file to be present in the output even wdl fails
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        tabix ~{VcfFile}        
        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/

        margin phase ~{bamFile} ~{refFile} ~{VcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_hp_vcf


        # gzip vcf and index bam
        bgzip -@ ~{threads} output/~{sampleName}_hp_vcf.phased.vcf
        tabix output/~{sampleName}_hp_vcf.phased.vcf.gz
        samtools index -@ ~{threads} output/~{sampleName}_hp_vcf.haplotagged.bam


    >>>
    output {
        File phasedVcf = "output/~{sampleName}_hp_vcf.phased.vcf.gz"
        File phasedVcfIdx = "output/~{sampleName}_hp_vcf.phased.vcf.gz.tbi"
        File phasedVCFPhaseSetBED = "output/~{sampleName}_hp_vcf.phaseset.bed"
        File haplotaggedBam = "output/~{sampleName}_hp_vcf.haplotagged.bam"
        File haplotaggedBamIdx = "output/~{sampleName}_hp_vcf.haplotagged.bam.bai"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        preemptible: preemptible_count
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
