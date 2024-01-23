version 1.0

task hapdup_t {
  input {
    Int threads
    File alignedBam
	File contigs
    String readType = "ont"
    Int memSizeGb = 256
    Int diskSizeGb = 1024
  }

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #Workaround for Python multiprocessing, fixes "AF_UNIX path too long"
    mkdir -p /tmp
    export TMPDIR="/tmp"

    flye-samtools index -@ 10 ~{alignedBam}
    hapdup --assembly ~{contigs} --bam ~{alignedBam} --out-dir hapdup -t ~{threads} --rtype ~{readType} --use-unphased 2>&1 | tee hapdup.log
  >>>

  output {
    File hapdupDual1 = "hapdup/hapdup_dual_1.fasta"
    File hapdupDual2 = "hapdup/hapdup_dual_2.fasta"
    File hapdupPhased1 = "hapdup/hapdup_phased_1.fasta"
    File hapdupPhased2 = "hapdup/hapdup_phased_2.fasta"
    File hapdupPhaseBed1 = "hapdup/phased_blocks_hp1.bed"
    File hapdupPhaseBed2 = "hapdup/phased_blocks_hp2.bed"
    File hapdupLog = "hapdup.log"
  }

  runtime {
    docker: "mkolmogo/hapdup:0.11"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
    runtime_minutes: 1440
  }
}
