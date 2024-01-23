version 1.0

task dipcall_t {
  input {
    File ctgsPat
    File ctgsMat
    File reference
    String prefix = "dipcall"
    Int threads = 16
    Int memSizeGb = 64
    Int diskSizeGb = 128
  }

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    /opt/dipcall.kit/samtools faidx ~{reference}
    /opt/dipcall.kit/run-dipcall -t ~{threads / 2} ~{prefix} ~{reference} ~{ctgsPat} ~{ctgsMat} > dipcall.mak
    make -j2 -f dipcall.mak
  >>>

  output {
    File dipcallVcf = "~{prefix}.dip.vcf.gz"
  }

  runtime {
    docker: "mkolmogo/card_dipcall:0.3"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
    runtime_minutes: 120
  }
}
