version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
    File bamAlignment
    String sample
    String mapMode = "ont"
    String extraArgs = ""
    Int memSizeGb = 256
    Int diskSizeGb = 1024
  }
  String pepperMode = if mapMode == "ont" then "--ont_r10_q20" else "--hifi"

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    samtools index -@ 10 ~{bamAlignment}
    run_pepper_margin_deepvariant call_variant -b ~{bamAlignment} -f ~{reference} -o `pwd` -t ~{threads} -s ~{sample} ~{pepperMode} --phased_output -p ~{sample}_PMDV_FINAL ~{extraArgs} 2>&1 | tee pmdv.log
    samtools index -@ 10 ~{sample}_PMDV_FINAL.haplotagged.bam
  >>>

  output {
    File pepperVcf = "~{sample}_PMDV_FINAL.phased.vcf.gz"
    File pepperLog = "pmdv.log"
    File haplotaggedBam = "~{sample}_PMDV_FINAL.haplotagged.bam"
    File haplotaggedBamIdx = "~{sample}_PMDV_FINAL.haplotagged.bam.bai"
  }

  runtime {
    docker: "kishwars/pepper_deepvariant:r0.8"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}
