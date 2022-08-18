version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
	File bamAlignment
	String mapMode = "ont"
	Int memSizeGb = 256
	Int diskSizeGb = 1024
  }

  String pepperMode = if mapMode == "ont" then "--ont_r9_guppy5_sup" else "--hifi"

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    samtools index -@ 10 ~{bamAlignment}
    run_pepper_margin_deepvariant call_variant -b ~{bamAlignment} -f ~{reference} -o `pwd` -t ~{threads} ~{pepperMode} -p PMDV_FINAL 2>&1 | tee pmdv.log
  >>>

  output {
	File pepperVcf = "PMDV_FINAL.vcf.gz"
	File pepperLog = "pmdv.log"
    #File haplotaggedBam = "intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam"
  }

  runtime {
    docker: "kishwars/pepper_deepvariant:r0.8"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
