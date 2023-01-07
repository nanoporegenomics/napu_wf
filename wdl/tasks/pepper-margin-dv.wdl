version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
	File bamAlignment
      String pepperMode = "--ont_r10_q20"  ## --ont_r9_guppy5_sup, --hifi, or --ont_r10_q20
	Int memSizeGb = 256
	  Int diskSizeGb = 1024
  }

  

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    samtools index -@ 10 ~{bamAlignment}
    run_pepper_margin_deepvariant call_variant -b ~{bamAlignment} -f ~{reference} -o `pwd` -t ~{threads} ~{pepperMode} --phased_output -p PMDV_FINAL 2>&1 | tee pmdv.log
    samtools index -@ 10 PMDV_FINAL.haplotagged.bam
  >>>

  output {
	File pepperVcf = "PMDV_FINAL.phased.vcf.gz"
	File pepperLog = "pmdv.log"
    	File haplotaggedBam = "PMDV_FINAL.haplotagged.bam"
    	File haplotaggedBamIdx = "PMDV_FINAL.haplotagged.bam.bai"
  }

  runtime {
    docker: "kishwars/pepper_deepvariant:r0.8"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
