version 1.0

task sniffles_t {
  input {
    Int threads
	File bamAlignment
	File bamAlignmentIndex
	File? vntrAnnotations
	Int minSvLen = 25
	Int memSizeGb = 128
	Int diskSizeGb = 256
  }

  String trfString = if defined(vntrAnnotations) then "--tandem-repeats ~{vntrAnnotations}" else ""
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    sniffles -i ~{bamAlignment} -v sniffles.vcf -t ~{threads} ~{trfString} --minsvlen ~{minSvLen} 2>&1 | tee sniffles.log
  >>>

  output {
	File snifflesVcf = "sniffles.vcf"
	File snifflesLog = "sniffles.log"
  }

  runtime {
    preemptible: 1
    docker: "mkolmogo/card_sniffles:2.0.3"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
