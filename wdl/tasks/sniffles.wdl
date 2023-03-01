version 1.0

task sniffles_t {
  input {
    Int threads
	File bamAlignment
	File bamAlignmentIndex
	File vntrAnnotations = ""
	Int minSvLen = 25
	Int memSizeGb = 128
	Int diskSizeGb = 256
  }
  
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    TRF_STRING=""
    if [ ! -z ~{vntrAnnotations} ]
    then
       TRF_STRING="--tandem-repeats ~{vntrAnnotations}"
    fi
    echo $TRF_STRING

    sniffles -i ~{bamAlignment} -v sniffles.vcf -t ~{threads} ${TRF_STRING} --minsvlen ~{minSvLen} 2>&1 | tee sniffles.log
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
