version 1.0

task minimap2_t {
  input {
    Int threads
    File reference
    File reads
	String mapMode = "map-ont"
	Boolean useMd = false
	Int memSizeGb = 128
	Int diskSizeGb = 1024
	Int kmerSize = 17
  }

  String mdString = if useMd then "--MD" else ""

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace
    
    minimap2 -ax ~{mapMode} ~{reference} ~{reads} -k ~{kmerSize} -K 5G -t ~{threads} ~{mdString} | samtools sort -@ 4 -m 4G > minimap2.bam
    samtools index -@ ~{threads} minimap2.bam
  >>>

  output {
    File bam = "minimap2.bam"
	File bamIndex = "minimap2.bam.bai"
  }

  runtime {
    docker: "mkolmogo/card_minimap2:2.23"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
