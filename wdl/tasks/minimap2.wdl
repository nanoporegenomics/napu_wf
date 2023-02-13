version 1.0

task minimap2_t {
  input {
    Int threads
    File reference
    File reads
	String mapMode = "map-ont"
	Boolean useMd = false
	Boolean useEqx = true
	Int memSizeGb = 128
	Int diskSizeGb = 1024
	  Int kmerSize = 17
      String minibatchSize = "5G"
      Int sortMemgb = "4"
  }

  String mdString = if useMd then "--MD" else ""
  String eqxString = if useEqx then "--eqx" else ""

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace
    
    MM_INPUT=~{reads}
    if [ "${MM_INPUT: -3}" == "bam" ]
    then
      samtools fastq -TMm,Ml ~{reads} | \
        minimap2 -ax ~{mapMode} ~{reference} - -k ~{kmerSize} -y -K ~{minibatchSize} -t ~{threads} ~{mdString} ~{eqxString} | samtools sort -@4 -m ~{sortMemgb}G > minimap2.bam
    else
      minimap2 -ax ~{mapMode} ~{reference} ~{reads} -k ~{kmerSize} -K ~{minibatchSize} -t ~{threads} ~{mdString} ~{eqxString} | samtools sort -@4 -m ~{sortMemgb}G > minimap2.bam
    fi

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
