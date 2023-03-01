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
      Int preemptible = 0
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
    preemptible: preemptible
  }
}

task splitReads {
    input {
        File reads
        Int readsPerChunk
        Int diskGb = 5 * round(size(reads, "G")) + 20
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CHUNK_LINES=$(( ~{readsPerChunk} * 4 ))

        MM_INPUT=~{reads}
        if [ "${MM_INPUT: -3}" == "bam" ]
        then
            samtools fastq -TMm,Ml ~{reads} | split -l $CHUNK_LINES --filter='gzip > ${FILE}.fq.gz' - "fq_chunk.part."
        else
            gzip -cd ~{reads} | split -l $CHUNK_LINES --filter='gzip > ${FILE}.fq.gz' - "fq_chunk.part."
        fi
    >>>
    output {
        Array[File] readChunks = glob("fq_chunk.part.*")
    }
    runtime {
        preemptible: 2
        time: 120
        cpu: 4
        memory: "8 GB"
        disks: "local-disk " + diskGb + " SSD"
        docker: "mkolmogo/card_minimap2:2.23"
    }
}

task mergeBAM {
    input {
        Array[File] bams
        String outname = "merged"
        Int threads = 16
        Int diskGb = round(5 * size(bams, 'G')) + 20
        Int memGb = 8
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools merge -f -p -c --threads ~{threads} ~{outname}.bam ~{sep=" " bams}
        samtools index ~{outname}.bam
    >>>
    output {
        File bam = "~{outname}.bam"
        File bamIndex = "~{outname}.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 240
        memory: memGb + " GB"
        cpu: threads
        disks: "local-disk " + diskGb + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}
