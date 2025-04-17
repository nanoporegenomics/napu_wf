version 1.0

workflow shasta {

    input {
        Array[File] readFiles = []
        String shastaArgs = ""
        Boolean inMemory = false
        Boolean qScoreCutoff = false
        Int preemptible = 2
        Int diskSizeGB = 1024
    }
    
    File readsFile = select_first(readFiles)

    if(qScoreCutoff){
      # convert to fastq to remove reads below QscoreMin
      call convertToFastq {
            input:
            readfiles=readFiles,
            preemptible=preemptible
        }
    }

    if(!qScoreCutoff){
      if ((basename(readsFile, ".fasta") == basename(readsFile)) && (basename(readsFile, ".fa") == basename(readsFile))){
        call convertToFasta {
            input:
            readfiles=readFiles,
            preemptible=preemptible
        }
    }
    }
    
    File readsFasta = select_first([convertToFastq.fastq, convertToFasta.fasta, readsFile])

    if(inMemory){
        call shasta_inmem_t {
            input:
            reads=readsFasta,
            shastaArgs=shastaArgs,
            diskSizeGb=diskSizeGB
        }
    }
    
    if(!inMemory){
        call shasta_t {
            input:
            reads=readsFasta,
            shastaArgs=shastaArgs,
            diskSizeGb=diskSizeGB
        }
    }
    
    File shastaFasta = select_first([shasta_t.shastaFasta, shasta_inmem_t.shastaFasta])
    File shastaGfa = select_first([shasta_t.shastaGfa, shasta_inmem_t.shastaGfa])
    File shastaLog = select_first([shasta_t.shastaLog, shasta_inmem_t.shastaLog])

  output {
        File fasta = shastaFasta
        File gfa = shastaGfa
        File log = shastaLog
        File readsToFasta = readsFasta
    }
}

task shasta_t {
  input {
    File reads
    String shastaArgs = ""
    Int threads = 96
    String shastaConfig = "/opt/shasta_config/Nanopore-R10-Fast-Nov2022.conf"
    Int memSizeGb = 624
    Int diskSizeGb = 1125
    String dockerImage = "meredith705/shasta@sha256:f0b2350446e5772232bbd027ad3f27414d20fd5b26d4a57ca73281593b1e21a2"
  }

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    SHASTA_INPUT=~{reads}
    if [ "${SHASTA_INPUT: -3}" == ".gz" ]
    then
      if [ "${SHASTA_INPUT: -4}" == "q.gz" ]
      then
        UNGZIPPED=${SHASTA_INPUT:0:-3}.fasta
        zcat $SHASTA_INPUT | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${UNGZIPPED}
      else
        UNGZIPPED=${SHASTA_INPUT:0:-3}
        zcat $SHASTA_INPUT > ${UNGZIPPED}
      fi
      rm $SHASTA_INPUT
      SHASTA_INPUT=${UNGZIPPED}
    fi

    if [ "${SHASTA_INPUT: -3}" == "bam" ]
    then
      UNGZIPPED=${SHASTA_INPUT:0:-4}.fasta
      samtools fasta $SHASTA_INPUT > $UNGZIPPED
      SHASTA_INPUT=${UNGZIPPED}
    fi

    #In-memory mode, in case you are shure all samples will fit
    #shasta --input $SHASTA_INPUT --config ~{shastaConfig} --threads ~{threads} 2>&1 | tee shasta.log
    
    #By default use disk caching
    shasta --input $SHASTA_INPUT --config ~{shastaConfig} --threads ~{threads} --memoryMode filesystem --memoryBacking disk ~{shastaArgs} 2>&1 | tee shasta.log

    tar -czvf shasta.log.tar.gz shasta.log ShastaRun/performance.log ShastaRun/stdout.log ShastaRun/AssemblySummary.html
  >>>

  output {
    File shastaFasta = "ShastaRun/Assembly.fasta"
    File shastaGfa = "ShastaRun/Assembly.gfa"
    File shastaHtml = "ShastaRun/AssemblySummary.html"
    File shastaLog = "shasta.log.tar.gz"
  }

  #This is optimized for GCP/Terra environemnt to get maximum available RAM. May need to adjust for other cloud environemnts or HPC
  runtime {
    docker: dockerImage
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " LOCAL"
    #cpuPlatform: "Intel Cascade Lake"
  }
}

task shasta_inmem_t {
  input {
    File reads
    String shastaArgs = ""
    Int threads = 80
    String shastaConfig = "/opt/shasta_config/Nanopore-R10-Fast-Nov2022.conf"
    Int memSizeGb = 768
    Int diskSizeGb = 1125
    String cpuPlatform = "AMD Rome"
  }

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    SHASTA_INPUT=~{reads}
    if [ "${SHASTA_INPUT: -3}" == ".gz" ]
    then
      if [ "${SHASTA_INPUT: -4}" == "q.gz" ]
      then
        UNGZIPPED=${SHASTA_INPUT:0:-3}.fasta
        zcat $SHASTA_INPUT | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${UNGZIPPED}
      else
        UNGZIPPED=${SHASTA_INPUT:0:-3}
        zcat $SHASTA_INPUT > ${UNGZIPPED}
      fi
      rm $SHASTA_INPUT
      SHASTA_INPUT=${UNGZIPPED}
    fi

    if [ "${SHASTA_INPUT: -3}" == "bam" ]
    then
      UNGZIPPED=${SHASTA_INPUT:0:-4}.fasta
      samtools fasta $SHASTA_INPUT > $UNGZIPPED
      SHASTA_INPUT=${UNGZIPPED}
    fi

    shasta --input $SHASTA_INPUT --config ~{shastaConfig} --threads ~{threads} ~{shastaArgs} 2>&1 | tee shasta.log
    
    tar -czvf shasta.log.tar.gz shasta.log ShastaRun/performance.log ShastaRun/stdout.log ShastaRun/AssemblySummary.html
  >>>

  output {
    File shastaFasta = "ShastaRun/Assembly.fasta"
    File shastaGfa = "ShastaRun/Assembly.gfa"
    File shastaHtml = "ShastaRun/AssemblySummary.html" 
    File shastaLog = "shasta.log.tar.gz"
  }

  #This is optimized for GCP/Terra environemnt to get maximum available RAM. May need to adjust for other cloud environemnts or HPC
  runtime {
    docker: "meredith705/shasta@sha256:f0b2350446e5772232bbd027ad3f27414d20fd5b26d4a57ca73281593b1e21a2"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " LOCAL"
    cpuPlatform: cpuPlatform
  }
}

task convertToFastq {
  input {
    Array[File] readfiles = []
    Int QscoreMin = 10
    Int threads = 10
    Int memSizeGb = 8
    Int diskSizeGb = 5 * round(size(readfiles, 'G')) + 50
    Int preemptible = 2
  }

  String outname = sub(sub(basename(select_first(readfiles)), ".gz$", ""), ".bam", "")
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    for READS in ~{sep=' ' readfiles}
    do
      if [ "${READS: -3}" == ".gz" ]
      then
        if [ "${READS: -4}" == "q.gz" ]
        then
          zcat $READS | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' >> ~{outname}.fasta
        else
          zcat $READS >> ~{outname}.fasta
        fi
      fi

      if [ "${READS: -3}" == "bam" ]
      then
        #samtools fasta -@ ~{threads} $READS >> ~{outname}.fasta
        # after converting to fastq and chopping, maybe pipe back? do we still need quality values?
        samtools fastq -@ ~{threads} $READS | chopper -t ~{threads} -q ~{QscoreMin} | bgzip >> ~{outname}.Q~{QscoreMin}.fastq.gz
      fi
    done;
  >>>

  output {
    File fastq = "~{outname}.Q~{QscoreMin}.fastq.gz"
  }

  runtime {
      docker: "meredith705/shasta@sha256:f0b2350446e5772232bbd027ad3f27414d20fd5b26d4a57ca73281593b1e21a2"
      preemptible: preemptible
      cpu: threads
      memory: memSizeGb + " GB"
      disks: "local-disk " + diskSizeGb + " SSD"
  }
}


task convertToFasta {
  input {
    Array[File] readfiles = []
    Int threads = 4
    Int memSizeGb = 8
    Int diskSizeGb = 5 * round(size(readfiles, 'G')) + 50
    Int preemptible = 2
  }

  String outname = sub(sub(basename(select_first(readfiles)), ".gz$", ""), ".bam", "")
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    for READS in ~{sep=' ' readfiles}
    do
      if [ "${READS: -3}" == ".gz" ]
      then
        if [ "${READS: -4}" == "q.gz" ]
        then
          zcat $READS | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' >> ~{outname}.fasta
        else
          zcat $READS >> ~{outname}.fasta
        fi
      fi

      if [ "${READS: -3}" == "bam" ]
      then
        samtools fasta -@ ~{threads} $READS >> ~{outname}.fasta
      fi
    done;
  >>>

  output {
    File fasta = "~{outname}.fasta"
  }

  runtime {
      docker: "meredith705/shasta@sha256:f0b2350446e5772232bbd027ad3f27414d20fd5b26d4a57ca73281593b1e21a2"
      preemptible: preemptible
      cpu: threads
      memory: memSizeGb + " GB"
      disks: "local-disk " + diskSizeGb + " SSD"
  }
}
