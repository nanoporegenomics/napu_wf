version 1.0

task dv_t {
  input {
      Int threads
      File reference
	  File bamAlignment
      String dvModel = "ONT_R104"
	  Int memSizeGb = 256
	  Int diskSizeGb = 1024
  }  

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa

    ln -s ~{bamAlignment} reads.bam
    samtools index -@ {threads} reads.bam
    
    /opt/deepvariant/bin/run_deepvariant \
        --model_type ~{dvModel} \
        --ref ref.fa \
        --reads reads.bam \
        --output_vcf dv.vcf.gz \
        --num_shards ~{threads}
  >>>

  output {
	File dvVcf = "dv.vcf.gz"
  }

  runtime {
    docker: "google/deepvariant:cl508467184"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}

task margin_t {
  input {
      File reference
      File vcfFile
	  File bamAlignment
      String sampleName
      Int threads
      String marginOtherArgs = ""
	  Int memSizeGb = 256
	  Int diskSizeGb = 1024
  }  

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ln -s ~{bamAlignment} reads.bam
    samtools index -@ ~{threads} reads.bam

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa
    
    mkdir output/
    margin phase reads.bam ref.fa ~{vcfFile} /opt/margin/params/phase/allParams.haplotag.ont-r104q20.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}

    bgzip output/~{sampleName}.phased.vcf

    samtools index -@ ~{threads} output/~{sampleName}.haplotagged.bam
  >>>

  output {
      File phasedVcf = "output/~{sampleName}.phased.vcf.gz"
      File haplotaggedBam = "output/~{sampleName}.haplotagged.bam"
      File haplotaggedBamIdx = "output/~{sampleName}.haplotagged.bam.bai"
  }

  runtime {
    docker: "mkolmogo/card_harmonize_vcf:0.1"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
