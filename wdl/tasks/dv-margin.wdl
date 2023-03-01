version 1.0

task dv_t {
  input {
      Int threads
      File reference
	  File bamAlignment
	  File bamAlignmentIndex
      String dvModel = "ONT_R104"
      String region = ""
	  Int memSizeGb = 256
	  Int diskSizeGb = 1024
  }  

  String regionArg = if region != "" then "--regions ~{region}" else ""
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa
    
    /opt/deepvariant/bin/run_deepvariant \
        --model_type ~{dvModel} \
        --ref ref.fa \
        --reads ~{bamAlignment} \
        --output_vcf dv.vcf.gz ~{regionArg} \
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
	  File bamAlignmentIndex
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

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa
    
    mkdir output/
    margin phase ~{bamAlignment} ref.fa ~{vcfFile} /opt/margin/params/phase/allParams.haplotag.ont-r104q20.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}

    bgzip output/~{sampleName}.phased.vcf

    samtools index -@ ~{threads} output/~{sampleName}.haplotagged.bam
  >>>

  output {
      File phasedVcf = "output/~{sampleName}.phased.vcf.gz"
      File haplotaggedBam = "output/~{sampleName}.haplotagged.bam"
      File haplotaggedBamIdx = "output/~{sampleName}.haplotagged.bam.bai"
  }

  runtime {
    preemptible: 1
    docker: "mkolmogo/card_harmonize_vcf:0.1"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}

task mergeVCFs {
  input {
      Array[File] vcfFiles
      String outname = "merged"
	  Int memSizeGb = 6
	  Int diskSizeGb = 5 * round(size(vcfFiles, 'G')) + 20
  }  

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    mkdir bcftools.tmp
    bcftools concat -n ~{sep=" " vcfFiles} | bcftools sort -T bcftools.tmp -O z -o ~{outname}.vcf.gz -
    bcftools index -t -o ~{outname}.vcf.gz.tbi ~{outname}.vcf.gz
  >>>

  output {
      File vcf = "~{outname}.vcf.gz"
      File vcfIndex = "~{outname}.vcf.gz.tbi"
  }

  runtime {
    preemptible: 1
    docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    cpu: 1
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}
