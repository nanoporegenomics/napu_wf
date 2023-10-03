version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
	File? bamAlignment
    File bamAlignmentIndex
    String sampleName
    String extraArguments = ""
	String mapMode = "ont"
	Boolean oneChr = false
	Int memSizeGb = 256
	Int diskSizeGb = 1024
	Int preemptible = 0
    File? resourceLogScript
  }

  String pepperMode = if mapMode == "ont" then "--ont_r9_guppy5_sup" else "--hifi"

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ## run a recurrent "top" in the background to monitor resource usage
    if [[ "~{resourceLogScript}" != "" ]]
    then
        bash ~{resourceLogScript} 20 top.log &
    fi

    # soft link the reference and create an index
    ln -s ~{reference} ref.fa
    samtools faidx ref.fa

    ## soft link the alignment and the index so they are reachable in this directory on Terra
    ln -s ~{bamAlignment} reads.bam
    ln -s ~{bamAlignmentIndex} reads.bam.bai

    ## if BAM has reads only for one chromosome
    ## figure out which one and add argument
    REGION_ARG=""
    if [ ~{oneChr} == true ]
    then
        CONTIG_ID=`head -1 < <(samtools view ~{bamAlignment}) | cut -f3`
        REGION_ARG="--region $CONTIG_ID"
    fi

    run_pepper_margin_deepvariant call_variant -b reads.bam -f ref.fa -o `pwd` -t ~{threads} ~{pepperMode} --phased_output -s ~{sampleName} -p ~{sampleName}_PMDV $REGION_ARG ~{extraArguments} 2>&1 | tee pmdv.log
    if [[ -f "~{sampleName}_PMDV.haplotagged.bam" ]]
    then
      samtools index -@ 10 ~{sampleName}_PMDV.haplotagged.bam
    fi

  >>>

  output {
	File pepperVcf = "~{sampleName}_PMDV.phased.vcf.gz"
	File pepperLog = "pmdv.log"
    File? haplotaggedBam = "~{sampleName}_PMDV.haplotagged.bam"
    File? haplotaggedBamIdx = "~{sampleName}_PMDV.haplotagged.bam.bai"
    File? pepperGVcf = "~{sampleName}_PMDV.g.vcf.gz"
  }

  runtime {
    preemptible: preemptible
    docker: "kishwars/pepper_deepvariant:r0.8"
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
    preemptible: 2
    docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    cpu: 1
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}