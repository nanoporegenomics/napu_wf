version 1.0

# cohort sniffles merge 

workflow run_sniffles_merge{
  input {
    File snf_tsv
    File reference
    String cohort 
    Int numSamples
    Boolean phaseVariants = true
    String extraArgs = ""
  }


  call snifflesMerge_t as snifflesMerge{
    input:
        snf_tsv = snf_tsv,
        referenceFa = reference,
        cohort = cohort,
        numSamples = numSamples,
        phaseVariants = phaseVariants,
        extraArgs = extraArgs
  }

  output{
    File snifflesVcf = snifflesMerge.snifflesVcf
    File snifflesVcfIdx = snifflesMerge.snifflesVcfIdx 
  }
}

task snifflesMerge_t {
  input {
    Int threads = 64
    File snf_tsv 
    File referenceFa
    String cohort 
    Int numSamples 
    Boolean phaseVariants 
    String extraArgs 
    Int memSizeGb = 512
    Int diskSizeGb = 256
    File? resourceLogScript
  }

  String phaseArg = if phaseVariants then "--phase " else ""
  String maxInMemArg = if numSamples>20 then "--combine-max-inmemory-result ~{numSamples}" else ""

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace
    
    # work on this command:
    sniffles ~{phaseArg} ~{maxInMemArg} --input ~{snf_tsv} --vcf ~{cohort}_multisample.vcf.gz --reference ~{referenceFa}

    tabix ~{cohort}_multisample.vcf.gz

  >>>

  output {
    File snifflesVcf = "~{cohort}_multisample.vcf.gz"
    File snifflesVcfIdx = "~{cohort}_multisample.vcf.gz.tbi"
    File? toplog = "top.log"
  }

  runtime {
    preemptible: 0
    docker: "meredith705/card_sniffles:2.7.5"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}
