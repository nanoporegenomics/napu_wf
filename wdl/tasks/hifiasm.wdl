version 1.0

workflow hifiasm {

    input {
        Array[File] readFiles = []
        String sample_name
        Boolean convert2fastq = false
        String hifiasmArgs = ""
        Int diskSizeGB = 1024
        Int preemptible = 2
    }
    
    File readsFile = select_first(readFiles)

    if(convert2fastq){
      # convert to fastq 
      call convertToFastq {
            input:
            readfiles=readFiles,
            preemptible=preemptible
        }
    }

    
    File readsFastq = select_first([convertToFastq.fastq, readsFile])

    call hifiasm_t {
        input:
        reads=readsFastq,
        sample_name=sample_name,
        hifiasmArgs=hifiasmArgs,
        diskSizeGb=diskSizeGB
    }


  output {
        File hap1_fa = hifiasm_t.asm_hap1_fa
        File hap2_fa = hifiasm_t.asm_hap2_fa
        File hap1_noseq_gfa = hifiasm_t.asm_hap1_gfa
        File hap2_noseq_gfa = hifiasm_t.asm_hap2_gfa
        File gfa = hifiasm_t.asm_gfa
        File hifiasm_log = hifiasm_t.hifiasm_log
    }
}



task hifiasm_t {
    input {
        File reads
        String sample_name
        String hifiasmArgs = ""
        Int threads = 96
        String hifiasmONToption = "--ont"
        Int memSizeGb = 360
        Int diskSizeGb = 1125
        String dockerImage = "meredith705/hifiasm@sha256:c86d7f2750a52a75814828462fe4d6dad5d3980c95ddadb6ed62af629a17aad7"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        hifiasm -t~{threads} ~{hifiasmONToption} ~{hifiasmArgs} -o ~{sample_name}.hifiasm.ont ~{reads} 2> ~{sample_name}.hifiasm.ont.log

        awk '/^S/{print ">"$2;print $3}' ~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.gfa > ~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.gfa > ~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.fa

        bgzip -@ ~{threads} ~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.fa
        bgzip -@ ~{threads} ~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.fa

    

    >>>

    output {
        File asm_gfa = "~{sample_name}.hifiasm.ont.bp.p_ctg.noseq.gfa"
        File asm_hap1_fa = "~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.fa.gz"
        File asm_hap2_fa = "~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.fa.gz"
        File asm_hap1_gfa = "~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.noseq.gfa"
        File asm_hap2_gfa = "~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.noseq.gfa"
        File hifiasm_log = "~{sample_name}.hifiasm.ont.log"

    }

    runtime {
      docker: dockerImage
      cpu: threads
      memory: memSizeGb + " GB"
      disks: "local-disk " + diskSizeGb + " SSD"
}

}


task convertToFastq {
  input {
    Array[File] readfiles = []
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
      samtools fastq -@ ~{threads} $READS  >> ~{outname}.fastq
    done;

    bgzip ~{outname}.fastq
  >>>

  output {
    File fastq = "~{outname}.fastq.gz"
  }

  runtime {
      docker: "meredith705/shasta@sha256:f0b2350446e5772232bbd027ad3f27414d20fd5b26d4a57ca73281593b1e21a2"
      preemptible: preemptible
      cpu: threads
      memory: memSizeGb + " GB"
      disks: "local-disk " + diskSizeGb + " SSD"
  }
}
