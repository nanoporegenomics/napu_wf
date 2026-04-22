version 1.0

workflow hifiasm_downsampled {

    input {
        Array[File] readFiles = []
        String sample_name
        Boolean convert2fastq = true
        String hifiasmArgs = ""
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
        hifiasmArgs=hifiasmArgs
    }


  output {
        File hap1_fa = hifiasm_t.asm_hap1_fa
        File hap2_fa = hifiasm_t.asm_hap2_fa
        File hap1_noseq_gfa = hifiasm_t.asm_hap1_gfa
        File hap2_noseq_gfa = hifiasm_t.asm_hap2_gfa
        File hap1_lowq_bed = hifiasm_t.asm_hap1_lowQ_bed
        File hap2_lowq_bed = hifiasm_t.asm_hap2_lowQ_bed
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
        String dockerImage = "meredith705/hifiasm@sha256:00be1f1b4ea950d3199d8c947ab7a62aa3e2d3f8d50e8186c07aa6f6056109fe"
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
        File asm_hap1_lowQ_bed = "~{sample_name}.hifiasm.ont.bp.hap1.p_ctg.lowQ.bed"
        File asm_hap2_lowQ_bed = "~{sample_name}.hifiasm.ont.bp.hap2.p_ctg.lowQ.bed"
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
    Float genomeSize = 3.1
    Int coverage = 50
    Int threads = 10
    Int memSizeGb = 50
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
      samtools fastq -@ ~{threads} $READS >> ~{outname}.fastq
    done;

    # add a counter to fastq reads to avoid name clashes 
    # could add counter ^ as they are combined, but might not avoid the clash if 
    # reads happen to be on the same line in each file..
    bash addCounterToFastqReadNames.sh ~{outname}.fastq | bgzip > ~{outname}.uniqName.fastq.gz

    # use rsusa to downsample to 150 Gbases, by using 50x coverage in 3.1gig refsize
    rasusa reads --genome-size ~{genomeSize}g -c ~{coverage} ~{outname}.uniqName.fastq.gz | bgzip > ~{outname}.uniqName.50x_rasusa_downsampled.fastq.gz


  >>>

  output {
    File fastq = "~{outname}.uniqName.50x_rasusa_downsampled.fastq.gz"
  }

  runtime {
      docker: 
      "meredith705/hifiasm@sha256:00be1f1b4ea950d3199d8c947ab7a62aa3e2d3f8d50e8186c07aa6f6056109fe"
      preemptible: preemptible
      cpu: threads
      memory: memSizeGb + " GB"
      disks: "local-disk " + diskSizeGb + " SSD"
  }
}
