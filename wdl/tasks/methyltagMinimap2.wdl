version 1.0

workflow run_methylTag_minimap2 {
    meta {
	author: "Melissa Meredith"
        email: "mmmeredi@ucsc.edu"
        description: "Align an ONT Guppy6.1.2 basecalled  unaligned bam with Remora methyl tags (Mm,Ml) and return a sorted bam file along with it's index (.bai) file. Task uses samtools to convert the unalinged bam into a fastq and piped into minimap2 with hardcoded arguments. Based off https://github.com/jmonlong/wdl-workflows/blob/master/minimap2/workflow.wdl " 
    }
    input {
        File UNALIGNED_METHYL_BAM
        File REF_FILE
        String REF_NAME
        String SAMPLE_NAME
        String CONTAINER = "meredith705/ont_methyl:latest"
        String MINIMAP2_ARGS="-y -x map-ont -a --eqx -k 17 -K 10g"
        String OUT_LABEL=""
        Int CORES = 64
        Int DISK = (4 * round(size(REF_FILE, 'G'))) + (10 * round(size(UNALIGNED_METHYL_BAM, 'G'))) + 1000
        Int MEM = 300
    }

    call fastqAlignAndSortBam {
        input:
        unaligned_methyl_bam  = UNALIGNED_METHYL_BAM,
        ref_file            = REF_FILE,
        ref_name            = REF_NAME,
        sample              = SAMPLE_NAME,
        in_container        = CONTAINER,
        in_args             = MINIMAP2_ARGS,
        in_cores            = CORES,
        in_disk             = DISK,
        in_mem              = MEM
    }
    
    output {
        File out_bam     = fastqAlignAndSortBam.out_bam
        File out_bam_idx = fastqAlignAndSortBam.out_bam_idx
    }
}

task fastqAlignAndSortBam {
    input {
        File unaligned_methyl_bam
        File ref_file
        String ref_name
        String sample
        String in_container = "meredith705/ont_methyl:latest"
        String in_args      = "-y -x map-ont -a --eqx -k 17 -K 10g"
        Int in_cores        = 64
        Int in_disk = 4 * round(size(ref_file, 'G')) + 10 * round(size(unaligned_methyl_bam, 'G')) + 100
        Int in_mem          = 300
    }
    command <<<
    set -eux -o pipefail
    
    samtools fastq -TMm,Ml ~{unaligned_methyl_bam} | minimap2 -t ~{in_cores} ~{in_args} ~{ref_file} - | samtools view -@ ~{in_cores} -bh - | samtools sort -@ ~{in_cores} - > ~{sample}.fastq.cpg.~{ref_name}.bam

    samtools index ~{sample}.fastq.cpg.~{ref_name}.bam

    >>>
    
    output {
        File out_bam     = "~{sample}.fastq.cpg.~{ref_name}.bam"
        File out_bam_idx = "~{sample}.fastq.cpg.~{ref_name}.bam.bai"
    }
    runtime {
        docker: in_container
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
    }
}
