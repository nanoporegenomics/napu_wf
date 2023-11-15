version 1.0

task modkit {
    input {
        File haplotaggedBam
        File haplotaggedBamBai
        File ref
        File? regional_bed
        String sample_name
        String ref_name = "ref"
        String out_type_filter = "--cpg"
        String partitionTag = "--partition-tag HP"
        String? extraArgs = ""
        Int memSizeGB = 64
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(haplotaggedBam, 'G')) + round(size(ref, 'G')) + 100
        String dockerImage = "meredith705/modkit:latest"
        File? resourceLogScript
    }

    parameter_meta {
        haplotaggedBam: "Guppy with Remora reads aligned to assembly. in BAM format."
        ref: "Assembly (reference) to that reads are aligned to."
        sample_name: "Sample name. Will be used in output bed file."
        ref_name: "Reference name. Will be used in output bed file."
        #modType: "Modified base of interest, one of: 5mC (default), 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao."
        out_type_filter: " Output records filtered to these sites: cpg (default), chg, chh"
    }

    String sample_name_ref = "~{sample_name}"+"_"+"~{ref_name}"
    Boolean regionaly =  defined(regional_bed)

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## run a recurrent "top" in the background to monitor resource usage
        if [ ~{resourceLogScript} != "" ]
        then
            bash ~{resourceLogScript} 20 top.log &
        fi

        ln -s ~{haplotaggedBam} reads.bam
        ln -s ~{haplotaggedBamBai} reads.bam.bai


        mkdir modkit_out

        # note output format: https://nanoporetech.github.io/modkit/intro_bedmethyl.html
        # filtering threshold default 10-th percentile of calls https://github.com/nanoporetech/modkit/blob/master/filtering.md
        if [ ~{regionaly} == true ]
        then
            modkit pileup ~{out_type_filter} ~{partitionTag} --prefix ~{sample_name_ref} --ref ~{ref} \
                   --threads ~{threadCount} --include-bed ~{regional_bed} ~{extraArgs} reads.bam modkit_out
        else
            # modkit command with reference on input reads
            modkit pileup ~{out_type_filter} ~{partitionTag} --prefix ~{sample_name_ref} --ref ~{ref} --threads ~{threadCount} ~{extraArgs} reads.bam modkit_out
        fi
        bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_1.bed
        bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_2.bed
        bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_ungrouped.bed

    >>>

        output {
            File hap1bedOut      = "modkit_out/~{sample_name_ref}_1.bed.gz"
            File hap2bedOut      = "modkit_out/~{sample_name_ref}_2.bed.gz"
            File? ungroupedBedOut      = "modkit_out/~{sample_name_ref}_ungrouped.bed.gz"
            File? toplog = "top.log"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}
