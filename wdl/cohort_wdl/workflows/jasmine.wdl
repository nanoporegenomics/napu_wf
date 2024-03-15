version 1.0

workflow run_jasmine{
    input {
        Array[File] vcfFiles1 = []
        Array[File] vcfFiles2 = []
        File? regional_bed
        String out_name
        Int preemptibles = 0
        String? extraArgs = ""
        String dockerImage = "meredith705/jasmine:latest"

    }


    call jasmine_merge as jasmine{
        input:
					vcfFiles1 = vcfFiles1,
                    vcfFiles2 = vcfFiles2,
					out_name = out_name,
					regional_bed = regional_bed,
                    dockerImage = dockerImage,
                    preemptibles = preemptibles
    }

    output{
        File? merged_sv_vcf = jasmine.merged_sv_vcf
        File? merged_sv_vcf_tbi = jasmine.merged_sv_vcf_tbi
    }
}

task jasmine_merge {
    input {
        Array[File] vcfFiles1 = []
        Array[File] vcfFiles2 = []
        File? regional_bed
        String out_name
        String? extraArgs = ""
        Int specReads = 10
        Int specLen = 30
        Int memSizeGB = 64
        Int threadCount = 64
        Int preemptibles = 0
        Int diskSizeGB = 4 * round(size(vcfFiles, 'G')) + 800
        String dockerImage = "meredith705/jasmine:latest"

    }

    Boolean regionaly =  defined(regional_bed)
    # combine vcfs from multiple SV callers
    Array[File] vcfFiles = flatten( [vcfFiles1, vcfFiles2] )


    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # pipeline for Jasmine
        # https://github.com/mkirsche/Jasmine/tree/master/pipeline

        # Normalize SV types in each sample
        jasmine --preprocess_only --pre_normalize file_list=~{write_lines(vcfFiles)}

        # Mark high-confidence callset (high-specificity callset) in each sample
        jasmine file_list=~{write_lines(vcfFiles)} --preprocess_only --mark_specific \
            spec_reads=~{specReads} spec_len=~{specLen}

        # make an empty file to store intra sample merged vcfs
        touch files.txt
        for vcf in ~{sep=" " select_all(vcfFiles)}
        do
            echo $vcf
            echo
            jasmine file_list=$vcf max_dist=200 \
                --allow_intrasample out_file=$vcf.intraSampleMerged.vcf \
                --nonlinear_dist --comma_filelist
            echo $vcf.intraSampleMerged.vcf >> files.txt
        done

        #### Merge SVs across samples ####
        echo starting merge
        jasmine --output_genotypes file_list=files.txt out_file=~{out_name}.cohort.merged.vcf

        # Convert insertions back to duplications
        echo convert insertions
        jasmine --dup_to_ins --postprocess_only out_file=~{out_name}.cohort.merged.vcf

        # Remove low-confidence or imprecise calls:
        cat ~{out_name}.cohort.merged.vcf | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' \
            | bcftools sort | bgzip -@ ~{threadCount} > ~{out_name}.cohort.merged.conf.vcf.gz

        tabix ~{out_name}.cohort.merged.conf.vcf.gz


    >>>

        output {
            File? merged_sv_vcf = "~{out_name}.cohort.merged.conf.vcf.gz"
            File? merged_sv_vcf_tbi = "~{out_name}.cohort.merged.conf.vcf.gz.tbi"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: preemptibles
    }
}
