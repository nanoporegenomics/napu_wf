version 1.0



workflow deleteStagedData {
    input {
        String sample 
        String cohortnum
        String staging_gs_bucket
    }

    if (cohortnum=="1"){
        call deletehaplotagged_coh1 {
        input:
        sample=sample,
        staging_gs_bucket=staging_gs_bucket
        }
    }

    if (cohortnum=="2"){
        call deletehaplotagged_coh2 {
        input:
        sample=sample,
        staging_gs_bucket=staging_gs_bucket

        }
    }

    File outfilet = select_first([deletehaplotagged_coh1.outfile, deletehaplotagged_coh2.outfile])

    output {
        File outfile = outfilet
    }    

}


task deletehaplotagged_coh2 {
    input {
        String sample
        String staging_gs_bucket
        Int memSizeGB = 2
        Int diskSizeGB = 2
        Int threads = 1
    }

    command <<<

        set -eux -o pipefail

        # 1: get haplotagged bam, and unmapped bam
        phasedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.haplotagged.bam" 
        unmappedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.unmappedGRCh38.bam"

        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/

        # delete the .haplotagged.bam
        gsutil rm ${phasedBAM}

        # delete the .unmappedGRCh38.bam
        gsutil rm ${unmappedBAM}

        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/ > ~{sample}.outfile.txt
    >>>

    output {
        File outfile = "~{sample}.outfile.txt"
        }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }

}    


task deletehaplotagged_coh1 {
    input {
        String sample
        String staging_gs_bucket
        String filesuffix = ""
        Int memSizeGB = 2
        Int diskSizeGB = 2
        Int threads = 1
    }

    command <<<

        set -eux -o pipefail

        # 1: get haplotagged bam, and unmapped bam : all samples will have this
        phasedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}~{filesuffix}.bam"
        unmappedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.unmappedGRCh38.bam"

        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/

        # delete the .haplotagged.bam
        gsutil rm ${phasedBAM}

        # delete the .unmappedGRCh38.bam
        gsutil rm ${unmappedBAM}

        # some samples will have chrM
        chrmBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.chrM_GRCh38.bam"

        # check if chrM exists: 
        if gsutil ls ${chrmBAM} > /dev/null 2>&1; then
            echo "ChrM file exists"
            # delete the .chrM.bam
            gsutil rm ${chrmBAM}
        else
            echo "chrM does not exist"
        fi

        gsutil ls -lh "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/ > ~{sample}.outfile.txt
    >>>

    output {
        File outfile = "~{sample}.outfile.txt"
        }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }

}    







