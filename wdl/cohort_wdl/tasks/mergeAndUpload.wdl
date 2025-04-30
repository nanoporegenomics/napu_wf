version 1.0



workflow combineAndUpload {
    input {
        Array[File] unphasedMappedBAMs = []
        File? unphasedMappedBAM
        File? unphasedMappedBAI
        File altchroms_file
        String sample 
        String cohortnum
        String staging_gs_bucket
        Boolean findUnmapped = false
        #Int diskSizeGB = 1024 # 5 * round(size(reads, "G")) + 20
        #Int memSizeGb = 128

    }

    if (cohortnum=="1"){
        call indexMergeUpload {
        input:
        unphasedMappedBAMs=unphasedMappedBAMs,
        altchroms_file=altchroms_file,
        sample=sample,
        staging_gs_bucket=staging_gs_bucket,
        findUnmapped=findUnmapped

        }
    }

    if (cohortnum=="2"){
        call indexMergeUpload_coh2 {
        input:
        unphasedMappedBAM=unphasedMappedBAM,
        unphasedMappedBAI=unphasedMappedBAI,
        altchroms_file=altchroms_file,
        sample=sample,
        staging_gs_bucket=staging_gs_bucket

        }
    }
    

    String tstring = select_first([indexMergeUpload.trackerString, indexMergeUpload_coh2.trackerString])
    File rcFile = select_first([indexMergeUpload.readcount, indexMergeUpload_coh2.readcount])
    File tfile = select_first([indexMergeUpload.trackerfile, indexMergeUpload_coh2.trackerfile])

    output {
        String trackerString = tstring
        File readcount = rcFile
        File trackerfile = tfile

    }
}


task indexMergeUpload_coh2 {


    input {
        File? unphasedMappedBAM
        File? unphasedMappedBAI
        File altchroms_file
        String sample
        String staging_gs_bucket
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 3 * round(size(unphasedMappedBAM, "GB")) + 40
    }

    String outname = "~{sample}"+".GRCh38.bam"
    String tracker_string = "~{sample}.GRCh38.bam\tWGS\tBAM\t#\t#\t~{sample}.GRCh38.bam\t~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.GRCh38.bam\tGRCh38\n"
    

    command <<<

        set -eux -o pipefail

        #1 : get the header from the first unphased bam into a tmp.sam to append all thre reads to
        samtools view -H ~{unphasedMappedBAM} > tmp.extracted_reads.sam

        # 2: get alts from unphasedMappedBAMs array append reads (no header) to the tmp sam
        echo "indexing and extracting "

        # append the alt reads to tmp sam for easy concatination 
        samtools view -@ ~{threads} ~{unphasedMappedBAM} $(cat ~{altchroms_file}) >> tmp.extracted_reads.sam


        # 3: convert the tmp.sam to a bam
        samtools view -b -@ ~{threads} tmp.extracted_reads.sam | samtools sort -@ ~{threads} - > tmp.alt_reads.bam

        # 4: get haplotagged bam, and unmapped bam
        phasedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.haplotagged.bam" 
        unmappedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.unmappedGRCh38.bam"

        # 5: merge the unmapped BAM with the alts
        gsutil cat ${unmappedBAM} | samtools merge -o tmp.~{sample}.unmapped.alts.bam - tmp.alt_reads.bam

        # 5: merge the haplotagged BAM with the unmapped.alts
        gsutil cat ${phasedBAM} | samtools merge -@ ~{threads} -o - - tmp.~{sample}.unmapped.alts.bam | samtools sort -@~{threads} - > ~{outname}


        # 5: index the merged BAM
        samtools index -@ ~{threads} ~{outname}

        # 6: move to staging workspace
        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/
        gsutil cp ~{outname} "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam
        gsutil cp ~{outname}.bai "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam.bai
        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/

        samtools view -@ ~{threads} -c ~{outname} > readcount.txt

        echo $(echo ~{tracker_string}) > trackerfile.txt


    >>>

    output {
        String trackerString = "~{tracker_string}"
        File readcount = "readcount.txt"
        File trackerfile = "trackerfile.txt"
    }

    runtime {
        #preemptible: 2
        #time: 240
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }
}

task indexMergeUpload {


    input {
        Array[File] unphasedMappedBAMs
        File altchroms_file
        String sample
        String staging_gs_bucket
        Boolean findUnmapped
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 2 * round(size(unphasedMappedBAMs, "GB")) + 40
    }

    String outname = "~{sample}"+".GRCh38.bam"
    String tracker_string = "~{sample}.GRCh38.bam\tWGS\tBAM\t#\t#\t~{sample}.GRCh38.bam\t~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.GRCh38.bam\tGRCh38\n"
    File firstReadFile = select_first(unphasedMappedBAMs)

    command <<<

        set -eux -o pipefail

        #1 : get the header from the first unphased bam into a tmp.sam to append all thre reads to
        samtools view -H ~{firstReadFile} > tmp.extracted_reads.sam

        # 2: get alts from unphasedMappedBAMs array append reads (no header) to the tmp sam
        for bam in ~{sep=' ' unphasedMappedBAMs}; do
            echo "indexing and extracting ${bam}"

            # index the unphaed BAM
            samtools index -@ ~{threads} ${bam}

            # append the alt reads to tmp sam for easy concatination 
            samtools view -@ ~{threads} ${bam} $(cat ~{altchroms_file}) >> tmp.extracted_reads.sam


            if [[ ~{findUnmapped} == "true" ]]
            then
                echo "find unmapped reads"
                samtools view -f 4 -@ ~{threads} ${bam} >> tmp.extracted_reads.sam
            fi

        done


        # 3: convert the tmp.sam to a bam
        samtools view -b -@ ~{threads} tmp.extracted_reads.sam | samtools sort -@ ~{threads} - > tmp.alt_reads.bam

        # 4: get haplotagged bam, and unmapped bam
        phasedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.bam"

        if [[ ~{findUnmapped} == "false" ]]
        then
            unmappedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.unmappedGRCh38.bam"

            # 5: merge the unmapped BAM with the alts
            gsutil cat ${unmappedBAM} | samtools merge -o tmp.~{sample}.unmapped.alts.bam - tmp.alt_reads.bam

            # 6: merge the haplotagged BAM with the unmapped.alts
            gsutil cat ${phasedBAM} | samtools merge -@ ~{threads} -o - - tmp.~{sample}.unmapped.alts.bam | samtools sort -@~{threads} - > ~{outname}

        else
            # 6: merge the haplotagged BAM with the alts
            gsutil cat ${phasedBAM} | samtools merge -@ ~{threads} -o - - tmp.alt_reads.bam | samtools sort -@~{threads} - > ~{outname}
        fi
        
        
        # 5: index the merged BAM
        samtools index -@ ~{threads} ~{outname}

        # 6: move to staging workspace
        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/
        gsutil cp ~{outname} "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam
        gsutil cp ~{outname}.bai "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam.bai
        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/

        samtools view -@ ~{threads} -c ~{outname} > readcount.txt

        echo $(echo ~{tracker_string}) > trackerfile.txt


    >>>

    output {
        String trackerString = "~{tracker_string}"
        File readcount = "readcount.txt"
        File trackerfile = "trackerfile.txt"
    }

    runtime {
        #preemptible: 2
        #time: 240
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }
}