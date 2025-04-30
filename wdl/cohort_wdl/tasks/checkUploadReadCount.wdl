version 1.0



workflow combineAndUpload {
    input {
        Array[File] unphasedMappedBAMs = []
        File? unphasedMappedBAM
        File? unphasedMappedBAI
        String sample 
        String cohortnum
        String staging_gs_bucket
        #Int diskSizeGB = 1024 # 5 * round(size(reads, "G")) + 20
        #Int memSizeGb = 128

    }

    if (cohortnum=="1"){
        call checkUploadedReads {
        input:
        unphasedMappedBAMs=unphasedMappedBAMs,
        sample=sample,
        staging_gs_bucket=staging_gs_bucket

        }
    }

    if (cohortnum=="2"){
        call checkUploadedReads_coh2 {
        input:
        unphasedMappedBAM=unphasedMappedBAM,
        unphasedMappedBAI=unphasedMappedBAI,
        sample=sample,
        staging_gs_bucket=staging_gs_bucket

        }
    }
    

    File origionalReadsPerChr = select_first([checkUploadedReads.origionalReadsPerChr, checkUploadedReads_coh2.origionalReadsPerChr])
    File uploadedReadsPerChr = select_first([checkUploadedReads.uploadedReadsPerChr, checkUploadedReads_coh2.uploadedReadsPerChr])
    File readcount = select_first([checkUploadedReads.readcount, checkUploadedReads_coh2.readcount])

    output {
        File origionalReadsPerChrf = origionalReadsPerChr
        File uploadedReadsPerChrf = uploadedReadsPerChr
        File readcountf = readcount

    }
}


task checkUploadedReads_coh2 {


    input {
        File? unphasedMappedBAM
        File? unphasedMappedBAI
        String sample
        String staging_gs_bucket
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = round(size(unphasedMappedBAM, "GB")) + 40
    }


    command <<<

        set -eux -o pipefail

        #1 : count reads in the unphased bam
        samtools view -@ ~{threads} ~{unphasedMappedBAM} | cut -f1 | sort | uniq -c > ~{sample}.readsInUnphasedBam.chr.txt

        #2 : number of reads in unphased bam
        echo "unphased_aln _unphased_read count" > ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUnphasedBam.chr.txt >> ~{sample}.numReads.txt

        #3 : count reads in the uploaded bam
        uploadedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.haplotagged.bam" 
        gsutil cat ${uploadedBAM} | samtools view -@ ~{threads} | cut -f1 | sort | uniq -c > ~{sample}.readsInUploadedBam.chr.txt

        #2 : number of reads in unphased bam
        echo "uploaded_aln uploaded_read count" >> ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUploadedBam.chr.txt >> ~{sample}.numReads.txt


    >>>

    output {
        File origionalReadsPerChr = "~{sample}.readsInUnphasedBam.chr.txt"
        File uploadedReadsPerChr = "~{sample}.readsInUploadedBam.chr.txt"
        File readcount = "~{sample}.numReads.txt"

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

task checkUploadedReads {


    input {
        Array[File] unphasedMappedBAMs
        String sample
        String staging_gs_bucket
        Boolean findUnmapped
        String filesuffix = ".GRCh38"
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 2 * round(size(unphasedMappedBAMs, "GB")) + 40
    }


    command <<<

        set -eux -o pipefail


        # 1: count reads in the unphased bams
        for bam in ~{sep=' ' unphasedMappedBAMs}; do

            # 2 : count reads in the unphased bam
            samtools view ${bam} | cut -f1 | sort | uniq -c >> ~{sample}.readsInUnphasedBam.chr.txt


        done

        #3 : number of reads in unphased bam
        echo "unphased_aln _unphased_read count" > ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUnphasedBam.chr.txt >> ~{sample}.numReads.txt

        #4 : count reads in the uploaded bam
        uploadedBAM="~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}~{filesuffix}.bam" 
        gsutil cat ${uploadedBAM} | samtools view -@ ~{threads} | cut -f1 | sort | uniq -c > ~{sample}.readsInUploadedBam.chr.txt

        #5 : number of reads in unphased bam
        echo "uploaded_aln uploaded_read count" >> ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUploadedBam.chr.txt >> ~{sample}.numReads.txt

    >>>

        output {
        File origionalReadsPerChr = "~{sample}.readsInUnphasedBam.chr.txt"
        File uploadedReadsPerChr = "~{sample}.readsInUploadedBam.chr.txt"
        File readcount = "~{sample}.numReads.txt"

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