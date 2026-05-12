version 1.0



workflow combineUnmapped {
    input {
        File phasedMappedBAM 
        File phasedMappedBAI
        File unphasedMappedBAM
        File unphasedMappedBAI
        String sample 
        String cohortnum
        #String staging_gs_bucket
        #Boolean findUnmapped = false
        #Int diskSizeGB = 1024 # 5 * round(size(reads, "G")) + 20
        #Int memSizeGb = 128

    }


    if (cohortnum=="3"){
        call combine_unmapped_coh3 {
            input:
                phasedBAM=phasedMappedBAM,
                phasedBAI=phasedMappedBAI,
                unphasedMappedBAM=unphasedMappedBAM,
                unphasedMappedBAI=unphasedMappedBAI,
                sample=sample
        }
    }
    

    String tstring = select_first([indexMergeUpload.trackerString, indexMergeUpload_coh2.trackerString])
    File rcFile = select_first([indexMergeUpload.readcount, indexMergeUpload_coh2.readcount])
    File tfile = select_first([indexMergeUpload.trackerfile, indexMergeUpload_coh2.trackerfile])

    output {
        File readcount = combine_unmapped_coh3.readcount
        File outPhasedBAM = combine_unmapped_coh3.outPhasedBAM

    }
}

task combine_unmapped_coh3 {
    
    input {
        File unphasedMappedBAM
        File unphasedMappedBAI
        File phasedBAM
        File phasedBAI
        String sample
        Int memSizeGB = 80
        Int threads = 16
        Int diskSizeGB = 2 * round(size(unphasedMappedBAM, "GB")) + 40
    }

    String outname = "~{sample}"+".GRCh38.bam"

    command <<<

        set -eux -o pipefail

        #1 : get the header from the first unphased bam into a tmp.sam to append all the reads to
        samtools view -H ~{unphasedMappedBAM} > tmp.extracted_reads.sam

                # 2: get alts from unphasedMappedBAMs array append reads (no header) to the tmp sam
                #echo "indexing and extracting "

                # append the alt reads to tmp sam for easy concatination 
                #samtools view -@ ~{threads} ~{unphasedMappedBAM} $(cat ~{altchroms_file}) >> tmp.extracted_reads.sam

        # 3: get unmapped reads from bam
        echo "find unmapped reads"
        samtools view -f 4 -@ ~{threads} ~{unphasedMappedBAM} >> tmp.extracted_reads.sam
        echo "unmapped reads ~{unphasedMappedBAM}" >> readcount.txt
        samtools view -c -f 4 -@ ~{threads} ~{unphasedMappedBAM} >> readcount.txt

        # 4: convert the tmp.sam to a bam
        samtools view -b -@ ~{threads} tmp.extracted_reads.sam | samtools sort -@ ~{threads} - > tmp.extracted_reads.bam

        # 5: merge unmapped and haplotagged bams
        samtools merge -@ ~{threads} -o ~{outname} ~{phasedBAM} tmp.extracted_reads.bam

        # 6: index the merged BAM
        samtools index -@ ~{threads} ~{outname}

    >>>

    output {
        File outPhasedBAM = "~{outname}"
        File outPhasedBAI = "~{outname}".bai
        File readcount = "readcount.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }

}

