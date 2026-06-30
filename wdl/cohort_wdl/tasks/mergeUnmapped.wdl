version 1.0



workflow combineUnmapped {
    input {
        File phasedMappedBAM 
        File phasedMappedBAI
        File unphasedMappedBAM
        File unphasedMappedBAI
        String sample 
        #Int diskSizeGB = 1024 # 5 * round(size(reads, "G")) + 20
        #Int memSizeGb = 128

    }


    call combine_unmapped {
        input:
            phasedBAM=phasedMappedBAM,
            phasedBAI=phasedMappedBAI,
            unphasedMappedBAM=unphasedMappedBAM,
            unphasedMappedBAI=unphasedMappedBAI,
            sample=sample
    }

    output {
        File readcount = combine_unmapped.readcount
        File outPhasedBAM = combine_unmapped.outPhasedBAM
        File outPhasedBAI = combine_unmapped.outPhasedBAI
        Int  unmapped_reads = combine_unmapped.unmapped_reads
    }
}

task combine_unmapped {
    
    input {
        File unphasedMappedBAM
        File unphasedMappedBAI
        File phasedBAM
        File phasedBAI
        String sample
        Int memSizeGB = 2 * round(size(unphasedMappedBAM, "GB"))
        Int threads = 16
        Int diskSizeGB = 2 * round(size(unphasedMappedBAM, "GB")) + 40
    }

    String outname = "~{sample}"+".haplotagged.GRCh38.bam"

    command <<<

        set -eux -o pipefail

        #1 : get the header from the first unphased bam into a tmp.sam to append all the reads to
        samtools view -H ~{unphasedMappedBAM} > tmp.extracted_reads.sam

        # 2: get unmapped reads from bam
        echo "find unmapped reads"
        samtools view -f 4 -@ ~{threads} ~{unphasedMappedBAM} >> tmp.extracted_reads.sam
                #echo "unmapped reads ~{unphasedMappedBAM}" >> readcount.txt
        UNMAPPED=$(samtools view -c -f 4 -@ ~{threads} ~{unphasedMappedBAM})
        echo "Unmapped reads: ${UNMAPPED}"
        echo "Unmapped reads: $UNMAPPED" >&2
        echo "$UNMAPPED" > readcount.txt

        # 3: convert the tmp.sam to a bam
        samtools view -b -@ ~{threads} tmp.extracted_reads.sam | samtools sort -@ ~{threads} - > tmp.extracted_reads.bam

        # 4: merge unmapped and haplotagged bams and sort
        # samtools merge -@ ~{threads} -o ~{outname} ~{phasedBAM} tmp.extracted_reads.bam
        samtools merge -@ ~{threads} -u - ~{phasedBAM} tmp.extracted_reads.bam | \
          samtools sort -@ ~{threads} -o ~{outname} -

        # 5: index the merged BAM
        samtools index -@ ~{threads} ~{outname}

    >>>

    output {
        File outPhasedBAM = "~{outname}"
        File outPhasedBAI = "~{outname}.bai"
        File readcount = "readcount.txt"
        Int  unmapped_reads  = read_int("readcount.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }

}

