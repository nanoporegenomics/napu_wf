version 1.0



workflow pilotUpload {
    input {
        #File modkitBEDunphased
        #File modkitBed_1
        #File modkitBed_2
        #File modkitBed_ungrouped

        #File snifflesvcf
        #File snifflessnf
        #File hapdiff09
        #File harmonizedVCF
        #File dvGvcf

        File haplotaggedBAM
        File mappedBam

        File altchroms_file
        String sample 
        String cohortnum
        String staging_gs_bucket


    }

    call indexMergeUpload {
        input:
            unphasedMappedBAM=mappedBam,
            haplotaggedBAM=haplotaggedBAM,
            altchroms_file=altchroms_file,
            sample=sample,
            staging_gs_bucket=staging_gs_bucket
    }

    #call uploadStagingData {
    #    input:
    #        modkitBEDunphased = modkitBEDunphased,
    #        modkitBed_1 = modkitBed_1,
    #        modkitBed_2 = modkitBed_2,
    #        modkitBed_ungrouped = modkitBed_ungrouped,
    #        harmonizedVCF = harmonizedVCF, 
    #        snifflesvcf = snifflesvcf,
    #        snifflessnf = snifflessnf,
    #        hapdiff09 = hapdiff09,
    #        dvGvcf = dvGvcf,
    #        sample=sample,
    #        staging_gs_bucket=staging_gs_bucket
    #}
    

    output {
        #File trackerfileAll = uploadStagingData.trackerfile
        File readcount = indexMergeUpload.readcount
        File trackerfileBAM = indexMergeUpload.trackerfile
        File origionalReadsPerChro = indexMergeUpload.origionalReadsPerChr
        File uploadedReadsPerChro = indexMergeUpload.uploadedReadsPerChr
        File numReadFileo = indexMergeUpload.numReadFile

    }
}

 
task uploadStagingData{
    input {
        File modkitBEDunphased
        File modkitBed_1
        File modkitBed_2
        File modkitBed_ungrouped

        File harmonizedVCF
        File snifflesvcf
        File snifflessnf
        File hapdiff09
        File dvGvcf 

        String sample 
        String staging_gs_bucket

        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 2 * (round(size(dvGvcf, "GB")) + round(size(hapdiff09, "GB")) + 80 )

    }

    command <<<

        set -eux -o pipefail

        # gsutil cp commands with file name changes if necessary; for example:
        # gsutil cp file "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam

        # methylBED unphased
        gsutil cp ~{modkitBEDunphased} "~{staging_gs_bucket}"/data_files/"~{sample}"/methylation/
        echo -e "~{sample}.modkit_unphased\tWGS\tbed\t#\t#\t~{sample}_GRCh38.bed.gz\t~{staging_gs_bucket}/data_files/~{sample}/methylation/~{sample}_GRCh38.bed.gz\tGRCh38\n" >> trackerfile.txt

        # methylBED hap1
        if [[ -f "~{modkitBed_1}" ]]; then
            gsutil cp ~{modkitBed_1} "~{staging_gs_bucket}"/data_files/"~{sample}"/methylation/
            echo -e "~{sample}.modkit_hap1\tWGS\tbed\t#\t#\t~{sample}_GRCh38_1.bed.gz\t~{staging_gs_bucket}/data_files/~{sample}/methylation/~{sample}_GRCh38_1.bed.gz\tGRCh38\n" >> trackerfile.txt
        fi

        # methylBED hap2
        if [[ -f "~{modkitBed_2}" ]]; then
            gsutil cp ~{modkitBed_2} "~{staging_gs_bucket}"/data_files/"~{sample}"/methylation/
            echo -e "~{sample}.modkit_hap2\tWGS\tbed\t#\t#\t~{sample}_GRCh38_2.bed.gz\t~{staging_gs_bucket}/data_files/~{sample}/methylation/~{sample}_GRCh38_2.bed.gz\tGRCh38\n" >> trackerfile.txt
        fi

        # methylBED ungrouped phased
        if [[ -f "~{modkitBed_ungrouped}" ]]; then
            gsutil cp ~{modkitBed_ungrouped} "~{staging_gs_bucket}"/data_files/"~{sample}"/methylation/
            echo -e "~{sample}.modkit_ungrouped\tWGS\tbed\t#\t#\t~{sample}_GRCh38_ungrouped.bed.gz\t~{staging_gs_bucket}/data_files/~{sample}/methylation/~{sample}_GRCh38_ungrouped.bed.gz\tGRCh38\n" >> trackerfile.txt
        fi

        #variant_calls
        # harmonized svs snvs
        gsutil cp ~{harmonizedVCF} "~{staging_gs_bucket}"/data_files/"~{sample}"/variant_calls/"~{sample}".harmonized.phased.vcf.gz
        echo -e "~{sample}_harmonizedGT_svs_snvs\tWGS\tvcf\t#\t#\t~{sample}.harmonized.phased.vcf.gz\t~{staging_gs_bucket}/data_files/~{sample}/variant_calls/~{sample}.harmonized.phased.vcf.gz\tGRCh38\n" >> trackerfile.txt

        # sniffles vcf
        bgzip ~{snifflesvcf} | gsutil cp - "~{staging_gs_bucket}"/data_files/"~{sample}"/variant_calls/"~{sample}".sniffles2.2.3.grch38.vcf.gz
        echo -e "~{sample}_sniffles_vcf\tWGS\tvcf\t#\t#\t~{sample}.sniffles2.2.3.grch38.vcf.gz\t~{staging_gs_bucket}/data_files/~{sample}/variant_calls/~{sample}.sniffles2.2.3.grch38.vcf.gz\tGRCh38\n" >> trackerfile.txt

        # sniffle snf
        gsutil cp ~{snifflessnf} "~{staging_gs_bucket}"/data_files/"~{sample}"/variant_calls/"~{sample}".sniffles2.2.3.grch38.snf
        echo -e "~{sample}_sniffles_snf\tWGS\tsnf\t#\t#\t~{sample}.sniffles2.2.3.grch38.snf\t~{staging_gs_bucket}/data_files/~{sample}/variant_calls/~{sample}.sniffles2.2.3.grch38.snf\tGRCh38\n" >> trackerfile.txt

        # hapdiff 0.9 : structuralVariantsAsmVcf_09_GRCh38
        gsutil cp ~{hapdiff09} "~{staging_gs_bucket}"/data_files/"~{sample}"/variant_calls/"~{sample}".hapdiff0.9_unphased.grch38.vcf.gz
        echo -e "~{sample}_hapdiff_vcf\tWGS\tvcf\t#\t#\t~{sample}.hapdiff0.9_unphased.grch38.vcf.gz\t~{staging_gs_bucket}/data_files/~{sample}/variant_calls/~{sample}.hapdiff0.9_unphased.grch38.vcf.gz\tGRCh38\n" >> trackerfile.txt

        # gvcf
        gsutil cp ~{dvGvcf} "~{staging_gs_bucket}"/data_files/"~{sample}"/variant_calls/
        echo -e "~{sample}_PMDV_gvcf\tWGS\tvcf\t#\t#\t~{sample}_PMDV.g.vcf.gz\t~{staging_gs_bucket}/data_files/~{sample}/variant_calls/~{sample}_PMDV.g.vcf.gz\tGRCh38\n" >> trackerfile.txt

        
    >>>



    output {
        File trackerfile = "trackerfile.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }
}

task indexMergeUpload {


    input {
        File unphasedMappedBAM
        File haplotaggedBAM
        File altchroms_file
        String sample
        String staging_gs_bucket
        String filesuffix = ""
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 3 * round(size(unphasedMappedBAM, "GB")) + 80
    }

    String outname = "~{sample}"+".GRCh38.bam"
    String tracker_string = "~{sample}.GRCh38.bam\tWGS\tBAM\t#\t#\t~{sample}.GRCh38.bam\t~{staging_gs_bucket}/data_files/~{sample}/reads/~{sample}.GRCh38.bam\tGRCh38\n"


    command <<<

        set -eux -o pipefail

        #1 : get the header from the unphased bam into a tmp.sam to append the reads to
        samtools view -H ~{unphasedMappedBAM} > tmp.extracted_reads.sam

        # 2: get alts from unphasedMappedBAM array append reads (no header) to the tmp sam
        echo "indexing and extracting "
        samtools index -@ ~{threads} ~{unphasedMappedBAM}

        # append the alt reads to tmp sam for easy concatination 
        samtools view -@ ~{threads} ~{unphasedMappedBAM} $(cat ~{altchroms_file}) >> tmp.extracted_reads.sam

        # 3: isolate and combine unmapped reads 
        echo "find unmapped reads"
        samtools view -f 4 -@ ~{threads} ~{unphasedMappedBAM} >> tmp.extracted_reads.sam
        echo "unmapped reads ~{unphasedMappedBAM}" >> readcount.txt
        samtools view -c -f 4 -@ ~{threads} ~{unphasedMappedBAM} >> readcount.txt

        # 4: convert the tmp.sam to a bam
        samtools view -b -@ ~{threads} tmp.extracted_reads.sam | samtools sort -@ ~{threads} - > tmp.alt_reads.bam

        # 5: merge the haplotagged BAM with the alts and unmapped reads
        samtools merge -@ ~{threads} -o ~{outname} ~{haplotaggedBAM} tmp.alt_reads.bam

        # 6: index the merged BAM
        samtools index -@ ~{threads} ~{outname}

        # 7: move to staging workspace
        #gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/
        gsutil cp ~{outname} "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam
        gsutil cp ~{outname}.bai "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/"~{sample}".GRCh38.bam.bai
        gsutil ls "~{staging_gs_bucket}"/data_files/"~{sample}"/reads/

        samtools view -@ ~{threads} -c ~{outname} >> readcount.txt

        echo $(echo ~{tracker_string}) > trackerfile.txt

        # check reads in 
        samtools view ~{unphasedMappedBAM} | cut -f1 | sort | uniq -c >> ~{sample}.readsInUnphasedBam.chr.txt

        echo "unphased_aln _unphased_read count" > ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUnphasedBam.chr.txt >> ~{sample}.numReads.txt

        samtools view -@ ~{threads} ~{outname} | cut -f1 | sort | uniq -c > ~{sample}.readsInUploadedBam.chr.txt

        echo "uploaded_aln uploaded_read count" >> ~{sample}.numReads.txt
        awk '{sum += $1; count++} END {print sum, count}' ~{sample}.readsInUploadedBam.chr.txt >> ~{sample}.numReads.txt

        echo "uploaded_chr_aln_count" >> ~{sample}.numReads.txt
        samtools view -@ ~{threads} ~{outname} | cut -f3 | sort | uniq -c >> ~{sample}.numReads.txt

    >>>

    output {
        String trackerString = "~{tracker_string}"
        File readcount = "readcount.txt"
        File trackerfile = "trackerfile.txt"
        File origionalReadsPerChr = "~{sample}.readsInUnphasedBam.chr.txt"
        File uploadedReadsPerChr = "~{sample}.readsInUploadedBam.chr.txt"
        File numReadFile = "~{sample}.numReads.txt"

    }

    runtime {
        #preemptible: 2
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/gsutilsamtools@sha256:d649c4eb695a32aa15e78bdc7a95ef81ae88b741f9ce10e40831566342e13480"
    }
}

