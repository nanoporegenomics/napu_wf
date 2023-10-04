version 1.0

import "../tasks/pepper-margin-dv.wdl" as pmdv_haplotag_t
import "../tasks/minimap2.wdl" as minimap_t


workflow pepper_margin_dv_chrs {

    input {
        Array[File] bamChrs
        Array[File] bamChrsIdx
        File referenceFile
        Boolean gvcf = false
        String sampleName = ""
        Int threads
        Int preemptible = 0
    }

    # run pmdv on chromosome chunks
    scatter (bamChr in zip(bamChrs, bamChrsIdx)){
        call pmdv_haplotag_t.pepper_margin_dv_t as pmdvHap_chrs{
            input:
                threads = threads,
                gvcf = gvcf,
                reference = referenceFile,
                bamAlignment = bamChr.left,
                bamAlignmentIndex = bamChr.right,
                sampleName = sampleName
            }
        }

    # merge the per chr bams
    call minimap_t.mergeBAM as merge_PEPPER_DV_BAMs {
        input:
            bams = select_all(pmdvHap_chrs.haplotaggedBam),
            outname = sampleName
        }

    # merge the per chr vcfs
    call pmdv_haplotag_t.mergeVCFs{
        input:
            vcfFiles = pmdvHap_chrs.pepperVcf,
            outname = sampleName
        }

    if (gvcf){
        # merge the per chr g.vcfs
        call pmdv_haplotag_t.mergeGVCFs as mergeGVCFs{
            input:
                vcfFiles = select_all(pmdvHap_chrs.pepperGVcf),
                outname = sampleName
            }
    }

    output {
        File pmdv_merged_vcf = mergeVCFs.vcf
        File pmdv_merged_vcf_Idx = mergeVCFs.vcfIndex
        File? pepper_merged_GVcf = mergeGVCFs.gvcf
        File? pepper_merged_GVcf_Idx = mergeGVCFs.gvcfIndex
        File? phasedBam = merge_PEPPER_DV_BAMs.bam
        File? phasedBamBai = merge_PEPPER_DV_BAMs.bamIndex
    }
}