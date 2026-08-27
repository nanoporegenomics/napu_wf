version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/dv-margin.wdl" as dv_margin_t

workflow dvMargin {

    input {
        File referenceFile
        File? bamAlignment
        File? bamAlignmentIndex
        File? fastqReads
        String sampleName = "sample"
        Array[String] chrs = []
        Boolean phaseVariants = false
        Int threads
    }

    if (defined(fastqReads)){
        call minimap_t.minimap2_t as mm_align {
            input:
                reads = select_first([fastqReads]),
                reference = referenceFile,
                threads = threads
        }        
    }

    File bamFile = select_first([bamAlignment, mm_align.bam])
    File bamIdxFile = select_first([bamAlignmentIndex, mm_align.bamIndex])

    # if the chr arr is provided, split the bams prior to running DV
    # adding this for running DV/chrom outside of the end2end pipeline
    if(length(chrs) > 0 {
        call minimap_t.indexBAM as indexSingleInputBam{
            input: 
                bam = inputBam,
                chrs = chrs
        }

        # don't actually need select_first as there is only one option for split bams
        Array[File] bamChrs = select_first([indexSingleInputBam.bamPerChrs])
        Array[File] bamChrsIndex = select_first([indexSingleInputBam.bamPerChrsIndex])

        # run dv on each chr bam
        scatter (bamChr in zip(bamChrs, bamChrsIndex)){
            call dv_margin_t.dv_t as chr_dv_t {
                input:
                    threads = threads,
                    reference = referenceFasta,
                    bamAlignment = bamChr.left,
                    bamAlignmentIndex = bamChr.right,
                    sampleName = sampleName,
                    oneChr = true,
                    preemptible = 2
            }
        }
        # combine the vcf's and gvcf's from each bam run
        call dv_margin_t.mergeVCFs {
            input:
                vcfFiles = chr_dv_t.dvVcf,
                gvcfFiles = chr_dv_t.dvgVcf,
                outname = sampleName
        }
    }
    
    # otherwise run dv whole genome
    if(length(chrs) == 0{
        call dv_margin_t.dv_t{
        input:
            threads = threads,
            reference = referenceFile,
            bamAlignment = bamFile,
            bamAlignmentIndex = bamIdxFile,
            sampleName = sampleName
        }
    }
    
    # phase the variants.
    if (phaseVariants){
        call dv_margin_t.margin_t{
            input:
                threads = threads,
                reference = referenceFile,
                bamAlignment = bamFile,
                bamAlignmentIndex = bamIdxFile,
                vcfFile = dv_t.dvVcf,
                sampleName = sampleName
        }
    }

    output {
        File? phasedVcf = margin_t.phasedVcf
        File  dvUnphasedVcf  = dv_t.dvVcf
        File  dvUnphasedgVcf = dv_t.dvVcf
        #File phasedgVcf = margin_t.phasedgVcf
        #File haplotaggedBam = margin_t.haplotaggedBam
        #File haplotaggedBamBai = margin_t.haplotaggedBamIdx
    }
}
