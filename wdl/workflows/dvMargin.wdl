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
    
    call dv_margin_t.dv_t{
        input:
            threads = threads,
            reference = referenceFile,
            bamAlignment = bamFile,
            bamAlignmentIndex = bamIdxFile,
            sampleName = sampleName
    }

    call dv_margin_t.margin_t{
        input:
            threads = threads,
            reference = referenceFile,
            bamAlignment = bamFile,
            bamAlignmentIndex = bamIdxFile,
            vcfFile = dv_t.dvVcf,
            gvcfFile = dv_t.dvgVcf,
            sampleName = sampleName
    }

    output {
        File phasedVcf = margin_t.phasedVcf
        File phasedgVcf = margin_t.phasedgVcf
        File haplotaggedBam = margin_t.haplotaggedBam
        File haplotaggedBamBai = margin_t.haplotaggedBamIdx
    }
}
