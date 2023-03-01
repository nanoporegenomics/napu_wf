version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/dv-margin.wdl" as dv_margin_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "shasta_hapdup_denovo.wdl" as denovo_asm_wf
import "marginPhase.wdl" as margin_phase_wf

workflow cardEndToEndVcf
{
	input {
		File  inputFastq	
		File  referenceFasta
		Int   threads
		File?  referenceVntrAnnotations
        Int   nbReadsPerChunk = 0
		String  sampleName = "sample"
        Array[String] chrs = []
	}

    if(nbReadsPerChunk == 0){
	    call minimap_t.minimap2_t as mm_align {
		    input:
			reads = inputFastq,
			reference = referenceFasta,
			threads = threads
	    }
    }
    if(nbReadsPerChunk > 0){
	    call minimap_t.splitReads {
		    input:
			reads = inputFastq,
            readsPerChunk = nbReadsPerChunk
	    }

        scatter (readChunk in splitReads.readChunks){
	        call minimap_t.minimap2_t as mm_align_chunk {
		        input:
			    reads = readChunk,
			    reference = referenceFasta,
                preemptible = 2,
			    threads = threads
	        }
        }

        call minimap_t.mergeBAM {
		    input:
			bams = mm_align_chunk.bam,
            outname = sampleName
	    }
    }

    File bamFile = select_first([mm_align.bam, mergeBAM.bam])
    File bamFileIndex = select_first([mm_align.bamIndex, mergeBAM.bamIndex])
    
    if(length(chrs) > 0){
        # shard by chromosome
        scatter (chrn in chrs){
	        call dv_margin_t.dv_t as chr_dv_t {
		        input:
		        threads = threads,
		        reference = referenceFasta,
		        bamAlignment = bamFile,
		        bamAlignmentIndex = bamFileIndex,
                region = chrn
	        }
        }

        call dv_margin_t.mergeVCFs {
            input:
            vcfFiles = chr_dv_t.dvVcf,
            outname = sampleName
        }
    }
    if(length(chrs) == 0){
	    call dv_margin_t.dv_t{
		    input:
		    threads = threads,
		    reference = referenceFasta,
		    bamAlignment = bamFile,
		    bamAlignmentIndex = bamFileIndex
	    }
    }
    File dvVCF = select_first([mergeVCFs.vcf, dv_t.dvVcf])
    
    call dv_margin_t.margin_t{
		input:
		threads = threads,
		reference = referenceFasta,
		bamAlignment = bamFile,
		bamAlignmentIndex = bamFileIndex,
        vcfFile = dvVCF,
        sampleName = sampleName
	}

	call sniffles_t.sniffles_t as sniffles {
		input:
			threads = threads,
			bamAlignment = margin_t.haplotaggedBam,
		    bamAlignmentIndex = margin_t.haplotaggedBamIdx,
			vntrAnnotations = referenceVntrAnnotations
	}

	call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
		input:
		readsFile = inputFastq,
        chunkedReadsFiles=select_first([splitReads.readChunks, []]),
		threads = threads,
	}

	call hapdiff_t.hapdiff_t as hapdiff {
		input:
			ctgsPat = asm.asmDual1,
			ctgsMat = asm.asmDual2,
			reference = referenceFasta,
			vntrAnnotations = referenceVntrAnnotations
	}

	call margin_phase_wf.runMarginPhase as margin_phase {
		input:
			smallVariantsFile = margin_t.phasedVcf,
			structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
			refFile = referenceFasta,
			bamFile = margin_t.haplotaggedBam,
			sampleName = sampleName
	}

	output {
		File phasedBam = margin_t.haplotaggedBam
		File smallVariantsVcf = margin_t.phasedVcf
		File snifflesVcf = sniffles.snifflesVcf
		File assemblyHap1 = asm.asmDual1
		File assemblyHap2 = asm.asmDual2
		File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
		File harmonizedVcf = margin_phase.out_margin_phase_svs
	}
}
