version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/pepper-margin-dv.wdl" as pmdv_haplotag_t
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
		File  referenceVntrAnnotations = ""
		String  sampleName = "sample"
	}

	call minimap_t.minimap2_t as mm_align {
		input:
			reads = inputFastq,
			reference = referenceFasta,
			threads = threads
	}

	call pmdv_haplotag_t.pepper_margin_dv_t as pmdvHap{
		input:
			threads = threads,
			reference = referenceFasta,
			bamAlignment = mm_align.bam,
	}
	
	call sniffles_t.sniffles_t as sniffles {
		input:
			threads = threads,
			bamAlignment = pmdvHap.haplotaggedBam,
			vntrAnnotations = referenceVntrAnnotations
	}

	call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
		input:
			readsFile = inputFastq,
			threads = threads
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
			smallVariantsFile = pmdvHap.pepperVcf,
			structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
			refFile = referenceFasta,
			bamFile = pmdvHap.haplotaggedBam,
			sampleName = sampleName
	}

	output {
		File phasedBam = pmdvHap.haplotaggedBam
		File smallVariantsVcf = pmdvHap.pepperVcf
		File snifflesVcf = sniffles.snifflesVcf
		File assemblyHap1 = asm.asmDual1
		File assemblyHap2 = asm.asmDual2
		File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
		File harmonizedVcf = margin_phase.out_margin_phase_svs
	}
}
