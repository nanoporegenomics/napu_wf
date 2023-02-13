version 1.0

import "../tasks/methyltagMinimap2.wdl" as minimap_methyl_t
import "../tasks/modbam2bed.wdl" as modbam2bed_t
import "../tasks/dv-margin.wdl" as dv_margin_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "shasta_hapdup_denovo.wdl" as denovo_asm_wf
import "marginPhase.wdl" as margin_phase_wf

workflow cardEndToEndVcfMethyl
{
	input {
		File  inputUnphasedMethylBam	
		File  referenceFasta
		Int   threads
		File  referenceVntrAnnotations = ""
		String  sampleName = "sample"
		String  referenceName = "ref"
	}

	call minimap_methyl_t.fastqAlignAndSortBam as mm_align {
		input:
			unaligned_methyl_bam = inputUnphasedMethylBam,
			ref_file = referenceFasta,
			ref_name = referenceName,
			sample = sampleName,
			in_cores = threads
	}

	call dv_margin_t.dv_t{
		input:
		threads = threads,
		reference = referenceFasta,
		bamAlignment = mm_align.out_bam,
	}

    call dv_margin_t.margin_t{
		input:
		threads = threads,
		reference = referenceFasta,
		bamAlignment = mm_align.out_bam,
        vcfFile = dv_t.dvVcf,
        sampleName = sampleName
	}

	call modbam2bed_t.modbam2bed as modbam2bed {
		input:
			haplotaggedBam = margin_t.haplotaggedBam,
			haplotaggedBamBai = margin_t.haplotaggedBamIdx,
			ref = referenceFasta,
			sample_name = sampleName,
			ref_name = referenceName
	}

	call sniffles_t.sniffles_t as sniffles {
		input:
			threads = threads,
			bamAlignment = margin_t.haplotaggedBam,
			vntrAnnotations = referenceVntrAnnotations
	}

	call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
		input:
			readsFile = inputUnphasedMethylBam,
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
			smallVariantsFile = margin_t.phasedVcf,
			structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
			refFile = referenceFasta,
			bamFile = margin_t.haplotaggedBam,
			sampleName = sampleName
	}

	output {
		File phasedMethylBam = margin_t.haplotaggedBam
		File smallVariantsVcf = margin_t.phasedVcf
		File methylationBed1 = modbam2bed.hap1bedOut
		File methylationBed2 = modbam2bed.hap2bedOut
		File snifflesVcf = sniffles.snifflesVcf
		File assemblyHap1 = asm.asmDual1
		File assemblyHap2 = asm.asmDual2
		File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
		File harmonizedVcf = margin_phase.out_margin_phase_svs
	}
}
