version 1.0

import "../tasks/minimap2.wdl" as minimap2_t
import "../tasks/shasta.wdl" as shasta_t
import "../tasks/hapdup.wdl" as hapdup_t

workflow structuralVariantsDenovoAssembly {

    input {
        File readsFile
        Int threads
        Int shastaDiskSizeGB = 1024
    }

    ### Shasta assembly ###
    call shasta_t.shasta_t as shasta_t {
        input:
        #threads=threads,
        reads=readsFile,
        diskSizeGb=shastaDiskSizeGB
    }

	### minimap2 alignent ###
    call minimap2_t.minimap2_t as minimap2 {
        input:
            threads=threads,
            reference=shasta_t.shastaFasta,
            reads=readsFile,
            useEqx=false
    }

	### hapdup
	call hapdup_t.hapdup_t as hapdup_t {
		input:
			threads=threads,
			alignedBam=minimap2.bam,
			contigs=shasta_t.shastaFasta
	}

	output {
        File asmDual1 = hapdup_t.hapdupDual1
        File asmDual2 = hapdup_t.hapdupDual2
        File asmPhased1 = hapdup_t.hapdupPhased1
        File asmPhased2 = hapdup_t.hapdupPhased2
        File phaseBed1 = hapdup_t.hapdupPhaseBed1
        File phaseBed2 = hapdup_t.hapdupPhaseBed2 
		File shastaHaploid = shasta_t.shastaFasta
		File shastaLog = shasta_t.shastaLog
	}
}
