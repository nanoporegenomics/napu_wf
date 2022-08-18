version 1.0

import "../tasks/minimap2.wdl" as minimap2_t

workflow minimapWorkflow {

    input {
        File referenceFile
        File readsFile
        Int threads
    }

    ### minimap2 alignent ###
    call minimap2_t.minimap2_t as minimap2 {
        input:
            threads=threads,
            reference=referenceFile,
            reads=readsFile,
    }

	output {
        File minimapBam = minimap2.bam
	}
}
