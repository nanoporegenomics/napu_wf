version 1.0

import "../tasks/sniffles.wdl" as sniffles_t

workflow snifflesWf {

    input {
        File bamAlignment
        File vntrAnnotations
        Int threads
    }

	### Sniffles
    call sniffles_t.sniffles_t as sniffles_t {
        input:
            threads=threads,
			bamAlignment=bamAlignment,
			vntrAnnotations=vntrAnnotations
    }

	output {
        File structuralVariantsSniffles = sniffles_t.snifflesVcf
	}
}
