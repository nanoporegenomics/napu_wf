version 1.0

import "../tasks/pepper-margin-dv.wdl" as pepper_t

workflow pepperWf {

    input {
        File referenceFile
        File bamAlignment
		Int threads
    }

    call pepper_t.pepper_margin_dv_t as pepper_t {
        input:
            threads=threads,
            reference=referenceFile,
			bamAlignment=bamAlignment
    }

	output {
		File pepperVcf = pepper_t.pepperVcf
		File haplotaggedBam = pepper_t.haplotaggedBam
		File haplotaggedBamBai = pepper_t.haplotaggedBamIdx
	}
}
