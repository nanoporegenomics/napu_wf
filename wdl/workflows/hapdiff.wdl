version 1.0

import "../tasks/hapdiff.wdl" as hapdiff_t

workflow hapdiffWf {

    input {
        File referenceFile
        File? vntrAnnotations
		String sample = "Sample"
        File assemblyHap1
        File assemblyHap2
        Int minSvSize = 30
        Int threads = 32
    }

	### hapdiff
	call hapdiff_t.hapdiff_t as hapdiff_t {
		input:
			threads=threads,
			reference=referenceFile,
			vntrAnnotations=vntrAnnotations,
			sample=sample,
			ctgsPat=assemblyHap1,
			ctgsMat=assemblyHap2,
			minSvSize=minSvSize
	}

	output {
		File vcfUnphased = hapdiff_t.hapdiffUnphasedVcf
		File vcfPhased = hapdiff_t.hapdiffPhasedVcf
	}
}
