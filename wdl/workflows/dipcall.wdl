version 1.0

import "../tasks/dipcall.wdl" as dipcall_t

workflow dipcallWf {

    input {
        File referenceFile
        File assemblyHap1
        File assemblyHap2
        Int threads = 32
    }

	call dipcall_t.dipcall_t as dipcall {
		input:
			threads=threads,
			reference=referenceFile,
			ctgsPat=assemblyHap1,
			ctgsMat=assemblyHap2
	}

	output {
		File vcfDipcall = dipcall.dipcallVcf
	}
}
