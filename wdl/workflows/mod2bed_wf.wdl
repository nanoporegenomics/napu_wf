version 1.0

import "../tasks/modbam2bed.wdl" as modbam2bed_t

workflow mod2bed_wf
{
	input {
		File  haplotaggedBam	
		File  haplotaggedBamIdx
		File  referenceFasta
		String  sampleName = "sample"
		String  referenceName = "ref"
	}

	call modbam2bed_t.modbam2bed as modbam2bed {
		input:
			haplotaggedBam = haplotaggedBam,
			haplotaggedBamBai = haplotaggedBamIdx,
			ref = referenceFasta,
			sample_name = sampleName,
			ref_name = referenceName
	}

	output {
		File methylationBed1 = modbam2bed.hap1bedOut
		File methylationBed2 = modbam2bed.hap2bedOut
	}
}
