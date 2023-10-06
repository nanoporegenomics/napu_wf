version 1.0

import "../tasks/modkit.wdl" as modkit_t

workflow modkit_wf
{
	input {
		File  haplotaggedBam	
		File  haplotaggedBamIdx
		File  referenceFasta
		String  sampleName = "sample"
		String  referenceName = "ref"
	}

	call modkit_t.modkit as modkit {
		input:
			haplotaggedBam = haplotaggedBam,
			haplotaggedBamBai = haplotaggedBamIdx,
			ref = referenceFasta,
			sample_name = sampleName,
			ref_name = referenceName
	}

	output {
		File methylationBed1 = modkit.hap1bedOut
		File methylationBed2 = modkit.hap2bedOut
		File? methylationBedUngrouped = modkit.ungroupedBedOut
	}
}
