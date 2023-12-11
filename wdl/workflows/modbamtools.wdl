version 1.0

import "../tasks/modkit.wdl" as modkit_t


workflow modbamtools_wf
{
	input {
		File  haplotaggedBam
		File  haplotaggedBamIdx
		File  referenceFasta
		File  regional_bed
		String  sampleName = "sample"
		String  regionName = "bed_region"
		String  referenceName = "ref"
	}

	call modkit_t.modbamtools as modbamtools {
		input:
			haplotaggedBam = haplotaggedBam,
			haplotaggedBamBai = haplotaggedBamIdx,
			ref = referenceFasta,
			sample_name = sampleName,
			ref_name = referenceName,
			regional_bed = regional_bed,
			region_type = regionName
	}

	output {
		File regionalBed = modbamtools.regionally_aggregated_bed
	}
}
