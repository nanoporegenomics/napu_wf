version 1.0

import "../tasks/minimap2.wdl" as minimap_t


workflow indexTheBAM
{
	input {
		File inputMappedBam 
		Boolean sortingInputBAM = false
		Array[String] chrs = []
	}

	call minimap_t.indexBAM as indexBAM_t {
		input: 
                bam = inputMappedBam,
                sortInputBAM = sortingInputBAM,
                chrs = chrs
	}

	output {
		File out_bam_index = indexBAM_t.bamIndex
		File? out_sorted_bam = indexBAM_t.sortedBam
		#Array[File]? chr_bams = bamChrs
        #Array[File]? chr_bams_idx = bamChrsIndex
	}
}