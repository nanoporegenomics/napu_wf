#!/bin/bash

input_vcf=$1
sample_id=$2
window=$3
min_cluster_size=$4
threshold_SD=$5

# compose the name of the filtered bed output
filteredbed="${sample_id}.merged_small_svs.filt${window}bp_${threshold_SD}_sds.bed"
# use the script to identify dense regions of variants
python3 /opt/filter_high_depth_variants.py $input_vcf $sample_id $filteredbed $window $min_cluster_size $threshold_SD

# if there any high density clusters identified remove them from the vcf
if [ -f "$filteredbed" ] && [ -s "$filteredbed" ]; then
	# compose the name of the bed after merging nearby dense clusters
	mergedFilteredBed="${sample_id}.merged_small_svs.filt${window}bp_${threshold_SD}_sds.100kbmerged.bed"
	bedtools merge -i $filteredbed -d 100000 > $mergedFilteredBed
	echo "written filtered regions to ${mergedFilteredBed}"

	# use bcftools to remove those sense clusters from the filtered vcf for harmonized phasing
	filtVcf="${sample_id}.merged_small_svs.${threshold_SD}_sd_depthFilt.vcf.gz"
	bcftools view -T ^$mergedFilteredBed $input_vcf -Oz -o $filtVcf
	echo "filtered VCF output: ${filtVcf}"

	tabix $filtVcf
fi
