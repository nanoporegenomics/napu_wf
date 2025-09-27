#!/bin/bash

input_vcf=$1
sample_id=$2
window=$3
min_cluster_size=$4
threshold_SD=$5

filteredbed="${sample_id}.merged_small_svs.filt${window}bp_${threshold_SD}_sds.bed"

python3 /opt/filter_high_depth_variants.py $input_vcf $sample_id $filteredbed $window $min_cluster_size $threshold_SD


mergedFilteredBed="${sample_id}.merged_small_svs.filt${window}bp_${threshold_SD}_sds.100kbmerged.bed"
bedtools merge -i $filteredbed -d 100000 > $mergedFilteredBed
echo "written filtered regions to ${mergedFilteredBed}"

filtVcf="${sample_id}.merged_small_svs.${threshold_SD}_sd_depthFilt.vcf.gz"
bcftools view -T ^$mergedFilteredBed $input_vcf -Oz -o $filtVcf
echo "filtered VCF output: ${filtVcf}"

tabix $filtVcf
