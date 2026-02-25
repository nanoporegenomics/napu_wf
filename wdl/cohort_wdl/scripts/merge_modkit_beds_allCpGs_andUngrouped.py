#!/usr/bin/env python3

import pandas as pd
import polars as pl
import os
import sys
import gzip
import gcsfs
from datetime import datetime
import argparse


"""
Script to merge individual cpg sites across CARD cohorts

inputs: a TSV with samples as rows, first column is sample ID
        haplotype 1 modkit bed files 

example: python3 merge_modkit_beds_allCpGs.py -i NABEC_cohort_methyl_012025.tsv -o out_HarmPhase_methylationBeds

        
This script reads in data using gs links so you must authenticate using: 
gcloud auth application-default login
or 
gcloud auth application-default login --no-launch-browser

Author: Melissa Meredith
2/2026
"""

def log_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def read_in_gslinks(tsvfile):
    """
    Function reads in a tsv of gcp links to phased modkit files, 
        columns must be named 
        
        methylationBed1
        methylationBed2
        methylationBedUngrouped

    """

    log_time(f"Reading in methylation gs links from {tsvfile}")
    delimiter="\t"
    if tsvfile[-4:]!=".tsv":
        print('is input file not a .tsv?')
        sys.exit(-1)

    try:
        input_df = pd.read_csv(tsvfile, sep=delimiter)
    except Exception as e:
        log_time(f"Error reading TSV: {e}")
        sys.exit(-1)

    hap1_links = input_df['methylationBed1'].to_list()
    hap2_links = input_df['methylationBed2'].to_list()
    ungrouped_links = input_df['methylationBedUngrouped'].to_list()

    log_time('done')
    return hap1_links, hap2_links, ungrouped_links

def merge_beds(gslinks, outputdir, haplotype):
    """
    Function reads in each modkit bed and stores the valid coverage, 
    number of reads with mods at each position and the modified fraction
    column. The individual dataframes are merged together based on chromosomal 
    position. Every 10 samples the dataframe is written to disk. 
    """

    log_time('Innitialize gcs sytem')
    # initialize GCS FileSystem
    fs = gcsfs.GCSFileSystem()

    # empty dataframe to fill with cpg data
    combined_df = None

    for i, file in enumerate(gslinks, 1):

        print('file', file)
        if not isinstance(file, str):
            log_time(f"Invalid file entry at line {i}. Skipping...")
            continue

        # get sample name from the file path ex: NABEC_KEN-1066_FTX_GRCh38_2.bed.gz
        sample_name = file.split('/')[-1].replace('.bed.gz', '')
        log_time(f'Processing file: {sample_name}')

        try:
            # open GCS file and read using gzip 
            # with fs.open(file, 'rb') as f:
                # with gzip.open(f, 'rt') as gz_file:
            # polars can just read the gs link
            df = pl.read_csv(
                file,
                separator="\t",
                has_header=False,
                columns=[0, 1, 2, 3, 9, 10, 11],
                new_columns=[
                    "#chrom",
                    "start",
                    "end",
                    "modtype",
                    f"{sample_name}_validCov",
                    f"{sample_name}_modFraction",
                    f"{sample_name}_modReads",
                ],
            )

        except Exception as e:
            log_time(f"Error reading {sample_name}: {e}")
            log_time(f"Make sure you authenticated your google account with 'gcloud auth application-default login --no-launch-browser'")
            continue


        # outer merge on genomic coordinates
        if combined_df is None:
            combined_df = df
        else:
            combined_df = combined_df.join(
                df,
                on=["#chrom", "start", "end", "modtype"],
                how="outer",
                coalesce=True,
            )

        log_time(f'Completed {sample_name}')

        # # write out the merged df every 100 samples to prevent loosing all merged data
        if i%100==0:
            log_time(f'writing out combined file with {i} samples')
            combined_df.to_csv(f'{outputdir}/combined_methylation_{haplotype}.tsv', sep="\t", index=False)

    log_time('Filling in zeros')
    # fill missing values with a 0 for zero coverage/measurements of that position
    combined_df = combined_df.fill_null(0)

    log_time('making output tsv')
    # save the combined data to CSV
    combined_df.write_csv(f'{outputdir}/combined_methylation_{haplotype}.tsv', separator="\t")


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Merge methylBeds from Modkit across cohorts.")
    parser.add_argument(
        "-i","--in_tsv_file",
        type=str,
        required=True,
        help="Path to the input gs link TSV file with header. The first column should be the sample IDs, subsequent columns are haplotype1 and hap2 file links."
    )

    parser.add_argument(
        "-o","--output_directory",
        type=str,
        required=True,
        help="name of output directory."
    )

    parser.add_argument(
        '--merge_groups',
        nargs='*',        # '+' = one or more, '*' = zero or more
        type=str,
        help='A space-separated list of groups to merge: hap1, hap2, ungrouped. if none are provided all three will be merged.'
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    # Read in the gs links from imput tsv
    hap1_gs, hap2_gs, ungrouped_gs = read_in_gslinks(args.in_tsv_file)
    
    # Create output directory
    output_dir = args.output_directory
    os.makedirs(output_dir, exist_ok=True)
    log_time(f"Output directory ensured at: {output_dir}")

    print(args.merge_groups)
    # Merge the data within each haplotype
    # methylation matches Variant genoytpes - not methylation status 

    # if a particular merge is input only run that
    if len(args.merge_groups) > 0:
        if 'hap1' in args.merge_groups:
            merge_beds(hap1_gs, output_dir, 'hap1')
        if 'hap2' in args.merge_groups:
            merge_beds(hap2_gs, output_dir, 'hap2')
        if 'ungrouped' in args.merge_groups:
            merge_beds(ungrouped_gs, output_dir, 'ungrouped')
    else:
        merge_beds(hap1_gs, output_dir, 'hap1')
        merge_beds(hap2_gs, output_dir, 'hap2')
        merge_beds(ungrouped_gs, output_dir, 'ungrouped')



