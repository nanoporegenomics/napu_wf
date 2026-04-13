#!/usr/bin/env python3

import pandas as pd
import os
import sys
import gzip
from datetime import datetime
import argparse


def log_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def merge_partial_cohort_tsvs(parital1, parital2):
	"""
	04/2026
	Perform an outter merge on '#chrom', 'start', 'end', 'modtype'
	Added modtype now that we call 5mC and 5hmC, the modtype column matters
	"""
	

	try:
		log_time(f'load in {parital1}: {parital2[-3:]}')
		if parital1[-3:] == ".gz":
			with gzip.open(parital1,'rb') as file1:
				p1 = pd.read_csv(file1, sep='\t')
		else:
			p1 = pd.read_csv(parital1, sep='\t')

		log_time(f'load in {parital2}')
		if parital2[-3:] == ".gz":
			with gzip.open(parital1,'rb') as file2:
				p2 = pd.read_csv(file2, sep='\t')
		else:
			p2 = pd.read_csv(parital2, sep='\t')
	except Exception as e:
		log_time(f"Error reading files {e}" )

	log_time(f'merge files')
	combined_df = pd.merge(p1, p2, on=['#chrom', 'start', 'end', 'modtype'], how='outer')

	log_time(f'replace NaNs with 0')
	combined_df.fillna(0, inplace=True)


	outfile = parital1.split("/")[-1][:-3]+"_combined_"+parital2.split("/")[-1][:-4]+".tsv"
	combined_df.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Merge partially combined methylbeds from Modkit across cohorts. Merged on ['chrom', 'start', 'end']")
    parser.add_argument(
        "-a","--in_a_tsv",
        type=str,
        required=True,
        help="Path to the first input TSV file with header ( can be .gz)."
    )

    parser.add_argument(
        "-b","--in_b_tsv",
        type=str,
        required=True,
        help="Path to the second input TSV file with header ( can be .gz)."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()


    merge_partial_cohort_tsvs(args.in_a_tsv, args.in_b_tsv)
