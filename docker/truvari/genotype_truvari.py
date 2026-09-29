import sys
import os
# import gzip

def parse_bed(bed_file):
    """
    Parses a BED file containing assembly coverages.
    Returns a dictionary where keys are genomic positions (0-based) and values are coverages.
    """
    coverage_dict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                parts = line.strip().split(":")
                sample = parts[1]
                print(sample)
            else:
                parts = line.strip().split("\t")
                # print(parts)
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                # contig = int(parts[3])
                # for pos in range(start, end):
                # coverage_dict[(chrom, start)] = end
                if chrom in coverage_dict:
                    coverage_dict[chrom].append((start,end))
                else:
                    coverage_dict[chrom] = [(start,end)]
    return sample, coverage_dict


def fill_vcf(vcf_file, sample_name, coverage_dict, output_file):
    """
    Fills in the genotype of a specified sample column in a multi-sample VCF with "0/0".
    if the genotype entry is './.'.
    Writes the modified VCF to the output file.
    Returns False (and writes nothing) if sample_name is not in the VCF.
    """
    # check if the sample is in the VCF first
    sample_index = None  # Initialize sample index
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                # check for sample in header
                if sample_name in header[header.index('FORMAT') + 1:]:
                    sample_index = header.index(sample_name, header.index('FORMAT') + 1)
                break

    # if the sample is not found, exit this function and write nothing
    if sample_index is None:
        print(f'{sample_name} is not in vcf. Exiting without writing output.')
        return False

    tmp_output_file = "tmp_"+output_file
    print('writing to tmp:', tmp_output_file)
    
    gtFillCounter=0
    with open(vcf_file, 'r') as f, open(tmp_output_file, 'w') as out:
        for line in f:
            # write out the header
            if line.startswith('#'):
                out.write(line)
                continue

            # for variant lines compare ./. genotype calls with 
            # sample coverage bed, and replace with 0/0 if covered
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1]) - 1  # VCF is 1-based, BED is 0-based
            # print('vcf', chrom)

            # look for missing genotypes and check for coverage
            genotype_info = parts[sample_index].split(':')
            if genotype_info[0] == './.':  # genotype_info[0][-1] != '1': #
                if (chrom) in coverage_dict.keys():
                    for start,end in coverage_dict[chrom]:
                        # only change genotype if sv is covered by a single block
                        if pos >= start and pos <= end:
                            if sample_index is not None:
                                genotype = "0/0"
                                parts[sample_index] = ":".join([genotype] + parts[sample_index].split(':')[1:])
                                gtFillCounter+=1
                                break

            out.write('\t'.join(parts) + '\n')


    # Rename the temporary output VCF file to the final output VCF file
    print('moving tmp:', tmp_output_file, 'to:', output_file, f'with {gtFillCounter} filled genotypes')
    os.rename(tmp_output_file, output_file)
    return True



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <bed_file> <multi_sample_vcf> <output_vcf>")
        sys.exit(1)

    bed_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_file = sys.argv[3]

    sampleName, coverage_dict = parse_bed(bed_file)
    fill_vcf(vcf_file, sampleName, coverage_dict, output_file)
    print("Genotypes filled successfully.")

