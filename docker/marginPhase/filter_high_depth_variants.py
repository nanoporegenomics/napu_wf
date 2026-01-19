#!/usr/bin/env python3
import sys
import statistics
import pysam

def compute_depth_stats(vcf_file, sample):
    """Compute genome-wide mean and stdev of per-sample depth."""
    dps = []
    with pysam.VariantFile(vcf_file) as vcf:
        if sample not in vcf.header.samples:
            raise ValueError(f"Sample {sample} not found in VCF")
        for chrom in vcf.header.contigs:
            for rec in vcf.fetch(chrom):
                if "DP" in rec.samples[sample] and rec.samples[sample]["DP"] is not None:
                    dps.append(rec.samples[sample]["DP"])
    if not dps:
        raise ValueError(f"No DP values found for sample {sample}")
    mean_dp = statistics.mean(dps)
    stdev_dp = statistics.pstdev(dps)
    return mean_dp, stdev_dp


def find_dense_variants(vcf_file, sample, mean_dp, stdev_dp,
                        window_size=100000, min_cluster_size=10, threshold=2.0):
    """Find clusters of densely packed variants with DP > mean + threshold*stdev."""
    clusters = []
    cutoff = mean_dp + threshold * stdev_dp

    with pysam.VariantFile(vcf_file) as vcf:
        for chrom in vcf.header.contigs:
            current_cluster = []
            for rec in vcf.fetch(chrom):
                dp = rec.samples[sample].get("DP")
                if dp is None:
                    continue
                pos = rec.pos
                if dp > cutoff:
                    if not current_cluster:
                        current_cluster.append((chrom, pos, dp))
                    else:
                        prev_chrom, prev_pos, _ = current_cluster[-1]
                        if chrom == prev_chrom and pos - prev_pos <= window_size:
                            current_cluster.append((chrom, pos, dp))
                        else:
                            if len(current_cluster) >= min_cluster_size:
                                clusters.append(current_cluster)
                            current_cluster = [(chrom, pos, dp)]
                else:
                    if len(current_cluster) >= min_cluster_size:
                        clusters.append(current_cluster)
                    current_cluster = []

            if len(current_cluster) >= min_cluster_size:
                clusters.append(current_cluster)

    return clusters, cutoff


def write_bed(clusters, bed_file):
    """Write clusters to a BED file with size and mean DP."""
    with open(bed_file, "w") as out:
        for cluster in clusters:
            chrom = cluster[0][0]
            start = cluster[0][1] - 1  # BED start is 0-based
            end = cluster[-1][1]       # BED end is exclusive
            size = len(cluster)
            mean_dp = statistics.mean(dp for _, _, dp in cluster)
            out.write(f"{chrom}\t{start}\t{end}\tcluster_size={size};meanDP={mean_dp:.2f}\n")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <vcf_file.gz> <sample_name> <output.bed> [window_size] [min_cluster_size] [threshold_SD]")
        sys.exit(1)

    vcf_file = sys.argv[1]
    sample = sys.argv[2]
    bed_file = sys.argv[3]
    window_size = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    min_cluster_size = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    threshold = float(sys.argv[6]) if len(sys.argv) > 2 else 2.0

    try:
        mean_dp, stdev_dp = compute_depth_stats(vcf_file, sample)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    clusters, cutoff = find_dense_variants(vcf_file, sample, mean_dp, stdev_dp,
                                           window_size, min_cluster_size, threshold)

    print(f"Sample: {sample}")
    print(f"Genome-wide mean DP: {mean_dp:.2f}")
    print(f"Genome-wide stdev DP: {stdev_dp:.2f}")
    print(f"Using cutoff: DP > {cutoff:.2f} ({threshold} SDs above mean)")
    print(f"Found {len(clusters)} dense clusters (written to {bed_file})")

    if len(clusters) > 0:
        write_bed(clusters, bed_file)
