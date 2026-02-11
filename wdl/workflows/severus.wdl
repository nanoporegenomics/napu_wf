version 1.0

workflow runSeverus {
	input {
		String sample
		File haplotaggedBAM 
		File? phasedVCF
		Int threads = 16
		String dockerImage = "meredith705/severus:1.6"
	}

	call severus {
		input:
			bam = haplotaggedBAM, 
			vcf = phasedVCF, 
			threads = threads,
			dockerImage = dockerImage
		}
	
	output {
		File severus_tar = severus.severus_tar
		}
	}



task severus {
	input {
		String sample
		File bam
		File? vcf
		Int threads
		String dockerImage
		Boolean includeVNTRS = true
		Boolean includePON = true
		String extraArgs = ""
		Int memSizeGb = 2 * round(size(bam, 'G'))
		Int disksizeGb = 2 * round(size(bam, 'G'))
	}

	String vntrsArg = if includeVNTRS then "--vntr-bed /opt/Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed" else ""
	String ponArg = if includePON then "--PON /opt/Severus/pon/PoN_1000G_hg38.tsv.gz" else ""
	String phasedVcf = if defined(vcf) then "--phasing-vcf" else ""

	command <<<
		set -eux -o pipefail

		# source the conda install
		source /opt/conda/etc/profile.d/conda.sh
		# activate the severus environment
		conda activate severus_env

		./opt/Severus/severus.py --target-bam {bam} --out-dir {sample}_severus -t {threads} \
		{vntrsArg} {ponArg} {phasedVcf}{vcf} {extraArgs}

		tar -czf {sample}_severus.tar.gz {sample}_severus

	>>>

	output {
		File severus_tar = "{sample}_severus.tar.gz"
	}

}
