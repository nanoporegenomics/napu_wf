version 1.0

workflow runSeverus {
	input {
		String sample
		File haplotaggedBAM 
		File? phasedVCF
		Int threads = 16
		Int preemptible_count = 0
		String dockerImage = "meredith705/severus:1.6"
	}

	call severus {
		input:
			bam = haplotaggedBAM, 
			sample = sample,
			vcf = phasedVCF, 
			threads = threads,
			dockerImage = dockerImage
		}
	
	output {
		File severus_tar = severus.severus_tar
		File severus_somatic_vcf = severus.severus_somatic_vcf
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
		Int preemptible_count
		Int memSizeGb = 2 * round(size(bam, 'G'))
		Int diskSizeGb = 2 * round(size(bam, 'G'))
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
		File severus_somatic_vcf = "{sample}_severus/somatic_SVs/{sample}_somatic.vcf"
	}

	runtime {
        preemptible: preemptible_count
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }

}
