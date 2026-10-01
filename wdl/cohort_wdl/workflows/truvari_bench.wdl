version 1.0

workflow run_truvari_collapse{
    input {
        File base_vcf
        File base_vcfIdxs 
        File comp_vcf
        File comp_vcfIdxs
        String out_prefix
        String sample

        # Matching parameters
        Float refdist = 1000
        Float pctsize = 0.75
        Float pctseq = 0.75
        Float pctovl = 0.0
        Int typeignore = 0 
        String passonly = "--passonly"

        # Collapse parameters
        String keep = "first"
        String? extraArgs = ""
        String dockerImage = "meredith705/truvari"

    }


    call truvari_bench {
        input:
                    base_vcf = base_vcf,
                    base_vcfIdxs = base_vcfIdxs,
                    comp_vcf = comp_vcf,
                    comp_vcfIdxs = comp_vcfIdxs,
                    refdist = refdist, 
                    pctsize = pctsize, 
                    pctseq = pctseq,
                    pctovl = pctovl,
                    typeignore = typeignore,
                    passonly = passonly,
                    keep = keep,
                    sample = sample,
                    extraArgs = extraArgs,
                    out_prefix = out_prefix,
                    dockerImage = dockerImage
    }

    output{
        File truvari_tarball = truvari_bench.truvari_tarball
    }
}

task truvari_bench {
    input {
        File base_vcf
        File base_vcfIdxs 
        File comp_vcf
        File comp_vcfIdxs 

        # Matching parameters
        Float refdist 
        Float pctsize 
        Float pctseq 
        Float pctovl 
        Int typeignore  
        String passonly 
        String sample
        String out_prefix

        # Collapse parameters
        String keep 
        String? extraArgs  


        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 3 * round(size(base_vcf, 'G')) + 300
        String dockerImage 

    }


    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        truvari bench -b ~{base_vcf} -c ~{comp_vcf} -o ~{out_prefix}_bench_out --pctsize ~{pctsize} --pctseq ~{pctseq} -r ~{refdist} --bSample ~{sample} --cSample ~{sample}

        # compress the whole output directory into one file
        tar -czf ~{out_prefix}_bench_out.tar.gz ~{out_prefix}_bench_out

    >>>

        output {
            File truvari_tarball = "~{out_prefix}_bench_out.tar.gz"

    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

