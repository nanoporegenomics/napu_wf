version 1.0



workflow phaseblocks {
    input {
        File phasedVCF 
        String sample 
        #Int diskSizeGB = 2 * round(size(phasedVCF, "G")) 
        #Int memSizeGb = 128

    }


    call phaseblocks_t {
        input:
            phasedVCF=phasedVCF,
            sample=sample
    }

    output {
        File bed = phaseblocks_t.bed
    }
}

task phaseblocks_t {
    
    input {
        File phasedVCF
        String sample
        Int memSizeGB = 20
        Int threads = 26
        Int diskSizeGB = 2 * round(size(phasedVCF, "GB")) 
    }

    command <<<

        set -eux -o pipefail

        # run the phase block calculator
        python3 /opt/vcf_phaseBlock_calculator.py {phasedVCF} --sample {sample}

    >>>

    output {
        File bed = "~{sample}.phaseset.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/card_harmonize_vcf:0.2"
    }

}

