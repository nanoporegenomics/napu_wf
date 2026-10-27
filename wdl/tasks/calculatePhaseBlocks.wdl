version 1.0



workflow phaseblocks {
    input {
        File phasedVCF 
        String sample 

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
        Int memSizeGB = 5
        Int threads = 1
        Int diskSizeGB = 12 
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
        docker: "meredith705/card_harmonize_vcf:0.3"
    }

}

