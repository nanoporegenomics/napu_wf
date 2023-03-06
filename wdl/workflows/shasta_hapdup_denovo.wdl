version 1.0

import "../tasks/minimap2.wdl" as minimap2_t
import "../tasks/shasta.wdl" as shasta_t
import "../tasks/hapdup.wdl" as hapdup_t

workflow structuralVariantsDenovoAssembly {

    input {
        File readsFile
        File? shastaFasta
        Array[File] chunkedReadsFiles = []
        Int threads
        Int shastaDiskSizeGB = 1024
    }

    
    ### Shasta assembly ###
    ## skip if a shastaFasta is provided
    if(!defined(shastaFasta)){
        if ((basename(readsFile, ".fasta") == basename(readsFile)) && (basename(readsFile, ".fa") == basename(readsFile))){
            call shasta_t.convertToFasta {
                input:
                reads=readsFile
            }
        }
        File readsFasta = select_first([convertToFasta.fasta, readsFile])
        call shasta_t.shasta_t as shasta_t {
            input:
            reads=readsFasta,
            diskSizeGb=shastaDiskSizeGB
        }
    }
    File ambFasta = select_first([shasta_t.shastaFasta, shastaFasta])
    
	### minimap2 alignment ###
    if(length(chunkedReadsFiles) == 0){
	    call minimap2_t.minimap2_t as minimap2 {
		    input:
			reads = readsFile,
            reference=ambFasta,
			useEqx=false,
            threads = threads
	    }
    }
    if(length(chunkedReadsFiles) > 0){
        scatter (readChunk in chunkedReadsFiles){
	        call minimap2_t.minimap2_t as minimap2_chunk {
		        input:
			    reads = readChunk,
                reference=ambFasta,
			    useEqx=false,
                preemptible=2,
			    threads = threads
	        }
        }

        call minimap2_t.mergeBAM as mergeBAMhapdup {
		    input:
			bams = minimap2_chunk.bam
	    }
    }

    File bamFile = select_first([minimap2.bam, mergeBAMhapdup.bam])

	### hapdup
	call hapdup_t.hapdup_t as hapdup_t {
		input:
			threads=threads,
			alignedBam=bamFile,
			contigs=ambFasta
	}

	output {
        File asmDual1 = hapdup_t.hapdupDual1
        File asmDual2 = hapdup_t.hapdupDual2
        File asmPhased1 = hapdup_t.hapdupPhased1
        File asmPhased2 = hapdup_t.hapdupPhased2
        File phaseBed1 = hapdup_t.hapdupPhaseBed1
        File phaseBed2 = hapdup_t.hapdupPhaseBed2 
		File shastaHaploid = ambFasta
		File? shastaLog = shasta_t.shastaLog
	}
}
