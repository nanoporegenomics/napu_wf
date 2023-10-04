version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/modbam2bed.wdl" as modbam2bed_t
import "../tasks/pepper-margin-dv.wdl" as pmdv_haplotag_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "../tasks/dipcall.wdl" as dipcall_t
import "shasta_hapdup_denovo.wdl" as denovo_asm_wf
import "marginPhase.wdl" as margin_phase_wf

workflow cardEndToEndVcfMethyl
{
    input {
        Array[File] inputReads  = []
        File        referenceFasta
        Int         threads
        File?       referenceVntrAnnotations
        File?       shastaFasta
        Array[File] inputMappedBams = []
        Int         nbReadsPerChunk = 0
        String      sampleName = "sample"
        Array[String] chrs = []
        Boolean runModbam2bed = true
        File? inputHaplotaggedBam
        File? inputHaplotaggedBamIdx
        File? inputPhasedVCF
        File? singleInputMappedBamIdx
    }

    parameter_meta {
        inputReads: "Array of Unmapped BAM/s or FASTQ file/s containing ONT R9 reads."
        referenceFasta: "Reference"
        threads: "Threads to pass to minimap2, DV, Margin, & Shasta"
        referenceVntrAnnotations: "Optional vntr annotation input for sniffles"
        shastaFasta: "Optional input Shasta assembly, assembly is skipped in workflow"
        inputMappedBams: "Array of input sorted BAMs aligned to the reference"
        sampleName: "Name of Sample"
        nbReadsPerChunk: "Number of reads to put into a chunk for parallel alignment chunking"
        chrs: "list of mapped chromosomes to run PMDV on, instead of calling variants whole genome"
        inputHaplotaggedBam: "haplotagged BAM from a previous run, skips PMDV"
        inputHaplotaggedBamIdx: "haplotagged BAM.bai file froma previous run"
        inputPhasedVCF: "small variant PMDV VCF from previous run"
        singleInputMappedBamIdx: "mapped BAM.bai from a previous run, skips indexing again"
        runModbam2bed: "modbam2bed boolean flag in preparation of switching to modkit"
    }

    ### Either align input, merge multiple mapped input, or reorganize single input
    ## If only one mapped bam is provided just index it (and split into chromosomes if chrs is provided input)
    if (length(inputMappedBams) == 1){
        File inputBam = select_first(inputMappedBams)
        if (!defined(singleInputMappedBamIdx)){
            call minimap_t.indexBAM as indexSingleInputBam{
                input:
                    bam = inputBam,
                    chrs = chrs
            }
        }
    }

    # If more than one mapped bams are provided as input, merge them into one
    if (length(inputMappedBams) > 1){
        call minimap_t.mergeBAM as mergeInputBams{
            input:
            bams = inputMappedBams,
            outname = sampleName,
            chrs = chrs
            }
    }

    ## If input ubam/fastq files are provided align the input reads
    if (length(inputReads) > 0){
        scatter (inputReadsFile in inputReads){

            ##### Aligning the reads to the reference genome
            ## Reads can be split into chunks to make shorter jobs that can be run on parallel preemptible instances
            # However, if no chunks are asked for align the whole file
            if(nbReadsPerChunk == 0){
                call minimap_t.minimap2_t as mm_align {
                    input:
                        reads = inputReadsFile,
                        reference = referenceFasta,
                        threads = threads
                }
            }
            # Split the reads and align them to the reference
            if(nbReadsPerChunk > 0){
                call minimap_t.splitReads {
                    input:
                        reads = inputReadsFile,
                        readsPerChunk = nbReadsPerChunk
                }
                scatter (readChunk in splitReads.readChunks){
                    call minimap_t.minimap2_t as mm_align_chunk {
                        input:
                            reads = readChunk,
                            reference = referenceFasta,
                            preemptible = 2,
                            threads = threads
                    }
                }
            }
        }
        if(nbReadsPerChunk == 0){
            call minimap_t.mergeBAM as mergeAlignedBAMs {
                input:
                    bams = select_all(mm_align.bam),
                    outname = sampleName,
                    chrs=chrs
            }
        }
        if(nbReadsPerChunk > 0){
            call minimap_t.mergeBAM as mergeScatteredBAMs {
                input:
                    bams = flatten(select_all(mm_align_chunk.bam)),
                    outname = sampleName,
                    chrs=chrs
            }
        }
        Array[File] chunkedReads = flatten(select_all(splitReads.readChunks))

    }


    ## Select aligned reads to the reference genome from the multiple possible inputs
    File bamFile = select_first([inputBam, mergeInputBams.bam, mergeAlignedBAMs.bam, mergeScatteredBAMs.bam])
    File bamFileIndex = select_first([singleInputMappedBamIdx, indexSingleInputBam.bamIndex, mergeInputBams.bamIndex, mergeAlignedBAMs.bamIndex, mergeScatteredBAMs.bamIndex])

    if (!defined(inputHaplotaggedBam)){
        ##### Reference-based variant calling with DeepVariant
        ## if the reads/BAMs were chunked by chromosomes, use directly those chunks
        if(length(chrs) > 0 && nbReadsPerChunk > 0){
            Array[File] bamChrs = select_first([mergeScatteredBAMs.bamPerChrs, mergeAlignedBAMs.bamPerChrs, indexSingleInputBam.bamPerChrs, mergeInputBams.bamPerChrs, mergeScatteredBAMs.bamPerChrs])
            Array[File] bamChrsIndex = select_first([mergeScatteredBAMs.bamPerChrsIndex, mergeAlignedBAMs.bamPerChrsIndex, indexSingleInputBam.bamPerChrsIndex, mergeInputBams.bamPerChrsIndex, mergeScatteredBAMs.bamPerChrsIndex])
            scatter (bamChr in zip(bamChrs, bamChrsIndex)){
                call pmdv_haplotag_t.pepper_margin_dv_t as pmdvHap_chrs{
                    input:
                        threads = threads,
                        reference = referenceFasta,
                        bamAlignment = bamChr.left,
                        bamAlignmentIndex = bamChr.right,
                        sampleName = sampleName,
                        oneChr = true,
                        preemptible = 2
                }
            }
            # merge the per chr haplotagged bams
            call minimap_t.mergeBAM as merge_PEPPER_DV_BAMs {
                input:
                    bams = select_all(pmdvHap_chrs.haplotaggedBam),
                    outname = sampleName
                }

            # merge the vcfs from each chromosome
            call pmdv_haplotag_t.mergeVCFs{
                input:
                    vcfFiles = pmdvHap_chrs.pepperVcf,
                    outname = sampleName
            }
        }

        ## otherwise, one job for the whole-genome on the entire BAM
        if(length(chrs) == 0 || nbReadsPerChunk == 0){
            call pmdv_haplotag_t.pepper_margin_dv_t as pmdvHap{
                input:
                    threads = threads,
                    reference = referenceFasta,
                    bamAlignment = bamFile,
                    bamAlignmentIndex = bamFileIndex,
                    sampleName = sampleName,
                    preemptible = 2
            }
        }
    }


    # isolate the haplotagged BAM from PEPPER_Margin_DeepVariant
	File haplotaggedBam = select_first([inputHaplotaggedBam, merge_PEPPER_DV_BAMs.bam, pmdvHap.haplotaggedBam])
	File haplotaggedBamIdx = select_first([inputHaplotaggedBamIdx, merge_PEPPER_DV_BAMs.bamIndex, pmdvHap.haplotaggedBamIdx])

    ##### estimate methylation at CpG sites
    ## only if the input reads file was a BAM
    if (length(inputMappedBams) > 0){
        File mappedBAM1 = select_first(inputMappedBams)
    }
    if (length(inputReads) > 0){
        File mappedRead1 = select_first(inputReads)
    }
    File inReadFile = select_first([mappedRead1, mappedBAM1])
    if(basename(inReadFile, ".bam") != basename(inReadFile) && runModbam2bed){
        call modbam2bed_t.modbam2bed as modbam2bed {
            input:
                haplotaggedBam = haplotaggedBam,
                haplotaggedBamBai = haplotaggedBamIdx,
                ref = referenceFasta,
                sample_name = sampleName
        }
    }

	call sniffles_t.sniffles_t as sniffles {
		input:
			#threads = threads,
			bamAlignment = haplotaggedBam,
			bamAlignmentIndex = haplotaggedBamIdx,
			vntrAnnotations = referenceVntrAnnotations
	}

    ##### De novo phased assembly

    ## if any fastq reads are suppled as input use those for shasta
    if(basename(inReadFile, ".bam") == basename(inReadFile)){
        ## If one fastq is provided as input read/s store as a File
        if (length(inputReads) == 1){
            File readFile = select_first(inputReads)
        }

        ## or merge multiple unaligned read fastqs into a single File
        if (length(inputReads) > 1){
            call minimap_t.mergeFASTQ as mergeInReadsFQs{
                input:
                    reads = inputReads,
                    outname = sampleName,
            }
        }
        File singleReadsFastq = select_first([mergeInReadsFQs.fq, readFile])
    }

    # if any non-BAM reads are supplied as input use those for shasta
    File shastaInputReads = select_first([singleReadsFastq, bamFile])

    ## Run assembly
	call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
		input:
            readsFile = shastaInputReads,
            chunkedReadsFiles=select_first([chunkedReads, []]),
            shastaFasta = shastaFasta,
            threads = threads
	}

	call hapdiff_t.hapdiff_t as hapdiff {
		input:
			ctgsPat = asm.asmDual1,
			ctgsMat = asm.asmDual2,
			reference = referenceFasta,
			vntrAnnotations = referenceVntrAnnotations,
            sample = sampleName
	}

	call dipcall_t.dipcall_t as dipcall {
		input:
			ctgsPat = asm.asmDual1,
			ctgsMat = asm.asmDual2,
			reference = referenceFasta
	}

    File phasedVCF = select_first([inputPhasedVCF, mergeVCFs.vcf, pmdvHap.pepperVcf])
	call margin_phase_wf.runMarginPhase as margin_phase {
		input:
			smallVariantsFile = phasedVCF,
			structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
			refFile = referenceFasta,
			bamFile = haplotaggedBam,
			sampleName = sampleName
	}

	output {
		File phasedMethylBam = haplotaggedBam
        File phasedMethylBamBai = haplotaggedBamIdx
		File smallVariantsVcf = phasedVCF
		File? methylationBed1 = modbam2bed.hap1bedOut
		File? methylationBed2 = modbam2bed.hap2bedOut
		File snifflesVcf = sniffles.snifflesVcf
        File shastaHaploid = asm.shastaHaploid
        File? shastaGFA = asm.shastaGFA
        File? shastaLog = asm.shastaLog
		File assemblyHap1 = asm.asmPhased1
		File assemblyHap2 = asm.asmPhased2
        File asmHap1PhaseBed = asm.phaseBed1
        File asmHap2PhaseBed = asm.phaseBed2
        File assemblyDual1 = asm.asmDual1
        File assemblyDual2 = asm.asmDual2
		File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
		File harmonizedVcf = margin_phase.out_margin_phase_svs
		File asmDipcallVcf = dipcall.dipcallVcf
	}
}
