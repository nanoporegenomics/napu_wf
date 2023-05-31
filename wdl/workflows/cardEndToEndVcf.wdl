version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/dv-margin.wdl" as dv_margin_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "../tasks/dipcall.wdl" as dipcall_t
import "../tasks/modbam2bed.wdl" as modbam2bed_t
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
    }

    parameter_meta {
        inputReads: "Array of Unmapped BAM/s or FASTQ file/s containing ONT R10 reads."
        referenceFasta: "Reference"
        threads: "Threads to pass to minimap2, DV, Margin, & Shasta"
        referenceVntrAnnotations: "Optional vntr annotation input"
        shastaFasta: "Optional input Shasta assembly, assembly is skipped in workflow"
        inputMappedBams: "Array of input sorted BAMs aligned to the reference"
        sampleName: "Name of Sample"
        nbReadsPerChunk: "Number of reads to put into a chunk for using preemptible instances"
    }

    ### Either align input, merge multiple mapped input, or reorganize single input
    ## If one or more mapped bams are provided as input, merge them into one
    if (length(inputMappedBams) == 1){
        File inputBam = select_first(inputMappedBams)
        call minimap_t.indexBAM as indexSingleInputBam{
            input: 
                bam = inputBam,
                chrs = chrs
        }
    }
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
            ## Reads can be split into chunks to make shorter jobs that can be run on cheaper preemptible instances 
            if(nbReadsPerChunk == 0){
                call minimap_t.minimap2_t as mm_align {
                    input:
                        reads = inputReadsFile, 
                        reference = referenceFasta,
                        threads = threads
                }
            }
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


    ## Aligned reads to the reference genome 
    File bamFile = select_first([indexSingleInputBam.outBam, mergeInputBams.bam, mergeAlignedBAMs.bam, mergeScatteredBAMs.bam])
    File bamFileIndex = select_first([indexSingleInputBam.bamIndex, mergeInputBams.bamIndex, mergeAlignedBAMs.bamIndex, mergeScatteredBAMs.bamIndex])
    

    ##### Reference-based variant calling with DeepVariant
    ## if the reads/BAMs were chunked by chromosomes, use directly those chunks
    if(length(chrs) > 0 && nbReadsPerChunk > 0){
        Array[File] bamChrs = select_first([mergeScatteredBAMs.bamPerChrs, mergeAlignedBAMs.bamPerChrs, indexSingleInputBam.bamPerChrs, mergeInputBams.bamPerChrs, mergeScatteredBAMs.bamPerChrs])
        Array[File] bamChrsIndex = select_first([mergeScatteredBAMs.bamPerChrsIndex, mergeAlignedBAMs.bamPerChrsIndex, indexSingleInputBam.bamPerChrsIndex, mergeInputBams.bamPerChrsIndex, mergeScatteredBAMs.bamPerChrsIndex])
        scatter (bamChr in zip(bamChrs, bamChrsIndex)){
            call dv_margin_t.dv_t as chr_dv_t {
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
        call dv_margin_t.mergeVCFs {
            input:
                vcfFiles = chr_dv_t.dvVcf,
                outname = sampleName
        }
    }
    ## otherwise, one job for the whole-genome on the entire BAM
    if(length(chrs) == 0 || nbReadsPerChunk == 0){
        call dv_margin_t.dv_t{
            input:
                threads = threads,
                reference = referenceFasta,
                bamAlignment = bamFile,
                bamAlignmentIndex = bamFileIndex,
                sampleName = sampleName
        }
    }
    ## Variant calls from DeepVariant
    File dvVCF = select_first([mergeVCFs.vcf, dv_t.dvVcf])

    ##### Haplotag the reads
    call dv_margin_t.margin_t{
        input:
            threads = threads,
            reference = referenceFasta,
            bamAlignment = bamFile,
            bamAlignmentIndex = bamFileIndex,
            vcfFile = dvVCF,
            sampleName = sampleName
    }


    ##### estimate methylation at CpG sites
    ## only if the input reads file was a BAM
    if (length(inputMappedBams) > 0){
        File mappedBAM1 = select_first(inputMappedBams)
    }
    if (length(inputReads) > 0){
        File mappedRead1 = select_first(inputReads)
    }
    File inReadFile = select_first([mappedRead1, mappedBAM1])
    if(basename(inReadFile, ".bam") != basename(inReadFile)){
        call modbam2bed_t.modbam2bed as modbam2bed {
            input:
                haplotaggedBam = margin_t.haplotaggedBam,
                haplotaggedBamBai = margin_t.haplotaggedBamIdx,
                ref = referenceFasta,
                sample_name = sampleName
        }
    }

    ##### Reference-based structural variant calling
    call sniffles_t.sniffles_t as sniffles {
        input:
            bamAlignment = margin_t.haplotaggedBam,
            bamAlignmentIndex = margin_t.haplotaggedBamIdx,
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

    # if any non-BAM reads are suppled as input use those for shasta
    File shastaInputReads = select_first([singleReadsFastq, bamFile])

    ## Run assembly
    call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
        input:
            readsFile = shastaInputReads, 
            chunkedReadsFiles=select_first([chunkedReads, []]),
            shastaFasta = shastaFasta,
            threads = threads
    }

    ##### Assembly-based structural variant calling
    call hapdiff_t.hapdiff_t as hapdiff {
        input:
            ctgsPat = asm.asmDual1,
            ctgsMat = asm.asmDual2,
            reference = referenceFasta,
            vntrAnnotations = referenceVntrAnnotations
    }

    call dipcall_t.dipcall_t as dipcall {
        input:
            ctgsPat = asm.asmDual1,
            ctgsMat = asm.asmDual2,
            reference = referenceFasta
    }

    ##### Phase short variants and structural variants
    call margin_phase_wf.runMarginPhase as margin_phase {
        input:
            smallVariantsFile = margin_t.phasedVcf,
            structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
            refFile = referenceFasta,
            bamFile = margin_t.haplotaggedBam,
            sampleName = sampleName
    }

    output {
        File phasedBam = margin_t.haplotaggedBam
        File smallVariantsVcf = margin_t.phasedVcf
        File snifflesVcf = sniffles.snifflesVcf
        File assemblyHap1 = asm.asmDual1
        File assemblyHap2 = asm.asmDual2
        File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
        File harmonizedVcf = margin_phase.out_margin_phase_svs
        File? methylationBed1 = modbam2bed.hap1bedOut
        File? methylationBed2 = modbam2bed.hap2bedOut
        File asmDipcallVcf = dipcall.dipcallVcf
    }
}
