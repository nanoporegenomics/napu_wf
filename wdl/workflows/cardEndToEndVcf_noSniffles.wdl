version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/dv-margin.wdl" as dv_margin_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "../tasks/dipcall.wdl" as dipcall_t
import "../tasks/modkit.wdl" as modkit_t
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
        Boolean     shastaInMem = false
        File?       hapdupFasta1
        File?       hapdupFasta2
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
    File bamFile = select_first([indexSingleInputBam.sortedBam, inputBam, mergeInputBams.bam, mergeAlignedBAMs.bam, mergeScatteredBAMs.bam])
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
                gvcfFiles = chr_dv_t.dvgVcf,
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
    File dvgVCF = select_first([mergeVCFs.gvcf, dv_t.dvgVcf])

    ##### Haplotag the reads  ?
    #call dv_margin_t.margin_t{
    #    input:
    #        threads = threads,
    #        reference = referenceFasta,
    #        bamAlignment = bamFile,
    #        bamAlignmentIndex = bamFileIndex,
    #        vcfFile = dvVCF,
    #        gvcfFile = dvgVCF,
    #        sampleName = sampleName
    #}

    

    ##### De novo phased assembly
    # if hapdup assembly already provided skip assembly
    if(!defined(hapdupFasta1)){

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
                shastaInMem = shastaInMem,
                threads = threads
        }

    }

    # Isolate the haplotype resolved assemblies
    File asmDual1 = select_first([hapdupFasta1, asm.asmDual1])
    File asmDual2 = select_first([hapdupFasta2, asm.asmDual2])

    ##### Assembly-based structural variant calling
    call hapdiff_t.hapdiff_t as hapdiff {
        input:
            ctgsPat = asmDual1,
            ctgsMat = asmDual2,
            reference = referenceFasta,
            vntrAnnotations = referenceVntrAnnotations,
			sample = sampleName
    }

    call dipcall_t.dipcall_t as dipcall {
        input:
            ctgsPat = asmDual1,
            ctgsMat = asmDual2,
            reference = referenceFasta
    }

    ##### Phase short variants and structural variants
    call margin_phase_wf.runMarginPhase as margin_phase {
        input:
            smallVariantsgVCFFile = dvgVCF,
            structuralVariantsFile = hapdiff.hapdiffUnphasedVcf,
            refFile = referenceFasta,
            bamFile = bamFile, #margin_t.haplotaggedBam,
            sampleName = sampleName
    }


    ##### Reference-based structural variant calling; using the harmonized bam
    #call sniffles_t.sniffles_t as sniffles {
    #    input:
    #        #bamAlignment = margin_t.haplotaggedBam,
    #        #bamAlignmentIndex = margin_t.haplotaggedBamIdx,
    #        bamAlignment = margin_phase.out_margin_phase_bam,
    #        bamAlignmentIndex = margin_phase.out_margin_phase_bam_bai,
    #        reference = referenceFasta,
    #        vntrAnnotations = referenceVntrAnnotations,
    #        sample = sampleName
    #}


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
        call modkit_t.modkit as modkit {
            input:
                haplotaggedBam = margin_phase.out_margin_phase_bam,
                haplotaggedBamBai = margin_phase.out_margin_phase_bam_bai,
                ref = referenceFasta,
                sample_name = sampleName
        }
    }

    output {
        File harmonizedPhasedBam = margin_phase.out_margin_phase_bam
        File harmonizedPhasedBamBai = margin_phase.out_margin_phase_bam_bai
        File harmonizedVcf = margin_phase.out_margin_phase_svs
        File harmonizedVcfIdx = margin_phase.out_phasedVcfIdx
        File harmonizedVcfPhaseset = margin_phase.out_phasedVCFPhaseSetBED
        File margin_out_monitor = margin_phase.margin_out_monitor
        File? harmonizedVcfDenseFilterBed = margin_phase.out_exclusionBed
        File smallVariantsVcf = dvVCF
        File smallVariantsgVcf = dvgVCF
        #File snifflesVcf = sniffles.snifflesVcf
        #File snifflesSnf = sniffles.snifflesSnf
        File? shastaHaploid = asm.shastaHaploid
        #File? shastaLog = asm.shastaLog
        #File? shastaGFA = asm.shastaGfa
        #File? shastaHtml = asm.shastaHtml
        File? assemblyHap1 = asm.asmPhased1
        File? assemblyHap2 = asm.asmPhased2
        File? asmHap1PhaseBed = asm.phaseBed1
        File? asmHap2PhaseBed = asm.phaseBed2
        File? assemblyDual1 = asm.asmDual1
        File? assemblyDual2 = asm.asmDual2
        File structuralVariantsVcf = hapdiff.hapdiffUnphasedVcf
        File alignmentBedHap1 = hapdiff.alignmentBedHap1
        File alignmentBedHap2 = hapdiff.alignmentBedHap2
        File alignmentConfidantBed = hapdiff.confidentBed
        File? methylationBed1 = modkit.hap1bedOut
        File? methylationBed2 = modkit.hap2bedOut
        File? methylationBedUngrouped = modkit.ungroupedBedOut
        File? methylationBedUnPhased = modkit.wholeGenomeOut
        File asmDipcallVcf = dipcall.dipcallVcf
        #Array[File]? chr_bams = bamChrs
        #Array[File]? chr_bams_idx = bamChrsIndex
    }
}
