version 1.0

import "../tasks/minimap2.wdl" as minimap_t
import "../tasks/dv-margin.wdl" as dv_margin_t
import "../tasks/sniffles.wdl" as sniffles_t
import "../tasks/hapdiff.wdl" as hapdiff_t
import "../tasks/modbam2bed.wdl" as modbam2bed_t
import "shasta_hapdup_denovo.wdl" as denovo_asm_wf
import "marginPhase.wdl" as margin_phase_wf

workflow cardEndToEndVcfMethyl
{
    input {
        File        inputReads
        File        referenceFasta
        Int         threads
        File?       referenceVntrAnnotations
        File?       shastaFasta
        Array[File] inputBams = []
        Int         nbReadsPerChunk = 0
        String      sampleName = "sample"
        Array[String] chrs = []
    }

    parameter_meta {
        inputReads: "Fasta file containing reads to be assembled & variant called."
        referenceFasta: "Reference"
        threads: "Threads to pass to minimap2, DV, Margin, & Shasta"
        referenceVntrAnnotations: "Optional vntr annotation input"
        shastaFasta: "Optional input Shasta assembly, assembly is skipped in workflow"
        inputBams: "Array of input sorted BAMs aligned to the reference"
        sampleName: "Name of Sample"
        nbReadsPerChunk: "Number of reads to put into a chunk for using preemptible instances"

    }

    ## If no input bam files are provided align the input reads
    if (length(inputBams) == 0){
        ##### Aligning the reads to the reference genome
        ## Reads can be split into chunks to make shorter jobs that can be run on cheaper preemptible instances 
        if(nbReadsPerChunk == 0){
            call minimap_t.minimap2_t as mm_align {
                input:
                    reads = inputReads,
                    reference = referenceFasta,
                    threads = threads
            }
        }
        if(nbReadsPerChunk > 0){
            call minimap_t.splitReads {
                input:
                    reads = inputReads,
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
            ## if "chrs" is not empty, BAMs for each specified chromosomes will also be output
            call minimap_t.mergeBAM {
                input:
                    bams = mm_align_chunk.bam,
                    outname = sampleName,
                    chrs=chrs
            }
        }
    }

    ## If one or more bams are provided as input, merge them into one
    if (length(inputBams) > 0){
        call minimap_t.mergeBAM as mergeInputBam{
                input:
                    bams = inputBams,
                    outname = sampleName,
                    chrs=chrs
            }
    }



    ## Aligned reads to the reference genome 
    File bamFile = select_first([mergeInputBam.bam, mm_align.bam, mergeBAM.bam])
    File bamFileIndex = select_first([mergeInputBam.bamIndex, mm_align.bamIndex, mergeBAM.bamIndex])

    ##### Reference-based variant calling with DeepVariant
    ## if the reads/BAMs were chunked by chromosomes, use directly those chunks
    if(length(chrs) > 0 && nbReadsPerChunk > 0){
        scatter (bamChr in zip(select_first([mergeBAM.bamPerChrs]), select_first([mergeBAM.bamPerChrsIndex]))){
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
    if(basename(inputReads, ".bam") != basename(inputReads)){
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
    call denovo_asm_wf.structuralVariantsDenovoAssembly as asm {
        input:
            readsFile = inputReads,
            chunkedReadsFiles=select_first([splitReads.readChunks, []]),
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
    }
}
