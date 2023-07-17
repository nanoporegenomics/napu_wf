Napu
====

### Pipelines version: R10

Napu (Nanopore Analysis Pipeline) is a collection of WDL workflows for variant calling and de novo assembly of ONT data,
optimized for [single-flowcell ONT sequencing protocol](https://dx.doi.org/10.17504/protocols.io.ewov1n93ygr2/v1).
The wet-lab/informatics protocol is now applied to sequence and characterize thousands of human brain genomes at 
the [Center for Alzheimer's and Related Dementias at NIH](https://card.nih.gov/).

Versions for R9/R10 data
------------------------

Please use either `r9` or `r10` branch for your corresponding data type.

Installation and Usage
---------------------

### Using Terra/DNAnexus/AnVIL

If you are using Terra, the workflows are already available at the 
[Dockstore collection](https://dockstore.org/organizations/NIHCARD/collections/NanoporeSequencing).
All you need to do is to import the relevant workflow to your Terra workspace.

### On a single local compute node

There are multiple existing WDL engine implementations. We performed our tests using
[Cromwell](https://cromwell.readthedocs.io/en/stable/), and the following instructions
assume this WDL implementation.

First, follow [these instructions](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) 
to download and install the latest version of Cromwell. Also make sure that [Docker](https://docs.docker.com/get-docker/)
is installed and running.

Then, you will need to prepare an inputs file. For example, for the end-to-end pipeline, 
you can generate blank input file as follows (`-o` option hides optional workflow inputs)

```
java -jar womtool-XX.jar inputs -o false wdl/workflows/cardEndToEndVcf.wdl > inputs.json
```

Then, edit the `inputs.json` file with your input parameters. Afterwards, you will be able to run the corresponding 
workflow (for example, `cardEndToEndVcf.wdl`) as follows:

```
java -jar cromwell-XY.jar run wdl/workflows/cardEndToEndVcf.wdl --inputs inputs.json
```

### On a custom HPC server or cloud environemnt

Cromwell could be configured to run on an HPC or cloud. This
configuration is more involving and requires optimization for a particular environemnt.
Please refer to the [corresponding manual](https://cromwell.readthedocs.io/en/stable/Configuring/) for details

### Running the assembly with Shasta separately

The Shasta assembler can be run in *in-memory* mode (faster/cheaper) or *disk-backed* mode (slower but doesn't require as much memory).
The end-to-end workflow uses the *disk-backed* mode by default because Terra can't (currently) launch instances with more than 768 Gb of memory which might not be enough for some samples.
However, because the *in-memory* mode is much cheaper (~$20 vs ~$70), it might be worth attempting to run Shasta in *in-memory* mode first.
Then, the rest of the workflow could be run on the samples that worked, or the full workflow with the *disk-backed* mode on the few that failed.

Shasta can be run separately with the workflow defined at `wdl/tasks/shasta.wdl` and deposited [on Dockstore](https://dockstore.org/workflows/github.com/jmonlong/card_nanopore_wf/shasta:r10?tab=info).
To use the *in-memory* mode on Terra, the suggested inputs are:
- `inMemory=true`
- `diskSizeGB=1500` (needed to get high-memory instances on Terra)

Then, the FASTA file produced by Shasta can be provided to the end-to-end workflow described above (defined at `wdl/workflows/cardEndToEndVcf.wdl` and deposited [on Dockstore](https://dockstore.org/workflows/github.com/jmonlong/card_nanopore_wf/cardEndToEndVcfMethyl:r10?tab=info)) using the optional input `shastaFasta`.

### Quick demo

```
wget https://zenodo.org/record/7532080/files/card_pipeline_small_test.tar.gz
tar -xvf card_pipeline_small_test.tar.gz
cd card_pipeline_small_test
java -jar PATH_TO/cromwell-XX.jar run PATH_TO/cardEndToEndVcf.wdl --inputs inputs_end2end.js
```

Make sure cromwell is installed and substitute the paths to cromwell and WDL workflow according to your setup.

Input data requirements
-----------------------

Napu was tested using 30-40x ONT sequencing using R9.4 pore with read N50 ~30kb.
Basecalling and mehtylation calls were done using Guppy 6.1.2. Napu should
work for similar or newer nanopore data. We are currently planning to release
a special version of this pipeline for R10 ONT data.

The input for end-to-end workflow is a single unmapped bam file with methylation tags
produced by Guppy. Other workflows can take either unmapped bam or fastq file as input.

Other kinds of input include reference genome and corresponding VNTR annotations (provided
in this repository).


bioRxiv manuscript data availability
------------------------------------

The cell line data (HG002, HG0073 and HG02723) and openly available through this [Terra workspace](https://anvil.terra.bio/#workspaces/anvil-datastorage/ ANVIL_NIA_CARD_Coriell_Cell_Lines_Open). 

Human brain sequencing datasets are under controlled access and require a dbGap application (phs001300.v4). Afterwards, the data will be available through the [restricted Terra workspace](https://anvil.terra.bio/#workspaces/anvil-datastorage/ANVIL_NIA_CARD_LR_WGS_NABEC_GRU).

Pipeline description
---------------------

Napu begins by generating a diploid de novo assembly using a combination of Shasta, 
which produces a haploid assembly and Hapdup, which generates locally phased diploid contigs. 
We then use the generated assemblies to call structural variants (at least 50 bp in size) 
against a given reference genome using a new assembly-to-reference pipeline called hapdiff.

Ideally, small variants could also be recovered from diploid contigs, as has been successfully done 
for HiFi-based assemblies. Our Shasta-Hapdup assemblies had mean substitution error rates of ~8 per 100 kb, which is  
higher than current contig assemblies produced with PacBio HiFi (<1 per 100kb). 
Reference-based small variant calling methods are less error-prone because they can better 
approximate the biases of ONT base calling errors via supervised learning. 
We therefore use  PEPPER-Margin-DeepVariant software to call small variants against a reference.

Given a set of structural variants produced with de novo assemblies, and reference-based small variant calls, 
our pipeline phases them into a harmonized variant call set using Margin. 
In addition, given the phased reference alignment with methylation tags (produced by Guppy), 
we produce haplotype-specific calls of hypo- and hyper-methylated regions of the genome.

The workflows are buit around the following tools:

* [Pepper-Margin-DeepVariant](https://github.com/kishwarshafin/pepper)
* [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
* [hapdiff](https://github.com/KolmogorovLab/hapdiff)
* [margin phase](https://github.com/UCSC-nanopore-cgl/margin)
* [Shasta](https://github.com/chanzuckerberg/shasta)
* [Hapdup](https://github.com/KolmogorovLab/hapdup)
* [modbam2bed](https://github.com/epi2me-labs/modbam2bed)

Credits
-------

Napu was developed at in collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/)
and the [Cancer Data Science Laboratory](https://ccr.cancer.gov/cancer-data-science-laboratory), National Cancer Institute.

Main code contributors.

* Mikhail Kolmogorov (NCI)
* Mira Mastoras (UCSC)
* Melissa Meredith (UCSC)
* Jean Monlong (UCSC)

Citation
--------
Kolmogorov, Billingsley et al, "Scalable Nanopore sequencing of human genomes provides a 
comprehensive view of haplotype-resolved variation and methylation". bioRxiv 2023
[doi.org/10.1101/2023.01.12.523790](https://doi.org/10.1101/2023.01.12.523790)


License 
--------
Workflows are distributed under a BSD license. See the LICENSE file for details.
