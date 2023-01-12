CARD Nanopore Varaint Calling and Assmebly Workflows
===================================================

This repository contains pieplines for variant calling and de novo assembly of ONT data,
optimized for [single-flowcell ONT sequencing protocol](https://dx.doi.org/10.17504/protocols.io.ewov1n93ygr2/v1).
The wet-lab/informatics protocol is now applied to sequence and characterize thousands of human brain genomes at 
the [Center for Alzheimer's and Related Dementias at NIH](https://card.nih.gov/).

CARD Data availability
---------------------

The cell line and brain sample sequencing data, assemblies and variant calls are currently
being uploaded to [Terra](https://terra.bio/). We will post the links as soon as the data is available.
If you want to access the data sooner, don't hesitate to contact us in the meantime.

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

Cromwell could be configures to run on an HPC or cloud. This
configuration is more involving and requires optimization for a particular environemnt.
Please refer to the [corresponding manual](https://cromwell.readthedocs.io/en/stable/Configuring/) for details

### Quick demo

XXX

Input data requirements
-----------------------

The pipeline was tested using 30-40x ONT sequencing using R9.4 pore with read N50 ~30kb.
Basecalling and mehtylation calls were done using Guppy 6.1.2. The pipeline should
work for similar or newer nanopore data. We are currently planning to release
a special version of this pipeline for R10 ONT data.

The input for end-to-end workflow is a single unmapped bam file with methylation tags
produced by Guppy. Other workflows can take either unmapped bam or fastq file as input.

Other kinds of input include reference genome and corresponding VNTR annotations (provided
in this repository).

Pipeline description
---------------------

The pipeline begins by generating a diploid de novo assembly using a combination of Shasta, 
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

License 
--------

Workflows are distributed under a BSD license. See the LICENSE file for details.

Credits
-------

The pipeline was developed at in collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/)
and the [Cancer Data Science Laboratory](https://ccr.cancer.gov/cancer-data-science-laboratory), National Cancer Institute.

Main code contributors.

* Mikhail Kolmogorov (NCI)
* Mira Mastoras (UCSC)
* Melissa Meredith (UCSC)

Citation
--------

Kolmogorov, Billingsley et al, "Scalable Nanopore sequencing of human genomes provides a 
comprehensive view of haplotype-resolved variation and methylation". bioRxiv 2023
