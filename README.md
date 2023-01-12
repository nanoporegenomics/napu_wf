CARD Nanopore Varaint Calling and Assmebly Workflows
===================================================

WDL workflows for processing nanopore sequencing of brain samples, generated at NIH CARD.

CARD Data availability
---------------------

The cell line and brain sample sequencing data, assemblies and variant calls are currently
being uploaded to Terra. We will post the links as soon as the data is available.
If you want to access the data sooner, don't hesitate to contact us in the meantime.

Installation and Usage
---------------------

# Using Terra

If you are using Terra, the workflows are already available at the 
[Dockstore collection](https://dockstore.org/organizations/NIHCARD/collections/NanoporeSequencing).
All you need to do is to import the relevant workflow to your Terra workspace.

# On a single local compute node

There are multiple existing WDL engine implementations. We performed our tests using
[Cromwell](https://cromwell.readthedocs.io/en/stable/), and the following instructions
assume this WDL implementation.

First, follow [these instructions](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) 
to download and install the latest version of Cromwell. Make sure that Docker is installed and running.

Then, you will need to prepare an inputs file. For example, for the endToEnd pipeline, 
you can generate blank input file as follows.

```
java womtool.jar XXX
```

Then, you will be able to run the corresponding pipeline (for example, endToEnd) as follows:

```
java -jar cromwell-XY.jar run wdl/workflows/cardEndToEndVcf.wdl --inputs inputs.json
```

# On a custom HPC server or cloud environemnt

# Small dataset example

Prepare a small dataset to illustrate how to run everythin locally.

Pipeline description
---------------------

The workflows are buit around the following tools:

1. [Pepper-Margin-DeepVariant](https://github.com/kishwarshafin/pepper)
2. [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
3. [hapdiff](https://github.com/KolmogorovLab/hapdiff)
4. [margin phase](https://github.com/UCSC-nanopore-cgl/margin)
5. [Shasta](https://github.com/chanzuckerberg/shasta)
6. [Hapdup](https://github.com/KolmogorovLab/hapdup)

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
