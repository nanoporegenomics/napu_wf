CARD Nanopore Varaint Calling and Assmebly Workflows
===================================================

1-2 paragraphs description

WDL workflows for processing nanopore sequencing of brain samples, generated at NIH CARD.

CARD Data availability
---------------------

Installation
------------

Dockstore collection: https://dockstore.org/organizations/NIHCARD/collections/NanoporeSequencing

Link to Cromwell installation


Quick example
--------------

Pipeline description
---------------------

Usage
-----

The workflows are buit around the followinf tools:

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
