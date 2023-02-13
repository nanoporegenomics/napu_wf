# CARD nanopore workflows

WDL workflows for processing nanopore sequencing of brain samples, generated at NIH CARD.

Dockstore collection: https://dockstore.org/organizations/NIHCARD/collections/NanoporeSequencing

1. [Pepper-Margin-DeepVariant](https://github.com/kishwarshafin/pepper)
2. [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
3. [hapdiff](https://github.com/KolmogorovLab/hapdiff)
4. [margin phase](https://github.com/UCSC-nanopore-cgl/margin)
5. [Shasta](https://github.com/chanzuckerberg/shasta)
6. [Hapdup](https://github.com/KolmogorovLab/hapdup)


## Test locally

```
## DV+margin workflow
miniwdl run --as-me  --copy-input-files wdl/workflows/dvMargin.wdl -i test/test.input.dvmargin.json

## End-to-end workflow (no methylation)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcfFastq.wdl -i test/test.input.endtoendfastq.json
```
