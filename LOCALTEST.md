## Test locally

Uses the tiny test dataset in [test](test).

Test locally, for example using [miniwdl](https://github.com/chanzuckerberg/miniwdl):

```
## End-to-end workflow (no methylation)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcfFastq.wdl -i test/test.input.endtoendfastq.json

## End-to-end workflow (no methylation) cheaper by chunking reads (and preempting)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcfFastq.wdl -i test/test.input.endtoendfastq.chunks.json

## DV+margin workflow
miniwdl run --as-me  --copy-input-files wdl/workflows/dvMargin.wdl -i test/test.input.dvmargin.json
```

### Next

- Simulate unmapped BAMs with methylation information
