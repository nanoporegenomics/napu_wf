## Tiny test dataset

Simulated tiny dataset to test locally that the workflow runs.

Prepare the dataset with:

```sh
python3 sim_tiny_test_data.py
gzip -f test_ont.fastq
```

This will create a reference FASTA (`ref.fa`) and some long reads (`test_ont.fastq.gz`).
