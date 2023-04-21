## Tiny test dataset

Simulated tiny dataset to test locally that the workflow runs.

Prepare the dataset with:

```sh
# simulate reference and reads
python3 sim_tiny_test_data.py

# gzip FASTQ
gzip -f test_ont.fastq

# make unmapped BAM
docker run -it -u `id -u $USER` -v `pwd`:/app quay.io/biocontainers/samtools:1.16.1--h6899075_1 samtools import -0 /app/test_ont.fastq.gz -o /app/test_ont.u.bam -O BAM
```

This will create a reference FASTA (`ref.fa`) and some long reads (`test_ont.fastq.gz`).
