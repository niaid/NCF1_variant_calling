# Snakemake workflow: dna-seq-gatk-variant-calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.1-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)](https://cdn.rawgit.com/snakemake-workflows/dna-seq-gatk-variant-calling/master/.test/report.html)

This Snakemake pipeline implements the [GATK best-practices workflow](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) for calling small genomic variants.

This workflow is adapted from this Snakemake pipeline: https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling

Many updates were made to the pipeline for calling variants in duplicated regions, like the NCF1 gene.

## Authors

* Eric Karlins


## Usage

You'll first need to create a fasta file with duplicated regions that match your region of interest masked with Ns.

#### Step 1: Create masked reference fasta file

1. For NCF1, using the hg38 reference fasta, I created a bed file with two lines to mask the regions of NCF1b and NCF1c. This bed looked like this:
```
chr7    73220639    73235945
chr7    75156639    75172044
```

2. Next I used `bedtools maskfasta` to create the masked reference file. The command was:

```
bedtools maskfasta -fi Homo_sapiens_assembly38_plus.fasta -bed NCF1_region_to_mask.bed -fo Homo_sapiens_assembly38_plus_NCF1_mask.fasta
```

3. Create the bwa index files for your new fasta reference:
   
```
bwa index Homo_sapiens_assembly38_plus_NCF1_mask.fasta
```

#### Step 2: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.
