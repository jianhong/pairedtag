# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

# Sample and unit setup

The sample and unit setup is specified via tab-separated tabular files (`.tsv`).
Missing values can be specified by empty columns or by writing `NA`.

## sample sheet

The default sample sheet is `config/samples.tsv` (as configured in `config/config.yaml`).
Each sample refers to an actual physical sample, and replicates (both biological and technical) may be specified as separate samples.
For each sample, you will always have to specify a `sample_name`.

## unit sheet

The default unit sheet is `config/units.tsv` (as configured in `config/config.yaml`).
For each sample, add one or more sequencing units (for example if you have several runs or lanes per sample).

### `.fastq` file source

For each unit, you will have to define a source for your `.fastq` files.
This can be done via the columns `fq1`, `fq2`, with two `.fastq` files for paired-end reads.

### adapter trimming

TODO
