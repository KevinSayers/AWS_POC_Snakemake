### Overview

RNAseq is one the most common NGS datatypes. It typically comes in two flavors:

- total RNAseq
- polyA RNAseq

The core of most RNAseq data analysis workflows are mostly identical involving the following steps:

- trim fastq reads
- align trimmed reads to genome using split aware aligner
- count reads-per-gene based on existing gene model
- collect QC metrics at various stages 

The above steps assume that one already has a pre-built genomic index to align the trimmed fastq reads against. If you do not have such index, the genomic sequence in FASTA format and gene annotations in GTF format are required to generate a index prior to alignment.

The follow graphic depicts the basic workflow that we will attempt to create on AWS:
<!-- ![diagram](https://i.imgur.com/Rh95JU9.png =250x) -->
<img src="https://i.imgur.com/Rh95JU9.png" alt="drawing" width="600"/>

The inputs will be read from a s3 bucket and the outputs will be written to the same s3 bucket. Individual dockers will be created for various steps of the pipeline as described below:

| Pipeline Step | Tool in Docker |
| ------------- | -------------- |
| Trim          | [CutAdapt<sup>1</sup>](https://cutadapt.readthedocs.io/)       |
| Align         | [STAR<sup>2</sup>](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)           |
| BuildIndex    | [STAR<sup>2</sup>](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)           |
| Count         | [RSEM<sup>3</sup>](https://github.com/deweylab/RSEM#readme-for-rsem)           |
| Report        | [MultiQC<sup>4</sup>](https://multiqc.info/)        |

### Inputs

The inputs for the pipeline can be broadly be classified into:

- inputs required to create the index for mapping, i.e., genomic FASTA file and its gene annotations file (GTF)
- raw FASTQ reads per sample

#### Fasta

This is the genome in FASTA format. We will be using hg38 version of the human genome for our purposes and it can be downloaded from [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz).

#### GTF

GTF or Gene Transfer Format is a file which includes all the gene annotations or gene models. This includes the information about the location of the genes and their splicing events. We will use the GENCODEs [release 38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz) for this.

#### FASTQ

Paired end (PE) FASTQ files per sample will be provided as input. _which dataset to use as test input is TBD_

### Outputs

The pipeline is expected to produce 2 primary outputs:

#### CountsMatrix

This is a tab-delimited file with _nsample + 1_ number of columns (first column is the gene identifier and it is followed by one column per sample of counts data) and _ngenes + 1_ number of rows (first row is the header containing sample names and every subsequent row gives counts per gene tab-delimited for each sample).

#### Report

Basic stats are collected to generated a HTML report using [MultiQC](https://multiqc.info/).


### Pipelining frameworks

There are 3 or 4 pipelining frameworks which are popular among bioinformatic community, specifically for NGS data analysis. [CCBR](https://github.com/CCBR) has successfully used [Snakemake<sup>5</sup>](https://snakemake.readthedocs.io/en/stable/) for the past few years to execute reproducible NGS data analysis on the [Biowulf](https://hpc.nih/gov) HPC cluster. We will be leveraging this extensive experience to build a Snakemake-based pipeline for the above outlined workflow to be run on the AWS cloud specifically using [AWS Genomics CLI](https://aws.amazon.com/genomics-cli/).

In addition to Snakemake, a [Nextflow](https://www.nextflow.io/)-based pipeline will also be built to mimic the same tasks achieved by the previously mentioned Snakemake-based workflow. The major reasons for this repetition of efforts are:

- AWS Genomic CLI adaptation of Nextflow seems more mature than Snakemake.
- Nextflow workflows seemed to be more widely accepted eg. [SBG](https://www.sevenbridges.com/press/releases/seven-bridges-expands-platform-to-support-execution-of-nextflow-and-wdl-workflows/), [DNAnexus](https://documentation.dnanexus.com/developer/workflows/importing-workflows), etc.
- Assess differences in building workflows in Snakemake vs Nextflow. This insight will be valuable in shaping future pipeline-development directions we take.

#### Snakemake

We will be using [CCBR's RNAseq repository](https://github.com/CCBR/RNA-seek) as reference to build our MVP workflow in Snakemake.

#### Nextflow

We will be using [Nfcore's RNAseq pipeline](https://nf-co.re/rnaseq) as a reference to build our MVP workflow in Nextflow


### References

<sup>1</sup> [http://dx.doi.org/10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)

<sup>2</sup> [http://dx.doi.org/10.1093/bioinformatics/bts635](http://dx.doi.org/10.1093/bioinformatics/bts635)

<sup>3</sup> [https://doi.org/10.1186/1471-2105-12-32](https://doi.org/10.1186/1471-2105-12-32)

<sup>4</sup> [http://dx.doi.org/10.1093/bioinformatics/btw354](http://dx.doi.org/10.1093/bioinformatics/btw354)

<sup>5</sup> [https://doi.org/10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)