## PIPELINE DETAILS

RNAseq is one the most common NGS datatypes. It typically comes in two flavors:

- total RNAseq
- polyA RNAseq

The core of most RNAseq data analysis workflows are mostly identical involving the following steps:

- trim fastq reads
- align trimmed reads to genome using split aware aligner
- count reads-per-gene based on existing gene model
- collect QC metrics at various stages 

The follow graphic depicts the basic workflow that we will attempt to create on AWS:
<!-- ![diagram](https://i.imgur.com/Rh95JU9.png =250x) -->
<img src="https://i.imgur.com/Rh95JU9.png" alt="drawing" width="600"/>

### INPUTS

#### FASTA

This is the genome in FASTA format. We will be using hg38 version of the human genome for our purposes and it can be downloaded from [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz).

#### GTF

GTF or Gene Transfer Format is a file which includes all the gene annotations or gene models. This includes the information about the location of the genes and their splicing events. We will use the GENCODEs [release 38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz) for this.

### OUTPUTS

The pipeline is expected to produce 2 primary outputs:

#### CountsMatrix

This is a tab-delimited file with _nsample + 1_ number of columns (first column is the gene identifier and it is followed by one column per sample of counts data) and _ngenes + 1_ number of rows (first row is the header containing sample names and every subsequent row gives counts per gene tab-delimited for each sample).