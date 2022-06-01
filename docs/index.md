## DESCRIPTION


This is a living document for the HPC component of the MVP on AWS. This is a code repository and a living document describing ideas/"how-tos"/slides etc. related to the Minimal Viable Product (MVP) [^1] to be built for the Genetics Branch (GB) on Cloud One or AWS.

> <ins>Disclaimer</ins>: There are two parallel efforts underway for the _"cloudification"_ of GB, namely, a. Moving the database management and its web-interface to AWS also referred to as **"the Database component"** and b. Orchestrating NGS analysis workflows on AWS which is also referred to as **"the HPC compotent"**. This repository solely focuses on _the HPC component_.

As MVP, we propose creating a basic RNAseq workflow on AWS using AWS CLI. We would like to use **both** the following pipelining frameworks:

- [Snakemake](https://snakemake.readthedocs.io/)
- [Nextflow](https://www.nextflow.io/)

The code [repository](https://github.com/CCRGeneticsBranch/AWS_MVP_HPC) will hold the following:

- Snakemake pipeline
- Nextflow pipeline
- Recipes for Docker containers used by the pipeline and
- This documentation


[^1]: Please send your comments/questions/suggestions to [Vishal Koparde](mailto:vishal.koparde@nih.gov).