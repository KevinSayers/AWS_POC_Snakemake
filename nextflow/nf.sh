#!/bin/bash

set -e
module load singularity nextflow
#nextflow run -profile biowulf main.nf -resume
nextflow run -profile biowulf 1_main.nf -resume
