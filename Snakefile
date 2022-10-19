
from os.path import join
from os import listdir
import os
import re


from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

SAMPLES = ["ggal_gut"]

#workpath="./"


rule all:
	input:
		expand(join("{sample}","fastqc","{sample}_1_fastqc.zip"),sample=SAMPLES),	


rule fastqc:
	input:
		R1 = S3.remote("myncitestbucket/inputs/{sample}_1.fastq",sample=SAMPLES),
		R2 = S3.remote("myncitestbucket/inputs/{sample}_2.fastq",sample=SAMPLES),
	params:
		out = join("{sample}","fastqc")
	output:
		join("{sample}","fastqc","{sample}_1_fastqc.zip"),
		join("{sample}","fastqc","{sample}_2_fastqc.zip")

	conda: 'env.yaml'
		
	shell: """

	fastqc {input} -o  {params.out}

	"""
