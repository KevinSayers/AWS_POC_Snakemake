
#def getstar_inputs(wildcards):
#    fastq =[]
#    s = wildcards.name
#    fastq.append(join(fastqdir,s + "_R1.fastq.gz"))
#    fastq.append(join(fastqdir,s + "_R2.fastq.gz"))
#    return fastq
#    print(wildcards)
rule star:
	input: 
		R1= join(fastqdir,"{name}_R1.fastq.gz"),
		R2= join(fastqdir,"{name}_R2.fastq.gz"),
	params:
		out = join(workpath,"{name}","star_out"),
		STARgenome = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/STAR_2.7.8a",
	output:
		G_bam = join(workpath,"{name}","star_out","{name}.star.bam"),
		T_bam = join(workpath,"{name}","star_out","{name}.star_transcriptome.bam")

	container: 'docker://nciccbr/ncigb_star_v2.7.10a:latest'

	shell: """
	
	cd {params.out}
	STAR --genomeDir {params.STARgenome} \
		--readFilesIn  {input} \
		--readFilesCommand zcat \
		--runThreadN {threads} \
		--twopassMode Basic \
		--outSAMunmapped Within \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 12 \
		--alignSJDBoverhangMin 10 \
		--alignMatesGapMax 100000 \
		--chimSegmentReadGapMax 3 \
		--outFilterMismatchNmax 2 \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--outBAMsortingThreadN 6 \
		--limitBAMsortRAM 80000000000
	mv *Aligned.sortedByCoord.out.bam {output.G_bam}
	mv *Aligned.toTranscriptome.out.bam {output.T_bam}
	"""

rule rsem:
	input: join(workpath,"{name}","star_out","{name}.star_transcriptome.bam")
	params:
		out = join(workpath,"{name}","rsem_out"),
		ref = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/rsem_1.3.2/rsem_1.3.2",
		name = "{name}"
	output:
		genes = join(workpath,"{name}","rsem_out","{name}.genes.results"),
		isoform = join(workpath,"{name}","rsem_out","{name}.isoforms.results"),
	container: 'docker://nciccbr/ccbr_rsem_1.3.3:v1.0'

	shell: """

	cd {params.out}
	rsem-calculate-expression --no-bam-output --paired-end -p {threads}  --estimate-rspd  --bam {input} {params.ref} {params.name}
	"""


