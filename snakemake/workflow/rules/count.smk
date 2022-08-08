rule star:
	input: getfastqc_inputs
	params:
		out = expand(join(workpath,"{name}","star_out"),name=samples),
		STARgenome = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/STAR_2.7.8a",
	output:
		G_bam = expand(join(workpath,"{name}","star_out","{name}.star.bam"),name=samples),
		T_bam = expand(join(workpath,"{name}","star_out","{name}.star_transcriptome.bam"),name=samples)

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
	input: expand(join(workpath,"{name}","star_out","{name}.star_transcriptome.bam"),name=samples)
	params:
		out = expand(join(workpath,"{name}","rsem_out"),name=samples),
		ref = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/rsem_1.3.2/rsem_1.3.2",
		name = expand("{name}",name=samples)
	output:
		genes = expand(join(workpath,"{name}","rsem_out","{name}.genes.results"),name=samples),
		isoform = expand(join(workpath,"{name}","rsem_out","{name}.isoforms.results"),name=samples)
	container: 'docker://nciccbr/ccbr_rsem_1.3.3:v1.0'

	shell: """

	cd {params.out}
	rsem-calculate-expression --no-bam-output --paired-end -p {threads}  --estimate-rspd  --bam {input} {params.ref} {params.name}
	"""


