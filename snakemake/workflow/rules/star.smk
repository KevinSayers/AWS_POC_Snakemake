rule star:
	input: getfastqc_inputs
	params:
		out = join(workpath,"star_out"),
		STARgenome = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/STAR_2.7.8a",
	output:
		G_bam = expand(join(workpath,"star_out","{name}.star.bam"),name=samples),
		T_bam = expand(join(workpath,"star_out","{name}.star_transcriptome.bam"),name=samples)

	container: 'docker://nciccbr/ncigb_star_v2.7.10a:latest'

	shell: """
	
	cd {params.out}
	STAR --genomeDir {params.STARgenome} \
		--readFilesIn  {input} \
		--readFilesCommand  zcat \
		--runThreadN ${{threads}} \
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


