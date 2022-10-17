rule rsem:
	input: expand(join(workpath,"star_out","{name}.star_transcriptome.bam"),name=samples)
	params:
		out = join(workpath,"rsem_out"),
		ref = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/rsem_1.3.2/rsem_1.3.2",
		name = expand("{name}",name=samples)
	output:
		genes = expand(join(workpath,"rsem_out","{name}.genes.results"),name=samples),
		isoform = expand(join(workpath,"rsem_out","{name}.isoforms.results"),name=samples)
	container: 'docker://nciccbr/ccbr_rsem_1.3.3:v1.0'

	shell: """

	cd {params.out}
	rsem-calculate-expression --no-bam-output --paired-end -p {threads}  --estimate-rspd  --bam {input} {params.ref} {params.name}
	"""

