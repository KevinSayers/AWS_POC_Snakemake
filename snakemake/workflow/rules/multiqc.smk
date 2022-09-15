
rule multiqc:
	input:
		expand(join(workpath,"fastqc","{name}_R1_fastqc.zip"),name=samples),
		expand(join(workpath,"fastqc","{name}_R2_fastqc.zip"),name=samples)

	params:
		out = join(workpath,"multiqc")
	output:
		expand(join(workpath,"multiqc","multiqc_report.html"),name=samples),

	envmodules:
	"multiqc/1.9"

	container: "docker://nciccbr/ccbr_multiqc_1.9:v0.0.1"
		
	shell: """
	mkdir -p {params.out}
	multiqc {input} -o {params.out}

	"""

