

def getfastqc_inputs(wildcards):
	fastq =[]
	for s in samples:
		fastq.append(join(fastqdir,s + "_R1.fastq.gz"))
		fastq.append(join(fastqdir,s + "_R2.fastq.gz"))
#		print(fastq)
	return fastq


rule fastqc:
	input:
		getfastqc_inputs
	params:
		out = expand(join(workpath,"{name}","fastqc"),name=samples)
	output:
		expand(join(workpath,"{name}","fastqc","{name}_R1_fastqc.zip"),name=samples),
		expand(join(workpath,"{name}","fastqc","{name}_R2_fastqc.zip"),name=samples)

	envmodules:
	"fastqc/0.11.9"

	container: "docker://nciccbr/ccbr_fastqc_0.11.9:v1.1"
		
	shell: """

	fastqc {input} -o  {params.out}

	"""

rule multiqc:
	input:
		join(workpath,"{name}","fastqc","{name}_R1_fastqc.zip"),
		join(workpath,"{name}","fastqc","{name}_R2_fastqc.zip"),

	params:
		out = join(workpath,"{name}","multiqc")
	output:
		join(workpath,"{name}","multiqc","multiqc_report.html"),

	envmodules:
	"multiqc/1.9"

	container: "docker://nciccbr/ccbr_multiqc_1.9:v0.0.1"
		
	shell: """
	mkdir -p {params.out}
	multiqc {input} -o {params.out}

	"""

