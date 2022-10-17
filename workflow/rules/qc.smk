

def getfastqc_inputs(wildcards):
	fastq =[]
	for s in samples:
		fastq.append(join(fastqdir,s + "_R1.fastq.gz"))
		fastq.append(join(fastqdir,s + "_R2.fastq.gz"))
#		print(fastq)
	return fastq


rule fastqc:
	input:
		R1= join(fastqdir,"{name}_R1.fastq.gz"),
		R2= join(fastqdir,"{name}_R2.fastq.gz"),
	params:
		out = join(workpath,"{name}","fastqc")
	output:
		join(workpath,"{name}","fastqc","{name}_R1_fastqc.zip"),
		join(workpath,"{name}","fastqc","{name}_R2_fastqc.zip")

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

