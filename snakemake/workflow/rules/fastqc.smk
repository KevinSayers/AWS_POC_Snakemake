

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
		out = join(workpath,"fastqc")
	output:
		expand(join(workpath,"fastqc","{name}_R1_fastqc.html"),name=samples),
		expand(join(workpath,"fastqc","{name}_R2_fastqc.html"),name=samples)

	envmodules:
	"fastqc/0.11.9"

	container: "docker://nciccbr/ccbr_fastqc_0.11.9:v1.1"
		
	shell: """

	fastqc {input} -o  {params.out}

	"""

