#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// define user inputs

params.reads = 'path to reads'
params.resultsdir = 'path to output folder'


// load params
params.genome_index = "/data/khanlab/projects/ngs_pipeline_testing/index-STAR_2.7.9a"
params.gtf = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/gencode.v37lift37.annotation_ERCC92.gtf"
params.rsem_index = "/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/rsem_1.3.2"


//Print out log
log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         star_index   : ${params.genome_index}
         rsem_index   : ${params.rsem_index}
         reads        : ${params.reads}
         outdir       : ${params.resultsdir}
         """
         .stripIndent()



// set up processes
process cutadapt {
        tag { dataset_id }

        input:
        tuple val(dataset_id), 
        path(forward), 
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}_R1.fastq"),
        path("trim_${dataset_id}_R2.fastq")

        container 'docker://nciccbr/ncigb_cutadapt_v1.18:latest'

        script:
        """
	cutadapt  -o trim_${dataset_id}_R1.fastq -p trim_${dataset_id}_R2.fastq $forward $reverse
        """

}


process fastqc {
        tag { dataset_id }
//	cache false
	input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

	output:
	tuple val("${dataset_id}"),
        path("fastqc_trim_${dataset_id}")

	container 'docker://nciccbr/ccbr_fastqc_0.11.9:v1.1'

	script:
	"""	
	mkdir fastqc_trim_${dataset_id}
	fastqc -o fastqc_trim_${dataset_id} -q $forward $reverse
	"""
}


process star {
	tag { dataset_id }
        input:
        tuple val(dataset_id),
        path(forward), 
        path(reverse),
    	path(genomeIndex),
        path(gtf)
        
        output:
        tuple val("${dataset_id}"),
            path("trim_${dataset_id}Aligned.toTranscriptome.out.bam")

        container 'docker://nciccbr/ncigb_star_v2.7.10a:latest'

        script:
        """
	STAR --genomeDir ${genomeIndex} \
                --readFilesIn $forward $reverse \
                --sjdbGTFfile ${gtf} \
                --runThreadN ${task.cpus} \
                --twopassMode Basic \
                --outSAMunmapped Within \
                --outFileNamePrefix trim_$dataset_id \
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

        """

}

process rsem {
        tag { dataset_id }

	publishDir "$params.resultsdir/$dataset_id", mode: 'move'

	input:
        tuple val(dataset_id),
        path(bam),
        path(genomeIndex)

	output:
        tuple val("${dataset_id}"), path("trim_${dataset_id}.genes.results")
	
	container 'docker://nciccbr/ccbr_rsem_1.3.3:v1.0'

	script:
	"""
	rsem-calculate-expression --no-bam-output --paired-end -p ${task.cpus}  --estimate-rspd  --bam $bam ${genomeIndex}/rsem_1.3.2 trim_$dataset_id
	"""

}

process multiqc {
        tag { dataset_id }
	publishDir "$params.resultsdir/$dataset_id", mode: 'move'
	input:
        tuple val(dataset_id),
        path(qc)
	output:
	path "multiqc_report.html"
//	tuple val("${dataset_id}"), path("trim_${dataset_id}_multiqc_report.html")
	container 'docker://nciccbr/ccbr_multiqc_1.9:v0.0.1'

	script:
	"""
	multiqc -m fastqc .

	"""

}



// invoke workflow:

workflow{
    read_pairs = Channel
        .fromFilePairs(params.reads, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
    gtf = Channel.of(file(params.gtf, checkIfExists:true))
    genomeIndex = Channel.of(file(params.genome_index, checkIfExists:true))
    rsemIndex = Channel.of(file(params.rsem_index, checkIfExists:true))

    cutadapt(read_pairs)
    fastqc(cutadapt.out)
    star(
        cutadapt.out
            .combine(genomeIndex)
            .combine(gtf)
    )
    rsem(
        star.out
            .combine(rsemIndex)
    )
    multiqc(fastqc.out)
}



