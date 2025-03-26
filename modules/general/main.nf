#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get params
ref = file(params.ref, type: 'file')
dbsnp = file(params.dbsnp, type: 'file') 
call_conf = params.call_conf
emit_conf = params.emit_conf

process align {
    tag { "${params.project_name}.${sample_id}.align" }
    memory { 128.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), path(fastq_r1_file), path(fastq_r2_file), val(flowcell), val(lane)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        tuple val("$sample_id"), file("${sample_id}.cram"), file("${sample_id}.cram.crai"), emit: raw_bam

    script:
        readgroup_info="@RG\\tID:$flowcell.$lane\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"

        if(lane == "0") {
            sample_id = "$sample_id"
        } else {
            sample_id = "$sample_id-${flowcell}.${lane}"
        }
                         
        """   
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}        
        sentieon bwa mem \
        -M \
        -R \"${readgroup_info}\" \
        -t ${sentieon_threads}  \
        -K 100000000 \
        -Y \
        ${ref} \
        ${fastq_r1_file} \
        ${fastq_r2_file} | \
        sentieon util sort \
        ${sentieon_bam_option} \
        -r ${ref} \
        -o ${sample_id}.cram \
        -t ${sentieon_threads} \
        --sam2bam \
        -i -
        """ 
}

process get_metrics {
    tag { "${params.project_name}.${sample_id}.get_metrics" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    cpus { "${sentieon_threads}" }
    label 'sentieon'

    input:
        tuple val(sample_id), file(bam_file), file(bam_file_index)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple val(sample_id), file("${sample_id}.*.txt"), file("${sample_id}.*.pdf"), emit: metrics

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc} 

        sentieon driver \
        -r $ref \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo MeanQualityByCycle ${sample_id}.mq_metrics.txt \
        --algo QualDistribution ${sample_id}.qd_metrics.txt \
        --algo GCBias --summary ${sample_id}.gc_summary.txt ${sample_id}.gc_metrics.txt \
        --algo AlignmentStat --adapter_seq '' ${sample_id}.aln_metrics.txt \
        --algo InsertSizeMetricAlgo ${sample_id}.is_metrics.txt
	
	    sentieon plot GCBias \
        -o ${sample_id}.gc-report.pdf ${sample_id}.gc_metrics.txt

	    sentieon plot QualDistribution \
        -o ${sample_id}.qd-report.pdf ${sample_id}.qd_metrics.txt
	
	    sentieon plot MeanQualityByCycle \
        -o ${sample_id}.mq-report.pdf ${sample_id}.mq_metrics.txt

	    sentieon plot InsertSizeMetricAlgo \
        -o ${sample_id}.is-report.pdf ${sample_id}.is_metrics.txt
        
        """ 
}

process locus_collector { 
    tag { "${params.project_name}.${sample_id}.locus_collector" } 
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    cpus { "${sentieon_threads}" }
    label 'sentieon' 
 
    input: 
        tuple val(sample_id), file(bam_file), file(bam_file_index) 
        val sentieon_license 
        val sentieon_libjemalloc 
        val sentieon_threads 
   
    output: 
        tuple val(sample_id), file("${sample_id}.score.txt"), file("${sample_id}.score.txt.idx"), emit: score_info 
 
    script: 
        """ 
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        -i ${bam_file} \
        --algo LocusCollector \
        --fun score_info ${sample_id}.score.txt 
        """  
} 

process dedup {
    tag { "${params.project_name}.${sample_id}.dedup" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    cpus { "${sentieon_threads}" }
    label 'sentieon'

    input:
        tuple val(sample_id), file(bam_file), file(bam_file_index), file(score_info), file(score_info_index)

        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple val("$sample_id"), file("${sample_id}.dedup.cram"), file("${sample_id}.dedup.cram.crai"), emit: dedup_bam 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        -i ${bam_file} \
        --algo Dedup \
        --rmdup \
        --score_info ${score_info} \
        --metrics ${sample_id}.dedup_metrics.txt \
        ${sample_id}.dedup.cram
        """ 
}

process cram_to_bam {
    tag { "${params.project_name}.${sample_id}.cram_to_bam" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    cpus { "${sentieon_threads}" }
    label 'samtools'

    input:
        tuple val(sample_id), file(cram_file), file(cram_file_index)
        val sentieon_threads

    output:
        tuple val("$sample_id"), file("${bam_file}"), file("${bam_file_index}"), emit: bam 

    script:
        bam_file = cram_file.getName().replaceFirst(/\.[^.]+$/, ".bam")  //
        bam_file_index = cram_file.getName().replaceFirst(/\.[^.]+$/, ".bam.bai")  //
        """
        samtools \
        view \
        -T ${ref} \
        -o ${bam_file} \
        -@ ${sentieon_threads} \
        ${cram_file}

        samtools \
        index \
        -@ ${sentieon_threads} \
        ${bam_file}
        """ 
}