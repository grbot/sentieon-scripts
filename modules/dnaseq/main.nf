#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get params
ref = file(params.ref, type: 'file')
dbsnp = file(params.dbsnp, type: 'file') 
known_indels_1 = file(params.known_indels_1, type: 'file')
known_indels_2 = file(params.known_indels_2, type: 'file')
call_conf = params.call_conf
emit_conf = params.emit_conf

process dnaseq_get_metrics {
    tag { "${params.project_name}.${sample_id}.dnaseq_get_metrics" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
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

process dnaseq_bqsr_table {
    tag { "${params.project_name}.${sample_id}.dnaseq_bqsr-table" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), file(bam_file), file(bam_file_index)

        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple val("$sample_id"), file("${sample_id}.dedup.bqsr.table"), emit: bqsr_table 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -r ${ref} \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo QualCal \
        -k ${dbsnp} -k ${known_indels_1} -k ${known_indels_2} \
        ${sample_id}.dedup.bqsr.table
       """ 
}

process dnaseq_bqsr_bam {
    tag { "${params.project_name}.${sample_id}.dnaseq_bqsr-bam" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), file(bam_file), file(bam_file_index), file(table_file) 

        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple val("$sample_id"), file("${sample_id}.dedup.bqsr.bam"), file("${sample_id}.dedup.bqsr.bam.bai"), emit: bqsr_bam 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -r ${ref} \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        -q ${table_file} \
        --algo ReadWriter \
        ${sample_id}.dedup.bqsr.bam 
       """ 
}

process dnaseq_call_variants {
    tag { "${params.project_name}.${sample_id}.dnaseq_call" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
    tuple val(sample_id), file(bam_file), file(bam_file_index)
    val sentieon_license
    val sentieon_libjemalloc
    val sentieon_threads
  
    output:
    tuple val("$sample_id"), file("${sample_id}.dedup.bqsr.vcf.gz"), file("${sample_id}.dedup.bqsr.vcf.gz.tbi"), emit: call_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -r ${ref} \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo Haplotyper \
        -d ${dbsnp} \
        --genotype_model multinomial \
        --emit_conf ${emit_conf} --call_conf ${call_conf} \
        ${sample_id}.dedup.bqsr.vcf.gz
        """ 
}

process dnaseq_genotype_gvcfs {
    tag { "${params.project_name}.dnaseq_genotype_gvcfs" }
    memory { 64.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        val(project_id)
        file(gvcf_list)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        tuple file("${project_id}.vcf.gz"), file("${project_id}.vcf.tbi"), emit: vcf

    script:                     
        """   
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc} 
        cat ${gvcf_list} | \
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        --algo GVCFtyper \
        ${project_id}.vcf.gz  -
        """ 
}

process vqsr_snps {
    tag { "${params.project_name}.${sample_id}.get_metrics" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
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

        
        
        """ 
}

process vqsr_indels { 
    tag { "${params.project_name}.${sample_id}.locus_collector" } 
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false 
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
        
        
        """  
} 