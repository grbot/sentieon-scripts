#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get params
ref = file(params.ref, type: 'file')
dbsnp = file(params.dbsnp, type: 'file') 
known_indels_1 = file(params.known_indels_1, type: 'file')
known_indels_2 = file(params.known_indels_2, type: 'file')
mills = file(params.mills, type: 'file')
omni = file(params.omni, type: 'file')
hapmap = file(params.hapmap, type: 'file')
phase1 = file(params.phase1, type: 'file')
phase1_indel = file(params.phase1_indel, type: 'file')
call_conf = params.call_conf
emit_conf = params.emit_conf
max_gaussians = params.max_gaussians

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
    tuple val("$sample_id"), file("${sample_id}.dedup.bqsr.gvcf.gz"), file("${sample_id}.dedup.bqsr.gvcf.gz.tbi"), emit: call_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -r ${ref} \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo Haplotyper \
        --emit_mode gvcf \
        -d ${dbsnp} \
        --genotype_model multinomial \
        --emit_conf ${emit_conf} --call_conf ${call_conf} \
        ${sample_id}.dedup.bqsr.gvcf.gz
        """ 
}

process dnaseq_genotype_gvcfs {
    tag { "${params.project_name}.dnaseq_genotype_gvcfs" }
    memory { 64.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        file(gvcf_list)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        tuple file("${params.project_name}.vcf.gz"), file("${params.project_name}.vcf.gz.tbi"), emit: vcf

    script:                     
        """   
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc} 
        cat ${gvcf_list} | \
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        --algo GVCFtyper \
        ${params.project_name}.vcf.gz -
        """ 
}

process dna_seq_vqsr_snps {
    tag { "${params.project_name}.dna_seq_vqsr_snps" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple file(vcf), file(vcf_index)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple file("${params.project_name}.vqsr.snps.vcf.gz"), file("${params.project_name}.vqsr.snps.vcf.gz.tbi"), emit: snps_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc} 

        resource_text="--resource ${phase1} \
        --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
        resource_text="\$resource_text --resource $omni \
        --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
        resource_text="\$resource_text --resource $dbsnp \
        --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
        resource_text="\$resource_text --resource $hapmap \
        --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"

        # Create the annotation argument
        annotate_text=""
        annotation_array="DP QD FS SOR MQRankSum ReadPosRankSum"
        for annotation in \$annotation_array; do
         annotate_text="\$annotate_text --annotation \$annotation"
        done

        sentieon driver \
        -r ${ref} \
        --algo VarCal \
        -v ${vcf} \
        \$resource_text \$annotate_text \
        --var_type SNP \
        --plot_file ${params.project_name}.vqsr.snps.plotfile \
        --tranche 90.0 \
        --tranche 99.0 \
        --tranche 99.5 \
        --tranche 99.9 \
        --tranche 100.0 \
        --max_gaussians ${max_gaussians} \
        --tranches_file ${params.project_name}.vqsr.snps.tranches \
        ${params.project_name}.vqsr.snps.recal

        # Apply VQSR
        sentieon driver \
        -r ${ref} \
        --algo ApplyVarCal \
        -v ${vcf} \
        --var_type SNP --tranches_file ${params.project_name}.vqsr.snps.tranches \
        --sensitivity 99.0 \
        --recal ${params.project_name}.vqsr.snps.recal ${params.project_name}.vqsr.snps.vcf.gz
        
        # Plot the report
        sentieon plot \
        vqsr -o ${params.project_name}.vqsr.snps.pdf ${params.project_name}.vqsr.snps.plotfile        
        """ 
}

process dna_seq_vqsr_indels {
    tag { "${params.project_name}.dna_seq_vqsr_indels" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple file(vcf), file(vcf_index)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple file("${params.project_name}.vqsr.snps.indels.vcf.gz"), file("${params.project_name}.vqsr.snps.indels.vcf.gz.tbi"), emit: indel_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}

        resource_text="--resource $mills \
        --resource_param mills,known=false,training=true,truth=false,prior=12.0 "
        resource_text="\$resource_text --resource $dbsnp \
        --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "

        # Create the annotation argument        
        annotation_array="QD DP FS SOR MQRankSum ReadPosRankSum"
        annotate_text=""
        for annotation in \$annotation_array; do
          annotate_text="\$annotate_text --annotation \$annotation"
        done

        # Run VQSR
        sentieon driver \
        -r ${ref} \
        --algo VarCal \
        -v ${vcf} \
        \$resource_text \$annotate_text --var_type INDEL \
        --plot_file ${params.project_name}.vqsr.snps.indels.plotfile \
        --tranche 90.0 \
        --tranche 99.0 \
        --tranche 99.5 \
        --tranche 99.9 \
        --tranche 100.0 \
        --max_gaussians ${max_gaussians} \
        --tranches_file ${params.project_name}.vqsr.snps.indels.tranches ${params.project_name}.vqsr.snps.indels.recal

        # Apply the VQSR
        sentieon driver \
        -r ${ref} \
        --algo ApplyVarCal \
        -v ${vcf} \
        --var_type INDEL \
        --recal ${params.project_name}.vqsr.snps.indels.recal \
        --tranches_file ${params.project_name}.vqsr.snps.indels.tranches \
        --sensitivity 99.5 ${params.project_name}.vqsr.snps.indels.vcf.gz
   
        # Plot the report
        sentieon plot \
        vqsr -o ${params.project_name}.vqsr.snps.indels.pdf \
        ${params.project_name}.vqsr.snps.indels.plotfile
        """
}
