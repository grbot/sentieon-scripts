#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get params
ref = file(params.ref, type: 'file')
dbsnp = file(params.dbsnp, type: 'file')

process dnascope_align {
    tag { "${params.project_name}.${sample_id}.dnascope_align" }
    memory { 128.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), path(fastq_r1_file), path(fastq_r2_file), val(flowcell), val(lane)
        path sentieon_bwa_model
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        tuple val("$sample_id"), file("${sample_id}.cram"), file("${sample_id}.cram.bai"), emit: raw_bam

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
        -x ${sentieon_bwa_model} \
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

process dnascope_call_variants {
    tag { "${params.project_name}.${sample_id}.dnascope_call" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
    tuple val(sample_id), path(bam_file), path(bam_file_index)
    path sentieon_dnascope_model
    val sentieon_license
    val sentieon_libjemalloc
    val sentieon_threads
    val pcr_free
  
    output:
    tuple val("$sample_id"), file("${sample_id}.dedup.gvcf.gz"), file("${sample_id}.dedup.gvcf.gz.tbi"), emit: call_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        
        if [ "${pcr_free}" = true ]; then
            sentieon driver \
            -t ${sentieon_threads} \
            -i ${bam_file} \
            -r ${ref} \
            --algo DNAscope \
            -d ${dbsnp} \
            --pcr_indel_model none \
            --model ${sentieon_dnascope_model} \
            --emit_mode gvcf \
            ${sample_id}.dedup.gvcf.gz
        else
            sentieon driver \
            -t ${sentieon_threads} \
            -i ${bam_file} \
            -r ${ref} \
            --algo DNAscope \
            -d ${dbsnp} \
            --model ${sentieon_dnascope_model} \
            --emit_mode gvcf \
            ${sample_id}.dedup.gvcf.gz
        fi
        """ 
}

process dnascope_model {
    tag { "${sample_id}.dnascope_model" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
        tuple val(sample_id), file(gvcf), file(gvcf_index)
        file sentieon_dnascope_model
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple file("${sample_id}.dedup.model.gvcf.gz"), file("${sample_id}.dedup.model.gvcf.gz.tbi"), emit: model_gvcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        --algo DNAModelApply \
        --model ${sentieon_dnascope_model} \
        -v ${gvcf} \
        ${sample_id}.dedup.model.gvcf.gz
        """ 
}


process dnascope_genotype_gvcfs {
    tag { "${params.project_name}.dnascope_genotype_gvcfs" }
    memory { 64.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        path gvcf_list
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

