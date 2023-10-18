#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Get params

ref = file(params.ref, type: 'file')
dbsnp = file(params.dbsnp, type: 'file') 
sentieon_license = params.sentieon_license
sentieon_libjemalloc = params.sentieon_libjemalloc
sentieon_bam_option = params.sentieon_bam_option
sentieon_threads = params.sentieon_threads
sentieon_dnascope_model = params.sentieon_dnascope_model
pcr_free = params.pcr_free


process align {
    tag { "${params.project_name}.${sample_id}.align" }
    memory { 64.GB * task.attempt }
    cpus { "${sentieon_threads}" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), path(fastq_r1_file), path(fastq_r2_file), val(flowcell), val(lane)
        //file ref
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        tuple val("$sample_id"), file("${sample_id}.bam"), file("${sample_id}.bam.bai"), emit: raw_bam

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
        -o ${sample_id}.bam \
        -t ${sentieon_threads} \
        --sam2bam \
        -i -

        """ 
}

process locus_collector {
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
        sentieon driver \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo LocusCollector \
        --fun score_info ${sample_id}.score.txt
        """ 
}

process dedup {
    tag { "${params.project_name}.${sample_id}.dedup" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        tuple val(sample_id), file(bam_file), file(bam_file_index), file(score_info), file(score_info_index)

        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        tuple val("$sample_id"), file("${sample_id}.dedup.bam"), file("${sample_id}.dedup.bam.bai"), emit: dedup_bam 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        sentieon driver \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        --algo Dedup \
        --rmdup \
        --score_info ${score_info} \
        --metrics ${sample_id}.dedup_metrics.txt \
        ${sample_id}.dedup.bam
        """ 
}

process call_variants {
    tag { "${params.project_name}.${sample_id}.call" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
    tuple val(sample_id), file(bam_file), file(bam_file_index)
    file sentieon_dnascope_model
    val sentieon_license
    val sentieon_libjemalloc
    val sentieon_threads
    val pcr_free
  
    output:
    tuple val("$sample_id"), file("${sample_id}.dedup.vcf.gz"), file("${sample_id}.dedup.vcf.gz.tbi"), emit: call_vcf 

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
            ${sample_id}.dedup.vcf.gz
        else
            sentieon driver \
            -t ${sentieon_threads} \
            -i ${bam_file} \
            -r ${ref} \
            --algo DNAscope \
            -d ${dbsnp} \
            --model ${sentieon_dnascope_model} \
            ${sample_id}.dedup.vcf.gz
        fi
        """ 
}

process model {
    tag { "${params.project_name}.${sample_id}.model" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
        tuple val(sample_id), file(vcf_file), file(vcf_file_index)
        file sentieon_dnascope_model
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads

  
    output:
        tuple val("$sample_id"), file("${sample_id}.dedup.model.vcf.gz"), emit: model_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        
        sentieon driver \
        -t ${sentieon_threads} \
        -r ${ref} \
        --algo DNAModelApply \
        --model ${sentieon_dnascope_model} \
        -v ${vcf_file} \
        ${sample_id}.dedup.model.vcf.gz       
        """ 
}

