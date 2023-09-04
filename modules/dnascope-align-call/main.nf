process align {
    tag { "${params.project_name}.${params.protocol}.align.${name}" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        set val(sample_id), file(fastq_r1_file), file(fastq_r2_file), val(flowcell), val(lane)
        file ref
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_bam_option
        val sentieon_threads

    output:
        set val("$sample_id"), file("${sample_id}.bam"), file("${sample_id}.bam.bai")  into raw_bam

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
        -R \"${readgroup_info}\" \
        -t ${sentieon_threads}  \
        -K 100000000 \
        -Y \
        ${ref} \
        ${fastq_r1_file} \
        ${fastq_r2_file} | \
        sentieon util sort \
        ${sentieon_bam_option} \
        -t ${sentieon_threads} \
        -r ${ref} \
        -o  ${sample_id}.bam \
        --sam2bam \
        -i -
        """ 
}

process locus_collector {
    tag { "${params.project_name}.${params.protocol}.locus_collector" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        set val(sample_id), file(bam_file), file(bam_file_index)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        set val(sample_id), file("${sample_id}.score.txt")  into score_info

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
    tag { "${params.project_name}.${params.protocol}.dedup" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
        set val(sample_id), file(bam_file), file(bam_file_index), file(score_info)
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads
  
    output:
        set val("$sample_id"), file("${sample_id}.dedup.bam") into dedup_bam 

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

process call {
    tag { "${params.project_name}.${params.protocol}.call" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'sentieon'

    input:
    set val(sample_id), file(bam_file)
    file ref
    file dbsnp
    file sentieon_dnascope_model
    val sentieon_license
    val sentieon_libjemalloc
    val sentieon_threads
    val pcr_free
  
    output:
    set val($sample_id), file("${sample_id}.dedup.vcf.gz") into call_vcf 

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
    input:
        set val(sample_id), file(vcf_file)
        file ref
        file sentieon_dnascope_model
        val sentieon_license
        val sentieon_libjemalloc
        val sentieon_threads

  
    output:
        set val($sample_id), file("${sample_id}.dedup.model.vcf.gz") into model_vcf 

    script:
        """
        export SENTIEON_LICENSE=${sentieon_license}
        export LD_PRELOAD=${sentieon_libjemalloc}
        
        sentieon driver \
        -t ${sentieon_threads} \
        -i ${bam_file} \
        -r ${ref} \
        --algo DNAModelApply \
        --model ${sentieon_dnascope_model} \
        -v ${vcf_file}
        ${sample_id}.dedup.model.vcf.gz       
        """ 
}

