#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// For general use
include { align; locus_collector; dedup } from './modules/general/main.nf'

// dnascope
include { dnascope_align; dnascope_call_variants; dnascope_model; dnascope_genotype_gvcfs  } from './modules/dnascope/main.nf'

// dnaseq
include { dnaseq_get_metrics; dnaseq_bqsr_table; dnaseq_bqsr_bam; dnaseq_call_variants; dnaseq_genotype_gvcfs; dna_seq_vqsr_snps; dna_seq_vqsr_indels } from './modules/dnaseq/main.nf'

// Get params
workflow = params.workflow
sentieon_license = params.sentieon_license
sentieon_libjemalloc = params.sentieon_libjemalloc
sentieon_bam_option = params.sentieon_bam_option
sentieon_threads = params.sentieon_threads
sentieon_bwa_model = file(params.sentieon_bwa_model, type: 'file')
sentieon_dnascope_model = file(params.sentieon_dnascope_model, type: 'file')

pcr_free = params.pcr_free

if (workflow == 'dnascope-align-call' || params.workflow == 'dnaseq-align-call' ) {
    Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ "${row.SampleID}",
                       "${row.Gender}",
                       "${row.FastqR1}",
                       "${row.FastqR2}",
                       "${row.Flowcell}",
                       "${row.Lane}",
                       "${row.BAM}",
                       "${row.gVCF}" ] }
        .set { samples }
    if (workflow == 'dnascope-align-call'){
                samples.map { [ it[0], it[2], it[3], it[4], it[5]] }
            .set { samples_dnascope_align_call }
    }
    if (workflow == 'dnaseq-align-call'){
                samples.map { [ it[0], it[2], it[3], it[4], it[5]] }
            .set { samples_dnaseq_align_call }
    }  
}

if (workflow == 'dnascope-call' || params.workflow == 'dnaseq-call' ) {
    Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ "${row.SampleID}",
                       "${row.Gender}",
                       "${row.FastqR1}",
                       "${row.FastqR2}",
                       "${row.Flowcell}",
                       "${row.Lane}",
                       "${row.BAM}",
                       "${row.gVCF}" ] }
        .set { samples }
    if (workflow == 'dnascope-call'){
                  samples.map { 
		    if(it[6].endsWith('cram')){
                      return [ it[0], it[6], it[6].concat(".crai") ]
		    }  else if (it[6].endsWith('bam')){
                     return [ it[0], it[6], it[6].concat(".bai") ]
                    }
	          }
            .set { samples_dnascope_call }
    }
    if (workflow == 'dnaseq-call'){
                samples.map { [ it[0], it[6] ] }
            .set { samples_dnaseq_call }
    }         
}

if (workflow == 'dnascope-genotype-gvcfs' || workflow == 'dnaseq-genotype-gvcfs'){
    gvcf_list = Channel.fromPath(params.gvcf_list).toList()
 }

workflow dnascope_align_call {
    take:
        samples

    main:
        dnascope_align(samples,sentieon_bwa_model,sentieon_license,sentieon_libjemalloc,sentieon_bam_option,sentieon_threads)
        locus_collector(dnascope_align.out.raw_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dedup(dnascope_align.out.raw_bam.join(locus_collector.out.score_info),sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dnascope_call_variants(dedup.out.dedup_bam,sentieon_dnascope_model,sentieon_license,sentieon_libjemalloc,sentieon_threads,pcr_free)
        dnascope_model(dnascope_call_variants.out.call_vcf,sentieon_dnascope_model,sentieon_license,sentieon_libjemalloc,sentieon_threads)       
}

workflow dnaseq_align_call {
    take:
        samples

    main:
        align(samples,sentieon_license,sentieon_libjemalloc,sentieon_bam_option,sentieon_threads)
        dnaseq_get_metrics(align.out.raw_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
        locus_collector(align.out.raw_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dedup(align.out.raw_bam.join(locus_collector.out.score_info),sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dnaseq_bqsr_table(dedup.out.dedup_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dnaseq_bqsr_bam(dedup.out.dedup_bam.join(dnaseq_bqsr_table.out.bqsr_table),sentieon_license,sentieon_libjemalloc,sentieon_threads)
        dnaseq_call_variants(dnaseq_bqsr_bam.out.bqsr_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
}

workflow dnascope_call {
    take:
        samples

    main:
        dnascope_call_variants(samples,sentieon_dnascope_model,sentieon_license,sentieon_libjemalloc,sentieon_threads,pcr_free)
        dnascope_model(dnascope_call_variants.out.call_vcf,sentieon_dnascope_model,sentieon_license,sentieon_libjemalloc,sentieon_threads)       
}

workflow dnaseq_call {
    take:
        samples

    main:
        dnaseq_call_variants(samples,sentieon_license,sentieon_libjemalloc,sentieon_threads)
}

workflow dnascope_gg {
    take:
        gvcf_list

    main:
        dnascope_genotype_gvcfs(gvcf_list,sentieon_license,sentieon_libjemalloc,sentieon_bam_option,sentieon_threads)        
}

workflow dnaseq_gg {
    take:
        gvcf_list

    main:
        dnaseq_genotype_gvcfs(gvcf_list,sentieon_license,sentieon_libjemalloc,sentieon_bam_option,sentieon_threads)
        dna_seq_vqsr_snps(dnaseq_genotype_gvcfs.out.vcf,sentieon_license,sentieon_libjemalloc,sentieon_threads)       
        dna_seq_vqsr_indels(dna_seq_vqsr_snps.out.snps_vcf,sentieon_license,sentieon_libjemalloc,sentieon_threads)       
}

workflow {
    switch (workflow) {
            
        case['dnascope-align-call']:
            dnascope_align_call(samples_dnascope_align_call)
            break

        case['dnaseq-align-call']:
            dnaseq_align_call(samples_dnaseq_align_call)
            break

        case['dnascope-call']:
            dnascope_call(samples_dnascope_call)
            break

        case['dnascope-genotype-gvcfs']:
            dnascope_gg(gvcf_list)
            break
        
        case['dnaseq-genotype-gvcfs']:
            dnaseq_gg(gvcf_list)
            break
          
        
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
        
        workflow.onComplete {
            println ( workflow.success ? """
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            """ : """
            Failed: ${workflow.errorReport}
            exit status : ${workflow.exitStatus}
            """
            )
        }
    }
}
