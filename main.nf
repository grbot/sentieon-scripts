#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// For alignment and calling
include { align } from './modules/dnascope-align-call/main.nf'
include { locus_collector } from './modules/dnascope-align-call/main.nf'
include { dedup } from './modules/dnascope-align-call/main.nf'
include { call } from './modules/dnascope-align-call/main.nf'
include { model } from './modules/dnascope-align-call/main.nf'

// Get params
ref = Channel.fromPath(params.ref).toList()
dbsnp = Channel.fromPath(params.dbsnp).toList()
workflow = params.workflow
sentieon_license = params.sentieon_license
sentieon_libjemalloc = params.sentieon_libjemalloc
sentieon_bam_option = params.sentieon_bam_option
sentieon_threads = params.sentieon_threads
sentieon_dnascope_model = params.sentieon_dnascope_model

pcr_free = params.pcr_free

// Get info from sample sheets for workflows (currently having only one workflow)
if (params.workflow == 'dnascope-align-call') {
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
    if (params.workflow == 'dnascope-align-call'){
                samples.map { [ it[0], it[2], it[3], it[4], it[5]] }
            .set { samples_dnascope_align_call }
    }
}

workflow dnascope_align_call {
    take:
        samples

    main:
        align(samples,ref,sentieon_license,sentieon_libjemalloc,sentieon_bam_option,sentieon_threads)
        //locus_collector(align.out.raw_bam,sentieon_license,sentieon_libjemalloc,sentieon_threads)
        //dedup(align.out.raw_bam.mix(locus_collector.out.score_info).groupTuple(),sentieon_license,sentieon_libjemalloc,sentieon_threads)
        
        //call(dedup.out.dedup_bam,ref,dbsnp,sentieon_dnascope_model,sentieon_license,sentieon_libjemalloc,sentieon_threads,pcr_free)
        //model(call.out.call_vcf,sentieon_license,sentieon_libjemalloc,sentieon_threads)
}

workflow {
    switch (workflow) {
            
        case['dnascope-align-call']:
            dnascope_align_call(samples_dnascope_align_call)
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
