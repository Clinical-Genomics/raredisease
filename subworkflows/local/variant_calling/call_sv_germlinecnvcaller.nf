#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PREPROCESSINTERVALS           } from '../../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_ANNOTATEINTERVALS             } from '../../../modules/nf-core/gatk4/annotateintervals/main.nf'
include { GATK4_FILTERINTERVALS               } from '../../../modules/nf-core/gatk4/filterintervals/main.nf'
include { GATK4_INTERVALLISTTOOLS             } from '../../../modules/nf-core/gatk4/intervallisttools/main.nf'
include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main.nf'

workflow CALL_SV_GERMLINECNVCALLER {
    take:
        ch_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
        ch_fasta          // channel: [ val(meta), path(ch_fasta_no_meta) ]
        ch_fai            // channel: [ val(meta), path(ch_fai) ]
        ch_target_bed     // channel: [ val(meta), path(bed), path(tbi) ]
        ch_blacklist_bed  // channel: [ val(meta), path(ch_blacklist_bed) ]
        ch_dict           // channel: [ val(meta), path(ch_dict) ]
        ch_ploidy_model   // channel: [ path(ch_ploidy_model) ]
        ch_gcnvcaller_model      // channel: [ path(ch_gcnvcaller_model) ]

    main:
        ch_versions = Channel.empty()

        input = ch_bam_bai.combine( ch_target_bed.collect{ it[1] } )

        GATK4_COLLECTREADCOUNTS ( input, ch_fasta, ch_fai, ch_dict )

        dgcp_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
                .groupTuple()
                .map({ meta, tsv -> return [meta, tsv, [], [] ]})
        GATK4_DETERMINEGERMLINECONTIGPLOIDY ( dgcp_case_input, [], ch_ploidy_model )

        gcnvc_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
                .groupTuple()
                .map({ meta, tsv -> return [meta, tsv, [] ]})
        GATK4_GERMLINECNVCALLER ( gcnvc_case_input, ch_gcnvcaller_model, GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls.collect{ it[1] } )

        GATK4_POSTPROCESSGERMLINECNVCALLS ( GATK4_GERMLINECNVCALLER.out.calls, ch_gcnvcaller_model, GATK4_GERMLINECNVCALLER.out.calls.collect{ it[1] } )

        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
        ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)
        ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions)
        ch_versions = ch_versions.mix(GATK4_POSTPROCESSGERMLINECNVCALLS.out.versions)

    emit:
        genotyped_intervals_vcf = GATK4_POSTPROCESSGERMLINECNVCALLS.out.intervals  // channel: [ val(meta), path(*.tar.gz) ]
        genotyped_segments_vcf  = GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments   // channel: [ val(meta), path(*.tar.gz) ]
        denoised_vcf            = GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised   // channel: [ val(meta), path(*.tar.gz) ]
        versions                = ch_versions                                      // channel: [ versions.yml ]
}
