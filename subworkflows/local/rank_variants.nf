//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE                                       } from '../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS                                         } from '../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE                                          } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_SCORE as GENMOD_SCORE_FOR_GICAM                } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND                                       } from '../../modules/nf-core/genmod/compound/main'
include { MIVMIR_INFER                                          } from '../../modules/local/mivmir/main'
include { GICAM_INFER                                           } from '../../modules/local/gicam/main'
include { BCFTOOLS_SORT                                         } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP                                           } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_GENMOD_GICAM               } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GICAM            } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX                                           } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_MERGE_GENMOD_GICAM      } from '../../modules/nf-core/bcftools/annotate/main'

workflow RANK_VARIANTS {

    take:
        ch_vcf                          // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_pedfile                      // channel: [mandatory] [ path(ped) ]
        ch_reduced_penetrance           // channel: [mandatory] [ path(pentrance) ]
        ch_score_config                 // channel: [mandatory] [ path(ini) ]
        ch_genmod_gicam_score_config    // channel: [mandatory] [ path(ini) ]
        rank_with_mivmir_gicam          // value

    main:
        ch_versions = Channel.empty()

        GENMOD_ANNOTATE(ch_vcf)

        ch_models_in = GENMOD_ANNOTATE.out.vcf.combine(ch_pedfile)

        GENMOD_MODELS(ch_models_in, ch_reduced_penetrance)

        ch_score_in = GENMOD_MODELS.out.vcf.combine(ch_pedfile)

        GENMOD_SCORE(ch_score_in, ch_score_config)

        // Run MIVMIR - GICAM scoring (not supported for MT SNVs and SVs)
        if (rank_with_mivmir_gicam) {
            GENMOD_SCORE_FOR_GICAM(ch_score_in, ch_genmod_gicam_score_config)
            MIVMIR_INFER(GENMOD_SCORE_FOR_GICAM.out.vcf)
            GICAM_INFER(MIVMIR_INFER.out.vcf)
            TABIX_BGZIPTABIX_GICAM(GICAM_INFER.out.vcf)
        }

        GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)
        TABIX_BGZIP(GENMOD_COMPOUND.out.vcf) //run only for SNVs

        // Merge Genmod and MIVMIR-GICAM scores
        if (rank_with_mivmir_gicam) {
            TABIX_BGZIP.out.output
                .join(TABIX_BGZIPTABIX_GICAM.out.gz_tbi, failOnMismatch: true)
                .map {meta, vcf_genmod, vcf_gicam, vcf_index_gicam -> return [ meta, vcf_genmod, [], vcf_gicam, vcf_index_gicam ]}
                .set {ch_merge_genmod_gicam}
            BCFTOOLS_MERGE_GENMOD_GICAM(ch_merge_genmod_gicam, [])
            TABIX_BGZIP_GENMOD_GICAM(BCFTOOLS_MERGE_GENMOD_GICAM.out.vcf)
        }

        BCFTOOLS_SORT(GENMOD_COMPOUND.out.vcf) // SV file needs to be sorted before indexing
        // Mix SNVs and SVs
        if (rank_with_mivmir_gicam) {
            ch_vcf = TABIX_BGZIP_GENMOD_GICAM.out.output.mix(BCFTOOLS_SORT.out.vcf)
        } else {
            ch_vcf = TABIX_BGZIP.out.output.mix(BCFTOOLS_SORT.out.vcf)
        }

        TABIX_TABIX (ch_vcf)

        ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
        ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
        vcf      = ch_vcf       // channel: [ val(meta), path(vcf) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
