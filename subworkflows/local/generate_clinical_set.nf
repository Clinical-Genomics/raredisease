//
// Generate clinical set of variants
//

include { ENSEMBLVEP_FILTERVEP } from '../../modules/nf-core/ensemblvep/filtervep'
include { TABIX_BGZIP          } from '../../modules/nf-core/tabix/bgzip'
include { BCFTOOLS_FILTER      } from '../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_PLUGINSETGT } from '../../modules/nf-core/bcftools/pluginsetgt/main'

workflow GENERATE_CLINICAL_SET {
    take:
        ch_vcf      // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_hgnc_ids // channel: [mandatory] [ val(hgnc_ids) ]
        val_ismt    // value: if mitochondria, set to true

    main:
        ch_versions = Channel.empty()

        ENSEMBLVEP_FILTERVEP(
            ch_vcf,
            ch_hgnc_ids
        )
        .output
        .set { ch_filtervep_out }

        if (val_ismt) {
            BCFTOOLS_FILTER (ch_clin_research_vcf.clinical.map { meta, vcf -> return [meta, vcf, []]})
            ch_clinical_filtered = BCFTOOLS_FILTER.out.vcf
            ch_versions = ch_versions.mix( BCFTOOLS_FILTER.out.versions )

            BCFTOOLS_PLUGINSETGT (
                ch_clinical_filtered.mix(ch_clin_research_vcf.research).map { meta, vcf -> return [meta, vcf, []] },
                Channel.value('q'),
                Channel.value("c:'1/1'"),
                [],
                []
            )
            ch_clinical = BCFTOOLS_PLUGINSETGT.out.vcf.filter { meta, vcf -> meta.set == "clinical" }
            ch_research = BCFTOOLS_PLUGINSETGT.out.vcf.filter { meta, vcf -> meta.set == "research" }
            ch_versions = ch_versions.mix( BCFTOOLS_PLUGINSETGT.out.versions )
        } else {
            ENSEMBLVEP_FILTERVEP(
                ch_clin_research_vcf.clinical,
                ch_hgnc_ids
            )
            .output
            .set { ch_filtervep_out }

            TABIX_BGZIP( ch_filtervep_out )
            ch_clinical = TABIX_BGZIP.out.output

            ch_research = ch_clin_research_vcf.research

            ch_versions = ch_versions.mix( ENSEMBLVEP_FILTERVEP.out.versions )
            ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )
        }

        ch_versions = ch_versions.mix( ENSEMBLVEP_FILTERVEP.out.versions )

    emit:
        vcf      = ch_clinical // channel: [ val(meta), path(vcf) ]
        versions = ch_versions // channel: [ path(versions.yml) ]
}
