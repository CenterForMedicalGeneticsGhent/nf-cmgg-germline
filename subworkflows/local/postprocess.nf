//
// GENOTYPE
//

include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS  } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_COMBINEGVCFS as COMBINEGVCFS     } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF       } from '../../modules/nf-core/modules/gatk4/reblockgvcf/main'
include { TABIX_TABIX as TABIX_GVCFS             } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS    } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS    } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_PED_VCFS          } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { RTGTOOLS_PEDFILTER as PEDFILTER        } from '../../modules/nf-core/modules/rtgtools/pedfilter/main'
include { MERGE_VCF_HEADERS                      } from '../../modules/local/merge_vcf_headers'
include { BCFTOOLS_FILTER as FILTER_SNPS         } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_INDELS       } from '../../modules/nf-core/modules/bcftools/filter/main'

workflow POST_PROCESS {
    take:
        gvcfs        // channel: [mandatory] gvcfs
        peds         // channel: [mandatory] peds
        fasta        // channel: [mandatory] fasta reference
        fasta_fai    // channel: [mandatory] fasta reference index
        dict         // channel: [mandatory] sequence dictionary

    main:

    post_processed_vcfs  = Channel.empty()
    ch_versions          = Channel.empty()

    //
    // Create indexes for all the GVCF files
    //

    TABIX_GVCFS(
        gvcfs
    )

    indexed_gvcfs = gvcfs
                    .combine(TABIX_GVCFS.out.tbi, by: 0)
                    .map({ meta, gvcf, tbi -> 
                        [ meta, gvcf, tbi, []]
                    })

    ch_versions = ch_versions.mix(TABIX_GVCFS.out.versions)

    //
    // Reblock the single sample GVCF files
    //

    REBLOCKGVCF(
        indexed_gvcfs,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    ch_versions = ch_versions.mix(REBLOCKGVCF.out.versions)

    //
    // Combine the GVCFs in each family
    //

    combine_gvcfs_input = REBLOCKGVCF.out.vcf.map(
    { meta, gvcf, tbi ->
        def new_meta = [:]
        new_meta.id = meta.family

        [ new_meta, gvcf, tbi ]
    })
    .groupTuple()

    COMBINEGVCFS(
        combine_gvcfs_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COMBINEGVCFS.out.versions)

    //
    // Create indexes for the combined GVCFs
    //

    combined_gvcfs = COMBINEGVCFS.out.combined_gvcf

    TABIX_COMBINED_GVCFS(
        combined_gvcfs
    )

    ch_versions = ch_versions.mix(TABIX_COMBINED_GVCFS.out.versions)

    indexed_combined_gvcfs = combined_gvcfs
                            .combine(TABIX_COMBINED_GVCFS.out.tbi, by:0)
                            
    //
    // Genotype the combined GVCFs
    //

    genotype_gvcfs_input = indexed_combined_gvcfs
                            .map({ meta, gvcf, tbi ->
                                [ meta, gvcf, tbi, [], [] ]
                            })

    GENOTYPE_GVCFS(
        genotype_gvcfs_input,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    BGZIP_GENOTYPED_VCFS(
        GENOTYPE_GVCFS.out.vcf
    )

    ch_versions = ch_versions.mix(BGZIP_GENOTYPED_VCFS.out.versions)

    //
    // Add pedigree information
    //

    PEDFILTER(
        peds
    )

    ch_versions = ch_versions.mix(PEDFILTER.out.versions)


    merge_vcf_headers_input = BGZIP_GENOTYPED_VCFS.out.output
                                .combine(PEDFILTER.out.vcf, by:0)

    MERGE_VCF_HEADERS(
        merge_vcf_headers_input
    )

    BGZIP_PED_VCFS( 
        MERGE_VCF_HEADERS.out.vcf
    )

    ch_versions = ch_versions.mix(MERGE_VCF_HEADERS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_PED_VCFS.out.versions)

    //
    // Filter the variants 
    //

    FILTER_SNPS(
        BGZIP_PED_VCFS.out.output
    )

    FILTER_INDELS(
        FILTER_SNPS.out.vcf
    )

    ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)
    ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)
    
    post_processed_vcfs = FILTER_INDELS.out.vcf

    emit:
    post_processed_vcfs  
    versions = ch_versions
}