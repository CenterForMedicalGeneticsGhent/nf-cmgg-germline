/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTva.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { GERMLINE_VARIANT_CALLING } from '../subworkflows/local/germline_variant_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { BEDTOOLS_SPLIT              } from '../modules/nf-core/modules/bedtools/split/main'
include { BCFTOOLS_CONCAT             } from '../modules/nf-core/modules/bcftools/concat/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TVA {

    ch_versions = Channel.empty()

    //
    // Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Split the BED files into multiple subsets
    //

    beds = INPUT_CHECK.out.crams.map(
    {meta, cram, crai, bed ->
        [meta, bed]
    })

    if(params.scatter_count > 1){

        BEDTOOLS_SPLIT(
            beds,
            params.scatter_count
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

        interval_beds = BEDTOOLS_SPLIT.out.beds
        .transpose()
        .map({ meta, bed ->
            [ meta, bed, params.scatter_count ]
        })
    }
    else{
        interval_beds = beds
        .map({ meta, bed ->
            [ meta, bed, params.scatter_count ]
        })
    }

    //
    // Perform the variant calling
    //

    germline_variant_calling_input_cram = INPUT_CHECK.out.crams.map(
    {meta, cram, crai, bed ->
        [meta, cram, crai]
    })

    GERMLINE_VARIANT_CALLING(
        germline_variant_calling_input_cram,
        interval_beds,
    )

    ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)

    //
    // Merge the VCFs if split BED files were used
    //

    if (params.scatter_count > 1){

        concat_input = GERMLINE_VARIANT_CALLING.out.vcfs
                    .map({meta, vcf, tbi -> 
                        new_meta = meta.clone()
                        new_meta.id = new_meta.sample
                        [ new_meta, vcf, tbi ]
                    })
                    .groupTuple()

        BCFTOOLS_CONCAT(concat_input)

        dummy_variable = BCFTOOLS_CONCAT.out.vcf

        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    }
    else {
        dummy_variable = GERMLINE_VARIANT_CALLING.out.vcfs
                        .map({ meta, vcf, tbi ->
                            [ meta, vcf ]
                        })
    }

    dummy_variable.view()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
