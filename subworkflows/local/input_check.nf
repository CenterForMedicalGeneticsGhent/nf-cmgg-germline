//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_cram_channel(it) }
        .set { crams }

    emit:
    crams                                     // channel: [ val(meta), cram, crai, bed ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, cram, crai, bed ]
def create_cram_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def cram_meta = []
    if (!file(row.cram).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> CRAM file does not exist!\n${row.cram}"
    }
    if (!file(row.crai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> CRAI file does not exist!\n${row.crai}"
    }
    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BED file does not exist!\n${row.bed}"
    }

    cram_meta = [ meta, file(row.cram), file(row.crai), file(row.bed)]

    return cram_meta
}
