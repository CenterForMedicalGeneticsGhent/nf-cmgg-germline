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
    crams                                     // channel: [ val(meta), cram, crai, bed, ped ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, cram, crai, bed, ped ]
def create_cram_channel(LinkedHashMap row) {

    // Read the PED file
    def ped = file(row.ped).text

    // Perform a validity check on the PED file since vcf2db is picky and not capable of giving good error messages
    comment_count = 0
    line_count = 0

    for( line : ped.readLines()) {
        line_count++
        if (line_count == 1 && line ==~ /^#.*$/) {
            continue
        }
        else if (line_count > 1 && line ==~ /^#.*$/) {
            exit 1, "A commented line was found on line $line_count in $row.ped, the only commented line allowed is an optional header on line 1."
        }
        else if (line_count == 1 && line ==~ /^#.* $/) {
            exit 1, "The header in $row.ped contains a trailing space, please remove this."
        }
        else if (line ==~ /^.+#.*$/) {
            exit 1, "A '#' has been found as a non-starting character on line $line_count in $row.ped, this is an illegal character and should be removed."
        }
        else if (line ==~ /^[^#].* .*$/) {
            exit 1, "A space has been found on line $line_count in $row.ped, please only use tabs to seperate the values (and change spaces in names to '_')."
        }
        else if ((line ==~ /^(\w+\t)+\w+$/) == false) {
            exit 1, "An illegal character has been found on line $line_count in $row.ped, only a-z; A-Z; 0-9 and '_' are allowed as column values."
        }
        else if ((line ==~ /^(\w+\t){5}\w+$/) == false) {
            exit 1, "$row.ped should contain exactly 6 tab-delimited columns (family_id    individual_id    paternal_id    maternal_id    sex    phenotype). This is not the case on line ${line_count}."
        }
    }
    if (ped =~ /\n$/) {
        exit 1, "An empty new line has been detected at the end of $row.ped, please remove this line."
    }

    // get family_id
    def family_id = (ped =~ /\n([^#]\w+)/)[0][1]

    // create meta map
    def meta = [:]
    meta.id = row.sample
    meta.family = family_id

    // add path(s) of the file(s) to the channel with meta
    def cram_meta = []
    if (!file(row.cram).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> CRAM file does not exist!\n${row.cram}"
    }
    if (row.crai != ""){
        if (!file(row.crai).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> CRAI file does not exist!\n${row.crai}"
        }
    }
    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BED file does not exist!\n${row.bed}"
    }
    if (!file(row.ped).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> PED file does not exist!\n${row.ped}"
    }

    def crai = row.crai == "" ? [] : file(row.crai)

    cram_meta = [ meta, file(row.cram), crai, file(row.bed), file(row.ped) ]

    return cram_meta
}
