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

    // Perform a validity check on the PED file since vcf2db is not capable of giving good error messages
    comment_count = 0
    line_count = 0
    for( line : ped.readLines()) {
        line_count++
        if (line =~ /^#.*$/) {
            comment_count++
            if(comment_count > 1){
                exit 1, "$row.ped contains more than 1 commented line (line starting with a '#'). Only the header is allowed as a commented line."
            }
            continue
        }
        if (line =~ /^.* .*$/) {
            exit 1, "A space has been found on line $line_count in $row.ped, please only use tabs to seperate the values (and change spaces in names to '_')."
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
    if (!file(row.ped).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> PED file does not exist!\n${row.ped}"
    }

    cram_meta = [ meta, file(row.cram), file(row.crai), file(row.bed), file(row.ped) ]

    return cram_meta
}
