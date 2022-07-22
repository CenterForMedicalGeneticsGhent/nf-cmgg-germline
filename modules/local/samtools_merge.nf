process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")
    path fasta
    path fai
    val cram_merge

    output:
    tuple val(meta), path("*.bam") , optional:true, emit: bam
    tuple val(meta), path("*.cram"), optional:true, emit: cram
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args   ?: ''
    def args2           = task.ext.args2   ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def reference       = fasta ? "--reference ${fasta}" : ""
    def convert_to_cram = cram_merge ? 
        "samtools view --threads ${task.cpus} --reference ${fasta} $args2 ${prefix}.bam -C -o ${prefix}.cram && rm ${prefix}.bam" : ""
    """
    samtools \\
        merge \\
        --threads ${task.cpus} \\
        $args \\
        ${reference} \\
        ${prefix}.bam \\
        $input_files

    $convert_to_cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
