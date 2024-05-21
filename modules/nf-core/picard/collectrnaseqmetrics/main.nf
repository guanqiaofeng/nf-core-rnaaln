process PICARD_COLLECTRNASEQMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path ref_flat
    path fasta
    path rrna_intervals

    output:
    tuple val(meta), path("*.rna_metrics")  , emit: metrics
    tuple val(meta), path("*.pdf")          , emit: pdf, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def rrna = rrna_intervals ? "--RIBOSOMAL_INTERVALS ${rrna_intervals}" : ""
    def strandedness = meta.strandedness == 'forward' ? '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND' : meta.strandedness == 'reverse' ? '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND' : '--STRAND_SPECIFICITY NONE'
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectRnaSeqMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectRnaSeqMetrics \\
        $args \\
        $reference \\
        $rrna \\
        $strandedness \\
        --REF_FLAT $ref_flat \\
        --INPUT ${bam[0]} \\
        --OUTPUT ${prefix}.rna_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectRnaSeqMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rna_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectRnaSeqMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
