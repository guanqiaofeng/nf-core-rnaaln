process NOVEL_SPLICE_MERGE {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(input) //reads may contain multiple pairs of fastq: [fastq1_1, fastq2_1, fastq1_2, fastq2_2], it can handle multiple pairs of fastq, but the read group information can not be added properly. When run with single pair of fastq, the read group info is added properly

    output:
    tuple val(meta), path("*txt"), emit: all_novel_splice
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.study_id}.${meta.patient}.${meta.sample}"

    """
    cat ${input.join(' ')} | sort -k1,1 -k2,3n -k4,4 -u > ${prefix}.${meta.aln}.novel_splicesites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | head -n1)
    END_VERSIONS
    """
}
