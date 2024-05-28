process TAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("out/*tgz"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"

"""
mkdir out
echo ${input} | tr ' ' '\n' | xargs -I {} bash -c 'base=\$(basename {}); base_no_ext=\$(echo \${base%.*}); tar -czvf out/\${base_no_ext}.tgz {}'
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tar: \$(echo \$(tar --version | cut -f3 -d' ' ))
END_VERSIONS
"""
}
