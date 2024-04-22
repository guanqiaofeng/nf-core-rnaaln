process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.log")                   , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def strandedness = '' // Strandedness is updated for ICGC ARGO, as only 4 entries are allowed "UNSTRANDED, FIRST_READ_SENSE_STRAND, FIRST_READ_ANTISENSE_STRAND, or empty(null)"
    if (meta.library_strandedness == 'FIRST_READ_SENSE_STRAND') {
    strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.library_strandedness == 'FIRST_READ_ANTISENSE_STRAND') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    } else {
        strandedness = ''  // Handle both UNSTRANDED and null the same way
    }

    ss = "$splicesites" ? "--known-splicesite-infile $splicesites" : ''

    //def seq_center = params.sequencing_center ? "--rg-id ${prefix} --rg SM:$prefix --rg CN:${params.sequencing_center.replaceAll('\\s','_')}" : "--rg-id ${prefix} --rg SM:$prefix"
    def seq_center = ''
    if (meta.containsKey('sequencing_center') && meta.sequencing_center) {
        // Replaces spaces with underscores in the sequencing center's name and constructs the read group string.
        seq_center = "--rg-id ${prefix} --rg SM:$prefix --rg CN:${meta.sequencing_center.replaceAll('\\s','_')}"
    } else {
        // Default read group string when no sequencing center is specified.
        seq_center = "--rg-id ${prefix} --rg SM:$prefix"
    }

    if (meta.single_end) {
        def unaligned = params.save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
        hisat2 \\
            -x \$INDEX \\
            -U $reads \\
            $strandedness \\
            $ss \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $seq_center \\
            $unaligned \\
            $args \\
            | samtools view -bS -F 4 -F 256 - > ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: $VERSION
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        def unaligned = params.save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $strandedness \\
            $ss \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $seq_center \\
            $unaligned \\
            --no-mixed \\
            --no-discordant \\
            $args \\
            | samtools view -bS -F 4 -F 8 -F 256 - \
            | samtools sort -o ${prefix}.bam

        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
        fi
        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: $VERSION
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}