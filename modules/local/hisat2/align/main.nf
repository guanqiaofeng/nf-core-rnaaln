process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads) //reads may contain multiple pairs of fastq: [fastq1_1, fastq2_1, fastq1_2, fastq2_2], it can handle multiple pairs of fastq, but the read group information can not be added properly. When run with single pair of fastq, the read group info is added properly
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)

    output:
    tuple val(meta), path("*.hisat2_Aligned.bam")    , emit: bam
    tuple val(meta), path("*_summary.txt")           , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val(meta), path("*_metrics.txt")           , emit: metrix
    tuple val(meta), path("*.novel_splicesites.txt") , emit: novel_splice
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.study_id}.${meta.patient}.${meta.sample}.${meta.id}"
    def VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def strandedness = meta.strandedness == 'forward' ? meta.single_end ? '--rna-strandness F' : '--rna-strandness FR' : meta.strandedness == 'reverse' ? meta.single_end ? '--rna-strandness R' : '--rna-strandness RF' : ''
    // Calculate the midpoint of the reads array
    // def mid = reads.size() / 2
    // def reads1 = reads[0..<mid].join(',')
    // def reads2 = reads[mid..<reads.size()].join(',')
    def reads1 = [], reads2 = [] // first half in reads1, second half in reads2
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> (ix < reads.size() / 2 ? reads1 : reads2) << v }

    """
    INDEX=\$(find -L ${index} -name "*.1.ht2" | sed 's/\\.1.ht2\$//')
    hisat2 \\
        -t -q \\
        -x \$INDEX \\
        -p $task.cpus \\
        ${meta.single_end ? "-1 ${reads.join(",")}" : "-1 ${reads1.join(",")} -2 ${reads2.join(",")}"} \\
        $strandedness \\
        --min-intronlen 20 \\
        --max-intronlen 500000 \\
        --met-file ${prefix}_metrics.txt \\
        --summary-file ${prefix}_summary.txt \\
        --new-summary \\
        --known-splicesite-infile $splicesites \\
        --novel-splicesite-outfile ${prefix}.novel_splicesites.txt \\
        -k 20 --secondary \\
        -I 50 -X 800 \\
        --rg-id $meta.id \\
        $args \\
        | samtools view -bS --no-PG - \\
        | samtools sort - \\
        | samtools reheader -P -c \'sed -e "s/^@RG.*/${meta.read_group.replaceAll(/'/,'')}/"\' - > ${prefix}.hisat2_Aligned.bam

    sort -k1,1 -k2,3n -k4,4 -u ${prefix}.novel_splicesites.txt > ${prefix}.hisat2.novel_splicesites.txt
    rm ${prefix}.novel_splicesites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
// echo "hisat2: $VERSION\nsamtools: \$(samtools --version | sed -n '1 s/^.*samtools //p')" > versions.yml
// echo "hisat2: $VERSION\nsamtools: \$(samtools --version | sed 's/^.*samtools //; s/Using.*\$//')" > versions.yml
