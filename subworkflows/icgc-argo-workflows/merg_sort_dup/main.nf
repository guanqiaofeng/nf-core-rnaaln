// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_MERGE                        } from '../../../modules/icgc-argo-workflows/samtools/merge/main'
include { BIOBAMBAM_BAMMARKDUPLICATES2          } from '../../../modules/icgc-argo-workflows/biobambam/bammarkduplicates2/main'
include { SAMTOOLS_INDEX                        } from '../../../modules/icgc-argo-workflows/samtools/index/main'
include { SAMTOOLS_CONVERT                      } from '../../../modules/icgc-argo-workflows/samtools/convert/main'
include { TAR                                   } from '../../../modules/icgc-argo-workflows/tar/main'

workflow MERG_SORT_DUP {

    take:
    bam
    reference_files

    main:

    ch_versions = Channel.empty()

    //Categorize reference_files ([meta, .fasta|.fa] [meta, fai]) into two separate channels based on file extension (reg_org.fasta, reg_org.fai)
    // reference_files.branch{
    //     fasta : it[1].name.endsWith(".fasta") || it[1].name.endsWith(".fa")
    //     fai : it[1].name.endsWith(".fai")
    // }.set{ref_org}

    reference_files.branch{
        fasta : it[1].toString().endsWith(".fasta") || it[1].toString().endsWith(".fa")
        fai : it[1].toString().endsWith(".fai")
    }.set{ref_org}

    //Collect channel (e.g. [metaA,bamA,metaB,bamB] and seperate back in channels of [meta,bam])
    //Simplfy metadata to group and collect BAMs : [meta, [bamA,bamB,bamC]] for merging
    bam.flatten().buffer( size: 2 )
    .map{
        meta,bam ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}",
            tool: "${meta.tool}"
            ],
            [
            read_group:"${meta.id}",
            data_type:"${meta.data_type}",
            size:"${meta.size}",
            ],
            bam
        ]
    }.groupTuple(by: 0).
    map{
        meta,info,bam ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}",
            read_group:"${info.read_group.collect()}",
            data_type:"${info.data_type.collect()}",
            size:"${info.size.collect()}",
            tool: "${meta.tool}"
            ],bam.collect()
        ]
    }.set{ch_bams}

    //Merge the bam files
    SAMTOOLS_MERGE(
        ch_bams,
        ref_org.fasta,
        ref_org.fai
    )

    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Prepare channel for markdup, id updates
    if (params.tools.split(',').contains('markdup')){
        SAMTOOLS_MERGE.out.bam
        .map{
            meta,file ->
            [
                [
                    id:"${meta.id}.csort.markdup",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    date:"${meta.date}",
                    numLanes:"${meta.numLanes}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    size:"${meta.size}",
                    experiment:"${meta.experiment}",
                    tool: "${meta.tool}"
                ],
                file
            ]
        }.set{ch_markdup}
    } else {
        SAMTOOLS_MERGE.out.bam
        .map{
            meta,file ->
            [
                [
                    id:"${meta.id}.csort",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    date:"${meta.date}",
                    numLanes:"${meta.numLanes}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    size:"${meta.size}",
                    experiment:"${meta.experiment}",
                    tool: "${meta.tool}"
                ],
                file
            ]
        }.set{ch_markdup}
    }

    //If markdup specified, markdup file else return as is
    if (params.tools.split(',').contains('markdup')){
        BIOBAMBAM_BAMMARKDUPLICATES2(
            ch_markdup
        )
        ch_versions = ch_versions.mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.versions)
        BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.set{markdup_bam}
    } else {
        ch_markdup.set{markdup_bam}//meta,bam
    }

    //Index Csort.Markdup.Bam
    SAMTOOLS_INDEX(markdup_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //Prepare channel [meta, bam, bai(new index)] for conversion: use new index and Bam
    markdup_bam.combine(SAMTOOLS_INDEX.out.bai)//meta,bai
    .map{
        metaA,bam,metaB,index ->
        [
            [
                id:"${metaA.id}",
                study_id:"${metaA.study_id}",
                patient:"${metaA.patient}",
                sex:"${metaA.sex}",
                sample:"${metaA.sample}",
                numLanes:"${metaA.numLanes}",
                date:"${metaA.date}",
                read_group:"${metaA.read_group}",
                data_type:"${metaA.data_type}",
                size:"${metaA.size}",
                experiment:"${metaA.experiment}",
                tool: "${metaA.tool}"
            ],
            bam,index
        ]
    }.set{ch_convert}

    //Convert bam, bai to cram crai
    SAMTOOLS_CONVERT(
        ch_convert,
        ref_org.fasta,
        ref_org.fai
    )

    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //Prepare output channel [meta, cram, crai]
    SAMTOOLS_CONVERT.out.cram.combine(SAMTOOLS_CONVERT.out.crai)
    .map{
        metaA,cram,metaB,index ->
        [
            [
                id:"${metaA.id}",
                study_id:"${metaA.study_id}",
                patient:"${metaA.patient}",
                sex:"${metaA.sex}",
                sample:"${metaA.sample}",
                numLanes:"${metaA.numLanes}",
                date:"${metaA.date}",
                read_group:"${metaA.read_group}",
                data_type:"${metaA.data_type}",
                size:"${metaA.size}",
                experiment:"${metaA.experiment}",
                tool: "${metaA.tool}"
            ],
            cram,index
        ]
    }.set{alignment_index}

    //If Markdup specified, TAR metrics file, set channels for cleanup and metrics
    if (params.tools.split(',').contains('markdup')){
        TAR(
            BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics
            .map{ meta,file->
            [
                [
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    date:"${meta.date}",
                    numLanes:"${meta.numLanes}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    size:"${meta.size}",
                    experiment:"${meta.experiment}",
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}.aln.cram.duplicates_metrics",
                    tools:"${meta.tool}"
                ],file
            ]
            }
        )
        TAR.out.stats.set{metrics}
        Channel.empty()
        .mix(ch_bams.map{meta,files -> files}.collect())
        .mix(SAMTOOLS_MERGE.out.bam.map{meta,file -> file}.collect())
        .mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.map{meta,file -> file}.collect())
        .mix(SAMTOOLS_INDEX.out.bai.map{meta,file -> file}.collect())
        .collect()
        .set{ch_cleanup}
    } else {
        Channel.empty()
        .mix(ch_bams.map{meta,files -> files}.collect())
        .mix(SAMTOOLS_MERGE.out.bam.map{meta,file -> file}.collect())
        .mix(SAMTOOLS_INDEX.out.bai.map{meta,file -> file}.collect())
        .collect()
        .set{ch_cleanup}

        Channel.empty().set{metrics}
    }

    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}

    emit:
    cram_alignment_index = alignment_index
    tmp_files = ch_cleanup
    metrics = metrics
    versions = ch_versions                     // channel: [ versions.yml ]
}
