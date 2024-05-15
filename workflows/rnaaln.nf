/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { HISAT2_ALIGN } from '../modules/local/hisat2/align/main'
include { STAR_ALIGN } from '../modules/local/star/align/main'
include { MERG_SORT_DUP as MERG_SORT_DUP_S } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'
include { MERG_SORT_DUP as MERG_SORT_DUP_H } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_S } from '../modules/local/payload/rnaseqalignment/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_H } from '../modules/local/payload/rnaseqalignment/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_S } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_H } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { NOVEL_SPLICE_MERGE as NOVEL_SPLICE_MERGE_S } from '../modules/local/novelsplice/main.nf'
include { NOVEL_SPLICE_MERGE as NOVEL_SPLICE_MERGE_H } from '../modules/local/novelsplice/main.nf'
include { PAYLOAD_NOVEL_SPLICE as PAYLOAD_NOVEL_SPLICE_S } from '../modules/local/payload/novel_splice/main'
include { PAYLOAD_NOVEL_SPLICE as PAYLOAD_NOVEL_SPLICE_H } from '../modules/local/payload/novel_splice/main'
include { SONG_SCORE_UPLOAD as UPLOAD_NOVEL_SPLICE_S } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_NOVEL_SPLICE_H } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAALN {

    ch_versions = Channel.empty()

    // Validate input, generate metadata, prepare fastq channel
    STAGE_INPUT(
        params.study_id,
        params.analysis_id,
        params.samplesheet
        )

    ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)
    // STAGE_INPUT.out.meta_files.subscribe { println("Meta Files Output: ${it}") }

    // Prepare reference file [meta fasta] [meta, fai]
    ch_ref = Channel.fromPath(params.reference_fasta)
                            .map{ path -> [ [id: 'fasta'], path ] }
                            .mix( Channel.fromPath(params.reference_fai)
                            .map{ path -> [ [id: 'fai'], path ] } )

    // HISAT2_ALIGN
    if (params.tools.split(',').contains('hisat2_aln')){

        index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

        // HISAT2_ALIGN in read group level
        HISAT2_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            splicesites
        )
        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

        // HISAT2_ALIGN.out.versions.subscribe {println("HISAT2 version output: ${it}")}

        HISAT2_ALIGN.out.bam.subscribe {println("Bam output: ${it}")}

        // MERG_SORT_DUP in sample level
        MERG_SORT_DUP_H( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            HISAT2_ALIGN.out.bam,
            ch_ref
        )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_H.out.versions)

        // MERG_SORT_DUP_H.out.cram_alignment_index.subscribe {println("Cram Output: ${it}")}

        // Combine channels to determine upload status and payload creation
        MERG_SORT_DUP_H.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            ch_ref.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],[cram,crai],analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_h_aln_payload}

        // ch_h_aln_payload.subscribe { println("Payload Prepare: ${it}") }

        // Make payload - alignment
        PAYLOAD_ALIGNMENT_H(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_h_aln_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(HISAT2_ALIGN.out.versions)
            .mix(MERG_SORT_DUP_H.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_H.out.versions)
        // PAYLOAD_ALIGNMENT_H.out.payload_files.subscribe { println("Payload Files: ${it}") }

        // Upload files - aligment
        UPLOAD_ALIGNMENT_H(PAYLOAD_ALIGNMENT_H.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_H.out.versions)
        // UPLOAD_ALIGNMENT_H.out.analysis_id.subscribe { println("Upload Analysis Id: ${it}") }

        // Collect Novel Splice sites
        HISAT2_ALIGN.out.novel_splice.flatten().buffer( size: 2 )
        .map{
            meta,txt ->
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
                txt
            ]
        }.groupTuple(by: 0)
        .map{
            meta,info,txt ->
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
                tool: "${meta.tool}",
                aln: "hisat2"
                ],txt.collect()
            ]
        }.set{ch_h_txts}

        // Merge novel splice sites
        NOVEL_SPLICE_MERGE_H(ch_h_txts)
        ch_versions = ch_versions.mix(NOVEL_SPLICE_MERGE_H.out.versions)

        NOVEL_SPLICE_MERGE_H.out.all_novel_splice.subscribe { println("Merged Novel Splice Files: ${it}") }

        // Combine channels to determine upload status and payload creation
        NOVEL_SPLICE_MERGE_H.out.all_novel_splice
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            ch_ref.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,txt,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],txt,analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_h_novel_splice_payload}

        ch_h_novel_splice_payload.subscribe { println("Novel Splice Payload Prepare: ${it}") }

        // Make payload - novel splice sites
        PAYLOAD_NOVEL_SPLICE_H(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_h_novel_splice_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(HISAT2_ALIGN.out.versions)
            .mix(MERG_SORT_DUP_H.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_H.out.versions)
        PAYLOAD_NOVEL_SPLICE_H.out.payload_files.subscribe { println("Novel Splice Payload Files: ${it}") }

        // Upload files - aligment
        UPLOAD_NOVEL_SPLICE_H(PAYLOAD_NOVEL_SPLICE_H.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_NOVEL_SPLICE_H.out.versions)
        UPLOAD_NOVEL_SPLICE_H.out.analysis_id.subscribe { println("Upload Analysis Id: ${it}") }

    }

    // STAR_ALIGN each read group
    if (params.tools.split(',').contains('star_aln')){

        index = Channel.fromPath(params.star_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        gtf = Channel.fromPath(params.reference_gtf).collect()
                    .map { path -> [ [id: 'gtf'], path ] }

        // STAR_ALIGN in read group level
        STAR_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            gtf
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        // STAR_ALIGN.out.bam.subscribe { println("bam Files: ${it}") }

        // MERG_SORT_DUP in sample level
        MERG_SORT_DUP_S( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            STAR_ALIGN.out.bam,
            ch_ref
        )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_S.out.versions)

        // MERG_SORT_DUP_S.out.cram_alignment_index.subscribe {println("Cram Output: ${it}")}

        // Combine channels to determine upload status and payload creation
        MERG_SORT_DUP_S.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            ch_ref.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],[cram,crai],analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_s_aln_payload}

        // ch_s_aln_payload.subscribe { println("Payload Prepare: ${it}") }

        // Make payload - alignment
        PAYLOAD_ALIGNMENT_S(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_s_aln_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .mix(MERG_SORT_DUP_S.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_S.out.versions)
        // PAYLOAD_ALIGNMENT_S.out.payload_files.subscribe { println("Payload Files: ${it}") }

        // Upload files - alignment
        UPLOAD_ALIGNMENT_S(PAYLOAD_ALIGNMENT_S.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_S.out.versions)
        // UPLOAD_ALIGNMENT_S.out.analysis_id.subscribe { println("Upload Analysis Id: ${it}") }

        // spl_junc_tab
        // Collect Novel Splice sites
        STAR_ALIGN.out.spl_junc_tab.flatten().buffer( size: 2 )
        .map{
            meta,txt ->
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
                txt
            ]
        }.groupTuple(by: 0)
        .map{
            meta,info,txt ->
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
                tool: "${meta.tool}",
                aln: "star"
                ],txt.collect()
            ]
        }.set{ch_s_txts}

        ch_s_txts.subscribe { println("star novel splice files: ${it}") }

        // Merge novel splice sites
        NOVEL_SPLICE_MERGE_S(ch_s_txts)
        ch_versions = ch_versions.mix(NOVEL_SPLICE_MERGE_S.out.versions)

        NOVEL_SPLICE_MERGE_S.out.all_novel_splice.subscribe { println("Merged Novel Splice Files: ${it}") }
        ch_versions.subscribe { println("Channel version after splice merge: ${it}") }

        // Combine channels to determine upload status and payload creation
        NOVEL_SPLICE_MERGE_S.out.all_novel_splice
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            ch_ref.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,txt,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],txt,analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_s_novel_splice_payload}

        ch_s_novel_splice_payload.subscribe { println("Novel Splice Payload Prepare: ${it}") }

        // Make payload - novel splice sites
        PAYLOAD_NOVEL_SPLICE_S(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_s_novel_splice_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .mix(MERG_SORT_DUP_S.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_S.out.versions)
        PAYLOAD_NOVEL_SPLICE_S.out.payload_files.subscribe { println("Novel Splice Payload Files: ${it}") }

        // Upload files - aligment
        UPLOAD_NOVEL_SPLICE_S(PAYLOAD_NOVEL_SPLICE_S.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_NOVEL_SPLICE_S.out.versions)
        UPLOAD_NOVEL_SPLICE_S.out.analysis_id.subscribe { println("Upload Analysis Id: ${it}") }


    }

    // MERGE_SORT_DUP


    //Markduplicate

    emit:
    meta_analysis = STAGE_INPUT.out.meta_analysis // channel: /path/to/multiqc_report.html
    meta_files = STAGE_INPUT.out.meta_files
    upRdpc = STAGE_INPUT.out.upRdpc

    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
