/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { HISAT2_ALIGN } from '../modules/local/hisat2/align/main'
include { STAR_ALIGN } from '../modules/local/star/align/main'
include { MERGE_DUP as MERG_DUP_S } from '../subworkflows/icgc-argo-workflows/merge_dup/main'
include { MERGE_DUP as MERG_DUP_ST } from '../subworkflows/icgc-argo-workflows/merge_dup/main'
include { MERGE_DUP as MERG_DUP_H } from '../subworkflows/icgc-argo-workflows/merge_dup/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_S } from '../modules/icgc-argo-workflows/payload/alignment/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_ST } from '../modules/icgc-argo-workflows/payload/alignment/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_H } from '../modules/icgc-argo-workflows/payload/alignment/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_S } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_ST } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_H } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { MERGE_SPLICE_JUNCTION as MERGE_SPLICE_JUNCTION_S } from '../modules/local/merge_splice_junction/main.nf'
include { MERGE_SPLICE_JUNCTION as MERGE_SPLICE_JUNCTION_H } from '../modules/local/merge_splice_junction/main.nf'
include { PAYLOAD_SPLICE_JUNCTION as PAYLOAD_SPLICE_JUNCTION_S } from '../modules/icgc-argo-workflows/payload/splicejunction/main'
include { PAYLOAD_SPLICE_JUNCTION as PAYLOAD_SPLICE_JUNCTION_H } from '../modules/icgc-argo-workflows/payload/splicejunction/main'
include { SONG_SCORE_UPLOAD as UPLOAD_NOVEL_SPLICE_S } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_NOVEL_SPLICE_H } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { PICARD_COLLECTRNASEQMETRICS as PICARD_COLLECTRNASEQMETRICS_S } from '../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTRNASEQMETRICS as PICARD_COLLECTRNASEQMETRICS_H } from '../modules/nf-core/picard/collectrnaseqmetrics/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_S } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_H } from '../modules/nf-core/samtools/stats/main'
include { MULTIQC as MULTIQC_S } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_H } from '../modules/nf-core/multiqc/main'
include { PREP_METRICS as PREP_METRICS_S } from '../modules/icgc-argo-workflows/prep/metrics/main'
include { PREP_METRICS as PREP_METRICS_H } from '../modules/icgc-argo-workflows/prep/metrics/main'
include { PAYLOAD_QCMETRICS as  PAYLOAD_QCMETRICS_S} from '../modules/icgc-argo-workflows/payload/qcmetrics/main'
include { PAYLOAD_QCMETRICS as  PAYLOAD_QCMETRICS_H} from '../modules/icgc-argo-workflows/payload/qcmetrics/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_S } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_H } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { PREP_REF_TRANS } from '../subworkflows/local/gen_transcript_ref'
include { CLEANUP as CLEAN_ALN_H} from '../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_ALN_S} from '../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_ALN_ST} from '../modules/icgc-argo-workflows/cleanup/main'


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

    // Prepare reference file [meta fasta] [meta, fai]
    ch_ref = Channel.fromPath(params.reference_fasta)
                            .map{ path -> [ [id: 'fasta'], path ] }
                            .mix( Channel.fromPath(params.reference_fai)
                            .map{ path -> [ [id: 'fai'], path ] } )

    // HISAT2 //
    if (params.tools.split(',').contains('hisat2_aln')){

        // HISAT2 - ALIGN //
        index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

        // ALN in read group level
        HISAT2_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            splicesites
        )
        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

        // MERG in sample level
        MERG_DUP_H( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            HISAT2_ALIGN.out.bam,
            ch_ref
        )
        ch_versions = ch_versions.mix(MERG_DUP_H.out.versions)

        // Combine channels to determine upload status and payload creation
        MERG_DUP_H.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta,cram,crai,upRdpc,metaB,analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
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

        // Make ALN payload
        PAYLOAD_ALIGNMENT_H(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_h_aln_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(HISAT2_ALIGN.out.versions)
            .mix(MERG_DUP_H.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_H.out.versions)

        // Upload files - aligment
        UPLOAD_ALIGNMENT_H(PAYLOAD_ALIGNMENT_H.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_H.out.versions)

        // HISAT2 - SPLICE JUNCTION //
        // Collect Splice Junctions
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
                data_type:"${info.data_type.collect()}",  // later check whether data type is correct **
                size:"${info.size.collect()}",
                tool: "${meta.tool}",
                aln: "hisat2"
                ],txt.collect()
            ]
        }.set{ch_h_txts}

        // MERG splice junctions in sample level
        MERGE_SPLICE_JUNCTION_H(ch_h_txts)
        ch_versions = ch_versions.mix(MERGE_SPLICE_JUNCTION_H.out.versions)

        // Combine channels to determine upload status and payload creation
        MERGE_SPLICE_JUNCTION_H.out.all_novel_splice
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta,txt,upRdpc,metaB,analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
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

        // Make payload - splice junctions
        PAYLOAD_SPLICE_JUNCTION_H(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_h_novel_splice_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(HISAT2_ALIGN.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_SPLICE_JUNCTION_H.out.versions)

        // Upload files - aligment
        UPLOAD_NOVEL_SPLICE_H(PAYLOAD_SPLICE_JUNCTION_H.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_NOVEL_SPLICE_H.out.versions)

        // QC Matrics \\
        // Samtools stats
        SAMTOOLS_STATS_H(
            MERG_DUP_H.out.cram_alignment_index,
            Channel.fromPath(params.reference_fasta).map{ it -> [ [ id:'fasta' ], it ] }
        )
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_H.out.versions)

        // Picard
        PICARD_COLLECTRNASEQMETRICS_H(
            MERG_DUP_H.out.cram_alignment_index,
            Channel.fromPath(params.ref_flat),
            Channel.fromPath(params.reference_fasta),
            Channel.fromPath(params.reference_fai),
            Channel.fromPath(params.rrna_intervals)
        )

        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS_H.out.versions)

        // MultiQC
        ch_reports = (
            Channel.empty()
            .mix(PICARD_COLLECTRNASEQMETRICS_H.out.metrics)
            .mix(HISAT2_ALIGN.out.summary)
            .mix(HISAT2_ALIGN.out.metrix)
            .mix(SAMTOOLS_STATS_H.out.stats)
        )

        ch_multiqc = Channel.empty()
        ch_multiqc = ch_multiqc.mix(ch_reports.collect{meta, report -> report}).ifEmpty([])

        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
        ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

        MULTIQC_H (
        ch_multiqc.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
        )
        ch_versions = ch_versions.mix(MULTIQC_H.out.versions)

        // metrics preparation
        PICARD_COLLECTRNASEQMETRICS_H.out.metrics
        .combine(MULTIQC_H.out.data)
        .map{
            meta, path, multiqc_dic ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    experiment:"${meta.experiment}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    read_groups_count: "${meta.read_groups_count}"
                ],multiqc_dic
            ]
        }
        .set{ch_h_prep_metrics}

        PREP_METRICS_H(
            ch_h_prep_metrics,
            []
        )

        PICARD_COLLECTRNASEQMETRICS_H.out.metrics
        .combine(MULTIQC_H.out.picard_multi)
        .combine(MULTIQC_H.out.hisat2_multi)
        .combine(MULTIQC_H.out.samtools_multi)
        .map{
            meta, path, picard, hisat2, samtools ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    experiment:"${meta.experiment}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    read_groups_count: "${meta.read_groups_count}"
                ],[picard, hisat2, samtools]
            ]
        }
        .set{ch_h_prep_metrics_files}


        // Payload generation - qc metrics
        // TAR.out.stats
        ch_h_prep_metrics_files
        .combine(PREP_METRICS_H.out.metrics_json)
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta, qcfiles_to_upload, metaB, multiqc, upRdpc, metaC, meta_analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
                    read_groups_count: "${meta.read_groups_count}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ], qcfiles_to_upload, meta_analysis, multiqc
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_h_qcmetrics_payload}

        Channel.empty()
        .mix(SAMTOOLS_STATS_H.out.stats.map{meta,files -> files}.collect())
        .mix(PICARD_COLLECTRNASEQMETRICS_H.out.metrics.map{meta,file -> file}.collect())
        .mix(MULTIQC_H.out.data.collect())
        .mix(PREP_METRICS_H.out.metrics_json.map{meta,file -> file}.collect())
        .collect()
        .set{hisat2_qc_cleanup}

        PAYLOAD_QCMETRICS_H( // [val(meta) path(json), [path(picard_multiQC), path(hisat2_multiQC)], path(multiQC)]
            ch_h_qcmetrics_payload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(HISAT2_ALIGN.out.versions)
            .mix(MERG_DUP_H.out.versions)
            .mix(PICARD_COLLECTRNASEQMETRICS_H.out.versions)
            .mix(MULTIQC_H.out.versions)
            .collectFile(name: 'collated_versions.yml')
            )
        ch_versions = ch_versions.mix(PAYLOAD_QCMETRICS_H.out.versions)

        // upload - qc metrics
        UPLOAD_QC_H(PAYLOAD_QCMETRICS_H.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_QC_H.out.versions)

        hisat2OutFlag_ch = UPLOAD_ALIGNMENT_H.out.analysis_id.concat(UPLOAD_NOVEL_SPLICE_H.out.analysis_id, UPLOAD_QC_H.out.analysis_id).collect()
    }

    // STAR //
    if (params.tools.split(',').contains('star_aln')){

        // STAR - alignment //
        index = Channel.fromPath(params.star_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        gtf = Channel.fromPath(params.reference_gtf).collect()
                    .map { path -> [ [id: 'gtf'], path ] }

        // ALN in readgroup level
        STAR_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            gtf
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        // MERG in sample level - alignment genome
        MERG_DUP_S( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            STAR_ALIGN.out.bam,
            ch_ref
        )
        ch_versions = ch_versions.mix(MERG_DUP_S.out.versions)

        // Combine channels to determine upload status and payload creation
        MERG_DUP_S.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta,cram,crai,upRdpc,metaB,analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
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

        // Make payload - alignment
        PAYLOAD_ALIGNMENT_S(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_s_aln_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .mix(MERG_DUP_S.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_S.out.versions)

        // Upload files - alignment
        UPLOAD_ALIGNMENT_S(PAYLOAD_ALIGNMENT_S.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_S.out.versions)

        // Transcriptome level alignment
        // prepare transcriptome reference
        if (params.reference_trans_fasta && params.reference_trans_fai) {
            // if transcript fasta and fai provided, use them
            ch_ref_trans = Channel.fromPath(params.reference_trans_fasta)
                            .map{ path -> [ [id: 'fasta'], path ] }
                            .mix( Channel.fromPath(params.reference_trans_fai)
                            .map{ path -> [ [id: 'fai'], path ] } )
        } else {
            // Prepare transcript fasta and fai
            PREP_REF_TRANS(
                Channel.fromPath(params.reference_fasta), // path(fasta)
                gtf // [meta, path(fasta)]
            )
            ch_ref_trans = PREP_REF_TRANS.out.trans_ref
        }

        // MERG in sample level - alignment transcriptome
        MERG_DUP_ST( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            STAR_ALIGN.out.bam_transcript,
            ch_ref_trans
        )
        ch_versions = ch_versions.mix(MERG_DUP_ST.out.versions)

        // Combine channels to determine upload status and payload creation - alignment transcriptome
        MERG_DUP_ST.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta,cram,crai,upRdpc,metaB,analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],[cram,crai],analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_st_aln_payload}

        // Make payload - alignment
        PAYLOAD_ALIGNMENT_ST(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_st_aln_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .mix(MERG_DUP_ST.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_ST.out.versions)

        // Upload files - alignment
        UPLOAD_ALIGNMENT_ST(PAYLOAD_ALIGNMENT_ST.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_ST.out.versions)

        // Collect Splice Junctions
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

        // Merge novel splice sites
        MERGE_SPLICE_JUNCTION_S(ch_s_txts)
        ch_versions = ch_versions.mix(MERGE_SPLICE_JUNCTION_S.out.versions)

        // Combine channels to determine upload status and payload creation
        MERGE_SPLICE_JUNCTION_S.out.all_novel_splice
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta,txt,upRdpc,metaB,analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
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

        // Make payload - splice junction
        PAYLOAD_SPLICE_JUNCTION_S(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_s_novel_splice_payload.upload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_S.out.versions)

        // Upload files - aligment
        UPLOAD_NOVEL_SPLICE_S(PAYLOAD_SPLICE_JUNCTION_S.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_NOVEL_SPLICE_S.out.versions)

        // QC Matrics \\
        // Samtools stats
        SAMTOOLS_STATS_S(
            MERG_DUP_S.out.cram_alignment_index,
            Channel.fromPath(params.reference_fasta).map{ it -> [ [ id:'fasta' ], it ] }
        )
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_S.out.versions)

        PICARD_COLLECTRNASEQMETRICS_S(
            MERG_DUP_S.out.cram_alignment_index,
            Channel.fromPath(params.ref_flat),
            Channel.fromPath(params.reference_fasta),
            Channel.fromPath(params.reference_fai),
            Channel.fromPath(params.rrna_intervals)
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS_S.out.versions)

        // MultiQC
        ch_reports_s = (
            Channel.empty()
            .mix(PICARD_COLLECTRNASEQMETRICS_S.out.metrics)
            .mix(STAR_ALIGN.out.log_final)
            .mix(SAMTOOLS_STATS_S.out.stats)
        )

        // Check if STAR_ALIGN.out.read_per_gene_tab exists
        if (STAR_ALIGN.out.read_per_gene_tab) {
            ch_reports_s = ch_reports_s.mix(STAR_ALIGN.out.read_per_gene_tab)
        }

        ch_multiqc_s = Channel.empty()
        ch_multiqc_s = ch_multiqc_s.mix(ch_reports_s.collect{meta, report -> report}).ifEmpty([])

        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
        ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

        MULTIQC_S (
        ch_multiqc_s.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
        )
        ch_versions = ch_versions.mix(MULTIQC_S.out.versions)

        // Prepare Metrics
        PICARD_COLLECTRNASEQMETRICS_S.out.metrics
        .combine(MULTIQC_S.out.data)
        .map{
            meta, path, multiqc_dic ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    experiment:"${meta.experiment}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    read_groups_count: "${meta.numLanes}"
                ],multiqc_dic
            ]
        }
        .set{ch_s_prep_metrics}

        PREP_METRICS_S(
            ch_s_prep_metrics,
            []
        )

        PICARD_COLLECTRNASEQMETRICS_S.out.metrics
        .combine(MULTIQC_S.out.picard_multi)
        .combine(MULTIQC_S.out.star_multi)
        .combine(MULTIQC_S.out.samtools_multi)
        .map{
            meta, path, picard, star, samtools ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    experiment:"${meta.experiment}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    read_groups_count: "${meta.read_groups_count}"
                ],[picard, star, samtools]
            ]
        }
        .set{ch_s_prep_metrics_files}

        // Payload generation - qc metrics
        // TAR.out.stats
        ch_s_prep_metrics_files
        .combine(PREP_METRICS_S.out.metrics_json)
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .map{
            meta, qcfiles_to_upload, metaB, multiqc, upRdpc, metaC, meta_analysis ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genome_build: "${params.genome_build}",
                    genome_annotation: "${params.genome_annotation}",
                    read_groups_count: "${meta.read_groups_count}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],qcfiles_to_upload, meta_analysis, multiqc
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_s_qcmetrics_payload}

        Channel.empty()
        .mix(SAMTOOLS_STATS_S.out.stats.map{meta,files -> files}.collect())
        .mix(PICARD_COLLECTRNASEQMETRICS_S.out.metrics.map{meta,file -> file}.collect())
        .mix(MULTIQC_S.out.data.collect())
        .mix(PREP_METRICS_S.out.metrics_json.map{meta,file -> file}.collect())
        .collect()
        .set{star_qc_cleanup}

        PAYLOAD_QCMETRICS_S( // [val(meta) path(json), [path(picard_multiQC), path(hisat2_multiQC)], path(multiQC)]
            ch_s_qcmetrics_payload,
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(STAR_ALIGN.out.versions)
            .mix(MERG_DUP_S.out.versions)
            .mix(PICARD_COLLECTRNASEQMETRICS_S.out.versions)
            .mix(MULTIQC_S.out.versions)
            .collectFile(name: 'collated_versions.yml')
            )
        ch_versions = ch_versions.mix(PAYLOAD_QCMETRICS_S.out.versions)

        // upload - qc metrics
        UPLOAD_QC_S(PAYLOAD_QCMETRICS_S.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_QC_S.out.versions)
    }

    if (params.tools.split(',').contains('cleanup')){
        if (params.samplesheet) {
            ch_cleanup_H = Channel.empty()
            ch_cleanup_S = Channel.empty()
        } else {
            if (params.tools.split(',').contains('hisat2_aln') && params.tools.split(',').contains('star_aln')){
                ch_cleanup_H=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
                ch_cleanup_S=Channel.empty()
            } else if (params.tools.split(',').contains('hisat2_aln')) {
                ch_cleanup_H=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
            } else if (params.tools.split(',').contains('star_aln')) {
                ch_cleanup_S=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
            }
        }
        if ( params.tools.split(',').contains('star_aln') ){

            merge_dup_ST = Channel.empty()
            .mix(MERG_DUP_ST.out.tmp_files.collect())
            .mix(MERG_DUP_ST.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())

            ch_cleanup_S=ch_cleanup_S
            .mix(MERG_DUP_S.out.tmp_files.collect())
            .mix(MERG_DUP_S.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
            .mix(MERGE_SPLICE_JUNCTION_S.out.all_novel_splice.map{meta,file -> file}.collect())
            .mix(star_qc_cleanup.collect())
            if (params.api_token){
                ch_cleanup_S=ch_cleanup_S
                .mix(PAYLOAD_ALIGNMENT_S.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(PAYLOAD_ALIGNMENT_ST.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(PAYLOAD_SPLICE_JUNCTION_S.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(PAYLOAD_QCMETRICS_S.out.payload_files.map{meta,analysis,files -> files}.collect())

                CLEAN_ALN_S(
                    ch_cleanup_S.unique().collect(),
                    UPLOAD_QC_S.out.analysis_id
                )
                CLEAN_ALN_ST(
                    merge_dup_ST.unique().collect(),
                    UPLOAD_ALIGNMENT_ST.out.analysis_id
                )
            } else {
                CLEAN_ALN_S(
                    ch_cleanup_S.unique().collect(),
                    PREP_METRICS_S.out.metrics_json
                )
            }
        }
        if ( params.tools.split(',').contains('hisat2_aln') ){
            ch_cleanup_H=ch_cleanup_H
            .mix(MERG_DUP_H.out.tmp_files.collect())
            .mix(MERG_DUP_H.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
            .mix(MERGE_SPLICE_JUNCTION_H.out.all_novel_splice.map{meta,file -> file}.collect())
            .mix(hisat2_qc_cleanup.collect())
            if (params.api_token){
                ch_cleanup_H=ch_cleanup_H
                .mix(PAYLOAD_ALIGNMENT_H.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(PAYLOAD_SPLICE_JUNCTION_H.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(PAYLOAD_QCMETRICS_H.out.payload_files.map{meta,analysis,files -> files}.collect())

                CLEAN_ALN_H(
                    ch_cleanup_H.unique().collect(),
                    hisat2OutFlag_ch
                )
            } else {
                CLEAN_ALN_H(
                    ch_cleanup_H.unique().collect(),
                    PREP_METRICS_H.out.metrics_json
                )
            }
        }
    }

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
