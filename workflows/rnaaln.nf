/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { HISAT2_ALIGN } from '../modules/nf-core/hisat2/align/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAALN {

    ch_versions = Channel.empty()

    STAGE_INPUT(
        params.study_id,
        params.analysis_id,
        params.samplesheet
        )

    // STAGE_INPUT.out.meta_files.subscribe { println("Meta Files Output: ${it}") }

    // STAGE_INPUT.out.meta_files
    //     .map{ meta, fastq -> [
    //         [
    //             analysis_type:"${meta.analysis_type}",
    //             id:"${meta.study_id}.${meta.patient}.${meta.sample}",
    //             study_id:"${meta.study_id}"
    //             patient:"${meta.patient}"
    //             sex:"${meta.sex}"
    //             status:"${meta.status}"
    //             sample:"${meta.sample}"
    //             read_group:"${meta1.id}.${meta2.id}"
    //             data_type:"${meta.data_type}"
    //             numLanes:"${meta.numLanes}"
    //             experiment:"${meta.experiment}"
    //             single_end:"${meta.single_end}"
    //             sequencing_center:"${meta.sequencing_center}"
    //             library_strandedness:"${meta.library_strandedness}"
    //         ]
    //         , fastq
    //     ]

    //     }

    // STAGE_INPUT.out.meta_files
    // .map { tuple ->
    //     def meta = tuple[0]
    //     def fastq = tuple[1]
    //     // Assuming the 'id' format always contains the part before '-' as you need
    //     def laneId = meta.id.split('-')[1]
    //     return [
    //         meta: [
    //             analysis_type: meta.analysis_type,
    //             id: "${meta.study_id}.${meta.patient}.${meta.sample}",
    //             study_id: meta.study_id,
    //             patient: meta.patient,
    //             sex: meta.sex,
    //             status: meta.status,
    //             sample: meta.sample,
    //             read_group: laneId,  // Extracted part of the ID for read_group
    //             data_type: meta.data_type,
    //             numLanes: meta.numLanes,
    //             experiment: meta.experiment,
    //             single_end: meta.single_end,
    //             sequencing_center: meta.sequencing_center,
    //             library_strandedness: meta.library_strandedness
    //         ],
    //         fastq: fastq
    //     ]
    // }
    // .collect()
    // .map { allItems ->
    //     def allMetadata = allItems.collect { it.meta }
    //     def allFastqs = allItems.collect { it.fastq }
    //     def combinedLaneIds = allMetadata.collect { it.read_group }.join(',')
    //     allMetadata.each { meta ->
    //         meta.read_group = combinedLaneIds
    //     }
    //     [allMetadata, allFastqs.transpose()]
    // }
    // .subscribe { println("Reformatted Output: ${it}") }

    STAGE_INPUT.out.meta_files
    .map { tuple ->
        def meta = tuple[0]
        def fastq = tuple[1]
        def laneId = meta.id.split('-')[1]
        return [
            meta: [
                analysis_type: meta.analysis_type,
                id: "${meta.study_id}.${meta.patient}.${meta.sample}",
                study_id: meta.study_id,
                patient: meta.patient,
                sex: meta.sex,
                status: meta.status,
                sample: meta.sample,
                read_group: laneId,
                data_type: meta.data_type,
                numLanes: meta.numLanes,
                experiment: meta.experiment,
                single_end: meta.single_end,
                sequencing_center: meta.sequencing_center,
                library_strandedness: meta.library_strandedness
            ],
            fastq: fastq
        ]
    }
    .collect()
    .map { allItems ->
        def allMetadata = allItems.collect { it.meta }
        def allFastqs = allItems.collect { it.fastq }
        def combinedLaneIds = allMetadata.collect { it.read_group }.join(',')
        def unifiedMetadata = allMetadata[0]
        unifiedMetadata.read_group = combinedLaneIds
        [unifiedMetadata, allFastqs.transpose().flatten()]
    }
    .map { meta, fastqPaths ->
        tuple(meta, fastqPaths)}
    .set{ch_fastq}


    ch_fastq.subscribe { println("Reformatted Output: ${it}") }

    index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

    splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

    ch_fasta = Channel.fromPath(params.reference_fasta).collect()
                .map { path -> [ [id: 'fasta'], path ] }

    HISAT2_ALIGN(
        ch_fastq,
        index,
        splicesites
    )

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
