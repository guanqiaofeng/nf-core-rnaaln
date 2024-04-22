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

    STAGE_INPUT.out.meta_files.subscribe { println("Meta Files Output: ${it}") }

    // index = Channel.fromPath(params.hisat2_index)
    //             .map { path -> [ [id: 'index'], path ] }

    // splicesites = Channel.fromPath(params.reference_splicesites)
    //                 .map { path -> [ [id: 'splicesites'], path ] }

    // ch_fasta = Channel.fromPath(params.reference_fasta)
    //             .map { path -> [ [id: 'fasta'], path ] }



    // // Create a channel that emits a list with one element (the fasta file path)
    // // ch_fasta = Channel.value([file(params.reference_fasta)])

    // // println("HISAT2 index path: ${params.hisat2_index}")
    // // println("splicesites path: ${params.reference_splicesites}")

    // index.subscribe { println("Index Files Output: ${it}") }
    // splicesites.subscribe { println("Splice Files Output: ${it}") }
    // ch_fasta.subscribe { println("Fasta Files Output: ${it}") }

    // Resolve the actual index files
    // indexFiles = Channel.value(params.hisat2_index)


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
