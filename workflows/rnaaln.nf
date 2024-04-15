/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
 include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
 include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
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
