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


    // //Perform HISAT2 alignment,make payload and upload
    // if (params.tools.split(',').contains('hisat2_aln')){
    //     //Perform Alignment per Read group
    //     BWAMEM2( //[val(meta), [path(file1),path(file2)]],[val(meta),[path(fileA),path(fileB)]]
    //         STAGE_INPUT.out.meta_files,
    //         reference_files_M2
    //     )
    //     ch_versions = ch_versions.mix(BWAMEM2.out.versions)

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
