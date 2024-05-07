/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { HISAT2_ALIGN } from '../modules/local/hisat2/align/main'
include { STAR_ALIGN } from '../modules/local/star/align/main'
include { MERG_SORT_DUP } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'

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
    STAGE_INPUT.out.meta_files.subscribe { println("Meta Files Output: ${it}") }

    // HISAT2_ALIGN each read group
    if (params.tools.split(',').contains('hisat2_aln')){

        index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

        // ch_fasta = Channel.fromPath(params.reference_fasta).collect()
        //         .map { path -> [ [id: 'fasta'], path ] }

        HISAT2_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            splicesites
        )

        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

        HISAT2_ALIGN.out.bam.subscribe { println("Bam Output: ${it}") }
    }

    // STAR_ALIGN each read group
    if (params.tools.split(',').contains('star_aln')){

        index = Channel.fromPath(params.star_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        gtf = Channel.fromPath(params.reference_gtf).collect()
                    .map { path -> [ [id: 'gtf'], path ] }

        // ch_fasta = Channel.fromPath(params.reference_fasta).collect()
        //         .map { path -> [ [id: 'fasta'], path ] }

        STAR_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            gtf
        )

        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    }



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
