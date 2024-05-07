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
include { MERG_SORT_DUP as MERG_SORT_DUP_S } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'
include { MERG_SORT_DUP as MERG_SORT_DUP_H } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'

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

    // Prepare reference file [meta fasta] [meta, fai]
    ch_ref = Channel.fromPath(params.reference_fasta)
                            .map{ path -> [ [id: 'fasta'], path ] }
                            .mix( Channel.fromPath(params.reference_fai)
                            .map{ path -> [ [id: 'fai'], path ] } )

    // HISAT2_ALIGN each read group
    if (params.tools.split(',').contains('hisat2_aln')){

        index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

        HISAT2_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            splicesites
        )

        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

        // HISAT2_ALIGN.out.bam.subscribe { println("Bam Output: ${it}") }

        MERG_SORT_DUP_H( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            HISAT2_ALIGN.out.bam,
            ch_ref
        )

        ch_versions = ch_versions.mix(MERG_SORT_DUP_H.out.versions)

    }

    // STAR_ALIGN each read group
    if (params.tools.split(',').contains('star_aln')){

        index = Channel.fromPath(params.star_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        gtf = Channel.fromPath(params.reference_gtf).collect()
                    .map { path -> [ [id: 'gtf'], path ] }

        STAR_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            gtf
        )

        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        // STAR_ALIGN.out.bam.subscribe { println("Star bam Input: ${it}") }

        //Merge read groups into one file follow by sort,indexing,optinal markDup and CRAM conversion

        // ch_ref = Channel.fromPath(params.reference_fasta)
        //                     .map{ path -> [ [id: 'fasta'], path ] }
        //                     .mix( Channel.fromPath(params.reference_fai)
        //                     .map{ path -> [ [id: 'fai'], path ] } )

        // ch_ref.subscribe { println("Reference Input: ${it}") }

        // Create an input channel with the metadata and the reference files path


        MERG_SORT_DUP_S( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            STAR_ALIGN.out.bam,
            ch_ref
        )

        ch_versions = ch_versions.mix(MERG_SORT_DUP_S.out.versions)

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
