/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { HISAT2_ALIGN } from '../modules/local/hisat2/align/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
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

    // HISAT2_ALIGN
    if (params.tools.split(',').contains('hisat2_aln')){

        index = Channel.fromPath(params.hisat2_index).collect()
                .map { path -> [ [id: 'index'], path ] }

        splicesites = Channel.fromPath(params.reference_splicesites).collect()
                    .map { path -> [ [id: 'splicesites'], path ] }

        ch_fasta = Channel.fromPath(params.reference_fasta).collect()
                .map { path -> [ [id: 'fasta'], path ] }

        HISAT2_ALIGN(
            STAGE_INPUT.out.meta_files,
            index,
            splicesites
        )

        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

        HISAT2_ALIGN.out.bam.subscribe { println("Bam Output: ${it}") }
    }


    // Reformat fastq file channel for input to HISAT2_ALIGN
    // Original: [[meta1, [fastq1-1, fastq1-2]], [meta2, [fastq2-1, fastq2-2]]]
    // Transformed: [meta, [fastq1-1, fastq2-1, fastq1-2, fastq2-2]]
    // future work: need to handle readgroup information in the bam files in HISAT2_ALIGN

    // STAR align seperately
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
    }


    // STAR_ALIGN align together
    // if (params.tools.split(',').contains('star_aln')){

    //     STAGE_INPUT.out.meta_files
    //     .map { tuple ->
    //         def meta = tuple[0]
    //         def fastq = tuple[1]
    //         def laneId = meta.id
    //         return [
    //             meta: meta,
    //             fastq: fastq
    //         ]
    //     }
    //     .collect()
    //     .map { allItems ->
    //         def allMetadata = allItems.collect { it.meta }
    //         def allFastqs = allItems.collect { it.fastq }
    //         def combinedLaneIds = allMetadata.collect { it.read_group }.join("\\n")
    //         def combinedid = allMetadata.collect {it.id}.join(',')
    //         def unifiedMetadata = allMetadata[0]
    //         unifiedMetadata.read_group = combinedLaneIds
    //         unifiedMetadata.id = combinedid
    //         [unifiedMetadata, allFastqs.transpose().flatten()]
    //     }
    //     .map { meta, fastqPaths ->
    //         tuple(meta, fastqPaths)}
    //     .set{ch_fastq}

    //     // ch_fastq.subscribe { println("Meta Files Output: ${it}") }

    //     ch_fastq.subscribe { println("Reformatted Output: ${it}") }

    //     index = Channel.fromPath(params.star_index).collect()
    //             .map { path -> [ [id: 'index'], path ] }

    //     gtf = Channel.fromPath(params.reference_gtf).collect()
    //                 .map { path -> [ [id: 'gtf'], path ] }

    //     // ch_fasta = Channel.fromPath(params.reference_fasta).collect()
    //     //         .map { path -> [ [id: 'fasta'], path ] }

    //     STAR_ALIGN(
    //         ch_fastq,
    //         index,
    //         gtf
    //     )

    //     ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //     STAR_ALIGN.out.bam.subscribe { println("STAR Bam Output: ${it}") }
    // }

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
