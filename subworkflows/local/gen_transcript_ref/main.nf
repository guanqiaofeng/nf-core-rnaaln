/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GFFREAD } from '../../../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREP_REF_TRANS {

    take:
    genome_fasta // path(genome_fasta)
    genome_gtf      // [meta, path(genome_gtf)]

    main:

    ch_versions = Channel.empty()

    GFFREAD(
        genome_gtf,
        genome_fasta
    )
    ch_versions = ch_versions.mix(GFFREAD.out.versions)

    SAMTOOLS_FAIDX(GFFREAD.out.gffread_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    trans_ref = GFFREAD.out.gffread_fasta.mix(SAMTOOLS_FAIDX.out.fai)

    emit:
    trans_ref = trans_ref
    // tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}
