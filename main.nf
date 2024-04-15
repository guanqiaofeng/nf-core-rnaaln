#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfcore/rnaaln
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfcore/rnaaln
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.study_id                    = WorkflowMain.getGenomeAttribute(params, 'study_id')
params.analysis_id                 = WorkflowMain.getGenomeAttribute(params, 'analysis_id')
params.samplesheet                  = WorkflowMain.getGenomeAttribute(params, 'samplesheet')

params.reference_fasta             = WorkflowMain.getGenomeAttribute(params, 'reference_fasta')
params.reference_fasta_secondary   = WorkflowMain.getGenomeAttribute(params, 'reference_fasta_secondary')

params.api_token                   = WorkflowMain.getGenomeAttribute(params, 'api_token')
params.score_url_upload            = WorkflowMain.getGenomeAttribute(params, 'score_url_upload')
params.song_url_upload             = WorkflowMain.getGenomeAttribute(params, 'song_url_upload')
params.score_url_download          = WorkflowMain.getGenomeAttribute(params, 'score_url_download')
params.song_url_download           = WorkflowMain.getGenomeAttribute(params, 'song_url_download')
params.score_url                   = WorkflowMain.getGenomeAttribute(params, 'score_url')
params.song_url                    = WorkflowMain.getGenomeAttribute(params, 'song_url')

params.tools                       = WorkflowMain.getGenomeAttribute(params, 'tools')
params.outdir                      = WorkflowMain.getGenomeAttribute(params, 'outdir')

params.local_mode                  = WorkflowMain.getGenomeAttribute(params, 'local_mode')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAALN } from './workflows/rnaaln'

//
// WORKFLOW: Run main nfcore/dnaseqaln analysis pipeline
//
workflow NFCORE_RNAALN {
    RNAALN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_RNAALN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
