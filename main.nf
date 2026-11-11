#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/ezflux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/ezflux
    Website: https://nf-co.re/ezflux
    Slack  : https://nfcore.slack.com/channels/ezflux
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EZFLUX  } from './workflows/ezflux'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_ezflux_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_ezflux_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_EZFLUX {

    take:
    reads       // channel: [ id, single_end, n_runs, fastqs ] per sample (--input)
    asv_samples // channel: [ id ] per sample (--seqtab_input)
    seqtab      // channel: the ASV table (--seqtab_input)
    diets       // channel: [ id, medium, diet ] per sample and diet (--media)
    growth      // channel: [ sample, growth_zip, build_manifest ] per growth result (--growth_input)

    main:

    //
    // WORKFLOW: Run pipeline
    //
    EZFLUX (
        reads,
        asv_samples,
        seqtab,
        diets,
        growth,
        params.growth_input,
        params.micom_db,
        params.sylph_db,
        params.model_metadata,
        params.blast_db,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = EZFLUX.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.seqtab_input,
        params.growth_input,
        params.media,
        params.medium,
        params.sylph_db,
        params.micom_db,
        params.model_metadata,
        params.blast_db,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_EZFLUX (
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.asv_samples,
        PIPELINE_INITIALISATION.out.seqtab,
        PIPELINE_INITIALISATION.out.diets,
        PIPELINE_INITIALISATION.out.growth
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_EZFLUX.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
