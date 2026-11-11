/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPROCESSING          } from '../subworkflows/local/preprocessing/main'
include { GROW                   } from '../subworkflows/local/grow/main'
include { MERGE                  } from '../subworkflows/local/merge/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ezflux_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EZFLUX {

    take:
    ch_reads       // channel: [ id, single_end, n_runs, fastqs ] per sample with reads
    ch_asv_samples // channel: [ id ] per sample of the ASV table
    ch_seqtab      // channel: the ASV table
    ch_diets       // channel: [ id, medium, diet ] per sample and diet
    ch_growth      // channel: [ sample, growth_zip, build_manifest ] per growth result (--growth_input)
    growth_input   // string: path to the growth results samplesheet (merge only)
    micom_db       // string: path to the micom model database
    sylph_db       // string: path to the sylph database (reads only)
    model_metadata // string: path to the model metadata (taxonomy -> model files)
    blast_db       // string: path to a BLAST 16S database directory (ASV tables only)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // PHASE 1 - PREPROCESSING: reads (sylph) or ASV table (blastn) -> micom taxonomy tables
    // PHASE 2 - GROWTH: every sample on every diet that applies to it
    // (both skipped with --growth_input, which brings the growth results of earlier runs)
    //
    def ch_zips      = ch_growth.map { row -> row.growth_zip }
    def ch_manifests = ch_growth.filter { row -> row.build_manifest }.map { row -> row.build_manifest }
    if (!growth_input) {
        def ch_micom_db       = channel.value(file(micom_db, checkIfExists: true))
        def ch_model_metadata = channel.value(file(model_metadata, checkIfExists: true))
        // optional databases: an empty channel when not given (channel.empty() breaks .map in the typed subworkflow)
        def ch_sylph_db       = sylph_db ? channel.value(file(sylph_db, checkIfExists: true)) : channel.fromList([])
        def ch_blast_db       = blast_db ? channel.value(file(blast_db, checkIfExists: true)) : channel.fromList([])

        PREPROCESSING(ch_reads, ch_asv_samples, ch_seqtab, ch_sylph_db, ch_blast_db, ch_model_metadata, ch_micom_db)
        GROW(PREPROCESSING.out.taxonomy, ch_diets, ch_micom_db)
        ch_zips      = GROW.out.growth.map { g -> g.zip }
        ch_manifests = GROW.out.growth.map { g -> g.build_manifest }
    }

    //
    // PHASE 3 - MERGING: combined results, exchanges, growth rates and plots
    //
    MERGE(ch_zips, ch_manifests)

    // the static plots (growth rates, found fraction, niches) go into the MultiQC report
    ch_multiqc_files = ch_multiqc_files.mix(MERGE.out.plots)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'ezflux_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'ezflux'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
