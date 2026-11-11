//
// Subworkflow with functionality specific to the nf-core/ezflux pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet (reads)
    seqtab_input      //  string: Path to a wide ASV sequence table (alternative to the samplesheet)
    growth_input      //  string: Path to a samplesheet of growth results of earlier runs (merge only)
    media             //  string: Path to media samplesheet
    medium            //  string: Default medium used when no media samplesheet is given
    sylph_db          //  string: Path to the sylph database (only needed for reads)
    micom_db          //  string: Path to the micom model database (not needed with --growth_input)
    model_metadata    //  string: Path to the model metadata (taxonomy -> model file)
    blast_db          //  string: Path to a BLAST 16S database directory (only needed with --seqtab_input)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/ezflux ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/ezflux/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --media media.csv --micom_db agora.qza --model_metadata agora.metadata.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Read and cross-validate the inputs: EXACTLY ONE of reads (--input), an ASV table (--seqtab_input)
    // or the growth results of earlier runs (--growth_input, merging only)
    //
    if ([ input, seqtab_input, growth_input ].count { path -> path } != 1) {
        error("Use exactly one of --input (reads), --seqtab_input (ASV table) or --growth_input (growth results of earlier runs)")
    }
    if (!growth_input) {
        if (!micom_db) {
            error("--micom_db is required with --input or --seqtab_input")
        }
        if (!model_metadata) {
            error("--model_metadata is required (the metadata file shipped with the model database)")
        }
    }
    if (input && !sylph_db) {
        error("--input (reads) needs --sylph_db")
    }
    if (seqtab_input && !blast_db) {
        error("--seqtab_input needs --blast_db (a BLAST 16S database directory)")
    }

    def growth_rows = growth_input ? growthSamples(samplesheetToList(growth_input, "${projectDir}/assets/schema_growth.json")) : []
    def reads_rows = input ? readsSamples(samplesheetToList(input, "${projectDir}/assets/schema_input.json")) : []
    def seqtab_file = seqtab_input ? file(seqtab_input, checkIfExists: true) : null
    def asv_rows = seqtab_file ? seqtabSamples(seqtab_file) : []
    def ids = (reads_rows + asv_rows).collect { sample -> sample.id }
    // diets are only needed when the growth is simulated
    def media_rows = growth_input ? [] : media
        ? samplesheetToList(media, "${projectDir}/assets/schema_media.json")
        : [ [ 'default', file(medium, checkIfExists: true), [] ] ]

    // one Map per sample (reads route / ASV route) and per sample x diet
    channel.fromList(reads_rows).set { ch_reads }
    channel.fromList(asv_rows).set { ch_asv_samples }
    channel.fromList(growth_input ? [] : expandDiets(media_rows, ids)).set { ch_diets }
    channel.fromList(seqtab_file ? [ seqtab_file ] : []).set { ch_seqtab }
    channel.fromList(growth_rows).set { ch_growth }

    emit:
    reads       = ch_reads         // channel: [ id, single_end, n_runs, fastqs ] (maps)
    asv_samples = ch_asv_samples   // channel: [ id ] (maps)
    seqtab      = ch_seqtab        // channel: path (empty unless --seqtab_input)
    diets       = ch_diets         // channel: [ id, medium, diet ] (maps), one per sample and diet
    growth      = ch_growth        // channel: [ sample, growth_zip, build_manifest ] (maps), only with --growth_input
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Growth samplesheet rows [ [ id ], growth_zip, build_manifest or [] ] -> one Map per row
//
def growthSamples(rows) {
    return rows.collect { row -> [ sample: row[0].id.toString(), growth_zip: row[1], build_manifest: row[2] ?: null ] }
}

//
// Samplesheet rows [ [ id, run_accession ], fastq_1, fastq_2 or [] ] -> one Map per sample.
// fastqs are ordered by run_accession as r1_run1, r2_run1, r1_run2, ... (the layout cat/fastq expects).
//
def readsSamples(rows) {
    return rows.groupBy { row -> row[0].id.toString() }.collect { id, runs ->
        if (runs.collect { run -> run[2] as boolean }.unique().size() > 1) {
            error("Sample '${id}' mixes single-end and paired-end runs")
        }
        def sorted = runs.sort(false) { run -> (run[0].run_accession ?: '').toString() }
        return [
            id: id,
            single_end: !sorted[0][2],
            n_runs: sorted.size(),
            fastqs: sorted.collectMany { run -> run[2] ? [ run[1], run[2] ] : [ run[1] ] },
        ]
    }
}

//
// The sample ids of a wide ASV table are its first column
//
def seqtabSamples(seqtab) {
    def lines = seqtab.readLines().findAll { line -> line.trim() }
    def sep = lines[0].contains('\t') ? '\t' : ','
    def ids = lines.drop(1).collect { line -> line.split(sep, 2)[0].trim().replaceAll(/^"|"$/, '') }
    if (!ids) {
        error("The ASV table has no samples (rows): ${seqtab}")
    }
    def dups = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (dups) {
        error("The ASV table has duplicated sample ids in its first column: ${dups.join(', ')}")
    }
    if (ids.any { id -> !id || id =~ /\s/ }) {
        error("The ASV table sample ids (first column) cannot be empty or contain spaces")
    }
    return ids.collect { id -> [ id: id ] }
}

//
// Media rows [ medium_name, medium_path, sample_id or [] ] -> one Map per sample and diet.
// A row without a sample_id applies to every sample.
//
def expandDiets(media_rows, ids) {
    def unknown = media_rows.collect { row -> row[2] }.findAll { id -> id && !(id.toString() in ids) }.unique()
    if (unknown) {
        error("--media lists sample_id(s) that are not in the input: ${unknown.join(', ')}")
    }
    return ids.collectMany { id ->
        def mine = media_rows.findAll { row -> !row[2] || row[2].toString() == id }
        if (!mine) {
            error("Sample '${id}' has no diet: add a --media row with sample_id '${id}' or one without a sample_id")
        }
        def names = mine.collect { row -> row[0].toString() }
        if (names.toSet().size() != names.size()) {
            error("Sample '${id}' gets the same medium_name more than once in --media")
        }
        return mine.collect { row -> [ id: id, medium: row[0].toString(), diet: row[1] ] }
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "MICOM (Diener et al. 2020),",
            "sylph (Shaw & Yu 2024),",
            "BLAST+ (Camacho et al. 2009),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Diener, C., Gibbons, S. M., & Resendis-Antonio, O. (2020). MICOM: Metagenome-scale modeling to infer metabolic interactions in the gut microbiota. mSystems, 5(1), e00606-19. doi: 10.1128/mSystems.00606-19</li>",
            "<li>Camacho, C., Coulouris, G., Avagyan, V., et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi: 10.1186/1471-2105-10-421</li>",
            "<li>Shaw, J., & Yu, Y. W. (2024). Rapid species-level metagenome profiling and containment estimation with sylph. Nature Biotechnology. doi: 10.1038/s41587-024-02412-y</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
