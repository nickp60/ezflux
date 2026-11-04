nextflow.enable.types = true

process MICOM_MERGE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    growth_zips: Bag<Path>        // one GrowthResults zip per sample and medium (sample_id = <sample>__<medium>)
    build_manifests: Bag<Path>    // the matching micom build manifests (found fractions)

    output:
    record(
        growth: file("*.growth.zip"),
        exchanges: file("*.exchanges.csv"),
        growth_rates: file("*.growth_rates.csv"),
        production_rates: file("*.production_rates.csv"),
        consumption_rates: file("*.consumption_rates.csv"),
        found_fraction: file("*.found_fraction.csv"),
        import_niche: file("*.import_niche.html"),
        export_niche: file("*.export_niche.html"),
        growth_rates_html: file("*.growth_rates.html"),
        plots: files("*_mqc.png"),       // static plots for MultiQC
    )

    topic:
    tuple(task.process, 'micom', eval("python3 -c 'import micom; print(micom.__version__)'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'ezflux'
    """
    micom_merge.py --prefix ${prefix} ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'ezflux'
    """
    touch ${prefix}.growth.zip ${prefix}.exchanges.csv ${prefix}.growth_rates.csv \\
        ${prefix}.production_rates.csv ${prefix}.consumption_rates.csv ${prefix}.found_fraction.csv \\
        ${prefix}.import_niche.html ${prefix}.export_niche.html ${prefix}.growth_rates.html \\
        growth_rates_mqc.png found_fraction_mqc.png import_niche_mqc.png export_niche_mqc.png
    """
}
