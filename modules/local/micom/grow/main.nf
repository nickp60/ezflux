nextflow.enable.types = true

process MICOM_GROW {
    tag "${id}__${medium}"
    label 'process_medium'

    // No container: Gurobi/CPLEX need a license that the host usually provides (module load,
    // GRB_LICENSE_FILE, token server). Pin micom to the version you validated in environment.yml.
    conda "${moduleDir}/environment.yml"

    input:
    record(id: String, medium: String, taxonomy: Path, diet: Path)   // taxonomy CSV + q2-micom medium (.qza)
    model_db: Path                                                   // micom model database, e.g. AGORA

    output:
    record(
        id: id,
        medium: medium,
        zip: file("${id}__${medium}.zip"),
        build_manifest: file("${id}__${medium}.build_manifest.csv"),
    )

    topic:
    tuple(task.process, 'micom', eval("python3 -c 'import micom; print(micom.__version__)'")) >> 'versions'

    script:
    // task.ext.args holds the micom_grow.py options (see conf/modules.config)
    def args = task.ext.args ?: ''
    """
    micom_grow.py \\
        --taxonomy ${taxonomy} \\
        --medium ${diet} \\
        --model-db ${model_db} \\
        --label ${id}__${medium} \\
        --out ${id}__${medium} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${id}__${medium}.zip
    echo 'sample_id,found_taxa,total_taxa,found_fraction,found_abundance_fraction' > ${id}__${medium}.build_manifest.csv
    """
}
