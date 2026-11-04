nextflow.enable.types = true

process MICOM_TAXONOMY {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    record(id: String, abundances: Path)   // annotated_abundances.csv: placeholder, sample_id, species, abundance
    model_db: Path                         // species-level micom model database (.qza / .zip)

    output:
    record(
        id: id,
        taxonomy: file("${id}.taxonomy.csv"),    // micom Batch taxonomy table
        report: file("${id}.annotation.tsv"),    // per species: number of models found
    )

    topic:
    tuple(task.process, 'micom', eval("python3 -c 'import micom; print(micom.__version__)'")) >> 'versions'

    script:
    """
    micom_taxonomy.py --abundances ${abundances} --model-db ${model_db} --sample ${id} --prefix ${id}
    """

    stub:
    """
    touch ${id}.taxonomy.csv ${id}.annotation.tsv
    """
}
