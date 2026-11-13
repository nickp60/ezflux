nextflow.enable.types = true

process MICOM_TAXONOMY {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    record(id: String, abundances: Path, source: String)   // source 'sylph': the raw sylph profile; 'species': a species, abundance table (BLAST route)
    model_metadata: Path                                   // Genome_file <tab> clade_name, shipped with the model database
    model_db: Path                                         // species-level micom model database (.qza / .zip)

    output:
    record(
        id: id,
        taxonomy: file("${id}.taxonomy.csv"),    // micom Batch taxonomy table
        report: file("${id}.annotation.tsv"),    // per sylph genome or species: the species it maps to and the number of models found
    )

    topic:
    tuple(task.process, 'micom', eval("python3 -c 'import micom; print(micom.__version__)'")) >> 'versions'

    script:
    """
    micom_taxonomy.py --abundances ${abundances} --source ${source} --model-metadata ${model_metadata} --model-db ${model_db} --sample ${id} --prefix ${id}
    """

    stub:
    """
    touch ${id}.taxonomy.csv ${id}.annotation.tsv
    """
}
