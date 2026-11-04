nextflow.enable.types = true

process BLAST_FORMAT {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0'
        : 'quay.io/biocontainers/r-tidyverse:2.0.0'}"

    input:
    record(id: String, counts: Path, blast: Path)   // counts of all samples (ASV_TOFASTA) + tabular BLAST output
    model_metadata: Path

    output:
    record(
        id: id,
        abundances: file("${id}.annotated_abundances.csv"),   // placeholder,sample_id,species,abundance
        detailed: file("${id}.blast_detailed.tsv"),           // best hit(s) per ASV of this sample
    )

    topic:
    tuple(task.process, 'R', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    """
    format_blast.R ${counts} ${blast} ${model_metadata} ${id} ${id}
    """

    stub:
    """
    echo "placeholder,sample_id,species,abundance" > ${id}.annotated_abundances.csv
    touch ${id}.blast_detailed.tsv
    """
}
