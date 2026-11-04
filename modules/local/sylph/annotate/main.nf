nextflow.enable.types = true

process SYLPH_ANNOTATE {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-tidyverse:2.0.0'
        : 'quay.io/biocontainers/r-tidyverse:2.0.0'}"

    input:
    record(id: String, profile: Path)   // raw `sylph profile` TSV
    model_metadata: Path                // Genome_file <tab> clade_name, shipped with the model database

    output:
    record(id: id, abundances: file("${id}.annotated_abundances.csv"))   // placeholder,sample_id,species,abundance

    topic:
    tuple(task.process, 'R', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    """
    annotate_sylph.R ${profile} ${model_metadata} ${id} ${id}.annotated_abundances.csv
    """

    stub:
    """
    echo "placeholder,sample_id,species,abundance" > ${id}.annotated_abundances.csv
    """
}
