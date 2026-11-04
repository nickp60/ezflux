nextflow.enable.types = true

process ASV_TOFASTA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-base:4.3.3'
        : 'quay.io/biocontainers/r-base:4.3.3'}"

    input:
    seqtab: Path   // wide table: sample ids in the first column, ASV sequences as the other columns, counts as values

    output:
    record(
        fasta: file("seqtab.fasta"),        // ASV1, ASV2, ... in column order
        counts: file("seqtab.counts.csv"),  // all samples: sample_id, asv (= FASTA id), count
    )

    topic:
    tuple(task.process, 'R', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    """
    asv_to_fasta.R ${seqtab} seqtab
    """

    stub:
    """
    printf ">ASV1\\nACGT\\n" > seqtab.fasta
    printf "sample_id,asv,count\\nS1,ASV1,1\\n" > seqtab.counts.csv
    """
}
