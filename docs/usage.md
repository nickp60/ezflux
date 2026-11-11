# nf-core/ezflux: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/ezflux/usage](https://nf-co.re/ezflux/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/ezflux predicts which metabolites the members of a gut microbiome take up and secrete. It profiles
your samples, builds a [MICOM](https://micom-dev.github.io/micom/) community model per sample, simulates growth
on one or more diets and merges everything into exchange tables and niche plots. It runs in three phases:

1. **Preprocessing**: samples with reads (runs of one sample are concatenated with `CAT_FASTQ`) are profiled with
   [sylph](https://github.com/bluenote-1577/sylph) (`SYLPH_PROFILE`) and re-annotated with the model metadata
   (`SYLPH_ANNOTATE`). ASV tables are written to FASTA (`ASV_TOFASTA`), searched against a 16S database with
   [blastn](https://blast.ncbi.nlm.nih.gov/) (`BLAST_BLASTN`, one search for the whole table) and the best hits are formatted
   and mapped to species per sample (`BLAST_FORMAT`). Both routes produce `annotated_abundances.csv`, which `MICOM_TAXONOMY` matches to the models in your
   MICOM database and writes as a MICOM taxonomy table.
2. **Growth**: `MICOM_GROW` builds each community and simulates growth with cooperative trade-off FBA, once for every
   sample and diet that applies to it.
3. **Merging**: `MICOM_MERGE` combines all growth results, exports the exchanges and draws the import/export niche plots.

You need EITHER a reads sample sheet (`--input`) OR a wide ASV table (`--seqtab_input`), a MICOM model database (`--micom_db`) and the metadata that ships with it (`--model_metadata`). Diets come from a diet sheet (`--media`) or, without one, from `--medium`. Reads need `--sylph_db`, `--seqtab_input` needs `--blast_db`. Alternatively, `--growth_input` merges growth results of earlier runs and needs none of these (see [Merging results of earlier runs](#merging-results-of-earlier-runs---growth_input)).

## Samplesheet input

```bash
--input '[path to samplesheet file]'
```

A CSV with one row per sample (or per sequencing run of a sample):

```csv title="samplesheet.csv"
sample,run_accession,fastq_1,fastq_2
2612,ERR5766176,/data/2612_run1_1.fastq.gz,/data/2612_run1_2.fastq.gz
2612,ERR5766176_B,/data/2612_run2_1.fastq.gz,/data/2612_run2_2.fastq.gz
2613,,/data/2613_1.fastq.gz,
```

| Column          | Description                                                                                                  |
| --------------- | ------------------------------------------------------------------------------------------------------------ |
| `sample`        | Sample name or number, no spaces. Numbers are kept as text (`007` stays `007`).                              |
| `run_accession` | Optional. Rows with the same `sample` and different `run_accession` are concatenated with `cat/fastq`.      |
| `fastq_1`       | Optional. Gzipped FastQ. Leave `fastq_2` empty for single-end data.                                          |
| `fastq_2`       | Optional. Gzipped FastQ for the second read. Needs `fastq_1`.                                                |

Every row needs `fastq_1`; all runs of a sample must be all single-end or all paired-end. Several rows for one sample need distinct `run_accession` values. Other columns (for
example `instrument_platform` and `fasta` in nf-core/taxprofiler sheets) are ignored. Reads need `--sylph_db`; `--seqtab_input` needs
`--blast_db`; both need `--model_metadata`.

> Relative paths in a samplesheet are resolved against the directory you launch Nextflow from. Use absolute paths or URLs.
### Sylph Database

Instructions for building the sylph database can be found  [here](https://github.com/vdblab/resources/tree/main/workflow/sylphflux). In short, genomes are downloaded when available, and short reads were assembled with skesa. the genomes/assemblies were then sketched into a database with Sylph. The model metadata file links the file path to the taxonomy.

### Model metadata

`--model_metadata` is the file distributed with the model database. It has no header and two columns, `Genome_file` and
`clade_name` (a GTDB-style lineage ending in `s__<species>`), separated by tabs:

```tsv
/db/genomes/GCF_000005845.2_genomic.fna.gz	d__Bacteria;p__Pseudomonadota;...;g__Escherichia;s__Escherichia coli
```

For reads, the genome sylph reports is looked up by its file name. For ASVs, the BLAST subject accession (`saccver`) is looked
up by the file name without the `.fna`/`.fa`/`.fasta` (and `.gz`) suffix, tolerating a different version suffix. Entries that
are not `d__Bacteria` are dropped.

### Annotating reads

The raw `sylph profile` table is joined to the model metadata, restricted to bacteria, summed per species and renormalised to
relative abundances (`<sample>.annotated_abundances.csv` with columns `placeholder,sample_id,species,abundance`).
`MICOM_TAXONOMY` then matches these species to the models of `--micom_db` (case, spaces and GTDB `_A` suffixes are ignored; on
a strain-level database a species' abundance is split over its models; the database's summary rank must be species, as in AGORA) and renormalises over the matched species. Species
without a model are dropped and listed in `preprocessing/taxonomy/<sample>.annotation.tsv`.

### ASV table and 16S annotation

```bash
--seqtab_input '[path to ASV sequence table]'
```

The standard (dada2-style) **wide** table: one row per sample, sample ids in the first column, one column per ASV named by its
sequence, counts as values. CSV or TSV (detected from the header):

```csv title="seqtab.csv"
"","GCTAAAGACAATTACATAACATACACTTCAGCACG","TCGGAATCACTGGGCGTAAAGCGCACGCAGGCGG"
"sample1",600,300
"sample2",500,0
```

The **sample ids are read from the first column** of the table, so no samplesheet is needed (and `--input` must not be given).
Ids must be unique and contain no spaces; numeric ids are kept as text. A transposed table (sequences as rows) is detected
and rejected with an error.


All ASVs are written to one FASTA (`ASV1`, `ASV2`, ... in column order; zero counts are dropped from the long count table)
and searched once with `blastn -max_target_seqs 500 -outfmt "6 qseqid staxids saccver stitle qlen length nident pident
bitscore evalue score"` against the database in `--blast_db` (a directory holding the `makeblastdb` output). Per sample, the
hit(s) with the highest bitscore are kept for each of its ASVs (`<sample>.blast_detailed.tsv`). The ASV counts are then
converted to species abundances: each ASV's count goes to the species of its best hits, split equally when several species
tie, and the species are normalised to relative abundances. Hits whose accession is not in `--model_metadata` are lost; the
share is printed in the task log.

### BLAST DB
Follow the instructions [here](https://github.com/nickp60/micom-snakemake/blob/main/workflow/rules/format16S.smk) for creating the blast database from the model genomes.
## Media input

```bash
--media '[path to media sheet]'
```

Optional. A CSV with one row per diet. Without `--media`, every sample is grown on the diet in `--medium` (a q2-micom `.qza`
file or URL; by default the western diet, `https://zenodo.org/record/3755182/files/western_diet_gut.qza?download=1`).

Optional. A CSV with one row per diet. Without `--media`, every sample is grown on the diet in `--medium` (a q2-micom `.qza` or a [CSV/TSV medium](#medium-formats)

```csv title="media.csv"
medium_name,medium_path,sample_id
western_diet,/data/western_diet_gut_agora.qza,
high_fiber,/data/high_fiber_gut_agora.qza,GUT_ASV
```

| Column        | Description                                                                                       |
| ------------- | ------------------------------------------------------------------------------------------------- |
| `medium_name` | Name of the diet, no spaces. Used in output file names.                                           |
| `medium_path` | A q2-micom medium artifact (`.qza`).                                                              |
| `sample_id`   | Optional. Restricts the diet to this sample. Empty means the diet applies to every sample.        |

Every sample must end up with at least one diet, either its own or one without a `sample_id`, and every `sample_id` here must exist
in `--input`; the pipeline stops with an error otherwise. A sample on several diets is simulated once per diet, and
`<sample>__<medium>` becomes the sample id in the merged results.

## Merging results of earlier runs (`--growth_input`)

A third input type skips the whole analysis (reading, annotation, model building and growth) and only runs the merge step and MultiQC on growth results that were produced earlier. Use it to combine runs, for example different batches of a project. `--input`, `--seqtab_input` and `--growth_input` are mutually exclusive, and `--micom_db`, `--model_metadata`, `--media` and the sylph/BLAST databases are not needed.

```csv title="growth_samplesheet.csv"
sample,growth_zip,build_manifest
batch1,/results/batch1/merged/batch1.growth.zip,/results/batch1/merged/batch1.found_fraction.csv
batch2,/results/batch2/merged/batch2.growth.zip,/results/batch2/merged/batch2.found_fraction.csv
```

| Column           | Description                                                                                                  |
| ---------------- | ------------------------------------------------------------------------------------------------------------ |
| `sample`         | Label for the row (no spaces). Must be unique together with `growth_zip`.                                    |
| `growth_zip`     | Required. A MICOM growth results zip: either the `*.growth.zip` of a whole run, or the per-sample zips in `growth/`. |
| `build_manifest` | Optional. Build manifest or `*.found_fraction.csv` matching the zip, used for the found-fraction plot.        |

```bash
nextflow run nf-core/ezflux --growth_input growth_samplesheet.csv --name combined --outdir results -profile docker
```

Each sample and medium may only be merged once: if the same `<sample>__<medium>` id appears in two zips, the merge fails and names the id. Without any build manifest the found-fraction table is empty and its plot reports that no manifests were provided.

## MICOM database and solver

`--micom_db` is a MICOM model database (`.qza` or `.zip`), for example AGORA. The taxonomy tables are matched against its
taxonomy up to its summary rank.

MICOM needs a quadratic-programming solver. The default `--micom_solver gurobi` needs a **Gurobi license visible to the
compute nodes** (for example `module load gurobi`, `GRB_LICENSE_FILE`, or a token server). Without it Gurobi falls back to a
restricted license that is limited to about 2000 variables and fails on real models with "Model too large for size-limited
license". `--micom_solver cplex` and `--micom_solver hybrid` (open source, much slower) are alternatives.

### Retries and numerical trouble

`MICOM_GROW` runs one sample on one diet. It fails if the sample cannot be solved, and Nextflow retries it once. The retry sets Gurobi's `NumericFocus` to
`--micom_numeric_focus_retry` (default 3, slower but more robust); the first attempt uses MICOM's default (2). Set
`--micom_numeric_focus` to use a different value from the start. Related options:

- `--micom_coupling false` switches off MICOM's resource coupling constraints, which are the most common cause of numerical
  trouble in large communities and are much cheaper to solve without.
- `--micom_presolve true|false|auto` overrides the solver presolve (MICOM's default is on).

MICOM only sets these Gurobi options internally, so they are applied by overriding its solver setup inside the task.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/ezflux --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/ezflux -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/ezflux
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/ezflux releases page](https://github.com/nf-core/ezflux/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
