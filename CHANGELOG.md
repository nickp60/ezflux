# nf-core/ezflux: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/ezflux, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Initial release of nf-core/ezflux: reads (sylph) or ASV tables (blastn) -> MICOM community models -> growth on one or more diets -> merged exchanges, rates and niche plots.
- `--input` reads samplesheet OR `--seqtab_input` wide ASV table (sample ids from its first column), exactly one of the two; `--input` columns `sample,run_accession,fastq_1,fastq_2` ( `sample` may be numeric; runs of a sample are merged with `cat/fastq`), and optional `--media` samplesheet (`medium_name,medium_path,sample_id`) with cross-validation. Without `--media`, `--medium` (default: western diet on Zenodo) applies to all samples.
- `--model_metadata` (taxonomy -> model files, shipped with the model database), `--sylph_db` and `--blast_db`.
- Local modules `SYLPH_ANNOTATE`, `ASV_TOFASTA`, `BLAST_FORMAT`, `MICOM_TAXONOMY`, `MICOM_GROW` (retries with a higher Gurobi `NumericFocus`) and `MICOM_MERGE`.
- Plots of growth rates (`plot_growth`), the found fraction of each community in the model database, and the import/export niches; static PNGs of all four are added to the MultiQC report.
- Test profiles `test` (nf-core/funcprofiler samplesheet, local AGORA2 files) and `test_local` (bundled toy data).

### `Changed`

- Local code is written in Nextflow typed syntax (records, typed processes and subworkflows, `topic:` versions, no `when:`). Template files and nf-core modules stay in the template's syntax.
- `MICOM_GROW` and `MICOM_MERGE` take their options from `task.ext.args`, built from the pipeline parameters in `conf/modules.config`; the Python code lives in the module `resources/usr/bin` scripts.
- `MICOM_GROW` handles exactly one sample on one diet; the sample id in all results is always `<sample>__<medium>`. `--micom_allow_partial` was removed.
- The `MICOM_MERGE` manifest was removed: merged samples are identified by the ids inside the growth zips.
- `--model_metadata` is tab-separated only; taxonomy matching is species level.

### `Fixed`

### `Dependencies`

### `Deprecated`
