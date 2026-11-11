//
// Phase 3 - merging: combine all growth results, export the tables and draw the plots
//
nextflow.enable.types = true

include { MICOM_MERGE } from '../../../modules/local/micom/merge/main'

workflow MERGE {

    take:
    growth_zips: Channel<Path>        // growth results (zips) of this run or of earlier runs
    build_manifests: Channel<Path>    // optional: found fractions of these results (may be empty)

    main:
    merged = MICOM_MERGE(growth_zips.collect(), build_manifests.collect())

    emit:
    plots: Channel<Path> = merged.flatMap { r -> r.plots }   // static PNGs for MultiQC
}
