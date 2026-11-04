//
// Phase 3 - merging: combine all growth results, export the tables and draw the plots
//
nextflow.enable.types = true

include { MICOM_MERGE } from '../../../modules/local/micom/merge/main'
include { Growth      } from '../../../workflows/types'

workflow MERGE {

    take:
    growth: Channel<Growth>

    main:
    merged = MICOM_MERGE(
        growth.map { g -> g.zip }.collect(),
        growth.map { g -> g.build_manifest }.collect()
    )

    emit:
    plots: Channel<Path> = merged.flatMap { r -> r.plots }   // static PNGs for MultiQC
}
