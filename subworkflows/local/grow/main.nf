//
// Phase 2 - growth: simulate every sample on every diet that applies to it
//
nextflow.enable.types = true

include { MICOM_GROW } from '../../../modules/local/micom/grow/main'
include { Growth     } from '../../../workflows/types'
include { Taxonomy   } from '../../../workflows/types'

workflow GROW {

    take:
    taxonomy: Channel<Taxonomy>
    diet_rows: Channel<Map>          // [ id, medium, diet ] per sample and diet
    micom_db: Value<Path>

    main:
    diets = diet_rows.map { row -> record(id: row.id as String, medium: row.medium as String, diet: row.diet as Path) }
    growth = MICOM_GROW(taxonomy.join(diets, by: 'id'), micom_db)

    emit:
    growth: Channel<Growth> = growth
}
