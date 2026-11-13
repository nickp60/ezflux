//
// Phase 1 - preprocessing: turn reads or an ASV table into micom taxonomy tables
//
//   reads     -> (CAT_FASTQ if >1 run) -> sylph profile ---------------------------------------\
//   ASV table -> FASTA -> blastn (16S db) -> best hit + model_metadata, per sample ------------+-> MICOM_TAXONOMY
//
// MICOM_TAXONOMY looks every sylph genome up in the model metadata by its file path and maps its species
// to a model of the micom database.
//
nextflow.enable.types = true

include { CAT_FASTQ        } from '../../../modules/nf-core/cat/fastq/main'
include { SYLPH_PROFILE    } from '../../../modules/nf-core/sylph/profile/main'
include { BLAST_BLASTN     } from '../../../modules/nf-core/blast/blastn/main'
include { ASV_TOFASTA      } from '../../../modules/local/asv/tofasta/main'
include { BLAST_FORMAT     } from '../../../modules/local/blast/format/main'
include { MICOM_TAXONOMY   } from '../../../modules/local/micom/taxonomy/main'
include { Taxonomy         } from '../../../workflows/types'

workflow PREPROCESSING {

    take:
    reads_rows: Channel<Map>         // [ id, single_end, n_runs, fastqs ] per sample with reads
    asv_rows: Channel<Map>           // [ id ] per sample of the ASV table
    seqtab: Channel<Path>            // the ASV table (empty without --seqtab_input)
    sylph_db: Channel<Path>          // empty without --input
    blast_db: Channel<Path>          // empty without --seqtab_input
    model_metadata: Value<Path>
    micom_db: Value<Path>

    main:

    //
    // Reads route: the runs of a sample are concatenated first
    //
    cat = CAT_FASTQ(
        reads_rows
            .filter { row -> row.n_runs > 1 }
            .map { row -> tuple([ id: row.id, single_end: row.single_end ], row.fastqs as List<Path>) }
    )
    merged = cat.reads.map { meta, fastq -> tuple(meta, fastq instanceof List ? fastq.sort() : [ fastq ]) }
    unmerged = reads_rows
        .filter { row -> row.n_runs == 1 }
        .map { row -> tuple([ id: row.id, single_end: row.single_end ], row.fastqs as List<Path>) }

    profile = SYLPH_PROFILE(unmerged.mix(merged), sylph_db)
    sylph = profile.profile_out.map { meta, tsv -> record(id: meta.id as String, abundances: tsv, source: 'sylph') }

    //
    // ASV route: one FASTA and one blastn search for the whole table, then best hits per sample
    //
    asv = ASV_TOFASTA(seqtab)
    blast = BLAST_BLASTN(
        asv.map { r -> tuple([ id: 'seqtab' ], r.fasta) },
        blast_db.map { db -> tuple([ id: 'blast_db' ], db) },
        [], [], false
    )
    hits = asv.map { r -> r.counts }.combine(blast.txt.map { _meta, txt -> txt })
    formatted = BLAST_FORMAT(
        asv_rows.combine(hits).map { row, counts, txt -> record(id: row.id as String, counts: counts, blast: txt) },
        model_metadata
    )

    //
    // Both routes end in one table per sample that is matched to the models of the micom database
    //
    taxonomy = MICOM_TAXONOMY(
        sylph.mix(formatted.map { r -> record(id: r.id, abundances: r.abundances, source: 'species') }),
        model_metadata,
        micom_db
    )

    emit:
    taxonomy: Channel<Taxonomy> = taxonomy.map { r -> record(id: r.id, taxonomy: r.taxonomy) }
}
