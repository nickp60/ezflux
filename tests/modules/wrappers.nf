//
// Thin workflows around the typed local modules so that nf-test can call them with plain paths
//
nextflow.enable.types = true

include { ASV_TOFASTA    } from '../../modules/local/asv/tofasta/main'
include { BLAST_FORMAT   } from '../../modules/local/blast/format/main'
include { SYLPH_ANNOTATE } from '../../modules/local/sylph/annotate/main'
include { MICOM_TAXONOMY } from '../../modules/local/micom/taxonomy/main'
include { MICOM_GROW     } from '../../modules/local/micom/grow/main'
include { MICOM_MERGE    } from '../../modules/local/micom/merge/main'

workflow TEST_ASV_TOFASTA {
    take:
    seqtab: Path

    main:
    out = ASV_TOFASTA(channel.value(seqtab))

    emit:
    fasta: Value<Path> = out.map { r -> r.fasta }
    counts: Value<Path> = out.map { r -> r.counts }
}

workflow TEST_BLAST_FORMAT {
    take:
    id: String
    counts: Path
    blast: Path
    model_metadata: Path

    main:
    out = BLAST_FORMAT(channel.value(record(id: id, counts: counts, blast: blast)), channel.value(model_metadata))

    emit:
    abundances: Value<Path> = out.map { r -> r.abundances }
    detailed: Value<Path> = out.map { r -> r.detailed }
}

workflow TEST_SYLPH_ANNOTATE {
    take:
    id: String
    profile: Path
    model_metadata: Path

    main:
    out = SYLPH_ANNOTATE(channel.value(record(id: id, profile: profile)), channel.value(model_metadata))

    emit:
    abundances: Value<Path> = out.map { r -> r.abundances }
}

workflow TEST_MICOM_TAXONOMY {
    take:
    id: String
    abundances: Path
    micom_db: Path

    main:
    out = MICOM_TAXONOMY(channel.value(record(id: id, abundances: abundances)), channel.value(micom_db))

    emit:
    taxonomy: Value<Path> = out.map { r -> r.taxonomy }
    report: Value<Path> = out.map { r -> r.report }
}

workflow TEST_MICOM_GROW {
    take:
    id: String
    medium: String
    taxonomy: Path
    diet: Path
    model_db: Path

    main:
    out = MICOM_GROW(channel.value(record(id: id, medium: medium, taxonomy: taxonomy, diet: diet)), channel.value(model_db))

    emit:
    zip: Value<Path> = out.map { r -> r.zip }
    build_manifest: Value<Path> = out.map { r -> r.build_manifest }
}

workflow TEST_MICOM_MERGE {
    take:
    zips: Channel<Path>
    manifests: Channel<Path>

    main:
    out = MICOM_MERGE(zips.collect(), manifests.collect())

    emit:
    result: Value<Record> = out
}
