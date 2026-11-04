#!/usr/bin/env Rscript
# Best BLAST hit(s) per ASV, then ASV counts -> species abundances through the model metadata.
# usage: format_blast.R <counts.csv> <blast.txt> <model_metadata> <sample_id> <prefix>
# writes <prefix>.blast_detailed.tsv and <prefix>.annotated_abundances.csv
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]; blast_file <- args[2]; metadata_file <- args[3]; sample <- args[4]; prefix <- args[5]

# counts and BLAST hits cover all samples: keep this sample's ASVs
counts <- read.csv(counts_file, colClasses = c(sample_id = "character", asv = "character")) %>%
  filter(sample_id == sample)

blast_names <- c("ASVId", "taxid", "accession", "species", "query_length", "align_length",
                 "nident", "pident", "bitscore", "evalue", "score")
blast_types <- c("character", "character", "character", "character", "integer", "integer", "integer",
                 "numeric", "numeric", "numeric", "numeric")
results <- if (file.size(blast_file) == 0) {
  read.csv(text = "", header = FALSE, col.names = blast_names, colClasses = blast_types)   # no hits at all
} else {
  read.csv(blast_file, sep = "\t", header = FALSE, col.names = blast_names, colClasses = blast_types)
}
results <- filter(results, ASVId %in% counts$asv)

# all hits with the highest bitscore of each ASV (ties stay)
top <- results %>% group_by(ASVId) %>% slice_max(bitscore, n = 1, with_ties = TRUE) %>% ungroup()

# one line per ASV: the tied accessions joined by ';', then the alignment length and identity
top %>%
  group_by(ASVId) %>%
  arrange(accession) %>%
  mutate(name = paste(paste(unique(accession), collapse = ";"), align_length[1], pident[1], sep = ";")) %>%
  slice(1) %>%
  ungroup() %>%
  write.table(paste0(prefix, ".blast_detailed.tsv"), sep = "\t", row.names = FALSE)

# Genome_file -> accession: file extension, version and a trailing _genomic are ignored
# (GCF_000005845.2_genomic.fna.gz matches the BLAST accession GCF_000005845.2)
key <- read.delim(metadata_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(key) != 2) stop("model_metadata must have two tab-separated columns: Genome_file, clade_name")
key <- key %>%
  transmute(acc = sub("(\\.[0-9]+)?(_genomic)?$", "", sub("\\.(fna|fa|fasta|ffn)(\\.gz)?$", "", basename(V1))),
            species = sub(".*s__", "", V2),
            bacteria = startsWith(V2, "d__Bacteria")) %>%
  filter(bacteria) %>%
  distinct(acc, species)

# each ASV count goes to the species of its best hits, split equally when several species tie
species_of_asv <- top %>%
  transmute(ASVId, acc = sub("\\.[0-9]+$", "", accession)) %>%   # BLAST's own species title is not used
  inner_join(key, by = "acc") %>%
  distinct(ASVId, species)
abundances <- counts %>%
  inner_join(species_of_asv, by = c("asv" = "ASVId")) %>%
  group_by(asv) %>% mutate(count = count / n()) %>% ungroup() %>%
  group_by(species) %>% summarize(abundance = sum(count), .groups = "drop")

message(sprintf("%.1f%% of the ASV counts could not be mapped to a species via the model metadata",
                100 * (1 - sum(abundances$abundance) / max(sum(counts$count), 1))))
abundances %>%
  mutate(abundance = abundance / sum(abundance), sample_id = sample, placeholder = "a") %>%
  select(placeholder, sample_id, species, abundance) %>%
  arrange(desc(abundance)) %>%
  write.table(paste0(prefix, ".annotated_abundances.csv"), row.names = FALSE, sep = ",", quote = FALSE)
