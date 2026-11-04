#!/usr/bin/env Rscript
# Re-annotate raw `sylph profile` output with the model metadata.
# usage: annotate_sylph.R <sylph.tsv> <model_metadata> <sample_id> <out.csv>
# model_metadata is tab-separated without header: Genome_file, clade_name (GTDB lineage ending in s__<species>)
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
sylph_file <- args[1]; metadata_file <- args[2]; sample <- args[3]; out_file <- args[4]

key <- read.delim(metadata_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(key) != 2) stop("model_metadata must have two tab-separated columns: Genome_file, clade_name")
colnames(key) <- c("Genome_file", "clade_name")
key$Genome_file <- basename(key$Genome_file)

profile <- read.delim(sylph_file, stringsAsFactors = FALSE) %>% mutate(Genome_file = basename(Genome_file))
unmatched <- sum(!(profile$Genome_file %in% key$Genome_file))
if (unmatched > 0) message(unmatched, " sylph genome(s) are not in the model metadata and were dropped")

profile %>%
  left_join(key, by = "Genome_file") %>%
  filter(startsWith(clade_name, "d__Bacteria")) %>%
  mutate(sample_id = sample,                      # one sylph profile per sample
         species = sub(".*s__", "", clade_name)) %>%
  group_by(sample_id, species) %>%
  summarize(abundance = sum(Taxonomic_abundance), .groups = "drop") %>%
  mutate(abundance = abundance / sum(abundance), placeholder = "a") %>%
  select(placeholder, sample_id, species, abundance) %>%
  filter(abundance > 0) %>%
  arrange(desc(abundance)) %>%
  write.table(out_file, row.names = FALSE, sep = ",", quote = FALSE)
