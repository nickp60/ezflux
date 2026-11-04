#!/usr/bin/env Rscript
# Convert a wide ASV sequence table (rows = samples, first column = sample id, other columns =
# ASV sequences, values = counts) to a FASTA file and a long count table.
# usage: asv_to_fasta.R <seqtab> <out_prefix>
# writes <prefix>.fasta (ids ASV1, ASV2, ... in column order) and <prefix>.counts.csv
# (sample_id, asv, count; zero counts dropped)
args <- commandArgs(trailingOnly = TRUE)
seqtab_file <- args[1]; prefix <- args[2]

header <- readLines(seqtab_file, n = 1)
sep <- if (grepl("\t", header)) "\t" else ","
n_cols <- length(strsplit(header, sep, fixed = TRUE)[[1]])
tab <- read.table(seqtab_file, sep = sep, header = TRUE, quote = "\"", check.names = FALSE,
                  comment.char = "", stringsAsFactors = FALSE,
                  colClasses = c("character", rep("numeric", n_cols - 1)))
sample_ids <- tab[[1]]
seqs <- colnames(tab)[-1]
is_dna <- function(x) all(grepl("^[ACGTUNacgtun]+$", x))
if (!is_dna(seqs)) {
  if (is_dna(sample_ids)) {
    stop("the table looks transposed: ASV sequences must be the columns and the sample ids the first column")
  }
  stop("the column names (after the first) must be ASV sequences; found e.g. '", seqs[!grepl("^[ACGTUNacgtun]+$", seqs)][1], "'")
}
if (anyDuplicated(sample_ids)) stop("duplicated sample ids in the first column: ", paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", "))
if (anyDuplicated(seqs)) stop("duplicated ASV sequences in the header")
if (any(is.na(as.matrix(tab[-1])))) stop("missing or non-numeric counts in the table")

ids <- paste0("ASV", seq_along(seqs))
writeLines(as.vector(rbind(paste0(">", ids), seqs)), paste0(prefix, ".fasta"))

m <- as.matrix(tab[-1])
idx <- which(m > 0, arr.ind = TRUE)
long <- data.frame(sample_id = sample_ids[idx[, 1]], asv = ids[idx[, 2]], count = m[idx], stringsAsFactors = FALSE)
long <- long[order(match(long$sample_id, sample_ids), long$asv), ]
write.csv(long, paste0(prefix, ".counts.csv"), row.names = FALSE, quote = FALSE)
message(length(seqs), " ASVs and ", length(sample_ids), " samples (", nrow(long), " non-zero counts)")
