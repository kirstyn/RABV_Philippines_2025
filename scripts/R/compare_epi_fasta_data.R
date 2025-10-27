library(seqinr)

# Load files
seq <- read.fasta("processed_data/processed_sequences/20251023_PHL_all.aln.fasta")
meta <- read.csv("processed_data/processed_metadata/gathered_metadata/27Oct25_gathered_metadata_n812_raddl_and_manual_Corrected.csv")


# Extract sequence names
seq_names <- names(seq)

# ---- 1. Sequences with no matching metadata ----
matched_to_sampleID <- seq_names %in% meta$Sample_ID
matched_to_accession <- seq_names %in% meta$Accession
matched_any <- matched_to_sampleID | matched_to_accession
unmatched_seqs <- seq_names[!matched_any]

# ---- 2. Metadata with no matching sequence ----
matched_to_seq <- meta$Sample_ID %in% seq_names | meta$Accession %in% seq_names
unmatched_meta <- meta[!matched_to_seq, c("Sample_ID", "Accession", "Source")]

# ---- Summary ----
cat("Sequences matched to metadata:", sum(matched_any), "of", length(seq_names), "\n")
cat("Unmatched sequences:", length(unmatched_seqs), "\n")
cat("Metadata entries matched to sequences:", sum(matched_to_seq), "of", nrow(meta), "\n")
cat("Unmatched metadata entries:", nrow(unmatched_meta), "\n")

# Optional: inspect results
unmatched_seqs        # names of sequences missing metadata
unmatched_meta        # metadata rows missing sequences #these are all low coverage sequences

# Remove the unmatched meta that is confirmed to be low cov
# Keep only metadata rows that have a matching sequence
meta_filtered <- meta[matched_to_seq, ]

# Optional: check
nrow(meta_filtered)  # should be 795
write.csv(meta_filtered,"processed_data/processed_metadata/gathered_metadata/27Oct25_gathered_metadata_n812_raddl_and_manual_Corrected_filteredTon795.csv")
