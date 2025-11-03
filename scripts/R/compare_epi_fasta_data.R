# ===============================================================
#  PHILIPPINES SEQUENCE–METADATA MATCHING AND CLEANING PIPELINE
#  Author: Kirstyn Brunker
#  Date created: 23 Oct 2025
#  Purpose: Match sequences and metadata, remove low coverage,
#           retain only records with valid dates, and export
#           cleaned metadata and FASTA files.
# ===============================================================

# ---- Load libraries ----
library(seqinr)
library(dplyr)
library(lubridate)
library(stringr)

# ---- Input files ----
seq_file  <- "processed_data/processed_sequences/241025_PHL_all_n794.fasta"
meta_file <- "processed_data/processed_metadata/gathered_metadata/final/29Oct25_gathered_metadata_n811_raddl_and_manual_Corrected_stdGeo.csv"

# ---- Output directories ----
out_meta_dir <- "processed_data/processed_metadata/gathered_metadata/final/"
out_seq_dir  <- "processed_data/processed_sequences/"

# Create directories if missing
dir.create(out_meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_seq_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading sequence and metadata files...\n")

seq  <- read.fasta(seq_file)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Clean whitespace from names
seq_names <- trimws(names(seq))
meta <- meta %>%
  mutate(across(c(Sample_ID, Accession), ~ trimws(.)))

cat("Sequences loaded:", length(seq_names), "\n")
cat("Metadata loaded:", nrow(meta), "\n\n")

# ---- Identify matches and mismatches ----
cat("Identifying sequence–metadata matches...\n")

matched_to_sampleID <- seq_names %in% meta$Sample_ID
matched_to_accession <- seq_names %in% meta$Accession
matched_any <- matched_to_sampleID | matched_to_accession
unmatched_seqs <- seq_names[!matched_any]

matched_to_seq <- meta$Sample_ID %in% seq_names | meta$Accession %in% seq_names
unmatched_meta <- meta[!matched_to_seq, c("Sample_ID", "Accession", "Source")]

cat("✅ Sequences matched to metadata:", sum(matched_any), "of", length(seq_names), "\n")
cat("❌ Unmatched sequences:", length(unmatched_seqs), "\n")
cat("✅ Metadata matched to sequences:", sum(matched_to_seq), "of", nrow(meta), "\n")
cat("❌ Unmatched metadata:", nrow(unmatched_meta), "\n\n")

# ---- Remove confirmed low-coverage (unmatched) metadata and Create a consistent 'seq_name' column ----
meta_filtered <- meta[matched_to_seq, ]
cat("Metadata retained after removing unmatched entries:", nrow(meta_filtered), "\n\n")

meta_filtered <- meta_filtered %>%
  mutate(
    seq_name = case_when(
      Sample_ID %in% seq_names ~ Sample_ID,
      Accession %in% seq_names ~ Accession,
      TRUE ~ NA_character_
    )
  )

meta_unmatched <- meta_filtered %>%
  filter(is.na(seq_name)) %>%
  select(Sample_ID, Accession)

cat(nrow(meta_unmatched), "metadata entries still didn’t match after filtering.\n\n")

# Save sequence-epi matched metadata
matched_meta_path <- file.path(
  out_meta_dir,
  paste0("PHL_metadata_matchedToSeq_", Sys.Date(), "_n", nrow(meta_filtered), ".csv")
)
write.csv(meta_filtered, matched_meta_path, row.names = FALSE)
cat("Matched metadata saved to:\n", matched_meta_path, "\n\n")

# ---- Filter for sequences with valid dates ----
has_date <- !is.na(meta_filtered$Preferred_date) & meta_filtered$Preferred_date != "NA"
seq_filtered <- seq[has_date]
matched_filtered <- meta_filtered[has_date, ]

cat("Sequences retained with valid dates:", length(seq_filtered), "\n\n")

# ---- Save filtered metadata ----
filtered_meta_path <- file.path(
  out_meta_dir,
  paste0("PHL_metadata_matchedWithDates_", Sys.Date(), "_n", nrow(matched_filtered), ".csv")
)
write.csv(matched_filtered, filtered_meta_path, row.names = FALSE)
cat("Filtered metadata (with valid dates) saved to:\n", filtered_meta_path, "\n\n")

# ---- Save filtered FASTA (unaltered names) ----
fasta_no_dates_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDates_", Sys.Date(), "_n", length(seq_filtered), ".fasta")
)
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_no_dates_path)
cat("Filtered FASTA saved (unaltered names):\n", fasta_no_dates_path, "\n\n")

# ----  Append dates to sequence names ----
#  prepare date strings for names
clean_dates <- gsub('["=]', '', matched_filtered$Preferred_date)   # remove quotes and equal signs
clean_dates <- trimws(clean_dates)
new_names_date <- paste0(matched_filtered$seq_name, "_", clean_dates)
names(seq_filtered) <- new_names_date

fasta_with_dates_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDatesLabeled_", Sys.Date(), "_n", length(seq_filtered), ".aln.fasta")
)
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_with_dates_path)
cat("FASTA with dates appended to sequence names saved to:\n", fasta_with_dates_path, "\n")

# ----  Append dates and region to sequence names ----
new_names_date_region <- paste0(matched_filtered$seq_name, "_", clean_dates, "_", matched_filtered$Region_std)
names(seq_filtered) <- new_names_date_region

fasta_with_dates_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDatesRegionsLabeled_", Sys.Date(), "_n", length(seq_filtered), ".aln.fasta")
)
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_with_dates_path)
cat("FASTA with dates appended to sequence names saved to:\n", fasta_with_dates_path, "\n")

# ---- 13. Summary ----
cat("\n================ SUMMARY ================\n")
cat("Total sequences:", length(seq_names), "\n")
cat("Sequences matched to metadata:", sum(matched_any), "\n")
cat("Sequences with valid dates:", length(seq_filtered), "\n")
cat("Final metadata records exported:", nrow(matched_filtered), "\n")
cat("=========================================\n")
