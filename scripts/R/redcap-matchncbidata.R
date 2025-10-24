# ============================================================
# Compare public NCBI data with in-house REDCap data
# ============================================================

# ---- Load packages ----
library(dplyr)
library(stringr)
library(seqinr)

# ---- 1. Input files ----
ncbi_file <- "raw_data/vgtk_ncbi_data/150725_metadata_coverage_90_country_Philippines.csv"
redcap_file <- "raw_data/redcap_sequences_and_metadata/redcap_download_20251023_134520251023_1345redcap_meta_phl.csv"
fasta_file <- "raw_data/vgtk_ncbi_data/150725_metadata_coverage_90_country_Philippines_sequences.fa"

# ---- 2. Output directories ----
meta_out_dir <- "processed_data/processed_metadata"
seq_out_dir  <- "processed_data/processed_sequences"

# ---- 3. Read data ----
ncbi <- read.csv(ncbi_file, stringsAsFactors = FALSE)
redcap <- read.csv(redcap_file, stringsAsFactors = FALSE)

# ============================================================
# PART A: Filter NCBI metadata (exclude team & exclusion_status)
# ============================================================

ncbi_additional <- ncbi %>%
  filter(
    !str_detect(authors, regex("brunker", ignore_case = TRUE)),
    is.na(exclusion_status) | exclusion_status != 1
  )

# ---- Write filtered NCBI metadata ----
ncbi_output <- file.path(meta_out_dir, "150725_metadata_coverage_90_country_Philippines_filtered.csv")
write.csv(ncbi_additional, ncbi_output, row.names = FALSE)

# ============================================================
# PART B: Filter NCBI FASTA sequences
# ============================================================

# ---- Read full FASTA ----
ncbi_seq <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)

# ---- Keep only sequences present in filtered metadata ----
seq_names_to_keep <- ncbi_additional$locus
ncbi_seq_filtered <- ncbi_seq[names(ncbi_seq) %in% seq_names_to_keep]

# ---- Write filtered FASTA ----
fasta_output <- file.path(seq_out_dir, "150725_metadata_coverage_90_country_Philippines_sequences_filtered.fa")
write.fasta(
  sequences = ncbi_seq_filtered,
  names = names(ncbi_seq_filtered),
  file.out = fasta_output
)
# ============================================================
# Optional summaries
# ============================================================
message("Filtered NCBI metadata records: ", nrow(ncbi_additional))
message("Filtered sequences: ", length(ncbi_seq_filtered))
