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
meta_file <- "processed_data/processed_metadata/gathered_metadata/final/11Nov25_gathered_metadata_n811_raddl_and_manual_Corrected_stdGeo.csv"
clade_file <- "processed_data/processed_metadata/241025_glue_clade_assignment.txt"
lineage_file <- "analysis/MADDOG/Outputs/sequence_data.csv"

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
clade <- read.table(clade_file, header=T, sep="\t")
clade <- clade %>%
  mutate(across(where(is.character), ~ trimws(.)))
lineage <- read.csv(lineage_file,stringsAsFactors = FALSE)

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

# Add clade info to metadata # 

meta_filtered <- meta_filtered %>%
  left_join(
    clade %>% 
      select(queryName, phylogenetic_clade = minor_cladeFinalClade),
    by = c("seq_name" = "queryName")
  )

cat("Added clade information for", sum(!is.na(meta_filtered$major_clade)), "records.\n\n")

# Add lineage info to metadata #

meta_filtered <- meta_filtered %>%
  left_join(
    lineage %>% 
      select(ID, lineage),
    by = c("seq_name" = "ID")
  )


# Add decimal date col to metadata # 

meta_filtered <- meta_filtered %>%
  mutate(
    decimal_date = Preferred_date %>%
      gsub('="|"', "", .) %>%   # remove =" "
      dmy() %>%                 # parse date
      decimal_date()            # convert to decimal
  )

cat("Added decimal date for", sum(!is.na(meta_filtered$Preferred_date)), "records.\n\n")

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
#write.table(meta_filtered, matched_meta_path, row.names = FALSE, sep="\t")
cat("Matched metadata saved to:\n", matched_meta_path, "\n\n")


## --- and filter metadata to ones just with dates ---

meta_with_dates <- meta_filtered %>%
  filter(!is.na(decimal_date))

cat("Sequences retained with valid dates:", nrow(meta_with_dates ), "\n\n")
names(seq) <- seq_names  
seq_filtered <- seq[names(seq) %in% meta_with_dates$seq_name]

cat("Sequences retained after filtering by valid dates:", length(seq_filtered), "of", length(seq), "\n")

# ---- Save filtered metadata ----
# Use the correct object (meta_with_dates) and safe file name counts
filtered_meta_path <- file.path(
  out_meta_dir,
  paste0("PHL_metadata_matchedWithDates_", Sys.Date(), "_n", nrow(meta_with_dates), ".csv")
)
write.csv(meta_with_dates, filtered_meta_path, row.names = FALSE)
cat("Filtered metadata (with valid dates) saved to:\n", filtered_meta_path, "\n\n")

# ---- Save filtered FASTA (unaltered names) ----
fasta_no_dates_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDates_", Sys.Date(), "_n", length(seq_filtered), ".fasta")
)
# Ensure seqinr::write.fasta gets sequences as a list of character vectors and corresponding names
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_no_dates_path)
cat("Filtered FASTA saved (unaltered names):\n", fasta_no_dates_path, "\n\n")

# ----  Append dates to sequence names (matched by seq_name, safe ordering + sanitise) ----
# helper to sanitise strings for FASTA headers (remove problematic characters, trim)
sanitise_header <- function(x) {
  x2 <- ifelse(is.na(x), "", x)
  x2 <- trimws(x2)
  
  # Replace problematic punctuation except the hyphen used in dates
  x2 <- gsub("[/:;,+\\\\]", "_", x2)       # replace unsafe punctuation
  x2 <- gsub("\\s+", "_", x2)             # collapse spaces to underscore
  
  # Preserve legitimate date hyphen format (e.g., 18-Jun-1998 or 18-Jun-98)
  # We'll temporarily protect hyphens that match date patterns:
  x2 <- gsub("(\\d{1,2})-([A-Za-z]{3})-(\\d{2,4})", "\\1§\\2§\\3", x2)
  
  # Now strip all other unsafe characters but leave hyphen intact
  x2 <- gsub("[^A-Za-z0-9_\\-\\.§]", "", x2)  # now hyphen is allowed
  
  # Restore protected date hyphens
  x2 <- gsub("§", "-", x2)
  
  ifelse(nchar(x2) == 0, NA_character_, x2)
}

# Build lookups (named vectors) keyed by seq_name
date_lookup   <- setNames(meta_with_dates$Preferred_date, meta_with_dates$seq_name)
region_lookup <- setNames(meta_with_dates$Region_std,    meta_with_dates$seq_name)

# original sequence names (these are the IDs that matched)
orig_seq_names <- names(seq_filtered)

# create cleaned date and region strings aligned to seq_filtered order
clean_dates_vec  <- sanitise_header(date_lookup[orig_seq_names])
clean_regions_vec <- sanitise_header(region_lookup[orig_seq_names])

# create new names (fall back to original ID if date missing)
new_names_date <- ifelse(
  !is.na(clean_dates_vec),
  paste0(orig_seq_names, "__", clean_dates_vec),
  orig_seq_names
)
names(seq_filtered) <- new_names_date

fasta_with_dates_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDatesLabeled_", Sys.Date(), "_n", length(seq_filtered), ".fasta")
)
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_with_dates_path)
cat("FASTA with dates appended to sequence names saved to:\n", fasta_with_dates_path, "\n\n")

# ----  Append dates and region to sequence names (maintaining correct order) ----
# Use the original IDs again (so we can re-build deterministically)
orig_seq_names <- names(seq_filtered)
# If names already have dates appended (they do now), extract base seq_name (split on "__")
base_seq_names <- sub("__.*$", "", orig_seq_names)

# build final label: seqname_date_region (skip parts that are NA)
final_labels <- vapply(seq_along(base_seq_names), function(i) {
  id <- base_seq_names[i]
  d  <- clean_dates_vec[id]
  r  <- clean_regions_vec[id]
  parts <- c(id)
  if (!is.na(d)) parts <- c(parts, d)
  if (!is.na(r)) parts <- c(parts, r)
  paste(parts, collapse = "_")
}, FUN.VALUE = character(1), USE.NAMES = FALSE)

names(seq_filtered) <- final_labels

fasta_with_dates_regions_path <- file.path(
  out_seq_dir,
  paste0("PHL_sequences_withDatesRegionsLabeled_", Sys.Date(), "_n", length(seq_filtered), ".fasta")
)
write.fasta(sequences = seq_filtered, names = names(seq_filtered), file.out = fasta_with_dates_regions_path)
cat("FASTA with dates and regions appended to sequence names saved to:\n", fasta_with_dates_regions_path, "\n\n")


head(data.frame(original = names(seq)[1:5], new = names(seq_filtered)[1:5]))
# ---- 13. Summary (fixed variable names) ----
cat("\n================ SUMMARY ================\n")
cat("Total sequences:", length(seq_names), "\n")
cat("Sequences matched to metadata:", sum(matched_any), "\n")
cat("Sequences with valid dates (exported):", length(seq_filtered), "\n")
cat("Final metadata records exported (with dates):", nrow(meta_with_dates), "\n")
cat("=========================================\n")
