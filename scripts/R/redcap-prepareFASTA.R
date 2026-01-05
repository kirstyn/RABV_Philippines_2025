library(dplyr)
library(stringr)

# ---- SETTINGS ----
metadata_path <- "/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/redcap/rabv_4b_2025_r10_run10-11/2025_run11/rabv_4b_2025_r10_run10_sequencing_form.csv" # CSV with columns: sampleID, ngs_runid, redcap_repeat_instance
fasta_dir <- "/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/redcap/rabv_4b_2025_r10_run10-11/2025_run11/split_fasta/"           # Folder containing your FASTA files
output_dir <- "/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/redcap/rabv_4b_2025_r10_run10-11/2025_run11/split_fasta_renamed/"  # Output folder for renamed FASTAs

# Create output dir if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)

# ---- LOAD METADATA ----
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Ensure columns are present
required_cols <- c("sample_id", "ngs_runid", "redcap_repeat_instance")
if (!all(required_cols %in% names(meta))) {
  stop("Metadata must contain: sample_id, ngs_runid, redcap_repeat_instance")
}

# ---- PROCESS FASTA FILES ----
fasta_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

for (file_path in fasta_files) {
  
  # Extract sample ID from filename
  this_sample_id <- str_remove(basename(file_path), "\\.fasta$")
  
  # Look up metadata
  row <- meta %>% filter(sample_id == this_sample_id)
  
  if (nrow(row) == 0) {
    message("⚠️ No metadata found for: ", this_sample_id, " — skipping")
    next
  }
  
  run_id  <- row$ngs_runid[1]
  instance <- row$redcap_repeat_instance[1]
  
  # New filename (FIXED)
  new_name <- sprintf(
    "%s__%s__instance%s.fasta",
    this_sample_id,
    run_id,
    instance
  )
  
  new_path <- file.path(output_dir, new_name)
  
  # Read FASTA
  fasta_lines <- readLines(file_path)
  
  # Update FASTA header
  if (length(fasta_lines) > 0 && startsWith(fasta_lines[1], ">")) {
    fasta_lines[1] <- paste0(
      ">",
      this_sample_id,
      "__",
      run_id,
      "__instance",
      instance
    )
  }
  
  # Write renamed FASTA
  writeLines(fasta_lines, new_path)
  
  message("✅ Renamed ", basename(file_path), " → ", new_name)
}

message("🎯 Done! Processed ", length(fasta_files), " files.")
