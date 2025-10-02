## Pull latest MIMIROPA seq and phylogenetics-relevant metadata

# redcap link
source("scripts/R/redcap-api.R")
source("scripts/R/redcap-functions.r")

# Get sequence data
# Export all records including repeat instances filtered by DAG
df <- exportRecordsTyped(
  rcon,
  dag = TRUE
)

# Merge the 2 forms into one row
df_merged <- redcap_form_merge(rcon, df, form1 = "diagnostic", form2 = "sequencing")

# Filter to just Philippines data (only necessary if you have access to more than one dag)
df_phl <- df_merged %>%
  dplyr::filter(country == "PHL: Philippines") %>%
  mutate(
    consensus_coverage = as.numeric(as.character(consensus_coverage)),  # convert to numeric
    consensus_coverage = consensus_coverage / 11932 * 100
  ) %>%
  filter(consensus_coverage>=90)

# Get metadata and filter to mimiropa
all <- read.csv("processed_data/processed_metadata/gathered_metadata_n794_20250922_155050.csv")
mimiropa_outbreak=all %>%
  filter(Region %in% c("MIMIROPA", "Southwestern Tagalog Region")) %>%
  mutate(Preferred_date = as.Date(Preferred_date, format = "%d-%b-%Y")) %>%
  filter(Preferred_date >= as.Date("2022-01-01")) %>%
  filter(Province %in% c("Marinduque","Romblon")) 
## Note: DON'T filter by genome cov because epi data only has this info if pulled from vgtk or mydata. Can filter using seq metadata instead. 

# Correct some common typos that are still present
mimiropa_outbreak <- mimiropa_outbreak %>%
  mutate(Sample_ID = if_else(
    Sample_ID == "RADDL4B-24-184-C",
    "RADDL4B-24-184",
    Sample_ID
  ))

# Join mimiropa_outbreak to df_phl, keeping only matches and adding df_phl columns
mimiropa_matched <- mimiropa_outbreak %>%
  inner_join(
    df_phl,
    by = c("Sample_ID" = "sample_id")  # adjust column names if needed
  )

# Rows in mimiropa_outbreak with no match in df_phl
mimiropa_unmatched <- mimiropa_outbreak %>%
  anti_join(
    df_phl,
    by = c("Sample_ID" = "sample_id")  # adjust if column names differ
  )

# # Download the associated fasta files

# Create an output folder named with the current date/time

output_dir <- file.path("raw_data/redcap_sequences_and_metadata/", paste0("redcap_download_", format(Sys.time(), "%Y%m%d_%H%M")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Write the associated metadata there too
write.csv(mimiropa_matched, paste0(output_dir, format(Sys.time(), "%Y%m%d_%H%M"),"redcap_meta_mimiropa.csv"), row.names=F)

# Export <-exportFilesMultiple(
#   rcon,
#   record=record_ids,
#   field="consensus_fasta",
#   repeat_instance = instances,
#   dir=output_dir,
#   file_prefix = F,
# )


# record IDs and repeat instances
record_ids <- mimiropa_matched$Sample_ID
instances <- as.integer(mimiropa_matched$redcap_repeat_instance)

total_records <- length(record_ids)
successful_downloads <- 0
failed_downloads <- 0
failed_records <- character()  # to store IDs with failures or no files

# Loop through both record IDs and their corresponding instances
for (i in seq_along(record_ids)) {
  rec <- record_ids[i]
  inst <- instances[i]

  message("Downloading files for record: ", rec, " (instance ", inst, ")")

  tryCatch({
    files_downloaded <- exportFiles(
      rcon,
      record = rec,
      field = "consensus_fasta",
      repeat_instance = inst,
      dir = output_dir
    )

    if (length(files_downloaded) > 0) {
      successful_downloads <- successful_downloads + 1
    } else {
      failed_downloads <- failed_downloads + 1
      failed_records <- c(failed_records, rec)
    }

  }, error = function(e) {
    message("Error downloading for record: ", rec, " - ", e$message)
    failed_downloads <- failed_downloads + 1
    failed_records <- c(failed_records, rec)
  })
}

message("✅ Downloads complete: ", successful_downloads, "/", total_records)
if (failed_downloads > 0) {
  message("⚠️ Failed downloads: ", failed_downloads)
  print(failed_records)
}

message("All files saved in: ", output_dir)



