library(redcapAPI)
library(dplyr)

exists("rcon")
# REDCap connection
exists("rcon")
#If this doesn't exist then need to run this to set up:
unlockREDCap(c(rcon= 'rage-redcap'),
             keyring= "login",
             envir= globalenv(),
             url= 'https://cvr-redcap.mvls.gla.ac.uk/redcap/redcap_v15.5.13/API/')


# Export all records including repeat instances filtered by DAG
df <- exportRecordsTyped(
  rcon
)

# Diagnostic form rows (non-repeating)
df_diag <- df %>%
  filter(is.na(redcap_repeat_instrument) | redcap_repeat_instrument == "")

# Sequencing form rows (repeating)
df_seq <- df %>%
  filter(redcap_repeat_instrument == "Sequencing")

# Get column names belonging to each form from metadata
meta <- rcon$metadata()

diagnostic_fields <- meta %>%
  filter(form_name == "diagnostic") %>%
  pull(field_name)

sequencing_fields <- meta %>%
  filter(form_name == "sequencing") %>%
  pull(field_name)

# Select only diagnostic columns + ID columns from diagnostic df
df_diag_clean <- df_diag %>%
  select(sample_id, any_of(diagnostic_fields))

# Select only sequencing columns + ID columns from sequencing df
df_seq_clean <- df_seq %>%
  select(sample_id, redcap_repeat_instance, any_of(sequencing_fields))

# Now join by sample_id 
df_merged <- df_seq_clean %>%
  left_join(df_diag_clean, by = "sample_id")
# Note this will contain all sequencing instances so there may be multiple rows per sample

## df_merged is a tidied version of the entire dataset from redcap (including different countries, neg control, repeat instances etc)

# Remove negative controls
df_positives <- df_merged %>%
  filter(negative_control_sample != "Yes")

# filter by latest instance
df_latest <- df_positives %>%
  group_by(sample_id) %>%
  filter(redcap_repeat_instance == max(redcap_repeat_instance, na.rm = TRUE)) %>%
  ungroup()

# filter by country
df_phl <- df_latest %>%
  filter(grepl("Philippines", country, ignore.case = TRUE))
dim(df_phl)
#write.csv(df_phl,"~/Downloads/philippines_sequencing.csv", row.names=F)

# filter by genome coverage >=90%
df_phl_90 <- df_phl %>%
  filter(!is.na(consensus_coverage)) %>% 
  mutate(consensus_coverage = as.numeric(consensus_coverage)) %>%
  filter((consensus_coverage / 11932) * 100 >= 90)

# # # Download the associated fasta files
# 
# # Create an output folder named with the current date/time
# 
# output_dir <- file.path("processed_data/redcap_sequences_and_metadata/", paste0("redcap_download_", format(Sys.time(), "%Y%m%d_%H%M")))
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# # Write the associated metadata there too
# write.csv(df_phl_90, paste0(output_dir, format(Sys.time(), "%Y%m%d_%H%M"),"redcap_meta_phl.csv"), row.names=F)
# 
# # Export <-exportFilesMultiple(
# #   rcon,
# #   record=record_ids,
# #   field="consensus_fasta",
# #   repeat_instance = instances,
# #   dir=output_dir,
# #   file_prefix = F,
# # )
# 
# 
# # record IDs and repeat instances
# record_ids <- df_phl_90$sample_id
# instances <- as.integer(df_phl_90$redcap_repeat_instance)
# 
# total_records <- length(record_ids)
# successful_downloads <- 0
# failed_downloads <- 0
# failed_records <- character()  # to store IDs with failures or no files
# 
# # Loop through both record IDs and their corresponding instances
# for (i in seq_along(record_ids)) {
#   rec <- record_ids[i]
#   inst <- instances[i]
# 
#   message("Downloading files for record: ", rec, " (instance ", inst, ")")
# 
#   tryCatch({
#     files_downloaded <- exportFiles(
#       rcon,
#       record = rec,
#       field = "consensus_fasta",
#       repeat_instance = inst,
#       dir = output_dir
#     )
# 
#     if (length(files_downloaded) > 0) {
#       successful_downloads <- successful_downloads + 1
#     } else {
#       failed_downloads <- failed_downloads + 1
#       failed_records <- c(failed_records, rec)
#     }
# 
#   }, error = function(e) {
#     message("Error downloading for record: ", rec, " - ", e$message)
#     failed_downloads <- failed_downloads + 1
#     failed_records <- c(failed_records, rec)
#   })
# }
# 
# message("✅ Downloads complete: ", successful_downloads, "/", total_records)
# if (failed_downloads > 0) {
#   message("⚠️ Failed downloads: ", failed_downloads)
#   print(failed_records)
# }
# 
# message("All files saved in: ", output_dir)
