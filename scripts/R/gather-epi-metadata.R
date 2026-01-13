# -----------------------------
# Libraries
# -----------------------------
library(dplyr)
library(lubridate)
library(gsheet)
library(tibble)
library(stringr)
library(purrr)

# -----------------------------
# Create empty phylo_meta template
# -----------------------------
phylo_meta <- data.frame(
  Sample_ID = character(),
  Case_no = character(),
  Accession = character(),
  Host = character(),
  Sample_type = character(),
  Preferred_date = character(),  # store as day-month-year string
  Date_source = character(),
  Barangay = character(),
  Municipality = character(),
  Province = character(),
  Region = character(),
  Latitude = numeric(),
  Longitude = numeric(),
  Source = character(),
  Author = character(),
  Pubmed= numeric(),
  Genome_coverage = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------
# Create directory (if doesn't already exist) for outputs
# -----------------------------
dir.create("processed_data/processed_metadata/gathered_metadata/intermediate", showWarnings = FALSE, recursive = TRUE)
dir.create("processed_data/processed_metadata/gathered_metadata/final", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper function: map, clean, align
# -----------------------------
map_and_clean <- function(df, col_map, template = phylo_meta) {
  
  # Keep only columns that exist in df
  existing_cols <- names(col_map)[names(col_map) %in% names(df)]
  col_map_filtered <- col_map[existing_cols]
  
  # Select and dplyr::rename
  df_clean <- df %>%
    select(all_of(names(col_map_filtered))) %>%
    dplyr::rename_with(~ col_map_filtered[.x], .cols = everything())
  
  # Parse and format date as dd-MMM-yyyy
  if ("Preferred_date" %in% names(df_clean)) {
    df_clean <- df_clean %>%
      mutate(
        Preferred_date = parse_date_time(
          Preferred_date,
          orders = c("d-b-y", "d b Y", "d-b-Y", "mdy", "ymd", "Y"), 
          locale = "C",
          quiet = TRUE
        ),
        Preferred_date = format(Preferred_date, "%d-%b-%Y")
      )
  }
  
  # Add missing columns as NA
  missing_cols <- setdiff(names(template), names(df_clean))
  if (length(missing_cols) > 0) {
    df_clean <- bind_cols(
      df_clean,
      as_tibble(setNames(replicate(length(missing_cols), NA, simplify = FALSE), missing_cols))
    )
  }
  
  # Reorder columns to match template
  df_clean <- df_clean %>% select(all_of(names(template)))
  
  return(df_clean)
}

# -----------------------------
# Load datasets
# -----------------------------
# Speedier (Google Sheet)
genomics_link <- "https://docs.google.com/spreadsheets/d/1o9Ykf__3YTs33tqczZahjwmcnOEve90uw-AkduzDWXU/edit?gid=469362849#gid=469362849"
speedier <- gsheet2tbl(genomics_link) %>%
  dplyr::rename(Host = Source) %>%
  # Remove duplicated sample first
  filter(sample_id != "RADDL4B-24-184") %>%
  mutate(Source = "speedier", Region = "MIMIROPA",
  # Correct typos in sample IDs
  sample_id = case_when(
    sample_id == "RADDL4B-24-184-C" ~ "RADDL4B-24-184",
    sample_id == "RADDL4B-24-162-D" ~ "RADDL4B-24-162",
    TRUE ~ sample_id
  )) %>%
  mutate(`% coverage (nonMasked)` = `% coverage (nonMasked)`/ 100) %>%
  filter(`% coverage (nonMasked)` >= 0.9) %>%
  mutate(Date_source = "date_collected") 

# vgtk (NCBI metadata)
vgtk <- read.csv("processed_data/processed_metadata/150725_metadata_coverage_90_country_Philippines_filtered.csv") %>%
  mutate(Source = "vgtk")%>%
  mutate(Date_source = "date_collected")

# Essel/REDCap
mydata <- read.csv("raw_data/gathered_epi_metadata/ph_redcap_2024.v1.csv") %>%
  mutate(Source = "phd") %>%
  mutate(Isolate.ID = str_replace(Isolate.ID, "H-23-011Sk12", "H-23-011Sk_12")) %>% #corrects common typo
mutate(Date_source = "date_collected")

# Zhang 2025 paper
zhang <- read.csv("raw_data/gathered_epi_metadata/zhang2025/Supplementary Table 3.csv")%>%
  mutate(Date_source = "date_collected")

# 2018 workshop trip data
workshop <- read.csv("raw_data/gathered_epi_metadata/2018_workshop/2018_sequenced_collated_epi.csv")%>%
  mutate(Date_source = "date_collected")%>%
  mutate(across(
    c(lab_date, Date_sorted, Date_tested),
    ~ format(parse_date_time(., orders = c("mdy", "dmy")), "%d-%b-%Y"),
    .names = "{.col}_std"
  ))%>%
  mutate(
    Preferred_date = coalesce(lab_date_std, Date_sorted_std, Date_tested_std),
    Date_source = case_when(
      !is.na(lab_date_std) ~ "lab_date",
      !is.na(Date_sorted_std) ~ "Date_sorted",
      !is.na(Date_tested_std) ~ "Date_tested",
      TRUE ~ NA_character_
    )
  )%>%
  mutate(Source="2018_workshop")

# -----------------------------
# Check for duplicates of mydata in vgtk
# -----------------------------
# Ensure both datasets have a common key for matching
common_key <- intersect(mydata$Isolate.ID, vgtk$isolate)

if(length(common_key) > 0) {
  
  # For matched samples, add Accession from vgtk to mydata
  mydata <- mydata %>%
    left_join(
      vgtk %>%
        select(isolate, primary_accession),
      by = c("Isolate.ID" = "isolate")
    ) %>%
    dplyr::rename(Accession = primary_accession)
  
  # Remove the duplicates from vgtk
  vgtk <- vgtk %>%
    filter(!isolate %in% common_key)
}
# -----------------------------
# Check for duplicates of workshop in vgtk
# -----------------------------

# Ensure both datasets have a common key for matching
common_key <- intersect(workshop$sample_id, vgtk$isolate)

if (length(common_key) > 0) {
  
  # For matched samples, add Accession from vgtk to workshop
  workshop <- workshop %>%
    left_join(
      vgtk %>%
        select(isolate, primary_accession),
      by = c("sample_id" = "isolate")
    ) %>%
    dplyr::rename(Accession = primary_accession)
  
  # Remove the duplicates from vgtk
  vgtk <- vgtk %>%
    filter(!isolate %in% common_key)
}

# -----------------------------
# Check for duplicates of speedier in myphd
# -----------------------------

# Identify common sample IDs
common_key <- intersect(mydata$Isolate.ID, speedier$sample_id)

if(length(common_key) > 0) {
  
  # Add any useful metadata from speedier into mydata
  mydata <- mydata %>%
    left_join(
      speedier %>%
        select(sample_id, case_number, SAMPLE_TYPE, Date_collected),
      by = c("Isolate.ID" = "sample_id")
    ) %>%
    mutate(
      # Keep mydata values if present, otherwise use speedier’s
      Case_no = coalesce(Case_no, case_number),
      Sample_type = coalesce(Sample_type, SAMPLE_TYPE),
      Preferred_date = coalesce(Preferred_date, Date_collected)
    ) %>%
    select(-case_number, -SAMPLE_TYPE, -Date_collected)
  
  # Remove duplicates from speedier
  speedier <- speedier %>%
    filter(!sample_id %in% common_key)
}

# -----------------------------
# Check for duplicates of speedier in workshop
# -----------------------------

# Identify common sample IDs
common_key <- intersect(workshop$sample_id, speedier$sample_id)

if (length(common_key) > 0) {
  
  # Add any useful metadata from speedier into workshop
  workshop <- workshop %>%
    left_join(
      speedier %>%
        select(sample_id, case_number, SAMPLE_TYPE, Date_collected),
      by = "sample_id"
    ) %>%
    mutate(
      # Keep workshop values if present, otherwise use speedier’s
      Case_no        = coalesce(Case_no, case_number),
      Sample_type    = coalesce(Sample_type, SAMPLE_TYPE),
      Preferred_date = coalesce(Preferred_date, Date_collected)
    ) %>%
    select(-case_number, -SAMPLE_TYPE, -Date_collected)
  
  # Remove duplicates from speedier
  speedier <- speedier %>%
    filter(!sample_id %in% common_key)
}

# -----------------------------
# Check for duplicates of workshop in mydata (phd)
# -----------------------------

# Identify common sample IDs
common_key <- intersect(workshop$sample_id, mydata$Isolate.ID)

if (length(common_key) > 0) {
  
  # Add any useful metadata from workshop into mydata
  mydata <- mydata %>%
    left_join(
      workshop %>% select(sample_id, Preferred_date),
      by = c("Isolate.ID" = "sample_id")
    ) %>%
    mutate(
      Preferred_date = coalesce(date, Preferred_date)  # use mydata$date first
    ) %>%
    select(-date)  # drop original date column if no longer needed
  
  # Remove duplicates from workshop
  workshop <- workshop %>%
    filter(!sample_id %in% common_key)
}
# -----------------------------
# Define column maps
# -----------------------------
speedier_col_map <- c(
  "sample_id" = "Sample_ID",
  "case_number" = "Case_no",
  "SAMPLE_TYPE" = "Sample_type",
  "Province" = "Province",
  "Municipality" = "Municipality",
  "Barangay" = "Barangay",
  "Region" = "Region",
  "Latitude" = "Latitude",
  "Longitude" = "Longitude",
  "Host" = "Host",
  "Date_collected" = "Preferred_date",
  "Date_source" = "Date_source",
  "Source" = "Source",
  "% coverage (nonMasked)" = "Genome_coverage"
)

vgtk_col_map <- c(
  "isolate" = "Sample_ID",
  "isolation_source" = "Sample_type",
  "geo_loc" = "Province",
  "host" = "Host",
  "collection_date" = "Preferred_date",
  "Source" = "Source",
  "primary_accession" = "Accession",
  "coverage" = "Genome_coverage", 
  "authors" = "Author",
  "pubmed_id" = "Pubmed"
)

mydata_col_map <- c(
  "Isolate.ID" = "Sample_ID",
  "isolation_source" = "Sample_type",
  "brgy" = "Barangay",
  "municipality" = "Municipality",
  "province" = "Province",
  "region" = "Region",
  "species" = "Host",
  "Preferred_date" = "Preferred_date",
  "Source" = "Source",
  "latitude" = "Latitude",
  "longitude" = "Longitude",
  "Accession" = "Accession"
)

workshop_col_map <- c(
  "sample_id" = "Sample_ID",
  "Barangay" = "Barangay",
  "Municipality" = "Municipality",
  "Province" = "Province",
  "Region" = "Region",
  "Specimen" = "Host",
  "Preferred_date" = "Preferred_date",
  "Source" = "Source"
)

# -----------------------------
# Apply Zhang updates to vgtk
# -----------------------------
zhang_in_vgtk <- zhang[zhang$Acceccsion.No. %in% vgtk$primary_accession, ]
vgtk$geo_loc <- ifelse(
  vgtk$primary_accession %in% zhang_in_vgtk$Acceccsion.No.,
  zhang$Location[match(vgtk$primary_accession, zhang$Acceccsion.No.)],
  vgtk$geo_loc
)

# -----------------------------
# Map and bind datasets into phylo_meta
# -----------------------------
phylo_meta <- bind_rows(
  phylo_meta,
  map_and_clean(speedier, speedier_col_map),
  map_and_clean(vgtk, vgtk_col_map),
  map_and_clean(mydata, mydata_col_map),
  map_and_clean(workshop,workshop_col_map)
)

# Strip whitespace from all character columns
phylo_meta  <- phylo_meta %>%
  mutate(across(where(is.character), ~ str_trim(.)))

# -----------------------------
# Check result
# -----------------------------
head(phylo_meta)
dim(phylo_meta)

# -----------------------------
# Check for duplicate Sample_id in phylo_meta
# -----------------------------

# Find duplicates
dup_ids <- phylo_meta %>%
  group_by(Sample_ID = .data[["Sample_ID"]]) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(Sample_ID)

# Extract rows with duplicates for inspection
dup_rows <- phylo_meta %>%
  filter(Sample_ID %in% dup_ids) %>%
  arrange(Sample_ID)

# View results
dup_ids      # just the duplicated IDs
dup_rows     # full rows to check

# -----------------------------
# Clean and apply RADDL corrections (Zhang updates)
# -----------------------------
raddl_link <- "https://docs.google.com/spreadsheets/d/1h9fdNvgfcmA02QDT3M9drNqsC2Zr9RlCYnklGqhiR4o/edit?usp=sharing"
raddl_corrections <- gsheet2tbl(raddl_link)

# Clean column names and trim whitespace
raddl_corrections <- raddl_corrections %>%
  dplyr::rename_with(~ str_trim(.x)) %>%
  mutate(across(everything(), ~ str_trim(as.character(.x)))) %>%
  dplyr::rename(
    Sample_ID = "Sample_id",
    Accession = "Accession No.",
    Barangay = "raddl_Brgy",
    Municipality = "raddl_Municipality",
    Province = "raddl_location",
    Preferred_date = "raddl_date",
    Host = "Host"
  )

# Standardize date format
raddl_corrections <- raddl_corrections %>%
  mutate(
    Preferred_date = parse_date_time(Preferred_date, orders = c("d-b-y","d-B-y"), quiet = TRUE),
    Preferred_date = format(Preferred_date, "%d-%b-%Y")) %>%
      mutate(across(where(is.character), ~ str_to_sentence(.x)))

# -----------------------------
# Columns to update
# -----------------------------
key_cols <- c("Preferred_date", "Barangay", "Municipality", "Province", "Host")

# Make a copy of original
phylo_meta_corrected <- phylo_meta

# -----------------------------
# Apply corrections safely column by column
# -----------------------------
for(col in key_cols){
  # Get matching new values
  new_vals <- raddl_corrections[[col]][match(phylo_meta_corrected$Sample_ID, raddl_corrections$Sample_ID)]
  
  # Replace only if new value is not NA or blank
  phylo_meta_corrected[[col]] <- ifelse(
    !is.na(new_vals) & new_vals != "",
    new_vals,
    phylo_meta_corrected[[col]]
  )
}

# -----------------------------
# Update Source ONLY for rows where at least one key column changed
# -----------------------------
changed_rows <- sapply(seq_len(nrow(phylo_meta)), function(i){
  sid <- phylo_meta$Sample_ID[i]
  if(sid %in% raddl_corrections$Sample_ID){
    orig_vals <- phylo_meta[i, key_cols]
    new_vals  <- phylo_meta_corrected[i, key_cols]
    
    any(mapply(function(ov, nv){
      # Treat NA/blank → value as a change
      is.na(ov) && !is.na(nv) ||
        ov == "" && nv != "" ||
        (!is.na(ov) && ov != "" && ov != nv)
    }, orig_vals, new_vals))
    
  } else {
    FALSE
  }
})

# Ensure changed_rows is a logical vector with no NAs
changed_rows_safe <- ifelse(is.na(changed_rows), FALSE, changed_rows)

# Append "raddl_correction" to existing Source where changes occurred
phylo_meta_corrected$Source[changed_rows_safe] <- phylo_meta_corrected$Source[changed_rows_safe] %>%
  ifelse(
    grepl("raddl_correction", ., ignore.case = TRUE),
    .,  # leave as-is if already marked
    paste0(., ifelse(. == "" | is.na(.), "", "; "), "raddl_correction")
  )
# -----------------------------
# Clean and apply manual corrections (to Bacus data)
# -----------------------------
# manual extraction/standardisation of location data all deposited in province col

manual_corrections <- read.csv("processed_data/processed_metadata/LocationManualCorrections_05Jan26.csv")

# Make a copy of original
phylo_meta_corrected2 <- phylo_meta_corrected

key_cols2 <- c("Barangay", "Municipality", "Province", "Region")

# -----------------------------
# Apply corrections safely column by column
# -----------------------------
for(col in key_cols2){
  # Get matching new values
  new_vals <- manual_corrections[[col]][match(phylo_meta_corrected2$Sample_ID, manual_corrections$Sample_ID)]
  
  # Replace only if new value is not NA or blank
  phylo_meta_corrected2[[col]] <- ifelse(
    !is.na(new_vals) & new_vals != "",
    new_vals,
    phylo_meta_corrected2[[col]]
  )
}

# -----------------------------
# Update Source ONLY for rows where key_cols2 changed
# -----------------------------
changed_rows <- sapply(seq_len(nrow(phylo_meta_corrected2)), function(i) {
  sid <- phylo_meta_corrected2$Sample_ID[i]
  if (sid %in% manual_corrections$Sample_ID) {
    # Compare the same row between pre- and post-manual correction for key_cols2 only
    orig_vals <- phylo_meta_corrected[i, key_cols2, drop = FALSE]
    new_vals  <- phylo_meta_corrected2[i, key_cols2, drop = FALSE]
    
    # TRUE if any of these three columns actually changed
    any(mapply(function(ov, nv) {
      (is.na(ov) && !is.na(nv)) ||
        (ov == "" && nv != "") ||
        (!is.na(ov) && ov != "" && ov != nv)
    }, orig_vals, new_vals))
  } else {
    FALSE
  }
})


# Ensure changed_rows is logical and no NAs
changed_rows_safe <- ifelse(is.na(changed_rows), FALSE, changed_rows)

# Append "manual_correction" to Source where changes occurred
phylo_meta_corrected2$Source[changed_rows_safe] <- phylo_meta_corrected2$Source[changed_rows_safe] %>%
  ifelse(
    grepl("manual_correction", ., ignore.case = TRUE),
    .,  # leave as-is if already marked
    paste0(., ifelse(. == "" | is.na(.), "", "; "), "manual_correction")
  )

# -----------------------------
# Final fixes
# -----------------------------

# -----------------------------
# Manual date entries for Zhang paper (raddl cases) based on sample ids 
# -----------------------------

# Update Preferred_date values
phylo_meta_corrected2$Preferred_date[phylo_meta_corrected2$Sample_ID == "R2-2019-7465"] <- "01-Jan-2019"
phylo_meta_corrected2$Preferred_date[phylo_meta_corrected2$Sample_ID == "R4B-2022-00011"] <- "01-Jan-2022"

# Identify rows you just changed
changed_rows <- phylo_meta_corrected2$Sample_ID %in% c("R2-2019-7465", "R4B-2022-00011")

# Append "manual_correction" to Source for those rows
phylo_meta_corrected2$Source[changed_rows] <- phylo_meta_corrected2$Source[changed_rows] %>%
  ifelse(
    grepl("manual_correction", ., ignore.case = TRUE),
    .,  # leave as-is if already marked
    paste0(., ifelse(. == "" | is.na(.), "", "; "), "manual_correction")
  )

# -----------------------------
# Convert Preferred_date to Excel-safe text (no visible tick)
# -----------------------------
# Manually fix Preferred_date for a single sample
phylo_meta_corrected2$Preferred_date[phylo_meta_corrected2$Sample_ID == "RADDL4B-23-092"] <- "07-Dec-2023"
for (df_name in c("phylo_meta", "phylo_meta_corrected", "phylo_meta_corrected2")) {
  df <- get(df_name)
  
  if ("Preferred_date" %in% names(df)) {
    df <- df %>%
      mutate(
        Preferred_date = ifelse(
          !is.na(Preferred_date) & Preferred_date != "",
          paste0("=\"", Preferred_date, "\""),
          Preferred_date
        )
      )
  }
  
  assign(df_name, df)
}

# -----------------------------
# Fill in date source col
# -----------------------------

# Assume date collected for all non-speedier sequences
phylo_meta_corrected2 <- phylo_meta_corrected2 %>%
  mutate(
    Date_source = if_else(is.na(Date_source), "date_collected", Date_source)
  )



# -----------------------------
# Standardise host column - add cols for scientific, common and type host
# -----------------------------
# Standardise host data
source("scripts/R/standardise_hosts.R")

# -----------------------------
# Write results to file
# -----------------------------
# Create timestamp
timestamp <- format(Sys.time(), "%d%b%y")
records <- nrow(phylo_meta)

# Build filename with timestamp
outfile <- paste0("processed_data/processed_metadata/gathered_metadata/intermediate/",timestamp,"_gathered_metadata_n", records, ".csv")
outfile2 <- paste0("processed_data/processed_metadata/gathered_metadata/intermediate/",timestamp,"_gathered_metadata_n",records,"_raddlCorrected.csv")
outfile3 <- paste0("processed_data/processed_metadata/gathered_metadata/final/",timestamp,"_gathered_metadata_n", records,"_raddl_and_manual_Corrected.csv")
# Write file
write.csv(phylo_meta, outfile, row.names = FALSE)
write.csv(phylo_meta_corrected, outfile2, row.names = FALSE)
write.csv(phylo_meta_corrected2, outfile3, row.names = FALSE)
