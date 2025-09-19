# -----------------------------
# Libraries
# -----------------------------
library(dplyr)
library(lubridate)
library(gsheet)
library(tibble)
library(stringr)
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
  Barangay = character(),
  Municipality = character(),
  Province = character(),
  Region = character(),
  Latitude = numeric(),
  Longitude = numeric(),
  Source = character(),
  Genome_coverage = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------
# Helper function: map, clean, align
# -----------------------------
map_and_clean <- function(df, col_map, template = phylo_meta) {
  
  # Keep only columns that exist in df
  existing_cols <- names(col_map)[names(col_map) %in% names(df)]
  col_map_filtered <- col_map[existing_cols]
  
  # Select and rename
  df_clean <- df %>%
    select(all_of(names(col_map_filtered))) %>%
    rename_with(~ col_map_filtered[.x], .cols = everything())
  
  # Parse and format date as dd-MMM-yyyy
  if ("Preferred_date" %in% names(df_clean)) {
    df_clean <- df_clean %>%
      mutate(
        Preferred_date = parse_date_time(
          Preferred_date,
          orders = c("dmy", "d-b-y", "d-B-y", "mdy", "ymd", "Y"), 
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
  rename(Host = Source) %>%
  mutate(Source = "speedier", Region = "MIMIROPA") %>%
  mutate(`% coverage (nonMasked)` = `% coverage (nonMasked)`/ 100) %>%
  filter(`% coverage (nonMasked)` >= 0.9)


# vgtk (NCBI metadata)
vgtk <- read.csv("processed_data/processed_metadata/20250917_filtered_ncbi.csv") %>%
  mutate(Source = "vgtk")

# Essel/REDCap
mydata <- read.csv("raw_data/gathered_epi_metadata/ph_redcap_2024.v1.csv") %>%
  mutate(Source = "phd") %>%
  mutate(Isolate.ID = str_replace(Isolate.ID, "H-23-011Sk12", "H-23-011Sk_12")) #corrects common typo

# Zhang 2025 paper
zhang <- read.csv("raw_data/gathered_epi_metadata/zhang2025/Supplementary Table 3.csv")

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
    rename(Accession = primary_accession)
  
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
  "coverage" = "Genome_coverage"
)

mydata_col_map <- c(
  "Isolate.ID" = "Sample_ID",
  "isolation_source" = "Sample_type",
  "brgy" = "Barangay",
  "municipality" = "Municipality",
  "province" = "Province",
  "region" = "Region",
  "species" = "Host",
  "date" = "Preferred_date",
  "Source" = "Source",
  "latitude" = "Latitude",
  "longitude" = "Longitude",
  "Accession" = "Accession"
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
  map_and_clean(mydata, mydata_col_map)
)

# -----------------------------
# Check result
# -----------------------------
head(phylo_meta)
dim(phylo_meta)

# -----------------------------
# Check for duplicate Sample_id in phylo_meta
# -----------------------------

# Find duplicates
phylo_meta %>%
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
# Write results to file
# -----------------------------
# Create timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
records <- nrow(phylo_meta)

# Build filename with timestamp
outfile <- paste0("processed_data/processed_metadata/gathered_metadata_n", records,"_",timestamp, ".csv")

# Write file
write.csv(phylo_meta, outfile, row.names = FALSE)
