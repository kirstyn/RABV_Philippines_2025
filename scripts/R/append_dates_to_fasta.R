library(seqinr)
library(dplyr)
library(stringr)
library(lubridate)

# ---- 1. Load sequences and metadata ----
seq <- read.fasta("processed_data/processed_sequences/241025_PHL_all_n794.fasta", as.string = TRUE)
meta <- read.csv("processed_data/processed_metadata/gathered_metadata/23Oct25_gathered_metadata_n797_raddl_and_manual_Corrected_stdGeo.csv")

seq_names <- names(seq)

# ---- 2. Manual date fixes ----
manual_dates <- tribble(
  ~Sample_ID,       ~Preferred_date,
  "19-30",          "2019/01/01",
  "19-40",          "2019/01/01",
  "CAR-22-018",     "2018/01/01",
  "CAR-22-020",     "2020/01/01",
  "COW",            "2020/01/01",
  "GOAT",           "2020/01/01",
  "R11-22-001",     "2022/01/01",
  "R2-22-5302",     "2022/01/01",
  "R2-22-5814",     "2022/01/01",
  "R2-22-84711",    "2022/01/01",
  "R2-22-8924",     "2022/01/01",
  "R3-0483",        "2022/01/01",
  "R6-21-4343",     "2021/01/01",
  "R7-21-149",      "2021/01/01"
)

# ---- 3. Clean and standardize Preferred_date ----
meta_lookup <- meta %>%
  select(Sample_ID, Accession, Preferred_date, Region_std) %>%
  mutate(
    across(c(Sample_ID, Accession, Region_std), as.character),
    Preferred_date = Preferred_date %>%
      str_replace_all("^=\\\"|\\\"$", "") %>%
      str_replace_all("\"", "") %>%
      str_trim() %>%
      parse_date_time(orders = c("d-b-Y", "d-B-Y", "d-m-Y", "Y-m-d")) %>%
      format("%Y-%m-%d")
  )

# ---- 4. Add manual cases ----
meta_combined <- meta_lookup %>%
  bind_rows(manual_dates %>%
              mutate(
                Preferred_date = parse_date_time(Preferred_date, orders = "Y/m/d") %>%
                  format("%Y-%m-%d"),
                Region_std = NA_character_
              )
  )

# ---- 5. Match sequences to metadata ----
matched_df <- data.frame(seq_name = seq_names, stringsAsFactors = FALSE) %>%
  left_join(meta_combined, by = c("seq_name" = "Sample_ID")) %>%
  mutate(
    Preferred_date = if_else(
      is.na(Preferred_date) | Preferred_date == "NA",
      meta_combined$Preferred_date[match(seq_name, meta_combined$Accession)],
      Preferred_date
    ),
    Region_std = if_else(
      is.na(Region_std) | Region_std == "NA",
      meta_combined$Region_std[match(seq_name, meta_combined$Accession)],
      Region_std
    )
  )

# ---- 6. Filter out sequences with no date ----
has_date <- !is.na(matched_df$Preferred_date) & matched_df$Preferred_date != "NA"

seq_filtered <- seq[has_date]
matched_filtered <- matched_df[has_date, ]

write.fasta(
  sequences = seq_filtered,
  names = names(seq_filtered),
  file.out = "processed_data/processed_sequences/20251023_PHL_all_filtered.fasta"
)

# ---- 7. Append dates only ----
new_names_date <- paste0(matched_filtered$seq_name, "_", matched_filtered$Preferred_date)
names(seq_filtered) <- new_names_date

write.fasta(
  sequences = seq_filtered,
  names = names(seq_filtered),
  file.out = "processed_data/processed_sequences/20251023_PHL_all_withDates_YYYYMMDD_filtered.aln.fasta"
)

# ---- 8. Append dates + Region_std ----
new_names_date_region <- paste0(
  matched_filtered$seq_name, "_", matched_filtered$Preferred_date,
  ifelse(!is.na(matched_filtered$Region_std) & matched_filtered$Region_std != "",
         paste0("_", matched_filtered$Region_std),
         "")
)

# Replace spaces or commas in region names with underscores
new_names_date_region <- str_replace_all(new_names_date_region, "[,\\s]+", "_")
names(seq_filtered) <- new_names_date_region

write.fasta(
  sequences = seq_filtered,
  names = names(seq_filtered),
  file.out = "processed_data/processed_sequences/20251023_PHL_all_withDates_Region_filtered.aln.fasta"
)
)