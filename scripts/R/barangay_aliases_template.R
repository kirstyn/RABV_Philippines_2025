library(dplyr)
library(tidyr)
library(stringr)

# -------------------------------
# 1. Load data
# -------------------------------
meta <- read.csv(
  "processed_data/processed_metadata/gathered_metadata/final/05Jan26_gathered_metadata_n811_raddl_and_manual_Corrected_stdGeo.csv",
  stringsAsFactors = FALSE
)
barangay_lookup <- read.csv(
  "processed_data/gis_data/psgc_barangay_municipality_province_lookup_09Jan26.csv",
  stringsAsFactors = FALSE
)

# -------------------------------
# 2. Keep original municipality and split for matching
# -------------------------------
meta_fixed <- meta %>%
  mutate(Municipality_orig = Municipality_std) %>%  # keep original
  separate(Municipality_std, 
           into = c("Province_from_meta", "Municipality_from_meta"), 
           sep = "-", extra = "merge", fill = "right") %>%
  mutate(
    Province_from_meta = str_trim(Province_from_meta),
    Municipality_from_meta = str_trim(Municipality_from_meta),
    # clean barangay for matching
    Barangay_clean = str_to_lower(Barangay),
    Barangay_clean = str_replace_all(Barangay_clean, "-", " "),
    Barangay_clean = str_replace_all(Barangay_clean, "  ", " "),
    Barangay_clean = str_replace_all(Barangay_clean, "barangay 1", "barangay i"),
    Barangay_clean = str_replace_all(Barangay_clean, "barangay 2", "barangay ii"),
    Barangay_clean = str_replace_all(Barangay_clean, "barangay 3", "barangay iii"),
    Barangay_clean = str_trim(Barangay_clean)
  )

# -------------------------------
# 3. Clean barangay lookup similarly
# -------------------------------
barangay_lookup <- barangay_lookup %>%
  mutate(
    barangay_std_clean = str_to_lower(barangay_std),
    barangay_std_clean = str_replace_all(barangay_std_clean, "-", " "),
    barangay_std_clean = str_replace_all(barangay_std_clean, "  ", " "),
    barangay_std_clean = str_trim(barangay_std_clean)
  )

# -------------------------------
# 4. Join with barangay lookup using cleaned names
# -------------------------------
meta_linked <- meta_fixed %>%
  left_join(
    barangay_lookup,
    by = c(
      "Barangay_clean" = "barangay_std_clean",
      "Municipality_from_meta" = "municipality_std",
      "Province_from_meta" = "province_std"
    )
  )

# -------------------------------
# 5. Separate matched and unmatched rows
# -------------------------------
matched <- meta_linked %>% filter(!is.na(barangay_id))
unmatched <- meta_linked %>% filter(is.na(barangay_id))
dim(matched);dim(unmatched)

# -------------------------------
# 6. Create alias table template for unmatched
# Keep original Municipality_std and Province for reference
# -------------------------------
barangay_aliases <- unmatched %>%
  select(
    raw_name = Barangay, 
    municipality_std = Municipality_orig,  
    province_std = Province
  ) %>%
  mutate(barangay_std = NA) %>%
  select(raw_name, barangay_std, municipality_std, province_std)

# -------------------------------
# 7. Inspect and save
# -------------------------------
head(barangay_aliases)

write.csv(barangay_aliases, "processed_data/gis_data/barangay_aliases.csv", row.names = FALSE)
