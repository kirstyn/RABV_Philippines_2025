library(dplyr)
library(stringr)

# Load raw PSGC data
psgc_raw <- read.csv("raw_data/gis_data/psa_data/PSGC-3Q-2025-Publication-Datafile.csv", stringsAsFactors = FALSE)

# Ensure PSGC is 10 digits
psgc_raw <- psgc_raw %>%
  mutate(psgc_10 = str_pad(as.character(X10.digit.PSGC), 10, pad = "0"))

# Provinces
provinces <- psgc_raw %>%
  filter(Geographic.Level == "Prov") %>%
  mutate(province_code = str_sub(psgc_10, 3, 5)) %>%
  select(province_code, province_std = Name) %>%
  mutate(province_std = str_trim(province_std))

# Municipalities
municipalities <- psgc_raw %>%
  filter(Geographic.Level %in% c("Mun", "City","SubMun")) %>%
  mutate(
    province_code = str_sub(psgc_10, 3, 5),
    municipality_code = str_sub(psgc_10, 6, 7)
  ) %>%
  select(province_code, municipality_code, municipality_std = Name) %>%
  mutate(municipality_std = str_trim(municipality_std))

# Barangays
barangays <- psgc_raw %>%
  filter(Geographic.Level == "Bgy") %>%
  mutate(
    province_code = str_sub(psgc_10, 3, 5),
    municipality_code = str_sub(psgc_10, 6, 7),
    barangay_code = str_sub(psgc_10, 8, 10),
    barangay_id = str_sub(psgc_10, 3, 10)
  ) %>%
  select(psgc_10, barangay_id, barangay_std = Name,
         municipality_code, province_code, barangay_code) %>%
  # join municipality names
  left_join(municipalities, by = c("province_code", "municipality_code")) %>%
  # join province names
  left_join(provinces, by = "province_code") %>%
  select(psgc_10, barangay_id, barangay_std, municipality_std, province_std) %>%
  # clean names
  mutate(
    barangay_std = str_to_lower(barangay_std) %>%
      str_replace_all("-", " ") %>%
      str_replace_all("  ", " ") %>%
      str_trim(),
    municipality_std = str_trim(municipality_std),
    province_std = str_trim(province_std)
  )

# Save cleaned lookup
write.csv(barangays, "processed_data/gis_data/psgc_barangay_municipality_province_lookup_09Jan26.csv", row.names=FALSE)

