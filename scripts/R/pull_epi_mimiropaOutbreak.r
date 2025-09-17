## Gathering from mimiropa outbreak sheets if necessary
# code needs to be incorporated into gather-epi-metadata, won't work standalone without phylo_meta 
library(dplyr)
## mimiropa animals
mimi_animal <- read.csv(
  "raw_data/gathered_epi_metadata/MIMAROPA outbreak 2022-24 - ANIMAL.csv",
  skip = 1
)
# add region before merging with other data
mimi_animal$REGION="MIMAROPA"
mimi_animal$Source="SPEEDIER"
mimi_animal <- mimi_animal %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

## mimiropa humans
mimi_humans=read.csv("raw_data/gathered_epi_metadata/MIMAROPA outbreak 2022-24 - HUMAN.csv")
mimi_humans$Region="MIMAROPA"
mimi_humans$Host="Human"
mimi_humans$Source="SPEEDIER"
mimi_humans <- mimi_humans %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

## -------
# Extract the sequenced samples from mimiropa datasets
mimi_animal %in% mimi_humans
# None of the col names are the same! arrrggghh
names(mimi_animal)
# sequenced mimiropa animal
mimi_animal_seq <- mimi_animal %>%
  filter(!is.na(Samples.sequenced..Seq.ID) & Samples.sequenced..Seq.ID != "")
dim(mimi_animal_seq)
# sequenced mimiropa human
mimi_humans_seq <- mimi_humans %>%
  filter(!is.na(Samples.sequenced..Seq.ID) & Samples.sequenced..Seq.ID != "")
dim(mimi_humans_seq)



## Map mimi_humans_seq to phylo_meta
mimi_humans_col_map <- c(
  "Samples.sequenced..Seq.ID" = "Sample_ID",
  "CASE_NO"                   = "Case_no",
  "Specimen_Collection"       = "Sample_type",
  "Province"                  = "Province",
  "Municipality"              = "Municipality",
  "Barangay"                  = "Barangay",
  "Region"                    = "Region",
  "Host"                      = "Host",
  "Lat"                       = "Latitude",
  "Lon"                       = "Longitude",
  "Date"                      = "Preferred_date",
  "Source"                    = "Source"
)

# Keep only those columns in mimi_humans_seq that exist in the map
existing_cols <- names(mimi_humans_col_map)[names(mimi_humans_col_map) %in% names(mimi_humans_seq)]

# Subset mapping to only existing cols
col_map_filtered <- mimi_humans_col_map[existing_cols]

# Select and rename
filled_dataset <- mimi_humans_seq %>%
  select(all_of(names(col_map_filtered))) %>%
  rename_with(~ col_map_filtered[.x], .cols = everything())

# Combine into phylo_meta structure
phylo_meta <- bind_rows(phylo_meta, filled_dataset)

# --------------
## Map mimi_humans_seq to phylo_meta
mimi_animal_col_map <- c(
  "Samples.sequenced..Seq.ID" = "Sample_ID",
  "CASE_NO"                   = "Case_no",
  "SAMPLE_TYPE"       = "Sample_type",
  "PROVINCE"                  = "Province",
  "MUNICIPALITY"              = "Municipality",
  "BARANGAY"                  = "Barangay",
  "REGION"                    = "Region",
  "BITING_ANIMAL"                      = "Host",
  "Date"                      = "Preferred_date",
  "Source"                    = "Source"
)

# Keep only those columns in mimi_humans_seq that exist in the map
existing_cols <- names(mimi_animal_col_map)[names(mimi_animal_col_map) %in% names(mimi_animal_seq)]

# Subset mapping to only existing cols
col_map_filtered <- mimi_animal_col_map[existing_cols]

# Select and rename
filled_dataset <- mimi_animal_seq %>%
  select(all_of(names(col_map_filtered))) %>%
  rename_with(~ col_map_filtered[.x], .cols = everything())

# Combine into phylo_meta structure
phylo_meta <- bind_rows(phylo_meta, filled_dataset)
