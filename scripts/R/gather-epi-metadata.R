## Gather epi metadata for the philippines to link to sequences

## essels data
mydata=read.csv("gathered_epi_metadata/ph_redcap_2024.v1.csv")

## ncbi (vgtk) data
vgtk=read.csv("gathered_epi_metadata/metadata_coverage_90_country_Philippines.csv")

## mimiropa animals
mimi_animal=read.csv("gathered_epi_metadata/MIMAROPA outbreak 2022-24 - ANIMAL.csv")

## mimiropa humans
mimi_humans=read.csv("gathered_epi_metadata/MIMAROPA outbreak 2022-24 - HUMAN.csv")

## zhang 2025 paper
zhang=read.csv("gathered_epi_metadata/zhang2025/Supplementary Table 3.csv")


## Start my adding additional data from the Zhang paper to the vgtk version
## match by accession id
# Filter Zhang samples that are in vgtk
zhang_in_vgtk <- zhang[zhang$Acceccsion.No. %in% vgtk$primary_accession, ]

# Merge the Location info into vgtk
# We'll use match to keep original vgtk order
vgtk$geo_loc <- ifelse(
  vgtk$primary_accession %in% zhang_in_vgtk$Acceccsion.No., 
  zhang$Location[match(vgtk$primary_accession, zhang$Acceccsion.No.)],
  vgtk$geo_loc
)

# Optional: check how many were updated
sum(vgtk$accession %in% zhang_in_vgtk$accession)
