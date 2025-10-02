# Tidy location data in Philippines metadata
library(dplyr)
library(stringdist)
library(stringr)
library(tidyr)

# Pull in adm centroids (generated prev by speedier, with adm1 added by kb) and use as reference table
adm <- read.csv("raw_data/gis_data/PHL_all_centroids.csv")
head(adm)

# Notes: barangays are named "village" in centroid file, municipality and village entries are hierarchical (contain province-municipality-village)

## Function to standardise Region information. Considers the different names used like Region I and Illocus. 

standardise_region <- function(region_col, adm) {
  # Reference regions
  adm_regions <- adm %>% filter(Type == "Region") %>% pull(Loc_ID)
  
  # Map known special cases (normalized lowercase)
  special_cases <- c(
    "mimiropa" = "Mimaropa Region",
    "mimaropa" = "Mimaropa Region",
    "ncr" = "National Capital Region (NCR)",
    "car" = "Cordillera Administrative Region (CAR)",
    "barmm" = "Bangsamoro Autonomous Region In Muslim Mindanao (BARMM)"
  )
  
  sapply(region_col, function(x) {
    if(is.na(x) || x == "") return(NA)
    
    x_lower <- tolower(str_squish(x))
    
    # Check special cases
    if(x_lower %in% names(special_cases)) return(special_cases[[x_lower]])
    
    # Remove Roman numerals in brackets and extra spaces for matching
    x_clean <- str_remove_all(x_lower, "\\([ivxlcdm]+\\)") %>% str_squish()
    
    # Try partial match anywhere in adm_regions
    matches <- adm_regions[str_detect(tolower(adm_regions), fixed(x_clean))]
    
    if(length(matches) > 0) {
      return(matches[1])
    } else {
      # fallback: return original
      return(x)
    }
  })
}

# Apply to data
data1 <- read.csv("/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/raw_data/gathered_epi_metadata/2018_workshop/2018_sequenced_collated_epi.csv") # subset of data to test
data1$Region_std <- standardise_region(data1$Region, adm)
table(is.na(data1$Region_std))  # check unmatched
data1$Region_std[1:50]          # inspect results

# function for the other adm levels, which have hierarchical structure

hierarchical_standardise_adm_simple <- function(df, adm, max_dist = 2) {
  
  df_std <- df
  
  # --- Province ---
  df_std$Province_std <- sapply(seq_len(nrow(df_std)), function(i) {
    prov <- df_std$Province[i]
    
    # Treat NA or empty string as NA
    if(is.na(prov) || prov == "") return(NA)
    
    candidates <- adm %>% 
      filter(Type == "Province") %>% 
      pull(Loc_ID)
    
    if(length(candidates) == 0) return(NA)
    
    candidates[which.min(stringdist(tolower(prov), tolower(candidates), method = "lv"))]
  })
  
  # --- Municipality ---
  df_std$Municipality_std <- sapply(seq_len(nrow(df_std)), function(i) {
    prov_std <- df_std$Province_std[i]
    mun      <- df_std$Municipality[i]
    
    # Treat NA or empty string as NA
    if(is.na(mun) || mun == "") return(NA)
    
    candidates <- adm %>% 
      filter(Type == "Municipality" & str_detect(Loc_ID, fixed(prov_std, ignore_case = TRUE))) %>% 
      pull(Loc_ID)
    
    if(length(candidates) == 0) return(NA)
    
    candidates[which.min(stringdist(tolower(mun), tolower(candidates), method = "lv"))]
  })
  
  # --- Optionally, leave Barangay untouched ---
  # df_std$Barangay_std <- df_std$Barangay
  
  return(df_std)
}

## Apply to data
data_std <- hierarchical_standardise_adm_simple(data1, adm)
head(data_std[, c("Province","Province_std","Municipality","Municipality_std")])
## seems to work pretty well! 

## Try with larger dataset
# load the map region to provine data
map_province=read.csv("raw_data/gis_data/PHL_provinceTo_region_mapping.csv")
data_all=read.csv("processed_data/processed_metadata/gathered_metadata_n794_20250922_155050_manuallyCorrected.csv")

# Identify rows to update: Source == "vgtk" and author contains Bacus or Cruz
rows_to_update <- which(data_all$Source == "vgtk" & grepl("Bacus|Cruz", data_all$Author, ignore.case = TRUE))

# Update Region from map_province
na_region_idx <- rows_to_update[is.na(data_all$Region[rows_to_update])]

data_all$Region[na_region_idx] <- map_province$ADM1_EN[
  match(data_all$Province[na_region_idx], map_province$ADM2_EN)
]

# Optional: check which provinces still have NA Region
remaining_unmatched <- data_all$Province[na_region_idx][
  is.na(data_all$Region[na_region_idx])
]
unique(remaining_unmatched)


# Apply functions
data_all$Region
data_all$Region_std <- standardise_region(data_all$Region, adm)
table(is.na(data_all$Region_std))  # check unmatched
# 3. Extract unmatched rows
unmatched_regions <- data_all %>%
  filter(is.na(Region_std)) %>%
  select(Province, Region, Region_std, everything()) 
# Optional: get unique unmatched regions
unique_unmatched_regions <- unique(unmatched_regions$Region)

data_all$Region_std[1:50]          # inspect results
data_std <- hierarchical_standardise_adm_simple(data_all, adm)
head(data_std[, c("Province","Province_std","Municipality","Municipality_std")])

## write the standardised data to file. 
## Note this is just a stepping stone code, will still need manually checked and enhanced
write.csv(data_std, "processed_data/processed_metadata/gathered_metadata_n794_20250922_155050_manuallyCorrected_Rstd.csv", row.names=F)
