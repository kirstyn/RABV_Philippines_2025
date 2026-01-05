# Tidy location data in Philippines metadata
library(dplyr)
library(stringdist)
library(stringr)
library(tidyr)

# Input data files
# Province to region mapping data:
map_province=read.csv("raw_data/gis_data/PHL_provinceTo_region_mapping.csv")
# Gathered metadata file to add geo standardisations to
input_file="processed_data/processed_metadata/gathered_metadata/final/05Jan26_gathered_metadata_n811_raddl_and_manual_Corrected.csv"
data_all=read.csv(input_file)


# Pull in adm centroids (generated prev by speedier, with adm1 added by kb) and use as reference table
adm <- read.csv("raw_data/gis_data/PHL_all_centroids.csv")
head(adm)

# Notes: barangays are named "village" in centroid file, municipality and village entries are hierarchical (contain province-municipality-village)

## Function to standardise Region information. Considers the different names used like Region I and Illocus. 

# improved standardise_region: handles "Central Luzon (III)", "Region III", "III (Central Luzon)", etc.
standardise_region <- function(region_col, adm, max_dist = 6) {
  library(stringr)
  library(stringdist)
  
  adm_regions <- adm %>%
    filter(Type == "Region") %>%
    mutate(Loc_ID = str_squish(Loc_ID)) %>%
    pull(Loc_ID)
  adm_regions_lower <- tolower(adm_regions)
  
  special_cases <- c(
    "mimiropa" = "Mimaropa Region",
    "mimaropa" = "Mimaropa Region",
    "ncr" = "National Capital Region (NCR)",
    "car" = "Cordillera Administrative Region (CAR)",
    "barmm" = "Bangsamoro Autonomous Region In Muslim Mindanao (BARMM)"
  )
  
  extract_roman <- function(s) {
    s2 <- gsub("[\\(\\)\\.]"," ", s) %>% str_squish()
    m <- str_match(s2, "\\b([IVXLCDM]+(?:-[A-Z])?)\\b")
    if (!is.na(m[1,2])) return(m[1,2]) else return(NA_character_)
  }
  
  sapply(region_col, function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    x_clean <- str_squish(as.character(x))
    x_lower <- tolower(x_clean)
    
    # special cases
    if (x_lower %in% names(special_cases)) return(special_cases[[x_lower]])
    
    # exact match
    idx_exact <- which(adm_regions_lower == x_lower)
    if (length(idx_exact) > 0) return(adm_regions[idx_exact[1]])
    
    # --- improved Roman numeral logic ---
    roman_token <- extract_roman(x_clean)
    if (!is.na(roman_token)) {
      idx_roman <- which(
        str_detect(tolower(adm_regions), paste0("\\(", tolower(roman_token), "\\)")) |
          str_detect(tolower(adm_regions), paste0("\\bregion\\s+", tolower(roman_token), "\\b"))
      )
      if (length(idx_roman) > 0) return(adm_regions[idx_roman[1]])
    }
    
    # partial substring matches (non-roman)
    idx_partial <- which(
      str_detect(adm_regions_lower, fixed(x_lower)) |
        str_detect(x_lower, fixed(adm_regions_lower))
    )
    if (length(idx_partial) > 0) return(adm_regions[idx_partial[1]])
    
    # fuzzy fallback (skip if string has Roman numeral — prevents Bicol(V)→II)
    if (is.na(roman_token)) {
      dists <- stringdist(tolower(x_clean), adm_regions_lower, method = "lv")
      best_idx <- which.min(dists)
      if (!is.infinite(dists[best_idx]) && dists[best_idx] <= max_dist)
        return(adm_regions[best_idx])
    }
    
    # remove roman token and retry partial
    if (!is.na(roman_token)) {
      no_roman <- gsub(roman_token, "", x_clean, ignore.case = TRUE)
      no_roman <- gsub("[\\(\\)\\-]", " ", no_roman) %>% str_squish()
      idx_partial2 <- which(
        str_detect(adm_regions_lower, fixed(tolower(no_roman))) |
          str_detect(tolower(no_roman), fixed(adm_regions_lower))
      )
      if (length(idx_partial2) > 0) return(adm_regions[idx_partial2[1]])
    }
    
    return(x_clean)
  }, USE.NAMES = FALSE)
}
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


# First, make corrections related to Metro Manila/NCR
# Technically metro manila has no provinces, but keeping Metropolitan Manila as province label 
data_all <- data_all %>%
  mutate(
    # Standardise Manila variants
    Province = if_else(
      Province %in% c("Metro Manila", "Metropolitan Manila"),
      "Metropolitan Manila",
      Province
    ),
    # Fix case for specific provinces
    Province = if_else(
      Province %in% c("Sultan kudarat", "South cotabato","Misamis oriental","La union","Ilocos sur"),
      str_to_title(Province),
      Province
    ),
    # Add NCR region label for Metropolitan Manila
    Region = if_else(
      Province == "Metropolitan Manila",
      "National Capital Region (NCR)",
      Region
    )
  )

# ---- Manual fixes for unusual data instances ----
data_all <- data_all %>%
  mutate(
    Province = as.character(Province),  # ensure character (prevents factor/NA coercion)
    Municipality = as.character(Municipality),
    
    # Fix specific sample
    Province = if_else(
      Sample_ID == "R11-21-61",
      "North Cotabato",
      Province,
      missing = Province  # keep existing values if NA encountered
    ),
    
    # Fix province when municipality is Zamboanga City
    Province = if_else(
      str_to_lower(Municipality) == "zamboanga city",
      "Zamboanga del Sur",
      Province,
      missing = Province
    )
  )


# Identify rows to update: 
rows_to_update <- which(is.na(data_all$Region))

# Update Region from map_province
na_region_idx <- rows_to_update[is.na(data_all$Region[rows_to_update])]

data_all$Region[na_region_idx] <- map_province$ADM1_EN[
  match(
    tolower(trimws(data_all$Province[na_region_idx])),
    tolower(trimws(map_province$ADM2_EN))
  )
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
# Check what Province is unmatched - all that is left should be "", which are for old ncbi sequences (Troupin paper)
 unique(unmatched_regions$Province)
 # Extract all unmatched rows
 unmatched_regions <- data_all %>%
   filter(is.na(Region_std)) %>%
   select(everything()); unmatched_regions
 

data_all$Region_std[1:50]          # inspect results

# now apply function to standardise provinces and municipality
data_std <- hierarchical_standardise_adm_simple(data_all, adm)

# Reorder columns so that *_std columns come next to originals
cols_order <- unlist(lapply(names(data_std), function(col) {
  if(grepl("_std$", col)) return(NULL)  # skip for now
  c(col, paste0(col, "_std")[paste0(col, "_std") %in% names(data_std)])  # add original then std if exists
}))

# Ensure we don’t lose any columns
cols_order <- unique(c(cols_order, names(data_std)[!names(data_std) %in% cols_order]))

data_std <- data_std[, cols_order]

# ---- Correct specific Region_std entries ----
data_std <- data_std %>%
  mutate(
    Region_std = case_when(
      Region_std == "Southwestern Tagalog Region" ~ "Mimaropa Region",
      Region_std == "Calabarzon (IV-A)" ~ "Region IV-A (Calabarzon)",
      TRUE ~ Region_std
    )
  )

length(unique(data_std$Region_std))

# ---- Assign Major Island Group (Adjusted for your Region_std entries) ----
data_std <- data_std %>%
  mutate(
    Major_Island = case_when(
      Region_std %in% c(
        "Region I (Ilocos Region)",
        "Region II (Cagayan Valley)",
        "Region III (Central Luzon)",
        "Region IV-A (Calabarzon)",
        "Calabarzon (IV-A)",
        "Mimaropa Region",
        "Southwestern Tagalog Region",
        "National Capital Region (NCR)",
        "Cordillera Administrative Region (CAR)",
        "Region V (Bicol Region)"
      ) ~ "Luzon",
      
      Region_std %in% c(
        "Region VI (Western Visayas)",
        "Region VII (Central Visayas)",
        "Region VIII (Eastern Visayas)"
      ) ~ "Visayas",
      
      Region_std %in% c(
        "Region IX (Zamboanga Peninsula)",
        "Region X (Northern Mindanao)",
        "Region XI (Davao Region)",
        "Region XII (Soccsksargen)",
        "Region XIII (Caraga)",
        "Bangsamoro Autonomous Region In Muslim Mindanao (BARMM)"
      ) ~ "Mindanao",
      
      TRUE ~ NA_character_
    )
  )

# Inspect
head(data_std[, grep("Province|Municipality", names(data_std))])

## write the standardised data to file. 
## Note this is just a stepping stone code, will still need manually checked and enhanced
write.csv(data_std, paste0(gsub(".csv","",input_file),"_stdGeo.csv"), row.names=F)
