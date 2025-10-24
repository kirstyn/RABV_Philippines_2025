# -----------------------------
# Standardise host data
# -----------------------------
library(dplyr)
library(stringr)

# Function to standardise host columns
standardise_hosts <- function(df) {
  
  df <- df %>%
    mutate(
      Host_scientific = case_when(
        str_to_lower(Host) %in% str_to_lower(c("Canis familiaris", "Canis lupus familiaris", "Dog", "Canine")) ~ "Canis lupus familiaris",
        str_to_lower(Host) %in% str_to_lower(c("Pig", "Sus scrofa")) ~ "Sus scrofa domesticus",
        str_to_lower(Host) %in% str_to_lower(c("Human")) ~ "Homo sapiens",
        str_to_lower(Host) %in% str_to_lower(c("Goat")) ~ "Capra hircus",
        str_to_lower(Host) %in% str_to_lower(c("Cow", "Bos taurus")) ~ "Bos taurus",
        str_to_lower(Host) %in% str_to_lower(c("Feline", "Felis catus")) ~ "Felis catus",
        TRUE ~ Host
      ),
      Host_common = case_when(
        Host_scientific == "Canis lupus familiaris" ~ "Dog",
        Host_scientific == "Sus scrofa domesticus"  ~ "Pig",
        Host_scientific == "Capra hircus"           ~ "Goat",
        Host_scientific == "Bos taurus"             ~ "Cow",
        Host_scientific == "Homo sapiens"           ~ "Human",
        Host_scientific == "Felis catus"            ~ "Cat",
        TRUE                                        ~ Host_scientific
      ),
      Host_type = case_when(
        Host_scientific %in% c("Canis lupus familiaris", "Felis catus") ~ "Domestic",   # pets
        Host_scientific %in% c("Sus scrofa domesticus", "Capra hircus", "Bos taurus") ~ "Livestock", # farm animals
        Host_scientific == "Homo sapiens" ~ "Human",
        TRUE ~ "Wildlife"
      )
    )
  
  return(df)
}

# -----------------------------
# Example: update main metadata objects
# -----------------------------
# If your metadata objects exist in the global environment
if(exists("phylo_meta")) phylo_meta <- standardise_hosts(phylo_meta)
if(exists("phylo_meta_corrected")) phylo_meta_corrected <- standardise_hosts(phylo_meta_corrected)
if(exists("phylo_meta_corrected2")) phylo_meta_corrected2 <- standardise_hosts(phylo_meta_corrected2)