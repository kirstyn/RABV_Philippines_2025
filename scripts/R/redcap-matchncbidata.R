## Compare public data with in-house data to:
#- identify any missing records in redcap
#- filter the ncbi data to only additional data (i.e not ours)
library(dplyr)
library(stringr)

## ncbi data
ncbi=read.csv("raw_data/vgtk_ncbi_data/150725_metadata_coverage_90_country_Philippines.csv")

# latest redcap data
redcap_seq=read.csv("raw_data/redcap_sequences_and_metadata/redcap_download_20251023_134520251023_1345redcap_meta_phl.csv")

head(redcap_seq)
head(ncbi)

# filter based on author contains brunker to get all records submitted by team
ncbi.team=ncbi %>%
  filter(str_detect(authors, regex("brunker", ignore_case = TRUE)))

# samples in ncbi that have a match in redcap
ncbi.redcap.matches=redcap_seq %>%
  filter(sample_id %in% ncbi.team$isolate)

# check to see if all have accession numbers
ncbi.redcap.matches$genbank_accession

# what ones don't
no.accession.redcap=ncbi.redcap.matches %>%
  filter(is.na(genbank_accession) | genbank_accession == "")

which(redcap_seq$sample_id %in% ncbi.team$isolate)
# Join redcap with ncbi.team accessions
redcap_updated <- redcap_seq %>%
  left_join(
    ncbi.team %>%
      select(sample_id = isolate, primary_accession),  # match names for join
    by = "sample_id"
  ) %>%
  mutate(
    genbank_accession = if_else(
      is.na(genbank_accession) | genbank_accession == "",
      primary_accession,  # fill from ncbi.team
      genbank_accession
    )
  ) %>%
  select(-primary_accession)  # remove helper column
write.csv(redcap_updated, "redcap_imports/130825_genbankAccessionUpdates.csv" ,row.names=F)



# samples in ncbi that have no match in redcap:
ncbi.only=ncbi.team %>%
  filter(!isolate %in% redcap$sample_id)
ncbi.only$isolate
#[1] "Z15-185" "Z14-142" "Z14-152" "Z12-012" : problem related to sample id typos (missing hypen)


# Filter out all records submitted by your team
ncbi.additional <- ncbi %>%
  filter(!str_detect(authors, regex("brunker", ignore_case = TRUE)))

# Write to CSV
write.csv(ncbi.additional, "redcap_imports/ncbi_additional_data_excl_brunker.csv", row.names = FALSE)

# Read the NCBI sequences
ncbi.seq <- read.fasta("raw_data/vgtk_ncbi_data/150725_metadata_coverage_90_country_Philippines_sequences.fa", seqtype = "DNA", as.string = TRUE)

# Keep only sequences that are in ncbi.additional
# Assuming sequence names in FASTA match `isolate` column
seq.names.to.keep <- ncbi.additional$isolate
ncbi.seq.filtered <- ncbi.seq[names(ncbi.seq) %in% seq.names.to.keep]

# Write filtered sequences to new FASTA
write.fasta(sequences = ncbi.seq.filtered, 
            names = names(ncbi.seq.filtered), 
            file.out = "redcap_imports/ncbi_additional_sequences_excl_brunker.fa")
