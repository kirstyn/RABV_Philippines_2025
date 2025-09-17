## Compare public data with in-house data to:
#- identify any missing records in redcap
#- filter the ncbi data to only additional data (i.e not ours)
library(dplyr)

## ncbi data
ncbi=read.csv("/Users/kirstyn.brunker/philippines/files_received/fasta_fromNCBI_15jul2025/metadata_coverage_90_country_Philippines.csv")

# latest redcap data
redcap=read.csv("/Users/kirstyn.brunker/GitHub/rage-redcap-developer/philippines/processed_data/redcap_sequences_and_metadata/redcap_download_20250814_114620250814_1146redcap_meta_phl.csv")

# filter to only sequencing rows
redcap_seq=redcap %>%
  filter(redcap_repeat_instrument=="sequencing")
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
write.csv(redcap_updated, "/Users/kirstyn.brunker/philippines/redcap_imports/130825_genbankAccessionUpdates.csv" ,row.names=F)



# samples in ncbi that have no match in redcap:
ncbi.only=ncbi.team %>%
  filter(!isolate %in% redcap$sample_id)
ncbi.only$isolate
#[1] "Z15-185" "Z14-142" "Z14-152" "Z12-012" : problem related to sample id typos (missing hypen)
