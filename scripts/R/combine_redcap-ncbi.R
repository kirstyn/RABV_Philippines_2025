library(redcapAPI)
library(dplyr)
library(seqinr)

source("scripts/redcap-api.R")

# METADATA
source("scripts/redcap-getMetadata.R")
# redcap data = 'df_phl_90'

# ncbi dataset (pulled from vgtk)
ncbi <- read.csv("raw_data/150725_metadata_coverage_90_country_Philippines.csv")

# Remove sequences excluded by VGTK
ncbi_filt <- ncbi %>%
  filter(exclusion_status==1)

# Filter ncbi data to remove records already in redcap
# Note that there are later sequencing instances for samples on ncbi that have improved coverage/quality, so prefer to use these new versions rather than ncbi version
# Therefore, search by isolate ID to remove, instead of accession id

# manually remove samples that won't match due to sampleid typos (known issue for 4 samples)
manual_remove <- c("Z15-185", "Z14-152", "Z14-142", "Z12-012")

# samples in ncbi that have a match in redcap
ncbi_otherseq <- ncbi %>%
  filter(
    !isolate %in% df_phl$sample_id,       # remove records already in redcap
    !isolate %in% manual_remove           # also remove manual list
  )

# 412 sequences
write.csv(ncbi_otherseq,paste0("philippines/processed_data/combined_redcap-ncbi/",format(Sys.time(), "%Y%m%d"),"_filtered_ncbi.csv"),row.names=F)

# SEQUENCES
# now filter the sequences to same set
ncbi_seq <- read.fasta("philippines/raw_data/150725_metadata_coverage_90_country_Philippines_sequences.fa")

# Keep only sequences in ncbi_otherseq and not in typo list
keep_samples <- ncbi_otherseq$primary_accession[!ncbi_otherseq$isolate %in% manual_remove]

# Filter FASTA list by names
ncbi_seq_filtered <- ncbi_seq[names(ncbi_seq) %in% keep_samples]

# Optional: write filtered FASTA
write.fasta(
  sequences = ncbi_seq_filtered,
  names = names(ncbi_seq_filtered),
  file.out = paste0("philippines/processed_data/combined_redcap-ncbi/",format(Sys.time(), "%Y%m%d"),"_filtered_ncbi_sequences.fasta")
)

# JOIN DATASETS
# Join redcap and ncbi data (metadata + sequences) 
redcap_seq=read.fasta("philippines/processed_data/redcap_sequences_and_metadata/redcap_download_20250814_1219_all_seq.fasta")

# Combine them
combined_seq <- c(redcap_seq, ncbi_seq_filtered)
# Write to new FASTA
write.fasta(
  sequences = combined_seq,
  names = names(combined_seq),
  file.out = paste0("philippines/processed_data/combined_redcap-ncbi/",format(Sys.time(), "%Y%m%d"),"_phl_combined_sequences.fasta")
)



