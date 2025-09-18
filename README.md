# README
phl_ritm_rabv_n359.aln.fasta  
Alignment manually edited to remove 3bp gap 11665-10667 created by one sequence R6-21-4511

## To do
- obtain metadata for sequences from Wellcome paper and upload to redcap
- upload the associated sequences
- update records that have accession numbers
- pull together ncbi, bea and redcap data and filter to remove duplicates
- merge the redcap ready data for Bea so forms are combined


# Preparing data

cat raw_data/redcap_sequences_and_metadata/redcap_download_20250917_1441/*fasta raw_data/vgtk_ncbi_data/150725_metadata_coverage_90_country_Philippines_sequences.fa >processed_data/processed_sequences/gathered_sequences/170925_all_phl_sequences.fasta
