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


# Preparing metadata   
This workflow gathers and standardises metadata from multiple sources into a unified dataset.

Step 1 – Gather and Standardise Metadata

Script: scripts/R/gather-epi-metadata.R
Purpose:
    •    Combine metadata from various sources.
    •    Standardise entries (e.g., parsing multiple date formats, trimming whitespace).
    •    Remove duplicate entries.

Outputs:
    1.    Complete, standardised metadata for all sequences.
    2.    Metadata with additional RADDL corrections applied.
    3.    Metadata with manual corrections applied (currently pending full validation).

⸻

Step 2 – Compare Metadata Versions

Script: scripts/R/compare_metadata_versions.R
Purpose:
    •    Compare successive metadata outputs to verify that corrections have been incorporated correctly.
    •    Identify additions, deletions, or changes between versions.

⸻

Step 3 – Standardise Administrative Levels

Script: scripts/R/standardise_adm_levels.R
Purpose:
    •    Populate missing Region values for VGTK sequences using the Province column.
    •    Standardise Municipality and Province entries to match the ADM levels defined in the shapefiles.
    •    Apply string-distance matching to correct typos or variant spellings.
    
Outputs:
     1.    Updated CSV file with both original and std columns, placing standardised columns next to their original counterparts for traceability.
    Standardised columns include:
    •    Region_std → Standardised Region names
    •    Province_std → Corrected Province names
    •    Municipality_std → Corrected Municipality names

