library(Biostrings)
library(DECIPHER)
library(dplyr)
library(stringr)

# Folder containing FASTA files
fasta_folder <- "/Users/kirstyn.brunker/philippines/redcap-troubleshooting/repeating-instances"

# List all FASTA files
fasta_files <- list.files(fasta_folder, pattern = "\\.fasta$", full.names = TRUE)

# Extract sample_id from filenames
sample_ids <- str_extract(basename(fasta_files), "^[^_]+")

# Create a data frame linking files to sample_ids
file_df <- data.frame(
  file = fasta_files,
  sample_id = sample_ids,
  stringsAsFactors = FALSE
)

# Function to compare sequences for a single sample
compare_sample_sequences <- function(files, sample_id) {
  
  # Read all sequences for this sample
  seqs <- do.call(c, lapply(files, readDNAStringSet))
  
  # Align sequences
  alignment <- AlignSeqs(seqs, anchor = NA)
  
  # Convert to character matrix
  aln_matrix <- as.matrix(alignment)
  
  # Initialize SNP list
  snp_list <- list()
  n_seqs <- nrow(aln_matrix)
  
  for (i in 1:(n_seqs-1)) {
    for (j in (i+1):n_seqs) {
      seq1 <- aln_matrix[i, ]
      seq2 <- aln_matrix[j, ]
      # Positions with differences, ignoring gaps and Ns
      positions <- which(seq1 != seq2 & !(seq1 %in% c("-", "N")) & !(seq2 %in% c("-", "N")))
      if (length(positions) > 0) {
        snp_df <- data.frame(
          sample_id = sample_id,
          seq1 = rownames(aln_matrix)[i],
          seq2 = rownames(aln_matrix)[j],
          position = positions,
          allele_seq1 = seq1[positions],
          allele_seq2 = seq2[positions],
          stringsAsFactors = FALSE
        )
        snp_list[[length(snp_list)+1]] <- snp_df
      }
    }
  }
  
  # Combine SNPs for this sample
  snp_table <- if (length(snp_list) > 0) do.call(rbind, snp_list) else NULL
  
  # Sample summary
  sample_summary <- data.frame(
    sample_id = sample_id,
    n_sequences = n_seqs,
    total_snps = ifelse(is.null(snp_table), 0, nrow(snp_table)),
    max_pairwise_snps = ifelse(is.null(snp_table), 0, 
                               max(table(snp_table$seq1, snp_table$seq2))),
    stringsAsFactors = FALSE
  )
  
  list(snp_table = snp_table, summary = sample_summary)
}

# Loop through all samples
results <- lapply(unique(file_df$sample_id), function(sid) {
  files <- file_df$file[file_df$sample_id == sid]
  compare_sample_sequences(files, sid)
})

# Combine SNP tables
all_snps_df <- do.call(rbind, lapply(results, function(x) x$snp_table))

# Combine summaries
all_summary_df <- do.call(rbind, lapply(results, function(x) x$summary))

# Save outputs
write.csv(all_snps_df, "/Users/kirstyn.brunker/philippines/redcap-troubleshooting/repeating-instances/all_samples_snp_table.csv", row.names = FALSE)
write.csv(all_summary_df, "/Users/kirstyn.brunker/philippines/redcap-troubleshooting/repeating-instances/all_samples_summary.csv", row.names = FALSE)

# Quick check
head(all_snps_df)
head(all_summary_df)
