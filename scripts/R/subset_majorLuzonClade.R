library(ape)        # For phylogenetic trees
library(Biostrings) # For fasta handling
library(treeio)
library(phytools)
# Subsetting Luzon major clade from the full tree

tree <- read.beast("analysis/BEAST/outputs/PHL_explore7and8_subsampled.hipstr.tre")

# extract phylo
phy <- tree@phylo

# now plot
plot(phy, cex=0.5)
nodelabels(cex=0.5)
tiplabels(cex=0.4)

# major luzon node label is 1024
# Replace with your node number
node_major <- 1024

# Get all descendants
descendants <- getDescendants(phy, node_major)

# Keep only tip numbers
tip_numbers <- descendants[descendants <= length(phy$tip.label)]

# Get tip labels
tip_labels_major <- phy$tip.label[tip_numbers]

length(tip_labels_major)  # check how many sequences
head(tip_labels_major)

fasta <- readDNAStringSet("processed_data/processed_sequences/PHL_sequences_withDatesLabeled_2026-01-05_n786.fasta")
length(fasta)  # check total number of sequences

# Keep only sequences that are in tip_labels_major
fasta_subset <- fasta[names(fasta) %in% tip_labels_major]

# Optional: reorder to match the tree tip order
fasta_subset <- fasta_subset[match(tip_labels_major, names(fasta_subset))]

length(fasta_subset)  # check how many sequences were kept
head(names(fasta_subset))

writeXStringSet(fasta_subset, "processed_data/processed_sequences/PHL_major_Luzon_subset.fasta")

## Metadata
meta <- read.csv("processed_data/processed_metadata/gathered_metadata/final/PHL_metadata_matchedWithDates_2026-01-05_n786.csv")
head(meta)

# Extract base sample ID from tree tips
tip_ids_base <- sapply(strsplit(tip_labels_major, "__"), `[`, 1)

head(tip_ids_base)
length(tip_ids_base)

meta_subset <- meta %>%
  dplyr::filter(seq_name  %in% tip_ids_base)

nrow(meta_subset)
head(meta_subset)

write.csv(meta_subset,"processed_data/processed_metadata/gathered_metadata/final/PHL_major_Luzon_subset.csv")

