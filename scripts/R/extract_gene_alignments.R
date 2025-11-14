#extract partitions for beast
library(seqinr)
#BiocManager::install("ORFik")
#BiocManager::install("DECIPHER")
library(ORFik)
library(DECIPHER)
library(devtools)
#install_github("thibautjombart/apex")
#install.packages("apex")
library(apex)
library(stringr)
library(seqinr)
library(Biostrings)

#input file
file="processed_data/processed_sequences/PHL_sequences_withDates_2025-11-03_n786_nextalign/PHL_sequences_withDates_2025-11-03_n786.nextalign.aligned.fasta"

#output files
dir.create(file.path(paste(dirname(file), "gene_alignments", sep="/")), showWarnings = FALSE)
# set output file location (same as input) and prefix
newfiles=paste(dirname(file), "gene_alignments", gsub(".fasta|.fst","",basename(file)), sep="/")

# read in sequences as string set for ORF function
string=readDNAStringSet(file)
writeXStringSet(string, paste0(newfiles, ".fasta"), format="fasta")

# find most complete genome and search for ORFs
chosen=which.max(str_count(as.character(string), "A|T|G|C"))
#chosen.seq=RemoveGaps(string[chosen], removeGaps = "common")
genes=as.data.frame(findORFs(string[chosen], startCodon = "ATG", minimumLength =200))

# it can't find the correct M gene start point so have to pull out manually
genes=genes[order(genes$start),]
genes=genes[-3,]
find.m=as.data.frame(findORFs(string[chosen], startCodon = "ATG", minimumLength =200, longestORF = F))
find.m=find.m[order(find.m$start),]
m=find.m[which(find.m$start>=2468 & find.m$width==609),]
#join with other genes
genes=rbind(genes, m)
genes=genes[order(genes$start),]
genes=genes[,-c(1,2)]
genes$gene=NA
genes$gene=c("n","p","m","g","l")

# read in sequences as fasta file
seq=read.fasta(paste0(newfiles, ".fasta"))

# split into coding partitions (5 genes)
# based on ORF positions
n=getFrag(seq, begin=genes$start[1],end=genes$end[1])
p=getFrag(seq, begin=genes$start[2],end=genes$end[2])
m=getFrag(seq, begin=genes$start[3],end=genes$end[3])
g=getFrag(seq, begin=genes$start[4],end=genes$end[4])
l=getFrag(seq, begin=genes$start[5],end=genes$end[5])

n_strip=n[sapply(n, function(x) any(x != "-"))]
p_strip=n[sapply(p, function(x) any(x != "-"))]
m_strip=n[sapply(m, function(x) any(x != "-"))]
g_strip=n[sapply(g, function(x) any(x != "-"))]
l_strip=n[sapply(l, function(x) any(x != "-"))]

# output partitions as fasta files
write.fasta(n_strip,names=sapply(n_strip, function(x) attr(x, "seqMother")), paste(newfiles,"n.fasta",sep="_"))
write.fasta(p_strip,names=seq_mothers <- sapply(p_strip, function(x) attr(x, "seqMother")), paste(newfiles,"p.fasta",sep="_"))
write.fasta(m_strip,names=seq_mothers <- sapply(m_strip, function(x) attr(x, "seqMother")), paste(newfiles,"m.fasta",sep="_"))
write.fasta(g_strip,names=seq_mothers <- sapply(g_strip, function(x) attr(x, "seqMother")), paste(newfiles,"g.fasta",sep="_"))
write.fasta(l_strip,names=seq_mothers <- sapply(l_strip, function(x) attr(x, "seqMother")), paste(newfiles,"l.fasta",sep="_"))


## Produce a concatenated gene file 
# Ensure all partitions have the same order of sequences
# Use the "seqMother" attribute as an identifier
get_names <- function(part) sapply(part, function(x) attr(x, "seqMother"))

seq_names <- get_names(n)  # assuming all partitions have same seqMother names
stopifnot(all(seq_names == get_names(p)))  # check order consistency
stopifnot(all(seq_names == get_names(m)))
stopifnot(all(seq_names == get_names(g)))
stopifnot(all(seq_names == get_names(l)))

# Concatenate sequences for each isolate
concat_seqs <- lapply(seq_along(seq_names), function(i) {
  c(n[[i]], p[[i]], m[[i]], g[[i]], l[[i]])
})

# Write concatenated sequences as FASTA
concat_file <- paste0(newfiles, "_concat.fasta")
write.fasta(sequences = concat_seqs,
            names = seq_names,
            file.out = concat_file)

#non-coding seq
nc <- function(x){
  a1=1
  a2=genes$start[1]-1
  b1=genes$end[1]+1
  b2=genes$start[2]-1
  c1=genes$end[2]+1
  c2=genes$start[3]-1
  d1=genes$end[3]+1
  d2=genes$start[4]-1
  e1=genes$end[4]+1
  e2=genes$start[5]-1
  f1=genes$end[5]+1
  f2=width(string[1])
  return(x[c(a1:a2, b1:b2,c1:c2,d1:d2,e1:e2,f1:f2)])
}

ncod=mapply(nc,seq)
ncod=t(ncod)
ncod=as.table(ncod)
new=rep(NA,length(seq))
new=data.frame(id=rownames(ncod), seq=NA)
for (i in 1:nrow(ncod)){
  new$seq[i]=paste(ncod[i,], collapse="")
}

#This creates an empty character vector, with length twice the length of your table; then puts the values from column1 in every second position starting at 1, and the values of column2 in every second position starting at 2.
Xfasta <- character(nrow(new) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", new$id)
Xfasta[c(FALSE, TRUE)] <- new$seq

#then write using writeLines:

writeLines(Xfasta,  paste(newfiles,"nc.fasta",sep="_"))

## ---- Codon-consistency checks ----
## Check for real stop codons in gene partitions

check_stops <- function(fasta_file) {
  dna <- readDNAStringSet(fasta_file)
  aa <- translate(dna, if.fuzzy.codon = "X")  # Ns/gaps -> "X"
  
  has_stop <- vcountPattern("*", aa) > 0
  
  stop_positions <- lapply(aa, function(x) {
    which(strsplit(as.character(x), "")[[1]] == "*")
  })
  
  data.frame(
    sequence = names(dna),
    has_stop = has_stop,
    stop_positions = I(stop_positions),   # keep as list column
    stringsAsFactors = FALSE
  )
}

# Run stop codon check for each gene partition
stop_reports <- list(
  n = check_stops(paste(newfiles, "n.fasta", sep="_")),
  p = check_stops(paste(newfiles, "p.fasta", sep="_")),
  m = check_stops(paste(newfiles, "m.fasta", sep="_")),
  g = check_stops(paste(newfiles, "g.fasta", sep="_")),
  l = check_stops(paste(newfiles, "l.fasta", sep="_"))
)

# Example: view all sequences with a stop codon in N gene
subset(stop_reports$n, has_stop)
subset(stop_reports$p, has_stop)
subset(stop_reports$m, has_stop)
subset(stop_reports$g, has_stop)
subset(stop_reports$l, has_stop)

###############################################
# CREATE IQ-TREE PARTITION FILE WITH INTERGENIC BLOCKS
###############################################

# Coding partitions
coding_names  <- genes$gene
coding_starts <- genes$start
coding_ends   <- genes$end

# Interspersed non-coding partitions
nc_blocks <- rbind(
  c(1, genes$start[1]-1),
  c(genes$end[1]+1, genes$start[2]-1),
  c(genes$end[2]+1, genes$start[3]-1),
  c(genes$end[3]+1, genes$start[4]-1),
  c(genes$end[4]+1, genes$start[5]-1),
  c(genes$end[5]+1, width(string[1]))
)
# Keep only valid blocks
nc_blocks <- nc_blocks[nc_blocks[,1] <= nc_blocks[,2], , drop=FALSE]

# Build lines for partition file
partition_file <- paste0(newfiles, "_partitions.txt")
lines <- character(0)

# Coding partitions
for(i in seq_along(coding_names)){
  lines <- c(lines, sprintf("DNA, %s = %d-%d", coding_names[i], coding_starts[i], coding_ends[i]))
}

# Non-coding partition, all blocks in one line separated by commas
nc_ranges <- apply(nc_blocks, 1, function(x) paste0(x[1], "-", x[2]))
lines <- c(lines, sprintf("DNA, nc = %s", paste(nc_ranges, collapse=",")))

# Write to file
writeLines(lines, partition_file)
