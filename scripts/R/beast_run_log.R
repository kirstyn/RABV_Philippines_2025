library(tibble)
library(readr)

# Define run metadata
run_log <- tibble(
  Run_ID = "PER_RABV_2024_SS_01",
  Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  # full date-time
  User = Sys.info()[["user"]],           # who ran it
  BEAST_version = "v1.10.5",             # update as needed
  Alignment = "150825_filt_seq.aln.fasta",
  XML_file = "peru/xmls/PER_RABV_2024_SS_01.xml",
  Clock_model = "Relaxed Lognormal",
  Clock_rate_prior = "Lognormal(real space mean=0.000216,sd=0.0000927,initial=0.0001993)",
  Tree_prior = "Coalescent Exponential",
  Chain_length = 2e8,
  Sampling_freq = 10000,
  Burn_in = "20%",
  SS_path_steps = 100,
  SS_chain_length_per_step = 2e6,
  SS_sampling_freq = 200,
  Location_trait = "District",
  Migration_model = "Asymmetric",
  BSSVS = TRUE,
  Location_strategy = "High-signal districts",
  Notes = "First SS run with updated rabies sequences"
)

# File to write log into
log_file <- "peru/BEAST_runs/BEAST_runs_log.csv"

# If file exists, append without header; otherwise, create with header
if (!file.exists(log_file)) {
  write_csv(run_log, log_file, append = FALSE)
} else {
  write_csv(run_log, log_file, append = TRUE)
}