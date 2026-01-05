library(tibble)
library(readr)

# Define run metadata
run_log <- tibble(
  Run_ID = "PHL_explore7",
  Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  # full date-time
  User = Sys.info()[["user"]],           # who ran it
  BEAST_version = "v1.10.5",             # update as needed
  Alignment = "PHL_sequences_withDatesLabeled_2025-12-02_n786.nextalign.aligned_concat/nc",
  XML_file = "analysis/xml/PHL_explore7.xml",
  Substitution_model = "coding=GTR+gamma, CP112, noncoding=GTR+gamma",
  Clock_model = "Relaxed Lognormal",
  Clock_rate_prior = "Lognormal(mean=0.00028113, sd=0.4)",
  Tree_prior = "Skygrid; 5 points since 86",
  Chain_length = 5e8,
  Sampling_freq = 50000,
  Burn_in = "10%",
  SS_path_steps = 200,
  SS_chain_length_per_step = 2e6,
  SS_sampling_freq = 1000,
  Location_trait = "na",
  Migration_model = "na",
  BSSVS = FALSE,
  Location_strategy = "na",
  Notes = "Reduced skyline intervals, imposed prior on clock rate"
)

# File to write log into
log_file <- "analysis/BEAST/BEAST_runs_log.csv"

# If file exists, append without header; otherwise, create with header
if (!file.exists(log_file)) {
  write_csv(run_log, log_file, append = FALSE)
} else {
  write_csv(run_log, log_file, append = TRUE)
}