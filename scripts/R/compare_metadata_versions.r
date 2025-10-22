# comparing metadata versions
# code will compare 2 dfs and identify discrepancies (additions/deletions/changes)

#devtools::install_github('alexsanjoseph/compareDF')

library(compareDF)


# import two datasets to compare (I call "new" the one with latest updates)
old=read.csv("processed_data/processed_metadata/gathered_metadata_n797_20251021_125349.csv")
new=read.csv("processed_data/processed_metadata/gathered_metadata_n797_raddl_and_manual_Corrected_20251021_125349.csv")


# Compare whole dataset: 

# do the comparison (group by Sample id, i.e. it compares the info for each sample id)
comparedf=compare_df(new, old, c("Sample_ID"))


# summarise

# for html (in R viewer):
#create_output_table(comparedf)

# for excel output
create_output_table(comparedf, output_type = 'xlsx', file_name = "processed_data/processed_metadata/gathered_metadata_n797_confirmManualCorrections.xlsx")


#Compare specific columns 

library(compareDF)

# Columns to compare
cols_to_compare <- c("Sample_ID", "Preferred_date")

# Subset the dataframes to only these columns
old_sub <- old[, cols_to_compare]
new_sub <- new[, cols_to_compare]

# Compare, matching on Sample_ID
comparedf <- compare_df(
  df_new = new_sub,
  df_old = old_sub,
  group_col = "Sample_ID"
)

# Export Excel
create_output_table(comparedf, output_type = "xlsx",
                    file_name = "processed_data/processed_metadata/gathered_metadata_n797_confirmManualCorrections.xlsx")
