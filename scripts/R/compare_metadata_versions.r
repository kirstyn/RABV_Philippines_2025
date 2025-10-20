# comparing metadata versions
# code will compare 2 dfs and identify discrepancies (additions/deletions/changes)

#devtools::install_github('alexsanjoseph/compareDF')

library(compareDF)
#https://github.com/alexsanjoseph/compareDF

#df1
old=read.csv("processed_data/processed_metadata/gathered_metadata_n797_20251020_154332.csv")


#df2
new=read.csv("processed_data/processed_metadata/gathered_metadata_n797_raddlCorrected_20251020_154332.csv")

# do the comparison
comparedf=compare_df(new, old, c("Sample_ID"))

# summarise
# for html (in viewer):
#create_output_table(comparedf)
# for excel output
#create_output_table(c1, output_type = 'xlsx', file_name = "outputs/discrepancies_15feb.xlsx")
create_output_table(comparedf, output_type = 'xlsx', file_name = "processed_data/processed_metadata/gathered_metadata_n797_compare.xls")

