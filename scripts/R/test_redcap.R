# install.packages("devtools") # if not already installed

devtools::install_github("RAGE-toolkit/rabvRedcapProcessing@devel", force = T)
library(rabvRedcapProcessing)

new_data=read_data(filepath = '/Users/kirstyn.brunker/philippines/metadata_fromBea_29jul2025/missing_metadata_wellcomePaper.xlsx - Redcap-v1.csv')

dict=read_and_parse_dict('~/Downloads/RABVlab_DataDictionary_2025-07-31.csv')

harmonised=compare_cols_to_dict(dayta=new_data, dictPath = '~/Downloads/RABVlab_DataDictionary_2025-07-31.csv')

scan_mismatched_levels(new_data, dict, col_to_check = "ngs_prep")
read_data
