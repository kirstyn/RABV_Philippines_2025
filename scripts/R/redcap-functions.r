## Useful redcap related functions 

# redcap data downloads as two forms, sequencing and diagnostic. This function will merge the two forms together by sample id

redcap_form_merge <- function(rcon, df, form1 = "diagnostic", form2 = "sequencing") {
  # Pull metadata directly from the project
  meta <- exportMetaData(rcon)
  
  # Split rows by repeating status
  df_form1 <- df %>%
    filter(is.na(redcap_repeat_instrument) | redcap_repeat_instrument == "")
  
  df_form2 <- df %>%
    filter(redcap_repeat_instrument == stringr::str_to_title(form2))
  
  # Get field lists from metadata
  form1_fields <- meta %>%
    filter(form_name == form1) %>%
    pull(field_name)
  
  form2_fields <- meta %>%
    filter(form_name == form2) %>%
    pull(field_name)
  
  # Select only relevant columns
  df_form1_clean <- df_form1 %>%
    select(any_of(c("sample_id", form1_fields)))
  
  df_form2_clean <- df_form2 %>%
    select(any_of(c("sample_id", "redcap_repeat_instance", form2_fields)))
  
  # Join on sample_id
  df_merged <- df_form2_clean %>%
    left_join(df_form1_clean, by = "sample_id")
  
  return(df_merged)
}