

# Functions that will be needed to run the global proteomics analyisis:


# Function 1: processed counts dataframe cleaner:

pcdf_cleaner <- function(processed_expression_df) {
  
  clean_expression_df <- processed_expression_df |> 
    
    # Step 1: Eliminate all Contaminant proteins from the df
    filter(!str_detect(protein_id, "Contaminant")) |> 
    
    # Step 2: Select columns of interest
    select(-protein_id,                     # We will keep entrez_ID as the identifier
           - gene_symbol, 
           - redundant_ids,
           -organism_name,
           -is_contaminant,
           -num_peptides,
           -percent_coverage,
           -protein_score) |> 
    
    # Step 3: Pivot longer Based on protein - Sample
    pivot_longer(
      cols = !starts_with("entrez"),
      names_to = "vial_label",
      values_to= "expression",
      values_drop_na = TRUE) |> 
    
    # Step 4: Add information of the sample
    inner_join(clean_metadata_table, by = "vial_label") |> 
    # Has to be an inner join, because the processed expression has samples taken 
    # at other moments beside 8 weeks, and the metadata has all the possible samples
    # we are interested in, (control or endurance exercise at 8 weeks).
  
  return(clean_expression_df)
}


