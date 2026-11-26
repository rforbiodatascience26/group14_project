

# In this script we will generate the input file to perform the diff prot expression analysis

# CORTEX DATA
cortex_expr_df <- pcdf_cleaner(cortex_processed_expression)
  
# Add info of the batch

cortex_expr_df <- cortex_expr_df |> 
  inner_join(Cortex_metadata,by = "vial_label")

colnames(cortex_expr_df)[6] <- "batch"
  
# GASTROCNEMIUS DATA

gastrocnemius_expr_df <- pcdf_cleaner(gastrocnemius_processed_expression)

# Add info of the batch

gastrocnemius_expr_df  <- gastrocnemius_expr_df  |> 
  inner_join(gastrocnemius_metadata,by = "vial_label")

colnames(gastrocnemius_expr_df)[6] <- "batch"

# Heart

heart_expr_df <- pcdf_cleaner(heart_processed_expression)

# Add info of the batch

heart_expr_df  <- heart_expr_df  |> 
  inner_join(heart_metadata,by = "vial_label")

colnames(heart_expr_df)[6] <- "batch"

# Kidney

kidney_expr_df <- pcdf_cleaner(kidney_processed_expression)

# Add info of the batch

kidney_expr_df  <- kidney_expr_df  |> 
  inner_join(kidney_metadata,by = "vial_label")

colnames(kidney_expr_df)[6] <- "batch"

# Liver

liver_expr_df <- pcdf_cleaner(liver_processed_expression)

# Add info of the batch

liver_expr_df  <- liver_expr_df  |> 
  inner_join(liver_metadata,by = "vial_label")

colnames(liver_expr_df)[6] <- "batch"

# Lung

lung_expr_df <- pcdf_cleaner(lung_processed_expression)

# Add info of the batch

lung_expr_df  <- lung_expr_df  |> 
  inner_join(lung_metadata,by = "vial_label")

colnames(lung_expr_df)[6] <- "batch"


# White adipose

white.a_expr_df <- pcdf_cleaner(white.a_processed_expression)

# Add info of the batch

white.a_expr_df  <- white.a_expr_df  |> 
  inner_join(white.a_metadata,by = "vial_label")

colnames(white.a_expr_df)[6] <- "batch"



