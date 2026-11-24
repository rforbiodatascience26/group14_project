

# In this script we will generate the input file to perform the diff prot expression analysis

# CORTEX DATA
cortex_expr_df <- pcdf_cleaner(cortex_processed_expression)
  
# Add info of the batch

cortex_expr_df <- cortex_expr_df |> 
  inner_join(Cortex_metadata,by = "vial_label")

colnames(cortex_expr_df)[6] <- "batch"

cortex_expr_df <- mutate(cortex_expr_df,
  condition = factor(condition, levels = c("Control", "Training")),
  sex = factor(sex, levels = c("M", "F")),
  batch = factor(batch, levels = c("S1","S2","S3","S4","S5","S6"))
)
  
# GASTROCNEMIUS DATA

gastrocnemius_expr_df <- pcdf_cleaner(gastrocnemius_processed_expression)

# Add info of the batch

gastrocnemius_expr_df  <- gastrocnemius_expr_df  |> 
  inner_join(gastrocnemius_metadata,by = "vial_label")

colnames(gastrocnemius_expr_df)[6] <- "batch"

# Heart

# Kidney

# Liver

# Lung

# White adipose