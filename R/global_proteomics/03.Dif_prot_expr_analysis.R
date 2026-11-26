

# Now we can move on to the differential protein expression analysis:

# Cortex :

DEP_cortex_results <- diff_prot_expr_calc(cortex_expr_df)

saveRDS(
  DEP_cortex_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_cortex_results.rds"
)

# Gastronemius :

DEP_gastronemius_results <- diff_prot_expr_calc(gastrocnemius_expr_df)

 saveRDS(
   DEP_gastronemius_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_gastronemius_results.rds"
)


# Heart :

DEP_heart_results <- diff_prot_expr_calc(heart_expr_df)
saveRDS(
  DEP_heart_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_heart_results.rds"
)

    # Remove obsolete objects now that we have a result:
rm(heart_expr_df, heart_metadata, heart_processed_expression)


# Kidney :

DEP_kidney_results <- diff_prot_expr_calc(kidney_expr_df)

saveRDS(
  DEP_kidney_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_kidney_results.rds"
)

# Remove obsolete objects now that we have a result:
rm(kidney_expr_df, kidney_metadata, kidney_processed_expression)


# Liver :

DEP_liver_results <- diff_prot_expr_calc(liver_expr_df)

saveRDS(
  DEP_liver_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_liver_results.rds"
)

# Remove obsolete objects now that we have a result:
rm(liver_expr_df, liver_metadata, liver_processed_expression)

# Lung :

DEP_lung_results <- diff_prot_expr_calc(lung_expr_df)

saveRDS(
  DEP_lung_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_lung_results.rds"
)

# Remove obsolete objects now that we have a result:
rm(lung_expr_df, lung_metadata, lung_processed_expression)

# White adipose :

DEP_white.a_results <- diff_prot_expr_calc(white.a_expr_df)

saveRDS(
  DEP_white.a_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_white.a_results.rds"
)

# Remove obsolete objects now that we have a result:
rm(white.a_expr_df, white.a_metadata, white.a_processed_expression)




