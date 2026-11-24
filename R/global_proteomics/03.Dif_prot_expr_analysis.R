

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
  DEP_cortex_results,
  file = "results/global_proteomics/Differential_expression_tables/DEP_gastronemius_results.rds"
)


# Heart :

DEP_heart_results <- diff_prot_expr_calc(heart_expr_df)

# Kidney :

DEP_kidney_results <- diff_prot_expr_calc(kidney_expr_df)

# Liver :

DEP_liver_results <- diff_prot_expr_calc(liver_expr_df)

# Lung :

DEP_lung_results <- diff_prot_expr_calc(lung_expr_df)

# White adipose :

DEP_wa_results <- diff_prot_expr_calc(wa_expr_df)






