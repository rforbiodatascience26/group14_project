# Cross-tissue comparison of differential proteins


# Load libraries

library("tidyverse")


# Set results directory

results_dir <- str_c(
  "results/global_proteomics/",
  "Differential_expression_tables/"
)


# Load differential expression results

cortex_results <- readRDS(
  str_c(
    results_dir,
    "DEP_cortex_results.rds"
  )
) |>
  mutate(
    tissue = "Cortex"
  )

gastro_results <- readRDS(
  str_c(
    results_dir,
    "DEP_gastronemius_results.rds"
  )
) |>
  mutate(
    tissue = "Gastrocnemius"
  )

heart_results <- readRDS(
  str_c(
    results_dir,
    "DEP_heart_results.rds"
  )
) |>
  mutate(
    tissue = "Heart"
  )

kidney_results <- readRDS(
  str_c(
    results_dir,
    "DEP_kidney_results.rds"
  )
) |>
  mutate(
    tissue = "Kidney"
  )

liver_results <- readRDS(
  str_c(
    results_dir,
    "DEP_liver_results.rds"
  )
) |>
  mutate(
    tissue = "Liver"
  )

lung_results <- readRDS(
  str_c(
    results_dir,
    "DEP_lung_results.rds"
  )
) |>
  mutate(
    tissue = "Lung"
  )

white_adipose_results <- readRDS(
  str_c(
    results_dir,
    "DEP_white.a_results.rds"
  )
) |>
  mutate(
    tissue = "White adipose"
  )


# Combine results from all tissues

all_tissue_results <- bind_rows(
  cortex_results,
  gastro_results,
  heart_results,
  kidney_results,
  liver_results,
  lung_results,
  white_adipose_results
)


# Select significant proteins

significant_results <- all_tissue_results |>
  filter(
    padj < 0.05
  )


# Count significant proteins in each tissue

significant_by_tissue <- significant_results |>
  group_by(
    tissue
  ) |>
  summarise(
    n_significant = n(),
    .groups = "drop"
  )

print(significant_by_tissue)


# Count how many tissues each significant protein appears in

protein_tissue_count <- significant_results |>
  distinct(
    entrez_id,
    tissue
  ) |>
  group_by(
    entrez_id
  ) |>
  summarise(
    n_tissues = n(),
    .groups = "drop"
  ) |>
  mutate(
    response_type = case_when(
      n_tissues == 1 ~ "Tissue-specific",
      n_tissues > 1 ~ "Shared"
    )
  )


# Summarise shared and tissue-specific proteins

response_summary <- protein_tissue_count |>
  group_by(
    response_type
  ) |>
  summarise(
    n_proteins = n(),
    .groups = "drop"
  )

print(response_summary)


# Add classification to significant results

significant_classified <- significant_results |>
  left_join(
    protein_tissue_count,
    by = "entrez_id"
  )


# Count tissue-specific proteins in each tissue

tissue_specific_summary <- significant_classified |>
  filter(
    response_type == "Tissue-specific"
  ) |>
  group_by(
    tissue
  ) |>
  summarise(
    n_proteins = n(),
    .groups = "drop"
  )

print(tissue_specific_summary)


# Compare direction of shared protein responses

shared_direction <- significant_classified |>
  filter(
    response_type == "Shared"
  ) |>
  group_by(
    entrez_id
  ) |>
  summarise(
    n_positive = sum(
      coef_condition > 0
    ),
    n_negative = sum(
      coef_condition < 0
    ),
    n_tissues = n(),
    .groups = "drop"
  ) |>
  mutate(
    direction = case_when(
      n_positive == n_tissues ~ "Positive in all tissues",
      n_negative == n_tissues ~ "Negative in all tissues",
      TRUE ~ "Different direction"
    )
  )


# Summarise direction of shared proteins

direction_summary <- shared_direction |>
  group_by(
    direction
  ) |>
  summarise(
    n_proteins = n(),
    .groups = "drop"
  )

print(direction_summary)


# Summarise number of tissues per significant protein

tissue_count_summary <- protein_tissue_count |>
  group_by(
    n_tissues
  ) |>
  summarise(
    n_proteins = n(),
    .groups = "drop"
  ) |>
  mutate(
    n_tissues = factor(
      n_tissues
    )
  )

print(tissue_count_summary)


# Plot tissue-specific proteins by tissue

tissue_specific_plot <- tissue_specific_summary |>
  ggplot(
    aes(
      x = tissue,
      y = n_proteins
    )
  ) +
  geom_col() +
  labs(
    x = "Tissue",
    y = "Number of tissue-specific proteins"
  )

print(tissue_specific_plot)


# Plot direction of shared proteins

direction_plot <- direction_summary |>
  ggplot(
    aes(
      x = direction,
      y = n_proteins
    )
  ) +
  geom_col() +
  labs(
    x = "Direction pattern",
    y = "Number of shared proteins"
  )

print(direction_plot)


# Plot number of tissues per significant protein

tissue_count_plot <- tissue_count_summary |>
  ggplot(
    aes(
      x = n_tissues,
      y = n_proteins
    )
  ) +
  geom_col() +
  labs(
    x = "Number of tissues",
    y = "Number of significant proteins"
  )

print(tissue_count_plot)

