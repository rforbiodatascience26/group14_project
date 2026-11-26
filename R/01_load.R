# Load data

## Load libraries

library("tidyverse")

## Gastro data

gastro_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t55-gastrocnemius_prot-pr_ratio-results_v2.0.txt"
)

gastro_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t55-gastrocnemius_prot-pr_vial-metadata_v2.0.txt"
)


## Cortex data

cortex_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t53-cortex_prot-pr_ratio-results_v2.0.txt"
)

cortex_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t53-cortex_prot-pr_vial-metadata_v2.0.txt"
)

## Heart data

heart_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t58-heart_prot-pr_ratio-results_v2.0.txt"
)

heart_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t58-heart_prot-pr_vial-metadata_v2.0.txt"
)


## Kidney data

kidney_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t59-kidney_prot-pr_ratio-results_v2.0.txt"
)

kidney_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t59-kidney_prot-pr_vial-metadata_v2.0.txt"
)


## Liver Data

liver_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t68-liver_prot-pr_ratio-results_v2.0.txt"
)

liver_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t68-liver_prot-pr_vial-metadata_v2.0.txt"
)

## Lung data


lung_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t66-lung_prot-pr_ratio-results_v2.0.txt"
)

lung_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t66-lung_prot-pr_vial-metadata_v2.0.txt"
)

## White Adipose Data


white_adipose_proteomics <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t70-white-adipose_prot-pr_ratio-results_v2.0.txt"
)

white_adipose_metadata <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_t70-white-adipose_prot-pr_vial-metadata_v2.0.txt"
)


## Load Phenotype Data

phenotype_data <- read_tsv(
  "../data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt"
)

proteomics_loaded <- bind_rows(
  cortex_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Cortex"),
  
  gastro_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Gastrocnemius"),
  
  heart_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Heart"),
  
  kidney_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Kidney"),
  
  liver_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Liver"),
  
  lung_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "Lung"),
  
  white_adipose_proteomics |>
    pivot_longer(
      cols = 10:69,
      names_to = "vial_label",
      values_to = "protein_abundance"
    ) |>
    mutate(tissue = "White adipose")
)

write_tsv(
  proteomics_loaded,
  "../data/01_dat_load.tsv"
)



# Metadata file Loading

meta_pheno_dic <- read_tsv("data/_raw/motrpac_pass1b-06_pheno_dictionary_merged_v4.0.txt")

meta_vial_label <- read_tsv("data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt")