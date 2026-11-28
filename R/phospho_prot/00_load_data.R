# Load data

###  For this step the existance of the phosphoproteomics counts is necessary, ç
# but if the file data/01_dat_phospho_load.tsv is already in the repository you 
# can skip the execution of this script (and the need of having the raw 
# files in the repository)

## Load libraries

library("tidyverse")

## Gastro data

gastro_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t55-gastrocnemius_prot-ph_ratio-results_v2.0.txt"
)

gastro_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t55-gastrocnemius_prot-ph_vial-metadata_v2.0.txt"
)

## Cortex data

cortex_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t53-cortex_prot-ph_ratio-results_v2.0.txt"
)

cortex_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t53-cortex_prot-ph_vial-metadata_v2.0.txt"
)

## Heart data

heart_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t58-heart_prot-ph_ratio-results_v2.0.txt"
)

heart_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t58-heart_prot-ph_vial-metadata_v2.0.txt"
)

## Kidney data

kidney_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t59-kidney_prot-ph_ratio-results_v2.0.txt"
)

kidney_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t59-kidney_prot-ph_vial-metadata_v2.0.txt"
)

## Liver Data

liver_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t68-liver_prot-ph_ratio-results_v2.0.txt"
)

liver_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t68-liver_prot-ph_vial-metadata_v2.0.txt"
)

## Lung data

lung_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t66-lung_prot-ph_ratio-results_v2.0.txt"
)

lung_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t66-lung_prot-ph_vial-metadata_v2.0.txt"
)

## White Adipose Data

white_adipose_phospho_proteomics <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t70-white-adipose_prot-ph_ratio-results_v2.0.txt"
)

white_adipose_phospho_metadata <- read_tsv(
  "data/_raw/motrpac_pass1b-06_t70-white-adipose_prot-ph_vial-metadata_v2.0.txt"
)

## Load Phenotype Data

phenotype_data <- read_tsv(
  "data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt"
)


phospho_loaded <- bind_rows(
  cortex_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Cortex"),
  
  gastro_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Gastrocnemius"),
  
  heart_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Heart"),
  
  kidney_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Kidney"),
  
  liver_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Liver"),
  
  lung_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "Lung"),
  
  white_adipose_phospho_proteomics |>
    pivot_longer(
      cols = 12:last_col(),
      names_to = "vial_label",
      values_to = "phosphosite_abundance"
    ) |>
    mutate(tissue = "White adipose")
)

write_tsv(
  phospho_loaded,
  "data/01_dat_phospho_load.tsv"
)

# Metadata file Loading

meta_pheno_dic <- read_tsv("data/_raw/motrpac_pass1b-06_pheno_dictionary_merged_v4.0.txt")

meta_vial_label <- read_tsv("data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt")
