# Clean data


# Load libraries

library("tidyverse")


# Load data

proteomics_loaded <- read_tsv(
  "data/01_dat_load.tsv"
)


# Inspect loaded data

glimpse(proteomics_loaded)


# Clean proteomics data

proteomics_clean <- proteomics_loaded |>
  mutate(
    vial_label = as.character(vial_label)
  ) |>
  filter(
    !is_contaminant
  )


# Inspect cleaned proteomics data

glimpse(proteomics_clean)


# Load phenotype data

phenotype_data <- read_tsv(
  "data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt"
)


# Select phenotype variables of interest

phenotype_selected <- phenotype_data |>
  select(
    viallabel,
    registration___sex,
    study_group_timepoint
  ) |>
  mutate(
    viallabel = as.character(viallabel)
  )


# Inspect phenotype data

glimpse(phenotype_selected)


# Join phenotype information to proteomics data

proteomics_clean <- proteomics_clean |>
  left_join(
    phenotype_selected,
    by = c("vial_label" = "viallabel")
  )


# Check for unmatched vial labels

proteomics_clean |>
  filter(
    is.na(study_group_timepoint)
  ) |>
  distinct(
    vial_label
  )


# Save clean data

write_tsv(
  proteomics_clean,
  "data/02_dat_clean.tsv"
)
