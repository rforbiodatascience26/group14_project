# Describe data


# Load library

library("tidyverse")


# Load augmented data

proteomics_aug <- read_tsv(
  "data/03_dat_aug.tsv"
)


# Sample overview

sample_overview <- proteomics_aug |>
  distinct(
    tissue,
    vial_label,
    study_group_timepoint
  ) |>
  group_by(
    tissue,
    study_group_timepoint
  ) |>
  summarise(
    n_samples = n(),
    .groups = "drop"
  )

sample_overview


# Number of proteins in each tissue

protein_overview <- proteomics_aug |>
  distinct(
    tissue,
    protein_id
  ) |>
  group_by(
    tissue
  ) |>
  summarise(
    n_proteins = n(),
    .groups = "drop"
  )

protein_overview


# Missing values in each tissue

missing_overview <- proteomics_aug |>
  group_by(
    tissue
  ) |>
  summarise(
    n_values = n(),
    n_missing = sum(
      is.na(protein_abundance)
    ),
    .groups = "drop"
  ) |>
  mutate(
    missing_percent =
      n_missing / n_values * 100
  )

missing_overview


# Missing values by group

missing_by_group <- proteomics_aug |>
  group_by(
    tissue,
    study_group_timepoint
  ) |>
  summarise(
    n_values = n(),
    n_missing = sum(
      is.na(protein_abundance)
    ),
    .groups = "drop"
  ) |>
  mutate(
    missing_percent =
      n_missing / n_values * 100
  )

missing_by_group

