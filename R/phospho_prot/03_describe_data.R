
# Describe Data

## Load library

library("tidyverse")

## Load augmented data (03)

phospho_aug <- read_tsv(
  "data/03_dat_phospho_aug.tsv"
)

## Sample overview

sample_overview <- phospho_aug |>
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

## see how many proteins are in each tissue

phospho_overview <- phospho_aug |>
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

phospho_overview

## check how many missing values (NA) for each tissue

missing_overview <- phospho_aug |>
  group_by(
    tissue
  ) |>
  summarise(
    n_values = n(),
    n_missing = sum(
      is.na(phosphosite_abundance)
    ),
    .groups = "drop"
  ) |>
  mutate(
    missing_percent =
      n_missing / n_values * 100
  )

missing_overview

## Check the missing values per group

missing_by_group <- phospho_aug |>
  group_by(
    tissue,
    study_group_timepoint
  ) |>
  summarise(
    n_values = n(),
    n_missing = sum(
      is.na(phosphosite_abundance)
    ),
    .groups = "drop"
  ) |>
  mutate(
    missing_percent =
      n_missing / n_values * 100
  )

missing_by_group
