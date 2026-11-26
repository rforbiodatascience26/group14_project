# Analysis 1: overview of training response across tissues


# Load libraries

library("tidyverse")


# Load augmented data

proteomics_aug <- read_tsv(
  "data/03_dat_aug.tsv"
)


# Select 8 week groups

proteomics_8weeks <- proteomics_aug |>
  filter(
    study_group_timepoint %in% c(
      "Control - 8 weeks",
      "Training - 8 weeks"
    )
  )


# Check sample numbers

proteomics_8weeks |>
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


# Calculate mean protein abundance

protein_means <- proteomics_8weeks |>
  group_by(
    tissue,
    protein_id,
    gene_symbol,
    study_group_timepoint
  ) |>
  summarise(
    mean_abundance = mean(
      protein_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Inspect protein means

glimpse(protein_means)


# Compare control and training

protein_changes <- protein_means |>
  pivot_wider(
    names_from = study_group_timepoint,
    values_from = mean_abundance
  )


# Calculate training effect

protein_changes <- protein_changes |>
  mutate(
    training_effect =
      `Training - 8 weeks` - `Control - 8 weeks`
  )


# Check missing effects

protein_changes |>
  summarise(
    n_proteins = n(),
    n_missing_effect = sum(
      is.na(training_effect)
    )
  )


# Remove proteins where the effect cannot be calculated

protein_changes_complete <- protein_changes |>
  drop_na(
    training_effect
  )


# Summarise training response by tissue

tissue_overview <- protein_changes_complete |>
  group_by(
    tissue
  ) |>
  summarise(
    n_proteins = n(),
    mean_abs_effect = mean(
      abs(training_effect)
    ),
    .groups = "drop"
  )

tissue_overview


# Plot training effect across tissues

tissue_overview |>
  ggplot(
    aes(
      x = tissue,
      y = mean_abs_effect
    )
  ) +
  geom_col() +
  labs(
    x = "Tissue",
    y = "Mean absolute training effect"
  )

