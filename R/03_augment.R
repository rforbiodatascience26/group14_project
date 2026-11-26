# Augment data


# Load libraries

library("tidyverse")


# Load clean data

proteomics_clean <- read_tsv(
  "data/02_dat_clean.tsv"
)


# Add analysis variables

proteomics_aug <- proteomics_clean |>
  mutate(
    training_group = case_when(
      study_group_timepoint == "Control - 8 weeks" ~ "Control",
      study_group_timepoint == "Training - 1 week" ~ "Training",
      study_group_timepoint == "Training - 2 weeks" ~ "Training",
      study_group_timepoint == "Training - 4 weeks" ~ "Training",
      study_group_timepoint == "Training - 8 weeks" ~ "Training"
    )
  )


# Check training groups

proteomics_aug |>
  distinct(
    study_group_timepoint,
    training_group
  )


# Save augmented data

write_tsv(
  proteomics_aug,
  "data/03_dat_aug.tsv"
)
