# Augment Data

## Load Libraries

library("tidyverse")

## Load Clean Data

phospho_clean <- read_tsv(
  "data/02_dat_phospho_clean.tsv"
)

## Add analysis variables

phospho_aug <- phospho_clean |>
  mutate(
    training_group = case_when(
      study_group_timepoint == "Control - 8 weeks" ~ "Control",
      study_group_timepoint == "Training - 1 week" ~ "Training",
      study_group_timepoint == "Training - 2 weeks" ~ "Training",
      study_group_timepoint == "Training - 4 weeks" ~ "Training",
      study_group_timepoint == "Training - 8 weeks" ~ "Training"
    )
  )

phospho_aug |>
  distinct(
    study_group_timepoint,
    training_group
  )

write_tsv(
  phospho_aug,
  "data/03_dat_phospho_aug.tsv"
)
