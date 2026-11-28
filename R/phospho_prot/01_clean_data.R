# Clean Data

## Load Libraries

library("tidyverse")


## Load Data

phospho_loaded <- read_tsv(
  "data/01_dat_phospho_load.tsv"
)

glimpse(phospho_loaded)

phospho_clean <- phospho_loaded |>
  mutate(
    vial_label = as.character(vial_label)
  ) |>
  filter(
    !is_contaminant
  )

glimpse(phospho_clean)

phenotype_data <- read_tsv(
  "data/_raw/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt"
)

phenotype_selected <- phenotype_data |>
  select(
    viallabel,
    registration___sex,
    study_group_timepoint
  ) |>
  mutate(
    viallabel = as.character(viallabel)
  )

glimpse(phenotype_selected)

phospho_clean <- phospho_clean |>
  left_join(
    phenotype_selected,
    by = c("vial_label" = "viallabel")
  )

phospho_clean |>
  filter(
    is.na(study_group_timepoint)
  ) |>
  distinct(
    vial_label
  )

## Save clean data
write_tsv(
  phospho_clean,
  "data/02_dat_phospho_clean.tsv"
)
