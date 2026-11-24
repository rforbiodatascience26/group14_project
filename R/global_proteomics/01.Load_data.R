
library(tidyverse)

# In this script we will load the data necessary for the global proteomics analysis:

# cortex data
Cortex_metadata <- read_tsv("data/_raw/global_proteomics/t53-cortex/prot-pr/motrpac_pass1b-06_t53-cortex_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

cortex_processed_expression <- read_tsv("data/_raw/global_proteomics/t53-cortex/prot-pr/motrpac_pass1b-06_t53-cortex_prot-pr_ratio-results_v2.0.txt")

# Gastrocnemius data

gastrocnemius_metadata <- read_tsv("data/_raw/global_proteomics/t55-gastrocnemius/prot-pr/motrpac_pass1b-06_t55-gastrocnemius_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

gastrocnemius_processed_expression <- read_tsv("data/_raw/global_proteomics/t55-gastrocnemius/prot-pr/motrpac_pass1b-06_t55-gastrocnemius_prot-pr_ratio-results_v2.0.txt")

# Heart

# Kidney

# Liver

# Lung

# White adipose
