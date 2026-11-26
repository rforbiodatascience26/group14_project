
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

heart_metadata <- read_tsv("data/_raw/global_proteomics/t58-heart/prot-pr/motrpac_pass1b-06_t58-heart_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

heart_processed_expression <- read_tsv("data/_raw/global_proteomics/t58-heart/prot-pr/motrpac_pass1b-06_t58-heart_prot-pr_ratio-results_v2.0.txt")

# Kidney

kidney_metadata <- read_tsv("data/_raw/global_proteomics/t59-kidney/prot-pr/motrpac_pass1b-06_t59-kidney_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

kidney_processed_expression <- read_tsv("data/_raw/global_proteomics/t59-kidney/prot-pr/motrpac_pass1b-06_t59-kidney_prot-pr_ratio-results_v2.0.txt")

# Liver

liver_metadata <- read_tsv("data/_raw/global_proteomics/t68-liver/prot-pr/motrpac_pass1b-06_t68-liver_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

liver_processed_expression <- read_tsv("data/_raw/global_proteomics/t68-liver/prot-pr/motrpac_pass1b-06_t68-liver_prot-pr_ratio-results_v2.0.txt")

# Lung

lung_metadata <- read_tsv("data/_raw/global_proteomics/t66-lung/prot-pr/motrpac_pass1b-06_t66-lung_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

lung_processed_expression <- read_tsv("data/_raw/global_proteomics/t66-lung/prot-pr/motrpac_pass1b-06_t66-lung_prot-pr_ratio-results_v2.0.txt")

# White adipose

white.a_metadata <- read_tsv("data/_raw/global_proteomics/t70-white-adipose/prot-pr/motrpac_pass1b-06_t70-white-adipose_prot-pr_vial-metadata_v2.0.txt") |> 
  select(-tmt11_channel)

white.a_processed_expression <- read_tsv("data/_raw/global_proteomics/t70-white-adipose/prot-pr/motrpac_pass1b-06_t70-white-adipose_prot-pr_ratio-results_v2.0.txt")

