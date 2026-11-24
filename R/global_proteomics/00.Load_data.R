
library(tidyverse)

# In this script we will load the data necessary for the global proteomics analysis:

# cortex data
Cortex_metadata <- read_tsv("data/_raw/global_proteomics/t53-cortex/prot-pr/motrpac_pass1b-06_t53-cortex_prot-pr_vial-metadata_v2.0.txt")

Cortex_dirty_counts <- read_tsv("data/_raw/global_proteomics/t53-cortex/prot-pr/motrpac_pass1b-06_t53-cortex_prot-pr_rii-results_v2.0.txt")

cortex_processed_counts <- read_tsv("data/_raw/global_proteomics/t53-cortex/prot-pr/motrpac_pass1b-06_t53-cortex_prot-pr_ratio-results_v2.0.txt")

