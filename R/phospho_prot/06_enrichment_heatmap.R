
## 1. Load packages
install.packages("BiocManager")
BiocManager::install("ReactomePA")

library(tidyverse)
library(ReactomePA)
library(pheatmap)

## 2. Create output folders

dir.create("results/phosphoproteomics/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/phosphoproteomics/heatmaps", recursive = TRUE, showWarnings = FALSE)

## 3. Load data

# Differential analysis results
diff_results_all <- read_tsv(
  "results/phosphoproteomics/06_differential_phospho_results.tsv",
  show_col_types = FALSE
)

# Augmented phosphoproteomics data
phospho_8weeks <- read_tsv(
  "data/03_dat_phospho_aug.tsv",
  col_types = cols(
    vial_label = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

## 4. Identify significant phosphosites

significant_phospho <- diff_results_all |>
  filter(
    !is.na(padj),
    padj < 0.05
  )

## 5. Functional enrichment

# Get unique proteins associated with significant sites
significant_proteins <- significant_phospho |>
  select(tissue, phosphosite_id) |>
  distinct() |>
  left_join(
    phospho_8weeks |>
      select(phosphosite_id, entrez_id) |>
      distinct(),
    by = "phosphosite_id"
  ) |>
  filter(!is.na(entrez_id))

enrichment_results <- list()

for (t in unique(significant_proteins$tissue)) {
  
  genes <- significant_proteins |>
    filter(tissue == t) |>
    pull(entrez_id) |>
    unique()
  
  enrichment_results[[t]] <- enrichPathway(
    gene = genes,
    organism = "rat",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # Save enrichment table
  write_tsv(
    as.data.frame(enrichment_results[[t]]),
    paste0("results/phosphoproteomics/enrichment/enrichment_", gsub(" ", "_", t), ".tsv")
  )
  
  # Save enrichment plot
  if (nrow(as.data.frame(enrichment_results[[t]])) > 0) {
    
    enrichment_plot <- dotplot(
      enrichment_results[[t]],
      showCategory = 15
    ) +
      ggtitle(paste("Functional enrichment -", t))
    
    ggsave(
      paste0(
        "results/phosphoproteomics/enrichment/enrichment_",
        gsub(" ", "_", t),
        ".png"
      ),
      plot = enrichment_plot,
      width = 9,
      height = 7,
      dpi = 300
    )
  }
}

## 6. Heatmaps

for (t in unique(significant_phospho$tissue)) {
  
  # Select the 50 most significant phosphosites
  top_sites <- significant_phospho |>
    filter(tissue == t) |>
    arrange(padj) |>
    slice_head(n = 50) |>
    pull(phosphosite_id)
  
  # Prepare data
  heatmap_data <- phospho_8weeks |>
    filter(
      tissue == t,
      study_group_timepoint %in% c(
        "Control - 8 weeks",
        "Training - 8 weeks"
      ),
      phosphosite_id %in% top_sites
    ) |>
    select(
      phosphosite_id,
      vial_label,
      study_group_timepoint,
      phosphosite_abundance
    ) |>
    group_by(phosphosite_id, vial_label) |>
    summarise(
      phosphosite_abundance = mean(
        phosphosite_abundance,
        na.rm = TRUE
      ),
      condition = first(study_group_timepoint),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = vial_label,
      values_from = phosphosite_abundance
    )
  
  # Create matrix
  heatmap_matrix <- heatmap_data |>
    select(-phosphosite_id, -condition) |>
    as.matrix()
  
  rownames(heatmap_matrix) <- heatmap_data$phosphosite_id
  
  # Remove rows without detected values
  heatmap_matrix <- heatmap_matrix[
    rowSums(!is.na(heatmap_matrix)) > 0,
    ,
    drop = FALSE
  ]
  
  # Impute missing values using row median
  for (i in seq_len(nrow(heatmap_matrix))) {
    
    missing <- is.na(heatmap_matrix[i, ])
    
    if (any(missing)) {
      heatmap_matrix[i, missing] <- median(
        heatmap_matrix[i, ],
        na.rm = TRUE
      )
    }
  }
  
  # Remove phosphosites with zero variance
  heatmap_matrix <- heatmap_matrix[
    apply(heatmap_matrix, 1, sd) > 0,
    ,
    drop = FALSE
  ]
  
  # Standardize values by phosphosite
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  
  # Sample annotation
  annotation_data <- heatmap_data |>
    select(vial_label, condition) |>
    distinct() |>
    as.data.frame()
  
  rownames(annotation_data) <- annotation_data$vial_label
  annotation_data$vial_label <- NULL
  
  annotation_data <- annotation_data[
    colnames(heatmap_matrix),
    ,
    drop = FALSE
  ]
  
  # Generate and save heatmap
  pheatmap(
    heatmap_matrix,
    annotation_col = annotation_data,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    scale = "none",
    na_col = "grey90",
    main = paste("Phosphorylation patterns -", t),
    filename = paste0(
      "results/phosphoproteomics/heatmaps/heatmap_",
      gsub(" ", "_", t),
      ".png"
    ),
    width = 10,
    height = 10
  )
}
