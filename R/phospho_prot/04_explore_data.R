# Explore phosphoproteomics data

library(tidyverse)
library(pheatmap)

# Read data
phospho_aug <- read_tsv(
  "data/03_dat_phospho_aug.tsv",
  col_types = cols(
    vial_label = col_character(),
    .default = col_guess()
  )
)


# 1. Distribution of phosphoproteomics values

ggplot(
  phospho_aug,
  aes(x = phosphosite_abundance)
) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "white"
  ) +
  facet_wrap(
    ~ tissue,
    scales = "free"
  ) +
  theme_bw() +
  labs(
    title = "Distribution of phosphoproteomics values",
    x = "Phosphosite abundance",
    y = "Frequency"
  )


# 2. Distribution by experimental group

ggplot(
  phospho_aug,
  aes(
    x = training_group,
    y = phosphosite_abundance,
    fill = training_group
  )
) +
  geom_boxplot(
    outlier.alpha = 0.2
  ) +
  facet_wrap(
    ~ tissue,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "Phosphoproteomics values by experimental group",
    x = "Experimental group",
    y = "Phosphosite abundance"
  )


# 3. Sample variability

sample_variability <- phospho_aug |>
  group_by(
    tissue,
    vial_label,
    training_group
  ) |>
  summarise(
    median_phospho = median(
      phosphosite_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ggplot(
  sample_variability,
  aes(
    x = vial_label,
    y = median_phospho,
    color = training_group
  )
) +
  geom_point(size = 2) +
  facet_wrap(
    ~ tissue,
    scales = "free"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 5
    )
  ) +
  labs(
    title = "Median phosphoproteomics values per sample",
    x = "Sample",
    y = "Median phosphosite abundance"
  )


# 4. Correlation between samples

tissues <- c(
  "Cortex",
  "Gastrocnemius",
  "Heart",
  "Kidney",
  "Liver",
  "Lung",
  "White adipose"
)

for (tissue_name in tissues) {
  
  tissue_matrix <- phospho_aug |>
    filter(tissue == tissue_name) |>
    group_by(
      ptm_id,
      vial_label
    ) |>
    summarise(
      phosphosite_abundance = median(
        phosphosite_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = vial_label,
      values_from = phosphosite_abundance
    ) |>
    column_to_rownames("ptm_id") |>
    as.matrix()
  
  tissue_cor <- cor(
    tissue_matrix,
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  
  pheatmap(
    tissue_cor,
    main = paste("Sample correlation -", tissue_name)
  )
}


# 5. Principal Component Analysis (PCA)

run_tissue_pca <- function(tissue_name) {
  
  # Filter tissue
  tissue_data <- phospho_aug |>
    filter(tissue == tissue_name) |>
    mutate(
      vial_label = as.character(vial_label)
    )
  
  # Create phosphosite x sample matrix
  tissue_matrix <- tissue_data |>
    group_by(
      ptm_id,
      vial_label
    ) |>
    summarise(
      phosphosite_abundance = median(
        phosphosite_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = vial_label,
      values_from = phosphosite_abundance
    ) |>
    column_to_rownames("ptm_id") |>
    as.matrix()
  
  # Keep phosphosites detected in at least 70% of samples
  detected <- rowMeans(!is.na(tissue_matrix)) >= 0.70
  
  tissue_matrix <- tissue_matrix[
    detected,
    ,
    drop = FALSE
  ]
  
  # Impute missing values using the median
  for (i in seq_len(nrow(tissue_matrix))) {
    
    missing <- is.na(tissue_matrix[i, ])
    
    if (any(missing)) {
      tissue_matrix[i, missing] <- median(
        tissue_matrix[i, ],
        na.rm = TRUE
      )
    }
  }
  
  # Remove phosphosites with zero variance
  site_sd <- apply(
    tissue_matrix,
    1,
    sd
  )
  
  tissue_matrix <- tissue_matrix[
    is.finite(site_sd) & site_sd > 0,
    ,
    drop = FALSE
  ]
  
  # PCA
  pca_result <- prcomp(
    t(tissue_matrix),
    center = TRUE,
    scale. = TRUE
  )
  
  # Percentage of variance
  variance <- summary(pca_result)$importance[2, ] * 100
  
  # PCA scores
  pca_scores <- as.data.frame(pca_result$x) |>
    rownames_to_column("vial_label") |>
    left_join(
      tissue_data |>
        distinct(
          vial_label,
          training_group
        ),
      by = "vial_label"
    )
  
  # PCA plot
  p <- ggplot(
    pca_scores,
    aes(
      x = PC1,
      y = PC2,
      color = training_group
    )
  ) +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      title = paste("PCA -", tissue_name),
      x = paste0("PC1 (", round(variance[1], 1), "%)"),
      y = paste0("PC2 (", round(variance[2], 1), "%)"),
      color = "Experimental group"
    )
  
  return(p)
}


# Run PCA for all tissues

for (tissue_name in tissues) {
  print(run_tissue_pca(tissue_name))
}
