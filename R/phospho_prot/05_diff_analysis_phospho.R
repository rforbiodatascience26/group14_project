## Load packages and data

library(tidyverse)

phospho_aug <- read_tsv(
  "data/03_dat_phospho_aug.tsv",
  col_types = cols(
    vial_label = col_character(),
    .default = col_guess()
  )
)

## Select Control and Training at 8 weeks
phospho_8weeks <- phospho_aug |>
  filter(
    study_group_timepoint %in% c(
      "Control - 8 weeks",
      "Training - 8 weeks"
    )
  ) |>
  mutate(
    condition = factor(
      study_group_timepoint,
      levels = c(
        "Control - 8 weeks",
        "Training - 8 weeks"
      )
    )
  )

## Creation of phosphosite identifiers 

phospho_8weeks <- phospho_8weeks |>
  mutate(
    phosphosite_id = paste(
      protein_id,
      ptm_id,
      flanking_sequence,
      sep = "_"
    )
  )

phospho_8weeks <- phospho_8weeks |>
  rename(
    sex = registration___sex
  )

## Function for differential analysis

diff_phospho_calc <- function(data) {
  
  diff_phospho <- data |>
    group_by(phosphosite_id) |>
    nest() |>
    mutate(
      
      model = map(
        data,
        possibly(
          ~ lm(
            phosphosite_abundance ~ condition + sex,
            data = .x
          ),
          otherwise = NULL
        )
      ),
      
      coef_condition = map_dbl(
        model,
        ~ {
          if (is.null(.x)) return(NA_real_)
          
          cf <- coef(.x)
          idx <- grep("^condition", names(cf))
          
          if (length(idx) == 0) return(NA_real_)
          
          unname(cf[idx[1]])
        }
      ),
      
      p_value = map_dbl(
        model,
        ~ {
          if (is.null(.x)) return(NA_real_)
          
          mat <- summary(.x)$coefficients
          idx <- grep("^condition", rownames(mat))
          
          if (length(idx) == 0) return(NA_real_)
          
          unname(mat[idx[1], "Pr(>|t|)"])
        }
      )
      
    ) |>
    mutate(
      padj = p.adjust(p_value, method = "BH")
    ) |>
    select(
      phosphosite_id,
      coef_condition,
      p_value,
      padj
    )
  
  return(diff_phospho)
}

## Differential analysis for the different tissues

tissues <- c(
  "Cortex",
  "Gastrocnemius",
  "Heart",
  "Kidney",
  "Liver",
  "Lung",
  "White adipose"
)

diff_results <- list()

for (tissue_name in tissues) {
  
  tissue_data <- phospho_8weeks |>
    filter(tissue == tissue_name)
  
  diff_results[[tissue_name]] <- diff_phospho_calc(
    tissue_data
  )
}

diff_results_all <- bind_rows(
  lapply(
    names(diff_results),
    function(tissue_name) {
      
      diff_results[[tissue_name]] |>
        mutate(tissue = tissue_name)
    }
  )
)



## Identify significant phosphosites

significant_phospho <- diff_results_all |>
  filter(
    !is.na(padj),
    padj < 0.05
  )

# Increased phosphorylation
increased_phospho <- significant_phospho |>
  filter(coef_condition > 0)

# Decreased phosphorylation
decreased_phospho <- significant_phospho |>
  filter(coef_condition < 0)


## Number of significant phosphosites per tissue

significant_summary <- significant_phospho |>
  group_by(tissue) |>
  summarise(
    significant_sites = n(),
    increased_sites = sum(coef_condition > 0),
    decreased_sites = sum(coef_condition < 0),
    .groups = "drop"
  )

significant_summary

########################################################

## Create folder to start saving plots
dir.create("results/phosphoproteomics", recursive = TRUE, showWarnings = FALSE)

write_tsv(
  diff_results_all,
  "results/phosphoproteomics/06_differential_phospho_results.tsv"
)
########################################################

## Volcano plots

volcano_plot <- ggplot(
  diff_results_all,
  aes(
    x = coef_condition,
    y = -log10(padj),
    color = padj < 0.05
  )
) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ tissue) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey",
      "TRUE" = "red"
    )
  ) +
  theme_bw() +
  labs(
    title = "Differential phosphoproteomics",
    x = "Training vs Control coefficient",
    y = "-log10(FDR)",
    color = "Significant"
  )

# Save volcano plot
ggsave(
  "results/phosphoproteomics/volcano_plots.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300
)

########################################################

## Number of significant proteins per tissue

significant_proteins <- significant_phospho |>
  select(tissue, phosphosite_id) |>
  distinct() |>
  left_join(
    phospho_8weeks |>
      select(phosphosite_id, protein_id) |>
      distinct(),
    by = "phosphosite_id"
  )

significant_summary <- significant_summary |>
  left_join(
    significant_proteins |>
      group_by(tissue) |>
      summarise(
        proteins = n_distinct(protein_id),
        .groups = "drop"
      ),
    by = "tissue"
  )

## Plot of significant phosphosites

plot_significant <- ggplot(
  significant_summary,
  aes(x = tissue, y = significant_sites)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Tissue",
    y = "Significant phosphosites",
    title = "Significant phosphosites across tissues"
  ) +
  theme_minimal()

## Save plot

ggsave(
  "results/phosphoproteomics/significant_phosphosites.png",
  plot = plot_significant,
  width = 8,
  height = 6,
  dpi = 300
)

## Save summary table

write_tsv(
  significant_summary,
  "results/phosphoproteomics/significant_summary.tsv"
)

########################################################

## Tissue comparison

# Number of significant phosphosites per tissue
significant_phospho |>
  count(tissue, name = "n_significant")

# Recover protein information

significant_proteins <- significant_phospho |>
  distinct(tissue, phosphosite_id) |>
  left_join(
    phospho_8weeks |>
      distinct(tissue, phosphosite_id, protein_id, gene_symbol),
    by = c("tissue", "phosphosite_id")
  )

# Shared proteins across tissues

shared_proteins <- significant_proteins |>
  filter(!is.na(protein_id)) |>
  distinct(tissue, protein_id, gene_symbol) |>
  count(protein_id, gene_symbol, name = "n_tissues") |>
  arrange(desc(n_tissues))


# Proteins altered in more than one tissue
shared_proteins |>
  filter(n_tissues > 1)

# Plotting 
plot_shared <- ggplot(
  shared_proteins,
  aes(x = reorder(gene_symbol, n_tissues), y = n_tissues)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Protein",
    y = "Number of tissues",
    title = "Shared proteins across tissues"
  ) +
  theme_minimal()

ggsave(
  "results/phosphoproteomics/shared_proteins.png",
  plot = plot_shared,
  width = 8,
  height = 7,
  dpi = 300
)

###############################################
