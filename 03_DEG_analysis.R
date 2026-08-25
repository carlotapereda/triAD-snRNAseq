# =============================================================================
# Step 3 of 6: Differential gene expression (APOE4/4 vs APOE3/3)
#
# Purpose: within each of the 6 sex-by-age groups and each major cell type,
#   test APOE4/4 vs APOE3/3 with FindMarkers (MAST, logfc.threshold = 0.01).
# Inputs:  annotated Seurat object "Emouse" from 02_cell_type_annotation.R
#   (uses metadata: sex, age, genotype, cell_type_ident).
# Outputs: combined DEG table DEG.list.all (all tested genes, per comparison
#   and cell type); filtering to the significant set (p_val_adj < 0.05 and
#   abs avg_log2FC >= 0.1) and export are in requested_code.R (DEG.csv,
#   sDEG.csv).
# Next:    04_cell_type_abundance.R
# =============================================================================

## identify DEGs and stratify by age and sex

Emouse$comparison <- paste(Emouse$sex, Emouse$age, sep = "_")

comparison <- c(
  "Female_06Mo", "Female_12Mo", "Female_18Mo",
  "Male_06Mo", "Male_12Mo", "Male_18Mo"
)

Celltype <- c(
  "Excitatory.Neurons",
  "Inhibitory.Neurons",
  "Astrocytes",
  "Microglia",
  "Oligodendrocytes",
  "Oligodendrocyte.Precursor"
)

library(Seurat)
library(stringr)

obj <- Emouse
DEG.list.all <- NULL

for (key in comparison) {

  # subset by comparison
  Idents(obj) <- "comparison"
  sub.comp <- subset(obj, idents = key)

  for (ct in Celltype) {

    # subset by cell type
    Idents(sub.comp) <- "cell_type_ident"
    sub.ct <- subset(sub.comp, idents = ct)

    # DEG
    DEG.list <- FindMarkers(
      sub.ct,
      ident.1 = "E44",
      ident.2 = "E33",
      verbose = TRUE,
      test.use = "MAST",
      logfc.threshold = 0.01
    )

    if (nrow(DEG.list) == 0) {
      message("No DEG found for ", key, " / ", ct)
      next
    }

    DEG.list$gene <- rownames(DEG.list)

    # Seurat version compatibility
    if ("avg_log2FC" %in% colnames(DEG.list)) {
      DEG.list$dir <- ifelse(DEG.list$avg_log2FC < 0, "neg", "pos")
    } else if ("avg_logFC" %in% colnames(DEG.list)) {
      DEG.list$dir <- ifelse(DEG.list$avg_logFC < 0, "neg", "pos")
    } else {
      DEG.list$dir <- NA
    }

    DEG.list$comparison <- key
    DEG.list$celltype <- ct

    DEG.list.all <- rbind(DEG.list.all, DEG.list)
  }
}
