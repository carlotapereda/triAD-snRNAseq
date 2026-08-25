#!/usr/bin/env Rscript
# Generates two review figures requested by the authors:
#  (1) a cell-type-labeled UMAP (with endothelial + mural populations labeled) -> Reviewer 1.3
#  (2) a PCA elbow plot justifying the choice of 15 PCs -> Reviewer 2.1 / Reviewer 3.1f
# Run AFTER the resolution sweep finishes (do not run two 19 GB loads at once on 24 GB RAM).
# Run:  Rscript code/new/umap_and_elbow.R    (from repo root)

mem.maxVSize(vsize = 60000)
options(future.globals.maxSize = 1e11)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

FIG  <- "results/figures/"
DEST <- "manuscript/figures/figures_reviews/"
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
save_both <- function(p, fig_name, supp_name, w, h) {
  ggsave(file.path(FIG, fig_name), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(paste0(DEST, supp_name), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  msg("saved", supp_name)
}

msg("loading object ...")
obj <- readRDS("data/processed/Emouse.rds")
DefaultAssay(obj) <- "RNA"
if ("SCT" %in% Assays(obj)) obj[["SCT"]] <- NULL
gc()

# readable labels for the celltype codes
lab_map <- c(ex.neu="Excitatory neurons", in.neu="Inhibitory neurons", ast="Astrocytes",
             oli="Oligodendrocytes", opc="OPC", mic="Microglia", end="Endothelial",
             per="Mural/Pericyte", negative="negative", hybrid="hybrid")
obj$celltype_lab <- factor(dplyr::recode(as.character(obj$celltype), !!!lab_map))
Idents(obj) <- "celltype_lab"
keep <- setdiff(levels(obj$celltype_lab), c("negative", "hybrid"))

## (1a) cell-type-labeled UMAP
p_umap <- DimPlot(obj, group.by = "celltype_lab", label = TRUE, repel = TRUE, label.size = 4,
                  raster = TRUE, raster.dpi = c(700, 700)) +
  ggtitle("Cell types (endothelial and mural populations labeled)") + NoLegend()
save_both(p_umap, "umap_celltypes_labeled-1.png", "SuppFig1_UMAP_cell_types_labeled.png", 8.5, 7)

## (1b) UMAP highlighting the vascular (endothelial + mural) nuclei
vasc <- WhichCells(obj, idents = c("Endothelial", "Mural/Pericyte"))
p_vasc <- DimPlot(obj, cells.highlight = vasc, cols.highlight = "firebrick",
                  sizes.highlight = 0.4, raster = TRUE, raster.dpi = c(700, 700)) +
  ggtitle(sprintf("Endothelial + mural nuclei on the UMAP (n = %d; %.1f%%)",
                  length(vasc), 100 * length(vasc) / ncol(obj))) + NoLegend()
save_both(p_vasc, "umap_vascular_highlight-1.png", "SuppFig1h_vascular_cells_on_UMAP.png", 7.5, 7)

## (2) PCA elbow plot to justify 15 PCs
sd  <- Stdev(obj, reduction = "pca")
ndm <- length(sd)
varexp <- sd^2 / sum(sd^2) * 100
edf <- data.frame(PC = seq_len(ndm), pct_var = varexp)
p_elbow <- ggplot(edf, aes(PC, pct_var)) +
  geom_line(color = "grey40") + geom_point(size = 1.6) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "firebrick") +
  annotate("text", x = 15.5, y = max(varexp) * 0.9, hjust = 0, color = "firebrick",
           label = "15 PCs used") +
  theme_bw() +
  labs(x = "principal component", y = "% variance explained",
       title = "PCA elbow: variance explained plateaus by ~15 PCs")
save_both(p_elbow, "pca_elbow-1.png", "SuppFig1o_PCA_elbow_15PCs.png", 7, 4.5)

# report the numbers so the rationale sentence can cite them
cat(sprintf("\nPCA dims available: %d\n", ndm))
cat(sprintf("cumulative %% variance at 10/15/20/30 PCs: %.1f / %.1f / %.1f / %.1f\n",
            sum(varexp[1:10]), sum(varexp[1:15]), sum(varexp[1:min(20,ndm)]),
            sum(varexp[1:min(30,ndm)])))
cat(sprintf("marginal %% variance at PC15 vs PC16 vs PC20: %.3f / %.3f / %.3f\n",
            varexp[15], varexp[min(16,ndm)], varexp[min(20,ndm)]))
msg("DONE")
