#!/usr/bin/env Rscript
# Reviewer 2.1 clustering-stability: recluster at several resolutions, cross-tabulate
# against the annotated cell-type labels, and quantify stability with:
#   (a) cross-resolution ARI  -> is the PARTITION stable as resolution changes (the real
#       stability metric; independent of the noisy per-cell module-score labels), and
#   (b) per-cell-type consolidation -> does each annotated cell type collapse into a
#       dedicated cluster at every resolution.
# NOTE: ARI of a clustering vs the per-cell celltype labels is intentionally NOT used as the
# headline metric: those labels are per-nucleus module-score argmax calls, so their noise (and
# the excitatory-neuron continuum) would masquerade as clustering instability.
# Run:  Rscript code/new/resolution_sweep.R   (from repo root)

mem.maxVSize(vsize = 60000)
options(future.globals.maxSize = 1e11)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })

FIG  <- "results/figures/"
DEST <- "manuscript/figures/figures_reviews/"
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)

ARI <- function(a, b) {                    # adjusted Rand index from a contingency table
  tab <- table(a, b); n <- sum(tab); ci <- function(x) sum(choose(x, 2))
  idx <- ci(tab); a_i <- ci(rowSums(tab)); b_j <- ci(colSums(tab))
  exp <- a_i * b_j / choose(n, 2); mx <- (a_i + b_j) / 2
  (idx - exp) / (mx - exp)
}

msg("loading object ...")
obj <- readRDS("data/processed/Emouse.rds")
DefaultAssay(obj) <- "RNA"
if ("SCT" %in% Assays(obj)) obj[["SCT"]] <- NULL
gc()
stopifnot("pca" %in% Reductions(obj))
ct <- as.character(obj$celltype)

msg("building neighbor graph on 15 PCs ...")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15,
                     graph.name = c("rev_nn", "rev_snn"), verbose = FALSE)

res_vec <- c(0.3, 0.5, 0.7, 1.0)
labs <- data.frame(celltype = ct)                    # collect per-cell cluster labels
summ <- data.frame()
sink(file.path(FIG, "resolution_sweep_crosstab.txt"))
for (r in res_vec) {
  obj <- FindClusters(obj, graph.name = "rev_snn", resolution = r, verbose = FALSE)
  cl  <- as.character(Idents(obj))
  labs[[paste0("res_", r)]] <- cl
  tab <- table(celltype = ct, cluster = cl)
  cat(sprintf("\n=================== resolution %.1f : %d clusters ===================\n", r, ncol(tab)))
  print(tab)
  summ <- rbind(summ, data.frame(resolution = r, n_clusters = ncol(tab)))
}
sink()
saveRDS(labs, file.path(FIG, "resolution_sweep_labels.rds"))   # so no reload is ever needed

## (a) cross-resolution ARI: each resolution vs the chosen 0.5, and adjacent pairs
ref <- "res_0.5"
summ$ARI_vs_0.5 <- sapply(paste0("res_", res_vec), function(c) ARI(labs[[c]], labs[[ref]]))
adj <- data.frame(pair = c("0.3 vs 0.5", "0.5 vs 0.7", "0.7 vs 1.0"),
                  ARI  = c(ARI(labs$res_0.3, labs$res_0.5),
                           ARI(labs$res_0.5, labs$res_0.7),
                           ARI(labs$res_0.7, labs$res_1)))

## (b) per-cell-type consolidation: fraction of each cell type's nuclei in its single
## dominant cluster at each resolution (discrete types -> near 1; neuronal continuum -> lower)
major <- c("oli","opc","mic","ast","in.neu","ex.neu","end","per")
cons <- do.call(rbind, lapply(res_vec, function(r) {
  cl <- labs[[paste0("res_", r)]]
  do.call(rbind, lapply(major, function(k) {
    idx <- ct == k; tt <- table(cl[idx])
    data.frame(resolution = r, celltype = k, consolidation = max(tt) / sum(tt))
  }))
}))

write.csv(summ, file.path(FIG, "resolution_sweep_summary.csv"), row.names = FALSE)
sink(file.path(FIG, "resolution_sweep_crosstab.txt"), append = TRUE)
cat("\n\n===== cross-resolution ARI (partition stability) =====\n"); print(summ); print(adj)
cat("\n===== per-cell-type consolidation (fraction in dominant cluster) =====\n")
print(tidyr::pivot_wider(cons, names_from = resolution, values_from = consolidation))
sink()
print(summ); print(adj)

## Figure
p1 <- ggplot(summ, aes(factor(resolution), n_clusters)) +
  geom_col(fill = "grey60") + geom_text(aes(label = n_clusters), vjust = -0.4) +
  theme_bw() + labs(x = "Louvain resolution", y = "# clusters", title = "a  Clusters per resolution")
p2 <- ggplot(adj, aes(pair, ARI, group = 1)) + geom_line() + geom_point(size = 3) +
  ylim(0, 1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = "adjacent-resolution ARI",
       title = "b  Partition stability across resolutions")
p3 <- ggplot(cons, aes(factor(resolution), consolidation, color = celltype, group = celltype)) +
  geom_line() + geom_point() + ylim(0, 1) + theme_bw() +
  labs(x = "Louvain resolution", y = "fraction in dominant cluster",
       title = "c  Each cell type maps to a dedicated cluster")
p <- p1 | p2 | p3
ggsave(file.path(FIG, "resolution_stability-1.png"), p, width = 15, height = 4.5, dpi = 300)
ggsave(paste0(DEST, "Supp Fig 1 _ clustering resolution stability.png"), p, width = 15, height = 4.5, dpi = 300)
msg("saved resolution_stability-1.png + Supp Fig 1 copy")
msg("DONE")
