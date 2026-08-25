#!/usr/bin/env Rscript
# Rebuild the stability figure from the saved cluster labels (no 19 GB reload).
# Replaces the misleading "consolidation" panel with NESTEDNESS: do finer-resolution
# clusters nest within coarser ones (i.e. higher resolution SUBDIVIDES rather than
# reshuffles) -- the correct evidence that the partition is stable.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
FIG  <- "results/figures/"
DEST <- "manuscript/figures/figures_reviews/"

labs <- readRDS(file.path(FIG, "resolution_sweep_labels.rds"))
res_vec <- c(0.3, 0.5, 0.7, 1.0)
col_of <- function(r) labs[[paste0("res_", sub("\\.?0+$", "", sprintf("%g", r)))]]
# robust column fetch (res_1 vs res_1.0)
getcol <- function(r) { nm <- paste0("res_", r); if (!nm %in% names(labs)) nm <- paste0("res_", sub("0$","",as.character(r))); labs[[paste0("res_", r)]] }
r_names <- c("0.3"="res_0.3","0.5"="res_0.5","0.7"="res_0.7","1"="res_1")

ARI <- function(a, b) { tab <- table(a, b); n <- sum(tab); ci <- function(x) sum(choose(x, 2))
  idx <- ci(tab); a_i <- ci(rowSums(tab)); b_j <- ci(colSums(tab))
  exp <- a_i*b_j/choose(n,2); mx <- (a_i+b_j)/2; (idx-exp)/(mx-exp) }

# weighted nestedness of a FINER partition within a COARSER one
nestedness <- function(fine, coarse) {
  tab <- table(fine, coarse); sizes <- rowSums(tab)
  frac <- apply(tab, 1, max) / sizes
  sum(frac * sizes) / sum(sizes)
}

n_clusters <- sapply(r_names, function(c) length(unique(labs[[c]])))
summ <- data.frame(resolution = factor(names(r_names), levels = names(r_names)),
                   n_clusters = as.integer(n_clusters))

adj <- data.frame(
  pair = factor(c("0.3 -> 0.5", "0.5 -> 0.7", "0.7 -> 1.0"),
                levels = c("0.3 -> 0.5","0.5 -> 0.7","0.7 -> 1.0")),
  ARI  = c(ARI(labs$res_0.3, labs$res_0.5), ARI(labs$res_0.5, labs$res_0.7), ARI(labs$res_0.7, labs$res_1)),
  nestedness = c(nestedness(labs$res_0.5, labs$res_0.3),
                 nestedness(labs$res_0.7, labs$res_0.5),
                 nestedness(labs$res_1,   labs$res_0.7)))
cat("=== stability metrics ===\n"); print(adj)
write.csv(adj, file.path(FIG, "resolution_stability_metrics.csv"), row.names = FALSE)

p1 <- ggplot(summ, aes(resolution, n_clusters)) +
  geom_col(fill = "grey60") + geom_text(aes(label = n_clusters), vjust = -0.4) +
  theme_bw() + labs(x = "Louvain resolution", y = "# clusters", title = "a  Clusters per resolution")
p2 <- ggplot(adj, aes(pair, ARI, group = 1)) + geom_line() + geom_point(size = 3) + ylim(0, 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = NULL, y = "adjacent-resolution ARI", title = "b  Partition stability (ARI)")
p3 <- ggplot(adj, aes(pair, nestedness, group = 1)) + geom_line(color = "steelblue") +
  geom_point(size = 3, color = "steelblue") + ylim(0, 1) +
  geom_text(aes(label = sprintf("%.2f", nestedness)), vjust = -1, size = 3.5) +
  theme_bw() + theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = NULL, y = "nestedness (finer within coarser)",
       title = "c  Higher resolution subdivides, not reshuffles")
p <- p1 | p2 | p3
ggsave(file.path(FIG, "resolution_stability-1.png"), p, width = 14, height = 4.5, dpi = 300)
ggsave(paste0(DEST, "SuppFig1n_clustering_resolution_stability.png"), p, width = 14, height = 4.5, dpi = 300)
cat("saved stability figure\n")
