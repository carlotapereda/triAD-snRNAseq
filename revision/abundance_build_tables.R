#!/usr/bin/env Rscript
# Build sample x cell-type count tables and a per-sample covariate table for the
# compositional abundance analyses (scCODA in Python; crumblr/propeller in R).
# Works from the saved metadata (no 19 GB reload).

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
OUT <- "results/abundance"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

md <- readRDS("results/figures/cell_metadata.rds")
md$age <- factor(md$age, levels = c("06Mo","12Mo","18Mo"))

# per-sample covariates: primary (genotype/sex/age) + nuisance (depth, complexity, %mt)
# NOTE on the mitochondrial covariate: in this single-nucleus data 58% of nuclei have percent.mt
# exactly 0, so the per-library MEDIAN is 0 for 37 of 47 libraries and takes only 11 distinct values.
# Z-scored, that is a sparse indicator of ten arbitrary libraries, not a measure of tissue quality,
# and it competed with sex for the dentate gyrus signal in the compositional model. Use the per-
# library MEAN instead: 47 distinct values, range 0.018 to 0.045. The median is kept for provenance.
cov <- md %>%
  group_by(orig.ident, genotype, sex, age) %>%
  summarise(n_cells = n(),
            med_nCount = median(nCount_RNA),
            med_nFeature = median(nFeature_RNA),
            mean_percent_mt = mean(percent.mt),
            med_percent_mt = median(percent.mt),
            med_log10GPU = median(log10GenesPerUMI), .groups = "drop")
write.csv(cov, file.path(OUT, "sample_covariates.csv"), row.names = FALSE)
cat("samples:", nrow(cov), "\n")

# count matrices at two annotation levels used by the manuscript (drop 'remove')
make_counts <- function(col) {
  d <- md[!md[[col]] %in% c("remove","negative","hybrid"), ]
  tab <- as.data.frame.matrix(table(d$orig.ident, d[[col]]))
  tab$orig.ident <- rownames(tab)
  tab[, c("orig.ident", setdiff(colnames(tab), "orig.ident"))]
}
# third level: the Louvain clusters behind Fig. 2c / Extended Data Fig. 2c. Prefixed so the column
# names are valid identifiers in the patsy formulas scCODA builds. Cluster 36 (848 nuclei, 0.15%) is
# annotated "remove" in cell_type_identity for every one of its nuclei, i.e. the authors excluded it,
# so it is dropped here too and 35 clusters are analyzed. Without this the cluster count table would
# disagree with percell_labels.csv, which filters on cell_type_identity.
md$cluster_id <- paste0("C", md$seurat_clusters)
md$cluster_id[md$cell_type_identity %in% "remove"] <- "remove"

major   <- make_counts("cell_type_identity")   # Ex/In neurons, Oligo, Astro, OPC, Microglia, ChoroidPlexus
fine    <- make_counts("hippocampus_cell_type") # CA1/CA3/DG/Subiculum/Interneurons + glia
cluster <- make_counts("cluster_id")            # 36 Louvain clusters at resolution 0.5
write.csv(major,   file.path(OUT, "counts_major.csv"),   row.names = FALSE)
write.csv(fine,    file.path(OUT, "counts_fine.csv"),    row.names = FALSE)
write.csv(cluster, file.path(OUT, "counts_cluster.csv"), row.names = FALSE)
cat("major cell types:", paste(setdiff(colnames(major),"orig.ident"), collapse=", "), "\n")
cat("fine cell types:", paste(setdiff(colnames(fine),"orig.ident"), collapse=", "), "\n")
cat("clusters:", ncol(cluster) - 1, "\n")

# also a long per-cell table (sample + celltype) for scCODA from_scanpy convenience
write.csv(md[!md$cell_type_identity %in% c("remove"),
             c("orig.ident","genotype","sex","age","cell_type_identity","hippocampus_cell_type",
               "cluster_id")],
          file.path(OUT, "percell_labels.csv"), row.names = FALSE)

# which compositional packages are available?
cat("\n=== R package availability ===\n")
for (p in c("crumblr","variancePartition","speckle","DirichletReg","edgeR","limma","lme4","glmmTMB"))
  cat(sprintf("%-18s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("DONE\n")
