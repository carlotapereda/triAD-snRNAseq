#!/usr/bin/env Rscript
# Export the per-nucleus metadata (small) so all abundance modeling can run without
# reloading the 19 GB object. Prints candidate cell-type columns for inspection.
mem.maxVSize(vsize = 60000)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
suppressPackageStartupMessages(library(Seurat))

msg("loading object ...")
obj <- readRDS("data/processed/Emouse.rds")
md <- obj@meta.data
msg("meta.data:", nrow(md), "rows x", ncol(md), "cols")
cat("\n=== ALL COLUMN NAMES ===\n"); print(colnames(md))

# show tables of every column that looks like an annotation / cluster / cell type
cand <- grep("type|cluster|celltype|ident|annot|class|subclass|seurat|leiden|louvain|hippo",
             colnames(md), ignore.case = TRUE, value = TRUE)
cat("\n=== candidate annotation columns ===\n"); print(cand)
for (cc in cand) {
  cat("\n---", cc, "--- (", length(unique(md[[cc]])), "levels )\n")
  if (length(unique(md[[cc]])) <= 60) print(sort(table(md[[cc]]), decreasing = TRUE))
}

# nuisance covariates present?
cat("\n=== numeric QC columns (nuisance covariates) ===\n")
num <- names(md)[sapply(md, is.numeric)]
print(num)

# drop any list columns, save compact
md2 <- md[, !sapply(md, is.list)]
saveRDS(md2, "results/figures/cell_metadata.rds")
data.table_ok <- requireNamespace("data.table", quietly = TRUE)
write.csv(md2, "results/figures/cell_metadata.csv", row.names = TRUE)
msg("saved results/figures/cell_metadata.rds + .csv (", ncol(md2), "cols )")
msg("DONE")
