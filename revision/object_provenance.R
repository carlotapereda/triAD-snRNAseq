# Extract package versions and pipeline parameters recorded inside Emouse.rds.
# Answers Reviewer 3 (software versions) and the 15 PCs / res 0.5 vs dims 1:30 / res 0.7 discrepancy.
# Run headless: Rscript code/new/object_provenance.R > results/csv/object_provenance.txt 2>&1

mem.maxVSize(vsize = 60000)
suppressPackageStartupMessages(library(Seurat))

obj <- readRDS("data/processed/Emouse.rds")
DefaultAssay(obj) <- "RNA"
if ("SCT" %in% Assays(obj)) obj[["SCT"]] <- NULL
gc()

cat("=== object version (Seurat version that created the object) ===\n")
print(obj@version)

cat("\n=== assays / reductions / graphs ===\n")
print(Assays(obj))
print(Reductions(obj))
print(Graphs(obj))
print(dim(obj))

cat("\n=== command log: names and timestamps ===\n")
cmds <- obj@commands
for (nm in names(cmds)) {
  cat(sprintf("%-45s %s\n", nm, format(cmds[[nm]]@time.stamp, "%Y-%m-%d %H:%M:%S")))
}

cat("\n=== command log: full parameters ===\n")
for (nm in names(cmds)) {
  cat("\n---", nm, "---\n")
  cat("call:", cmds[[nm]]@call.string, "\n")
  cat("assay:", cmds[[nm]]@assay.used, "\n")
  p <- cmds[[nm]]@params
  for (pn in names(p)) {
    v <- p[[pn]]
    if (is.atomic(v) && length(v) <= 40) {
      cat(sprintf("  %s = %s\n", pn, paste(format(v), collapse = ", ")))
    } else {
      cat(sprintf("  %s = <%s, length %d>\n", pn, class(v)[1], length(v)))
    }
  }
}

cat("\n=== dims actually used by FindNeighbors / RunUMAP / RunTSNE ===\n")
for (nm in grep("FindNeighbors|RunUMAP|RunTSNE|FindClusters|RunPCA", names(cmds), value = TRUE)) {
  p <- cmds[[nm]]@params
  d <- p[["dims"]]
  r <- p[["resolution"]]
  cat(sprintf("%-45s dims=%s resolution=%s\n", nm,
              if (is.null(d)) "NA" else paste0(min(d), ":", max(d), " (n=", length(d), ")"),
              if (is.null(r)) "NA" else paste(r, collapse = ",")))
}

cat("\n=== clustering-resolution columns present in metadata ===\n")
print(grep("res\\.|snn|clusters", colnames(obj@meta.data), value = TRUE))
cat("\nnlevels of seurat_clusters:", nlevels(factor(obj$seurat_clusters)), "\n")

cat("\n=== sessionInfo of THIS run (not the original) ===\n")
print(sessionInfo())
