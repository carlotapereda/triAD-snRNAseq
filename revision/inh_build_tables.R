#!/usr/bin/env Rscript
# Build inhibitory-subtype count + covariate tables for compositional analysis (Reviewer 2:
# does inhibitory-neuron subtype composition differ by sex/age/genotype?). Uses the 12 subclusters.
suppressPackageStartupMessages({ library(dplyr) })
OUT <- "results/abundance"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

md <- readRDS(file.path(OUT, "inh_metadata.rds"))
md$age <- factor(md$age, levels = c("06Mo","12Mo","18Mo"))
md$subtype <- paste0("In", md$seurat_clusters)   # In0..In11

# which SCT column? check 'cluster' vs seurat_clusters agreement, print markers hint
cat("subclusters (seurat_clusters):\n"); print(sort(table(md$subtype), decreasing = TRUE))
if ("cluster" %in% colnames(md)) { cat("\n'cluster' column levels:\n"); print(sort(table(md$cluster), decreasing = TRUE)) }

# per-sample covariates computed WITHIN the inhibitory compartment
cov <- md %>% group_by(orig.ident, genotype, sex, age) %>%
  summarise(n_cells = n(), med_nCount = median(nCount_RNA),
            mean_percent_mt = mean(percent.mt),   # median is degenerate, see abundance_build_tables.R
            med_percent_mt = median(percent.mt), .groups = "drop")
write.csv(cov, file.path(OUT, "sample_covariates_inh.csv"), row.names = FALSE)

counts <- as.data.frame.matrix(table(md$orig.ident, md$subtype))
counts$orig.ident <- rownames(counts)
counts <- counts[, c("orig.ident", setdiff(colnames(counts), "orig.ident"))]
write.csv(counts, file.path(OUT, "counts_inh.csv"), row.names = FALSE)

write.csv(md[, c("orig.ident","genotype","sex","age","subtype")],
          file.path(OUT, "percell_labels_inh.csv"), row.names = FALSE)
cat("\nsamples:", nrow(cov), "| subtypes:", ncol(counts)-1, "\n")
cat("DONE\n")
