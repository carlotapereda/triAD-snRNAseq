# =============================================================================
# Step 4 of 6: Cell type abundance
#
# Purpose: compute per-library cell type proportions (hybrid and negative
#   nuclei removed) and test genotype x sex x age effects per cell type with
#   a three-way Type III ANOVA (car::Anova), BH-adjusted within cell type.
# Inputs:  annotated Seurat object "Emouse" (metadata: celltype,
#   cell_type_gen, orig.ident, sex, genotype, age).
# Outputs: per-library proportion table (df3);
#   2025oct_celltype_3anova_results.csv (ANOVA results).
# Next:    05_cluster_sample_abundance.R
# =============================================================================

### cell type abundance analysis

## remove doublets and unresolved cells
Idents(Emouse) <- "cell_type_gen"
Emouse <- subset(Emouse, idents = c("hybrid", "negative"), invert = TRUE)

## extract data needed
# get total count of per sample
Emouse$age_sex_geno <- paste(Emouse$sex, Emouse$genotype, Emouse$age, sep = "_")
df <- as.data.frame(table(Emouse$orig.ident, Emouse$age_sex_geno ))
df <- df[df$Freq != 0, ]
colnames(df) <- c("id", "group", "sum")
# get count per sample per cell type
df2 <- as.data.frame(table(Emouse$celltype, Emouse$orig.ident ))
colnames(df2) <- c("celltype", "id", "freq")
# combine count
df3 <- left_join(df2, df, by = "id")
# calculate prop
df3$prop <- df3$freq/df3$sum


## anova test
# Packages
library(dplyr)
library(tidyr)
library(car)      # for Type III Anova
library(purrr)

# 0) Prep: split "group" into age / sex / genotype and factorize
df_anova <- df3 %>%
  separate(group, into = c( "sex", "genotype", "age"), sep = "_", remove = FALSE) %>%
  mutate(
    celltype = factor(celltype),
    age      = factor(age),
    sex      = factor(sex),
    genotype = factor(genotype)
  )

# Type III needs sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

# 1) Run 3-way ANOVA PER celltype
celltypes <- levels(df_anova$celltype)

res_list <- lapply(celltypes, function(ct) {
  dat <- filter(df_anova, celltype == ct)
  # model: prop ~ age * sex * genotype
  m <- lm(prop ~ age * sex * genotype, data = dat)
  a3 <- car::Anova(m, type = 3)
  out <- as.data.frame(a3)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out$celltype <- ct
  out[, c("celltype", "term", colnames(out)[1:4])]  # keep key columns
})

results_per_celltype <- bind_rows(res_list)
results_per_celltype

#filter for plotting
results_per_celltype <- na.omit(results_per_celltype)

results_per_celltype <- results_per_celltype %>%
  filter(term != "(Intercept)")
results_per_celltype <- results_per_celltype %>%
    group_by(celltype) %>%
    mutate(p_adj = p.adjust(`Pr(>F)`, method = "BH"))
write.csv(results_per_celltype, "2025oct_celltype_3anova_results.csv")



