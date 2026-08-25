# =============================================================================
# Step 5 of 6: Cluster-level abundance
#
# Purpose: same design as step 4 but at the level of the 36 Louvain clusters:
#   per-library cluster proportions, then a three-way Type III ANOVA
#   (genotype x sex x age) per cluster, BH-adjusted within cluster.
# Inputs:  annotated Seurat object "Emouse" (metadata: cluster, orig.ident,
#   sex, genotype, age).
# Outputs: per-library cluster proportion table (df3);
#   cluster_3anova_results.csv (ANOVA results).
# Next:    06_APOE_expression_level_comparison.R
# =============================================================================

Emouse$age_sex_geno <- paste(Emouse$sex, Emouse$genotype, Emouse$age, sep = "_")
df <- as.data.frame(table(Emouse$orig.ident, Emouse$age_sex_geno ))
df <- df[df$Freq != 0, ]
colnames(df) <- c("id", "group", "sum")

df2 <- as.data.frame(table(Emouse$cluster, Emouse$orig.ident ))
colnames(df2) <- c("celltype", "id", "freq")

df3 <- left_join(df2, df, by = "id")

df3$prop <- df3$freq/df3$sum




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

write.csv(results_per_celltype, "cluster_3anova_results.csv")
