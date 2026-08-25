# =============================================================================
# Step 6 of 6: APOE expression level comparison
#
# Purpose: summarize hApoE transgene expression (SCT "data" slot, Pearson
#   residuals) per genotype-sex-age group: number and percentage of positive
#   cells (residual > 0), mean over all cells, and mean over positive cells;
#   pairwise group tests (Wilcoxon for expression, proportion test for
#   percent positive), BH-adjusted.
# Inputs:  annotated Seurat object "obj"/"Emouse" with SCT assay and
#   age_sex_geno metadata.
# Outputs: per-group summary (avg_pos_df);
#   APOEexpr_pairwise_pval.csv (pairwise adjusted p-values).
# =============================================================================

DefaultAssay(obj) <- "SCT"

library(dplyr)

group_var  <- "age_sex_geno"
gene       <- "hapoE-transgene"
assay_use  <- "SCT"     # or "RNA"
slot_use   <- "data"    # if "SCT", residuals can be negative
thresh     <- 0         # positivity threshold for the chosen slot

# ---- Build per-cell data.frame ---------------------------------------------
expr_vec <- as.numeric(GetAssayData(obj, assay = assay_use, slot = slot_use)[gene, ])
pos_vec  <- expr_vec > thresh

grp <- obj[[group_var]][, 1]

df  <- data.frame(
  group = grp,
  expr  = expr_vec,
  pos   = pos_vec
)

# ---- Per-group summaries (your original logic) -----------------------------
avg_pos_df <- df %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n_total = sum(!is.na(expr)),                     # count rows without NA
    n_pos   = sum(pos %in% TRUE, na.rm = TRUE),      # number of positives
    pct_pos = 100 * n_pos / n_total,
    avg_all = mean(expr, na.rm = TRUE),              # mean over all cells
    avg_pos = ifelse(
      n_pos > 0,
      mean(expr[which(pos)], na.rm = TRUE),          # mean among positives only
      NA_real_
    ),
    percent_pos = n_pos / n_total,
    .groups = "drop"
  )

# --------------------------------------------------------------------
# Pairwise p-values for:
#   1) avg_all      -> Wilcoxon on expr (all cells)
#   2) avg_pos      -> Wilcoxon on expr (only pos cells)
#   3) percent_pos  -> prop.test on pos vs total
# --------------------------------------------------------------------

# Clean per-cell data (no NA expr, group as factor)
df_non_na <- df %>%
  filter(!is.na(expr)) %>%
  mutate(group = factor(group))

## 1) avg_all: pairwise Wilcoxon tests on expr across all cells ----------
pw_all <- pairwise.wilcox.test(
  x = df_non_na$expr,
  g = df_non_na$group,
  p.adjust.method = "BH"
)

pairwise_avg_all <- as.data.frame(as.table(pw_all$p.value)) %>%
  filter(!is.na(Freq)) %>%
  rename(
    group1 = Var1,
    group2 = Var2,
    p_adj  = Freq
  ) %>%
  mutate(metric = "avg_all")

## 2) avg_pos: pairwise Wilcoxon only among positive cells ---------------
df_pos <- df_non_na %>% filter(pos)

# only run if there are at least 2 groups with positive cells
if (length(unique(df_pos$group)) >= 2) {
  pw_pos <- pairwise.wilcox.test(
    x = df_pos$expr,
    g = df_pos$group,
    p.adjust.method = "BH"
  )

  pairwise_avg_pos <- as.data.frame(as.table(pw_pos$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(
      group1 = Var1,
      group2 = Var2,
      p_adj  = Freq
    ) %>%
    mutate(metric = "avg_pos")
} else {
  pairwise_avg_pos <- NULL
}

## 3) percent_pos: pairwise proportion test on n_pos / n_total ----------
x <- avg_pos_df$n_pos
n <- avg_pos_df$n_total
names(x) <- avg_pos_df$group  # optional but helps readability

pw_prop <- pairwise.prop.test(
  x = x,
  n = n,
  p.adjust.method = "BH"
)

pairwise_percent_pos <- as.data.frame(as.table(pw_prop$p.value)) %>%
  filter(!is.na(Freq)) %>%
  rename(
    group1 = Var1,
    group2 = Var2,
    p_adj  = Freq
  ) %>%
  mutate(metric = "percent_pos")

# ---- Combined table of all metrics -------------------------------------
pairwise_all_metrics <- bind_rows(
  pairwise_avg_all,
  pairwise_avg_pos,
  pairwise_percent_pos
)

# Objects you probably care about:
avg_pos_df             # per-group summary
pairwise_avg_all       # pairwise p for avg_all
pairwise_avg_pos       # pairwise p for avg_pos (may be NULL if no pos cells)
pairwise_percent_pos   # pairwise p for percent_pos
pairwise_all_metrics   # all three metrics stacked


write.csv(pairwise_all_metrics, "APOEexpr_pairwise_pval.csv")
