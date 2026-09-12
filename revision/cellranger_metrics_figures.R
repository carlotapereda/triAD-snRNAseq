# Cell Ranger per-library sequencing and mapping metrics across conditions
# (Reviewer 3, Major point 1: "plots of RIN, number of reads, ... across conditions").
#
# RIN was not measured: nuclei were loaded intact into the 10x Chromium and no bulk RNA was
# extracted, so there is no RNA integrity number for these libraries. The Cell Ranger count
# web summaries (cellranger-7.0.0, one per library, results/cellranger/ and
# results/csv/cellranger_metrics_per_library.csv) are the library-level record of sequencing
# and library quality and are shown here instead.
#
# Inputs : results/csv/cellranger_metrics_per_library.csv (parsed from the 47 web_summary.html)
#          results/csv/QC_per_library_summary.csv       (from qc_condition_figures.R)
# Outputs: results/csv/cellranger_metrics_anova.csv
#          results/csv/cellranger_metrics_condition_summary.csv
#          SuppFig1q_cellranger_metrics_by_condition.png (Extended Data Fig. 1q after the
#          2026-09-04 panel cuts; seven facets incl. total read pairs)
# Statistical unit = library (n = 47): three-way ANOVA (genotype x sex x age).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

cr <- read.csv("results/csv/cellranger_metrics_per_library.csv", check.names = FALSE)
qc <- read.csv("results/csv/QC_per_library_summary.csv")

lib <- cr %>%
  inner_join(qc %>% select(orig.ident, nuclei_retained = nuclei), by = "orig.ident") %>%
  mutate(
    age      = factor(age, levels = c("06Mo", "12Mo", "18Mo"),
                      labels = c("6 months", "12 months", "18 months")),
    genotype = factor(genotype, levels = c("E33", "E44"),
                      labels = c("APOE3/3", "APOE4/4")),
    sex      = factor(sex, levels = c("Female", "Male")),
    reads_M  = Number_of_Reads / 1e6,
    pct_retained = 100 * nuclei_retained / Estimated_Number_of_Cells
  )
stopifnot(nrow(lib) == 47)

geno_cols <- c(`APOE3/3` = "#0072B2", `APOE4/4` = "#D55E00")
out_rev <- "manuscript/figures/figures_reviews"; out_res <- "results/figures"
save_fig <- function(p, name, w, h) for (d in c(out_rev, out_res))
  ggsave(file.path(d, name), p, width = w, height = h, dpi = 300, bg = "white")

metrics <- c(
  reads_M                                        = "Total read pairs (millions)",
  Mean_Reads_per_Cell                            = "Mean reads per nucleus",
  Sequencing_Saturation                          = "Sequencing saturation (%)",
  Reads_Mapped_Confidently_to_Genome             = "Reads mapped confidently to genome (%)",
  Reads_Mapped_Confidently_to_Intronic_Regions   = "Reads mapped to intronic regions (%)",
  Fraction_Reads_in_Cells                        = "Fraction of reads in nuclei (%)",
  Q30_Bases_in_RNA_Read                          = "Q30 bases in RNA read (%)",
  Estimated_Number_of_Cells                      = "Nuclei called by Cell Ranger",
  Median_UMI_Counts_per_Cell                     = "Median UMIs per nucleus (Cell Ranger)"
)

## ---- three-way ANOVA per metric ------------------------------------------
anova_one <- function(metric) {
  fit <- lm(reformulate("genotype * sex * age", metric), data = lib)
  tab <- as.data.frame(anova(fit)); tab$term <- rownames(tab); tab$metric <- metric
  tab %>% filter(term != "Residuals") %>% transmute(metric, term, F = `F value`, p = `Pr(>F)`)
}
aov_tab <- bind_rows(lapply(names(metrics), anova_one))
write.csv(aov_tab, "results/csv/cellranger_metrics_anova.csv", row.names = FALSE)
cat("\nThree-way ANOVA p-values (library = unit, n = 47):\n")
print(aov_tab %>% select(-F) %>% pivot_wider(names_from = term, values_from = p) %>% as.data.frame(), digits = 2)

fp <- function(p) if (p < 0.01) sprintf("%.3f", p) else sprintf("%.2f", p)
p_line <- function(metric) {
  a <- aov_tab %>% filter(metric == !!metric)
  sprintf("genotype p = %s, sex p = %s, age p = %s; interactions p >= %s",
          fp(a$p[a$term == "genotype"]), fp(a$p[a$term == "sex"]), fp(a$p[a$term == "age"]),
          fp(min(a$p[grepl(":", a$term)])))
}

## ---- condition summary table ---------------------------------------------
cond <- lib %>% group_by(genotype, sex, age) %>%
  summarise(n_libraries = n(),
            mean_reads_M = round(mean(reads_M), 1),
            mean_reads_per_nucleus = round(mean(Mean_Reads_per_Cell)),
            mean_saturation = round(mean(Sequencing_Saturation), 1),
            mean_pct_mapped_confidently = round(mean(Reads_Mapped_Confidently_to_Genome), 1),
            mean_pct_intronic = round(mean(Reads_Mapped_Confidently_to_Intronic_Regions), 1),
            mean_fraction_reads_in_cells = round(mean(Fraction_Reads_in_Cells), 1),
            mean_cells_called = round(mean(Estimated_Number_of_Cells)),
            .groups = "drop")
write.csv(cond, "results/csv/cellranger_metrics_condition_summary.csv", row.names = FALSE)

base_thm <- theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, color = "grey30"),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

## ---- q (final letter): seven library-level metrics by condition ---------------------------
show <- c("reads_M", "Mean_Reads_per_Cell", "Sequencing_Saturation", "Reads_Mapped_Confidently_to_Genome",
          "Reads_Mapped_Confidently_to_Intronic_Regions", "Fraction_Reads_in_Cells",
          "Q30_Bases_in_RNA_Read")
long <- lib %>% select(orig.ident, genotype, sex, age, all_of(show)) %>%
  pivot_longer(all_of(show), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = show,
                         labels = sapply(show, function(m) {
                           paste0(metrics[[m]], "\n", p_line(m))
                         })))
pw <- ggplot(long, aes(sex, value, fill = genotype)) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA,
               linewidth = 0.3, width = 0.7) +
  geom_point(aes(shape = age),
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.15),
             size = 1.5, stroke = 0.4) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = geno_cols, name = NULL) +
  scale_shape_manual(values = c(`6 months` = 21, `12 months` = 24, `18 months` = 22), name = NULL) +
  guides(shape = guide_legend(override.aes = list(fill = "grey60"))) +
  base_thm + theme(strip.text = element_text(size = 6.8)) +
  labs(x = NULL, y = NULL,
       title = "Cell Ranger library metrics across conditions (each point = one library, n = 47)")
save_fig(pw, "SuppFig1q_cellranger_metrics_by_condition.png", 14, 6.4)

cat("\nOverall numbers for the response letter:\n")
rng <- function(x, d = 0) sprintf("median %s (range %s to %s)", format(round(median(x), d), big.mark = ","),
                                  format(round(min(x), d), big.mark = ","), format(round(max(x), d), big.mark = ","))
cat("read pairs (M):", rng(lib$reads_M, 1), "\n")
cat("mean reads/nucleus:", rng(lib$Mean_Reads_per_Cell), "\n")
cat("saturation %:", rng(lib$Sequencing_Saturation, 1), "\n")
cat("mapped confidently to genome %:", rng(lib$Reads_Mapped_Confidently_to_Genome, 1), "\n")
cat("intronic %:", rng(lib$Reads_Mapped_Confidently_to_Intronic_Regions, 1), "\n")
cat("exonic %:", rng(lib$Reads_Mapped_Confidently_to_Exonic_Regions, 1), "\n")
cat("fraction reads in cells %:", rng(lib$Fraction_Reads_in_Cells, 1), "\n")
cat("Q30 RNA read %:", rng(lib$Q30_Bases_in_RNA_Read, 1), "\n")
cat("cells called:", rng(lib$Estimated_Number_of_Cells), "\n")
cat("nuclei retained after QC, % of called:", rng(lib$pct_retained, 1), "\n")
cat("libraries with 'Low Fraction Reads in Cells' flag (< 70%):", sum(lib$Fraction_Reads_in_Cells < 70), "\n")
cat("antisense %:", rng(lib$Reads_Mapped_Antisense_to_Gene, 1), "\n")
