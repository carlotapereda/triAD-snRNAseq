#!/usr/bin/env Rscript
# Memory-frugal runnable version of subtype_and_oligo_validation.Rmd
# Run headless:  Rscript code/new/subtype_and_oligo_validation.R
# Key fixes for the 24 GB machine:
#   - raise R's vector-memory ceiling at RUNTIME (env var R_MAX_VSIZE is ignored once R is running)
#   - drop the SCT assay's dense scale.data right after load (saves several GB)
#   - process one analysis at a time, gc() between, save each figure to disk

mem.maxVSize(vsize = 60000)          # MB; allow spill to swap if needed
options(future.globals.maxSize = 1e11)

t0 <- Sys.time()
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})

FIG <- "results/figures/"                       # run from repo root
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  msg("saved", name)
}

## ---- load & slim ---------------------------------------------------------
msg("reading Emouse.rds (~19 GB, be patient) ...")
obj <- readRDS("data/processed/Emouse.rds")
DefaultAssay(obj) <- "RNA"
msg("loaded. assays:", paste(Assays(obj), collapse = ", "),
    "| cells:", ncol(obj), "| RNA genes:", nrow(obj))

# free the heavy SCT scale.data (and SCT entirely; we only need RNA + reductions)
if ("SCT" %in% Assays(obj)) { obj[["SCT"]] <- NULL; msg("dropped SCT assay") }
gc()

obj$age   <- factor(obj$age, levels = c("06Mo", "12Mo", "18Mo"))
obj$group <- paste(obj$genotype, obj$sex, obj$age, sep = "_")
Idents(obj) <- "celltype"
gs   <- rownames(obj)
have <- function(g) g[g %in% gs]
real_ct <- setdiff(levels(factor(obj$celltype)), c("negative", "hybrid"))

# make sure RNA has log-normalized data for scoring/plots
if (max(GetAssayData(obj, "RNA", "data")[1:50, 1:50]) == 0 ||
    identical(GetAssayData(obj,"RNA","data"), GetAssayData(obj,"RNA","counts"))) {
  msg("normalizing RNA ..."); obj <- NormalizeData(obj, verbose = FALSE)
}
gc()

## ---- 1. Oligodendrocyte identity (Fig 4c genes) --------------------------
tryCatch({
  msg("== 1. oligodendrocyte identity ==")
  oli_markers <- have(c("Plp1","Mbp","Mog","Mobp","Cnp"))
  neu_markers <- have(c("Snap25","Syt1","Rbfox3","Meg3"))
  obj <- AddModuleScore(obj, list(oli_markers), name = "oliScore", assay = "RNA")
  obj <- AddModuleScore(obj, list(neu_markers), name = "neuScore", assay = "RNA")
  md <- obj@meta.data
  oli_cells <- rownames(md)[md$celltype == "oli"]
  oli_thr <- quantile(md[oli_cells,"oliScore1"], 2/3)
  neu_med <- median(md[oli_cells,"neuScore1"])
  hi_oli  <- oli_cells[md[oli_cells,"oliScore1"] >= oli_thr & md[oli_cells,"neuScore1"] < neu_med]
  msg(sprintf("high-confidence oli: %d of %d", length(hi_oli), length(oli_cells)))
  obj$oli_conf <- ifelse(colnames(obj) %in% hi_oli, "hi_conf_oli", as.character(obj$celltype))

  fig4c <- have(c("Syt1","Snap25","Rims1","Erc2","Ppfia2","Cntnap2","Adgrl3","Gabrb3","Gabrb1",
                  "Nbea","Kcnma1","Cacna1e","Cacnb2","Slc8a1","Grm5","Kalrn","Negr1","Lrrc7",
                  "Ctnnd2","Epha5","Epha6","Lrfn5","Dlgap2","Plcb1","Pde4d","Ryr2","Robo1","Slit3",
                  "Dock4","Nav3","Dcc","Dpp6","Kctd16","Hs6st3","Lingo2"))

  p <- DotPlot(obj, features = fig4c, group.by = "celltype", idents = real_ct) +
       RotatedAxis() + labs(x=NULL,y=NULL,title="Fig 4c genes across cell types (oli vs neurons)")
  save_fig(p, "oli_identity_dotplot-1.png", 14, 6)

  pct_expr <- function(cells, genes) {
    m <- GetAssayData(obj,"RNA","counts")[genes, cells, drop=FALSE]; Matrix::rowMeans(m>0)*100 }
  ex_cells <- rownames(md)[md$celltype=="ex.neu"]
  pe <- data.frame(gene=fig4c, hi_conf_oli=pct_expr(hi_oli,fig4c),
                   all_oli=pct_expr(oli_cells,fig4c), ex_neu=pct_expr(ex_cells,fig4c)) |>
        pivot_longer(-gene, names_to="group", values_to="pct")
  p <- ggplot(pe, aes(reorder(gene,pct), pct, fill=group)) +
       geom_col(position="dodge") + coord_flip() + theme_bw() +
       labs(x=NULL, y="% nuclei expressing (>0 counts)",
            title="Fig 4c genes: bona fide oligodendrocytes vs neuronal ambient")
  save_fig(p, "oli_ambient_pct-1.png", 12, 8)
  write.csv(tidyr::pivot_wider(pe, names_from=group, values_from=pct),
            file.path(FIG, "oli_fig4c_pct_expr.csv"), row.names = FALSE)

  p <- VlnPlot(obj, features = have(c("Cntnap2","Grm5","Kcnma1","Snap25","Syt1","Plp1")),
               group.by = "oli_conf", pt.size = 0, ncol = 3) & theme(axis.title.x=element_blank())
  save_fig(p, "oli_conf_violin-1.png", 12, 5)
  gc()
}, error = function(e) msg("!! section 1 failed:", conditionMessage(e)))

## ---- 2. DAM / DAA / DAO scoring by condition -----------------------------
tryCatch({
  msg("== 2. glial-state scoring ==")
  DAM <- have(c("Trem2","Tyrobp","Apoe","Cst7","Lpl","Itgax","Clec7a","Cd9","Spp1","Gpnmb","Axl",
                "Ctsb","Ctsd","Cd63","Csf1"))
  DAA <- have(c("Gfap","Serpina3n","C3","Vim","Cd44","Osmr","Ggta1","C4b","S1pr3"))
  DAO <- have(c("Hspa1a","Hspa1b","Hsp90aa1","Hsp90ab1","Dnaja1","Calm1","Calm2","Ubb",
                "Serpina3n","C4b","B2m","Klk6"))
  obj <- AddModuleScore(obj, list(DAM), name="DAM", assay="RNA")
  obj <- AddModuleScore(obj, list(DAA), name="DAA", assay="RNA")
  obj <- AddModuleScore(obj, list(DAO), name="DAO", assay="RNA")

  plot_state <- function(ct, score) {
    d <- obj@meta.data |> filter(celltype==ct)
    ggplot(d, aes(age, .data[[score]], fill=interaction(genotype,sex))) +
      geom_boxplot(outlier.shape=NA, alpha=0.8) +
      scale_fill_brewer(palette="Set2", name="genotype.sex") + theme_bw() +
      labs(title=paste0(score," in ",ct), x=NULL, y="module score") }
  p <- plot_state("mic","DAM1") | plot_state("ast","DAA1") | plot_state("oli","DAO1")
  save_fig(p, "glial_state_scores-1.png", 13, 5)

  # per-animal pseudobulk + 3-way ANOVA on each state score
  sink(file.path(FIG, "glial_state_anova.txt"))
  for (cs in list(c("mic","DAM1"), c("ast","DAA1"), c("oli","DAO1"))) {
    pb <- obj@meta.data |> filter(celltype==cs[1]) |>
      group_by(orig.ident, genotype, sex, age) |>
      summarise(score=mean(.data[[cs[2]]]), .groups="drop")
    cat("\n=====", cs[2], "in", cs[1], "=====\n")
    print(summary(aov(score ~ genotype*sex*age, data=pb)))
  }
  sink()
  saveRDS(obj@meta.data[, grep("celltype|genotype|sex|age|orig.ident|DAM1|DAA1|DAO1|oliScore1|neuScore1",
          colnames(obj@meta.data))], file.path(FIG, "glial_state_metadata.rds"))
  msg("saved glial_state_anova.txt + metadata")
  gc()
}, error = function(e) msg("!! section 2 failed:", conditionMessage(e)))

## ---- 4. Vascular compartment --------------------------------------------
tryCatch({
  msg("== 4. vascular ==")
  print(round(100*prop.table(table(obj$celltype)),2))
  p <- DotPlot(obj, features = have(c("Cldn5","Flt1","Pecam1","Pdgfrb","Rgs5","Kcnj8","Acta2","Vtn","Notch3")),
               idents = c("end","per")) + RotatedAxis() +
       labs(x=NULL,y=NULL,title="Endothelial & mural markers (verify 'per')")
  save_fig(p, "vascular_dotplot-1.png", 9, 4)
  gc()
}, error = function(e) msg("!! section 4 failed:", conditionMessage(e)))

## ---- 3. Microglia subclustering (last; RNA-based to avoid SCT) -----------
tryCatch({
  msg("== 3. microglia subclustering ==")
  mic <- subset(obj, idents = "mic")
  gc()
  mic <- NormalizeData(mic, verbose=FALSE) |>
         FindVariableFeatures(nfeatures=2000, verbose=FALSE) |>
         ScaleData(verbose=FALSE) |>
         RunPCA(npcs=20, verbose=FALSE) |>
         FindNeighbors(dims=1:15, verbose=FALSE) |>
         FindClusters(resolution=0.4, verbose=FALSE) |>
         RunUMAP(dims=1:15, verbose=FALSE)
  p1 <- DimPlot(mic, label=TRUE) + NoLegend() + ggtitle("microglia subclusters")
  p2 <- DotPlot(mic, features=have(c("P2ry12","Tmem119","Cx3cr1","Csf1r",
                "Mrc1","Lyve1","Cd163","Ms4a7","Ly6c2","Ccr2"))) + RotatedAxis() +
        labs(x=NULL,y=NULL,title="homeostatic vs PVM vs monocyte")
  save_fig(p1|p2, "microglia_subcluster-1.png", 12, 5)
  saveRDS(mic, file.path(FIG, "microglia_subcluster.rds"))
  gc()
}, error = function(e) msg("!! section 3 failed:", conditionMessage(e)))

msg("DONE in", round(difftime(Sys.time(), t0, units="mins"),1), "min")
