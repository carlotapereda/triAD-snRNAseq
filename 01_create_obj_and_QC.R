# =============================================================================
# Step 1 of 6: Create Seurat object and quality control
#
# Purpose: read the 47 per-library Cell Ranger count matrices, merge them into
#   one Seurat object with sample metadata (age, sex, genotype), filter nuclei
#   and genes on QC metrics, normalize with SCTransform (v2), then run PCA,
#   UMAP, and Louvain clustering.
# Inputs:  47 Cell Ranger filtered count matrix directories, one per library
#   (Read10X paths under ~/data/ApoEmouse/ApoEcount_only/).
# Outputs: merged, QC-filtered, clustered Seurat object "Emouse"
#   (555,008 nuclei x 27,153 genes; saved as Emouse.rds); QC violin plots.
# Next:    02_cell_type_annotation.R
# =============================================================================

#read counts and creat seruat object

library(SeuratData)
library(Seurat)
library(patchwork)
library(dittoSeq)
library(ggplot2)
library(Seurat) 
library(dplyr)
library(sctransform)
library(cowplot)
library(BRETIGEA)
library(knitr)
library(stringr)
library(patchwork)

## E3 data process

E3F12Mo_275 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F12Mo_275/")
E3F12Mo_275 <- CreateSeuratObject(counts = E3F12Mo_275, project = "E3F12Mo_275")
E3F12Mo_288 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F12Mo_288/")
E3F12Mo_288 <- CreateSeuratObject(counts = E3F12Mo_288, project = "E3F12Mo_288")
E3F12Mo_304 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F12Mo_304/")
E3F12Mo_304 <- CreateSeuratObject(counts = E3F12Mo_304, project = "E3F12Mo_304")
E3F12Mo_305 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F12Mo_305/")
E3F12Mo_305 <- CreateSeuratObject(counts = E3F12Mo_305, project = "E3F12Mo_305")
E3F12Mo <- merge(E3F12Mo_275, y = c(E3F12Mo_288, E3F12Mo_304, E3F12Mo_305), add.cell.ids = c("E3F12Mo_275", "E3F12Mo_288", "E3F12Mo_304", "E3F12Mo_305"), project = "E3F12Mo")
E3F12Mo@meta.data$age <- "12 Mo"

E3F18Mo_123 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F18Mo_123/")
E3F18Mo_123 <- CreateSeuratObject(counts = E3F18Mo_123, project = "E3F18Mo_123")
E3F18Mo_162 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F18Mo_162/")
E3F18Mo_162 <- CreateSeuratObject(counts = E3F18Mo_162, project = "E3F18Mo_162")
E3F18Mo_208 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F18Mo_208/")
E3F18Mo_208 <- CreateSeuratObject(counts = E3F18Mo_208, project = "E3F18Mo_208")
E3F18Mo_180 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F18Mo_180/")
E3F18Mo_180 <- CreateSeuratObject(counts = E3F18Mo_180, project = "E3F18Mo_180")
E3F18Mo <- merge(E3F18Mo_123, y = c(E3F18Mo_162, E3F18Mo_208, E3F18Mo_180), add.cell.ids = c("E3F18Mo_123", "E3F18Mo_162", "E3F18Mo_208", "E3F18Mo_180"), project = "E3F18Mo")
E3F18Mo@meta.data$age <- "18 Mo"

E3F6Mo_348 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F6Mo_348/")
E3F6Mo_348 <- CreateSeuratObject(counts = E3F6Mo_348, project = "E3F6Mo_348")
E3F6Mo_349 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F6Mo_349/")
E3F6Mo_349 <- CreateSeuratObject(counts = E3F6Mo_349, project = "E3F6Mo_349")
E3F6Mo_360 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F6Mo_360/")
E3F6Mo_360 <- CreateSeuratObject(counts = E3F6Mo_360, project = "E3F6Mo_360")
E3F6Mo_375 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3F6Mo_375/")
E3F6Mo_375 <- CreateSeuratObject(counts = E3F6Mo_375, project = "E3F6Mo_375")
E3F6Mo <- merge(E3F6Mo_348, y = c(E3F6Mo_349, E3F6Mo_360, E3F6Mo_375), add.cell.ids = c("E3F6Mo_348", "E3F6Mo_349", "E3F6Mo_360", "E3F6Mo_375"), project = "E3F6Mo")
E3F6Mo@meta.data$age <- "6 Mo"

E3F <- merge(E3F6Mo, y = c(E3F12Mo, E3F18Mo), project = "E3F")
E3F@meta.data$sex <- "Female"

rm(list=setdiff(ls(), "E3F"))

E3M12Mo_280 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M12Mo_280/")
E3M12Mo_280 <- CreateSeuratObject(counts = E3M12Mo_280, project = "E3M12Mo_280")
E3M12Mo_285 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M12Mo_285/")
E3M12Mo_285 <- CreateSeuratObject(counts = E3M12Mo_285, project = "E3M12Mo_285")
E3M12Mo_286 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M12Mo_286/")
E3M12Mo_286 <- CreateSeuratObject(counts = E3M12Mo_286, project = "E3M12Mo_286")
E3M12Mo_306 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M12Mo_306/")
E3M12Mo_306 <- CreateSeuratObject(counts = E3M12Mo_306, project = "E3M12Mo_306")
E3M12Mo <- merge(E3M12Mo_280, y = c(E3M12Mo_285, E3M12Mo_286, E3M12Mo_306), add.cell.ids = c("E3M12Mo_280", "E3M12Mo_285", "E3M12Mo_286", "E3M12Mo_306"), project = "E3M12Mo")
E3M12Mo@meta.data$age <- "12 Mo"

E3M18Mo_147 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M18Mo_147/")
E3M18Mo_147 <- CreateSeuratObject(counts = E3M18Mo_147, project = "E3M18Mo_147")
E3M18Mo_149 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M18Mo_149/")
E3M18Mo_149 <- CreateSeuratObject(counts = E3M18Mo_149, project = "E3M18Mo_149")
E3M18Mo_169 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M18Mo_169/")
E3M18Mo_169 <- CreateSeuratObject(counts = E3M18Mo_169, project = "E3M18Mo_169")
E3M18Mo_185 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M18Mo_185/")
E3M18Mo_185 <- CreateSeuratObject(counts = E3M18Mo_185, project = "E3M18Mo_185")
E3M18Mo <- merge(E3M18Mo_147, y = c(E3M18Mo_149, E3M18Mo_169, E3M18Mo_185), add.cell.ids = c("E3M18Mo_147", "E3M18Mo_149", "E3M18Mo_169", "E3M18Mo_185"), project = "E3M18Mo")
E3M18Mo@meta.data$age <- "18 Mo"

E3M6Mo_4545 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M6Mo_4545/")
E3M6Mo_4545 <- CreateSeuratObject(counts = E3M6Mo_4545, project = "E3M6Mo_4545")
E3M6Mo_4546 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M6Mo_4546/")
E3M6Mo_4546 <- CreateSeuratObject(counts = E3M6Mo_4546, project = "E3M6Mo_4546")
E3M6Mo_4547 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E3M6Mo_4547/")
E3M6Mo_4547 <- CreateSeuratObject(counts = E3M6Mo_4547, project = "E3M6Mo_4547")

E3M6Mo <- merge(E3M6Mo_4545, y = c(E3M6Mo_4546, E3M6Mo_4547), add.cell.ids = c("E3M6Mo_4545", "E3M6Mo_4546", "E3M6Mo_4547"), project = "E3M6Mo")
E3M6Mo@meta.data$age <- "6 Mo"

E3M <- merge(E3M6Mo, y = c(E3M12Mo, E3M18Mo), project = "E3M")
E3M@meta.data$sex <- "Male"


E3 <- merge(E3M, y = c(E3F), project = "E3")
E3@meta.data$genotype <- "E3/3"

rm(list=setdiff(ls(), "E3"))





## E4 data process
E4F12Mo_271 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F12Mo_271/")
E4F12Mo_271 <- CreateSeuratObject(counts = E4F12Mo_271, project = "E4F12Mo_271")
E4F12Mo_272 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F12Mo_272/")
E4F12Mo_272 <- CreateSeuratObject(counts = E4F12Mo_272, project = "E4F12Mo_272")
E4F12Mo_273 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F12Mo_273/")
E4F12Mo_273 <- CreateSeuratObject(counts = E4F12Mo_273, project = "E4F12Mo_273")
E4F12Mo_274 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F12Mo_274/")
E4F12Mo_274 <- CreateSeuratObject(counts = E4F12Mo_274, project = "E4F12Mo_274")
E4F12Mo <- merge(E4F12Mo_271, y = c(E4F12Mo_272, E4F12Mo_273, E4F12Mo_274), add.cell.ids = c("E4F12Mo_271", "E4F12Mo_272", "E4F12Mo_273", "E4F12Mo_274"), project = "E4F12Mo")
E4F12Mo@meta.data$age <- "12 Mo"

E4F18Mo_137 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F18Mo_137/")
E4F18Mo_137 <- CreateSeuratObject(counts = E4F18Mo_137, project = "E4F18Mo_137")
E4F18Mo_138 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F18Mo_138/")
E4F18Mo_138 <- CreateSeuratObject(counts = E4F18Mo_138, project = "E4F18Mo_138")
E4F18Mo_159 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F18Mo_159/")
E4F18Mo_159 <- CreateSeuratObject(counts = E4F18Mo_159, project = "E4F18Mo_159")
E4F18Mo_210 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F18Mo_210/")
E4F18Mo_210 <- CreateSeuratObject(counts = E4F18Mo_210, project = "E4F18Mo_210")
E4F18Mo <- merge(E4F18Mo_137, y = c(E4F18Mo_138, E4F18Mo_159, E4F18Mo_210), add.cell.ids = c("E4F18Mo_137", "E4F18Mo_138", "E4F18Mo_159", "E4F18Mo_210"), project = "E4F18Mo")
E4F18Mo@meta.data$age <- "18 Mo"

E4F6Mo_31 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F6Mo_31/")
E4F6Mo_31 <- CreateSeuratObject(counts = E4F6Mo_31, project = "E4F6Mo_31")
E4F6Mo_32 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F6Mo_32/")
E4F6Mo_32 <- CreateSeuratObject(counts = E4F6Mo_32, project = "E4F6Mo_32")
E4F6Mo_378 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F6Mo_378/")
E4F6Mo_378 <- CreateSeuratObject(counts = E4F6Mo_378, project = "E4F6Mo_378")
E4F6Mo_NA <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4F6Mo_NA/")
E4F6Mo_NA <- CreateSeuratObject(counts = E4F6Mo_NA, project = "E4F6Mo_NA")

E4F6Mo <- merge(E4F6Mo_31, y = c(E4F6Mo_32, E4F6Mo_378, E4F6Mo_NA), add.cell.ids = c("E4F6Mo_31", "E4F6Mo_32", "E4F6Mo_378", "E4F6Mo_NA"), project = "E4F6Mo")
E4F6Mo@meta.data$age <- "6 Mo"

E4F <- merge(E4F6Mo, y = c(E4F12Mo, E4F18Mo), project = "E4F")
E4F@meta.data$sex <- "Female"

rm(list=ls()[! ls() %in% c("E4F","E3")])


E4M12Mo_276 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M12Mo_276/")
E4M12Mo_276 <- CreateSeuratObject(counts = E4M12Mo_276, project = "E4M12Mo_276")
E4M12Mo_277 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M12Mo_277/")
E4M12Mo_277 <- CreateSeuratObject(counts = E4M12Mo_277, project = "E4M12Mo_277")
E4M12Mo_293 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M12Mo_293/")
E4M12Mo_293 <- CreateSeuratObject(counts = E4M12Mo_293, project = "E4M12Mo_293")
E4M12Mo_295 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M12Mo_295/")
E4M12Mo_295 <- CreateSeuratObject(counts = E4M12Mo_295, project = "E4M12Mo_295")
E4M12Mo <- merge(E4M12Mo_276, y = c(E4M12Mo_277, E4M12Mo_293, E4M12Mo_295), add.cell.ids = c("E4M12Mo_276", "E4M12Mo_277", "E4M12Mo_293", "E4M12Mo_295"), project = "E4M12Mo")
E4M12Mo@meta.data$age <- "12 Mo"

E4M18Mo_151 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M18Mo_151/")
E4M18Mo_151 <- CreateSeuratObject(counts = E4M18Mo_151, project = "E4M18Mo_151")
E4M18Mo_152 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M18Mo_152/")
E4M18Mo_152 <- CreateSeuratObject(counts = E4M18Mo_152, project = "E4M18Mo_152")
E4M18Mo_168 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M18Mo_168/")
E4M18Mo_168 <- CreateSeuratObject(counts = E4M18Mo_168, project = "E4M18Mo_168")
E4M18Mo_238 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M18Mo_238/")
E4M18Mo_238 <- CreateSeuratObject(counts = E4M18Mo_238, project = "E4M18Mo_238")
E4M18Mo <- merge(E4M18Mo_151, y = c(E4M18Mo_152, E4M18Mo_168, E4M18Mo_238), add.cell.ids = c("E4M18Mo_151", "E4M18Mo_152", "E4M18Mo_168", "E4M18Mo_238"), project = "E4M18Mo")
E4M18Mo@meta.data$age <- "18 Mo"

E4M6Mo_34 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M6Mo_34/")
E4M6Mo_34 <- CreateSeuratObject(counts = E4M6Mo_34, project = "E4M6Mo_34")
E4M6Mo_35 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M6Mo_35/")
E4M6Mo_35 <- CreateSeuratObject(counts = E4M6Mo_35, project = "E4M6Mo_35")
E4M6Mo_36 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M6Mo_36/")
E4M6Mo_36 <- CreateSeuratObject(counts = E4M6Mo_36, project = "E4M6Mo_36")
E4M6Mo_37 <- Read10X(data.dir = "~/data/ApoEmouse/ApoEcount_only/E4M6Mo_37/")
E4M6Mo_37 <- CreateSeuratObject(counts = E4M6Mo_37, project = "E4M6Mo_37")

E4M6Mo <- merge(E4M6Mo_34, y = c(E4M6Mo_35, E4M6Mo_36, E4M6Mo_37), add.cell.ids = c("E4M6Mo_34", "E4M6Mo_35", "E4M6Mo_36", "E4M6Mo_37"), project = "E4M6Mo")
E4M6Mo@meta.data$age <- "6 Mo"

E4M <- merge(E4M6Mo, y = c(E4M12Mo, E4M18Mo), project = "E4M")
E4M@meta.data$sex <- "Male"

E4 <- merge(E4M, y = c(E4F), project = "E4")
E4@meta.data$genotype <- "E4/4"

rm(list=ls()[! ls() %in% c("E4","E3")])

Emouse <- merge(E3, y = c(E4), project = "ApoEmouse")





#QC: remove poor quality cells and low expressing features

Emouse <- PercentageFeatureSet(Emouse, pattern = "^mt-", col.name = "percent.mt")

# Add number of genes per UMI for each cell to metadata
Emouse$log10GenesPerUMI <- log10(Emouse$nFeature_RNA) / log10(Emouse$nCount_RNA)

VlnPlot(Emouse, features = c("log10GenesPerUMI", "percent.mt"), ncol = 3)
VlnPlot(Emouse, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
Emouse <- subset(Emouse, subset = nFeature_RNA > 250 & log10GenesPerUMI > 0.85 & nCount_RNA > 500)

#Gene-level filtering
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = Emouse, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
Emouse <- CreateSeuratObject(filtered_counts, meta.data = Emouse@meta.data)

rm(counts, filtered_counts, nonzero, individualID)
#rm(list=setdiff(ls(), "Emouse"))


## Normalize and find clusters
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)

Emouse <- SCTransform(Emouse, vst.flavor = "v2")

Emouse <- RunPCA(Emouse)
Emouse <- RunUMAP(Emouse, reduction = "pca", dims = 1:15)

Emouse <- FindNeighbors(Emouse, reduction = "pca", dims = 1:15) %>%
    FindClusters(resolution = 0.5)
save.image("~/data/ApoEmouse/sctransformed.RData")
saveRDS(Emouse, "Emouse.rds")

