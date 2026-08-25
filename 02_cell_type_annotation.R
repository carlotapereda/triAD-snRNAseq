# =============================================================================
# Step 2 of 6: Cell type annotation
#
# Purpose: score curated marker gene sets (PanglaoDB / HippoSeq derived) with
#   AddModuleScore and assign each nucleus the top-scoring cell type; nuclei
#   whose top two scores are within 10% are labeled "hybrid" and nuclei with
#   no positive score are labeled "negative". Cluster identities are then
#   assigned per cluster from the most abundant per-nucleus label.
# Inputs:  clustered Seurat object "Emouse" from 01_create_obj_and_QC.R.
# Outputs: metadata columns cell_type_gen (per-nucleus label incl. hybrid/
#   negative) and Seurat_Clusters, added to the object in place.
# Note:    the marker list named "Endothelial" scores Csf1r; Csf1r is a
#   myeloid marker, so this module captures myeloid rather than endothelial
#   identity (kept as run; see the response to reviewers).
# Next:    03_DEG_analysis.R
# =============================================================================

## Cell type annotation

#cell type markers
Dentate.Gyrus<- list(c("Ntng1",  "Plk5",  "Tanc1", "Ahcyl2", "Slc4a4", "Rfx3", "Stxbp6", "Dgkh"))
CA1.Pyramidal.Neurons<- list(c("Ccdc88c", "Galntl6", "Man1a", "Man1"))
CA23.Pyramidal.Neurons<- list(c("Slit2"))
CA3.Pyramidal.Neurons<- list(c("Trhde", "Rnf182", "Nectin3"))
Interneurons<- list(c("Nxph1", "Grik1", "Sst", "Cnr1", "Vip", "Npas1"))
Subiculum<- list(c("Tshz2", "Sgcz", "Meis2"))
Oligodendrocytes<- list(c("Mbp",  "Plp1", "Mag"))
OPC<- list(c("Vcan", "Pdgfra", "Arhgap31"))
Microglia<- list(c("Arhgap45", "Ly86",  "Cd37"))
Astrocytes<-list(c("Gfap", "Slc1a2",  "Atp1a2",  "Phkg1",   "Slc7a11", "Aldh1l1",  "Ranbp3l"))
Fibroblast.like<- list(c("Slc38a2", "Mgp", "Sned1", "Slc6a13"))
Choroid.Plexus<- list(c("Col9a3","Ttr","Htr2c","Folr1","Otx2"))
Endothelial <- list(c("Csf1r"))
# Define celltype for neuronal subtype identification in hippocampal samples
Celltype <- c( "Astrocytes", "Microglia", "Fibroblast.like", "Oligodendrocytes", "OPC", "Choroid.Plexus","Dentate.Gyrus","CA1.Pyramidal.Neurons", "CA23.Pyramidal.Neurons", "CA3.Pyramidal.Neurons", "Interneurons", "Subiculum", "Endothelial")

#calculate scores
DefaultAssay(obj) <- "SCT"
dataset2 <- Emouse

library(stringr)
for(key in Celltype){
    dataset2<-AddModuleScore(object=dataset2, features = get(key), pool = NULL,  nbin = 24,  ctrl = 100,  k = FALSE,  assay = NULL, name = key,  seed = 1,  search = FALSE)
  }
new_meta <- dataset2@meta.data
colnames(new_meta)
new_meta[is.na(new_meta)] <- 0
new_meta$cell_type<- names(new_meta[, 11:23])[apply(new_meta[, 11:23],1,which.max)]
table(new_meta$cell_type)
# Find max and second max
new_meta$x1<- apply(new_meta[, 11:23], 1, max)
new_meta$x2<- apply(new_meta[, 11:23], 1, function(x) x[order(x)[12]]) #in decreasing order
# Find difference
new_meta$x1_x2<- (new_meta$x1 - new_meta$x2)/new_meta$x1
new_meta$x1_x2 <- ifelse(
new_meta$x1 == 0 & new_meta$x2 == 0,     # both zero?
    0,                                         # yes → 0
ifelse(new_meta$x1_x2 < 0.1, 1, 0)        # otherwise original rule
) #if the difference between top two scores is less than 10%, 1, if not 0
table(new_meta$x1_x2) 

new_meta$cell_type2<- new_meta$cell_type
new_meta$cell_type2<- ifelse(new_meta$x1_x2 == 1, "hybrid", new_meta$cell_type2 ) #doublets
new_meta$cell_type2<- ifelse(new_meta$x1 <= 0 , "negative", new_meta$cell_type2 ) #cell do not express markers of any included cell types
table(new_meta$cell_type2)
new_meta$cell_type2<- gsub("1", "", new_meta$cell_type2)

#add cell type annotation as metadata to seurat object
Emouse <-AddMetaData(Emouse, new_meta$cell_type2, col.name = 'cell_type_gen')
table(Emouse$cell_type_gen)

obj$Seurat_Clusters <- as.numeric(obj$seurat_clusters)

# cluster identity is assigned per cluster based on most abundant cell type in Emouse$cell_type_gen
# cells labeled as hybrid and negative were removed for downstream analysis
