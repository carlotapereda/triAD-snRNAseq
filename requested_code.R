# =============================================================================
# Figure regeneration code (DEG tables, UpSet plots, pathway heatmaps)
#
# Purpose: regenerate the published DEG tables and figure panels:
#   (a) inhibitory-neuron all-gene DEG table (MAST, SCT-subset object),
#   (b) all-cell-type DEG tables and the significant set (p_val_adj < 0.05,
#       abs avg_log2FC >= 0.1), (c) UpSet plots of DEG overlap, (d) pathway
#       heatmaps (Supp Fig 4).
# Inputs:  annotated Seurat object "Emouse"/"obj"; pathway result CSVs under
#   data/original_pathways/ (InNeu_pathway/, Oli/).
# Outputs: Inh_all_gene_SCT_DEGs.csv, DEG.csv (all genes), sDEG.csv
#   (significant set), UpSet and heatmap panels.
# Note: this copy contains four corrections relative to the as-received file
#   (sDEG filtering, column name, cell type label, heatmap NA filter); see
#   the repository README.
# =============================================================================

## for Inh_all_gene_SCT_DEGs.csv
## setsub in.neurons
DefaultAssay(obj) <- "SCT"
Idents(obj) <- "celltype"
obj <- subset(obj, idents= "In.Neurons")
## remove doublets and negative cells
Idents(obj) <- "cell_type_gen"
obj <- subset(obj, idents= c("negative", "hybrid"), invert = TRUE)
DefaultAssay(obj) <- "RNA"
obj$sex_age <- paste(obj$sex, obj$age, sep = "_")
celltypes <- c("Female_06Mo", "Female_12Mo", "Female_18Mo", "Male_06Mo", "Male_12Mo", "Male_18Mo")
## DEG identification
library(stringr)
DEG.list.FDR2 <- NULL
DEG.list.FDR<- NULL
for (key in celltypes)  {
  Idents(obj) <- "sex_age"
  subcelltype <- subset(obj, idents= key)
  Idents(subcelltype) <- "genotype"
  DEG.list <- FindMarkers(subcelltype, ident.1= "E44", ident.2 ="E33", verbose = FALSE, test.use = "MAST", logfc.threshold = 0.01)
  DEG.list$gene <- row.names(DEG.list)
  #DEG.list.FDR <- subset(DEG.list, DEG.list$p_val_adj< 0.05)
  DEG.list.FDR <- DEG.list
  DEG.list.FDR$dir<- ifelse(DEG.list.FDR$avg_log2FC < 0, "neg","pos")
  DEG.list.FDR$celltype<- "Inhibitory.Neurons"
  DEG.list.FDR$group<- key
  DEG.list.FDR2<- rbind(DEG.list.FDR2, DEG.list.FDR)
}


write.csv(DEG.list.FDR2, "Inh_all_gene_SCT_DEGs.csv")


## DEG.csv all DEGs
## identify DEGs and stratify by age and sex

Emouse$comparison <- paste(Emouse$sex, Emouse$age, sep = "_")

comparison <- c(
  "Female_06Mo", "Female_12Mo", "Female_18Mo",
  "Male_06Mo", "Male_12Mo", "Male_18Mo"
)

Celltype <- c(
  "Excitatory.Neurons",
  "Inhibitory.Neurons",
  "Astrocytes",
  "Microglia",
  "Oligodendrocytes",
  "Oligodendrocyte.Precursor"
)

library(Seurat)
library(stringr)

obj <- Emouse
Idents(obj) <- "cell_type_gen"
obj <- subset(obj, idents= c("negative", "hybrid"), invert = TRUE)
DefaultAssay(obj) <- "SCT"

DEG.list.all <- NULL

for (key in comparison) {
  
  # subset by comparison
  Idents(obj) <- "comparison"
  sub.comp <- subset(obj, idents = key)
  
  for (ct in Celltype) {
    
    # subset by cell type
    Idents(sub.comp) <- "cell_type_ident"
    sub.ct <- subset(sub.comp, idents = ct)
    
    # DEG
    DEG.list <- FindMarkers(
      sub.ct,
      ident.1 = "E44",
      ident.2 = "E33",
      verbose = TRUE,
      test.use = "MAST",
      logfc.threshold = 0.01
    )
    
    if (nrow(DEG.list) == 0) {
      message("No DEG found for ", key, " / ", ct)
      next
    }
    
    DEG.list$gene <- rownames(DEG.list)
    
    # Seurat version compatibility
    if ("avg_log2FC" %in% colnames(DEG.list)) {
      DEG.list$dir <- ifelse(DEG.list$avg_log2FC < 0, "neg", "pos")
    } else if ("avg_logFC" %in% colnames(DEG.list)) {
      DEG.list$dir <- ifelse(DEG.list$avg_logFC < 0, "neg", "pos")
    } else {
      DEG.list$dir <- NA
    }
    
    DEG.list$comparison <- key
    DEG.list$celltype <- ct
    
    DEG.list.all <- rbind(DEG.list.all, DEG.list)
  }
}

write.csv(DEG.list.all, "DEG.csv")

## sDEG.csv

DEG.list <- subset(DEG.list.all, DEG.list.all$p_val_adj< 0.05)
DEG.list <- subset(DEG.list, abs(DEG.list$avg_log2FC) >=0.1)
write.csv(DEG.list, "sDEG.csv")

## Supp Figure 2a
#get cell numbers from seurat object
Emouse$group <- paste(Emouse$sex, Emouse$age, Emouse$genotype, sep = "_")
Emouse$group_celltype <- paste(Emouse$group, Emouse$cell_type_ident, sep = "_")
Emouse$individual_group_celltype <- paste(Emouse$orig.ident, Emouse$group_celltype, sep = "_")
#get cell number per animal
orig.ident <- as.data.frame(table(Emouse$orig.ident))
colnames(orig.ident) <- c("name", "total.cell")
group_celltype <- table(Emouse$individual_group_celltype)
group_celltype <- as.data.frame(group_celltype)
group_celltype <- group_celltype %>%
  separate(Var1, c("group","individual", "sex", "age", "genotype", "cell_type"), "_")
group_celltype$name <- paste(group_celltype$group, group_celltype$individual, sep = "_")
group_celltype <- inner_join(group_celltype, orig.ident, by ="name")
group_celltype$celltype.ratio <- (group_celltype$Freq / group_celltype$total.cell) * 100

# ---- plotting ----
# set display order 
group_celltype$age      <- factor(group_celltype$age,      levels = c("6mos", "12mos", "18mos"))
group_celltype$sex      <- factor(group_celltype$sex,      levels = c("Female", "Male"))
group_celltype$genotype <- factor(group_celltype$genotype, levels = c("E33", "E44"))
group_celltype$cell_type <- factor(group_celltype$cell_type,
                                   levels = c("Astrocytes", "Excitatory.Neurons", "Inhibitory.Neurons",
                                              "Microglia", "Oligodendrocyte.Precursor", "Oligodendrocytes"))

# order bars: genotype -> sex -> age -> sample
ord <- group_celltype %>%
  distinct(name, genotype, sex, age) %>%
  arrange(genotype, sex, age, name)
group_celltype$name <- factor(group_celltype$name, levels = ord$name)

# black lines between the 4 blocks
blk <- ord %>% mutate(b = paste(genotype, sex))
vlines <- which(diff(as.integer(factor(blk$b, levels = unique(blk$b)))) != 0) + 0.5

cols <- c("Astrocytes"                = "#B3DE69",
          "Excitatory.Neurons"        = "#BEBADA",
          "Inhibitory.Neurons"        = "#FDB462",
          "Microglia"                 = "#FFFFB3",
          "Oligodendrocyte.Precursor" = "#80B1D3",
          "Oligodendrocytes"          = "#FB8072")

ggplot(group_celltype, aes(x = name, y = Freq, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  geom_vline(xintercept = vlines, linewidth = 0.4) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Freq") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.position = "left")




## Figure 5c 

#pathway overlaps - Inhibitory neurons

#old.male.up <- read.csv("data/original_pathways/InNeu_pathway/")
old.female.up <- read.csv("data/original_pathways/InNeu_pathway/F18_up.csv")
old.male.down <- read.csv("data/original_pathways/InNeu_pathway/M18_down.csv")
old.female.down <- read.csv("data/original_pathways/InNeu_pathway/F18_down.csv")
#old.male.up$male_18mos <- 1
old.female.up$female_18mos <- 1
old.male.down$male_18mos <- -1
old.female.down$female_18mos <- -1


adult.male.up <- read.csv("data/original_pathways/InNeu_pathway/M12_up.csv")
#adult.female.up <- read.csv("data/original_pathways/InNeu_pathway/F12_up.csv")
adult.male.down <- read.csv("data/original_pathways/InNeu_pathway/M12_down.csv")
adult.female.down <- read.csv("data/original_pathways/InNeu_pathway/F12_down.csv")
adult.male.up$male_12mos <- 1
#adult.female.up$female_12mos <- 1
adult.male.down$male_12mos <- -1
adult.female.down$female_12mos <- -1
adult <- join_all(list( adult.male.up, adult.male.down, adult.female.down), by = 'term_id', type = 'full')


young.male.up <- read.csv("data/original_pathways/InNeu_pathway/M06_up.csv")
young.female.up <- read.csv("data/original_pathways/InNeu_pathway/F06_up.csv")
young.male.down <- read.csv("data/original_pathways/InNeu_pathway/M06_down.csv")
young.female.down <- read.csv("data/original_pathways/InNeu_pathway/F06_down.csv")

young.male.up$male_06mos <- 1
young.female.up$female_06mos <- 1
young.male.down$male_06mos <- -1
young.female.down$female_06mos <- -1


males <- join_all(list( young.male.up, young.male.down, adult.male.down, old.male.down, adult.male.up), by = 'term_id', type = 'full')
females <- join_all(list(young.female.up, young.female.down, adult.female.down, old.female.up, old.female.down), by = 'term_id', type = 'full')

combined1 <- join_all(list(males, females), by = 'term_id', type = 'full')
combined1 <- combined1[, c(1, 2, 12:17)]
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]}

combined_overlaps <- delete.na(combined1, 4) #116

combined_overlaps$term_name <- paste(combined_overlaps$source, combined_overlaps$term_name, sep = " : ")


row.names(combined_overlaps) <- combined_overlaps$term_name
combined_overlaps$term_name <-NULL
combined_overlaps$source <-NULL
combined_overlaps[is.na(combined_overlaps)] <- 0
row.names(combined_overlaps) <- ifelse(row.names(combined_overlaps) == "GO:BP : adenylate cyclase-modulating G protein-coupled receptor signaling pathway", "GO:BP : adenylate cyclase-modulating G protein-coupled receptor", row.names(combined_overlaps))

my.breaks <- c(seq(-3, -0.1, by=0.1), -0.01, 0.01, seq(0.1, 3, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white",  "red"))(length(my.breaks)/2))
pdf("InN_overlap_pathways3.pdf", width = 9, height = 15)
pheatmap(combined_overlaps, color = my.colors, breaks = my.breaks, cluster_rows = TRUE, cluster_cols = FALSE,
         labels_col = c( "Male- 06 mos", "Male- 12mos", "Male- 18mos", "Female- 06mos",  "Female- 12mos", "Female- 18mos"), angle_col = "315", angle_row = "0",
         fontsize_row = 14, fontsize_col = 14, border_color = TRUE, treeheight_row = 0, legend = FALSE)
dev.off()



# code for oli pathway Supp Fig 4


#pathway overlaps -oli


old.male.up <- read.csv("data/original_pathways/Oli/M18_up.csv")
#old.female.up <- read.csv("data/original_pathways/Oli/F18_up.csv")
old.male.down <- read.csv("data/original_pathways/Oli/M18_down.csv")
old.female.down <- read.csv("data/original_pathways/Oli/F18_down.csv")
old.male.up$male_18mos <- 1
#old.female.up$female_18mos <- 1
old.male.down$male_18mos <- -1
old.female.down$female_18mos <- -1

adult.male.up <- read.csv("data/original_pathways/Oli/M12_up.csv")
adult.female.up <- read.csv("data/original_pathways/Oli/F12_up.csv")
adult.male.down <- read.csv("data/original_pathways/Oli/M12_down.csv")
adult.female.down <- read.csv("data/original_pathways/Oli/F12_down.csv")
adult.male.up$male_12mos <- 1
adult.female.up$female_12mos <- 1
adult.male.down$male_12mos <- -1
adult.female.down$female_12mos <- -1

young.male.up <- read.csv("data/original_pathways/Oli/M06_up.csv")
#young.female.up <- read.csv("data/original_pathways/Oli/F06_up.csv")
young.male.down <- read.csv("data/original_pathways/Oli/M06_down.csv")
young.female.down <- read.csv("data/original_pathways/Oli/F06_down.csv")
young.male.up$male_06mos <- 1
#young.female.up$female_06mos <- 1
young.male.down$male_06mos <- -1
young.female.down$female_06mos <- -1


males <- join_all(list( young.male.up, young.male.down, adult.male.down,adult.male.up, old.male.down, old.male.up), by = 'term_id', type = 'full')
females <- join_all(list( young.female.down,adult.female.up, adult.female.down, old.female.down), by = 'term_id', type = 'full')

combined1 <- join_all(list(males, females), by = 'term_id', type = 'full')
combined1 <- combined1[, c(1, 2, 12:17)]
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]}

combined_overlaps <- delete.na(combined1, 4) #23 terms

combined_overlaps$term_name <- paste(combined_overlaps$source, combined_overlaps$term_name, sep = " : ")


row.names(combined_overlaps) <- combined_overlaps$term_name
combined_overlaps$term_name <-NULL
combined_overlaps$source <-NULL
combined_overlaps[is.na(combined_overlaps)] <- 0
row.names(combined_overlaps) <- ifelse(row.names(combined_overlaps) == "GO:BP : adenylate cyclase-modulating G protein-coupled receptor signaling pathway", "GO:BP : adenylate cyclase-modulating G protein-coupled receptor", row.names(combined_overlaps))

my.breaks <- c(seq(-3, -0.1, by=0.1), -0.01, 0.01, seq(0.1, 3, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white",  "red"))(length(my.breaks)/2))
pdf("Oli_overlap_pathways.pdf", width = 9, height = 8)
pheatmap(combined_overlaps, color = my.colors, breaks = my.breaks, cluster_rows = TRUE, cluster_cols = FALSE,
         labels_col = c( "Male- 06 mos", "Male- 12mos", "Male- 18mos", "Female- 06mos",  "Female- 12mos", "Female- 18mos"), angle_col = "315", angle_row = "0",
         fontsize_row = 14, fontsize_col = 14, border_color = TRUE, treeheight_row = 0, legend = FALSE)
dev.off()




## doublets
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










#upset plot for each cell type
```{r}

library(UpSetR, lib.loc = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
DEGall <- read.csv( file = "sDEG.csv" )
Celltype <- c( "Excitatory.Neurons", "Inhibitory.Neurons", "Astrocytes", "Microglia", "Oligodendrocytes", "Oligodendrocyte.Precursor")

CCA06 <- subset(DEGall, DEGall$age == "06Mo" )
CCA12 <- subset(DEGall, DEGall$age == "12Mo" )
CCA18 <- subset(DEGall, DEGall$age == "18Mo")

DEGall$gene[DEGall$group == "Female_12Mo_Astrocytes" & DEGall$dir == "pos"]
DEGall$gene[DEGall$group == "Male_12Mo_Astrocytes" & DEGall$dir == "neg"]

Celltype <- c("Astrocytes", "Microglia", "Excitatory.Neurons", "Inhibitory.Neurons", "Oligodendrocytes", "Oligodendrocyte.Precursor")
Celltype <- "Oligodendrocyte.Precursor"

de.list1<-(list(M.down.18mo= CCA18$gene[CCA18$celltype == Celltype & CCA18$dir == "neg" & CCA18$sex == "Male" ],
                M.down.12mo= CCA12$gene[CCA12$celltype == Celltype & CCA12$dir == "neg" & CCA12$sex == "Male"],
                M.down.06mo= CCA06$gene[CCA06$celltype == Celltype & CCA06$dir == "neg" & CCA06$sex == "Male"],
                M.up.18mo= CCA18$gene[CCA18$celltype == Celltype & CCA18$dir == "pos" & CCA18$sex == "Male"],
                M.up.12mo= CCA12$gene[CCA12$celltype == Celltype & CCA12$dir == "pos" & CCA12$sex == "Male"],
                M.up.06mo= CCA06$gene[CCA06$celltype == Celltype & CCA06$dir == "pos" & CCA06$sex == "Male"],
                F.down.18mo= CCA18$gene[CCA18$celltype == Celltype & CCA18$dir == "neg" & CCA18$sex == "Female"],
                F.down.12mo= CCA12$gene[CCA12$celltype == Celltype & CCA12$dir == "neg" & CCA12$sex == "Female"],
                F.down.06mo= CCA06$gene[CCA06$celltype == Celltype & CCA06$dir == "neg" & CCA06$sex == "Female"],
                F.up.18mo= CCA18$gene[CCA18$celltype == Celltype & CCA18$dir == "pos" & CCA18$sex == "Female"],
                F.up.12mo= CCA12$gene[CCA12$celltype == Celltype & CCA12$dir == "pos" & CCA12$sex == "Female"],
                F.up.06mo= CCA06$gene[CCA06$celltype == Celltype & CCA06$dir == "pos" & CCA06$sex == "Female"]))
x <- upset(fromList(de.list1), order.by = "freq", sets = names(de.list1), keep.order = T, set_size.show = TRUE, set_size.scale_max = 350, mainbar.y.label = "Oli DEG Counts", nsets=12, sets.x.label = "Set Size", text.scale = c(1.75, 1.75, 1.75, 1.75, 1.75, 1.5), mb.ratio = c(0.50, 0.50), nintersects = 80 ) #changed nsets to 2
x <- upset(fromList(de.list1), order.by = "freq", sets = names(de.list1), keep.order = T, set_size.show = TRUE, set_size.scale_max = 350, mainbar.y.label = "In.Neu DEG Counts", nsets=12, sets.x.label = "Set Size", text.scale = c(1.25, 1.25, 1.25, 1.25, 1.25, 1.5), mb.ratio = c(0.50, 0.50), nintersects = 80 ) #changed nsets to 2


# for Figure 3d
CCA <- subset(CCA18, CCA18$sex == "Male")

CCA <- subset(CCA18, CCA18$sex == "Female")

de.list1<-(list(down.Opc= CCA$gene[CCA$celltype == "Oligodendrocyte.Precursor" & CCA$dir == "neg"],
                down.Oli= CCA$gene[CCA$celltype == "Oligodendrocytes" & CCA$dir == "neg"],
                down.In.Neu= CCA$gene[CCA$celltype == "Inhibitory.Neurons" & CCA$dir == "neg"],
                down.Ex.Neu= CCA$gene[CCA$celltype == "Excitatory.Neurons" & CCA$dir == "neg"],
                down.Ast= CCA$gene[CCA$celltype == "Astrocytes" & CCA$dir == "neg"],
                down.Mic= CCA$gene[CCA$celltype == "Microglia" & CCA$dir == "neg"],
                up.Opc= CCA$gene[CCA$celltype == "Oligodendrocyte.Precursor" & CCA$dir == "pos"],
                up.Oli= CCA$gene[CCA$celltype == "Oligodendrocytes" & CCA$dir == "pos"],
                up.In.Neu= CCA$gene[CCA$celltype == "Inhibitory.Neurons" & CCA$dir == "pos"], 
                up.Ex.Neu= CCA$gene[CCA$celltype == "Excitatory.Neurons" & CCA$dir == "pos"],
                up.Ast= CCA$gene[CCA$celltype == "Astrocytes" & CCA$dir == "pos"],
                up.Mic= CCA$gene[CCA$celltype == "Microglia" & CCA$dir == "pos"]))
plot1 <- upset(fromList(de.list1), order.by = "freq", sets = names(de.list1), keep.order = T, set_size.show = TRUE,set_size.scale_max = 180, mainbar.y.label = "Gene Intersections", nsets=12, sets.x.label = "Set Size", text.scale =  c(1.75, 1.75, 1.75, 1.75, 1.75, 1.5), mb.ratio = c(0.50, 0.50), nintersects = 30 )
