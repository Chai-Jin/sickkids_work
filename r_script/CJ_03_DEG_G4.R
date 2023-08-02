library(dplyr)
library(Seurat)
library(reticulate)
library(harmony)
library(ggplot2)
library(Matrix)
library(ape)
library(future)
library(ggrepel)

my_path <- "/Users/cj/Desktop/PUGH/"

ID <- c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126")

custom_colors <- list()
colors_dutch <- c("#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67", "#F79F1F", "#A3CB38", "#1289A7", "#D980FA", "#B53471", "#EE5A24", "#009432", "#0652DD", "#9980FA", "#833471", "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51")
colors_spanish <- c("#40407a", "#706fd3", "#f7f1e3", "#34ace0", "#33d9b2", "#2c2c54", "#474787", "#aaa69d", "#227093", "#218c74", "#ff5252", "#ff793f", "#d1ccc0", "#ffb142", "#ffda79", "#b33939", "#cd6133", "#84817a", "#cc8e35", "#ccae62")
custom_colors$discrete <- c(colors_dutch, colors_spanish)

G4_whole_final <- readRDS(paste(my_path, "Integrated_G4_MBs_scRNAseq.rds", sep = ""))

# without harmony
G4_whole_final <- FindNeighbors(G4_whole_final, assay = "SCT", dims = 1:20)
G4_whole_final <- FindClusters(G4_whole_final, resolution = 0.4, algorithm = 2, assay = "SCT")
G4_whole_final <- RunUMAP(G4_whole_final, assay = "SCT", dims = 1:20)

pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_PCA_norm.pdf", sep = ""))
p <- DimPlot(G4_whole_final, group.by = "orig.ident", label = TRUE)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, label = TRUE, group.by = "Phase")
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, label = TRUE)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, label = TRUE, split.by = "orig.ident")
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_UMAP_norm.pdf", sep = ""))
p <- DimPlot(G4_whole_final, reduction = "umap", label = TRUE, group.by = "orig.ident", cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "umap", label = TRUE, group.by = "Phase", cols = c("#45aaf2", "#f1c40f", "#e74c3c"))
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "umap", label = TRUE, cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "umap", label = TRUE, split.by = "orig.ident", cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

G4_whole_final.markers <- FindAllMarkers(G4_whole_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
G4_whole_final.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(G4_whole_final, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

allmarkers <- FindAllMarkers(G4_whole_final, only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
library(writexl)
write_xlsx(allmarkers, path = paste(my_path, "All_markers.xlsx", sep = ""))

Idents(G4_whole_final) <- factor(Idents(G4_whole_final), levels = c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126"))
Idents(G4_whole_final) <- G4_whole_final$orig.ident
# up-regulated in primary
DEG1 <- FindMarkers(G4_whole_final, ident.1 = "G4_A_P_0527", ident.2 = "G4_A_R_0527", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
DEG2 <- FindMarkers(G4_whole_final, ident.1 = "G4_D_P_1303", ident.2 = "G4_D_R_1303", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
DEG3 <- FindMarkers(G4_whole_final, ident.1 = "G4_E_P_2126", ident.2 = "G4_E_R_2126", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")

# up-regulated in recurrent
DEG1 <- FindMarkers(G4_whole_final, ident.1 = "G4_A_R_0527", ident.2 = "G4_A_P_0527", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
DEG2 <- FindMarkers(G4_whole_final, ident.1 = "G4_D_R_1303", ident.2 = "G4_D_P_1303", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
DEG3 <- FindMarkers(G4_whole_final, ident.1 = "G4_E_R_2126", ident.2 = "G4_E_P_2126", only.pos = TRUE, min.pct = 0.25, test.use = "MAST")

VlnPlot(G4_whole_final, features = c("PTPRD", "KHDRBS2", "ZNF804A", "XACT", "TMSB10", "CACNA2D1", "HDAC9", "TNIK", "PTPRO", "SLC26A3", "AUTS2", "PEX5L"), pt.size = 0)
FeaturePlot(G4_whole_final, features = c("PTPRD", "KHDRBS2", "ZNF804A", "XACT", "TMSB10", "CACNA2D1", "HDAC9", "TNIK", "PTPRO", "SLC26A3", "AUTS2", "PEX5L"))