library(dplyr)
library(Seurat)
library(reticulate)
library(harmony)
library(ggplot2)
library(Matrix)
library(ape)
library(future)
library(ggrepel)
options(future.globals.maxSize = 50000 * 1024^2)
plan("multiprocess", workers = 4)

my_path <- "/Users/chaijinlee/Documents/Research/PUGH/"
# my_path <- "/hpf/largeprojects/mdtaylor/chaijin/Pugh_singlecell/"

ID <- c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126")
sample_number <- length(ID)

#read raw files (non normalized)
for (i in 1:sample_number) {
  assign(ID[i], readRDS(paste(my_path, ID[i], "_scRNAseq_raw.rds", sep = "")))
}

custom_colors <- list()
colors_dutch <- c("#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67", "#F79F1F", "#A3CB38", "#1289A7", "#D980FA", "#B53471", "#EE5A24", "#009432", "#0652DD", "#9980FA", "#833471", "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51")
colors_spanish <- c("#40407a", "#706fd3", "#f7f1e3", "#34ace0", "#33d9b2", "#2c2c54", "#474787", "#aaa69d", "#227093", "#218c74", "#ff5252", "#ff793f", "#d1ccc0", "#ffb142", "#ffda79", "#b33939", "#cd6133", "#84817a", "#cc8e35", "#ccae62")
custom_colors$discrete <- c(colors_dutch, colors_spanish)


#MERGE datasets
G4_whole <- merge(G4_A_P_0527, y = c(G4_A_R_0527, G4_D_P_1303, G4_D_R_1303, G4_E_P_2126, G4_E_R_2126), add.cell.ids = c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126"), project = "Integrated_G4_MBs_sc_PR")
saveRDS(G4_whole, file = (paste(my_path, "Integrated_G4_MBs_scRNAseq_merge_raw.rds", sep = "")))
G4_whole <- readRDS(paste(my_path, "Integrated_G4_MBs_scRNAseq_merge_raw.rds", sep = ""))


#normalization
G4_whole_final <- SCTransform(G4_whole, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = TRUE, variable.features.n = 3000, return.only.var.genes = FALSE)
saveRDS(G4_whole_final, file = paste(my_path, "Integrated_G4_MBs_scRNAseq.rds", sep = ""))
G4_whole_final <- readRDS(paste(my_path, "Integrated_G4_MBs_scRNAseq.rds", sep = ""))
G4_whole_final <- RunPCA(G4_whole_final, assay = "SCT")
pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_elbow_plot.pdf", sep = ""))
ElbowPlot(G4_whole_final)
dev.off()


G4_whole_final <- RunHarmony(G4_whole_final, group.by.vars = "orig.ident", assay.use = "SCT", max.iter.harmony = 15, theta = 2)
G4_whole_final <- RunHarmony(G4_whole_final, group.by.vars = "orig.ident", assay.use = "SCT", max.iter.harmony = 15, theta = 0)

### calculate stdv and select for highly variable components only
G4_whole_final@reductions$harmony@stdev <- apply(G4_whole_final@reductions$harmony@cell.embeddings, 2, sd)
pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_elbow_plot_harmony_theta0.pdf", sep = ""))
plot(sort(G4_whole_final@reductions$harmony@stdev, decreasing = TRUE))
dev.off()


## threshold stdv >5 but adapt it depending on the "elbow"
G4_whole_final <- FindNeighbors(G4_whole_final, reduction = "harmony", assay = "SCT", dims = which(G4_whole_final@reductions$harmony@stdev > 5))
G4_whole_final <- FindClusters(G4_whole_final, resolution = 0.4, algorithm = 2, assay = "SCT")
G4_whole_final <- RunUMAP(G4_whole_final, reduction = "harmony", assay = "SCT", dims = which(G4_whole_final@reductions$harmony@stdev > 5))


pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_PCA_norm_harmony_theta0.pdf", sep = ""))
p <- DimPlot(G4_whole_final, reduction = "harmony", group.by = "orig.ident", label = TRUE)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "harmony", label = TRUE, group.by = "Phase")
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "harmony", label = TRUE)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(G4_whole_final, reduction = "harmony", label = TRUE, split.by = "orig.ident")
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

pdf(paste(my_path, "Integrated_G4_MBs_scRNAseq_UMAP_norm_harmony_theta0.pdf", sep = ""))
p <- DimPlot(SHH_whole_final, reduction = "umap", label = TRUE, group.by = "orig.ident", cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(SHH_whole_final, reduction = "umap", label = TRUE, group.by = "Phase", cols = c("#45aaf2", "#f1c40f", "#e74c3c"))
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(SHH_whole_final, reduction = "umap", label = TRUE, cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
p <- DimPlot(SHH_whole_final, reduction = "umap", label = TRUE, split.by = "orig.ident", cols = custom_colors$discrete)
p + theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()


saveRDS(G4_whole_final, file = (paste(my_path, "Integrated_G4_MBs_scRNAseq.rds", sep = "")))
write.table(G4_whole_final@assays$SCT@scale.data, paste(my_path, "Integrated_G4_MBs_scRNAseq_scale_data.txt", sep = ""), quote = FALSE)
write.table(G4_whole_final@meta.data$seurat_clusters, paste(my_path, "Integrated_G4_MBs_scRNAseq_seuratclusters.txt", sep = ""), quote = FALSE)
saveRDS(G4_whole_final@reductions$umap@cell.embeddings, file = paste(my_path, "Integrated_G4_MBs_scRNAseq_umap_embeddings.rds", sep = ""))
mat <- as.data.frame(as.matrix(G4_whole_final@assays$SCT@data))
saveRDS(mat, file = paste(my_path, "Integrated_G4_MBs_scRNAseq_norm_data.rds", sep = ""))
