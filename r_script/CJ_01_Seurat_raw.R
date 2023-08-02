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

my_path <- "/Users/chaijinlee/Documents/Research/"
#server_path <- "/hpf/largeprojects/mdtaylor/chaijin/Pugh_singlecell/"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


repeated_path <- "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/"
# Load the G4 dataset
A_P_0527 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/A_P_MDT-AP-0527-0527_200218/filtered_feature_bc_matrix/")
A_R_0527 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/A_R-MET_MDT-AP-0527-0529_200218/filtered_feature_bc_matrix")
D_P_1303 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/D_P_MDT-AP-1303-1303_200206/filtered_feature_bc_matrix")
D_R_1303 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/D_R_MDT-AP-1303-0560_200206/filtered_feature_bc_matrix")
E_P_2126 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/E_P_MDT-AP-2126-2126_190913/filtered_feature_bc_matrix")
E_R_2126 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/E_R-MET_MDT-AP-2126-2127_190913/filtered_feature_bc_matrix")

# Load the SHH dataset
G_P_2622 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/G_P_MDT-AP-2622-2622_200115/filtered_feature_bc_matrix")
G_R_2622 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/G_R_MDT-AP-2622-2956_200115/filtered_feature_bc_matrix")
L_P_3661 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/L_P_MDT-AP-3661-3661_190913/filtered_feature_bc_matrix")
L_R_3661 <- Read10X(data.dir = "/Users/chaijinlee/Documents/Research/PUGH/CellRanger/L_R_MDT-AP-3661-3662_191017/filtered_feature_bc_matrix")


G4.list <- list(A_P_0527, A_R_0527, D_P_1303, D_R_1303, E_P_2126, E_R_2126)
ID <- c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126")

SHH.list <- list(G_P_2622, G_R_2622, L_P_3661, L_R_3661)
ID <- c("SHH_G_P_2622", "SHH_G_R_2622", "SHH_L_P_3661", "SHH_L_R_3661")

#Filter out outliers cells and compute cell cycle scoring
for (i in 1:length(G4.list)) {
	G4.list[[i]] <- CreateSeuratObject(counts = G4.list[[i]], project = ID[i], min.cells = 10, min.features = 500)
    G4.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(G4.list[[i]], pattern = "^MT-")
	#VlnPlot(G4.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	median_nCount <- median(G4.list[[i]]$nCount_RNA)
	median_nFeature <- median(G4.list[[i]]$nFeature_RNA)
	median_percent_MT <- median(G4.list[[i]]$percent.mt)
	mad_nCount <- mad(G4.list[[i]]$nCount_RNA)
	mad_nFeature <- mad(G4.list[[i]]$nFeature_RNA)
	mad_percent_MT <- mad(G4.list[[i]]$percent.mt)
	thresholds_nCount <- median_nCount + 5*mad_nCount
	thresholds_nFeature <- median_nFeature + 5*mad_nFeature
	thresholds_percent_MT <- median_percent_MT + 5*mad_percent_MT
	G4.list[[i]] <- subset(G4.list[[i]], subset = nFeature_RNA < thresholds_nFeature & nCount_RNA < thresholds_nCount & percent.mt < thresholds_percent_MT)
	#perform scoring of cell cycle genes in order to do linear regression during normalization steps
	G4.list[[i]] <- CellCycleScoring(G4.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
	saveRDS(G4.list[[i]], file = paste(my_path,ID[i],"_scRNAseq_raw.rds",sep=""))
}

for (i in 1:length(SHH.list)) {
	SHH.list[[i]] <- CreateSeuratObject(counts = SHH.list[[i]], project = ID[i], min.cells = 10, min.features = 500)
    SHH.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(SHH.list[[i]], pattern = "^MT-")
	#VlnPlot(SHH.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

	#perform scoring of cell cycle genes in order to do linear regression during normalization steps
	SHH.list[[i]] <- CellCycleScoring(SHH.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
	saveRDS(SHH.list[[i]], file = paste(my_path,ID[i],"_scRNAseq_raw.rds",sep=""))
}

#read raw files (non normalized)
for (i in 1:length(ID)) {
	assign(ID[i],readRDS(paste(my_path,ID[i],"_scRNAseq_raw.rds",sep="")))
}

SCT.list <- list(G4_A_P_0527,G4_A_R_0527,G4_D_P_1303,G4_D_R_1303,G4_E_P_2126,G4_E_R_2126)
SCT.list <- list(SHH_G_P_2622,SHH_G_R_2622,SHH_L_P_3661,SHH_L_R_3661)

#compute normalization on each file and run PCA to determine number of significant components
for (i in 1:length(SCT.list)) {
	SCT.list[[i]] <- SCTransform(SCT.list[[i]], vars.to.regress = c("percent.mt","S.Score","G2M.Score"),variable.features.n=3000,return.only.var.genes=FALSE)

	SCT.list[[i]] <- RunPCA(SCT.list[[i]],assay="SCT")
	pdf(paste(my_path,ID[i],"_scRNAseq_elbow_plot.pdf",sep=""))
	print(ElbowPlot(SCT.list[[i]]))
	dev.off()
}

# elbow plot을 확인해서 elbowvalue를 정한다. G4 SHH 각각
#modify with the appropriate elbow values
elbowvalue <- list(c(1:14), c(1:17), c(1:16), c(1:16), c(1:11), c(1:16))
elbowvalue <- list(c(1:11), c(1:16), c(1:8), c(1:12))

# compute clustering and dimension reduction and some plots
for (i in 1:length(SCT.list)) {
	SCT.list[[i]] <- FindNeighbors(SCT.list[[i]], reduction = "pca", dims = elbowvalue[[i]])
	SCT.list[[i]] <- FindClusters(SCT.list[[i]], resolution = 0.3, algorithm=2)
	SCT.list[[i]] <- RunUMAP(SCT.list[[i]], reduction = "pca", dims = elbowvalue[[i]])

	pdf(paste(my_path,ID[i],"_scRNAseq_PCA.pdf",sep=""))
	p <- DimPlot(SCT.list[[i]],reduction = "pca",label=TRUE, group.by="Phase")
	print(p+ theme(axis.ticks=element_blank(),axis.text=element_blank()))
	p <- DimPlot(SCT.list[[i]],reduction = "pca",label=TRUE)
	print(p+ theme(axis.ticks=element_blank(),axis.text=element_blank()))
	dev.off()

	pdf(paste(my_path,ID[i],"_scRNAseq_UMAP.pdf",sep=""))
	p <- DimPlot(SCT.list[[i]],reduction = "umap",label=TRUE, group.by="Phase")
	print(p+ theme(axis.ticks=element_blank(),axis.text=element_blank()))
	p <- DimPlot(SCT.list[[i]],reduction = "umap",label=TRUE)
	print(p+ theme(axis.ticks=element_blank(),axis.text=element_blank()))
	dev.off()

	pdf(paste(my_path,ID[i],"_scRNAseq_Tree.pdf",sep=""))
	SCT.list[[i]]=BuildClusterTree(SCT.list[[i]],dims=elbowvalue[[i]])
	print(PlotClusterTree(SCT.list[[i]]))
	dev.off()

	# plot some markers
	#markers=c("TAGLN","ACTA2","IGFBP3","ENG","PDE6H","CRX","CD74","CXCR4","TOP2A","CENPF","CALB1","CADM2","CYP26A1","GREM1")
	#markers=c("CHRNA3","CHRNA4")
	#plots <- VlnPlot(SCT.list[[i]],features=markers, pt.size = 0,combine=FALSE)
	#for(j in 1:length(plots)) {
	#	plots[[j]] <- plots[[j]] +geom_boxplot(width=0.1,fill="white",outlier.size=0) + theme(legend.position = "none",axis.text.y = element_text(size = 8), axis.title.y=element_text(size=8),axis.text.x = element_text(size = 8, angle=45)) + labs(x="")
	#}

	#pdf(paste(my_path,ID[i],"_scRNAseq_markers_plot.pdf",sep=""),width=15)
	#print(CombinePlots(plots))
	#dev.off()
}

# compute differential analysis
for (i in 1:length(SCT.list)) {
allmarkers <- FindAllMarkers(SCT.list[[i]], only.pos = TRUE, min.pct = 0.25, test.use = "MAST")
write.table(allmarkers, paste(my_path, ID[i], "_scRNAseq_markers_clusters.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}


# save files
for (i in 1:length(SCT.list)) {
	saveRDS(SCT.list[[i]], file = paste(my_path, ID[i], "_scRNAseq.rds", sep = ""))
	write.table(SCT.list[[i]]@assays$SCT@scale.data, paste(my_path, ID[i], "_scRNAseq_scale_data.txt", sep = ""), quote = FALSE)
	write.table(SCT.list[[i]]@meta.data$seurat_clusters, paste(my_path, ID[i], "_scRNAseq_seuratclusters.txt", sep=""), quote = FALSE)
	saveRDS(SCT.list[[i]]@reductions$umap@cell.embeddings, file = paste(my_path, ID[i], "_scRNAseq_umap_embeddings.rds", sep = ""))
	mat <- as.data.frame(as.matrix(SCT.list[[i]]@assays$SCT@data))
	saveRDS(mat, file = paste(my_path, ID[i], "_scRNAseq_norm_data.rds", sep = ""))
}
