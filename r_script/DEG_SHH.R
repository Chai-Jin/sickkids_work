library(dplyr)
library(Seurat)
library(reticulate)
library(harmony)
library(ggplot2)
library(Matrix)
library(ape)
library(future)
library(ggrepel)
library(writexl)

my_path <- "/Users/cj/Desktop/PUGH/"

ID=c("SHH_G_P_2622","SHH_G_R_2622","SHH_L_P_3661","SHH_L_R_3661")

custom_colors <- list()
colors_dutch <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471','#EE5A24','#009432','#0652DD','#9980FA','#833471','#EA2027','#006266','#1B1464','#5758BB','#6F1E51')
colors_spanish <- c('#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2','#2c2c54','#474787','#aaa69d','#227093','#218c74','#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79','#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
custom_colors$discrete <- c(colors_dutch, colors_spanish)

SHH_whole_final=readRDS(paste(my_path,"Integrated_SHH_MBs_scRNAseq.rds",sep=""))

# without harmony
SHH_whole_final <- FindNeighbors(SHH_whole_final, assay="SCT", dims= 1:20)
SHH_whole_final <- FindClusters(SHH_whole_final,resolution = 0.4, algorithm=2,assay="SCT")
SHH_whole_final <- RunUMAP(SHH_whole_final, assay="SCT",dims=1:20)

pdf(paste(my_path,"Integrated_SHH_MBs_scRNAseq_PCA_norm.pdf",sep=""))
p=DimPlot(SHH_whole_final, group.by="orig.ident",label=TRUE)
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,label=TRUE, group.by="Phase")
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,label=TRUE)
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,label=TRUE,split.by="orig.ident")
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
dev.off()

pdf(paste(my_path,"Integrated_SHH_MBs_scRNAseq_UMAP_norm.pdf",sep=""))
p=DimPlot(SHH_whole_final,reduction = "umap",label=TRUE, group.by="orig.ident",cols=custom_colors$discrete)
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,reduction = "umap",label=TRUE, group.by="Phase",cols=c('#45aaf2', '#f1c40f', '#e74c3c'))
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,reduction = "umap",label=TRUE,cols=custom_colors$discrete)
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
p=DimPlot(SHH_whole_final,reduction = "umap",label=TRUE,split.by="orig.ident",cols=custom_colors$discrete)
p+ theme(axis.ticks=element_blank(),axis.text=element_blank())
dev.off()

SHH_whole_final.markers <- FindAllMarkers(SHH_whole_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SHH_whole_final.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(SHH_whole_final, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

allmarkers <- FindAllMarkers(SHH_whole_final, only.pos=TRUE, min.pct = 0.25, test.use="MAST")

write_xlsx(allmarkers, path = paste(my_path,"SHH_All_markers.xlsx",sep=""))


Idents(SHH_whole_final)=factor(Idents(SHH_whole_final),levels=c("SHH_G_P_2622","SHH_G_R_2622","SHH_L_P_3661","SHH_L_R_3661"))
Idents(SHH_whole_final)=SHH_whole_final$orig.ident
# up-regulated in primary
DEG1 = FindMarkers(SHH_whole_final, ident.1="SHH_G_P_2622", ident.2="SHH_G_R_2622", only.pos=TRUE, min.pct = 0.25,test.use="MAST")
DEG2 = FindMarkers(SHH_whole_final, ident.1="SHH_L_P_3661", ident.2="SHH_L_R_3661", only.pos=TRUE, min.pct = 0.25,test.use="MAST")


# up-regulated in recurrent
DEG3 = FindMarkers(SHH_whole_final, ident.1="SHH_G_R_2622", ident.2="SHH_G_P_2622", only.pos=TRUE, min.pct = 0.25,test.use="MAST")
DEG4 = FindMarkers(SHH_whole_final, ident.1="SHH_L_R_3661", ident.2="SHH_L_P_3661", only.pos=TRUE, min.pct = 0.25,test.use="MAST")


VlnPlot(SHH_whole_final, features = c("MAML3", "FAM155A", "TMEM108", "CTNNA2", "DACH1", "NDST3"), pt.size=0)
FeaturePlot(SHH_whole_final, features = c("MAML3", "FAM155A", "TMEM108", "CTNNA2", "DACH1", "NDST3"))

saveRDS(SHH_whole_final, file = (paste(my_path,"Integrated_SHH_MBs_scRNAseq_DEG.rds",sep="")))
