# SingleR 돌리기

# Human Primary Cell Atlas data for reference
BiocManager::install("celldex")
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

# install singleR
BiocManager::install("SingleR")
library(SingleR)

# SingleR 용 data로 변환하기
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(DietSeurat(G4_whole_final))
sce

# singleR 돌리기
pred.G4 <- SingleR(test = sce, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)
table(pred.G4$pruned.labels)
G4_whole_final@meta.data$pred <- pred.G4$pruned.labels
G4_whole_final <- SetIdent(G4_whole_final, value = "pred")
DimPlot(G4_whole_final, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()