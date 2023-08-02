
top100 <- G4_whole_final.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

#ClusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
library("clusterProfiler")
library("org.Hs.eg.db")
library("AnnotationHub")

df <- top100pval[, 7:6]
dfsample <- split(df$gene, df$cluster)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#Here there at 17 clusters (0, 1, 2, 3, ~, 16)
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$`0` <- bitr(dfsample$`0`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`1` <- bitr(dfsample$`1`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`2` <- bitr(dfsample$`2`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`3` <- bitr(dfsample$`3`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`4` <- bitr(dfsample$`4`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`5` <- bitr(dfsample$`5`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`6` <- bitr(dfsample$`6`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`7` <- bitr(dfsample$`7`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`8` <- bitr(dfsample$`8`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`9` <- bitr(dfsample$`9`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`10` <- bitr(dfsample$`10`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`11` <- bitr(dfsample$`11`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`12` <- bitr(dfsample$`12`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`13` <- bitr(dfsample$`13`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`14` <- bitr(dfsample$`14`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`15` <- bitr(dfsample$`15`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dfsample$`16` <- bitr(dfsample$`16`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#do the same here, a line like below for each cluster

genelist <- list("0" = dfsample$`0`$ENTREZID, 
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID,
                 "7" = dfsample$`7`$ENTREZID,
                 "8" = dfsample$`8`$ENTREZID,
                 "9" = dfsample$`9`$ENTREZID, 
                 "10" = dfsample$`10`$ENTREZID,
                 "11" = dfsample$`11`$ENTREZID,
                 "12" = dfsample$`12`$ENTREZID,
                 "13" = dfsample$`13`$ENTREZID,
                 "14" = dfsample$`14`$ENTREZID,
                 "15" = dfsample$`15`$ENTREZID,
                 "16" = dfsample$`16`$ENTREZID
)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
dotplot(GOclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)
