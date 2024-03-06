# anders module

salloc -N 1 -c 1 --mem 48G -t 10:00:00 #ask for 1 node, 1 core, 16G RAM, 10 hours
srun --pty bash # log in using that allocation
conda activate /hpf/largeprojects/mdtaylor/aerickson/data/clones/G4MB/src/envs/r_general
R
library(Seurat)
library(Signac)
library(magrittr)
library(RColorBrewer)

custom_colors <- list()
colors_dutch <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471','#EE5A24','#009432','#0652DD','#9980FA','#833471','#EA2027','#006266','#1B1464','#5758BB','#6F1E51')
colors_spanish <- c('#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2','#2c2c54','#474787','#aaa69d','#227093','#218c74','#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79','#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
custom_colors$discrete <- c(colors_dutch, colors_spanish)

# SHH ---
shh <- readRDS("/hpf/largeprojects/mdtaylor/chaijin/Pugh_SingleCell/SHH/Integrated_SHH/Integrated_SHH_MBs_scRNAseq.rds")
shh

# G4again
g4again <- readRDS("/hpf/largeprojects/mdtaylor/chaijin/Pugh_SingleCell/G4again/harmony_integrated_G4_MBs_scRNAseq.rds")

setwd("/hpf/largeprojects/mdtaylor/aerickson/data/clones/G4MB")
# protein coding genes ---
prot <- data.table::fread("src/metadata/protein_coding_genes_updated.tsv")
prot <- data.table::fread("/Users/chaijinlee/Research/G4again/IntegratedG4/protein_coding_genes_updated.tsv")

# get aldinger modules ---
ald_mks <- qs::qread("out/symphony_refs/markers_aldinger.qs")
ald_mks <- qs::qread("/Users/chaijinlee/Research/G4again/IntegratedG4/markers_aldinger.qs")

ald_mods <- ald_mks %>% 
  dplyr::filter(gene %in% prot$use & !grepl("^RP[SL]", gene)) %>%
  split(.$cluster) %>% 
  purrr::map(\(x) {x %>% dplyr::pull(gene)}) %>% 
  magrittr::set_names(snakecase::to_snake_case(names(.)))
ald_mod_names <- names(ald_mods) %>% paste0(1:length(ald_mods))

# assign aldinger module scores ---
g4again %<>% AddModuleScore(ald_mods, name = names(ald_mods))

# plot
pdf(glue::glue("/hpf/largeprojects/mdtaylor/chaijin/Pugh_SingleCell/G4again/aldinger_modules_RB_g4.pdf"), w = 20, h = 10)
FeaturePlot(g4again, ald_mod_names, cols = rev(RColorBrewer::brewer.pal(9, "RdBu")), raster = T, ncol = 4) & NoAxes()
dev.off()

pdf(glue::glue("/hpf/largeprojects/mdtaylor/chaijin/Pugh_SingleCell/G4again/aldinger_modules_vln3_g4.pdf"), w = 20, h = 10)
VlnPlot(g4again, ald_mod_names, cols = custom_colors$discrete, pt.size = 0, ncol = 4, group.by="orig.ident") # & NoAxes()
dev.off()

pdf(glue::glue("/hpf/largeprojects/mdtaylor/chaijin/Pugh_SingleCell/G4again/aldinger_modules_vln4_g4.pdf"), w = 20, h = 10)
VlnPlot(g4again, ald_mod_names, cols = custom_colors$discrete, pt.size = 0, ncol = 4, split.by = 'group') # & NoAxes()
dev.off()