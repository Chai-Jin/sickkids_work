# chr arm modules ---------------------------------------------------------------
# conda activate r_general
# R
library(Seurat)
library(magrittr)
library(glue)
library(qs)
library(ggplot2)
source("src/scripts/utils.R")

# inputs
input_stamp <- snakemake@input[[1]]
sample <- input_stamp %>%
  basename() %>%
  stringr::str_replace(".stamp", "") %>%
  stringr::str_replace("initial_object_", "")

# directories
outdir <- "out/chr_arm_modules"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# get chrom arm info
cytoband_dest <- "src/metadata/cytoband_hg38.tsv"
if (!file.exists(cytoband_dest)) {
  cytoband_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
  download.file(cytoband_url, glue::glue("{cytoband_dest}.gz"))
  R.utils::gunzip(glue::glue("{cytoband_dest}.gz"))
}

cytob <- data.table::fread(cytoband_dest, col.names = c("chr", "start", "end", "band_name", "stain")) %>%
  dplyr::filter(chr %in% glue::glue("chr{c(1:22, 'X', 'Y')}")) %>%
  dplyr::mutate(arm = stringr::str_sub(band_name, 1, 1))

# get gtf
gtf <- readr::read_tsv("src/metadata/hg38_gencode_v27.txt", col_names = c("gene", "chr", "start", "end")) # from InferCNV

# declare iteration groups --> chroms and arms
chr_iter_groups <- c(glue::glue("chr{c(1:22, 'X')}"), glue::glue("chr{c(1:22, 'X')}p"), glue::glue("chr{c(1:22, 'X')}q"))


# tweak from loop
i <- sample

# plot name
plot_name <- glue::glue("{outdir}/{i}.pdf")

if (!file.exists(plot_name)) {
  # import
  so <- qs::qread(glue::glue("out/initial_objects/{i}.qs"))
  DefaultAssay(so) <- "SCT"
  message(glue::glue("Read in sample: {i}"))


  # plot loop
  pl <- list()
  for (x in chr_iter_groups) {
    if (!grepl("p|q", x)) {
      # whole chromosomes
      gn_set <- gtf %>%
        dplyr::filter(chr == x) %>%
        dplyr::pull(gene) %>%
        .[. %in% rownames(so)]
      so %<>% AddModuleScore(list(gn_set), name = glue::glue("{x}.module"))
      message(glue::glue("{i}: Gene set module added for: {x}"))
    } else {
      # chrom arms
      iter_chrom <- stringr::str_sub(x, 1, -2)
      iter_arm <- stringr::str_sub(x, -1, -1)
      iter_start <- cytob %>%
        dplyr::filter(chr == iter_chrom & arm == iter_arm) %>%
        dplyr::slice_min(start) %>%
        dplyr::pull(start)
      iter_end <- cytob %>%
        dplyr::filter(chr == iter_chrom & arm == iter_arm) %>%
        dplyr::slice_max(end) %>%
        dplyr::pull(end)

      gn_set <- gtf %>%
        dplyr::filter(chr == iter_chrom & start > iter_start & end < iter_end) %>%
        dplyr::pull(gene) %>%
        .[. %in% rownames(so)]

      tryCatch(
        {
          so %<>% AddModuleScore(list(gn_set), name = glue::glue("{x}.module"))
          message(glue::glue("{i}: Gene set module added for: {x}"))
        },
        error = function(e) {
          print(glue::glue("ERROR: {conditionMessage(e)}, skipping plot for {x}."))
        }
      )
    }

    # return featureplot
    fp_label <- x %>% stringr::str_replace("chr", "Chr. ")
    tryCatch(
      {
        pl[[x]] <- FeaturePlot(so, glue::glue("{x}.module1"), reduction = "umap", cols = rev(RColorBrewer::brewer.pal(9, "RdBu")), raster = TRUE) + NoLegend() + NoAxes() + ggtitle(glue::glue("{fp_label}")) # max.cutoff = "q95", min.cutoff = "q5"
        message(glue::glue("{i}: Plot added for: {x}"))
      },
      error = function(e) {
        print(glue::glue("ERROR: {conditionMessage(e)}, skipping plot for {x}."))
      }
    )
  }

  # test plot
  # pdf(glue::glue("{outdir}/test_{i}_{x}.pdf"), w = 5, h = 5)
  # print(pl[[1]])
  # dev.off()


  # print
  pdf(plot_name, w = 20, h = 56)
  print(cowplot::plot_grid(plotlist = pl, ncol = 5))
  dev.off()
}


# Stamp -------------------------------------------------------------------
file.create(snakemake@output[[1]])
message(glue::glue("@Stamped!"))
