## Load library
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggvenn)
library(readr)

so.data <- ReadMGATK(dir = "./data/atac_out_bj_mkn45_1pct/final/")
write_rds(so.data, "./tmp/signac_mt_bj_mkn45_1pct.rds")

so_mito <- CreateSeuratObject(
  counts = so.data$counts,
  meta.data = so.data$depth,
  assay = "mito"
)

variable.sites <- IdentifyVariants(so_mito, assay = "mito", refallele = so.data$refallele)

nrow(variable.sites)

## NOTE: it is 5, not 2 anymore
VariantPlot(variants = variable.sites, min.cells = 5)

## NOTE: There are 14 sites with high confidence
high.conf <- subset(
    variable.sites, subset = n_cells_conf_detected >= 5 & 
        strand_correlation >= 0.65 &
        vmr > 0.01
)

d = data.table(high.conf, keep.rownames = T)
write_tsv(d, "./tmp/table_signac_high_confidence.tsv")


