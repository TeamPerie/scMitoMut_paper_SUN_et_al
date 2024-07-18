## Load library
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggvenn)
library(readr)

## Processed so atac data
#+ eval = FALSE
so.data <- ReadMGATK(dir = "./data/mgatk_out/final/")
write_rds(so.data, "./tmp/raw_mgatk_10x_brain.rds")

so_mito <- CreateSeuratObject(
  counts = so.data$counts,
  meta.data = so.data$depth,
  assay = "mito"
)

so = read_rds("./tmp/so.rds")
so[["mito"]] = so_mito[, colnames(so)][["mito"]]
write_rds(so, "./tmp/raw_mito_so.rds")


#+ eval = TRUE
so = read_rds("./tmp/raw_mito_so.rds")
so.data <- read_rds("./tmp/raw_mgatk_10x_brain.rds")
# Calculate VMR and strand concordance using loci with coverage >= 10
variable.sites <- IdentifyVariants(so, assay = "mito", refallele = so.data$refallele)
# number of variable sites
nrow(variable.sites)
# plot the VMR and strand concordance distribution of the variable sites
VariantPlot(variants = variable.sites)
# min.cells = 5, concordance = 0.65, vmr.threshold = 0.01)
# filtering the variable sites by strand concordance and VMR
high.conf <- subset(
    variable.sites, subset = n_cells_conf_detected >= 10 & 
        strand_correlation >= 0.65 &
        vmr > 0.01
)

# NOTE: No high confidence sites
print(high.conf)

## Save high confidence variable sites
# d = data.table(high.conf, keep.rownames = T)
# write_tsv(d, "./tmp/table_signac_high_confidence.tsv")

## Mutant cells
# DefaultAssay(so) <- "mito"
# so <- AlleleFreq(
  # object = so,
  # variants = high.conf$variant,
  # assay = "mito"
# )
# DefaultAssay(so) <- "alleles"

# DoHeatmap(so, features = rownames(so), slot = "data", disp.max = 1) +
  # scale_fill_viridis_c()

## Save mutant cells
# so@assays
# x = so@assays$@data[rownames(so), ] %>% as.matrix
# x = melt(x)
# names(x) = c("mutant", "cell_barcode", "af")

# write_tsv(x, "./tmp/raw_signac_vaf.tsv")


