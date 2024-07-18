## Load library
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggvenn)
library(readr)

## Processed crc atac data
crc = read_rds("./data/raw_signac_crc.rds")
mito.data = read_rds("./data/raw_mito_data.rds")

# Calculate VMR and strand concordance using loci with coverage >= 10
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
# number of variable sites
nrow(variable.sites)
# plot the VMR and strand concordance distribution of the variable sites
VariantPlot(variants = variable.sites)
write_rds(variable.sites, "./tmp/variable_sites.rds")
# In variable.sites, the n_cells_conf_detected is defined by:
# https://github.com/stuart-lab/signac/blob/8ecdde291676102bb3b503f48926c993354b5471/R/mito.R#L686
# At least 2 reads in fwd and rev strand
# (fwd.counts >= 2) + (rev.counts >= 2)
# min.cells = 5, concordance = 0.65, vmr.threshold = 0.01)
# filtering the variable sites by strand concordance and VMR
high.conf <- subset(
    variable.sites, subset = n_cells_conf_detected >= 5 &
        strand_correlation >= 0.65 &
        vmr > 0.01
)

## Save high confidence variable sites
d = data.table(high.conf, keep.rownames = T)
write_tsv(d, "./tmp/table_signac_high_confidence.tsv")

## Mutant cells
crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
DefaultAssay(crc) <- "alleles"

DoHeatmap(crc, features = rownames(crc), slot = "data", disp.max = 1) +
  scale_fill_viridis_c()

## Save mutant cells
x = crc@assays$alleles@data[rownames(crc), ] %>% as.matrix
x = melt(x)
names(x) = c("mutant", "cell_barcode", "af")

write_tsv(x, "./tmp/raw_signac_vaf.tsv")
