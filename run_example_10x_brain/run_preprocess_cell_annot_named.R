library(data.table)
library(knitr)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Signac)
library(Seurat)
library(stringr)

set.seed(1234)

so = read_rds("./tmp/so.rds")

FeaturePlot(so, "PLP1")

## Neuron markers
FeaturePlot(so, c("MAP2", "SYN1"))

## Astrocyte markers
FeaturePlot(so, c("GFAP", "ALDH1L1"))

## Microglia markers: IBA1, CX3CR1
FeaturePlot(so, c("P2RY12", "CX3CR1"))

## Oligodendrocyte markers
FeaturePlot(so, c("MBP", "OLIG2"))

## Endothelial Cells
FeaturePlot(so, c("VWF"))

## Pericytes
FeaturePlot(so, c("PDGFRB", "RGS5"))

Idents(so) <- "seurat_clusters"

DimPlot(so, label = TRUE)


# 11, 4, 8, 9: Neuron
# 6,1,3,5, 0, 7: Oligodendrocyte
# 12: Microglia
# 10, 2, 13: Astrocyte
# 14: ?
revalue(so[['seurat_clusters']][[1]], c(
    "11" = "Neuron_1", "4" = "Neuron_2", "8" = "Neuron_3", "9" = "Neuron_4",
    "6" = "Oligodendrocyte_1", "1" = "Oligodendrocyte_2", "3" = "Oligodendrocyte_3", "5" = "Oligodendrocyte_4", "0" = "Oligodendrocyte_5", "7" = "Oligodendrocyte_6",
    "12" = "Microglia",
    "10" = "Astrocyte_1", "2" = "Astrocyte_2", "13" = "Astrocyte_3",
    "14" = "Unknown"
)) -> so[['seurat_clusters_named']]

Idents(so) <- "seurat_clusters_named"
DimPlot(so, label = TRUE)

write_rds(so, "./tmp/so_named.rds")
