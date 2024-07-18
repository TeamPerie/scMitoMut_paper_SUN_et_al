library(data.table)
library(knitr)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

set.seed(1234)

#+ eval=F
counts <- Read10X_h5("./data/cellranger_out/human_brain_3k_filtered_feature_bc_matrix.h5")
fragpath <- "./data/cellranger_out/human_brain_3k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels
# seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
so <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
so[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

so

## Quality control
DefaultAssay(so) <- "ATAC"

so <- NucleosomeSignal(so)
so <- TSSEnrichment(so)

VlnPlot(
  object = so,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
so <- subset(
  x = so,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
so

## call peaks using MACS2
peaks <- CallPeaks(so)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(so),
  features = peaks,
  cells = colnames(so)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
so[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(so) <- "RNA"
so <- SCTransform(so)
so <- RunPCA(so)
so <- FindNeighbors(so)
so <- FindClusters(so)
so <- RunUMAP(so, reduction = "pca", dims = 1:30)

DefaultAssay(so) <- "peaks"
so <- NucleosomeSignal(so)
so <- TSSEnrichment(so)
so <- FindTopFeatures(so, min.cutoff = 5)
so <- RunTFIDF(so)
so <- RunSVD(so)

gene.activities <- GeneActivity(so)

# add to the Seurat object as a new assay
so[['gene']] <- CreateAssayObject(counts = gene.activities)

so <- NormalizeData(
  object = so,
  assay = 'gene',
  normalization.method = 'LogNormalize',
  scale.factor = median(so$nCount_RNA)
)

write_rds(so, "./tmp/so.rds")


# write_rds(so, "./tmp/so_so.rds")
#+ eval=T
so = read_rds("./tmp/so.rds")

g = DimPlot(so, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
g

d_marker = FindAllMarkers(so, min.pct = 0.25, logfc.threshold = 2)
d_marker = data.table(d_marker)
d_marker[, .(rank = rank(p_val), gene), by = cluster][rank < 5]
gene_list = data.table(d_marker)[p_val_adj < 0.01 & abs(avg_log2FC) > 2, gene]
DoHeatmap(so, features = gene_list)

gene_list = c("PLP1")
FeaturePlot(so, features = gene_list)
VlnPlot(so, features = gene_list)

## Assign PLP1+ to subcluster 0, 1, 3, 5, 6, 7, 13
write_rds(so, "./tmp/so.rds")
        

## Rename other cells as "PLP1-" cells

# ggsave("./tmp/fig_umap.png", g)

cell_type = data.table(barcode = names(Idents(so)), cell_type = Idents(so))

write_tsv(cell_type, "./tmp/cell_type.tsv")


## Reduce cell type into 3 groups
# selected intreated cell types
# so = read_rds("./tmp/so.rds")
#
# cell_annot = fread("./tmp/cell_type.tsv")
#
# cell_annot_sub = cell_annot[cell_type %in% c("CD14 Mono", "CD16 Mono", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "Treg", "CD4 CTL", "gdT", "MAIT", "B memory", "B intermediate", "B naive")]
# kable(cell_annot_sub$cell_type %>% table) %>% write("./tmp/table_so_sub_cell_type_count_3g.txt")
#
# cell_id_high_depth = fread("./tmp/mgatk_out/final/mgatk_atac.depthTable.txt")[V2 > 5, V1]
# cell_id_selected = intersect(cell_id_high_depth, cell_annot_sub$barcode)
# cell_annot_sub = cell_annot_sub[barcode %in% cell_id_selected]
# # kable(cell_annot_sub$cell_type %>% table) %>% write("./tmp/table_so_sub2_cell_type_count_3g.txt")
# write_tsv(cell_annot_sub, "./tmp/table_cell_type_sub_3g.tsv")
#
# so_sub = so[, cell_id_selected]
# DimPlot(so_sub, reduction = "umap", label = TRUE, repel = TRUE)
# write_rds(so_sub, "./tmp/so_sub_3g.rds")

########
#  Explore the cell annotaiton #
########


