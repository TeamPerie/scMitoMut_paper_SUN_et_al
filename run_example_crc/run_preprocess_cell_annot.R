library(Signac)
library(readr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)

# Ref:
# https://github.com/stuart-lab/signac/blob/6cff061d29a0b2231d16ef49f6f5f1b9647884c2/vignettes/mito.Rmd

counts <- Read10X_h5(filename = "./data/cellranger_output/CRC_v12-mtMask_mgatk.filtered_peak_bc_matrix.h5")
metadata <- read.csv(
    file = "./data/cellranger_output/CRC_v12-mtMask_mgatk.singlecell.csv",
    header = TRUE,
    row.names = 1
)

# load gene annotations from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
# write_rds(annotations, "./data/annotations.rds")
annotations = read_rds("./data/annotations.rds")

# create object
crc_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  annotation = annotations,
  min.cells = 10,
  genome = "hg19",
  fragments = './data/cellranger_output/CRC_v12-mtMask_mgatk.fragments.tsv.gz'
)
crc <- CreateSeuratObject(
  counts = crc_assay,
  assay = 'peaks',
  meta.data = metadata
)

# Augment QC metrics that were computed by cellranger-atac
crc$pct_reads_in_peaks <- crc$peak_region_fragments / crc$passed_filters * 100
crc$pct_reads_in_DNase <- crc$DNase_sensitive_region_fragments / crc$passed_filters * 100
crc$blacklist_ratio <- crc$blacklist_region_fragments / crc$peak_region_fragments

# compute TSS enrichment score and nucleosome banding pattern
crc <- TSSEnrichment(crc)
crc <- NucleosomeSignal(crc)

# visualize QC metrics for each cell
VlnPlot(crc, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks", "pct_reads_in_DNase", "blacklist_ratio"), pt.size = 0, ncol = 3)

# remove low-quality cells
crc <- subset(
  x = crc,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 50000 &
    pct_reads_in_DNase > 40 &
    blacklist_ratio < 0.05 &
    TSS.enrichment > 3 & 
    nucleosome_signal < 4
)
crc

# load mgatk output
mito.data <- ReadMGATK(dir = "./data/mgatk_output/final/")
write_rds(mito.data, "./data/raw_mito_data.rds")

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assat
mito <- subset(mito, cells = colnames(crc))

# add assay and metadata to the seurat object
crc[["mito"]] <- mito
crc <- AddMetaData(crc, metadata = mito.data$depth, col.name = "mtDNA_depth")

VlnPlot(crc, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()

# filter cells based on mitochondrial depth
crc <- subset(crc, mtDNA_depth >= 10)
crc

set.seed(123)
crc <- RunTFIDF(crc)
crc <- FindTopFeatures(crc, min.cutoff = 10)
crc <- RunSVD(crc)
crc <- RunUMAP(crc, reduction = "lsi", dims = 2:50)
crc <- FindNeighbors(crc, reduction = "lsi", dims = 2:50)
crc <- FindClusters(crc, resolution = 0.5, algorithm = 3)

DimPlot(crc, label = TRUE) + NoLegend()

# compute gene accessibility
gene.activities <- GeneActivity(crc)

# add to the Seurat object as a new assay
crc[['RNA']] <- CreateAssayObject(counts = gene.activities)

crc <- NormalizeData(
  object = crc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(crc$nCount_RNA)
)

DefaultAssay(crc) <- 'RNA'

FeaturePlot(
  object = crc,
  features = c('TREM1', 'EPCAM', "PTPRC", "IL1RL1","GATA3", "KIT"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)

crc <- RenameIdents(
  object = crc,
  '0' = 'Epithelial',
  '1' = 'Epithelial',
  '2' = 'Basophil',
  '3' = 'Myeloid_1',
  '4' = 'Myeloid_2',
  '5' = 'Tcell'
)

colors = c(
    Tcell = "#f6bd60", 
    Basophil = "#9b948d",
    Myeloid_1 = "#f5cac3",
    Myeloid_2 = "#84a59d",
    Epithelial = "#f28482")
DimPlot(crc, reduction = "umap", cols = colors)

p1 <- FeatureScatter(crc, "mtDNA_depth", "pct_reads_in_peaks") + ggtitle("") + scale_x_log10()
p2 <- FeatureScatter(crc, "mtDNA_depth", "nucleosome_signal") + ggtitle("") + scale_x_log10()

p1 + p2 + plot_layout(guides = 'collect')

variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites)

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

high.conf[,c(1,2,5)]
write_rds(crc, "./data/raw_signac_crc.rds")
