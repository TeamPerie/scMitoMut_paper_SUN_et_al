library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Seurat)
library(Signac)
library(ggvenn)
library(ggpubr)
library(devtools)

install.packages(c("data.table", "ggplot2", "plyr", "magrittr", "readr", "stringr", "Seurat", "Signac", "ggvenn", "ggpubr"))
BiocManager::install("Signac")

load_all("~/projects/scMitoMut/scMitoMut-devel/")



## Function
## BB, Signac, MQuad
colors_methods = c("BB" = "#CE5374", "Signac" = "grey", "MQuad" = "#4F9D9D")
## BI, BM, BB
colors_models = c("BI" = "#0073C2FF", "BM" = "#EFC000FF", "BB" = "#CE5374")

colors_binary = c(
    "Epithelial" = "#f28482",
    "Blood" = "#f6bd60")


## Preapre data
# BB scMitoMut obj
so <- open_h5_file("./tmp/mut.h5")

# Cell annotation
crc = read_rds("./data/raw_signac_crc.rds")
x = Idents(crc) 
# Access the UMAP coordinates
DimPlot(crc)
umap_coords <- crc@reductions$umap@cell.embeddings %>% data.table(keep.rownames = TRUE)
d_cell_type = data.table(cell_id = names(x), cell_type = x)
d_cell_ann = merge(d_cell_type, umap_coords, by.x = "cell_id", by.y = "rn", all = T)
unique(d_cell_ann$cell_type)
d_cell_ann[, cell_type_bi := revalue(cell_type, c(
        "Epithelial" = "Epithelial",
        "Basophil" = "Blood",
        "Myeloid_2" = "Blood",
        "Myeloid_1" = "Blood",
        "Tcell" = "Blood"
        ))]


d_bb = fread("tmp/raw_crc_scmitomut_bb.tsv")
d_bm = fread("tmp/raw_crc_scmitomut_bm.tsv")
d_bi = fread("tmp/raw_crc_scmitomut_bi.tsv")

## Heatmap for mutations BB
m_count = filter_loc(
    mtmutObj = so, 
    min_cell = 5, 
    model = "bb",
    p_threshold = 0.01, 
    # alt_count_threshold = 3,
    p_adj_method = "fdr"
)
cell_ann_binary = data.frame(d_cell_ann[, .(CellType = as.character(cell_type_bi))])
rownames(cell_ann_binary) = d_cell_ann$cell_id
ann_colors_binary = list("CellType" = colors_binary)
plot_heatmap(m_count, cell_ann = cell_ann_binary, ann_colors = ann_colors_binary, type = "af", percent_interp= 0.2, n_interp = 3)
plot_heatmap(m_count, cell_ann = cell_ann_binary, ann_colors = ann_colors_binary, type = "p", percent_interp= 0.2, n_interp = 3)

## Heatmap for mutations BI
# m_count = filter_loc(
#     mtmutObj = so, 
#     min_cell = 5, 
#     model = "bi",
#     p_threshold = 0.01, 
#     alt_count_threshold = 3,
#     p_adj_method = "fdr"
# )
# cell_ann_binary = data.frame(d_cell_ann[, .(CellType = as.character(cell_type_bi))])
# rownames(cell_ann_binary) = d_cell_ann$cell_id
# ann_colors_binary = list("CellType" = colors_binary)
# plot_heatmap(m_count, cell_ann = cell_ann_binary, ann_colors = ann_colors_binary, type = "af", percent_interp= 0.2, n_interp = 3)
# plot_heatmap(m_count, cell_ann = cell_ann_binary, ann_colors = ann_colors_binary, type = "p", percent_interp= 0.2, n_interp = 3)


## Comparing loci
l_plot = list(
    # BI = unique(d_bi$loc),
    BM = unique(d_bm$loc),
    BB = unique(d_bb$loc)
    )

ggvenn(
  l_plot, fill_color = as.character(colors_models),
  stroke_size = 1.5, set_name_size = 5,
  show_outside = "none", show_percent = F
  )

## Comparing lineage precision
calc_precision = function(cell_type, mut_status) {
    tab = table(cell_type == "Epithelial", mut_status)
    max(tab[, 2]) / sum(tab[, 2])
}
setkey(d_cell_ann, cell_id)
d_pres_bb = d_bb[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
d_pres_bi = d_bi[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
d_pres_bm = d_bm[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]

# d_pres_bi[, method := "BI"]
d_pres_bm[, method := "BM"]
d_pres_bb[, method := "BB"]
# d_plot = rbind(d_pres_bi, d_pres_bm, d_pres_bb)
d_plot = rbind(d_pres_bm, d_pres_bb)
# d_plot$method = factor(d_plot$method, levels = c("BI", "BM", "BB"))
d_plot$method = factor(d_plot$method, levels = c("BM", "BB"))

comparison = list(
    # c("BI", "BM"),
    # c("BI", "BB"),
    c("BM", "BB")
    )

ggplot(d_plot, aes(x = method, y = precision, color = method)) +
    stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_models) +
    labs(y = "Lineage Precision", x = "Model") +
    theme_bw() + ylim(0.5, 1.1) 

ggplot(d_plot, aes(x = method, y = precision, color = method)) +
    geom_boxplot() +
    # stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    labs(y = "Lineage Precision", x = "Model") +
    scale_color_manual(values = colors_models) +
    theme_bw()



## Comparing clone size
d_clone_size_bi = d_bi[mut_status == TRUE, .(clone_size = .N), by = .(loc)]
d_clone_size_bm = d_bm[mut_status == TRUE, .(clone_size = .N), by = .(loc)]
d_clone_size_bb = d_bb[mut_status == TRUE, .(clone_size = .N), by = .(loc)]

# d_clone_size_bi[, method := "BI"]
d_clone_size_bm[, method := "BM"]
d_clone_size_bb[, method := "BB"]

# d_plot = rbind(d_clone_size_bi, d_clone_size_bm, d_clone_size_bb)
d_plot = rbind(d_clone_size_bm, d_clone_size_bb)
# d_plot$method = factor(d_plot$method, levels = c("BI", "BM", "BB"))
d_plot$method = factor(d_plot$method, levels = c("BM", "BB"))

comparison = list(
    # c("BI", "BM"),
    # c("BI", "BB"),
    c("BM", "BB")
    )

ggplot(d_plot, aes(x = method, y = log10(clone_size), color = method)) +
    stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_models) +
    theme_bw() + ylim(0, 3.1)

ggplot(d_plot, aes(x = method, y = clone_size, color = method)) +
    geom_boxplot() +
    # stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_models) +
    theme_bw() + scale_y_log10() + labs(y = "Clone Size", x = "Model")
