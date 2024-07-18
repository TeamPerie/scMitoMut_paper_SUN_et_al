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
library(pheatmap)
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
so <- open_h5_file("./tmp/bj_mkn45_1pct.h5")


## BJ and MKN45 cell line identification using copy number data
x = fread("./data/cnv/bj_mkn45_1pct_dloupe-group_2110-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv")
node_id = x[[1]]
cell_id = x[[2]] %>% str_split(";")
cell_rep = rep(1:length(cell_id), sapply(cell_id, length))
cell_id_v = unlist(cell_id)
node_id_v = node_id[cell_rep]
d_cell_ann = data.table(
    cell_id = cell_id_v,
    cell_type_bi = ifelse(node_id_v == "2105", "MKN45", "BJ")
    )


## Loading data from bb, bm and bi model
d_bb = fread("./tmp/bj_mkn45_1pct_bb.tsv")
d_bm = fread("./tmp/bj_mkn45_1pct_bm.tsv")
d_bi = fread("./tmp/bj_mkn45_1pct_bi.tsv")

## Loading signac data and mquad data
d_s = fread("./tmp/table_signac_high_confidence.tsv")
names(d_s)[1] = "mutant"
mut_list_s = unique(d_s$mutant) %>% str_split_fixed(".>", 2) %>% data.table
colnames(mut_list_s) = c("loc", "base")
mut_list_s$mutant = unique(d_s$mutant)
x = lapply(1:nrow(mut_list_s), function(i) {
    x = read_locus(so, paste0("chrM.", mut_list_s[i, 1]), mut_list_s[i, 2])
    x$mutant = mut_list_s[i, 3]
    x
}) %>% rbindlist
d_s = merge(d_s, x, by = c("mutant"), all = T)
d_s[, mut_status := fwd_depth > 1 & rev_depth > 1]
d_s[, loc := paste0("chrM.", gsub("...$", "", mutant))]
d_s = d_s[!is.na(mut_status)]

d_m = fread("./tmp/table_mquad_high_confidence.tsv")

## Comparing the mutation loci
l_plot = list(
    # BI = unique(d_bi$loc),
    BB = unique(d_bb$loc),
    MQuad = unique(d_m$loc),
    Signac = unique(d_s$loc)
    )

ggvenn(
  l_plot, fill_color = as.character(colors_methods),
  stroke_size = 1.5, set_name_size = 5,
  show_outside = "none", show_percent = F
  )

ggsave("./tmp/fig_venn_loci.png", width = 5, height = 5)

## count the mutation at least one cell in MKN45
v_loc_bb = d_bb[mut_status == T & cell_barcode %in% cell_id_v[node_id_v == "2105"], unique(loc)]
v_loc_signac = d_s[mut_status == T & cell_barcode %in% cell_id_v[node_id_v == "2105"], unique(loc)]
v_loc_mquad = d_m$loc

## Comparing the mutation loci
l_plot = list(
    # BI = unique(d_bi$loc),
    BB = unique(v_loc_bb),
    Signac = unique(v_loc_signac),
    MQuad = unique(v_loc_mquad)
    )

ggvenn(
  l_plot, fill_color = as.character(colors_methods),
  stroke_size = 1.5, set_name_size = 5,
  show_outside = "none", show_percent = F
  )

ggsave("./tmp/fig_venn_loci.png", width = 5, height = 5)



## Comparing lineage precision
calc_precision = function(cell_type, mut_status) {
    tab = table(cell_type == "BJ", mut_status)
    max(tab[, 2]) / sum(tab[, 2])
}

setkey(d_cell_ann, cell_id)
d_pres_bb = d_bb[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
d_pres_bm = d_bm[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
d_pres_s = d_s[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]

d_pres_bb[, method := "BB"]
d_pres_bm[, method := "BM"]
d_pres_s[, method := "Signac"]
# d_pres_m[, method := "MQ"]
# d_plot = rbind(d_pres_c, d_pres_s, d_pres_m)
d_plot = rbind(d_pres_c, d_pres_s)
d_plot$method = factor(d_plot$method, levels = c("scMitoMut", "Signac"))

comparison = list(
    c("scMitoMut", "Signac")
    # c("scMitoMut", "MQ"),
    # c("Signac", "MQ")
    )

ggplot(d_plot, aes(x = method, y = precision, color = method)) +
    geom_boxplot() +
    stat_compare_means(method = "wilcox", comparisons = comparison) +
    scale_color_manual(values = colors_methods) +
    geom_jitter(stat = "identity") +
    ylim(0.5, 1.1) +
    theme_bw()

ggplot(d_plot, aes(x = method, y = precision, color = method)) +
    geom_boxplot() +
    # stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_methods) +
    labs(y = "Lineage Precision", x = "Method") +
    theme_bw()

d_plot = dcast(d_plot, loc ~ method, value.var = "precision", fill = 0)
ggplot(d_plot) + aes(x = scMitoMut, y = Signac) + geom_jitter(width = 0.01, height = 0.01) + geom_abline(intercept = 0, slope = 1) + theme_bw()

## Comparing clone size
d_clone_size_c = d_c[mut_status == TRUE, .(clone_size = .N), by = .(loc)]
d_clone_size_s = d_s[mut_status == TRUE, .(clone_size = .N), by = .(loc)]
# d_clone_size_m = d_m[mut_status == TRUE, .(clone_size = .N), by = .(loc)]

d_clone_size_c[, method := "scMitoMut"]
d_clone_size_s[, method := "Signac"]
# d_clone_size_m[, method := "MQ"]

# d_plot = rbind(d_clone_size_c, d_clone_size_s, d_clone_size_m)
d_plot = rbind(d_clone_size_c, d_clone_size_s)
# d_plot$method = factor(d_plot$method, levels = c("scMitoMut", "Signac", "MQ"))
d_plot$method = factor(d_plot$method, levels = c("scMitoMut", "Signac"))

comparison = list(
    c("scMitoMut", "Signac")
    )

ggplot(d_plot, aes(x = method, y = clone_size, color = method)) +
    geom_boxplot() +
    stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_methods) +
    theme_bw() + scale_y_log10(limits = c(1, 2000))

ggplot(d_plot, aes(x = method, y = clone_size, color = method)) +
    geom_boxplot() +
    # stat_compare_means(method = "wilcox", comparisons = comparison) +
    geom_jitter(stat = "identity") +
    scale_color_manual(values = colors_methods) +
    theme_bw() + scale_y_log10() +
    labs(y = "Clone Size", x = "Method")

d_plot = dcast(d_plot, loc ~ method, value.var = "clone_size", fill = 0)
d_plot[Signac == 0][order(scMitoMut)]
ggplot(d_plot) + aes(x = scMitoMut + 1, y = Signac + 1) + geom_jitter(width = 0.01, height = 0.01) + geom_abline(intercept = 0, slope = 1) + 
    theme_bw() + scale_x_log10() + scale_y_log10()



## Compare VMR and strand concordance
variable.sites = read_rds("./tmp/variable_sites.rds")

# get the concordance and vmr of the variable sites
d_plot = data.table(variable.sites)[n_cells_conf_detected >= 2
    , .(vmr = vmr[which.max(strand_correlation)], strand_correlation = max(strand_correlation)), 
    by = position] %>% na.omit
d_plot[, loc := paste0("chrM.", position)]

ggplot() + aes(x = strand_correlation, y = vmr) +
    geom_point(data = d_plot, alpha = 0.2) +
    # geom_point(data = d_plot[class != "negative"], aes(color = class)) +
    scale_y_log10() +
    theme_bw() + 
    geom_hline(yintercept = 0.01, linetype = "dashed", col = 'red') +
    geom_vline(xintercept = 0.65, linetype = "dashed", col = 'red') +
    labs(x = "Strand concordance", y = "VMR")

# VMR of BB model called mutations
loc_v = l_plot$scMitoMut

d_plot[, called_scMitoMut := "No"]
d_plot[loc %in% loc_v, called_scMitoMut := "Yes"]

ggplot(d_plot) + aes(x = called_scMitoMut, y = vmr, color = called_scMitoMut) +
    # geom_boxplot() +
    geom_jitter() +
    theme_bw() + 
    scale_y_log10() +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0.01, linetype = "dashed", col = 'red') +
    labs(x = "Identified by BB", y = "VMR")

ggplot(d_plot) + aes(x = called_scMitoMut, y = strand_correlation, color = called_scMitoMut) +
    # geom_boxplot() +
    geom_jitter() +
    theme_bw() + 
    scale_y_log10() +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0.65, linetype = "dashed", col = 'red') +
    labs(x = "Identified by BB", y = "Strand concordance")



## Comparing lineage precision

## Comparing clone size



# x_m = fread("./tmp/bj_mkn45_1pct_bb.tsv")
# x_m[mut_status == T, .N, loc]
# x_m[mut_status == T, mean(af), loc]
#
#
#
#
#
# y_cell = x_m[loc == "chrM.9903" & mut_status == T, cell_barcode]
# y_cell = x_m[loc == "chrM.11226" & mut_status == T, cell_barcode]
# y_cell = x_m[loc == "chrM.60" & mut_status == T, cell_barcode]
#
# node_id_v[cell_id_v %in% y_cell]
#
# cell_id_v_sub = cell_id_v[node_id_v == "2105"]
#
#
# # x_c = fread("./tmp/BAF_bin_phased_adata_obs.csv")
# # x_cell_c = x_c[sub24_2_r1 == "Cluster 1", index]
# # x_c = fread("./tmp/mkn45_5k_750k_cnv_cluster_id.tsv")
#
# # x_cell_c = x_c[cluster_id == 1, cell_id]
# # x_cell_m = x_m[loc %in% c("chrM.2393") & mut_status == T, cell_barcode]
#
# # x_af[loc == "chrM.2393" & mut_status == T, ]
#
#
#
#
h5_obj = open_h5_file("./tmp/bj_mkn45_1pct.h5")

m_count = filter_loc(
    mtmutObj = h5_obj, 
    min_cell = 5, 
    model = "bb",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
m_count

d_bb = export_dt(m_count, all_cell = T)

# choose mutation with at least one mutation in MKN45 cells
v_loc = d_bb[mut_status == T & cell_barcode %in% cell_id_v[node_id_v == "2105"], unique(loc)]
m_count$loc_pass = v_loc
d_bb_sub = d_bb[loc %in% v_loc]

v_loc_shared = c(v_loc_signac, v_loc_mquad)

ann_colors = list(node_id = c("BJ" = "#f28482", "MKN45" = "#f6bd60"))
cell_ann = data.frame(`Cell Type` = ifelse(node_id_v == "2105", "MKN45", "BJ"))
rownames(cell_ann) = cell_id_v

## Plot the p value heatmap
m = export_pval(m_count, all_cell = T)
mut_ann = data.frame(`Mutation` = ifelse(rownames(m) %in% v_loc_shared, "Shared", "Unique"))
rownames(mut_ann) = rownames(m)
i_mut_order = c(v_loc_shared, setdiff(v_loc, v_loc_shared))
m = m[i_mut_order, ]
pheatmap::pheatmap(m,
    color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "GnBu")))(100)),
    show_colnames = FALSE, annotation_col = cell_ann, cluster_cols = FALSE,
    cluster_rows = FALSE, annotation_colors = ann_colors, na_col = "white",
    annotation_row = mut_ann
)

## Plot the vaf heatmap
m = export_af(m_count, all_cell = T)
mut_ann = data.frame(`Mutation` = ifelse(rownames(m) %in% v_loc_shared, "Shared", "Unique"))
rownames(mut_ann) = rownames(m)
i_mut_order = c(v_loc_shared, setdiff(v_loc, v_loc_shared))
m = m[i_mut_order, ]
pheatmap::pheatmap(m,
    color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "GnBu")))(100)),
    show_colnames = FALSE, annotation_col = cell_ann, cluster_cols = FALSE,
    cluster_rows = FALSE, annotation_colors = ann_colors, na_col = "white",
    annotation_row = mut_ann
)


