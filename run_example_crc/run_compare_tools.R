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
load_all("~/projects/scMitoMut/scMitoMut-devel/")

mtmutobj = open_h5_file("./tmp/mut.h5")

## Function
plot_qc_dot = function(umap_coords, d_bb, loc_x) {
    d_plot = merge(umap_coords, d_bb[loc == loc_x], by.x = "rn", by.y = "cell_barcode", all = F)
    d_plot = na.omit(d_plot)
    d_plot = d_plot[order(mut_status)]
    g_umap = ggplot(d_plot) + aes(x = UMAP_1, y = UMAP_2, color = mut_status) + geom_point() +
        scale_color_manual(values = c("grey", "red")) + theme_bw() + theme(legend.position = "none")
    # g_scatter = plot_af_coverage(m_count, loc_x, p_threshold = 0.01, alt_count_threshold = 4, p_adj_method = "bonferroni", model = "bb")
    # g_umap = g_umap + ggtitle(loc_x)
    # egg::ggarrange(g_umap, g_scatter, ncol = 2)
    g_umap
}

## Preapre data

# Cell annotation
crc = read_rds("./data/raw_signac_crc.rds")
x = Idents(crc) 
# Access the UMAP coordinates
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

d_c = fread("tmp/raw_crc_scmitomut_bb.tsv")

## Process Signac daat
d_s = fread("tmp/raw_signac_vaf.tsv")
mut_list_s = unique(d_s$mutant) %>% str_split_fixed(".>", 2) %>% data.table
colnames(mut_list_s) = c("loc", "base")
mut_list_s$mutant = unique(d_s$mutant)
x = lapply(1:nrow(mut_list_s), function(i) {
    x = read_locus(mtmutobj, paste0("chrM.", mut_list_s[i, 1]), mut_list_s[i, 2])
    x$mutant = mut_list_s[i, 3]
    x
}) %>% rbindlist
d_s = merge(d_s, x, by = c("mutant", "cell_barcode"), all = T)
d_s[, mut_status := fwd_depth > 1 & rev_depth > 1]
d_s[, loc := paste0("chrM.", gsub("...$", "", mutant))]
d_s = d_s[!is.na(mut_status)]

## Process MQuad data
# d_m = fread("tmp/raw_mq_vaf.tsv")
d_m = fread("./tmp/table_mquad_high_confidence.tsv")
# d_m[, loc := paste0("chrM.", gsub("...$", "", mutant))]

## Comparing loci
l_plot = list(
    scMitoMut = unique(d_c$loc),
    Signac = unique(d_s$loc),
    MQ = unique(d_m$loc)
    )

## scMitoMut specific loci
setdiff(l_plot[[1]], c(l_plot[[2]], l_plot[[3]]))

## MQuad specific loci
setdiff(l_plot[[3]], c(l_plot[[1]], l_plot[[2]]))

## BB, Signac, MQuad
colors_methods = c("#CE5374", "grey", "#4F9D9D")
## BI, BM, BB
colors_models = c("#0073C2FF", "#EFC000FF", "#CE5374")

ggvenn(
  l_plot, fill_color = as.character(colors_methods),
  stroke_size = 1.5, set_name_size = 5,
  show_outside = "none", show_percent = F
  )
ggsave("tmp/fig_venn_loci.pdf", width = 5, height = 5)

## Comparing lineage precision
calc_precision = function(cell_type, mut_status) {
    tab = table(cell_type == "Epithelial", mut_status)
    max(tab[, 2]) / sum(tab[, 2])
}
setkey(d_cell_ann, cell_id)
d_pres_c = d_c[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
d_pres_s = d_s[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]
# d_pres_m = d_m[, .(precision = calc_precision(d_cell_ann[cell_barcode, cell_type_bi], mut_status)), by = loc]

d_pres_c[, method := "scMitoMut"]
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

# Show the lineage precision of the identified mutations
# loc_x = d_plot[called_scMitoMut == "Yes" & vmr < 0.01, loc]
# plot_qc_dot(umap_coords,g d_bb, loc_x)

# loc_x = d_plot[called_scMitoMut == "Yes" & strand_correlation < 0.65, loc]
# plot_qc_dot(umap_coords, d_bb, loc_x)

# Example of 
# plot_qc_dot(umap_coords, d_bb, "chrM.1227")

# setdiff(l_plot[[1]], c(l_plot[[2]], l_plot[[3]]))
# plot_qc_dot(umap_coords, d_bb, "chrM.200")
# plot_qc_dot(umap_coords, d_bb, "chrM.204")
# plot_qc_dot(umap_coords, d_bb, "chrM.310")
plot_qc_dot(umap_coords, d_bb, "chrM.4330")
plot_qc_dot(umap_coords, d_c, "chrM.4330")


# setdiff(l_plot[[3]], c(l_plot[[1]], l_plot[[2]]))
# plot_qc_dot(umap_coords, d_bb, "chrM.12731")


