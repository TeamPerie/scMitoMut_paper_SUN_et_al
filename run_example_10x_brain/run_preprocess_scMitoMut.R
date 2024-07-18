library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Signac)
library(Seurat)
library(rhdf5)
library(paletteer) 

library(devtools)
load_all("~/projects/scMitoMut/scMitoMut-devel/")

# pbmc = read_rds("./tmp/raw_mito_pbmc.rds")

# d_c[loc == "chrM.16499" & mut_status == T, cell_barcode %in% cell_id_high_depth]
# cell_annot = fread("./tmp/cell_type.tsv")
# cell_annot[barcode %in% d_c[loc == "chrM.16499" & mut_status == T, cell_barcode]]

#+ eval=F
# Load mgatkdata
f_h5 <- parse_mgatk("./data/mgatk_out/final/", prefix = "mgatk_atac", h5_file = "./tmp/mut.h5")
f_h5
f_h5 = "./tmp/mut.h5"
x <- open_h5_file(f_h5)
# Filter good quality cells
cell_id_high_depth = fread("./data/mgatk_out/final/mgatk_atac.depthTable.txt")[V2 > 5, V1]
# x = subset_cell(x, colnames(pbmc))
x = subset_cell(x, cell_id_high_depth)
# x = subset_loc(x, unique(d_bb$loc))
x = subset_loc(x, x$loc_list)
# Fit model
# run_model_fit(x, mc.cores = 6, bb_over_bm = T, bb_over_bm_p = 0.001, bb_over_bm_adj = "fdr")
run_model_fit(x, mc.cores = 6, bb_over_bm = T)

#+ eval=T
x <- open_h5_file("./tmp/mut.h5")

# Filter and export mutation

## beta-binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 11, 
    model = "bb",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
d_bb = export_dt(m_count, all_cell = T)
# d_bb = d_bb[cell_barcode %in% colnames(so)]
# d_bb[mut_status == T & (coverage - alt_depth) >= 3, .N, by = loc]
d_bb[mut_status == T, .N, by = loc]
write_tsv(d_bb, "./tmp/raw_brain_scmitomut_bb.tsv")

#+ eval=F
##### Test
loc_i_v_new = d_bb$loc %>% unique
loc_i_v_new %>% sort

pbmc = read_rds("./tmp/raw_mito_pbmc.rds")
umap_coords <- pbmc@reductions$umap@cell.embeddings %>% data.table(keep.rownames = TRUE)

pdf("tmp/fig_loci_af_coverage_c_specific.pdf", width = 10, height = 5)
for (loc_x in loc_i_v_new) {
    d_plot = merge(umap_coords, d_bb[loc == loc_x], by.x = "rn", by.y = "cell_barcode", all = F)
    d_plot = na.omit(d_plot)
    d_plot = d_plot[order(mut_status)]
    g_umap = ggplot(d_plot) + aes(x = UMAP_1, y = UMAP_2, color = mut_status) + geom_point() +
        scale_color_manual(values = c("grey", "#CE5374")) + theme_bw() + theme(legend.position = "none")
    g_scatter = plot_af_coverage(m_count, loc_x, p_threshold = 0.01, alt_count_threshold = 4, p_adj_method = "bonferroni", model = "bb")
    g_umap = g_umap + ggtitle(loc_x)
    egg::ggarrange(g_umap, g_scatter, ncol = 2)
}
dev.off()
##### Test

## binomial-mixture
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 10, 
    model = "bm",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
d_bm = export_df(m_count)
write_tsv(d_bm, "./tmp/raw_pbmc_scmitomut_bm.tsv")

## binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 10, 
    model = "bi",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
d_bi = export_df(m_count)
write_tsv(d_bi, "./tmp/raw_pbmc_scmitomut_bi.tsv")


