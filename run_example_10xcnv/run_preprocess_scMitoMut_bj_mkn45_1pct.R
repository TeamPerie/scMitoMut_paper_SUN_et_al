library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(rhdf5)
library(paletteer) 

library(devtools)
load_all("~/projects/scMitoMut/scMitoMut-devel/")

#+ eval=F
dataset_dir = './data/atac_out_bj_mkn45_1pct/'
dataset_name = basename(dataset_dir) %>% str_replace("atac_out_", "")

h5_file = paste0("./tmp/", dataset_name, ".h5")
mgatk_out_dir = paste0(dataset_dir, "/final/")
print(h5_file)
print(mgatk_out_dir)
f_h5 <- parse_mgatk(mgatk_out_dir, prefix = "mgatk_atac", h5_file = h5_file)
depth_table_file = paste0(mgatk_out_dir, "mgatk_atac.depthTable.txt")
cell_id_high_depth = fread(depth_table_file)[V2 > 5, V1]
x <- open_h5_file(h5_file)
x = subset_cell(x, cell_id_high_depth)
x = subset_loc(x, x$loc_list)
run_model_fit(x, mc.cores = 14, bb_over_bm = T)
H5Fclose(x$h5f)

#+ eval=T
bb_file = paste0("./tmp/", dataset_name, "_bb.tsv")
bm_file = paste0("./tmp/", dataset_name, "_bm.tsv")
bi_file = paste0("./tmp/", dataset_name, "_bi.tsv")

x <- open_h5_file(f_h5)

# Filter and export mutation

## beta-binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bb",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
m_count
d_bb = export_dt(m_count, all_cell = T)
write_tsv(d_bb, bb_file)

## binomial-mixture
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bm",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
m_count
d_bm = export_dt(m_count)
d_bm[mut_status == T, .N, by = loc]
write_tsv(d_bm, bm_file)

## binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bi",
    p_threshold = 0.01, 
    p_adj_method = "fdr"
)
d_bi = export_dt(m_count)
write_tsv(d_bi, bi_file)



