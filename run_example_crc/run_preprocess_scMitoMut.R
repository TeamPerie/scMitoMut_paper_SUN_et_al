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

crc = read_rds("./data/raw_signac_crc.rds")

colors = c(
    Tcell = "#f6bd60", 
    Basophil = "#9b948d",
    Myeloid_1 = "#f5cac3",
    Myeloid_2 = "#84a59d",
    Epithelial = "#f28482")
DimPlot(crc, reduction = "umap", cols = colors) 

#+ eval=F
# Load mgatkdata
f_h5 <- parse_mgatk("./data/mgatk_output/final/", prefix = "CRC_v12-mtMask_mgatk", h5_file = "./tmp/mut.h5")
f_h5
x <- open_h5_file(f_h5)
# Filter good quality cells
x = subset_cell(x, colnames(crc))
# Fit model
run_model_fit(x, mc.cores = 6)

#+ eval=T
x <- open_h5_file("./tmp/mut.h5")

# Filter and export mutation

## beta-binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bb",
    p_threshold = 0.01, 
    # alt_count_threshold = 3,
    p_adj_method = "fdr"
)
d_bb = export_df(m_count)
write_tsv(d_bb, "./tmp/raw_crc_scmitomut_bb.tsv")

m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bb",
    p_threshold = 0.01, 
    # alt_count_threshold = 3,
    p_adj_method = "fdr"
)
d_bb = export_df(m_count)
write_tsv(d_bb, "./tmp/raw_crc_scmitomut_bb_0.01.tsv")



## binomial-mixture
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bm",
    p_threshold = 0.01, 
    # alt_count_threshold = 3,
    p_adj_method = "fdr"
)
d_bm = export_df(m_count)
write_tsv(d_bm, "./tmp/raw_crc_scmitomut_bm.tsv")

## binomial
m_count = filter_loc(
    mtmutObj = x, 
    min_cell = 5, 
    model = "bi",
    p_threshold = 0.01, 
    # alt_count_threshold = 3,
    p_adj_method = "fdr"
)
d_bi = export_df(m_count)
write_tsv(d_bi, "./tmp/raw_crc_scmitomut_bi.tsv")

