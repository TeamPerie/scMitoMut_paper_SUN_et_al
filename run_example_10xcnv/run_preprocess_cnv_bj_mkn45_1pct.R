library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

# Load data
d_cluster = fread("./data/cnv/bj_mkn45_1pct_dloupe-group_2110-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv")
v_cluster_id = d_cluster[[1]]
l_cell_list = d_cluster[[2]] %>% str_split(";")
v_cell_count = l_cell_list %>% lapply(length) %>% unlist

v_cluster_id[v_cell_count == 11]  # NOTE cluster 2105 is MKN45 cells

d_cell_ph = data.table(
    cluster_id = rep(v_cluster_id, v_cell_count),
    cell_barcode = unlist(l_cell_list)
)

d_cell_ph[, cell_type := ifelse(cluster_id == 2105, "MKN45", "BJ")]


## NOTE: Cell number of BJ and MKN45 is 1045 and 11, respectively.
length(unlist(l_cell_list))
1056 - 11

## The high quality cells
mgatk_out_dir = "./data/atac_out_bj_mkn45_1pct/final/"
depth_table_file = paste0(mgatk_out_dir, "mgatk_atac.depthTable.txt")
cell_id_high_depth = fread(depth_table_file)[V2 > 5, V1]
length(cell_id_high_depth)  # NOTE 962 cell totally
d_cell_ph[, high_quality := cell_barcode %in% cell_id_high_depth]

d_cell_ph[high_quality == TRUE, .N, by = cell_type] # NOTE 962 cells in total, 11 MKN45 cells

write_tsv(d_cell_ph, "./tmp/table_cell_phenotype.tsv")


# Heatmap of CNV variation
