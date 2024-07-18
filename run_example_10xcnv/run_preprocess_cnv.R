library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(pheatmap)

x = fread("./data/cnv/mkn45_5k_750k_rpc_dloupe-group_10436-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv")
node_id = x[[1]]
cell_id = x[[2]] %>% str_split(";")
cell_rep = rep(1:length(cell_id), sapply(cell_id, length))
x_m = x[, -c(1, 2)]
x_m = as.matrix(x_m)
x_m = x_m[cell_rep, ]
rownames(x_m) = cell_id %>% unlist

x_select = names(x) %>% grep("^4:", ., value = TRUE)
x_m_sub = x_m[, x_select]
x_m_sub = t(scale(t(x_m_sub)))
color = colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(x_m_sub, cluster_cols = F, cluster_rows = T, color = color, scale = 'none', show_rownames = F, show_colnames = F)

x_cluster = hclust(dist(x_m_sub))
plot(x_cluster)
x_cluster_id = cutree(x_cluster, h = 7)
table(x_cluster_id)
d_cluster_cnv = data.table(cell_id = rownames(x_m_sub), cluster_id = x_cluster_id)

write_tsv(d_cluster_cnv, "./tmp/mkn45_5k_750k_cnv_cluster_id.tsv")


baf = fread("./tmp/BAF_bin_phased_adata_obs.csv")
cnv = fread("./tmp/mkn45_5k_750k_cnv_cluster_id.tsv")
sum(baf[sub24_2_r1 == "Cluster 1", index] %in% cnv[cluster_id == 5, cell_id])

node_id_v = node_id[cell_rep]
cell_id_v = unlist(cell_id)

node_id_v[cell_id_v %in% baf[sub24_2_r1 == "Cluster 1", index]] %>% table
node_id_v[cell_id_v %in% cnv[cluster_id == 5, cell_id]] %>% table

table(node_id_v)

m = fread("./tmp/mkn45_5k_750k_rpc_bb.tsv")
node_id_v[cell_id_v %in% m[loc == "chrM.2393" & mut_status == T, cell_barcode]] %>% table
node_id_v[cell_id_v %in% m[loc == "chrM.8368" & mut_status == T, cell_barcode]] %>% table
