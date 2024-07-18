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


##########
#  UMAP  #
##########

# mtMutObj
# mtmutobj = open_h5_file("./tmp/mut.h5")
# Cell annotation
so = read_rds("./tmp/so_named.rds")
so
x = Idents(so) 
# Access the UMAP coordinates
umap_coords <- so@reductions$umap@cell.embeddings %>% data.table(keep.rownames = TRUE)


## cell type
d_cell_type = data.table(cell_id = names(x), cell_type = x)
FeaturePlot(so, "PLP1")

d_cell_ann = merge(d_cell_type, umap_coords, by.x = "cell_id", by.y = "rn", all = T)
table(d_cell_ann$cell_type)
d_cell_ann = d_cell_ann[cell_type != "Unknown"]
ggplot(d_cell_ann, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.5) +
    theme_classic() + theme(legend.position = "top")
unique(d_cell_ann$cell_type)

d_cell_ann$cell_type = sub("_\\d+$", "", d_cell_ann$cell_type)
d_cell_ann$cell_type = factor(d_cell_ann$cell_type)
ggplot(d_cell_ann, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.5) +
    theme_classic() + theme(legend.position = "top")
unique(d_cell_ann$cell_type)


## Function
## BI, BM, BB
# colors_models = c("BI" = "#0073C2FF", "BM" = "#EFC000FF", "BB" = "#CE5374")
# colors_binary = c(
    # "PLP1+" = "#f28482",
    # "PLP1-" = "#f6bd60")

####################
#  Mutation count  #
####################

## Preapre data
# BB scMitoMut obj
x <- open_h5_file("./tmp/mut.h5")


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

d_bb = fread("tmp/raw_brain_scmitomut_bb.tsv")
# d_bb = d_bb[N > 10]


## Process MQuad data
d_m = fread("tmp/table_mquad_high_confidence.tsv")

## Process Signac data
d_s = fread("tmp/table_signac_high_confidence.tsv")
d_s[, loc := paste0("chrM.", sub("...$", "", variant))]
d_s[, .(loc, n_cells_conf_detected)]


## Comparing loci
l_plot = list(
    scMitoMut = unique(d_bb$loc),
    Signac = c(), # NOTE: No mutation with clone size > 10
    MQ = unique(d_m$loc)
    )

## BB, Signac, MQuad
colors_methods = c("#CE5374", "grey", "#4F9D9D")
## BI, BM, BB
colors_models = c("#0073C2FF", "#EFC000FF", "#CE5374")

ggvenn(
  l_plot, fill_color = as.character(colors_methods),
  stroke_size = 1.5, set_name_size = 5,
  show_outside = "none", show_percent = F
  )

#####################
#  Cell type trees  #
#####################

## Cluster cell clusters
library(ggtree)
loc_i_v = unique(d_bb$loc)

d_plot = merge(d_bb, d_cell_type, by.x = "cell_barcode", by.y = "cell_id", all = T) %>% na.omit
d_plot = d_plot[cell_type != "Unknown"]
d_plot[, cell_type := sub("_\\d+$", "", cell_type)]

d_plot = d_plot[mut_status == T, .N, by = c("loc", "cell_type")]

m = dcast(d_plot, loc ~ cell_type, value.var = "N", fill = 0)
m = m[, -1]
m = t(t(m) / colSums(m))

d_m = dist(t(m))
plot(hclust(d_m), hang = -1)

# Draw tree with ggtree
library(ape)
tree = hclust(d_m)
tree = as.phylo(tree)
ggtree(tree, branch.length = "none") + geom_tiplab(size = 3) + theme_tree2() 

############################
#  Cell type distribution  #
############################
## Cell type distribution
d_mut = d_bb
d = merge(d_mut, d_cell_type, by.x = "cell_barcode", by.y = "cell_id", all = T)
d = d %>% na.omit
d = d[cell_type != "Unknown"]
d$cell_type = sub("_\\d+$", "", d$cell_type)
cell_count = d[, .N, .(cell_type, cell_barcode)][, table(cell_type)]

d_cell_type_x = d[mut_status == T, .N, by = c("loc", "cell_type")]
d_cell_type_order = d_cell_type_x[, .(N_al = sum(N), N_Mono = sum(N[cell_type == "PLP1+"])), by = loc]
d_cell_type_order = d_cell_type_order[order(-N_al, -N_Mono)]
d_cell_type_x$loc %<>% factor(levels = d_cell_type_order$loc)
ggplot(d_cell_type_x, aes(x = loc, y = N, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() + 
    theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(d_cell_type_x, aes(x = loc, y = N, fill = cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() + 
    theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1))
d_plot_ref = data.table(loc = "All cells", cell_type = names(cell_count), N = as.vector(cell_count))
d_plot_ref$cell_type 
ggplot(d_plot_ref) + aes(x = loc, y = N, fill = cell_type) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() + 
    theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1))

# d[loc == loc_i][is.na(cell_type)]
loc_i_v = unique(d_mut$loc)
d_fisher = lapply(loc_i_v, function(loc_i) {
    x = c(
        d[loc == loc_i & mut_status == T & cell_type == "PLP1+", .N],
        d[loc == loc_i & mut_status == F & cell_type == "PLP1+", .N],
        d[loc == loc_i & mut_status == T & cell_type == "PLP1-", .N],
        d[loc == loc_i & mut_status == F & cell_type == "PLP1-", .N]
        ) %>% matrix(ncol = 2, byrow = T)
    x
    p = fisher.test(x)$p.value
    data.table(
        loc = loc_i,
        PLP1pos_mut = x[1, 1],
        PLP1pos_wt = x[1, 2],
        PLP1neg_mut = x[2, 1],
        PLP1neg_wt = x[2, 2],
        p = p
    )
}) %>% rbindlist()

d_fisher$p_adj = p.adjust(d_fisher$p, method = "fdr")

d_fisher$p %<>% round(3)
d_fisher$p_adj %<>% round(3)
d_fisher[p < 0.05]
d_fisher
# write_tsv(d_fisher, "./tmp/table_fisher_test.tsv")

###################################
#  UMAP mutant cell distribution  #
###################################

# umap
d_bb$loc %>% unique
d_plot = merge(umap_coords, d_bb[loc == "chrM.14001"], by.x = 'rn', by.y = 'cell_barcode')
d_plot$mut_status %<>% as.integer
d_plot$mut_status %<>% factor(levels = c(1, 0), labels = c("Mutant", "WT"))
d_plot %<>% na.omit
ggplot(d_plot[order(mut_status, decreasing = T)]) + aes(x = UMAP_1, y = UMAP_2, color = mut_status) + geom_point() +
    scale_color_manual(values = c(Mutant = "red", WT = "grey")) + theme_bw()


