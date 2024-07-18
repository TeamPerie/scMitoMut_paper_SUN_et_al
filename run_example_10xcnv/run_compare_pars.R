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
library(ggpubr)
library(parallel)

library(devtools)
load_all("~/projects/scMitoMut/scMitoMut-devel/")

library(ggfigdone)

# fd_load("../ggfigdone_db")


d_cell_ph = fread("./tmp/table_cell_phenotype.tsv")


f_h5 = "./tmp/bj_mkn45_1pct.h5"
x <- open_h5_file(f_h5)

# Filter and export mutation

## beta-binomial
get_benchmark = function(obj = x, model = "bb", p_threshold = 0.01, cell_num = 1) {
    # obj = x
    # model = "bb"
    # p_threshold = 0.01
    # cell_num = 5

    m_count = filter_loc(
        mtmutObj = obj, 
        min_cell = cell_num, 
        model = model,
        p_threshold = p_threshold, 
        p_adj_method = "fdr"
    )

    d_out = export_dt(m_count, all_cell = T)
    d_out = merge(d_cell_ph, d_out, by = "cell_barcode")
    d_summ = d_out[mut_status == T, 
        .(
        lineage_precision = sum(cell_type == "MKN45") / .N,
        clone_size = sum(cell_type == "MKN45")
        ), 
    ## at least 1 cells in MKN45
    by = loc][clone_size >= 1]
}

###############################
#  Compare p value threshold  #
###############################
v_p = c(0.001, 0.005, 0.01, 0.05, 0.1)

#+ eval=F
l_summary = mclapply(v_p, function(p) {
    print(p)
    d_bb = get_benchmark(x, "bb", p, 1)
    d_bb$p_threshold = p
    d_bb
}, mc.cores = 12)
d_summary = rbindlist(l_summary)
write_tsv(d_summary, "./tmp/table_p_threshold.tsv")

#+ eval=T
d_summary = fread("./tmp/table_p_threshold.tsv")
d_summary$p_threshold = factor(d_summary$p_threshold, levels = rev(as.character(v_p)))
d_summary[, .N, by = p_threshold]

ggplot(d_summary[p_threshold %in% c(0.05, 0.01, 0.005), .N, by = p_threshold]) + aes(x = factor(p_threshold, levels = rev(as.character(v_p))), y = N) + geom_bar(stat = "identity") + theme_classic()

comparison = list(
    # c("0.001", "0.005"),
    c("0.005", "0.01"),
    # c("0.001", "0.01"),
    c("0.01", "0.05")
    # c("0.05", "0.1")
    # c("0.01", "0.1")
)

ggplot(d_summary[p_threshold %in% c(0.05, 0.01, 0.005)]) + aes(x = p_threshold, y = clone_size) + 
    geom_jitter(alpha = 0.5, width = 0.2) + stat_compare_means(comparisons = comparison) + ylim(0, 15) +
    theme_classic()

ggplot(d_summary[p_threshold %in% c(0.05, 0.01, 0.005)]) + aes(x = p_threshold, y = lineage_precision) + 
    geom_jitter(alpha = 0.5, width = 0.2) + stat_compare_means(comparisons = comparison) +
    theme_classic() + scale_y_continuous(limits = c(0, 1.2))


#################################
#  Compare Clone size threshold #
#################################

#+ eval=F
v_cell_n = c(1, 5, 9)
l_summary = mclapply(v_cell_n, function(n) {
    d_bb = get_benchmark(x, "bb", 0.01, n)
    d_bb$cell_n_threshold = n
    d_bb
}, mc.cores = 12)
d_summary = rbindlist(l_summary)
write_tsv(d_summary, "./tmp/table_cell_n_threshold.tsv")

#+ eval=T
d_summary = fread("./tmp/table_cell_n_threshold.tsv")
d_summary$cell_n_threshold = factor(d_summary$cell_n_threshold, levels = (as.character(v_cell_n)))
d_summary[, .N, by = cell_n_threshold]

ggplot(d_summary[cell_n_threshold %in% v_cell_n, .N, by = cell_n_threshold]) + aes(x = factor(cell_n_threshold, levels = (as.character(v_cell_n))), y = N) + geom_bar(stat = "identity") + theme_classic()

comparison = list(
    c("1", "5"),
    c("5", "9")
    # c("3", "5"),
    # c("5", "7"),
    # c("7", "9"),
    # c("1", "9")
)

ggplot(d_summary[cell_n_threshold %in% v_cell_n]) + aes(x = cell_n_threshold, y = clone_size) + 
    geom_jitter(alpha = 0.5, width = 0.2) + stat_compare_means(comparisons = comparison) + ylim(0, 15) +
    theme_classic()

ggplot(d_summary[cell_n_threshold %in% v_cell_n]) + aes(x = cell_n_threshold, y = lineage_precision) + 
    geom_jitter(alpha = 0.5, width = 0.2) + stat_compare_means(comparisons = comparison) + ylim(0, 1.2) +
    theme_classic()


