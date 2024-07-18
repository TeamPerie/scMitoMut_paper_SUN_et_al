library(ggplot2)
library(data.table)
library(magrittr)
 
# coverage_file, the *.coverage.txt.gz file in mgatk output
# exp: "./data/mgatk_out/CD34plus_MO_souris_allo/final/lib_id.coverage.txt.gz"
 
draw_depth_plot <- function(coverage_file, save_tsv = "", title = "", subset_cells = NULL) {
 
    coverage = fread(coverage_file)

    if (!is.null(subset)) {
        coverage = coverage[V2 %in% subset]
    }
 
    names(coverage) = c("loci", "cell_id", "depth")
 
    d = coverage[, .(
        "depth 95%" = quantile(depth, 0.95) %>% as.numeric,
        "depth median" = quantile(depth, 0.50) %>% as.numeric,
        "depth 5%" = quantile(depth, 0.05) %>% as.numeric
        ), by = loci]
 
    d
    d_long = melt(d, id.vars="loci", measure.vars=c("depth 95%", "depth median", "depth 5%"), variable.name = "Quantile", value.name = "depth")
 
    g2 = ggplot(d_long) + aes(x = loci, y = depth, color = Quantile) + geom_line() + 
        coord_polar(theta = "x", start = 0, direction = 1, clip = "on") + 
        scale_y_log10() + theme_classic()
        labs(x = "chrM loc", y = "Depth", title = title)
 
    if (save_tsv != "") {
        write_tsv(d_long, save_tsv)
    }
 
    g2
}

coverage_file = "./data/atac_out_bj_mkn45_1pct/final/mgatk_atac.coverage.txt.gz"

g = draw_depth_plot(coverage_file)
g

# ggsave("./tmp/fig_mtDNA_coverage.pdf", width = 12, height = 12, units = "cm")


mgatk_out_dir = "./data/atac_out_bj_mkn45_1pct/final/"
depth_table_file = paste0(mgatk_out_dir, "mgatk_atac.depthTable.txt")
d_plot = fread(depth_table_file)
names(d_plot) = c("cell_barcode", "depth")
d_plot[order(d_plot$depth), rank := 1:.N]

d_plot[depth > 5][, min(rank)]

g = ggplot(d_plot) + aes(x = rank, y = depth) + geom_bar(stat = "identity") +
    theme_classic() + geom_hline(yintercept = 5, linetype = "dashed") +
    labs(x = "Cell rank", y = "Depth") 
g
