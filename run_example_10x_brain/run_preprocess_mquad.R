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
library(ggpubr)
library(ggvenn)
library(readr)


d = fread("./tmp/mquad_output/BIC_params.csv")
d = d[PASS_KP == T & PASS_MINCELLS==T]
loc = d$variant_name %>% sub("SNP", "", .) %>% paste0("chrM.", .)
d$loc = loc

write_tsv(d, "./tmp/table_mquad_high_confidence.tsv")
