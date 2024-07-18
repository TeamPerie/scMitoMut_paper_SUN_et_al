library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(microbenchmark)

library(ggfigdone)
# fo = fd_load("../ggfigdone_db")


mbm_bm = read_rds("./tmp/fit_bm.rds")
mbm_bb = read_rds("./tmp/fit_bb.rds")

autoplot(mbm_bm)
autoplot(mbm_bb)

#######################################################################
#                                 BM                                  #
#######################################################################
d_plot = data.table(mbm_bm)
## normlize the time to speed
d_plot[, speed := 1 / (time/ 1000000)]
## normalize with the mixtools as 1
mbm_mean = d_plot[expr == "mixtools", median(speed)]
d_plot[, speed := speed / mbm_mean]
d_plot[, .(median = median(speed), Q75vQ25 = quantile(speed, 0.75) - quantile(speed, 0.25)), by = expr]

g = ggplot(d_plot) + aes(x = speed, y = expr) + 
    geom_jitter(alpha = 0.5) + 
    scale_x_log10() + theme_classic()
g

# fd_add(g, "bm_speed", fo)

#######################################################################
#                                 BB                                  #
#######################################################################
d_plot = data.table(mbm_bb)
## normlize the time to speed
d_plot[, speed := 1 / (time/ 1000000)]
## normalize with the mixtools as 1
mbm_mean = d_plot[expr == "VGAM", median(speed)]
d_plot[, speed := speed / mbm_mean]

d_plot[, .(median = median(speed), Q75vQ25 = quantile(speed, 0.75) - quantile(speed, 0.25)), by = expr]

g = ggplot(d_plot) + aes(x = speed, y = expr) + 
    geom_jitter(alpha = 0.5) + 
    scale_x_log10() + theme_classic()
g

# fd_add(g, "bb_speed", fo)

# fd_save(fo)
