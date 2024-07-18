library(devtools)
library(microbenchmark)
library(ggplot2)
library(data.table)
library(readr)

library(ggfigdone)

# fo = fd_load("../ggfigdone_db")

load_all("~/projects/scMitoMut/scMitoMut-devel/")

#+ eval=F
## Simulate 
simu_bb = function(p, t, obs_n, mean, sd) {
    n = round(rlnorm(obs_n, mean, sd))
    (s1 = p * t)
    (s2 = t - p * t)
    x = sapply(n, function(x) {
        VGAM::rbetabinom.ab(n = 1, size = x, shape1 = s1, shape2 = s2)
    })
    list(x = x, n = n)
}

## Fit by VGAM
fit_vgam = function(x, n) {
    fit = VGAM::vglm(cbind(x, n - x) ~ 1, VGAM::betabinomialff, trace = T)
    y = VGAM::Coef(fit)
    c(y[[1]] / (y[[1]] + y[[2]]), y[[1]] + y[[2]])
}

fit_bbmle_mle2 = function(x, n) {
    prob_start = 0.50; theta_start = 10; counts = x; total = n
    mtmp <- function(prob, size, theta) {
        -sum(emdbook::dbetabinom(counts, prob, size, theta, log = TRUE))
    }
    m0 = bbmle::mle2(
        mtmp, 
        start = list(prob = prob_start, theta = theta_start), 
        data = list(size = total)
    )
    bbmle::coef(m0)
}

fit_scMitoMut = function(x, n) {
    y = mle_bb(x, n)
    c(y$a / (y$a + y$b), y$a + y$b)
}


#######################################################################
#                           Benchmark speed                           #
#######################################################################

set.seed(2023)
simu_res = simu_bb(p = 0.5, t = 20, obs_n = 5000, mean = 2, sd = 1)
mbm <- microbenchmark(
    "bbmle_mle2" = { 
        b = fit_bbmle_mle2(simu_res$x, simu_res$n)
    },
    "VGAM" = {
        b = fit_vgam(simu_res$x, simu_res$n)
    },
    "scMitoMut" = {
        b = fit_scMitoMut(simu_res$x, simu_res$n)
    }
)

write_rds(mbm, "./tmp/fit_bb.rds")


#######################################################################
#                         Benchmark accurate                          #
#######################################################################

set.seed(2023)
simu_pars = data.table(
    p = c(0.1, 0.5, 0.9, 0.99),
    t = c(20, 40, 80, 160),
    obs_n = c(5000, 5000, 5000, 5000),
    mean = c(2, 2, 2, 2),
    sd = c(1, 1, 1, 1)
)
simu_pars = simu_pars[rep(1:nrow(simu_pars), 10), ]

res = parallel::mclapply(1:nrow(simu_pars), function(i) {
    simu_res = simu_bb(
        p = simu_pars$p[i], 
        t = simu_pars$t[i], 
        obs_n = simu_pars$obs_n[i], 
        mean = simu_pars$mean[i], 
        sd = simu_pars$sd[i]
    )
    vgam = try(fit_vgam(simu_res$x, simu_res$n))
    bbmle_mle2 = try(fit_bbmle_mle2(simu_res$x, simu_res$n))
    scMitoMut = try(fit_scMitoMut(simu_res$x, simu_res$n))

    data.table(
        fit_method = c("VGAM", "bbmle_mle2", "scMitoMut"),
        p = c(vgam[1], bbmle_mle2[1], scMitoMut[1]),
        t = c(vgam[2], bbmle_mle2[2], scMitoMut[2]),
        p_true = simu_pars$p[i],
        t_true = simu_pars$t[i],
        obs_n = simu_pars$obs_n[i]
        )
}, mc.cores = 6)

d_plot = res %>% rbindlist %>% na.omit
d_plot$p %<>% as.numeric
d_plot$t %<>% as.numeric
d_plot$fit_method = factor(d_plot$fit_method, levels = c("VGAM", "bbmle_mle2", "scMitoMut"))

d_plot2 = d_plot[, .(t_se = sd(t), t_mean = mean(t), p_se = sd(p), p_mean = mean(p)),
    by = .(fit_method, t_true, p_true)]

write_tsv(d_plot2, "./tmp/fit_bb_accurate.tsv")

#+ eval=T
d_plot2 = fread("./tmp/fit_bb_accurate.tsv")
d_plot2$fit_method = factor(d_plot2$fit_method, levels = c("VGAM", "bbmle_mle2", "scMitoMut"))

g = ggplot(d_plot2, aes(x = p_true, y = p_mean)) +
    # geom_point(aes(x = p_true, y = t_true), col = "red", size = 3) + 
    geom_errorbar(aes(ymin = p_mean - p_se, ymax = p_mean + p_se), width = 0.1) +
    geom_point(size = 3, alpha = 0.5) + 
    geom_abline(slope = 1, col = 'red') +
    facet_wrap(~fit_method) +
    # geom_abline(intercept = d_plot$t_true, slope = d_plot$p_true, color = "red") + 
    theme_bw() + theme(legend.position = "top") +
    labs(x = "True p", y = "Estimated p")
g

# fd_add(g, "bb_p", fo)

g = ggplot(d_plot2, aes(x = t_true, y = t_mean)) +
    # geom_point(aes(x = p_true, y = t_true), col = "red", size = 3) + 
    geom_errorbar(aes(ymin = t_mean - t_se, ymax = t_mean + t_se), width = 0.1) +
    geom_point(size = 3, alpha = 0.5) + 
    geom_abline(slope = 1, col = 'red') +
    facet_wrap(~fit_method) +
    # geom_abline(intercept = d_plot$t_true, slope = d_plot$p_true, color = "red") + 
    theme_bw() + theme(legend.position = "top") +
    labs(x = "True t", y = "Estimated t")

g

# fd_add(g, "bb_t", fo)

# fd_save(fo)

sessionInfo()
