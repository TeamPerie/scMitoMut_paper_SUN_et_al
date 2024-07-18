library(devtools)
library(microbenchmark)
library(readr)
library(data.table)

library(ggfigdone)

# fo = fd_load("../ggfigdone_db")

load_all("~/projects/scMitoMut/scMitoMut-devel")

#+ eval=F
## simulate data
simu_bm = function(observe_n, mean, sd, theta, p) {
    n = rlnorm(observe_n, mean, sd) %>% round(0) %>% abs
    p = sample(p, observe_n, replace = TRUE, prob = theta) %>% as.numeric()
    x = sapply(seq_along(n), function(i) { rbinom(1, n[i], p[i]) })
    return(list(x = x, n = n))
}

#######################################################################
#                           Benchmark speed                           #
#######################################################################
set.seed(123)
simu_res = simu_bm(observe_n = 1000, mean = 2, sd = 1, 
    theta = c(0.5, 0.5), p = c(0.992, 1))

mbm = microbenchmark(
    "mixtools" = {
        m = data.matrix(data.frame(simu_res$x, simu_res$n-simu_res$x))
        b = mixtools::multmixEM(m)
    },
    "scMitoMut" = {
        ave_p = sum(simu_res$x) / sum(simu_res$n)
        b = em_bm(simu_res$x, simu_res$n, p1 = ave_p, p2 = ave_p / 2, theta1 = 0.95, max_iter = 100, tol = 1e-6)
    }
)


write_rds(mbm, "./tmp/fit_bm.rds")

mbm
autoplot(mbm, log = F)


#######################################################################
#                         Benchmark accurate                          #
#######################################################################

set.seed(2023)

simu_pars = data.table(
    observe_n = 5000,
    mean = 2,
    sd = 1,
    theta1 = c(0.01, 0.1, 0.2, 0.4),
    theta2 = c(0.99, 0.9, 0.8, 0.6),
    p1 = 0.5,
    p2 = 1
    )

simu_pars = simu_pars[rep(1:nrow(simu_pars), 10), ]
# simu_res = simu_bm(observe_n = 1000, mean = 2, sd = 1, 

res = parallel::mclapply(1:nrow(simu_pars), function(i) {
    simu_res = simu_bm(observe_n = simu_pars$observe_n[i], mean = simu_pars$mean[i], 
        sd = simu_pars$sd[i], theta = c(simu_pars$theta1[i], simu_pars$theta2[i]), 
        p = c(simu_pars$p1[i], simu_pars$p2[i]))

    m = data.matrix(data.frame(simu_res$x, simu_res$n-simu_res$x))
    mixture_tools_out <- mixtools::multmixEM(m)

    ave_p = sum(simu_res$x) / sum(simu_res$n)
    smm_out = em_bm(simu_res$x, simu_res$n, p1 = ave_p, p2 = ave_p / 2, theta1 = 0.95, max_iter = 100, tol = 1e-6)

    smm_theta1 = smm_out$theta[1]
    smm_p1 = smm_out$p1
    smm_p2 = smm_out$p2
    if (xor(simu_pars$theta1[i] > 0.5 , smm_theta1 > 0.5)) {
        smm_theta1 = 1 - smm_theta1
        smm_p1 = smm_out$p2
        smm_p2 = smm_out$p1
    }

    mt_theta1 = mixture_tools_out$lambda[1]
    mt_p1 = mixture_tools_out$theta[1]
    mt_p2 = mixture_tools_out$theta[2]
    if (xor(simu_pars$theta1[i] > 0.5,  mt_theta1 > 0.5)) {
        mt_theta1 = 1 - mt_theta1
        mt_p1 = mixture_tools_out$theta[2]
        mt_p2 = mixture_tools_out$theta[1]
    }

    data.table(
        fit_method = c("mixtools", "scMitoMut"),
        theta1 = c(mt_theta1, smm_theta1),
        p1 = c(mt_p1, smm_p1),
        p2 = c(mt_p2, smm_p2),
        theta1_true = simu_pars$theta1[i],
        p1_true = simu_pars$p1[i],
        p2_true = simu_pars$p2[i]
        )
})


d_plot = res %>% rbindlist %>% na.omit
d_plot$fit_method = factor(d_plot$fit_method, levels = c("mixtools", "scMitoMut"))

d_plot2 = d_plot[, .(
    theta1_se   = sd(theta1),
    theta1_mean = mean(theta1),
    p1_se       = sd(p1),
    p1_mean     = mean(p1),
    p2_se       = sd(p2),
    p2_mean     = mean(p2)
    ),
    by = .(fit_method, theta1_true, p1_true, p2_true)
    ]
write_tsv(d_plot2, "./tmp/fit_bm_theta_accurate.tsv")


d_plot3 = d_plot[, .(
    p1_se       = sd(p1),
    p1_mean     = mean(p1),
    p2_se       = sd(p2),
    p2_mean     = mean(p2)
    ),
    by = .(fit_method, p1_true, p2_true)
    ]
write_tsv(d_plot3, "./tmp/fit_bm_p_accurate.tsv")

#+ eval=T
d_plot2 = fread("./tmp/fit_bm_theta_accurate.tsv")
d_plot3 = fread("./tmp/fit_bm_p_accurate.tsv")

g = ggplot(d_plot3, aes(x = fit_method, y = p1_mean)) +
    geom_errorbar(aes(ymin = p1_mean - p1_se, ymax = p1_mean + p1_se), width = 0.1) +
    geom_point(size = 3, alpha = 0.5) + 
    geom_abline(intercept = 0.5, slope = 0, color = "red") +
    theme_bw() + theme(legend.position = "top") +
    labs(x = "True p1", y = "Estimated p1") + ylim(0, 1)
g

# fd_add(g, "bm_p1", fo)

g = ggplot(d_plot3, aes(x = fit_method, y = p2_mean)) +
    geom_errorbar(aes(ymin = p2_mean - p2_se, ymax = p2_mean + p2_se), width = 0.1) +
    geom_point(size = 3, alpha = 0.5) + 
    geom_abline(intercept = 1, slope = 0, color = "red") +
    theme_bw() + theme(legend.position = "top") +
    labs(x = "True p2", y = "Estimated p2") + ylim(0, 1)

g

# fd_add(g, "bm_p2", fo)


g = ggplot(d_plot2, aes(x = theta1_true, y = theta1_mean)) +
    geom_errorbar(aes(ymin = theta1_mean - theta1_se, ymax = theta1_mean + theta1_se), width = 0.1) +
    geom_point(size = 3, alpha = 0.5) + 
    geom_abline(slope = 1, col = 'red') +
    facet_wrap(~fit_method) +
    # geom_abline(intercept = d_plot$t_true, slope = d_plot$p_true, color = "red") + 
    theme_bw() + theme(legend.position = "top") +
    labs(x = "True theta", y = "Estimated theta")

g

# fd_add(g, "bm_theta1", fo)


#######################################################################
#                                mics                                 #
#######################################################################

## mixture tool
simu_res = simu_bm(observe_n = 1000, mean = 2, sd = 1, theta = c(0.01, 0.99), p = c(0.7, 0.9))

time.start = Sys.time()
m = data.matrix(data.frame(simu_res$x, simu_res$n-simu_res$x))
mixture_tools_out <- mixtools::multmixEM(m, k = 2)
Sys.time() - time.start

mixture_tools_out[c("lambda", "theta", "loglik")]

time.start = Sys.time()
smm_out = em_bm(simu_res$x, simu_res$n, p1 = 0.30, p2 = 0.70, theta1 = 0.2, tol = 1e-3)
Sys.time() - time.start
str(smm_out)

smm_out

# fd_save(fo)
sessionInfo()
