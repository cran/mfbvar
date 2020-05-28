## ----preliminaries, echo=FALSE--------------------------------------
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
run_mod <- FALSE

## ----setup, include=FALSE, cache=FALSE------------------------------
library(knitr)
render_sweave()
local({
  hook_error = knit_hooks$get('warning')
  knit_hooks$set(warning = function(x, options) {
    x <- gsub("Warning:", "Warning:\n  ", x)
    hook_error(x, options)
  })
})
set.seed(100)

## ----message = FALSE------------------------------------------------
library("dplyr")
library("ggplot2")
library("alfred")

variables <- c("CPIAUCSL", "UNRATE", "GDPC1")
out <- lapply(variables, get_alfred_series,
           observation_start = "1980-01-01",
           observation_end = "2018-11-01",
           realtime_start = "2018-12-10",
           realtime_end = "2018-12-10")


## -------------------------------------------------------------------
alfred_to_ts <- function(x, freq) {
  ts(x[, 3],
     start = c(1980, 1),
     frequency = freq)
}

mf_list <- mapply(alfred_to_ts, x = out, freq = c(12, 12, 4))
names(mf_list) <- variables

## -------------------------------------------------------------------
log_diff <- function(x) {
  freq <- frequency(x)
  100 * freq * diff(log(x))
}

mf_list[c("CPIAUCSL", "GDPC1")] <- 
  lapply(mf_list[c("CPIAUCSL", "GDPC1")], log_diff)

## -------------------------------------------------------------------
mf_list <- mapply(window, x = mf_list, 
                  start = list(c(1980, 4), c(1980, 4), c(1980, 2)))

## -------------------------------------------------------------------
str(mf_list, vec.len = 2)

## -------------------------------------------------------------------
library("mfbvar")
prior <- set_prior(Y = mf_list, n_lags = 4, n_reps = 1000)

## -------------------------------------------------------------------
prior

## -------------------------------------------------------------------
prior_intervals <- matrix(c(1, 3,
                            4, 8,
                            1, 3), ncol = 2, byrow = TRUE)
moments <- interval_to_moments(prior_intervals)
prior <- update_prior(prior,
                      d = "intercept",
                      prior_psi_mean = moments$prior_psi_mean,
                      prior_psi_Omega = moments$prior_psi_Omega)

## ----ss_plot, fig.asp=0.5, fig.cap = "Prior steady-state intervals"----
plot(prior)

## -------------------------------------------------------------------
prior <- update_prior(prior, n_fcst = 24)

## -------------------------------------------------------------------
summary(prior)

## ----comp_chunk,include=FALSE---------------------------------------
if (run_mod) {
  mod_ss_iw  <- estimate_mfbvar(prior, prior = "ss", variance = "iw")
  mod_ssng_iw <- estimate_mfbvar(prior, prior = "ssng", variance = "iw")
  mod_ss_csv  <- estimate_mfbvar(prior, prior = "ss", variance = "csv")
  mod_ss_fsv <- estimate_mfbvar(prior, prior = "ss", variance = "fsv",
                              n_fac = 1)

  predict_example <- predict(mod_ss_iw, pred_bands = 0.8)
  p1 <- plot(mod_ss_iw, plot_start = "2010-01-01", nrow_facet = 3)
  p2 <- plot(mod_ssng_iw, plot_start = "2010-01-01", nrow_facet = 3)
  pred_df <- bind_rows("Inverse Wishart" = predict(mod_ss_iw, pred_bands = NULL),
                       "Common stochastic volatility" = predict(mod_ss_csv, pred_bands = NULL),
                       "Factor stochastic volatility" = predict(mod_ss_fsv, pred_bands = NULL),
                       .id = "Variance") %>%
    filter(variable == "GDPC1")
  p3 <- ggplot(pred_df, aes(y = factor(fcst_date), x = fcst, fill = Variance)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE,
                                quantiles = 2, alpha = 0.5) +
  labs(x = "US GDP Growth",
       y = "Date of Forecast") +
  coord_cartesian(xlim = c(-5, 10)) +
  theme_minimal() +
  scale_fill_brewer(palette = "YlGnBu")

  const_vol <- median(sqrt(mod_ss_iw$Sigma[3, 3, ]))

  p4 <- varplot(mod_ss_fsv, variables = "GDPC1") +
    geom_hline(yintercept = const_vol ,
               color = "red", linetype = "dashed") +
    coord_cartesian(ylim = c(0, 20))

  p5 <- varplot(mod_ss_csv, variables = "GDPC1") +
    geom_hline(yintercept = const_vol,
               color = "red", linetype = "dashed") +
    coord_cartesian(ylim = c(0, 20))

  mdd_example <- mdd(mod_ss_iw)

  library("parallel")
  par_fun <- function(lambda1, prior) {
    set.seed(2019)
    mod_par <- estimate_mfbvar(prior, prior = "ss", variance = "iw",
                               lambda1 = lambda1, lambda3 = 1)
    mdd(mod_par)
  }

  cl <- makeCluster(2)
  clusterEvalQ(cl, library("mfbvar"))
  lambda1_seq <- seq(0.05, 1, by = 0.05)
  result <- parSapply(cl, lambda1_seq,
                      par_fun, prior = prior)
  stopCluster(cl)
  max_val <- tibble(lambda1 = lambda1_seq,
         mdd = result) %>%
    filter(mdd == max(mdd)) %>% .$lambda1


  save(predict_example, mdd_example, lambda1_seq, result, max_val,
       file = "vignettes/vignette_data.RData",
       compress = "xz")
  ggsave("vignettes/figures/ss_plots-1.pdf", p1, "pdf", width = 5.5, height = 5.5*1.5, units = "in")
  ggsave("vignettes/figures/ss_plots-2.pdf", p2, "pdf", width = 5.5, height = 5.5*1.5, units = "in")
  ggsave("vignettes/figures/ridges-1.pdf", p3, "pdf", width = 7, height = 7*0.65, units = "in")
  ggsave("vignettes/figures/varplot-1.pdf", p4, "pdf", width = 5, height = 5*0.5, units = "in")
  ggsave("vignettes/figures/varplot-2.pdf", p5, "pdf", width = 5, height = 5*0.5, units = "in")
} else {
  load("vignette_data.RData")
}

## ----mod_chunk,eval=FALSE-------------------------------------------
#  mod_ss_iw  <- estimate_mfbvar(prior, prior = "ss", variance = "iw")
#  mod_ssng_iw <- estimate_mfbvar(prior, prior = "ssng", variance = "iw")

## ----mod_chunk2,eval=FALSE------------------------------------------
#  mod_ss_csv  <- estimate_mfbvar(prior, prior = "ss", variance = "csv")
#  mod_ss_fsv <- estimate_mfbvar(prior, prior = "ss", variance = "fsv",
#                                n_fac = 1)

## ----eval=FALSE-----------------------------------------------------
#  predict(mod_ss_iw, pred_bands = 0.8)

## ----echo=FALSE-----------------------------------------------------
predict_example

## ----eval=FALSE-----------------------------------------------------
#  plot(mod_ss_iw, plot_start = "2010-01-01", nrow_facet = 3)
#  plot(mod_ssng_iw, plot_start = "2010-01-01", nrow_facet = 3)

## ----eval=FALSE-----------------------------------------------------
#  pred_iw <- predict(mod_ss_iw, pred_bands = NULL)
#  pred_csv <- predict(mod_ss_csv, pred_bands = NULL)
#  pred_fsv <- predict(mod_ss_fsv, pred_bands = NULL)
#  pred_df <- bind_rows("Inverse Wishart" = pred_iw,
#                       "Common stochastic volatility" = pred_csv,
#                       "Factor stochastic volatility" = pred_fsv,
#                       .id = "Variance") %>%
#    filter(variable == "GDPC1")
#  ggplot(pred_df, aes(y = factor(fcst_date), x = fcst, fill = Variance)) +
#    ggridges::stat_density_ridges(quantile_lines = TRUE,
#                                  quantiles = 2, alpha = 0.5) +
#    labs(x = "US GDP Growth",
#         y = "Date of Forecast") +
#    coord_cartesian(xlim = c(-5, 10)) +
#    theme_minimal() +
#    scale_fill_brewer(palette = "YlGnBu")

## ----eval=FALSE-----------------------------------------------------
#  const_vol <- median(sqrt(mod_ss_iw$Sigma[3, 3, ]))
#  
#  varplot(mod_ss_fsv, variables = "GDPC1") +
#    geom_hline(yintercept = const_vol ,
#               color = "red", linetype = "dashed") +
#    coord_cartesian(ylim = c(0, 20))
#  
#  varplot(mod_ss_csv, variables = "GDPC1") +
#    geom_hline(yintercept = const_vol,
#               color = "red", linetype = "dashed") +
#    coord_cartesian(ylim = c(0, 20))

## ----eval=FALSE-----------------------------------------------------
#  mdd(mod_ss_iw)

## ----echo=FALSE-----------------------------------------------------
mdd_example

## ----par_chunk2, eval = FALSE---------------------------------------
#  library("parallel")
#  par_fun <- function(lambda1, prior) {
#    set.seed(2019)
#    mod_par <- estimate_mfbvar(prior, prior = "ss", variance = "iw",
#                               lambda1 = lambda1, lambda3 = 1)
#    mdd(mod_par)
#  }
#  
#  cl <- makeCluster(4)
#  clusterEvalQ(cl, library("mfbvar"))
#  lambda1_seq <- seq(0.05, 1, by = 0.05)
#  result <- parSapply(cl, lambda1_seq,
#                      par_fun, prior = prior)
#  stopCluster(cl)

## ----mdd_plot,fig.asp = 0.4, fig.cap = sprintf("Logarithm of the marginal data density as a function of $\\lambda_1$ with $\\lambda_3=1$. The point shows the maximum point at $\\lambda_1=%4.1f$.", max_val)----
plot_df <- tibble(lambda1 = lambda1_seq,
       mdd = result)
ggplot(plot_df, aes(x = lambda1, y = mdd)) +
  geom_line() +
  geom_point(data = filter(plot_df, mdd == max(mdd))) +
  labs(y = "Marginal data density (log)",
       x = bquote(lambda[1])) +
  theme_minimal()

