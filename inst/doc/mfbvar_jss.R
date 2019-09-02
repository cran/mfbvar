## ----preliminaries, echo=FALSE--------------------------------------
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

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
library("tidyverse")
library("alfred")

variables <- c("CPIAUCSL", "UNRATE", "GDPC1")
out <- map(variables, get_alfred_series,
           observation_start = "1980-01-01",
           observation_end = "2018-11-01",
           realtime_start = "2018-12-10",
           realtime_end = "2018-12-10")

## -------------------------------------------------------------------
out[[3]]$date <- out[[3]]$date + months(2)

## -------------------------------------------------------------------
log_diff <- function(x, lag = 1) {
  c(rep(NA, lag), 1200/lag * diff(log(x), lag = lag))
}

mf_df <- reduce(out, full_join, by = c("date", "realtime_period")) %>%
  mutate(CPIAUCSL = log_diff(CPIAUCSL),
         GDPC1 = log_diff(GDPC1, lag = 3)) %>%
  filter(date >= "1980-04-01") %>%
  select(-realtime_period)

tail(mf_df)

## -------------------------------------------------------------------
library("mfbvar")
prior <- set_prior(Y = mf_df, freq = c("m", "m", "q"),
                   n_lags = 4, n_burnin = 1000, n_reps = 1000)

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

## ----mod_chunk------------------------------------------------------
mod_ss_iw  <- estimate_mfbvar(prior, prior = "ss", variance = "iw")
mod_ssng_iw <- estimate_mfbvar(prior, prior = "ssng", variance = "iw")


## ----mod_chunk2-----------------------------------------------------
mod_ss_csv  <- estimate_mfbvar(prior, prior = "ss", variance = "csv")
mod_ss_fsv <- estimate_mfbvar(prior, prior = "ss", variance = "fsv",
                              n_fac = 1)


## -------------------------------------------------------------------
predict(mod_ss_iw, pred_bands = 0.8)

## ----ss_plots, fig.cap = "Forecasts and posterior steady-state intervals", fig.subcap= c("Steady-state prior", "Hierarchical steady-state prior"), out.width='0.49\\linewidth', fig.width = 5.5, fig.asp = 1.5----
plot(mod_ss_iw, plot_start = "2010-01-01", nrow_facet = 3)
plot(mod_ssng_iw, plot_start = "2010-01-01", nrow_facet = 3)

## ----ridges, fig.asp = 0.65, message = FALSE, fig.cap = "Distributions of forecasts produced using the steady-state prior with constant or time-varying error covariance. The vertical lines represent the medians."----
pred_iw <- predict(mod_ss_iw, pred_bands = NULL)
pred_csv <- predict(mod_ss_csv, pred_bands = NULL)
pred_fsv <- predict(mod_ss_fsv, pred_bands = NULL)
pred_df <- bind_rows("Inverse Wishart" = pred_iw,
                     "Common stochastic volatility" = pred_csv,
                     "Factor stochastic volatility" = pred_fsv,
                     .id = "Variance") %>%
  filter(variable == "GDPC1")
ggplot(pred_df, aes(y = factor(fcst_date), x = fcst, fill = Variance)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE,
                                quantiles = 2, alpha = 0.5) +
  labs(x = "US GDP Growth",
       y = "Date of Forecast") +
  coord_cartesian(xlim = c(-5, 10)) +
  theme_minimal() +
  scale_fill_brewer(palette = "YlGnBu")

## ----varplot, fig.asp = 0.5, fig.cap = "Standard deviation of the error term in the equation for GDP growth. The black solid and red dashed lines are the medians from the models with factor stochastic and constant volatility, respectively. The bands are obtained from the 95 \\% posterior point-wise intervals.", cache = TRUE, fig.subcap= c("Factor stochastic volatility", "Common stochastic volatility"), out.width='0.49\\linewidth', fig.width = 5----
const_vol <- median(sqrt(mod_ss_iw$Sigma[3, 3, ]))

varplot(mod_ss_fsv, variables = "GDPC1") +
  geom_hline(yintercept = const_vol ,
             color = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 20))

varplot(mod_ss_csv, variables = "GDPC1") +
  geom_hline(yintercept = const_vol,
             color = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 20))

## -------------------------------------------------------------------
mdd(mod_ss_iw)

## ----par_chunk, include = FALSE-------------------------------------
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

## ----par_chunk2, eval = FALSE---------------------------------------
#  
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

## ----include=FALSE--------------------------------------------------
max_val <- tibble(lambda1 = lambda1_seq,
       mdd = result) %>%
  filter(mdd == max(mdd)) %>% .$lambda1

## ----mdd_plot,fig.asp = 0.4, fig.cap = sprintf("Logarithm of the marginal data density as a function of $\\lambda_1$ with $\\lambda_3=1$. The point shows the maximum point at $\\lambda_1=%4.1f$.", max_val)----
plot_df <- tibble(lambda1 = lambda1_seq,
       mdd = result)
ggplot(plot_df, aes(x = lambda1, y = mdd)) +
  geom_line() +
  geom_point(data = filter(plot_df, mdd == max(mdd))) +
  labs(y = "Marginal data density (log)",
       x = bquote(lambda[1])) +
  theme_minimal()

