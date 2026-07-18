## ---------------------------------------------------------------------------
## CHIP-BCIS3: Bayesian re-analysis of all-cause mortality (HR 1.54, 95% CI 0.99-2.41)
##
## Four priors on log(HR), all fit as single-trial Bayesian re-analyses:
##   1) Vague:                  log(HR) ~ Normal(0, 2)                 [genuinely flat]
##   2) Empirical Bayes (2025): van Zwet, Wiecek & Gelman (2025) "Meta-analysis
##      with a single study" -- Student-t/lognormal two-level model fit on
##      1,635 Cochrane meta-analyses, applied AS-IS to the log(HR) scale
##   3) Empirical Bayes (2025, rescaled): same, rescaled from the probit/SMD
##      scale on which it was estimated to an approximate logit/log-HR scale
##   4) Empirical Bayes (2021 mixture): van Zwet, Schwab & Senn (2021) --
##      4-component zero-mean normal mixture on the z/SNR scale. This is the
##      prior your brophyj.com Shiny app actually implements as "Empirical
##      Bayes" (see R/eb_mixture.R in the shiny_light repo); it is NOT the
##      same prior as (2) and (3) above.
##
## Also included: a stand-alone validation check that reproduces your NEJM
## letter's reported "16% posterior probability of improved survival," using
## your app's actual default "Weak" prior on the RISK DIFFERENCE scale
## (Normal(0, pp SD backed out of a +/-6 percentage-point 95% interval)).
## This prior is NOT flat: its SD (~3.1pp) is tighter than the trial's own
## SE(RD) (~5.1pp), so it meaningfully shrinks the ~10pp observed excess
## mortality toward zero. That shrinkage -- not a scale difference between
## RD and HR -- is why it gives ~16% versus the ~3-5% produced by the
## genuinely flat log(HR) ~ Normal(0,2) prior in models 1-4 above.
##
## Models 2 and 3 reproduce, in brms/Stan, the same two-level model used by
## van Zwet et al. (2025) via their `baggr` package:
##     b_j       ~ Normal(beta_j, se_j)          [se_j known/fixed]
##     beta_j    =  mu + u_j,   u_j ~ Normal(0, tau)
##     mu        ~  prior on population-level effect
##     tau       ~  prior on between-study heterogeneity (SD)
## Model 4 (2021 mixture) instead shrinks the observed z-value directly via
## an analytic empirical-Bayes mixture posterior (no MCMC needed) -- ported
## directly from your app's R/eb_mixture.R so the numbers match.
##
## Author: J. Brophy, NEJM letter re: CHIP-BCIS3
## ---------------------------------------------------------------------------

## 0. Packages ----------------------------------------------------------------
# install.packages(c("brms", "bayesplot", "ggplot2", "posterior", "dplyr"))
library(brms)
library(bayesplot)
library(ggplot2)
library(posterior)
library(dplyr)

set.seed(20260716)
theme_set(bayesplot::theme_default(base_size = 13))

## 1. Trial data: death from any cause ----------------------------------------
# CHIP-BCIS3, Perera et al. NEJM 2026;394:1779-89, Table 2.
# HR = 1.54, 95% CI 0.99-2.41 (microaxial flow pump vs. standard care)
hr_point <- 1.54
hr_lo    <- 0.99
hr_hi    <- 2.41

log_hr    <- log(hr_point)
log_hr_se <- (log(hr_hi) - log(hr_lo)) / (2 * 1.96)   # back out SE from the 95% CI

cat(sprintf("log(HR) = %.4f, SE = %.4f  (HR = %.2f, 95%% CI %.2f-%.2f)\n",
            log_hr, log_hr_se, hr_point, hr_lo, hr_hi))

trial_dat <- data.frame(
  study  = "CHIP_BCIS3",
  y      = log_hr,
  se     = log_hr_se
)

# Raw event counts (Table 2), used only for the RD-scale validation check.
e_t <- 47; n_t <- 148   # microaxial flow pump
e_c <- 33; n_c <- 152   # standard care

## 2. Empirical Bayes prior parameters, 2025 (van Zwet/Wiecek/Gelman) ---------
# Table 1, "zero mean" row: population effect mu ~ generalized Student-t,
# heterogeneity: log(tau) ~ Normal(mean, sd)  [tau lognormal]
eb_mu_df    <- 5.18
eb_mu_loc   <- 0
eb_mu_scale <- 0.37        # as originally estimated, on probit/SMD scale

eb_tau_mean <- -1.82       # mean of log(tau)
eb_tau_sd   <- 0.90        # sd of log(tau)

# Rescaling to an approximate logit/log-HR scale. The paper's prior was fit on
# the probit scale for binary outcomes (and SMD for continuous outcomes), not
# on log(OR)/log(HR). A commonly used approximate conversion is
# logit ~= k * probit with k in ~1.6-1.8 (k = pi/sqrt(3) ~ 1.81 is one common
# value; we use 1.7 as a round middle estimate). This is an approximation,
# not an exact scale transformation -- shown here alongside the "as-is"
# version so the sensitivity to this assumption is transparent.
k_scale <- 1.7

eb_mu_scale_rescaled  <- eb_mu_scale * k_scale
eb_tau_mean_rescaled  <- eb_tau_mean + log(k_scale)   # log(tau*k) = log(tau)+log(k)
eb_tau_sd_rescaled    <- eb_tau_sd                    # sd unaffected by linear rescaling

## 2b. Empirical Bayes prior, 2021 mixture (van Zwet, Schwab & Senn) ----------
# Ported directly from R/eb_mixture.R in github.com/brophyj/shiny_light, so
# this reproduces exactly what your Shiny app currently computes as its
# "Empirical Bayes" prior. Prior is a 4-component zero-mean normal mixture
# on the signal-to-noise ratio SNR = beta/SE (z-scale); NOT the same model
# as the 2025 t/lognormal prior above.
EB2021_sd_z  <- c(1.17, 1.74, 2.38, 5.73)
EB2021_p_mix <- c(0.32, 0.31, 0.30, 0.07)

eb_posterior_snr <- function(z, sd_z = EB2021_sd_z, p_mix = EB2021_p_mix, n = NULL) {
  p_mix  <- p_mix / sum(p_mix)
  sd_snr <- sqrt(sd_z^2 - 1)
  like   <- dnorm(z, 0, sqrt(sd_snr^2 + 1))
  post_w <- p_mix * like
  post_w <- post_w / sum(post_w)
  Bk      <- sd_snr^2 / (sd_snr^2 + 1)
  post_mu <- Bk * z
  post_sd <- sqrt(Bk)
  sampler <- function(n) {
    comp <- sample.int(length(post_w), n, replace = TRUE, prob = post_w)
    rnorm(n, post_mu[comp], post_sd[comp])
  }
  list(post_w = post_w, post_mu = post_mu, post_sd = post_sd, sampler = sampler)
}

## 3. Model 1: vague prior, log(HR) ~ Normal(0, 2) ----------------------------
# Genuinely (near-)flat prior -- kept as the reference case. Because it is
# so wide relative to the trial's own SE (0.227), the posterior essentially
# reproduces the frequentist result.

prior_vague <- prior(normal(0, 2), class = "Intercept")

fit_vague <- brm(
  bf(y | se(se, sigma = FALSE) ~ 1),
  data   = trial_dat,
  prior  = prior_vague,
  family = gaussian(),
  chains = 4, iter = 4000, warmup = 1000,
  seed   = 20260716,
  control = list(adapt_delta = 0.995),
  refresh = 0
)

## 4. Model 2: empirical Bayes prior, 2025, applied as-is ---------------------
prior_eb_raw <- c(
  prior_string(sprintf("student_t(%f, %f, %f)", eb_mu_df, eb_mu_loc, eb_mu_scale),
               class = "Intercept"),
  prior_string(sprintf("lognormal(%f, %f)", eb_tau_mean, eb_tau_sd),
               class = "sd", group = "study")
)

fit_eb_raw <- brm(
  bf(y | se(se, sigma = FALSE) ~ 1 + (1 | study)),
  data   = trial_dat,
  prior  = prior_eb_raw,
  family = gaussian(),
  chains = 4, iter = 4000, warmup = 1000,
  seed   = 20260716,
  control = list(adapt_delta = 0.995),
  refresh = 0
)

## 5. Model 3: empirical Bayes prior, 2025, rescaled to logit/log-HR scale ----
prior_eb_rescaled <- c(
  prior_string(sprintf("student_t(%f, %f, %f)", eb_mu_df, eb_mu_loc, eb_mu_scale_rescaled),
               class = "Intercept"),
  prior_string(sprintf("lognormal(%f, %f)", eb_tau_mean_rescaled, eb_tau_sd_rescaled),
               class = "sd", group = "study")
)

fit_eb_rescaled <- brm(
  bf(y | se(se, sigma = FALSE) ~ 1 + (1 | study)),
  data   = trial_dat,
  prior  = prior_eb_rescaled,
  family = gaussian(),
  chains = 4, iter = 4000, warmup = 1000,
  seed   = 20260716,
  control = list(adapt_delta = 0.995),
  refresh = 0
)

## 6. Model 4: empirical Bayes prior, 2021 mixture (matches your app) --------
# Analytic -- no MCMC needed. Shrinks the observed z = log(HR)/SE directly;
# there is no separate population- vs. trial-level effect in this model
# (beta == mu here), unlike the hierarchical 2025-based models above.
z_obs        <- log_hr / log_hr_se
eb2021_fit   <- eb_posterior_snr(z_obs)
n_draws_eb21 <- 12000   # matches the post-warmup draw count of the brms fits
snr_draws_eb21 <- eb2021_fit$sampler(n_draws_eb21)
beta_draws_eb21 <- snr_draws_eb21 * log_hr_se   # back to log(HR) scale

## 7. Extract posterior draws (all models, log-HR / HR scale) -----------------
# beta = study-level (trial-specific) log(HR); mu = population-level log(HR).
# Models 1 and 4 have no hierarchy, so beta == mu == Intercept / EB draw.

draws_vague <- as_draws_df(fit_vague) %>%
  transmute(model = "Vague N(0,2)",
            beta  = b_Intercept,
            mu    = b_Intercept)

draws_eb_raw <- as_draws_df(fit_eb_raw) %>%
  transmute(model = "Empirical Bayes 2025 (as-is)",
            beta  = b_Intercept + `r_study[CHIP_BCIS3,Intercept]`,
            mu    = b_Intercept)

draws_eb_rescaled <- as_draws_df(fit_eb_rescaled) %>%
  transmute(model = "Empirical Bayes 2025 (rescaled)",
            beta  = b_Intercept + `r_study[CHIP_BCIS3,Intercept]`,
            mu    = b_Intercept)

draws_eb2021 <- tibble(
  model = "Empirical Bayes 2021 (mixture, matches app)",
  beta  = beta_draws_eb21,
  mu    = beta_draws_eb21
)

all_draws <- bind_rows(draws_vague, draws_eb_raw, draws_eb_rescaled, draws_eb2021) %>%
  mutate(model = factor(model, levels = c(
    "Vague N(0,2)",
    "Empirical Bayes 2025 (as-is)",
    "Empirical Bayes 2025 (rescaled)",
    "Empirical Bayes 2021 (mixture, matches app)"
  )))

## 8. Summary table: posterior HR + probability of benefit --------------------
# Set your own clinically-meaningful HR threshold below (placeholder = 0.80).
hr_benefit_threshold <- 0.80

summary_tbl <- all_draws %>%
  group_by(model) %>%
  summarise(
    beta_HR_median   = exp(median(beta)),
    beta_HR_lo95     = exp(quantile(beta, 0.025)),
    beta_HR_hi95     = exp(quantile(beta, 0.975)),
    mu_HR_median     = exp(median(mu)),
    mu_HR_lo95       = exp(quantile(mu, 0.025)),
    mu_HR_hi95       = exp(quantile(mu, 0.975)),
    P_beta_HR_lt_1   = mean(beta < 0),
    P_beta_HR_lt_thr = mean(beta < log(hr_benefit_threshold)),
    .groups = "drop"
  )

print(summary_tbl, width = Inf)
write.csv(summary_tbl, "chip_bcis3_bayesian_summary.csv", row.names = FALSE)

## 9. Plots ---------------------------------------------------------------------

## 9a. bayesplot: overlaid posterior densities of the trial-level HR
# All brms models are fit with the same chains/iter/warmup, and the 2021
# mixture draws are generated with a matching draw count, so columns align.
hr_draws_wide <- data.frame(
  `Vague N(0,2)`                                 = exp(draws_vague$beta),
  `Empirical Bayes 2025 (as-is)`                 = exp(draws_eb_raw$beta),
  `Empirical Bayes 2025 (rescaled)`              = exp(draws_eb_rescaled$beta),
  `Empirical Bayes 2021 (mixture, matches app)`  = exp(draws_eb2021$beta),
  check.names = FALSE
)

p_bayesplot <- mcmc_areas(
  hr_draws_wide,
  prob = 0.95
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "Posterior distribution of the trial-level hazard ratio",
       subtitle = "Death from any cause, CHIP-BCIS3 (mAFP vs. standard care)",
       x = "Hazard ratio") +
  scale_x_continuous(limits = c(0.3, 4))

ggsave("chip_bcis3_posterior_bayesplot.png", p_bayesplot,
       width = 9, height = 5.5, dpi = 300)

## 9b. ggplot2: shrinkage plot -- MLE vs. each posterior for beta (trial-level)
mle_df <- data.frame(
  model = "MLE (unadjusted)",
  median = exp(log_hr),
  lo95   = exp(log_hr - 1.96 * log_hr_se),
  hi95   = exp(log_hr + 1.96 * log_hr_se)
)

post_df <- summary_tbl %>%
  transmute(model, median = beta_HR_median, lo95 = beta_HR_lo95, hi95 = beta_HR_hi95)

shrinkage_df <- bind_rows(mle_df, post_df) %>%
  mutate(model = factor(model, levels = c(
    "MLE (unadjusted)",
    "Vague N(0,2)",
    "Empirical Bayes 2025 (as-is)",
    "Empirical Bayes 2025 (rescaled)",
    "Empirical Bayes 2021 (mixture, matches app)"
  )))

p_shrinkage <- ggplot(shrinkage_df, aes(x = model, y = median, ymin = lo95, ymax = hi95)) +
  geom_pointrange(size = 0.8, color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Shrinkage of the mortality hazard ratio estimate",
       subtitle = "CHIP-BCIS3, death from any cause",
       x = NULL, y = "Hazard ratio (log scale, 95% interval)") +
  theme_minimal(base_size = 13)

ggsave("chip_bcis3_shrinkage_ggplot.png", p_shrinkage,
       width = 9, height = 5.5, dpi = 300)

## 10. Console summary -----------------------------------------------------------
cat("\n--- Posterior probability HR < 1 (any survival benefit) ---\n")
print(summary_tbl %>% select(model, P_beta_HR_lt_1))

cat(sprintf("\n--- Posterior probability HR < %.2f (pre-specified benefit threshold) ---\n",
            hr_benefit_threshold))
print(summary_tbl %>% select(model, P_beta_HR_lt_thr))

## ---------------------------------------------------------------------------
## 11. VALIDATION CHECK: reproduce the letter's reported "16%" figure
## ---------------------------------------------------------------------------
## This uses your Shiny app's ACTUAL default "Weak" prior on the RISK
## DIFFERENCE scale: centered at 0, with a 95% prior interval of +/-6
## percentage points (i.e. Normal(0, 0.06/1.96) on the proportion scale),
## combined analytically (conjugate normal-normal) with the observed RD
## from the raw event counts. This is NOT flat: prior SD (~3.1pp) is
## tighter than the trial's own SE(RD) (~5.1pp) -- that mismatch, not a
## scale difference between RD and HR, is what produces ~16% rather than
## the ~3-5% seen with the genuinely flat log(HR) ~ Normal(0,2) prior above.

p_t <- e_t / n_t
p_c <- e_c / n_c
rd_obs   <- p_t - p_c
se_rd    <- sqrt(p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c)

weak_rd_prior_mean <- 0
weak_rd_prior_sd   <- 0.06 / 1.96     # +/-6pp at 95% -> SD in proportion units

post_var <- 1 / (1 / se_rd^2 + 1 / weak_rd_prior_sd^2)
post_mean <- post_var * (rd_obs / se_rd^2 + weak_rd_prior_mean / weak_rd_prior_sd^2)
post_sd   <- sqrt(post_var)

p_benefit_rd <- pnorm(0, mean = post_mean, sd = post_sd)

cat("\n--- Validation: reproducing the letter's RD-scale Bayesian reanalysis ---\n")
cat(sprintf("Observed RD = %.2f pp (SE %.2f pp)\n", rd_obs * 100, se_rd * 100))
cat(sprintf("Weak prior:  RD ~ Normal(0, %.2f pp)  [app default, +/-6pp at 95%%]\n",
            weak_rd_prior_sd * 100))
cat(sprintf("Posterior:   RD ~ Normal(%.2f pp, %.2f pp)\n",
            post_mean * 100, post_sd * 100))
cat(sprintf("P(RD < 0, i.e. improved survival) = %.1f%%\n", p_benefit_rd * 100))
