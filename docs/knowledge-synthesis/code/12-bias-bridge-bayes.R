# R code for chapter: 12-bias-bridge-bayes.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch12
#| include: false
# Shared setup, sourced at the top of every chapter --------------------------
knitr::opts_chunk$set(
  echo      = TRUE,      # code is part of the text: show it
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  fig.retina = 2,
  out.width = "85%"
)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
})
theme_set(theme_minimal(base_size = 13))
options(dplyr.summarise.inform = FALSE)
library(tidyr)
suppressPackageStartupMessages({
  library(metafor)
  library(survival)
  library(purrr)
})


## -----------------------------------------------------------------------------
#| label: avr-data
avr <- tibble::tribble(
  ~study,       ~treat,              ~n_int, ~n_cs, ~fu_yr, ~xover12, ~hr,  ~lo,  ~hi,
  "RECOVERY",   "SAVR",                  73,    72,   6.2,      9,   0.09, 0.01, 0.67,
  "AVATAR",     "SAVR",                  78,    79,   5.2,     12,   0.46, 0.23, 0.90,
  "EVOLVED",    "SAVR 75%, TAVR 25%",   113,   111,   3.5,     31,   0.79, 0.44, 1.43,
  "EARLY TAVR", "TAVR",                 455,   446,   3.8,    208,   0.50, 0.40, 0.63) |>
  mutate(yi = log(hr), sei = (log(hi) - log(lo)) / (2 * 1.96))   # log HR and SE from the 95% CI

avr |>
  transmute(Trial = study, `Early replacement` = treat, `n (early / surveillance)` = paste(n_int, "/", n_cs),
            `Median follow-up (years)` = fu_yr,
            `Surveillance patients replaced in year 1` = sprintf("%d (%.0f%%)", xover12, 100 * xover12 / n_cs),
            `Published HR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", hr, lo, hi)) |>
  kable(align = "llcccc", caption = "Four trials of early valve replacement versus surveillance in asymptomatic severe aortic stenosis.")


## -----------------------------------------------------------------------------
#| label: avr-ma
m_avr <- rma(yi, sei = sei, data = avr, slab = study, method = "REML", test = "knha")


## -----------------------------------------------------------------------------
#| label: fig-avr-forest
#| fig-cap: "Conventional random-effects meta-analysis of the four trials (REML, Hartung-Knapp interval)."
#| fig-width: 8
#| fig-height: 3.4
forest(m_avr, atransf = exp, at = log(c(0.01, 0.1, 0.5, 1, 2)), header = c("Trial", "HR [95% CI]"),
       xlab = "Hazard ratio, early replacement vs surveillance")


## -----------------------------------------------------------------------------
#| label: avr-crossover
xo   <- escalc(measure = "PLO", xi = xover12, ni = n_cs, data = avr, slab = study)  # logit of a proportion
m_xo <- rma(yi, vi, data = xo, method = "REML", subset = study != "EARLY TAVR")      # the three surgical trials
xo_surg <- predict(m_xo, transf = transf.ilogit)


## -----------------------------------------------------------------------------
#| label: avr-ipd
ipd <- readRDS("data/ch12/dat_all_km.rds") |>                # reconstructed IPD, 1427 patients
  mutate(time  = if_else(study == "AVATAR", time / 12, time), # AVATAR digitised in months
         trt   = as.integer(group != "Clinical Surveillance"),
         study = recode(study, EARLY = "EARLY TAVR"))

cox_hr <- function(d) {                                      # log HR and its SE for one trial
  f <- coxph(Surv(time, event) ~ trt, data = d)
  tibble::tibble(yi = unname(coef(f)), sei = sqrt(vcov(f)[1, 1]))
}
yr1 <- ipd |>                                                # events in year 1 only
  mutate(event = if_else(time > 1, 0, event), time = pmin(time, 1)) |>
  filter(study != "RECOVERY") |>                             # 2 events in year 1: no estimate
  group_by(study) |> group_modify(\(d, k) cox_hr(d)) |> ungroup()
lmk <- ipd |>                                                # landmark: from year 1 onwards
  filter(time > 1) |> mutate(time = time - 1) |>
  group_by(study) |> group_modify(\(d, k) cox_hr(d)) |> ungroup()


## -----------------------------------------------------------------------------
#| label: fig-early-tavr-km
#| fig-cap: "EARLY TAVR, reconstructed data: cumulative incidence of the primary outcome. The vertical line marks the one-year landmark."
#| fig-width: 7.5
#| fig-height: 4
km <- survfit(Surv(time, event) ~ trt, data = filter(ipd, study == "EARLY TAVR"))
tibble::tibble(time = km$time, inc = 1 - km$surv,
               arm = rep(c("Surveillance", "Early TAVR"), km$strata)) |>
  ggplot(aes(time, inc, colour = arm)) +
  geom_step(linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Early TAVR" = "#041E42", Surveillance = "#B22222"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Years from randomisation", y = "Cumulative incidence") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: tbl-avr-landmark
#| tbl-cap: "Hazard ratios (95% CI), early replacement vs surveillance, from reconstructed individual patient data."
hr_txt <- \(d) sprintf("%.2f (%.2f to %.2f)", exp(d$yi), exp(d$yi - 1.96 * d$sei), exp(d$yi + 1.96 * d$sei))
tibble::tibble(Trial = avr$study) |>
  left_join(tibble::tibble(Trial = yr1$study, `Year 1` = hr_txt(yr1)), by = "Trial") |>
  left_join(tibble::tibble(Trial = lmk$study, `After year 1 (landmark)` = hr_txt(lmk)), by = "Trial") |>
  mutate(across(everything(), \(x) tidyr::replace_na(x, "too few events"))) |>
  kable(align = "lcc")
m_lmk <- rma(yi, sei = sei, data = lmk, method = "REML", test = "knha")


## -----------------------------------------------------------------------------
#| label: avr-landmark-numbers
#| include: false
et1 <- filter(yr1, study == "EARLY TAVR"); etl <- filter(lmk, study == "EARLY TAVR")


## -----------------------------------------------------------------------------
#| label: conjugate
update_normal <- function(m0, s0, est, se) {
  prec <- 1 / s0^2 + 1 / se^2
  c(mean = (m0 / s0^2 + est / se^2) / prec, sd = sqrt(1 / prec))
}
update_normal(m0 = 0, s0 = 1.5, est = log(0.50), se = 0.20)   # vague prior: the data dominate
update_normal(m0 = 0, s0 = 0.20, est = log(0.50), se = 0.20)  # sceptical prior as precise as the data


## -----------------------------------------------------------------------------
#| label: bayes-grid-fn
bayes_nn <- function(yi, sei, mu_mean = 0, mu_sd = 1.5,
                     tau_prior = \(t) 2 * dt(t, df = 3),      # half-t(3) on tau, scale 1
                     mu_range = c(-4, 3), tau_max = 4, grid = 400, n_draws = 40000) {
  g  <- expand.grid(mu = seq(mu_range[1], mu_range[2], length.out = grid),
                    tau = seq(0, tau_max, length.out = grid))
  lp <- dnorm(g$mu, mu_mean, mu_sd, log = TRUE) + log(tau_prior(g$tau))   # log prior
  for (i in seq_along(yi))                                                # + log likelihood
    lp <- lp + dnorm(yi[i], g$mu, sqrt(sei[i]^2 + g$tau^2), log = TRUE)
  s  <- g[sample(nrow(g), n_draws, replace = TRUE, prob = exp(lp - max(lp))), ]
  dm <- diff(mu_range) / (grid - 1); dt_ <- tau_max / (grid - 1)          # jitter within a cell
  s$mu  <- s$mu + runif(n_draws, -dm / 2, dm / 2)
  s$tau <- abs(s$tau + runif(n_draws, -dt_ / 2, dt_ / 2))
  tibble::tibble(mu = s$mu, tau = s$tau, new = rnorm(n_draws, s$mu, s$tau))
}


## -----------------------------------------------------------------------------
#| label: avr-brms
#| eval: false
# library(brms)
# pri <- c(prior(normal(0, 1.5), class = Intercept),     # mu, log HR
#          prior(student_t(3, 0, 1), class = sd))        # tau; brms truncates sd priors at zero
# fit_pub <- brm(yi | se(sei) ~ 1 + (1 | study), data = avr, prior = pri,
#                iter = 6000, warmup = 1000, control = list(adapt_delta = 0.99), seed = 2026)
# fit_lmk <- update(fit_pub, newdata = lmk, seed = 2026)


## -----------------------------------------------------------------------------
#| label: avr-bayes
set.seed(2026)
post_avr <- bind_rows(
  bayes_nn(avr$yi, avr$sei) |> mutate(analysis = "Published HRs"),
  bayes_nn(lmk$yi, lmk$sei) |> mutate(analysis = "One-year landmark")) |>
  mutate(analysis = factor(analysis, levels = c("Published HRs", "One-year landmark")))
q_txt <- \(x) sprintf("%.2f (%.2f to %.2f)", exp(median(x)), exp(quantile(x, 0.025)), exp(quantile(x, 0.975)))
avr_tab <- post_avr |>
  summarise(.by = analysis, `HR, median (95% CrI)` = q_txt(mu),
            `New trial (95%)` = sprintf("%.2f to %.2f", exp(quantile(new, 0.025)), exp(quantile(new, 0.975))),
            `P(HR < 1)` = sprintf("%.2f", mean(mu < 0)), `P(HR < 0.9)` = sprintf("%.2f", mean(mu < log(0.9))),
            `τ, median` = sprintf("%.2f", median(tau)))
avr_tab |> kable(align = "lccccc")


## -----------------------------------------------------------------------------
#| label: fig-avr-posterior
#| fig-cap: "Posterior distributions of the pooled hazard ratio, from the published hazard ratios and from the one-year landmark analysis."
#| fig-width: 7.5
#| fig-height: 3.6
ggplot(post_avr, aes(exp(mu), fill = analysis)) +
  geom_density(alpha = 0.55, colour = NA, adjust = 1.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2), limits = c(0.08, 3)) +
  scale_fill_manual(values = c("#041E42", "#E69F00"), name = NULL) +
  labs(x = "Pooled hazard ratio (log scale)", y = "Posterior density") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: tbl-avr-compare
#| tbl-cap: "Landmark analysis: frequentist and Bayesian summaries of the same four estimates."
f_z  <- rma(yi, sei = sei, data = lmk, method = "REML")
f_hk <- rma(yi, sei = sei, data = lmk, method = "REML", test = "knha")
pl   <- filter(post_avr, analysis == "One-year landmark")
fmt3 <- \(m, l, u) sprintf("%.2f (%.2f to %.2f)", exp(m), exp(l), exp(u))
tibble::tibble(
  Analysis = c("REML, Wald interval", "REML, Hartung-Knapp interval", "Bayesian, median and 95% credible interval"),
  `τ` = c(sprintf("%.2f (point estimate)", sqrt(f_z$tau2)), sprintf("%.2f (point estimate)", sqrt(f_hk$tau2)),
          sprintf("posterior median %.2f", median(pl$tau))),
  `Pooled HR` = c(fmt3(f_z$b[1], f_z$ci.lb, f_z$ci.ub), fmt3(f_hk$b[1], f_hk$ci.lb, f_hk$ci.ub),
                  fmt3(median(pl$mu), quantile(pl$mu, 0.025), quantile(pl$mu, 0.975)))) |>
  kable(align = "llc")


## -----------------------------------------------------------------------------
#| label: lead-data
lead <- tibble::tribble(        # harmonised regression coefficients and SEs, Sell et al., Table 2
  ~study,                     ~yi,   ~sei,
  "Bauer 2020",              -0.10,  0.87,
  "Crump 2013",              -3.31,  0.63,
  "Dantzer 2020",            -5.71,  1.61,
  "Desrochers-Couture 2018",  0.02,  0.06,
  "Earl 2016",               -5.86,  2.22,
  "Hong 2015",               -3.14,  1.37,
  "Lee 2021",                -3.83,  2.25,
  "Lucchini 2019",           -1.10,  0.84,
  "Lucchini 2012",           -3.48,  1.48,
  "Menezes-Filho 2018",      -3.74,  1.29,
  "Roy 2013",                -4.70,  1.66,
  "Schnaas 2006",             0.17,  0.81)
m_lead_ee <- rma(yi, sei = sei, data = lead, method = "EE")


## -----------------------------------------------------------------------------
#| label: lead-fits
m_lead  <- rma(yi, sei = sei, data = lead, slab = study, method = "REML", test = "knha")
pr_lead <- predict(m_lead)
ci_lead <- confint(m_lead)

set.seed(2026)
post_lead <- bayes_nn(lead$yi, lead$sei, mu_mean = -1, mu_sd = 2,    # the paper's priors
                      tau_prior = \(t) 2 * dnorm(t, 0, 2),          # half-normal(0, 2) on tau
                      mu_range = c(-10, 6), tau_max = 10, grid = 500)


## -----------------------------------------------------------------------------
#| label: tbl-lead-compare
#| tbl-cap: "Blood lead and IQ: frequentist (REML, Hartung-Knapp) and Bayesian (posterior median, 95% credible interval) results."
fmt <- \(m, l, u) sprintf("%.2f (%.2f to %.2f)", m, l, u)
tibble::tibble(
  Quantity = c("Pooled effect μ", "τ²", "95% prediction interval, new study"),
  Frequentist = c(fmt(m_lead$b[1], m_lead$ci.lb, m_lead$ci.ub),
                  fmt(m_lead$tau2, ci_lead$random["tau^2", "ci.lb"], ci_lead$random["tau^2", "ci.ub"]),
                  sprintf("%.2f to %.2f", pr_lead$pi.lb, pr_lead$pi.ub)),
  Bayesian = c(fmt(median(post_lead$mu), quantile(post_lead$mu, 0.025), quantile(post_lead$mu, 0.975)),
               fmt(median(post_lead$tau^2), quantile(post_lead$tau^2, 0.025), quantile(post_lead$tau^2, 0.975)),
               sprintf("%.2f to %.2f", quantile(post_lead$new, 0.025), quantile(post_lead$new, 0.975)))) |>
  kable(align = "lcc")


## -----------------------------------------------------------------------------
#| label: tbl-lead-probs
#| tbl-cap: "Posterior probabilities for questions a regulator might ask."
tibble::tibble(
  Question = c("Is lead harmful at all?  P(μ < 0)",
               "Loss of more than 1 IQ point per log unit?  P(μ < −1)",
               "Loss of more than 2 IQ points?  P(μ < −2)",
               "Loss of more than 3 IQ points?  P(μ < −3)",
               "Would a new study's true effect be harmful?  P(θ_new < 0)",
               "Is heterogeneity substantial?  P(τ > 1)"),
  `Posterior probability` = sprintf("%.2f", c(mean(post_lead$mu < 0), mean(post_lead$mu < -1),
                                              mean(post_lead$mu < -2), mean(post_lead$mu < -3),
                                              mean(post_lead$new < 0), mean(post_lead$tau > 1)))) |>
  kable(align = "lc")


## -----------------------------------------------------------------------------
#| label: smallk-sim
sei_pool <- lead$sei[lead$sei > 0.1]                  # SEs of the real studies, minus the outlier
one_sim <- function(k, mu = -2.4, tau = 1.8) {
  sei <- sample(sei_pool, k, replace = TRUE)
  yi  <- rnorm(k, mu + rnorm(k, 0, tau), sei)         # observed = true + sampling error
  new <- rnorm(1, mu, tau)                            # true effect in a new study
  m   <- tryCatch(rma(yi, sei = sei, method = "REML", test = "knha",
                      control = list(stepadj = 0.5, maxiter = 500)), error = \(e) NULL)
  if (is.null(m)) return(NULL)
  pz <- predict(update(m, test = "z")); pk <- predict(m)
  pb <- quantile(bayes_nn(yi, sei, mu_mean = -1, mu_sd = 2, tau_prior = \(t) 2 * dnorm(t, 0, 2),
                          mu_range = c(-12, 8), tau_max = 8, grid = 150, n_draws = 4000)$new,
                 c(0.025, 0.975))
  tibble::tibble(k = k, tau_zero = m$tau2 < 1e-4,
                 `Frequentist (Wald)` = new > pz$pi.lb & new < pz$pi.ub,
                 `Frequentist (Hartung-Knapp)` = new > pk$pi.lb & new < pk$pi.ub,
                 Bayesian = new > pb[1] & new < pb[2])
}
set.seed(629)
smallk <- map_dfr(rep(c(4, 6, 12), each = 600), one_sim) |>
  summarise(.by = k, tau0 = mean(tau_zero), across(c(starts_with("Freq"), Bayesian), mean))


## -----------------------------------------------------------------------------
#| label: fig-smallk
#| fig-cap: "Coverage of 95% prediction intervals for the true effect in a new study, by number of studies (600 simulated meta-analyses at each size; true μ = −2.4, τ = 1.8)."
#| fig-width: 7.5
#| fig-height: 3.8
smallk |>
  pivot_longer(-c(k, tau0), names_to = "Method", values_to = "coverage") |>
  mutate(Method = factor(Method, levels = c("Frequentist (Wald)", "Frequentist (Hartung-Knapp)", "Bayesian"))) |>
  ggplot(aes(factor(k), coverage, fill = Method)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#0072B2"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "Number of studies", y = "Coverage of 95% prediction interval") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: tbl-lead-prior
#| tbl-cap: "Blood lead and IQ: sensitivity of the posterior to the prior on τ (the prior on μ is essentially flat, N(0, 100²))."
set.seed(1)
tau_sds <- c("Strong shrinkage, half-normal(0, 0.67)" = 0.67, "Half-normal(0, 1)" = 1,
             "Default, half-normal(0, 2)" = 2, "Weak shrinkage, half-normal(0, 6)" = 6)
imap_dfr(tau_sds, \(sd, nm) {
  p <- bayes_nn(lead$yi, lead$sei, mu_mean = 0, mu_sd = 100, tau_prior = \(t) 2 * dnorm(t, 0, sd),
                mu_range = c(-10, 6), tau_max = 10)
  tibble::tibble(`Prior on τ` = nm,
                 `μ, median (95% CrI)` = fmt(median(p$mu), quantile(p$mu, 0.025), quantile(p$mu, 0.975)),
                 `τ, median` = sprintf("%.2f", median(p$tau)),
                 `95% prediction, new study` = sprintf("%.1f to %.1f", quantile(p$new, 0.025), quantile(p$new, 0.975)))
}) |> kable(align = "lccc")

