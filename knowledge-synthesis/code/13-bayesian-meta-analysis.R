# R code for chapter: 13-bayesian-meta-analysis.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch13
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
  library(purrr)
})


## -----------------------------------------------------------------------------
#| label: fig-same-p
#| fig-cap: "Two trials with P = 0.03. With a flat prior, the shaded area is the probability of a clinically important effect (RR < 0.95)."
#| fig-width: 10
#| fig-height: 3.4
#| out-width: "100%"
z <- qnorm(1 - 0.03 / 2)
trl <- tibble::tibble(trial = c("Trial 1: 25 000 patients, RR 0.97", "Trial 2: 400 patients, RR 0.72"),
                      rr = c(0.97, 0.72)) |>
  mutate(se = abs(log(rr)) / z, p_imp = pnorm((log(0.95) - log(rr)) / se))
curves <- trl |>
  crossing(x = seq(log(0.4), log(1.4), length.out = 600)) |>
  mutate(dens = dnorm(x, log(rr), se)) |>
  group_by(trial) |> mutate(dens = dens / max(dens)) |> ungroup()
ggplot(curves, aes(exp(x), dens)) +
  geom_area(data = \(d) filter(d, exp(x) < 0.95), fill = "#B22222", alpha = 0.25) +
  geom_line(linewidth = 1, colour = "#041E42") +
  geom_vline(xintercept = 0.95, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(data = trl, aes(x = 1.03, y = 0.85, label = sprintf("P = 0.03\nP(RR < 0.95) = %.2f", p_imp)),
            hjust = 0, size = 3.6, colour = "#B22222") +
  facet_wrap(~ trial) +
  scale_x_log10(breaks = c(0.5, 0.7, 0.95, 1.2)) +
  labs(x = "Risk ratio (log scale)", y = "Relative likelihood")


## -----------------------------------------------------------------------------
#| label: fig-absence
#| fig-cap: "Two trials with P ≈ 0.3. The shaded area is the probability of an important benefit (HR < 0.9), with a flat prior."
#| fig-width: 10
#| fig-height: 3.4
#| out-width: "100%"
trl2 <- tibble::tibble(trial = c("Trial A: 900 patients, HR 0.93", "Trial B: 4400 patients, HR 0.97"),
                       hr = c(0.93, 0.97), se = c(0.0675, 0.03)) |>
  mutate(p = 2 * pnorm(-abs(log(hr)) / se), p_imp = pnorm((log(0.9) - log(hr)) / se))
curves2 <- trl2 |>
  crossing(x = seq(log(0.55), log(1.35), length.out = 600)) |>
  mutate(dens = dnorm(x, log(hr), se)) |>
  group_by(trial) |> mutate(dens = dens / max(dens)) |> ungroup()
ggplot(curves2, aes(exp(x), dens)) +
  geom_area(data = \(d) filter(d, exp(x) < 0.9), fill = "#B22222", alpha = 0.25) +
  geom_line(linewidth = 1, colour = "#041E42") +
  geom_vline(xintercept = 0.9, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(data = trl2, aes(x = 0.57, y = 0.85, label = sprintf("P = %.2f\nP(HR < 0.9) = %.3f", p, p_imp)),
            hjust = 0, size = 3.6, colour = "#B22222") +
  facet_wrap(~ trial) +
  scale_x_log10(breaks = c(0.6, 0.8, 0.9, 1, 1.2)) +
  labs(x = "Hazard ratio (log scale)", y = "Relative likelihood")


## -----------------------------------------------------------------------------
#| label: fig-prior-update
#| fig-cap: "Prior, likelihood and posterior for a pooled log RR of −0.10 (RR 0.90, SE 0.05), under a weakly informative and a sceptical prior. Curves scaled to the same height."
#| fig-width: 10
#| fig-height: 3.6
#| out-width: "100%"
y <- -0.10; s <- 0.05
pri <- tibble::tibble(prior = factor(c("Weakly informative N(0, 0.5²)", "Sceptical N(0, 0.05²)"),
                                     levels = c("Weakly informative N(0, 0.5²)", "Sceptical N(0, 0.05²)")),
                      m0 = 0, s0 = c(0.5, 0.05)) |>
  mutate(prec = 1 / s0^2 + 1 / s^2, pm = (m0 / s0^2 + y / s^2) / prec, ps = sqrt(1 / prec),
         p_ben = pnorm(0, pm, ps))
pri |>
  crossing(x = seq(-0.45, 0.25, length.out = 500)) |>
  mutate(Prior = dnorm(x, m0, s0), Likelihood = dnorm(x, y, s), Posterior = dnorm(x, pm, ps)) |>
  pivot_longer(Prior:Posterior, names_to = "curve") |>
  group_by(prior, curve) |> mutate(value = value / max(value)) |> ungroup() |>
  ggplot(aes(exp(x), value, colour = curve, linetype = curve)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_text(data = pri, aes(x = 0.68, y = 1.08, label = sprintf("P(RR < 1) = %.2f", p_ben)),
            inherit.aes = FALSE, hjust = 0, size = 3.6) +
  facet_wrap(~ prior) +
  scale_colour_manual(values = c(Likelihood = "grey50", Prior = "#E69F00", Posterior = "#041E42"), name = NULL) +
  scale_linetype_manual(values = c(Likelihood = "dashed", Prior = "dotted", Posterior = "solid"), name = NULL) +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.2)) +
  labs(x = "Risk ratio (log scale)", y = "Scaled density") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: fig-prior-types
#| fig-cap: "Four kinds of prior for a log risk ratio, shown on the risk ratio scale."
#| fig-width: 10
#| fig-height: 3
#| out-width: "100%"
lv <- c("Weakly informative\nN(0, 0.5²)", "Sceptical\nN(0, 0.1²)", "Enthusiastic\nN(−0.3, 0.3²)",
        "Empirical\nN(−0.18, 0.12²)")
tibble::tibble(type = factor(lv, levels = lv), m = c(0, 0, -0.3, -0.18), s = c(0.5, 0.1, 0.3, 0.12)) |>
  crossing(x = seq(log(0.25), log(2.5), length.out = 400)) |>
  mutate(d = dnorm(x, m, s)) |>
  ggplot(aes(exp(x), d)) +
  geom_line(linewidth = 1.1, colour = "#041E42") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  facet_wrap(~ type, nrow = 1, scales = "free_y") +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  labs(x = "Risk ratio (log scale)", y = "Prior density") +
  theme(axis.text.y = element_blank())


## -----------------------------------------------------------------------------
#| label: colch-data
colch <- tibble::tribble(          # events / patients from the published forest plot
  ~study,           ~ev_col, ~n_col, ~ev_pbo, ~n_pbo,
  "CISR",                 1,     30,       8,     60,
  "COLCHICINE-PCI",      24,    206,      25,    194,
  "COLCOT",             131,   2366,     170,   2379,
  "COPS",                24,    396,      41,    399,
  "Deftereos 2013",      16,    100,      32,     96,
  "LoDoCo-MI",            0,    119,       2,    118,
  "O'Keefe 1992",        41,    130,      25,     67)
dat <- escalc(measure = "RR", ai = ev_col, n1i = n_col, ci = ev_pbo, n2i = n_pbo,
              data = colch, slab = study)                 # 0.5 added to cells of LoDoCo-MI
m_dl <- rma(yi, vi, data = dat, method = "DL")            # as published
m_hk <- rma(yi, vi, data = dat, method = "REML", test = "knha")
pr_dl <- predict(m_dl, transf = exp)
f_rr <- \(m) sprintf("%.2f (%.2f to %.2f)", exp(m$b[1]), exp(m$ci.lb), exp(m$ci.ub))
tibble::tibble(Method = c("DerSimonian-Laird, Wald (as published)", "REML, Hartung-Knapp"),
               `RR (95% CI)` = c(f_rr(m_dl), f_rr(m_hk)),
               `95% prediction interval` = c(sprintf("%.2f to %.2f", pr_dl$pi.lb, pr_dl$pi.ub),
                                             sprintf("%.2f to %.2f", predict(m_hk, transf = exp)$pi.lb,
                                                     predict(m_hk, transf = exp)$pi.ub)),
               P = sprintf("%.4f", c(m_dl$pval, m_hk$pval)),
               `τ²` = sprintf("%.4f", c(m_dl$tau2, m_hk$tau2))) |>
  kable(align = "lcccc", caption = "The 2022 colchicine meta-analysis, reproduced.")


## -----------------------------------------------------------------------------
#| label: bayes-grid-fn13
bayes_nn <- function(yi, sei, mu_mean = 0, mu_sd = 1, tau_prior = \(t) 2 * dnorm(t, 0, 0.5),
                     mu_range = c(-3, 2), tau_max = 2.5, grid = 400, n_draws = 40000) {
  g  <- expand.grid(mu = seq(mu_range[1], mu_range[2], length.out = grid),
                    tau = seq(0, tau_max, length.out = grid))
  lp <- dnorm(g$mu, mu_mean, mu_sd, log = TRUE) + log(tau_prior(g$tau))
  for (i in seq_along(yi)) lp <- lp + dnorm(yi[i], g$mu, sqrt(sei[i]^2 + g$tau^2), log = TRUE)
  s  <- g[sample(nrow(g), n_draws, replace = TRUE, prob = exp(lp - max(lp))), ]
  dm <- diff(mu_range) / (grid - 1); dt_ <- tau_max / (grid - 1)
  s$mu  <- s$mu + runif(n_draws, -dm / 2, dm / 2)
  s$tau <- abs(s$tau + runif(n_draws, -dt_ / 2, dt_ / 2))
  tibble::tibble(mu = s$mu, tau = s$tau, new = rnorm(n_draws, s$mu, s$tau))
}
set.seed(2025)
post_c <- bayes_nn(dat$yi, sqrt(dat$vi))
q3 <- \(x) sprintf("%.2f (%.2f to %.2f)", exp(median(x)), exp(quantile(x, 0.025)), exp(quantile(x, 0.975)))


## -----------------------------------------------------------------------------
#| label: colch-bayesmeta
#| eval: false
# # The same model with bayesmeta, which integrates numerically rather than on a grid
# library(bayesmeta)
# bma <- bayesmeta(dat, mu.prior.mean = 0, mu.prior.sd = 1,
#                  tau.prior = \(tau) dhalfnormal(tau, scale = 0.5))
# bma$summary; bma$pposterior(mu = log(0.85))


## -----------------------------------------------------------------------------
#| label: fig-colch-post
#| fig-cap: "Seven trials in the 2022 meta-analysis: posterior for the mean risk ratio (left; shaded area RR < 0.85) and prior and posterior for τ (right)."
#| fig-width: 10
#| fig-height: 3.6
#| out-width: "100%"
p_mu <- ggplot(post_c, aes(exp(mu))) +
  geom_density(fill = "#56B4E9", alpha = 0.5, colour = "#041E42", adjust = 1.5) +
  geom_vline(xintercept = 0.85, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_log10(breaks = c(0.5, 0.7, 0.85, 1, 1.2), limits = c(0.4, 1.3)) +
  labs(x = "Mean RR, colchicine vs control (log scale)", y = "Posterior density",
       title = sprintf("μ: P(RR < 0.85) = %.2f", mean(post_c$mu < log(0.85))))
p_tau <- ggplot(post_c, aes(tau)) +
  geom_density(fill = "#B22222", alpha = 0.35, colour = "#B22222", adjust = 1.5) +
  stat_function(fun = \(t) 2 * dnorm(t, 0, 0.5), linetype = "dashed", colour = "grey40") +
  coord_cartesian(xlim = c(0, 1.5)) +
  labs(x = "τ (between-trial SD of the log RR)", y = "Density",
       title = "τ: posterior (shaded) and prior (dashed)")
patchwork::wrap_plots(p_mu, p_tau)


## -----------------------------------------------------------------------------
#| label: fig-colch-shrink
#| fig-cap: "Observed risk ratios (open circles, 95% CI) and Bayesian shrunken estimates of each trial's true effect (filled, 95% credible interval). The dashed line is the posterior median of μ."
#| fig-width: 8
#| fig-height: 4
shrink <- map_dfr(seq_len(nrow(dat)), \(i) {
  v <- dat$vi[i]; t2 <- pmax(post_c$tau^2, 1e-8)
  prec  <- 1 / v + 1 / t2
  theta <- rnorm(nrow(post_c), (dat$yi[i] / v + post_c$mu / t2) / prec, sqrt(1 / prec))
  tibble::tibble(study = dat$study[i], type = c("Observed", "Shrunken"),
                 est = c(dat$yi[i], median(theta)),
                 lo  = c(dat$yi[i] - 1.96 * sqrt(v), quantile(theta, 0.025)),
                 hi  = c(dat$yi[i] + 1.96 * sqrt(v), quantile(theta, 0.975)))
}) |>
  mutate(study = factor(study, levels = rev(dat$study)))
ggplot(shrink, aes(exp(est), study, colour = type, shape = type)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_vline(xintercept = exp(median(post_c$mu)), linetype = "dashed", colour = "#041E42") +
  geom_pointrange(aes(xmin = exp(lo), xmax = exp(hi)), position = position_dodge(width = 0.55)) +
  scale_x_log10(breaks = c(0.05, 0.2, 0.5, 1, 2)) +
  scale_colour_manual(values = c(Observed = "grey50", Shrunken = "#041E42"), name = NULL) +
  scale_shape_manual(values = c(Observed = 1, Shrunken = 16), name = NULL) +
  labs(x = "Risk ratio (log scale)", y = NULL) +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: clear-update
y_cl <- log(0.99); s_cl <- (log(1.16) - log(0.85)) / (2 * 1.96)
post_cl <- tibble::tibble(
  Prior = factor(c("Vague", "COLCOT", "Earlier meta-analysis"),
                 levels = c("Vague", "COLCOT", "Earlier meta-analysis")),
  m0 = c(0, log(0.77), log(0.73)),
  s0 = c(2, (log(0.96) - log(0.61)) / 3.92, (log(0.87) - log(0.61)) / 3.92)) |>
  mutate(prec = 1 / s0^2 + 1 / s_cl^2,
         pm = (m0 / s0^2 + y_cl / s_cl^2) / prec, ps = sqrt(1 / prec))
post_cl |>
  transmute(Prior,
            `Prior HR (95%)` = sprintf("%.2f (%.2f to %.2f)", exp(m0), exp(m0 - 1.96 * s0), exp(m0 + 1.96 * s0)),
            `Posterior HR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", exp(pm), exp(pm - 1.96 * ps), exp(pm + 1.96 * ps)),
            `P(HR < 1)` = sprintf("%.2f", pnorm(0, pm, ps)),
            `P(HR < 0.90)` = sprintf("%.2f", pnorm(log(0.90), pm, ps)),
            `P(HR < 0.85)` = sprintf("%.2f", pnorm(log(0.85), pm, ps)),
            `P(HR < 0.80)` = sprintf("%.2f", pnorm(log(0.80), pm, ps))) |>
  kable(align = "lcccccc", caption = "Posterior after CLEAR under three priors (normal approximation).")


## -----------------------------------------------------------------------------
#| label: fig-clear-post
#| fig-cap: "Prior (dotted) and posterior (solid) distributions of the hazard ratio after CLEAR, under three priors. The dashed line is the CLEAR likelihood."
#| fig-width: 10
#| fig-height: 3.4
#| out-width: "100%"
post_cl |>
  crossing(x = seq(log(0.5), log(1.4), length.out = 500)) |>
  mutate(Prior_d = dnorm(x, m0, s0), Posterior = dnorm(x, pm, ps), Likelihood = dnorm(x, y_cl, s_cl)) |>
  pivot_longer(c(Prior_d, Posterior, Likelihood), names_to = "curve") |>
  mutate(curve = recode(curve, Prior_d = "Prior")) |>
  ggplot(aes(exp(x), value, colour = curve, linetype = curve)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 1, colour = "grey60") +
  facet_wrap(~ Prior) +
  scale_colour_manual(values = c(Likelihood = "grey50", Prior = "#E69F00", Posterior = "#041E42"), name = NULL) +
  scale_linetype_manual(values = c(Likelihood = "dashed", Prior = "dotted", Posterior = "solid"), name = NULL) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.2)) +
  labs(x = "Hazard ratio (log scale)", y = "Density") +
  theme(legend.position = "top", axis.text.y = element_blank())


## -----------------------------------------------------------------------------
#| label: clear-brms
#| eval: false
# library(brms)
# arms <- tibble::tibble(study = "CLEAR", treat = c(1, 0), r = c(322, 327), n = c(3528, 3534))
# priors <- list(vague  = prior(normal(0, 2), class = b),          # log OR for treat
#                colcot = prior(normal(-0.26, 0.12), class = b),
#                meta   = prior(normal(-0.31, 0.09), class = b))
# fit_vague <- brm(r | trials(n) ~ 1 + treat, family = binomial(), data = arms,
#                  prior = priors$vague, seed = 2025)
# fit_colcot <- update(fit_vague, prior = priors$colcot)          # Stan recompiles for a new prior
# fit_meta   <- update(fit_vague, prior = priors$meta)


## -----------------------------------------------------------------------------
#| label: colch-eight
dat8 <- escalc(measure = "RR", ai = ev_col, n1i = n_col, ci = ev_pbo, n2i = n_pbo, slab = study,
               data = bind_rows(colch, tibble::tibble(study = "CLEAR", ev_col = 322, n_col = 3528,
                                                      ev_pbo = 327, n_pbo = 3534)))
set.seed(2025)
post_c8 <- bayes_nn(dat8$yi, sqrt(dat8$vi))
pr8 <- predict(rma(yi, vi, data = dat8, method = "REML"), transf = exp)
pr7 <- predict(rma(yi, vi, data = dat, method = "REML"), transf = exp)
pi_txt <- \(x) sprintf("%.2f to %.2f", exp(quantile(x, 0.025)), exp(quantile(x, 0.975)))
tibble::tibble(Evidence = c("2022 meta-analysis, 7 trials", "7 trials + CLEAR"),
               `Mean RR, Bayesian (95% CrI)` = c(q3(post_c$mu), q3(post_c8$mu)),
               `τ, posterior median` = sprintf("%.2f", c(median(post_c$tau), median(post_c8$tau))),
               `New trial, frequentist (REML)` = c(sprintf("%.2f to %.2f", pr7$pi.lb, pr7$pi.ub),
                                                   sprintf("%.2f to %.2f", pr8$pi.lb, pr8$pi.ub)),
               `New trial, Bayesian` = c(pi_txt(post_c$new), pi_txt(post_c8$new))) |>
  kable(align = "lcccc", caption = "Before and after adding CLEAR as an eighth trial.")


## -----------------------------------------------------------------------------
#| label: fig-decision
#| fig-cap: "Posterior probability after CLEAR of little or no benefit (HR ≥ 0.95), a grey zone (0.85 to 0.95) and important benefit (HR < 0.85), under three priors."
#| fig-width: 8
#| fig-height: 3.4
post_cl |>
  mutate(`Important benefit (HR < 0.85)` = pnorm(log(0.85), pm, ps),
         `Grey zone (0.85 to 0.95)` = pnorm(log(0.95), pm, ps) - pnorm(log(0.85), pm, ps),
         `Little or none (HR ≥ 0.95)` = 1 - pnorm(log(0.95), pm, ps)) |>
  pivot_longer(c(starts_with("Important"), starts_with("Grey"), starts_with("Little")),
               names_to = "Region", values_to = "p") |>
  mutate(Region = factor(Region, levels = c("Little or none (HR ≥ 0.95)", "Grey zone (0.85 to 0.95)",
                                            "Important benefit (HR < 0.85)"))) |>
  ggplot(aes(p, Prior, fill = Region)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(p > 0.06, sprintf("%.0f%%", 100 * p), "")),
            position = position_stack(vjust = 0.5), colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("#B22222", "#E69F00", "#2E7D32"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Posterior probability after CLEAR", y = NULL) +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: tbl-colch-tau
#| tbl-cap: "The seven trials of the 2022 meta-analysis: sensitivity to the prior on τ."
set.seed(7)
tau_priors <- list(`Half-normal(0.1)` = \(t) 2 * dnorm(t, 0, 0.1),
                   `Half-normal(0.5)` = \(t) 2 * dnorm(t, 0, 0.5),
                   `Half-normal(1)`   = \(t) 2 * dnorm(t, 0, 1),
                   `Half-Cauchy(0.5)` = \(t) 2 * dcauchy(t, 0, 0.5))
imap_dfr(tau_priors, \(pr, nm) {
  p <- bayes_nn(dat$yi, sqrt(dat$vi), tau_prior = pr)
  tibble::tibble(`Prior on τ` = nm, `RR (95% CrI)` = q3(p$mu),
                 `New trial (95%)` = sprintf("%.2f to %.2f", exp(quantile(p$new, 0.025)), exp(quantile(p$new, 0.975))),
                 `P(RR < 0.85)` = sprintf("%.2f", mean(p$mu < log(0.85))))
}) |> kable(align = "lccc")

