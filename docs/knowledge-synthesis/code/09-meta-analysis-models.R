# R code for chapter: 09-meta-analysis-models.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch09
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
ctl <- list(stepadj = 0.5, maxiter = 1000)   # gentler optimiser for iterative tau^2 estimators


## -----------------------------------------------------------------------------
#| label: vote-sim
set.seed(629)
k      <- 33
p0     <- 0.10                                # control risk
rr     <- 0.80                                # true risk ratio
n_arm  <- round(runif(k, 100, 800))           # patients per arm
ev_c   <- rbinom(k, n_arm, p0)
ev_t   <- rbinom(k, n_arm, p0 * rr)

vote <- escalc(measure = "RR", ai = ev_t, n1i = n_arm, ci = ev_c, n2i = n_arm,
               data = data.frame(trial = 1:k, n_arm, ev_t, ev_c))
vote <- summary(vote) |>                      # adds ci.lb, ci.ub and the z test
  mutate(signif = ci.ub < 0)                  # CI excludes 1, in favour of treatment
fit_vote <- rma(yi, vi, data = vote, method = "REML")

n_sig  <- sum(vote$signif)
n_harm <- sum(vote$ci.lb > 0)
power  <- mean(map_lgl(n_arm, \(n) power.prop.test(n = n, p1 = p0, p2 = p0 * rr)$power > 0.8))


## -----------------------------------------------------------------------------
#| label: fig-vote
#| fig-cap: "Thirty-three simulated trials of a treatment with a true risk ratio of 0.80. Filled points: confidence interval excludes 1. The diamond is the random-effects pooled estimate."
#| fig-width: 7.5
#| fig-height: 7.5
vote_plot <- vote |>
  mutate(rr = exp(yi), lo = exp(ci.lb), hi = exp(ci.ub),
         trial = reorder(factor(trial), n_arm))
pooled <- with(fit_vote, c(est = exp(b[1]), lo = exp(ci.lb), hi = exp(ci.ub)))

ggplot(vote_plot, aes(y = trial)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_vline(xintercept = rr, colour = "#B22222", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "grey45") +
  geom_point(aes(x = rr, shape = signif, size = n_arm), colour = "#041E42") +
  annotate("polygon", x = pooled[c("lo", "est", "hi", "est")],
           y = c(-1, -0.4, -1, -1.6), fill = "#E69F00") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("Not significant", "Significant"), name = NULL) +
  scale_size_area(max_size = 3.5, guide = "none") +
  scale_x_log10(breaks = c(0.25, 0.5, 0.8, 1, 2, 4)) +
  scale_y_discrete(labels = NULL, expand = expansion(add = c(2.2, 0.6))) +
  labs(x = "Risk ratio (log scale); trials ordered by size", y = NULL) +
  theme(legend.position = "top", panel.grid.minor = element_blank())


## -----------------------------------------------------------------------------
#| label: fe-hand
logRR <- c(-0.35, -0.28, -0.41)               # three studies, log risk ratios
se    <- c(0.12, 0.08, 0.15)
w     <- 1 / se^2
theta <- sum(w * logRR) / sum(w)
se_th <- sqrt(1 / sum(w))

fe_check <- rma(yi = logRR, sei = se, method = "EE")   # metafor's "equal effects" model
tibble::tibble(
  Study  = c("1", "2", "3", "Pooled (by hand)", "Pooled (rma)"),
  Weight = c(sprintf("%.0f%%", 100 * w / sum(w)), "", ""),
  RR     = sprintf("%.2f", exp(c(logRR, theta, fe_check$b[1]))),
  `95% CI` = c(sprintf("%.2f to %.2f", exp(logRR - 1.96 * se), exp(logRR + 1.96 * se)),
               sprintf("%.2f to %.2f", exp(theta - 1.96 * se_th), exp(theta + 1.96 * se_th)),
               sprintf("%.2f to %.2f", exp(fe_check$ci.lb), exp(fe_check$ci.ub)))) |>
  kable(align = "lrrr")


## -----------------------------------------------------------------------------
#| label: bcg-het
bcg <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg, slab = paste(author, year))
m_bcg  <- rma(yi, vi, data = bcg, method = "REML")
ci_het <- confint(m_bcg)$random                 # Q-profile intervals for tau^2, I^2
pr_bcg <- predict(m_bcg, transf = exp)


## -----------------------------------------------------------------------------
#| label: fig-ci-pi
#| fig-cap: "Width of the 95% confidence interval for the mean and of the 95% prediction interval, as studies are added. Simulated trials with true log risk ratios drawn from N(log 0.5, 0.56²) and within-study standard errors of 0.15 to 0.6; median of 500 simulations at each size. The dashed line is the width of the true 95% range of effects."
#| fig-width: 7.5
#| fig-height: 4.2
sim_width <- function(k, tau = 0.56, mu = log(0.5)) {
  se  <- runif(k, 0.15, 0.6)
  yi  <- rnorm(k, mu, tau) + rnorm(k, 0, se)
  fit <- rma(yi, sei = se, method = "REML", test = "knha", control = ctl)
  p   <- predict(fit)
  c(ci = fit$ci.ub - fit$ci.lb, pi = p$pi.ub - p$pi.lb)
}
set.seed(2026)
widths <- map_dfr(c(5, 10, 20, 40, 80), \(k) {
  w <- replicate(500, tryCatch(sim_width(k), error = \(e) c(ci = NA, pi = NA)))
  tibble::tibble(k = k, `Confidence interval for the mean` = median(w["ci", ], na.rm = TRUE),
                 `Prediction interval` = median(w["pi", ], na.rm = TRUE))
}) |>
  pivot_longer(-k, names_to = "interval", values_to = "width")

ggplot(widths, aes(k, width, colour = interval)) +
  geom_hline(yintercept = 2 * 1.96 * 0.56, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_log10(breaks = c(5, 10, 20, 40, 80)) +
  scale_colour_manual(values = c("#041E42", "#B22222"), name = NULL) +
  labs(x = "Number of studies (log scale)", y = "Width on the log RR scale") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: fig-mcid
#| fig-cap: "Five hypothetical risk ratios (95% CI) compared with no effect (RR 1.0) and a minimal clinically important difference (RR 0.85)."
#| fig-width: 8.5
#| fig-height: 3.4
scen <- tibble::tibble(
  id    = LETTERS[1:5],
  label = c("Not significant; important effect unlikely",
            "Not significant; important effect possible",
            "Significant; not clinically important",
            "Significant; important effect possible",
            "Significant and clinically important"),
  lo = c(0.884, 0.712, 0.865, 0.757, 0.609),
  hi = c(1.118, 1.074, 0.989, 0.926, 0.777)) |>
  mutate(est = sqrt(lo * hi), y = factor(id, levels = rev(LETTERS[1:5])))

ggplot(scen, aes(y = y)) +
  geom_vline(xintercept = 1, colour = "#041E42") +
  geom_vline(xintercept = 0.85, colour = "#E69F00", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "#041E42") +
  geom_point(aes(x = est), size = 2.4, colour = "#B22222") +
  geom_text(aes(x = 1.17, label = label), hjust = 0, size = 3.5) +
  annotate("text", x = 0.84, y = 5.6, label = "MCID", hjust = 1, colour = "#b07800") +
  annotate("text", x = 1.01, y = 5.6, label = "No effect", hjust = 0) +
  scale_x_log10(limits = c(0.58, 2.3), breaks = c(0.6, 0.7, 0.85, 1, 1.15)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.9))) +
  labs(x = "Risk ratio (95% CI), log scale", y = NULL) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())


## -----------------------------------------------------------------------------
#| label: heparin-data
heparin <- tibble::tribble(
  ~study,           ~ev_asa, ~n_asa, ~ev_hep, ~n_hep,
  "Theroux 1988",        4,    121,       2,    122,
  "RISC 1990",           7,    189,       3,    210,
  "Cohen 1990",          1,     32,       0,     37,
  "Cohen 1994",          9,    109,       4,    105,
  "Holdright 1994",     40,    131,      42,    154,
  "Gurfinkel 1995",      7,     73,       4,     70)

# escalc() adds 0.5 to all four cells of a trial with a zero cell (Cohen 1990)
hep_rr <- escalc(measure = "RR", ai = ev_hep, n1i = n_hep, ci = ev_asa, n2i = n_asa,
                 data = heparin, slab = study)
hep_or <- escalc(measure = "OR", ai = ev_hep, n1i = n_hep, ci = ev_asa, n2i = n_asa,
                 data = heparin, slab = study)

heparin |>
  mutate(Aspirin = sprintf("%d/%d (%.1f%%)", ev_asa, n_asa, 100 * ev_asa / n_asa),
         `Aspirin + heparin` = sprintf("%d/%d (%.1f%%)", ev_hep, n_hep, 100 * ev_hep / n_hep)) |>
  select(Trial = study, Aspirin, `Aspirin + heparin`) |>
  kable(align = "lrr", caption = "Myocardial infarction or death during treatment, from @oler1996.")


## -----------------------------------------------------------------------------
#| label: heparin-reproduce
dl_rr <- rma(yi, vi, data = hep_rr, method = "DL")
dl_or <- rma(yi, vi, data = hep_or, method = "DL")
f_ci  <- \(m) sprintf("%.2f (%.2f to %.2f)", exp(m$b[1]), exp(m$ci.lb), exp(m$ci.ub))
tibble::tibble(Analysis = c("As published", "DerSimonian-Laird, risk ratio",
                            "DerSimonian-Laird, odds ratio"),
               `Summary (95% CI)` = c("0.67 (0.44 to 1.02)", f_ci(dl_rr), f_ci(dl_or))) |>
  kable()


## -----------------------------------------------------------------------------
#| label: heparin-grid
ests <- c("DL", "HE", "PM", "ML", "REML", "EB", "SJ")
fit_one <- function(method, dat) {
  z  <- rma(yi, vi, data = dat, method = method, control = ctl)
  hk <- rma(yi, vi, data = dat, method = method, test = "adhoc", control = ctl)
  pk <- predict(hk, transf = exp)
  tibble::tibble(method, tau2 = z$tau2, est = exp(z$b[1]),
                 wald_lo = exp(z$ci.lb), wald_hi = exp(z$ci.ub),
                 hk_lo = exp(hk$ci.lb), hk_hi = exp(hk$ci.ub),
                 pi_lo = pk$pi.lb, pi_hi = pk$pi.ub)
}
grid <- bind_rows(
  map_dfr(ests, fit_one, dat = hep_or) |> mutate(scale = "Odds ratio"),
  map_dfr(ests, fit_one, dat = hep_rr) |> mutate(scale = "Risk ratio"))
g <- \(m, s = "Odds ratio") filter(grid, method == m, scale == s)


## -----------------------------------------------------------------------------
#| label: fig-heparin-choices
#| fig-cap: "The heparin meta-analysis repeated with seven estimators of τ², on the odds ratio and risk ratio scales. Thick bar: Wald 95% CI. Thin bar: HKSJ 95% CI (never narrower than Wald). Red: 95% prediction interval. Labels give the estimated τ²."
#| fig-width: 8.5
#| fig-height: 5
grid |>
  mutate(method = factor(method, levels = rev(ests))) |>
  ggplot(aes(y = method)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_errorbarh(aes(xmin = pi_lo, xmax = pi_hi), height = 0, colour = "#B22222",
                 linewidth = 0.6, position = position_nudge(y = -0.25)) +
  geom_errorbarh(aes(xmin = hk_lo, xmax = hk_hi), height = 0, colour = "#041E42",
                 linewidth = 0.5) +
  geom_errorbarh(aes(xmin = wald_lo, xmax = wald_hi), height = 0, colour = "#041E42",
                 linewidth = 2.2) +
  geom_point(aes(x = est), shape = 21, fill = "#E69F00", size = 2.5) +
  geom_text(aes(x = 4.2, label = sprintf("tau^2 == '%.3f'", tau2)), parse = TRUE,
            hjust = 1, size = 3.2) +
  facet_wrap(~ scale) +
  scale_x_log10(limits = c(0.25, 4.2), breaks = c(0.25, 0.5, 1, 2)) +
  labs(x = "Summary (aspirin + heparin vs aspirin), log scale",
       y = "Estimator of the between-study variance")


## -----------------------------------------------------------------------------
#| label: heparin-mh-peto
iv   <- rma(yi, vi, data = hep_or, method = "EE")
mh   <- rma.mh(measure = "OR", ai = ev_hep, n1i = n_hep, ci = ev_asa, n2i = n_asa, data = heparin)
peto <- rma.peto(ai = ev_hep, n1i = n_hep, ci = ev_asa, n2i = n_asa, data = heparin)
tibble::tibble(Method = c("Inverse variance (0.5 added to Cohen 1990)", "Mantel-Haenszel", "Peto"),
               `Odds ratio (95% CI)` = c(f_ci(iv), f_ci(mh), f_ci(peto))) |>
  kable(caption = "Fixed-effect odds ratio for the heparin trials, three ways.")


## -----------------------------------------------------------------------------
#| label: simpson-toy
toy <- tibble::tibble(trial = c("A (high risk, 4:1 allocation)", "B (low risk, 1:4 allocation)"),
                      ev_t = c(160, 8), n_t = c(400, 100), ev_c = c(50, 40), n_c = c(100, 400))
toy_or <- escalc(measure = "OR", ai = ev_t, n1i = n_t, ci = ev_c, n2i = n_c, data = toy)
odds   <- \(e, n) e / (n - e)
added  <- with(toy, odds(sum(ev_t), sum(n_t)) / odds(sum(ev_c), sum(n_c)))
mh_toy <- rma.mh(measure = "OR", ai = ev_t, n1i = n_t, ci = ev_c, n2i = n_c, data = toy)
tibble::tibble(Analysis = c(toy$trial, "Both trials added into one table", "Mantel-Haenszel, by trial"),
               `Odds ratio` = sprintf("%.2f", c(exp(toy_or$yi), added, exp(mh_toy$b[1])))) |>
  kable(align = "lr")

