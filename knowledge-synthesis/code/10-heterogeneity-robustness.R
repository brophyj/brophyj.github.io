# R code for chapter: 10-heterogeneity-robustness.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch10
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
ctl <- list(stepadj = 0.5, maxiter = 1000)


## -----------------------------------------------------------------------------
#| label: bcg-setup
bcg   <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                data = dat.bcg, slab = paste(author, year))
m_bcg <- rma(yi, vi, data = bcg, method = "REML", test = "knha")


## -----------------------------------------------------------------------------
#| label: bcg-metareg
m_reg <- rma(yi, vi, mods = ~ ablat, data = bcg, method = "REML", test = "knha")
m_reg


## -----------------------------------------------------------------------------
#| label: bcg-metareg-numbers
#| include: false
b      <- coef(m_reg)
r2     <- m_reg$R2
rr_deg <- exp(b[2])
v_typ  <- with(m_reg, tau2 * (100 - I2) / I2)       # "typical" within-trial variance
pr_eq  <- exp(b[1] + 0 * b[2]); pr_50 <- exp(b[1] + 50 * b[2])
m_reg_z <- rma(yi, vi, mods = ~ ablat, data = bcg, method = "REML")


## -----------------------------------------------------------------------------
#| label: fig-bcg-regplot
#| fig-cap: "Meta-regression of the BCG trials on absolute latitude. Point area is proportional to the trial's weight. Inner band: 95% confidence interval for the mean log RR at each latitude. Outer band: 95% prediction interval for a new trial."
#| fig-width: 7.5
#| fig-height: 5
regplot(m_reg, mod = "ablat", pi = TRUE, shade = c("grey80", "grey92"),
        xlab = "Absolute latitude (degrees from the equator)",
        ylab = "Log risk ratio (BCG vs control)", bg = "#B22222")
abline(h = 0, lty = 2, col = "grey40")


## -----------------------------------------------------------------------------
#| label: fig-bcg-dag
#| fig-cap: "Latitude affects the baseline risk of tuberculosis (grey), which randomisation balances between arms, and modifies the effect of BCG (red, an arrow into the BCG-to-TB arrow). There is no arrow from latitude to BCG within a trial."
#| fig-width: 6.5
#| fig-height: 4
#| out-width: "70%"
nodes <- tibble::tibble(x = c(0, 3, 1.5), y = c(0, 0, 2),
                        lab = c("BCG\n(randomised)", "Tuberculosis", "Latitude\n(mycobacterial\nexposure)"))
ggplot() +
  annotate("segment", x = 0.6, y = 0, xend = 2.35, yend = 0, linewidth = 1.2,
           arrow = arrow(length = unit(0.3, "cm")), colour = "#041E42") +
  annotate("segment", x = 1.95, y = 1.6, xend = 2.8, yend = 0.35, linewidth = 0.8,
           arrow = arrow(length = unit(0.25, "cm")), colour = "grey50") +
  annotate("segment", x = 1.5, y = 1.45, xend = 1.5, yend = 0.1, linewidth = 1,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"), colour = "#B22222") +
  annotate("text", x = 1.6, y = 0.8, label = "effect\nmodification", hjust = 0,
           size = 3.6, colour = "#B22222", lineheight = 0.9) +
  annotate("text", x = 2.55, y = 1.15, label = "baseline\nrisk", size = 3.4,
           colour = "grey35", lineheight = 0.9) +
  geom_label(data = nodes, aes(x, y, label = lab), size = 4, fill = "#F4F6F8", lineheight = 0.9) +
  coord_cartesian(xlim = c(-0.5, 3.5), ylim = c(-0.4, 2.5)) +
  theme_void()


## -----------------------------------------------------------------------------
#| label: bcg-confounded
fit_mod <- function(f) {
  m <- rma(yi, vi, mods = f, data = bcg, method = "REML", test = "knha")
  tibble::tibble(Moderators = deparse(f[[2]]) |> sub(pattern = "ablat", replacement = "latitude") |>
                   sub(pattern = "alloc", replacement = "allocation method"), `R² (%)` = sprintf("%.0f", m$R2),
                 `Residual I² (%)` = sprintf("%.0f", m$I2),
                 `Moderator test P` = sprintf("%.3f", m$QMp))
}
map_dfr(list(~ ablat, ~ year, ~ alloc, ~ ablat + year), fit_mod) |>
  kable(align = "lrrr",
        caption = sprintf("Meta-regressions of the BCG trials (REML, Knapp-Hartung tests). Latitude and year are correlated (r = %.2f).",
                          cor(bcg$ablat, bcg$year)))


## -----------------------------------------------------------------------------
#| label: dredge-sim
dredge_once <- function(k = 10, n_mods = 5, tau2 = 0.1) {
  vi <- runif(k, 0.02, 0.2)
  yi <- rnorm(k, log(0.8), sqrt(tau2 + vi))           # no moderator affects yi
  X  <- matrix(rnorm(k * n_mods), k)
  p  <- map(seq_len(n_mods), \(j) {
    c(z    = rma(yi, vi, mods = ~ X[, j], method = "REML", control = ctl)$QMp,
      knha = rma(yi, vi, mods = ~ X[, j], method = "REML", test = "knha", control = ctl)$QMp)
  }) |> do.call(what = rbind)
  c(z_one = unname(p[1, "z"] < 0.05), knha_one = unname(p[1, "knha"] < 0.05),
    z_any = any(p[, "z"] < 0.05), knha_any = any(p[, "knha"] < 0.05))
}
set.seed(10)
dredge <- replicate(1000, tryCatch(dredge_once(), error = \(e) rep(NA, 4))) |>
  rowMeans(na.rm = TRUE)


## -----------------------------------------------------------------------------
#| label: fig-dredge
#| fig-cap: "Probability of at least one 'significant' moderator (P < 0.05) when none is real. 1000 simulated meta-analyses of 10 trials with τ² = 0.1; default (z) and Knapp-Hartung tests."
#| fig-width: 7
#| fig-height: 3.4
tibble::tibble(test = rep(c("Default (z)", "Knapp-Hartung"), 2),
               scope = rep(c("One pre-specified moderator", "Best of five moderators"), each = 2),
               p = dredge[c("z_one", "knha_one", "z_any", "knha_any")]) |>
  mutate(scope = factor(scope, levels = c("One pre-specified moderator", "Best of five moderators"))) |>
  ggplot(aes(p, scope, fill = test)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * p)), position = position_dodge(width = 0.7),
            hjust = -0.15, size = 3.6) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("grey60", "#041E42"), name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.45)) +
  labs(x = "False-positive rate", y = NULL) +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: fig-simpson-eco
#| fig-cap: "Ecological Simpson's paradox (simulated data). Solid lines: within each region, more circumcision goes with less HIV. Dashed line: ignoring region, the association reverses."
#| fig-width: 7
#| fig-height: 4.6
set.seed(42)
eco <- tibble::tibble(
  region = rep(c("Region A: high circumcision, high HIV",
                 "Region B: low circumcision, low HIV"), each = 12),
  circ   = c(runif(12, 0.70, 0.95), runif(12, 0.10, 0.35))) |>
  mutate(hiv = if_else(startsWith(region, "Region A"), 0.40, 0.12) - 0.15 * circ +
           rnorm(n(), 0, 0.015))
b_all    <- coef(lm(hiv ~ circ, data = eco))[["circ"]]
b_within <- coef(lm(hiv ~ circ + region, data = eco))[["circ"]]

ggplot(eco, aes(circ, hiv)) +
  geom_point(aes(colour = region), size = 2.6, alpha = 0.85) +
  geom_smooth(aes(colour = region), method = "lm", se = FALSE, formula = y ~ x) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, colour = "#B22222",
              linetype = "dashed") +
  scale_colour_manual(values = c("#E69F00", "#041E42"), name = NULL) +
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Circumcision prevalence (country)", y = "HIV prevalence (country)") +
  theme(legend.position = "bottom", legend.direction = "vertical")


## -----------------------------------------------------------------------------
#| label: childbirth
cb <- data.frame(study = paste("Trial", 1:9),
                 ai  = c(10, 18, 22, 13, 40, 25,  8, 30, 15),   # epidurals, support
                 n1i = c(55, 80, 95, 60, 135, 110, 45, 120, 70),
                 ci  = c(22, 30, 35, 18, 55, 38, 15, 48, 24),   # epidurals, usual care
                 n2i = c(50, 75, 90, 58, 130, 105, 42, 118, 68))
cb_rr <- escalc("RR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = cb, slab = study)

fits <- list(
  "Equal effects (inverse variance)" = rma(yi, vi, data = cb_rr, method = "EE"),
  "Random effects, DL, Wald"         = rma(yi, vi, data = cb_rr, method = "DL"),
  "Random effects, REML, Wald"       = rma(yi, vi, data = cb_rr, method = "REML"),
  "Random effects, REML, HKSJ"       = rma(yi, vi, data = cb_rr, method = "REML", test = "knha"),
  "Random effects, REML, HKSJ (ad hoc)" = rma(yi, vi, data = cb_rr, method = "REML", test = "adhoc"),
  "Mantel-Haenszel"                  = rma.mh(measure = "RR", ai = ai, n1i = n1i, ci = ci,
                                              n2i = n2i, data = cb))
cb_tab <- imap_dfr(fits, \(m, nm) tibble::tibble(
  Method = nm, RR = exp(m$b[1]), lo = exp(m$ci.lb), hi = exp(m$ci.ub)))
cb_tab |>
  mutate(`RR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", RR, lo, hi),
         `CI width (log scale)` = sprintf("%.3f", log(hi) - log(lo))) |>
  select(Method, `RR (95% CI)`, `CI width (log scale)`) |>
  kable(align = "lrr")


## -----------------------------------------------------------------------------
#| label: childbirth-numbers
#| include: false
m_dl  <- fits[["Random effects, DL, Wald"]]
q_df  <- m_dl$QE / (m_dl$k - 1)
w_rat <- with(cb_tab, (log(hi[4]) - log(lo[4])) / (log(hi[3]) - log(lo[3])))


## -----------------------------------------------------------------------------
#| label: fig-bcg-loo
#| fig-cap: "Leave-one-out analysis of the BCG trials (REML, HKSJ). Each row is the pooled risk ratio with that trial omitted. Dashed line: the pooled risk ratio from all 13 trials."
#| fig-width: 7.5
#| fig-height: 5
loo <- leave1out(m_bcg, transf = exp)
forest(loo$estimate, ci.lb = loo$ci.lb, ci.ub = loo$ci.ub, slab = loo$slab,
       refline = exp(m_bcg$b[1]), header = c("Trial omitted", "RR [95% CI]"),
       xlab = "Pooled risk ratio with one trial omitted")


## -----------------------------------------------------------------------------
#| label: lau-cumul
lau   <- dat.lau1992
m_lau <- rma(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = lau,
             method = "EE", slab = paste(trial, year))
cum_lau <- cumul(m_lau)
n_cum   <- cumsum(lau$n1i + lau$n2i)
k_lau   <- length(cum_lau$estimate)
sig     <- cum_lau$ci.ub < 0                 # cumulative CI excludes OR = 1
first_sig <- max(which(!sig)) + 1            # significant from here on


## -----------------------------------------------------------------------------
#| label: fig-lau-cumul
#| fig-cap: "Cumulative fixed-effect meta-analysis of 33 trials of intravenous streptokinase for acute myocardial infarction (mortality), adding trials in order of publication. Data from Lau et al. (1992), as distributed with metafor."
#| fig-width: 8
#| fig-height: 8.5
forest(cum_lau, atransf = exp, at = log(c(0.25, 0.5, 1, 2, 4)),
       xlim = c(-7.5, 4.5), ilab = n_cum, ilab.xpos = -2.6, cex = 0.75,
       header = c("Trial added", "Cumulative OR [95% CI]"),
       xlab = "Cumulative pooled odds ratio for death")
text(-2.6, k_lau + 2, "Patients", font = 2, cex = 0.75)


## -----------------------------------------------------------------------------
#| label: fragility-fn
p_trial <- function(e1, n1, e0, n0)
  fisher.test(matrix(c(e1, n1 - e1, e0, n0 - e0), nrow = 2, byrow = TRUE))$p.value

fragility <- function(e1, n1, e0, n0) {
  k <- 0
  while (p_trial(e1, n1, e0, n0) < 0.05) {   # add events to the arm with fewer events
    if (e1 / n1 < e0 / n0) e1 <- e1 + 1 else e0 <- e0 + 1
    k <- k + 1
  }
  k
}
c(p = p_trial(30, 1000, 50, 1000), fragility_index = fragility(30, 1000, 50, 1000))


## -----------------------------------------------------------------------------
#| label: fragility-sim
set.seed(629)
sims <- tibble::tibble(e1 = rbinom(2000, 500, 0.07), e0 = rbinom(2000, 500, 0.10)) |>
  mutate(rr  = e1 / e0,
         p   = map2_dbl(e1, e0, \(a, b) p_trial(a, 500, b, 500)),
         sig = p < 0.05 & e1 < e0)
sig_trials <- sims |> filter(sig) |>
  mutate(fi = map2_dbl(e1, e0, \(a, b) fragility(a, 500, b, 500)))


## -----------------------------------------------------------------------------
#| label: fig-fragility
#| fig-cap: "Two thousand simulated trials of a drug with a true risk ratio of 0.70 (500 patients per arm, control risk 10%). Left: observed risk ratios; navy trials reached P < 0.05. Right: fragility indices of the significant trials."
#| fig-width: 10
#| fig-height: 3.8
#| out-width: "100%"
p1 <- ggplot(sims, aes(rr, fill = sig)) +
  geom_histogram(binwidth = 0.05, boundary = 1) +
  geom_vline(xintercept = 0.70, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#041E42"),
                    labels = c("Not significant", "P < 0.05"), name = NULL) +
  labs(x = "Observed risk ratio", y = "Trials") + theme(legend.position = "top")
p2 <- ggplot(sig_trials, aes(fi)) +
  geom_bar(fill = "#041E42") +
  geom_vline(xintercept = 5.5, colour = "#B22222", linetype = "dashed") +
  labs(x = "Fragility index (patients)", y = "Significant trials")
patchwork::wrap_plots(p1, p2)

