# R code for chapter: 15-publication-bias.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch15
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
#| label: fig-pub-mechanism
#| fig-cap: "Two hundred simulated studies of a small true effect (log OR 0.15, dashed green line). Significant studies are published with probability 0.95, others with probability 0.30. Left: all studies, with the unpublished in red. Right: the published studies alone, whose mean (red line) overstates the truth."
#| fig-width: 10
#| fig-height: 4
#| out-width: "100%"
set.seed(1)
n_st <- 200; true_effect <- 0.15
study_se <- runif(n_st, 0.05, 0.40)
observed <- rnorm(n_st, true_effect, study_se)
p_value  <- 2 * pnorm(-abs(observed / study_se))
published <- rbinom(n_st, 1, ifelse(p_value < 0.05, 0.95, 0.30)) == 1
df_pub <- tibble::tibble(observed, study_se, status = ifelse(published, "Published", "Unpublished"))
p1 <- ggplot(df_pub, aes(observed, 1 / study_se, colour = status)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = true_effect, colour = "#2E7D32", linetype = "dashed") +
  scale_colour_manual(values = c(Published = "#041E42", Unpublished = "#B22222"), name = NULL) +
  labs(x = "Observed log OR", y = "Precision (1/SE)", title = "All 200 studies") +
  theme(legend.position = "bottom")
p2 <- ggplot(filter(df_pub, status == "Published"), aes(observed, 1 / study_se)) +
  geom_point(colour = "#041E42", alpha = 0.7) +
  geom_vline(xintercept = true_effect, colour = "#2E7D32", linetype = "dashed") +
  geom_vline(xintercept = mean(observed[published]), colour = "#B22222") +
  labs(x = "Observed log OR", y = NULL, title = "Published studies only")
patchwork::wrap_plots(p1, p2)


## -----------------------------------------------------------------------------
#| label: fig-funnel-theory
#| fig-cap: "Simulated funnel plots of 30 studies with a true effect of 0.3. Left: all studies. Right: small studies (SE > 0.25) with negative results suppressed. Red line: inverse-variance pooled estimate."
#| fig-width: 10
#| fig-height: 4
#| out-width: "100%"
set.seed(42)
se_f <- runif(30, 0.05, 0.50)
y_f  <- rnorm(30, 0.3, se_f)
fun_df <- bind_rows(tibble::tibble(y = y_f, se = se_f, data = "No bias"),
                    tibble::tibble(y = y_f, se = se_f, data = "Small negative studies missing") |>
                      filter(!(y < 0 & se > 0.25)))
pools <- fun_df |> summarise(.by = data, pool = weighted.mean(y, 1 / se^2))
cones <- pools |> crossing(se = seq(0, 0.5, length.out = 50)) |>
  mutate(lo = pool - 1.96 * se, hi = pool + 1.96 * se)
ggplot(fun_df, aes(y, se)) +
  geom_line(data = cones, aes(lo, se), linetype = "dashed", colour = "grey60") +
  geom_line(data = cones, aes(hi, se), linetype = "dashed", colour = "grey60") +
  geom_point(colour = "#041E42", size = 2) +
  geom_vline(data = pools, aes(xintercept = pool), colour = "#B22222") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ data) +
  scale_y_reverse() +
  labs(x = "Effect size", y = "Standard error")


## -----------------------------------------------------------------------------
#| label: bcg-asym
bcg   <- escalc("RR", ai = tpos, n1i = tpos + tneg, ci = cpos, n2i = cpos + cneg, data = dat.bcg)
m_bcg <- rma(yi, vi, data = bcg, method = "REML")
eg <- regtest(m_bcg, model = "lm", predictor = "sei")      # Egger
bg <- ranktest(m_bcg)                                      # Begg
pt <- regtest(m_bcg, model = "lm", predictor = "ninv")     # Peters: effect on 1/n
pt_wrong <- regtest(m_bcg, model = "lm", predictor = "ni") # a common error: effect on n
round(c(Egger = eg$pval, Begg = bg$pval, Peters = pt$pval, `Peters on n (wrong)` = pt_wrong$pval), 3)


## -----------------------------------------------------------------------------
#| label: fig-trimfill
#| fig-cap: "Trim-and-fill applied to the BCG trials. Filled points would be imputed studies."
#| fig-width: 7
#| fig-height: 4.5
tf <- trimfill(m_bcg)
funnel(tf, xlab = "Log risk ratio (BCG vs control)")


## -----------------------------------------------------------------------------
#| label: pet-peese
set.seed(42)
se_i <- runif(25, 0.08, 0.45)
yi   <- rnorm(25, 0.20, se_i)                      # true effect 0.20
pub  <- !(yi < 0.05 & se_i > 0.25)                 # small, disappointing trials stay in the drawer
d_pp <- tibble::tibble(yi = yi[pub], sei = se_i[pub], vi = sei^2)
m_naive <- rma(yi, vi, data = d_pp)
m_pet   <- rma(yi, vi, mods = ~ sei, data = d_pp)
m_peese <- rma(yi, vi, mods = ~ vi, data = d_pp)
p_pet   <- coef(summary(m_pet))["intrcpt", "pval"]


## -----------------------------------------------------------------------------
#| label: fig-pet-peese
#| fig-cap: "PET (red, straight line in SE) and PEESE (orange, curve in SE²) fitted to the published trials and extrapolated to SE = 0. Green dashed: the true effect (0.20); grey dotted: naive pooled estimate."
#| fig-width: 7.5
#| fig-height: 4.2
b_pet <- coef(m_pet); b_peese <- coef(m_peese)
ggplot(d_pp, aes(sei, yi)) +
  geom_function(fun = \(x) b_pet[1] + b_pet[2] * x, xlim = c(0, 0.45), colour = "#B22222", linewidth = 1) +
  geom_function(fun = \(x) b_peese[1] + b_peese[2] * x^2, xlim = c(0, 0.45), colour = "#E69F00", linewidth = 1) +
  geom_hline(yintercept = 0.20, colour = "#2E7D32", linetype = "dashed") +
  geom_hline(yintercept = m_naive$b[1], colour = "grey50", linetype = "dotted") +
  geom_point(colour = "#041E42", size = 2.3) +
  annotate("point", x = 0, y = c(b_pet[1], b_peese[1]), size = 3.5, colour = c("#B22222", "#E69F00")) +
  labs(x = "Standard error (larger = smaller study)", y = "Observed effect")


## -----------------------------------------------------------------------------
#| label: lit-sim
make_lit <- function(delta, k_all = 150, seed = 629) {
  set.seed(seed)
  n     <- round(runif(k_all, 15, 120))                    # patients per arm
  sei   <- sqrt(2 / n)                                     # approximate SE of an SMD
  theta <- delta + rnorm(k_all, 0, 0.05)                   # each trial's true effect
  yi    <- rnorm(k_all, theta, sei)
  p     <- 2 * pnorm(-abs(yi / sei))
  sig   <- p < 0.05 & yi > 0
  keep  <- sample(which(!sig), round(0.15 * sum(!sig)))    # 15% of the rest are published
  all   <- data.frame(yi, sei, vi = sei^2, z = yi / sei, p, theta,
                      published = sig | seq_len(k_all) %in% keep)
  structure(all[all$published, ], all = all)
}
lit  <- make_lit(0.30)
lit0 <- make_lit(0)
naive  <- rma(yi, vi, data = lit)
naive0 <- rma(yi, vi, data = lit0)


## -----------------------------------------------------------------------------
#| label: selmodel
m_ml <- rma(yi, vi, data = lit, method = "ML")
p_one <- pnorm(lit$yi / lit$sei, lower.tail = FALSE)            # one-sided P values
table(cut(p_one, c(0, 0.025, 0.05, 0.10, 0.50, 1)))               # published studies per band
sel  <- selmodel(m_ml, type = "stepfun", steps = c(0.025, 0.05, 0.10, 0.50))
sel3 <- selmodel(m_ml, type = "stepfun", steps = 0.025)           # one cut-point: "3PSM"
tibble::tibble(Model = c("Naive (ML)", "Step function, four cut-points", "Step function, one cut-point (3PSM)"),
               `Pooled SMD (95% CI)` = sprintf("%.2f (%.2f to %.2f)",
                                               c(m_ml$b[1], sel$b[1], sel3$b[1]),
                                               c(m_ml$ci.lb, sel$ci.lb, sel3$ci.lb),
                                               c(m_ml$ci.ub, sel$ci.ub, sel3$ci.ub)),
               `Relative publication weights` = c("", paste(sprintf("%.2f", sel$delta), collapse = ", "),
                                                  paste(sprintf("%.2f", sel3$delta), collapse = ", ")),
               `LRT for selection, P` = c("", format.pval(sel$LRTp, digits = 2, eps = 0.001),
                                          format.pval(sel3$LRTp, digits = 2, eps = 0.001))) |>
  kable(align = "lccc")


## -----------------------------------------------------------------------------
#| label: fig-pcurves
#| fig-cap: "Distribution of significant P values in 0.01 bins. No effect: flat. A real effect with 80% power: right-skewed. No effect plus optional stopping (start with 20 per arm, add 10 up to five times until P < 0.05): left-skewed. Dashed line: 20% per bin."
#| fig-width: 10
#| fig-height: 3.5
#| out-width: "100%"
set.seed(42)
p_null <- 2 * pnorm(-abs(rnorm(20000)))
p_true <- 2 * pnorm(-abs(rnorm(20000, 2.8, 1)))
p_hack <- replicate(4000, {
  x <- rnorm(70); y <- rnorm(70)
  for (n in seq(20, 70, by = 10)) {
    p <- t.test(x[1:n], y[1:n], var.equal = TRUE)$p.value
    if (p < 0.05) break
  }
  p
})
tibble::tibble(p = c(p_null, p_true, p_hack),
               scenario = rep(c("No effect", "Real effect, 80% power", "No effect + p-hacking"),
                              c(20000, 20000, 4000))) |>
  filter(p < 0.05) |>
  mutate(bin = cut(p, seq(0, 0.05, 0.01), labels = c(".01", ".02", ".03", ".04", ".05")),
         scenario = factor(scenario, c("No effect", "Real effect, 80% power", "No effect + p-hacking"))) |>
  count(scenario, bin) |>
  mutate(.by = scenario, pct = 100 * n / sum(n)) |>
  ggplot(aes(bin, pct, group = scenario)) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1, colour = "#041E42") + geom_point(colour = "#041E42") +
  facet_wrap(~ scenario) +
  labs(x = "P value bin (significant results only)", y = "% of significant results")


## -----------------------------------------------------------------------------
#| label: pcurve-fn
pcurve_tests <- function(d) {
  s   <- subset(d, p < 0.05 & yi > 0)
  z   <- s$yi / s$sei
  pp0 <- s$p / 0.05                                         # uniform under no effect
  ncp <- qnorm(0.975) - qnorm(2 / 3)                        # z-shift giving 33% power
  pp33 <- (pnorm(z, ncp, lower.tail = FALSE) + pnorm(-z, ncp)) / (1 / 3)  # uniform under 33% power
  Z_rs <- sum(qnorm(pp0)) / sqrt(nrow(s))
  Z_fl <- sum(qnorm(pmin(pp33, 1 - 1e-12))) / sqrt(nrow(s))
  c(k_significant = nrow(s), `share P < 0.025` = mean(s$p < 0.025),
    `right-skew test P` = pnorm(Z_rs), `flatness (33% power) test P` = pnorm(Z_fl, lower.tail = FALSE))
}
round(rbind(`True SMD 0.30` = pcurve_tests(lit), `True SMD 0` = pcurve_tests(lit0)), 3)


## -----------------------------------------------------------------------------
#| label: puniform-fn
puniform_simple <- function(d) {
  s <- subset(d, p < 0.05 & yi > 0)
  crit <- qnorm(0.975) * s$sei                             # smallest significant effect
  cond_p <- \(delta) pnorm((s$yi - delta) / s$sei, lower.tail = FALSE) /
                     pnorm((crit - delta) / s$sei, lower.tail = FALSE)
  k <- nrow(s); half_w <- qnorm(0.975) * sqrt(k / 12)      # normal approximation to Irwin-Hall
  f <- \(delta, target) sum(cond_p(delta)) - target
  c(estimate = uniroot(f, c(-2, 2), target = k / 2)$root,
    lower = uniroot(f, c(-2, 2), target = k / 2 - half_w)$root,
    upper = uniroot(f, c(-2, 2), target = k / 2 + half_w)$root)
}
round(rbind(`True SMD 0.30` = puniform_simple(lit), `True SMD 0` = puniform_simple(lit0)), 2)


## -----------------------------------------------------------------------------
#| label: puniform-pkg
#| eval: !expr requireNamespace("puniform", quietly = TRUE)
library(puniform)
puniform(yi = lit$yi, vi = lit$vi, side = "right", method = "P")
puni_star(yi = lit$yi, vi = lit$vi, side = "right")


## -----------------------------------------------------------------------------
#| label: fig-eta
#| fig-cap: "Corrected pooled SMD as the assumed selection ratio η increases (fixed-effect version: non-affirmative studies up-weighted by η). Green dashed: the truth (0.30). Red dotted: the selection ratio used in the simulation (6.7)."
#| fig-width: 7.5
#| fig-height: 4
eta_fit <- function(d, eta) {
  aff <- d$p < 0.05 & d$yi > 0
  w   <- ifelse(aff, 1, eta) / d$vi
  est <- sum(w * d$yi) / sum(w)
  se  <- sqrt(sum(w^2 * d$vi)) / sum(w)
  c(eta = eta, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}
etas <- c(1, 1.5, 2, 3, 4, 6.7, 10, 20, 50, 200)
sens <- map_dfr(etas, \(e) eta_fit(lit, e))
ggplot(sens, aes(eta, est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#041E42", alpha = 0.15) +
  geom_line(colour = "#041E42", linewidth = 1) + geom_point(colour = "#041E42") +
  geom_hline(yintercept = 0.30, linetype = "dashed", colour = "#2E7D32") +
  geom_vline(xintercept = 1 / 0.15, linetype = "dotted", colour = "#B22222") +
  scale_x_log10(breaks = c(1, 2, 4, 10, 50, 200)) +
  labs(x = "Selection ratio η (log scale)", y = "Corrected pooled SMD (95% CI)")


## -----------------------------------------------------------------------------
#| label: robma
#| eval: !expr requireNamespace("RoBMA", quietly = TRUE)
library(RoBMA)                     # requires JAGS
fit1 <- RoBMA(d = lit$yi, se = lit$sei, algorithm = "ss", seed = 1)
summary(fit1)


## -----------------------------------------------------------------------------
#| label: fig-all-methods
#| fig-cap: "Several methods applied to the same published literature (true SMD 0.30, green dashed line)."
#| fig-width: 8
#| fig-height: 4
pet_l   <- rma(yi, vi, mods = ~ sei, data = lit)
peese_l <- rma(yi, vi, mods = ~ vi,  data = lit)
pp_l    <- if (coef(summary(pet_l))["intrcpt", "pval"] < 0.10) peese_l else pet_l
sel_l   <- selmodel(rma(yi, vi, data = lit, method = "ML"), type = "stepfun", steps = 0.025)
tf_l    <- trimfill(naive)
pu_l    <- puniform_simple(lit)
eta_l   <- eta_fit(lit, 1 / 0.15)
tibble::tribble(
  ~method,                               ~est,          ~lo,              ~hi,
  "Naive random effects",                naive$b[1],    naive$ci.lb,      naive$ci.ub,
  "Trim-and-fill",                       tf_l$b[1],     tf_l$ci.lb,       tf_l$ci.ub,
  "PET-PEESE (two-step)",                pp_l$b[1],     pp_l$ci.lb[1],    pp_l$ci.ub[1],
  "Selection model (3PSM)",              sel_l$b[1],    sel_l$ci.lb[1],   sel_l$ci.ub[1],
  "p-uniform (simple)",                  pu_l[1],       pu_l[2],          pu_l[3],
  "Selection ratio η = 6.7 (fixed effect)", eta_l["est"], eta_l["lo"],   eta_l["hi"]) |>
  mutate(method = factor(method, levels = rev(method))) |>
  ggplot(aes(est, method)) +
  geom_vline(xintercept = 0.30, colour = "#2E7D32", linetype = "dashed") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#041E42") +
  labs(x = "Estimated SMD (95% interval)", y = NULL)

