# R code for chapter: 19-observational-studies.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch19
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
  library(rms)          # rcs(): restricted cubic splines
  library(dosresmeta)   # dose-response meta-analysis
})


## -----------------------------------------------------------------------------
#| label: fig-obs-pool
#| fig-cap: "Eight simulated studies of an exposure with no effect, all sharing the same healthy-user bias (about 22% apparent protection). The pooled estimate is precise, consistent and wrong."
#| fig-width: 7
#| fig-height: 4
set.seed(42)
k  <- 8
lr <- rnorm(k, log(1.0) - 0.25, 0.08)        # truth RR = 1; shared bias of about -0.25 on the log scale
se <- runif(k, 0.06, 0.18)
pooled <- rma(lr, sei = se, method = "REML")
tibble::tibble(label = c(paste("Study", 1:k), "Pooled"), rr = exp(c(lr, pooled$b)),
               lo = exp(c(lr - 1.96 * se, pooled$ci.lb)), hi = exp(c(lr + 1.96 * se, pooled$ci.ub)),
               type = c(rep("Study", k), "Pooled")) |>
  mutate(label = factor(label, levels = rev(label))) |>
  ggplot(aes(rr, label, colour = type)) +
  geom_vline(xintercept = 1, colour = "#2E7D32", linewidth = 1) +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = type)) +
  scale_colour_manual(values = c(Study = "#041E42", Pooled = "#B22222"), guide = "none") +
  scale_shape_manual(values = c(Study = 16, Pooled = 18), guide = "none") +
  labs(x = "Risk ratio (truth = 1)", y = NULL)


## -----------------------------------------------------------------------------
#| label: fig-or-rr
#| fig-cap: "The odds ratio corresponding to a constant risk ratio of 0.70, as baseline risk increases."
#| fig-width: 6
#| fig-height: 3.5
tibble::tibble(p0 = seq(0.01, 0.5, 0.005)) |>
  mutate(p1 = 0.70 * p0, or = (p1 / (1 - p1)) / (p0 / (1 - p0))) |>
  ggplot(aes(p0, or)) +
  geom_hline(yintercept = 0.70, linetype = "dashed", colour = "#041E42") +
  geom_line(colour = "#B22222", linewidth = 1.1) +
  annotate("text", x = 0.5, y = 0.72, hjust = 1, label = "risk ratio = 0.70", colour = "#041E42") +
  annotate("text", x = 0.5, y = 0.52, hjust = 1, label = "odds ratio", colour = "#B22222") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Baseline risk", y = "Effect")


## -----------------------------------------------------------------------------
#| label: evalue
evalue_rr <- \(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev <- evalue_rr(0.78); ev_ci <- evalue_rr(0.85)          # an illustrative pooled HR 0.78 (0.72 to 0.85)
round(c(`E-value, estimate` = ev, `E-value, CI limit` = ev_ci), 2)


## -----------------------------------------------------------------------------
#| label: fig-evalue
#| fig-cap: "Combinations of confounder-exposure and confounder-outcome associations strong enough to explain away an observed HR of 0.78 (on and above the curve). Red point: the E-value. Dashed box: a plausible range for healthy-user bias."
#| fig-width: 6
#| fig-height: 5
rr_obs <- 1 / 0.78
tibble::tibble(rr_EU = seq(rr_obs + 0.01, 4, 0.005)) |>
  mutate(rr_UD = rr_obs * (rr_EU - 1) / (rr_EU - rr_obs)) |>
  filter(rr_UD <= 4) |>
  ggplot(aes(rr_EU, rr_UD)) +
  geom_ribbon(aes(ymin = rr_UD, ymax = 4), fill = "#B22222", alpha = 0.12) +
  geom_line(colour = "#1565C0", linewidth = 1.2) +
  annotate("rect", xmin = 1.5, xmax = 2, ymin = 1.5, ymax = 2, fill = NA, colour = "grey30", linetype = "dashed") +
  annotate("point", x = ev, y = ev, colour = "#B22222", size = 3) +
  annotate("text", x = ev, y = ev, hjust = -0.15, label = sprintf("E-value = %.2f", ev), colour = "#B22222") +
  annotate("text", x = 3.2, y = 3.4, label = "could explain\nthe association", colour = "#B22222") +
  coord_cartesian(xlim = c(1, 4), ylim = c(1, 4)) +
  labs(x = "Confounder-exposure association (RR)", y = "Confounder-outcome association (RR)")


## -----------------------------------------------------------------------------
#| label: statin-studies
studies <- tibble::tribble(
  ~study, ~hr,  ~se,  ~robins,
  "A", 0.76, 0.07, "Moderate", "B", 0.94, 0.06, "Serious",  "C", 0.73, 0.08, "Serious",
  "D", 0.84, 0.07, "Critical", "E", 0.97, 0.08, "Moderate", "F", 0.79, 0.06, "Serious",
  "G", 0.69, 0.09, "Critical", "H", 0.87, 0.05, "Moderate", "I", 0.75, 0.07, "Serious",
  "J", 0.90, 0.06, "Serious") |>
  mutate(robins = factor(robins, levels = c("Moderate", "Serious", "Critical")), label = paste("Study", study))
res     <- rma(yi = log(hr), sei = se, data = studies, slab = label, method = "REML")
res_mod <- rma(yi = log(hr), sei = se, data = studies, subset = robins == "Moderate", method = "REML")
res_nc  <- rma(yi = log(hr), sei = se, data = studies, subset = robins != "Critical", method = "REML")


## -----------------------------------------------------------------------------
#| label: fig-statins
#| fig-cap: "Ten simulated observational studies of statins and all-cause mortality, coloured by ROBINS-I rating, with pooled estimates for all studies and for moderate-risk studies only."
#| fig-width: 8
#| fig-height: 5.5
rob_col <- c(Moderate = "#2E7D32", Serious = "#E07B00", Critical = "#B22222")
fp <- forest(res, atransf = exp, refline = 0, at = log(c(0.5, 0.7, 1, 1.4)), colout = rob_col[studies$robins],
             ilab = as.character(studies$robins), ylim = c(-3, res$k + 3), mlab = "All studies (REML)",
             header = c("Study (simulated)", "HR [95% CI]"), xlab = "Hazard ratio, all-cause mortality")
text(fp$ilab.xpos, res$k + 2, "ROBINS-I", font = 2, cex = fp$cex)
addpoly(res_mod, row = -2, atransf = exp, col = rob_col["Moderate"], mlab = "Moderate risk of bias only")


## -----------------------------------------------------------------------------
#| label: statin-trend
studies$rob <- as.numeric(studies$robins)                       # 1 Moderate, 2 Serious, 3 Critical
mr <- rma(yi = log(hr), sei = se, mods = ~ rob, data = studies, method = "REML", test = "knha")
sub_hr <- sapply(levels(studies$robins), \(r)
  exp(coef(rma(yi = log(hr), sei = se, data = studies, subset = robins == r, method = "REML"))))


## -----------------------------------------------------------------------------
#| label: fig-cutoffs
#| fig-cap: "A hypothetical J-shaped dose-response curve (navy), and the 'exposed vs unexposed' risk ratio obtained when four studies define exposure with different cut-offs (coloured lines)."
#| fig-width: 8
#| fig-height: 3.8
f_true <- \(d) exp(-0.04 * d + 0.0012 * d^2)
thr <- c(10, 20, 30, 40)
binary <- sapply(thr, \(t) mean(f_true(seq(t, 50, 5))))
ggplot(tibble::tibble(d = seq(0, 50, 0.5)), aes(d, f_true(d))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#041E42", linewidth = 1.2) +
  geom_segment(data = tibble::tibble(thr, binary), aes(x = thr, xend = 50, y = binary, yend = binary,
                                                       colour = factor(thr)), linewidth = 1) +
  geom_text(data = tibble::tibble(thr, binary), aes(x = 51, y = binary, colour = factor(thr),
            label = sprintf("exposed = >= %d g/day: RR %.2f", thr, binary)), hjust = 0, size = 3.2) +
  scale_colour_manual(values = c("#B22222", "#E07B00", "#2E7D32", "#1565C0"), guide = "none") +
  scale_x_continuous(limits = c(0, 85), breaks = seq(0, 50, 10)) +
  labs(x = "Dose (g/day)", y = "Risk ratio")


## -----------------------------------------------------------------------------
#| label: alcohol-data
data(alcohol_cvd)
alcohol_cvd |> filter(id <= 2) |> select(author, type, dose, cases, n, logrr)


## -----------------------------------------------------------------------------
#| label: tbl-alcohol-units
#| tbl-cap: "Spline-based RR at 20 g/day vs non-drinkers, with the wine and beer analyses of the same studies treated in three ways."
knots <- quantile(alcohol_cvd$dose, c(0.1, 0.5, 0.9))
rr20 <- function(d) {
  m <- dosresmeta(logrr ~ rcs(dose, knots), id = id, type = type, cases = cases, n = n, se = se,
                  data = d, method = "reml")
  p <- predict(m, newdata = data.frame(dose = c(0, 20)), xref = 0, expo = TRUE)[2, ]
  tibble::tibble(Studies = length(unique(d$id)), `RR at 20 g/day (95% CI)` = sprintf("%.2f (%.2f to %.2f)", p$pred, p$ci.lb, p$ci.ub))
}
bind_rows(`All six, as if independent` = rr20(alcohol_cvd),
          `One per study: wine` = rr20(subset(alcohol_cvd, !id %in% c(4, 6))),
          `One per study: beer` = rr20(subset(alcohol_cvd, !id %in% c(3, 5))), .id = "Data used") |>
  kable(align = "lcc")


## -----------------------------------------------------------------------------
#| label: alcohol-fits
fit_dr <- \(f) dosresmeta(f, id = id, type = type, cases = cases, n = n, se = se, data = alcohol_cvd, method = "reml")
m_lin  <- fit_dr(logrr ~ dose)
m_quad <- fit_dr(logrr ~ dose + I(dose^2))
m_rcs  <- fit_dr(logrr ~ rcs(dose, knots))
grid   <- data.frame(dose = seq(0, 75, 0.5))
curve_of <- function(m, label) {
  p <- predict(m, newdata = grid, xref = 0, expo = TRUE)
  tibble::tibble(dose = grid$dose, est = p$pred, lo = p$ci.lb, hi = p$ci.ub, model = label)
}
p_quad <- waldtest(Sigma = vcov(m_quad), b = coef(m_quad), Terms = 2)$chitest["P"]
p_rcs  <- waldtest(Sigma = vcov(m_rcs), b = coef(m_rcs), Terms = 2)$chitest["P"]


## -----------------------------------------------------------------------------
#| label: alcohol-gp
nr <- subset(alcohol_cvd, dose > 0); y <- nr$logrr; dd <- nr$dose
V  <- as.matrix(Matrix::bdiag(m_lin$Slist))                     # within-study (G-L) covariances
kfun <- \(a, b, s, l) s^2 * exp(-outer(a, b, "-")^2 / (2 * l^2))
kt   <- \(a, b, s, l) kfun(a, b, s, l) - kfun(a, 0, s, l) %*% t(rep(1, length(b))) -
                      rep(1, length(a)) %*% kfun(0, b, s, l) + s^2   # covariance of f(d) - f(0)
slope_cov <- outer(nr$id, nr$id, "==") * outer(dd, dd)          # random slope per study
hyp <- expand.grid(s = seq(0.05, 1.5, length.out = 30), l = seq(5, 80, length.out = 30),
                   tau = c(0, 0.005, 0.01, 0.02, 0.03))
logpost <- apply(hyp, 1, \(h) {
  S <- kt(dd, dd, h["s"], h["l"]) + V + h["tau"]^2 * slope_cov
  -0.5 * (as.numeric(determinant(S)$modulus) + sum(y * solve(S, y))) +
    dnorm(h["s"], 0, 0.5, log = TRUE) + dgamma(1 / h["l"], 3, 50, log = TRUE) - 2 * log(h["l"]) +
    dnorm(h["tau"], 0, 0.02, log = TRUE)
})
set.seed(2026)
gp_draws <- t(sapply(sample(nrow(hyp), 2000, replace = TRUE, prob = exp(logpost - max(logpost))), \(i) {
  h  <- unlist(hyp[i, ]); S <- kt(dd, dd, h["s"], h["l"]) + V + h["tau"]^2 * slope_cov
  Ks <- kt(grid$dose, dd, h["s"], h["l"])
  C  <- kt(grid$dose, grid$dose, h["s"], h["l"]) - Ks %*% solve(S, t(Ks))
  e  <- eigen((C + t(C)) / 2, symmetric = TRUE)
  as.numeric(Ks %*% solve(S, y) + e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(nrow(grid))))
}))
gp_q <- apply(exp(gp_draws), 2, quantile, c(0.5, 0.025, 0.975))
gp_curve <- tibble::tibble(dose = grid$dose, est = gp_q[1, ], lo = gp_q[2, ], hi = gp_q[3, ], model = "Gaussian process")


## -----------------------------------------------------------------------------
#| label: fig-alcohol
#| fig-cap: "Alcohol and cardiovascular disease: four pooled dose-response curves (line and 95% interval). Grey points: study-specific RRs."
#| fig-width: 10
#| fig-height: 3.6
#| out-width: "100%"
curves <- bind_rows(curve_of(m_lin, "Linear"), curve_of(m_quad, "Quadratic"),
                    curve_of(m_rcs, "Restricted cubic spline"), gp_curve) |>
  mutate(model = factor(model, levels = c("Linear", "Quadratic", "Restricted cubic spline", "Gaussian process")))
ggplot(curves, aes(dose)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(data = mutate(nr, rr = exp(logrr)), aes(dose, rr), colour = "grey60", size = 1.5, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#041E42", alpha = 0.15) +
  geom_line(aes(y = est), colour = "#041E42", linewidth = 1) +
  facet_wrap(~ model, nrow = 1) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) + coord_cartesian(ylim = c(0.25, 4)) +
  labs(x = "Alcohol (g/day)", y = "RR vs non-drinkers (log scale)")


## -----------------------------------------------------------------------------
#| label: alcohol-summary
#| include: false
nadir <- curves |> slice_min(est, n = 1, by = model)
at60  <- curves |> filter(dose == 60)


## -----------------------------------------------------------------------------
#| label: fig-coffee-slopes
#| fig-cap: "Coffee and all-cause mortality: study-specific linear slopes (RR per additional cup/day) and the random-effects pooled slope with its 95% prediction interval."
#| fig-width: 7.5
#| fig-height: 6.5
data(coffee_mort)
cf_fit <- \(f, d = coffee_mort) dosresmeta(f, id = id, type = type, cases = cases, n = n, se = se, data = d, method = "reml")
cf_lin <- cf_fit(logrr ~ dose)
cf_slopes <- coffee_mort |> distinct(id, author, gender) |> arrange(id) |>
  mutate(bi = unname(cf_lin$bi[, 1]), vi = unlist(cf_lin$Si))      # Si holds variances
cf_re <- rma(yi = bi, vi = vi, data = cf_slopes, method = "REML", slab = paste0(author, " (", gender, ")"))
forest(cf_re, atransf = exp, at = log(c(0.6, 0.8, 1, 1.2)), addpred = TRUE, cex = 0.75,
       header = c("Study (sex)", "RR per cup/day [95% CI]"), xlab = "RR per additional cup/day")
cf_pr <- predict(cf_re, transf = exp)


## -----------------------------------------------------------------------------
#| label: fig-coffee-curves
#| fig-cap: "Coffee and all-cause mortality: linear and restricted cubic spline dose-response curves (knots at the quartiles of dose)."
#| fig-width: 8
#| fig-height: 3.6
cf_knots <- quantile(coffee_mort$dose, c(0.25, 0.5, 0.75))
cf_rcs   <- cf_fit(logrr ~ rcs(dose, cf_knots))
cf_kl    <- cf_fit(logrr ~ rcs(dose, cf_knots), subset(coffee_mort, id == 3))   # one study alone
cf_wald  <- \(m) waldtest(Sigma = vcov(m), b = coef(m), Terms = 2)$chitest
cf_grid  <- data.frame(dose = seq(0, 8, 0.1))
cf_curve <- \(m, label) { p <- predict(m, newdata = cf_grid, xref = 0, expo = TRUE)
  tibble::tibble(dose = cf_grid$dose, est = p$pred, lo = p$ci.lb, hi = p$ci.ub, model = label) }
cfc <- bind_rows(cf_curve(cf_lin, "Linear"), cf_curve(cf_rcs, "Restricted cubic spline"))
ggplot(cfc, aes(dose)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#041E42", alpha = 0.15) +
  geom_line(aes(y = est), colour = "#041E42", linewidth = 1) +
  facet_wrap(~ model) + scale_y_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1, 1.1)) +
  coord_cartesian(ylim = c(0.6, 1.15)) + labs(x = "Coffee (cups/day)", y = "RR of death vs 0 cups/day")
cf_at <- \(m, d) with(filter(cfc, model == m, abs(dose - d) < 1e-8), sprintf("%.2f (%.2f to %.2f)", est, lo, hi))

