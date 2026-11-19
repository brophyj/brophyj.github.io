# R code for chapter: 14-missing-data.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch14
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
  library(nlme)
  library(mice)
})


## -----------------------------------------------------------------------------
#| label: fig-mechanisms
#| fig-cap: "The same 2000 simulated patients under three missingness mechanisms. Bars: observed (navy) and missing (red) outcomes, stacked. Lines: the true mean of everyone (black), the mean of the observed (orange, dotted) and an estimate adjusted for the observed baseline value (green, dashed)."
#| fig-width: 10
#| fig-height: 4
#| out-width: "100%"
set.seed(42)
n <- 2000
baseline <- rnorm(n)                                        # observed for everyone
y <- 0.7 * baseline + rnorm(n, 0, sqrt(1 - 0.7^2))          # final outcome (higher = worse)
miss <- list(MCAR = rbinom(n, 1, 0.30),                     # a coin flip
             MAR  = rbinom(n, 1, plogis(-1 + 1.5 * baseline)),  # depends on observed baseline
             MNAR = rbinom(n, 1, plogis(-1 + 1.5 * y)))         # depends on the outcome itself
est <- imap_dfr(miss, \(m, nm) {
  obs <- m == 0
  fit <- lm(y ~ baseline, subset = obs)
  tibble::tibble(Mechanism = nm, pct = 100 * mean(m),
                 `True mean (everyone)` = mean(y), `Mean of observed only` = mean(y[obs]),
                 `Baseline-adjusted estimate` = mean(predict(fit, newdata = data.frame(baseline))))
}) |> mutate(Mechanism = factor(Mechanism, c("MCAR", "MAR", "MNAR")))
lines <- est |> select(-pct) |>
  pivot_longer(-Mechanism, names_to = "Estimate", values_to = "x") |>
  mutate(Estimate = factor(Estimate, c("True mean (everyone)", "Mean of observed only",
                                       "Baseline-adjusted estimate")))
lab <- est |> mutate(label = sprintf("%.0f%% missing\nbias, observed: %+.2f\nbias, adjusted: %+.2f", pct,
                                     `Mean of observed only` - `True mean (everyone)`,
                                     `Baseline-adjusted estimate` - `True mean (everyone)`))
imap_dfr(miss, \(m, nm) tibble::tibble(Mechanism = nm, y = y,
                                       Status = ifelse(m == 1, "Missing", "Observed"))) |>
  mutate(Status = factor(Status, c("Observed", "Missing")),
         Mechanism = factor(Mechanism, c("MCAR", "MAR", "MNAR"))) |>
  ggplot(aes(y)) +
  geom_histogram(aes(fill = Status), bins = 35, colour = "white") +
  geom_vline(data = lines, aes(xintercept = x, colour = Estimate, linetype = Estimate), linewidth = 1.1) +
  geom_label(data = lab, aes(x = -2.5, y = Inf, label = label), inherit.aes = FALSE,
             hjust = 0, vjust = 1.1, size = 3, fill = "white", label.size = 0, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  facet_wrap(~ Mechanism) +
  scale_fill_manual(values = c(Observed = "#041E42", Missing = "#B22222"), name = NULL) +
  scale_colour_manual(values = c("black", "#E69F00", "#2E7D32"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = NULL) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(x = "Final outcome (higher = worse)", y = "Patients") +
  theme(legend.position = "bottom", legend.box = "vertical")


## -----------------------------------------------------------------------------
#| label: missing-sim
sim_trial <- function(n = 150, mech = "MCAR", effect = -6) {
  trt <- rep(0:1, each = n)
  u   <- rnorm(2 * n, 0, 4)                                   # patient level
  s   <- rnorm(2 * n, 0, 6)                                   # patient response
  Y   <- sapply(0:3, \(t) 30 + u + (-3 + trt * effect / 3 + s / 3) * t + rnorm(2 * n, 0, 4))
  R   <- matrix(TRUE, 2 * n, 4)                               # TRUE = observed
  for (t in 2:4) {
    worse_prev <- Y[, t - 1] - Y[, 1] + 3 * (t - 2)           # lack of improvement, previous visit
    worse_now  <- Y[, t]     - Y[, 1] + 3 * (t - 1)           # ... current visit (unseen if missing)
    lp <- switch(mech, MCAR = rep(-2.2, 2 * n),               # plogis(-2.2) = 10% per visit
                 MAR = -2.2 + 0.3 * worse_prev, MNAR = -2.2 + 0.3 * worse_now)
    R[, t] <- R[, t - 1] & rbinom(2 * n, 1, plogis(lp)) == 0  # monotone dropout
  }
  Yobs <- Y; Yobs[!R] <- NA
  list(Y = Y, Yobs = Yobs, trt = trt)
}
coefse <- function(fit, term = "trt") {
  s <- if (inherits(fit, "gls")) summary(fit)$tTable else summary(fit)$coefficients
  c(est = s[term, 1], se = s[term, 2])
}
analyse <- function(d) {
  Yo <- d$Yobs; trt <- d$trt; y0 <- Yo[, 1]; y3 <- Yo[, 4]
  last <- apply(Yo, 1, \(r) tail(na.omit(r), 1))
  out <- list(
    `No missing data` = coefse(lm(d$Y[, 4] ~ trt + y0)),     # benchmark, before dropout
    `Complete case`   = coefse(lm(y3 ~ trt + y0)),
    LOCF              = coefse(lm(last ~ trt + y0)),
    BOCF              = coefse(lm(ifelse(is.na(y3), y0, y3) ~ trt + y0)),
    `Mean imputation` = coefse(lm(ifelse(is.na(y3), ave(y3, trt, FUN = \(x) mean(x, na.rm = TRUE)), y3) ~ trt + y0)))
  imp <- mice(data.frame(trt, y0, y1 = Yo[, 2], y2 = Yo[, 3], y3), m = 5, method = "norm", printFlag = FALSE)
  pm  <- summary(pool(with(imp, lm(y3 ~ trt + y0))))
  out$`Multiple imputation` <- c(est = pm$estimate[2], se = pm$std.error[2])
  N    <- length(trt)
  long <- tibble::tibble(id = rep(1:N, 3), trt = rep(trt, 3), y0 = rep(y0, 3),
                         visit = factor(rep(1:3, each = N)), y = c(Yo[, 2:4])) |> filter(!is.na(y))
  mm <- gls(y ~ visit * trt + visit * y0, data = long,           # mixed model for repeated measures
            correlation = corSymm(form = ~ as.integer(visit) | id),
            weights = varIdent(form = ~ 1 | visit))
  V  <- vcov(mm)[c("trt", "visit3:trt"), c("trt", "visit3:trt")]
  out$`Mixed model (MMRM)` <- c(est = sum(coef(mm)[c("trt", "visit3:trt")]), se = sqrt(sum(V)))
  bind_rows(out, .id = "method")
}
set.seed(629)
sims <- map_dfr(c("MCAR", "MAR", "MNAR"), \(m) map_dfr(1:150, \(r) {
  d <- sim_trial(mech = m)
  analyse(d) |> mutate(mech = m, miss_pbo = mean(is.na(d$Yobs[d$trt == 0, 4])),
                       miss_act = mean(is.na(d$Yobs[d$trt == 1, 4])))
}))
sim_summary <- sims |>
  mutate(mech = factor(mech, c("MCAR", "MAR", "MNAR")),
         method = factor(method, c("No missing data", "Complete case", "LOCF", "BOCF",
                                   "Mean imputation", "Multiple imputation", "Mixed model (MMRM)"))) |>
  group_by(mech, method) |>
  summarise(est_mean = mean(est), lo = mean(est - 1.96 * se), hi = mean(est + 1.96 * se),
            coverage = mean(abs(est + 6) < 1.96 * se),
            miss_pbo = mean(miss_pbo), miss_act = mean(miss_act), .groups = "drop")
miss_rates <- sim_summary |> distinct(mech, miss_pbo, miss_act)


## -----------------------------------------------------------------------------
#| label: fig-missing-sim
#| fig-cap: "Six analyses of 150 simulated trials under each mechanism (true effect −6, dashed line). Point: average estimate; bar: a typical single trial's 95% CI. Percentages at right: coverage, the proportion of trials whose 95% CI contains −6."
#| fig-width: 9
#| fig-height: 6.5
sim_summary |>
  mutate(method = forcats::fct_rev(method), cov_lab = sprintf("%.0f%%", 100 * coverage)) |>
  ggplot(aes(est_mean, method, colour = mech)) +
  geom_vline(xintercept = -6, linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi, group = forcats::fct_rev(mech)),
                  position = position_dodge(width = 0.7)) +
  geom_text(aes(x = 1.8, label = cov_lab, group = forcats::fct_rev(mech)), size = 3.3,
            position = position_dodge(width = 0.7), show.legend = FALSE) +
  coord_cartesian(xlim = c(-9, 2.3)) +
  scale_colour_manual(values = c(MCAR = "#2E7D32", MAR = "#E69F00", MNAR = "#B22222"), name = NULL) +
  labs(x = "Estimated treatment effect (true = −6)", y = NULL) +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: fig-locf
#| fig-cap: "A single patient who drops out at week 4. LOCF freezes the outcome at the last observed value; the true trajectory continues to change."
#| fig-width: 7
#| fig-height: 3.8
set.seed(99)
true_traj <- function(t) 50 + 8 * t - 0.5 * t^2
obs4 <- true_traj(0:4) + rnorm(5, 0, 4)
tibble::tibble(time = c(0:4, 4:12, 0:12),
               value = c(obs4, rep(last(obs4), 9), true_traj(0:12) + rnorm(13, 0, 4)),
               type = c(rep("Observed (dropout at week 4)", 5), rep("LOCF imputation", 9),
                        rep("True trajectory (if followed)", 13))) |>
  ggplot(aes(time, value, colour = type, linetype = type)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 4, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = c("#B22222", "#041E42", "#E69F00"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid", "solid"), name = NULL) +
  labs(x = "Week", y = "Outcome score") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: worst-best
rr <- function(et, nt, ec, nc) {
  e <- escalc("RR", ai = et, n1i = nt, ci = ec, n2i = nc)
  x <- exp(e$yi + c(0, -1.96, 1.96) * sqrt(e$vi))
  sprintf("%.2f (%.2f to %.2f)", x[1], x[2], x[3])
}
tibble::tibble(
  Scenario = c("Truth: everyone counted", "Complete case", "Missing = no event",
               "Worst case for treatment", "Best case for treatment"),
  `RR (95% CI)` = c(rr(12 + 4, 63, 22 + 2, 60), rr(12, 55, 22, 55), rr(12, 63, 22, 60),
                    rr(12 + 8, 63, 22, 60), rr(12, 63, 22 + 5, 60))) |>
  kable(align = "lc")


## -----------------------------------------------------------------------------
#| label: p-backcalc
md <- 5.0; p <- 0.03; n1 <- n2 <- 40                    # "MD 5.0 points, p = 0.03, 40 per arm"
t_val <- qt(1 - p / 2, df = n1 + n2 - 2)
se_md <- md / t_val
round(c(t = t_val, SE = se_md, SD = se_md / sqrt(1 / n1 + 1 / n2)), 2)


## -----------------------------------------------------------------------------
#| label: tbl-r-sensitivity
#| tbl-cap: "Change-score SMD for one trial, by the assumed baseline–final correlation r. The numerator (difference in mean change, 8.4 points) is the same in every row."
sd_change <- \(sb, sf, r) sqrt(sb^2 + sf^2 - 2 * r * sb * sf)
diff_change <- (28.9 - 17.9) - (28.4 - 25.8)
tibble::tibble(r = c(0, 0.25, 0.5, 0.75, 0.9, 0.97)) |>
  mutate(`SD change, fluoxetine` = sd_change(8.4, 8.1, r),
         `SD change, nortriptyline` = sd_change(8.7, 8.6, r),
         `Pooled SD` = sqrt((`SD change, fluoxetine`^2 + `SD change, nortriptyline`^2) / 2),
         `Change-score SMD` = diff_change / `Pooled SD`) |>
  kable(digits = 2, align = "ccccc")


## -----------------------------------------------------------------------------
#| label: tbl-wan
#| tbl-cap: "Mean and SD estimated from quartiles (Wan et al.) for three hypothetical studies of 6-minute walk distance (m)."
eta    <- \(n) 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))     # expected IQR of n standard normals
wan_s2 <- \(q1, m, q3, n) c(mean = (q1 + m + q3) / 3, sd = (q3 - q1) / eta(n))
tibble::tribble(~study, ~n, ~q1, ~med, ~q3,
                "Study 1", 45, 285, 320, 360,
                "Study 2", 112, 260, 298, 345,
                "Study 3", 28, 300, 335, 375) |>
  mutate(`η(n)` = eta(n), est = pmap(list(q1, med, q3, n), wan_s2),
         `Estimated mean` = map_dbl(est, "mean"), `Estimated SD` = map_dbl(est, "sd")) |>
  select(-est) |>
  kable(digits = 2)


## -----------------------------------------------------------------------------
#| label: metamedian
#| eval: !expr requireNamespace("metamedian", quietly = TRUE)
library(metamedian)     # dat.age: 51 COVID-19 studies, age of non-survivors vs survivors
set.seed(629)
fits <- list(
  `Means only (22 studies)` = rma(measure = "MD", data = dat.age, m1i = mean.g1, sd1i = sd.g1, n1i = n.g1,
                                  m2i = mean.g2, sd2i = sd.g2, n2i = n.g2),
  `Hozo / Wan / Bland`  = metamean(dat.age, mean_method = "hozo/wan/bland"),
  `Luo mean, Wan SD`    = metamean(dat.age, mean_method = "luo", sd_method = "wan"),
  `Quantile estimation` = metamean(dat.age, mean_method = "qe", nboot = 100))
imap_dfr(fits, \(f, nm) tibble::tibble(Method = nm, k = f$k,
  `Age difference, years (95% CI)` = sprintf("%.1f (%.1f to %.1f)", f$b, f$ci.lb, f$ci.ub))) |>
  kable()


## -----------------------------------------------------------------------------
#| label: haloperidol-data
haloperidol <- tibble::tribble(   # improved, not improved, missing; haloperidol (H) and placebo (P)
  ~study, ~rH, ~fH, ~mH, ~rP, ~fP, ~mP,
  "Arvanitis 1997", 25, 25, 2, 18, 33, 0,   "Beasley 1996", 29, 18, 22, 20, 14, 34,
  "Bechelli 1983", 12, 17, 1,  2, 28, 1,    "Borison 1992",  3,  9,  0,  0, 12,  0,
  "Chouinard 1993", 10, 11, 0, 3, 19, 0,    "Durost 1964",  11,  8,  0,  1, 14,  0,
  "Garry 1962",  7, 18, 1,  4, 21, 1,       "Howard 1974",   8,  9,  0,  3, 10,  0,
  "Marder 1994", 19, 45, 2, 14, 50, 2,      "Nishikawa 1982", 1, 9,  0,  0, 10,  0,
  "Nishikawa 1984", 11, 23, 3, 0, 13, 0,    "Reschke 1974", 20,  9,  0,  2,  9,  0,
  "Selman 1976", 17,  1, 11, 7,  4, 18,     "Serafetinides 1972", 4, 10, 0, 0, 13, 1,
  "Simpson 1967",  2, 14, 0, 0,  7, 1,      "Spencer 1992", 11,  1,  0,  1, 11,  0,
  "Vichaiya 1971", 9, 20, 1, 0, 29, 1)
aca <- rma(measure = "RR", ai = rH, bi = fH, ci = rP, di = fP, data = haloperidol,
           method = "EE", add = 0.5, to = "only0")          # complete case, common effect
fmt_rr <- \(f) sprintf("%.2f (%.2f to %.2f)", exp(f$b[1]), exp(f$ci.lb), exp(f$ci.ub))
haloperidol |>
  mutate(Haloperidol = sprintf("%d/%d (%d)", rH, rH + fH, mH),
         Placebo = sprintf("%d/%d (%d)", rP, rP + fP, mP),
         `Complete-case weight` = sprintf("%.0f%%", weights(aca))) |>
  select(Study = study, Haloperidol, Placebo, `Complete-case weight`) |>
  kable(align = "lccc", caption = "Haloperidol vs placebo: improved / analysed (missing), and each trial's weight in the complete-case meta-analysis.")


## -----------------------------------------------------------------------------
#| label: imor-functions
adj_arm <- function(r, f, m, imor, sd_log = 0, cc = 0) {    # one arm (Higgins, White & Wood 2008)
  N  <- r + f + m; pi <- m / N                               # randomised; proportion missing
  n  <- r + f + 2 * cc; p <- (r + cc) / n                    # analysed; observed risk
  pm <- imor * p / (imor * p + 1 - p)                        # risk among the missing
  ps <- (1 - pi) * p + pi * pm                               # risk among everyone randomised
  dp <- (1 - pi) + pi * imor / (imor * p + 1 - p)^2          # d ps / d p
  v  <- dp^2 * p * (1 - p) / n +                             # 1. sampling error in p
        (pm - p)^2 * pi * (1 - pi) / N +                     # 2. sampling error in pi
        (pi * pm * (1 - pm))^2 * sd_log^2                    # 3. uncertainty in log IMOR
  tibble::tibble(ps, v)
}
imor_ma <- function(d, imorH, imorP = imorH, sd_log = 0) {   # pooled log RR, common effect
  cc <- ifelse(with(d, rH * fH * rP * fP) == 0, 0.5, 0)
  h  <- with(d, adj_arm(rH, fH, mH, imorH, sd_log, cc))
  pl <- with(d, adj_arm(rP, fP, mP, imorP, sd_log, cc))
  rma(yi = log(h$ps / pl$ps), vi = h$v / h$ps^2 + pl$v / pl$ps^2, method = "EE")
}
scen <- list(`Complete case (IMOR = 1)` = imor_ma(haloperidol, 1),
             `Missing did worse (IMOR = 0.5)` = imor_ma(haloperidol, 0.5),
             `All missing failed (IMOR = 0)` = imor_ma(haloperidol, 0),
             `Missing did better (IMOR = 2)` = imor_ma(haloperidol, 2),
             `Worst case for haloperidol` = imor_ma(haloperidol, 0, 1e6),
             `Best case for haloperidol` = imor_ma(haloperidol, 1e6, 0))
imap_dfr(scen, \(f, nm) tibble::tibble(Assumption = nm, `RR (95% CI)` = fmt_rr(f),
                                       `Weight: Beasley` = sprintf("%.0f%%", weights(f)[2]),
                                       `Weight: Selman` = sprintf("%.0f%%", weights(f)[13]))) |>
  kable(align = "lccc", caption = "IMOR-adjusted meta-analyses of haloperidol vs placebo.")


## -----------------------------------------------------------------------------
#| label: metamiss
#| eval: !expr requireNamespace("metasens", quietly = TRUE)
mb <- meta::metabin(rH, rH + fH, rP, rP + fP, data = haloperidol, studlab = study,
                    sm = "RR", method = "Inverse", common = TRUE, random = FALSE)
mm <- metasens::metamiss(mb, miss.e = mH, miss.c = mP, IMOR.e = 0.5, IMOR.c = 0.5)
sprintf("%.2f (%.2f to %.2f)", exp(mm$TE.common), exp(mm$lower.common), exp(mm$upper.common))


## -----------------------------------------------------------------------------
#| label: fig-tipping
#| fig-cap: "Pooled RR over IMORs from 0.05 to 20 in each arm. Inside the dashed contour (top left) the 95% CI includes 1. The diamond marks complete-case analysis."
#| fig-width: 6.5
#| fig-height: 5.5
imor_brk <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20)
grid_df <- expand.grid(dH = seq(-3, 3, by = 0.2), dP = seq(-3, 3, by = 0.2)) |>
  mutate(fit = map2(dH, dP, \(a, b) imor_ma(haloperidol, exp(a), exp(b))),
         RR = map_dbl(fit, \(f) exp(f$b[1])), lo = map_dbl(fit, \(f) exp(f$ci.lb))) |>
  select(-fit)
ggplot(grid_df, aes(dH, dP)) +
  geom_raster(aes(fill = RR)) +
  geom_contour(aes(z = lo), breaks = 1, colour = "white", linewidth = 1.2, linetype = "dashed") +
  annotate("point", x = 0, y = 0, shape = 18, size = 5, colour = "white") +
  scale_fill_viridis_c(name = "Pooled RR", trans = "log", breaks = c(1, 1.25, 1.5, 2)) +
  scale_x_continuous(breaks = log(imor_brk), labels = imor_brk) +
  scale_y_continuous(breaks = log(imor_brk), labels = imor_brk) +
  labs(x = "IMOR, haloperidol arm (log scale)", y = "IMOR, placebo arm (log scale)") +
  coord_fixed()


## -----------------------------------------------------------------------------
#| label: tbl-imor-uncertain
#| tbl-cap: "Haloperidol vs placebo with an uncertain IMOR (log IMOR ~ N(−0.35, 0.35²) in each arm)."
set.seed(629)
draws <- map_dbl(1:2000, \(i) {
  f <- imor_ma(haloperidol, exp(rnorm(1, -0.35, 0.35)), exp(rnorm(1, -0.35, 0.35)))
  rnorm(1, f$b[1], f$se)
})
an <- imor_ma(haloperidol, exp(-0.35), sd_log = 0.35)
tibble::tibble(Analysis = c("Complete case", "Uncertain IMOR, analytic", "Uncertain IMOR, simulation"),
               `Pooled RR (95% interval)` = c(fmt_rr(aca), fmt_rr(an),
                 sprintf("%.2f (%.2f to %.2f)", exp(median(draws)), exp(quantile(draws, 0.025)),
                         exp(quantile(draws, 0.975)))),
               `Weight: Beasley` = sprintf("%.0f%%", c(weights(aca)[2], weights(an)[2], NA))) |>
  mutate(`Weight: Beasley` = replace(`Weight: Beasley`, 3, "")) |>
  kable(align = "lcc")

