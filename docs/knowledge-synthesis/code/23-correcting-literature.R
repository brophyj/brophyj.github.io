# R code for chapter: 23-correcting-literature.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch23
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
  library(patchwork)
})
f2   <- function(x) sprintf("%.2f", x)
f3   <- function(x) sprintf("%.3f", x)
ci_t <- function(est, lo, hi) sprintf("%.2f (%.2f to %.2f)", est, lo, hi)
pct  <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")
big  <- function(x) formatC(signif(x, 2), format = "d", big.mark = ",")
navy <- "#041E42"; red <- "#C8102E"; green <- "#2E7D32"


## -----------------------------------------------------------------------------
#| label: audit-data
dat <- data.frame(
  study = c("Zhang 2015", "PRAGUE-13", "FIRE", "PRAMI", "DANAMI-3-PRIMULTI",
            "Politi", "CvLPRIT", "COMPARE-ACUTE", "COMPLETE", "Hamza", "Ghani",
            "Zhang 2023", "FULL REVASC", "Zhao", "HELP-AMI"),
  year  = c(2015, 2015, 2023, 2013, 2015, 2010, 2019, 2017, 2019, 2016, 2012,
            2023, 2024, 2016, 2004),
  ai    = c(13,  6, 66, 12, 15, 10,  9,  4,   96, 1, 4, 3, 76, 3, 1),   # deaths, complete
  n1i   = c(215, 106, 720, 234, 313, 130, 150, 295, 2016, 50, 79, 115, 764, 148, 52),
  ci    = c(15,  7, 93, 16, 11, 13, 15, 10,  106, 4, 0, 4, 72, 4, 0),  # deaths, culprit-only
  n2i   = c(213, 108, 725, 231, 314, 84, 146, 590, 2025, 50, 40, 114, 778, 158, 17))
# log RRs; 0.5 added only to the two trials with a zero cell, as the authors did
es <- escalc(measure = "RR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
             data = dat, add = 0.5, to = "only0", slab = study)
es$sei <- sqrt(es$vi)
m_z     <- rma(yi, vi, data = es, method = "REML")                    # Wald z
m_kh    <- rma(yi, vi, data = es, method = "REML", test = "knha")     # HKSJ
m_adhoc <- rma(yi, vi, data = es, method = "REML", test = "adhoc")    # truncated HKSJ
p_z <- predict(m_z, transf = exp); p_kh <- predict(m_kh, transf = exp)
p_ad <- predict(m_adhoc, transf = exp)
w <- setNames(weights(m_z), dat$study)
big3 <- c("FIRE", "COMPLETE", "FULL REVASC")


## -----------------------------------------------------------------------------
#| label: fig-audit-forest
#| fig-cap: "All-cause mortality, complete vs culprit-only revascularization: the 15 trials, the pooled estimate reproduced with the authors' method (REML with HKSJ, navy diamond), and the correctly computed 95% prediction interval (red)."
#| fig-width: 11
#| fig-height: 6.4
#| out-width: "100%"
fd <- as.data.frame(es) |>
  mutate(rr = exp(yi), lo = exp(yi - 1.96 * sei), hi = exp(yi + 1.96 * sei), wt = w,
         y = rev(seq_len(n())) + 2, counts = sprintf("%d/%d", ai, n1i), counts0 = sprintf("%d/%d", ci, n2i),
         rr_txt = sprintf("%.2f (%.2f-%.2f)", rr, lo, hi))
pool <- tibble::tibble(y = c(1, 0), lab = c("Pooled (HKSJ, as published)", "Prediction interval (correct)"),
                       est = p_kh$pred, lo = c(p_kh$ci.lb, p_kh$pi.lb), hi = c(p_kh$ci.ub, p_kh$pi.ub),
                       col = c(navy, red)) |>
  mutate(rr_txt = sprintf("%.2f (%.2f-%.2f)", est, lo, hi))
yl <- c(-0.6, max(fd$y) + 1.2)
left <- ggplot() +
  geom_text(data = fd, aes(0, y, label = study), hjust = 0, size = 3.8) +
  geom_text(data = fd, aes(0.55, y, label = counts), hjust = 0, size = 3.6) +
  geom_text(data = fd, aes(0.80, y, label = counts0), hjust = 0, size = 3.6) +
  geom_text(data = pool, aes(0, y, label = lab), hjust = 0, size = 3.8, fontface = "bold", colour = pool$col) +
  annotate("text", x = c(0, 0.55, 0.80), y = max(fd$y) + 1.1, label = c("Trial", "Complete", "Culprit"),
           hjust = 0, fontface = "bold", size = 3.8) +
  scale_y_continuous(limits = yl) + scale_x_continuous(limits = c(0, 1)) + theme_void()
mid <- ggplot(fd, aes(y = y)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_vline(xintercept = p_kh$pred, colour = navy, linetype = "dashed") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey30") +
  geom_point(aes(x = rr, size = wt), shape = 15, colour = navy) +
  geom_polygon(data = tibble::tibble(x = c(p_kh$ci.lb, p_kh$pred, p_kh$ci.ub, p_kh$pred), y = c(1, 1.35, 1, 0.65)),
               aes(x = x, y = y), fill = navy) +
  geom_linerange(data = pool[2, ], aes(y = y, xmin = lo, xmax = hi), colour = red, linewidth = 1.2) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10), labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10"),
                oob = scales::squish) +
  coord_cartesian(xlim = c(0.08, 12)) + scale_size_area(max_size = 6, guide = "none") +
  scale_y_continuous(limits = yl, breaks = NULL) +
  labs(x = "Relative risk of death (log scale)", y = NULL) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
right <- ggplot() +
  geom_text(data = fd, aes(0, y, label = rr_txt), hjust = 0, size = 3.6) +
  geom_text(data = fd, aes(1, y, label = sprintf("%.1f%%", wt)), hjust = 1, size = 3.6) +
  geom_text(data = pool, aes(0, y, label = rr_txt), hjust = 0, size = 3.6, fontface = "bold", colour = pool$col) +
  annotate("text", x = c(0, 1), y = max(fd$y) + 1.1, label = c("RR (95% CI)", "Weight"),
           hjust = c(0, 1), fontface = "bold", size = 3.8) +
  scale_y_continuous(limits = yl) + scale_x_continuous(limits = c(0, 1)) + theme_void()
left + mid + right + plot_layout(widths = c(1.35, 1.1, 0.75))


## -----------------------------------------------------------------------------
#| label: tbl-audit-ci
#| tbl-cap: "The same REML model with three ways of building the confidence interval."
tibble::tibble(`CI method` = c("Wald z (τ treated as known)", "HKSJ", "Truncated HKSJ (q ≥ 1)"),
               `RR (95% CI)` = c(ci_t(p_z$pred, p_z$ci.lb, p_z$ci.ub), ci_t(p_kh$pred, p_kh$ci.lb, p_kh$ci.ub),
                                 ci_t(p_ad$pred, p_ad$ci.lb, p_ad$ci.ub)),
               P = f3(c(m_z$pval, m_kh$pval, m_adhoc$pval)),
               q = c("not used", f2(m_kh$s2w), sprintf("%s (truncated)", f2(m_adhoc$s2w))),
               ` ` = c("", "matches the paper", "")) |>
  kable(align = "lccll")


## -----------------------------------------------------------------------------
#| label: fig-audit-many
#| fig-cap: "The pooled relative risk of death and its 95% CI under six defensible analyses of the same 15 trials. The published analysis is in red."
#| fig-width: 9
#| fig-height: 4.2
m_ee   <- rma(yi, vi, data = es, method = "EE")
m_mh   <- rma.mh(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat, measure = "RR")
m_glmm <- rma.glmm(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat, model = "UM.FS")
many <- tibble::tibble(
  method = c("Random effects, HKSJ (published)", "Random effects, Wald z", "Random effects, truncated HKSJ",
             "Common effect, inverse variance", "Common effect, Mantel-Haenszel", "Binomial GLMM (odds ratio)"),
  est = exp(c(coef(m_kh), coef(m_z), coef(m_adhoc), coef(m_ee), coef(m_mh), coef(m_glmm))),
  lo  = exp(c(m_kh$ci.lb, m_z$ci.lb, m_adhoc$ci.lb, m_ee$ci.lb, m_mh$ci.lb, m_glmm$ci.lb)),
  hi  = exp(c(m_kh$ci.ub, m_z$ci.ub, m_adhoc$ci.ub, m_ee$ci.ub, m_mh$ci.ub, m_glmm$ci.ub)),
  p   = c(m_kh$pval, m_z$pval, m_adhoc$pval, m_ee$pval, m_mh$pval, m_glmm$pval)) |>
  mutate(method = forcats::fct_reorder(method, hi), pub = grepl("published", method),
         lab = sprintf("%.2f (%.2f-%.2f)  P = %.3f", est, lo, hi, p))
ggplot(many, aes(est, method, colour = pub)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1) + geom_point(size = 3) +
  geom_text(aes(x = 1.07, label = lab), hjust = 0, size = 3.5) +
  scale_colour_manual(values = c(navy, red), guide = "none") +
  scale_x_continuous(limits = c(0.68, 1.4), breaks = seq(0.7, 1.05, 0.05)) +
  labs(x = "Relative risk of death (95% CI)", y = NULL)


## -----------------------------------------------------------------------------
#| label: audit-tau
tau_ci <- confint(m_z)$random
est_tau2 <- sapply(c(DL = "DL", REML = "REML", PM = "PM", ML = "ML", SJ = "SJ"),
                   \(x) rma(yi, vi, data = es, method = x)$tau2)
tau_hi <- tau_ci["tau", "ci.ub"]
wi <- 1 / es$vi
s2_typ <- (m_z$k - 1) * sum(wi) / (sum(wi)^2 - sum(wi^2))   # "typical" within-trial variance
i2_of  <- \(tau2) tau2 / (tau2 + s2_typ)
round(est_tau2, 3)
round(tau_ci[c("tau", "I^2(%)"), ], 2)


## -----------------------------------------------------------------------------
#| label: tbl-audit-pi
#| tbl-cap: "Confidence and prediction intervals from the authors' model (REML with HKSJ), on the log and RR scales."
k <- m_z$k
pi_t <- exp(coef(m_z) + c(-1, 1) * qt(0.975, k - 2) * sqrt(m_z$tau2 + m_z$se^2))   # Higgins, t(k-2)
tibble::tibble(
  Interval = c("CI", "PI (HKSJ, t with k − 1 df, as metafor)", "PI (Higgins, t with k − 2 df)", "PI as published"),
  `Log scale` = c(sprintf("%.2f to %+.2f", log(p_kh$ci.lb), log(p_kh$ci.ub)),
                  sprintf("%.2f to %+.2f", log(p_kh$pi.lb), log(p_kh$pi.ub)),
                  sprintf("%.2f to %+.2f", log(pi_t[1]), log(pi_t[2])), "-0.39 to -0.08"),
  `RR scale` = c(sprintf("%.2f to %.2f", p_kh$ci.lb, p_kh$ci.ub), sprintf("%.2f to %.2f", p_kh$pi.lb, p_kh$pi.ub),
                 sprintf("%.2f to %.2f", pi_t[1], pi_t[2]), sprintf("%.2f to %.2f", exp(-0.39), exp(-0.08)))) |>
  kable(align = "lcc")


## -----------------------------------------------------------------------------
#| label: audit-influence
rr_of   <- \(s) exp(es$yi[dat$study == s])
no_fire <- predict(rma(yi, vi, data = es, subset = study != "FIRE"), transf = exp)
b3      <- predict(rma(yi, vi, data = es, subset = study %in% big3, method = "EE"), transf = exp)
tibble::tibble(Trial = big3, Population = c("Age ≥ 75, STEMI or NSTEMI", "STEMI", "Mostly STEMI"),
               Guidance = c("Physiology", "Angiography", "FFR"),
               RR = f2(sapply(big3, rr_of)), Weight = pct(w[big3] / 100)) |>
  kable(align = "lllcc")


## -----------------------------------------------------------------------------
#| label: fig-audit-l1o
#| fig-cap: "Pooled RR omitting each trial in turn, with Wald (left) and HKSJ (right, as published) 95% CIs. The point estimates are identical; only the intervals differ. Red: interval excludes 1."
#| fig-width: 11
#| fig-height: 5.2
l1o    <- leave1out(m_z,  transf = exp)        # Wald intervals
l1o_kh <- leave1out(m_kh, transf = exp)        # HKSJ intervals, same omissions
l1o_both <- bind_rows(
  tibble::tibble(method = "Wald z 95% CI", study = dat$study, est = l1o$estimate,
                 lo = l1o$ci.lb, hi = l1o$ci.ub, p = l1o$pval),
  tibble::tibble(method = "HKSJ 95% CI (as published)", study = dat$study, est = l1o_kh$estimate,
                 lo = l1o_kh$ci.lb, hi = l1o_kh$ci.ub, p = l1o_kh$pval)) |>
  mutate(method = forcats::fct_inorder(method),
         study  = forcats::fct_reorder(paste("without", study), est),   # same order in both panels
         sig    = hi < 1)
ggplot(l1o_both, aes(est, study, colour = sig)) +
  geom_vline(xintercept = 1, colour = "grey40") + geom_vline(xintercept = p_z$pred, colour = navy, linetype = 2) +
  geom_linerange(aes(xmin = lo, xmax = hi)) + geom_point(size = 2) +
  facet_wrap(~ method) +
  scale_colour_manual(values = c(`FALSE` = navy, `TRUE` = red), guide = "none") +
  labs(x = "Pooled RR omitting one trial (REML)", y = NULL) +
  theme(strip.text = element_text(face = "bold"))
n_sig <- tapply(l1o_both$p < 0.05, l1o_both$method, sum)


## -----------------------------------------------------------------------------
#| label: fig-audit-cumul
#| fig-cap: "Cumulative meta-analysis in order of publication (REML, HKSJ, as published). Red: P < 0.05."
#| fig-width: 7
#| fig-height: 4.6
o   <- order(dat$year, dat$n1i + dat$n2i)
pos <- integer(length(o)); pos[o] <- seq_along(o)   # each trial's place in publication order
cmk <- cumul(m_kh, order = pos, transf = exp)        # metafor >= 3: 'order' is a sort key, not a permutation
cum_df <- tibble::tibble(step = seq_along(o), study = dat$study[o], year = dat$year[o],
                         est = cmk$estimate, lo = cmk$ci.lb, hi = cmk$ci.ub, p = cmk$pval) |> filter(step >= 4)
cum_df |>
  mutate(lab = factor(paste0("+ ", study, " (", year, ")"), levels = rev(paste0("+ ", study, " (", year, ")")))) |>
  ggplot(aes(est, lab, colour = p < 0.05)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_linerange(aes(xmin = pmax(lo, 0.38), xmax = pmin(hi, 1.6))) + geom_point(size = 2) +   # pmax/pmin clip wide early CIs instead of dropping them
  scale_colour_manual(values = c(`FALSE` = navy, `TRUE` = red), guide = "none") +
  scale_x_continuous(limits = c(0.38, 1.6)) +
  labs(x = "Cumulative pooled RR", y = NULL)
cs <- \(s) cum_df[cum_df$study == s, ]


## -----------------------------------------------------------------------------
#| label: audit-brms
#| eval: false
# library(brms)
# fit_nn <- brm(yi | se(sei) ~ 1 + (1 | study), data = es,
#               prior = c(prior(normal(0, 1), class = Intercept), prior(normal(0, 0.5), class = sd)),
#               iter = 6000, warmup = 1000, chains = 4, seed = 2026, control = list(adapt_delta = 0.99))
# d <- as_draws_df(fit_nn)          # columns b_Intercept (mu) and sd_study__Intercept (tau)


## -----------------------------------------------------------------------------
#| label: audit-post
draws <- read.csv("data/ch23/revasc_draws.csv")
main  <- filter(draws, model == "main"); mu <- main$mu; tau <- main$tau
set.seed(2026)
theta_new <- rnorm(length(mu), mu, tau)                     # true log RR in a new setting
q_rr <- \(x) exp(quantile(x, c(0.5, 0.025, 0.975)))
mu_q <- q_rr(mu); new_q <- q_rr(theta_new)
round(c(P_benefit = mean(mu < 0), P_RR_below_0.9 = mean(mu < log(0.9)),
        P_RR_below_0.8 = mean(mu < log(0.8)), P_new_benefit = mean(theta_new < 0)), 2)


## -----------------------------------------------------------------------------
#| label: fig-audit-post
#| fig-cap: "Left: posterior for the pooled mean RR of death. Right: predictive distribution of the true RR in a new setting. Shading: reduction of at least 10% (dark), smaller reduction (light), harm (red)."
#| fig-width: 10
#| fig-height: 4
panel_lab <- c(mu = "A. Pooled mean effect", new = "B. True effect in a new setting")
dens_df <- function(x, from, to, which) {
  d <- density(exp(x), from = from, to = to, n = 1024)
  tibble::tibble(x = d$x, y = d$y, panel = panel_lab[[which]])
}
dd <- bind_rows(dens_df(mu, 0.55, 1.35, "mu"), dens_df(theta_new, 0.4, 1.7, "new")) |>
  mutate(zone = case_when(x < 0.9 ~ "a", x < 1 ~ "b", TRUE ~ "c"), panel = factor(panel, levels = panel_lab))
probs <- tibble::tibble(panel = factor(rep(panel_lab, each = 3), levels = panel_lab), zone = rep(c("a", "b", "c"), 2),
  p = c(mean(mu < log(0.9)), mean(mu >= log(0.9) & mu < 0), mean(mu >= 0),
        mean(theta_new < log(0.9)), mean(theta_new >= log(0.9) & theta_new < 0), mean(theta_new >= 0)),
  x = rep(c(0.6, 0.95, 1.3), 2)) |>
  mutate(lab = paste0(c("RR < 0.90", "0.90 to 1", "RR > 1"), "\n", sprintf("P = %.2f", p)))
ymax <- max(dd$y)
ggplot(dd, aes(x, y)) +
  geom_area(aes(fill = zone), alpha = 0.85) + geom_vline(xintercept = 1, colour = "grey30") +
  geom_label(data = probs, aes(x = x, y = ymax * 1.18, label = lab, colour = zone), size = 3.3,
             label.size = 0, fill = "white") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c(a = navy, b = "#9DB6D9", c = red), guide = "none") +
  scale_colour_manual(values = c(a = navy, b = "#4A7AB8", c = red), guide = "none") +
  scale_y_continuous(limits = c(0, ymax * 1.3)) +
  labs(x = "Relative risk of death, complete vs culprit-only", y = "Posterior density")


## -----------------------------------------------------------------------------
#| label: audit-arm
#| eval: false
# long <- dat |> transmute(study, e_1 = ai, n_1 = n1i, e_0 = ci, n_0 = n2i) |>
#   pivot_longer(-study, names_to = c(".value", "trt"), names_sep = "_") |>
#   mutate(trt = as.numeric(trt) - 0.5)                         # +/- 0.5 coding
# fit_arm <- brm(e ~ 0 + factor(study) + trt + (0 + trt | study) + offset(log(n)), family = poisson(),
#                data = long, prior = c(prior(normal(0, 1), class = b, coef = trt),
#                                       prior(normal(-3, 2), class = b), prior(normal(0, 0.5), class = sd)),
#                iter = 6000, warmup = 1000, chains = 4, seed = 2026, control = list(adapt_delta = 0.99))


## -----------------------------------------------------------------------------
#| label: tbl-audit-sens
#| tbl-cap: "Pooled mean and new-setting relative risks of death under four Bayesian analyses."
lab_m <- c(main = "Main: τ ~ half-N(0.5)", tau025 = "τ ~ half-N(0.25)", tau1 = "τ ~ half-N(1)",
           arm = "Arm-based Poisson, raw counts")
ci_x <- \(v) ci_t(exp(median(v)), exp(quantile(v, .025)), exp(quantile(v, .975)))
draws |> group_by(model) |>
  group_modify(\(d, ...) { set.seed(2026); xn <- rnorm(nrow(d), d$mu, d$tau)
    tibble::tibble(`Pooled mean RR (95% CrI)` = ci_x(d$mu), `P(RR < 1)` = f2(mean(d$mu < 0)),
                   `New-setting RR (95% PrI)` = ci_x(xn), `P(new RR < 1)` = f2(mean(xn < 0))) }) |>
  ungroup() |> mutate(Analysis = lab_m[model], model = match(model, names(lab_m))) |>
  arrange(model) |> select(Analysis, everything(), -model) |> kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: audit-absolute
cer <- sum(dat$ci) / sum(dat$n2i); r0 <- dat$ci / dat$n2i
arr <- cer * (1 - exp(mu))                                   # one ARR per posterior draw
arr_q <- quantile(arr, c(0.5, 0.025, 0.975))
rd <- escalc(measure = "RD", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat)
m_rd <- rma(yi, vi, data = rd, method = "REML")
c(baseline = cer, min_r0 = min(r0), max_r0 = max(r0)) |> round(3)
# pooled risk difference per 100 patients, with 95% CI and P
round(c(RD_per_100 = 100 * coef(m_rd)[[1]], lower = 100 * m_rd$ci.lb,
        upper = 100 * m_rd$ci.ub, P = m_rd$pval), 3)


## -----------------------------------------------------------------------------
#| label: audit-cost
set.seed(2026); N <- length(mu)
p0     <- rbeta(N, sum(dat$ci), sum(dat$n2i) - sum(dat$ci))
arr_h  <- p0 * (1 - exp(mu))                                  # deaths prevented per patient
cost   <- rgamma(N, shape = (8000 / 3000)^2, rate = 8000 / 3000^2)
saving <- runif(N, 0, 1500)
ly     <- runif(N, 5, 12)
cpd    <- ifelse(arr_h > 0, (cost - saving) / arr_h, Inf)     # $ per death prevented
cpl    <- cpd / ly                                            # $ per life-year gained
# medians and 95% intervals; Inf = draws in which no deaths are prevented
cq <- \(x) { v <- quantile(x, c(0.025, 0.5, 0.975)); ifelse(is.finite(v), big(v), "Inf") }
print(data.frame(centile = c("2.5th", "50th", "97.5th"),
                 deaths_per_100 = round(100 * quantile(arr_h, c(0.025, 0.5, 0.975)), 1),
                 cost_per_death = cq(cpd), cost_per_LY = cq(cpl)), row.names = FALSE)


## -----------------------------------------------------------------------------
#| label: fig-audit-cost
#| fig-cap: "Illustrative distribution of the cost per life-year gained with complete revascularization (log scale), combining the posterior for the relative risk with assumed costs and life-years. Draws with no deaths prevented are not shown."
#| fig-width: 9
#| fig-height: 4
thr <- c(50000, 100000); p_thr <- sapply(thr, \(k) mean(cpl < k))
fin <- cpl[is.finite(cpl) & cpl > 0]
y_top <- max(hist(log10(fin), breaks = 80, plot = FALSE)$counts)
ggplot(tibble::tibble(cpl = fin), aes(cpl)) +
  geom_histogram(bins = 80, fill = navy, alpha = 0.8) +
  geom_vline(xintercept = thr, colour = red, linetype = "dashed") +
  annotate("label", x = thr * c(1 / 1.08, 1.08), y = y_top * 1.12, hjust = c(1, 0), size = 3.5,
           label.size = 0, fill = "white", lineheight = 0.95,
           label = sprintf("$%sk per LY\nP(below) = %.2f", thr / 1000, p_thr)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, y_top * 1.25)) +
  scale_x_log10(labels = scales::label_dollar(), breaks = c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6)) +
  coord_cartesian(xlim = c(1e4, 5e6)) +
  labs(x = "Cost per life-year gained (CA$, log scale)", y = "Draws")

