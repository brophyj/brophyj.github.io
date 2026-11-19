# R code for chapter: 16-ipd-meta-analysis.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch16
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
  library(survival)
  library(lme4)
})


## -----------------------------------------------------------------------------
#| label: fig-ad-vs-ipd
#| fig-cap: "Six simulated trials in which the benefit grows with age (green dashed line: 4 units at age 65, 0.4 units per year). Left: aggregate data, one point per trial at its mean age, with a meta-regression. Right: IPD from the same patients, with the age interaction estimated within trials; points are age-band estimates."
#| fig-width: 10
#| fig-height: 4.3
#| out-width: "100%"
set.seed(42)
trials <- tibble::tibble(study = paste("Trial", 1:6), n = c(150, 220, 120, 300, 180, 160),
                         mean_age = c(63, 66, 64, 68, 65, 67), base = c(50, 55, 48, 62, 45, 52))
true_benefit <- \(age) 4 + 0.4 * (age - 65)
ipd_age <- trials |>
  rowwise() |>
  reframe(study = study, age = rnorm(n, mean_age, 10), trt = rbinom(n, 1, 0.5),
          y = base - 0.2 * (age - 65) - trt * true_benefit(age) + rnorm(n, 0, 12)) |>
  mutate(.by = study, age_c = age - mean(age))
ad <- ipd_age |>                                           # one benefit per trial
  summarise(.by = study, mean_age = mean(age),
            yi = mean(y[trt == 0]) - mean(y[trt == 1]),
            vi = var(y[trt == 1]) / sum(trt == 1) + var(y[trt == 0]) / sum(trt == 0))
mr <- rma(yi, vi, mods = ~ I(mean_age - 65), data = ad)
grid <- tibble::tibble(age = 40:90)
ad_line <- as.data.frame(predict(mr, newmods = grid$age - 65)) |> bind_cols(grid)
# IPD: trial-specific intercepts and age slopes; interaction from within-trial (centred) age
fit_age <- lm(y ~ study + study:age_c + trt + trt:age_c, data = ipd_age)
ctr <- mean(trials$mean_age)
X <- cbind(1, grid$age - ctr)
b <- -coef(fit_age)[c("trt", "age_c:trt")]                 # benefit = minus the effect on y
V <- vcov(fit_age)[c("trt", "age_c:trt"), c("trt", "age_c:trt")]
ipd_line <- grid |> mutate(pred = drop(X %*% b), se = sqrt(rowSums((X %*% V) * X)),
                           ci.lb = pred - 1.96 * se, ci.ub = pred + 1.96 * se)
bands <- ipd_age |>
  mutate(band = cut(age, c(-Inf, 50, 60, 70, 80, Inf))) |>
  group_by(band) |>
  group_modify(\(d, k) { m <- lm(y ~ study + trt, data = d)
    tibble::tibble(age = median(d$age), est = -coef(m)[["trt"]], se = sqrt(vcov(m)["trt", "trt"])) }) |>
  ungroup()
slope_ad  <- coef(summary(mr))[2, ]
slope_ipd <- c(b[2], -rev(confint(fit_age)["age_c:trt", ]))
panel <- function(line, pts, x, title) {
  ggplot() +
    geom_ribbon(data = line, aes(age, ymin = ci.lb, ymax = ci.ub), fill = "#041E42", alpha = 0.15) +
    geom_line(data = line, aes(age, pred), colour = "#041E42", linewidth = 1) +
    geom_function(fun = true_benefit, colour = "#2E7D32", linetype = "dashed", linewidth = 1) +
    geom_pointrange(data = pts, aes(.data[[x]], est, ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                    colour = "#B22222", size = 0.3) +
    coord_cartesian(xlim = c(40, 90), ylim = c(-12, 16)) +
    labs(x = "Age (years)", y = "Benefit (control − treated)", title = title)
}
patchwork::wrap_plots(
  panel(ad_line, mutate(ad, est = yi, se = sqrt(vi)), "mean_age",
        sprintf("Aggregate data: slope %.2f (%.2f to %.2f)", slope_ad["estimate"], slope_ad["ci.lb"], slope_ad["ci.ub"])),
  panel(ipd_line, bands, "age",
        sprintf("IPD: slope %.2f (%.2f to %.2f)", slope_ipd[1], slope_ipd[2], slope_ipd[3])))


## -----------------------------------------------------------------------------
#| label: one-two-stage
sim_study <- function(n, log_or, id, shift = 0) {
  age   <- rnorm(n, 68, 10)
  chd   <- rbinom(n, 1, plogis(-2 + 0.04 * age))
  nsaid <- rbinom(n, 1, 0.30)
  mi    <- rbinom(n, 1, plogis(-3 + log_or * nsaid + 0.04 * age + 0.5 * chd + shift))
  data.frame(study = id, mi, nsaid, age = age - 68, chd)
}
set.seed(42)
studies <- bind_rows(sim_study(6000, log(1.48), "RAMQ", 0.10),
                     sim_study(4500, log(1.55), "Finland", -0.10),
                     sim_study(450, log(1.40), "GPRD", 0),
                     sim_study(580, log(1.52), "Saskatchewan", 0.05))
# Two-stage: the same logistic model in each study, then a random-effects meta-analysis
stage1 <- studies |>
  group_by(study) |>
  group_modify(\(d, k) { m <- glm(mi ~ nsaid + age + chd, family = binomial, data = d)
    tibble::tibble(yi = coef(m)[["nsaid"]], vi = vcov(m)["nsaid", "nsaid"]) }) |>
  ungroup()
two <- rma(yi, vi, data = stage1, method = "REML")
# One-stage: one model; study-specific intercepts and covariate effects, random exposure effect
one <- glmer(mi ~ 0 + study + study:age + study:chd + nsaid + (0 + nsaid | study),
             family = binomial, data = studies, control = glmerControl(optimizer = "bobyqa"))
b1 <- fixef(one)[["nsaid"]]; se1 <- sqrt(vcov(one)["nsaid", "nsaid"])
tibble::tibble(Analysis = c("Two-stage (REML)", "One-stage (glmer)"),
               `Pooled OR (95% CI)` = c(sprintf("%.2f (%.2f to %.2f)", exp(two$b[1]), exp(two$ci.lb), exp(two$ci.ub)),
                                        sprintf("%.2f (%.2f to %.2f)", exp(b1), exp(b1 - 1.96 * se1), exp(b1 + 1.96 * se1))),
               τ = c(sqrt(two$tau2), attr(VarCorr(one)$study, "stddev"))) |>
  kable(digits = 2, align = "lcc")


## -----------------------------------------------------------------------------
#| label: tbl-nsaid
#| tbl-cap: "Selected results of the NSAID IPD meta-analysis: pooled adjusted odds ratios of acute myocardial infarction versus non-use, with 95% credible intervals and posterior probabilities [@bally2017nsaid]."
nsaid <- tibble::tribble(
  ~Drug, ~Category, ~or, ~lo, ~hi, ~p1, ~p15,
  "Celecoxib", "Any dose, 1–7 days", 1.24, 0.91, 1.82, 92.4, 12.1,
  "Celecoxib", "> 200 mg/day, 8–30 days", 1.23, 0.78, 1.80, 85.8, 12.5,
  "Diclofenac", "Any dose, 1–7 days", 1.50, 1.06, 2.04, 98.6, 50.1,
  "Diclofenac", "> 100 mg/day, > 30 days", 1.48, 1.08, 1.95, 98.8, 45.6,
  "Ibuprofen", "Any dose, 1–7 days", 1.48, 1.00, 2.26, 97.3, 47.1,
  "Ibuprofen", "> 1200 mg/day, 8–30 days", 1.75, 1.00, 2.93, 97.5, 70.8,
  "Naproxen", "Any dose, 1–7 days", 1.53, 1.07, 2.33, 98.8, 54.9,
  "Naproxen", "> 750 mg/day, 8–30 days", 1.76, 1.14, 2.65, 99.2, 77.2,
  "Rofecoxib", "Any dose, 1–7 days", 1.58, 1.07, 2.17, 98.8, 63.6,
  "Rofecoxib", "> 25 mg/day, 8–30 days", 2.65, 1.46, 4.67, 99.8, 96.9)
nsaid |>
  transmute(Drug, Category, `OR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", or, lo, hi),
            `P(OR > 1)` = sprintf("%.1f%%", p1), `P(OR > 1.5)` = sprintf("%.1f%%", p15)) |>
  kable(align = "llccc")


## -----------------------------------------------------------------------------
#| label: fig-nsaid-prob
#| fig-cap: "Approximate posterior probability that the odds ratio for 1 to 7 days of use at any dose exceeds each threshold (normal approximation on the log scale, from the published 95% credible intervals)."
#| fig-width: 7.5
#| fig-height: 4.2
nsaid |>
  filter(Category == "Any dose, 1–7 days") |>
  mutate(m = log(or), s = (log(hi) - log(lo)) / (2 * 1.96)) |>
  crossing(threshold = seq(0.8, 2.5, by = 0.01)) |>
  mutate(prob = 100 * pnorm(log(threshold), m, s, lower.tail = FALSE)) |>
  ggplot(aes(threshold, prob, colour = Drug)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dotted") +
  scale_x_log10(breaks = c(0.8, 1, 1.25, 1.5, 2, 2.5)) +
  labs(x = "Odds ratio threshold (log scale)", y = "P(OR > threshold), %", colour = NULL)


## -----------------------------------------------------------------------------
#| label: km-digitised
read_digitised <- function(file) {
  read.csv(file, header = FALSE, col.names = c("time", "cuminc")) |>
    bind_rows(tibble::tibble(time = 0, cuminc = 0)) |>
    arrange(time) |>
    mutate(cuminc = cummax(pmax(cuminc, 0)))                 # never negative, never decreasing
}
dig_tavr <- read_digitised("data/ch16/early_tavr_int.csv")
dig_cs   <- read_digitised("data/ch16/early_tavr_cs.csv")
at <- \(d, t) approx(d$time, d$cuminc, xout = t, rule = 2, ties = max)$y
tibble::tibble(Arm = rep(c("TAVR", "Surveillance"), each = 2), Year = c(2, 5, 2, 5),
               `Printed (%)` = c(17.7, 35.1, 37.5, 51.2),
               `Digitised (%)` = c(at(dig_tavr, c(2, 5)), at(dig_cs, c(2, 5)))) |>
  kable(digits = 1, align = "lccc")


## -----------------------------------------------------------------------------
#| label: km-reconstruct
t_risk    <- 0:5
risk_tavr <- c(455, 390, 363, 285, 142, 103)                 # numbers at risk, from the figure
risk_cs   <- c(446, 305, 266, 187, 117, 46)
reconstruct_ipd <- function(dig, at_risk, t_risk, arm) {
  out <- map(seq_len(length(t_risk) - 1), \(i) {
    a <- t_risk[i]; b <- t_risk[i + 1]
    ci <- at(dig, c(a, b)) / 100
    leaving <- at_risk[i] - at_risk[i + 1]
    p <- 1 - (1 - ci[2]) / (1 - ci[1])                       # risk in (a, b], given event-free at a
    events <- min(round((at_risk[i] - leaving / 2) * p / (1 - p / 2)), leaving)
    tibble::tibble(time = runif(leaving, a, b), event = rep(c(1, 0), c(events, leaving - events)))
  })
  bind_rows(out, tibble::tibble(time = last(t_risk), event = rep(0, last(at_risk)))) |>
    mutate(arm = arm)
}
make_ipd <- \() bind_rows(reconstruct_ipd(dig_tavr, risk_tavr, t_risk, "TAVR"),
                          reconstruct_ipd(dig_cs, risk_cs, t_risk, "Surveillance")) |>
  mutate(arm = factor(arm, levels = c("Surveillance", "TAVR")))
set.seed(2026)
ipd <- make_ipd()


## -----------------------------------------------------------------------------
#| label: fig-km-validate
#| fig-cap: "EARLY TAVR: digitised published curves (points) and Kaplan-Meier curves of the reconstructed data (lines)."
#| fig-width: 7
#| fig-height: 4.2
km <- survfit(Surv(time, event) ~ arm, data = ipd)
km_df <- tibble::tibble(time = km$time, cuminc = 100 * (1 - km$surv),
                        arm = rep(c("Surveillance", "TAVR"), km$strata))
bind_rows(mutate(dig_tavr, arm = "TAVR"), mutate(dig_cs, arm = "Surveillance")) |>
  ggplot(aes(time, cuminc, colour = arm)) +
  geom_point(size = 0.6, alpha = 0.35) +
  geom_step(data = km_df, linewidth = 0.9) +
  scale_colour_manual(values = c(Surveillance = "#E69F00", TAVR = "#041E42"), name = NULL) +
  labs(x = "Years since randomisation", y = "Cumulative incidence (%)") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: tbl-km-validate
#| tbl-cap: "EARLY TAVR: published values and the reconstructed data."
hr_txt <- function(fit) { x <- exp(c(coef(fit), confint(fit))); sprintf("%.2f (%.2f to %.2f)", x[1], x[2], x[3]) }
cox_full <- coxph(Surv(time, event) ~ arm, data = ipd)
s <- summary(km, times = c(2, 4.99))$surv
tibble::tibble(
  Quantity = c("Events, TAVR", "Events, surveillance", "Cumulative incidence at 2 years, TAVR / surveillance (%)",
               "Cumulative incidence at 5 years, TAVR / surveillance (%)", "Hazard ratio (95% CI)"),
  Published = c("122", "202", "17.7 / 37.5", "35.1 / 51.2", "0.50 (0.40 to 0.63)"),
  Reconstructed = c(sum(ipd$event[ipd$arm == "TAVR"]), sum(ipd$event[ipd$arm == "Surveillance"]),
                    sprintf("%.1f / %.1f", 100 * (1 - s[3]), 100 * (1 - s[1])),
                    sprintf("%.1f / %.1f", 100 * (1 - s[4]), 100 * (1 - s[2])), hr_txt(cox_full))) |>
  kable(align = "lcc")


## -----------------------------------------------------------------------------
#| label: km-repeat
hr_rep <- replicate(200, {
  d <- make_ipd()
  l <- d |> filter(time > 1) |> mutate(time = time - 1)
  c(whole = exp(coef(coxph(Surv(time, event) ~ arm, d))),
    landmark = exp(coef(coxph(Surv(time, event) ~ arm, l))))
})
round(apply(hr_rep, 1, quantile, c(0.025, 0.5, 0.975)), 2)


## -----------------------------------------------------------------------------
#| label: tbl-km-landmark
#| tbl-cap: "EARLY TAVR, reconstructed data: hazard ratios by period, and restricted mean event-free time (RMST) differences, TAVR minus surveillance, in months."
yr1 <- ipd |> mutate(event = event * (time <= 1), time = pmin(time, 1))
lmk <- ipd |> filter(time > 1) |> mutate(time = time - 1)        # event-free at 1 year; clock restarted
rmst_diff <- function(data, horizon) {
  tab <- summary(survfit(Surv(time, event) ~ arm, data = data), rmean = horizon)$table
  d <- unname(diff(tab[, "rmean"])); se <- sqrt(sum(tab[, "se(rmean)"]^2))
  sprintf("%.1f (%.1f to %.1f)", 12 * d, 12 * (d - 1.96 * se), 12 * (d + 1.96 * se))
}
tibble::tibble(Period = c("Whole follow-up (5 years)", "Year 1 only", "After year 1 (landmark, 4 years)"),
               `HR (95% CI)` = c(hr_txt(cox_full), hr_txt(coxph(Surv(time, event) ~ arm, yr1)),
                                 hr_txt(coxph(Surv(time, event) ~ arm, lmk))),
               `RMST difference, months (95% CI)` = c(rmst_diff(ipd, 5), rmst_diff(ipd, 1), rmst_diff(lmk, 4))) |>
  kable(align = "lcc")

