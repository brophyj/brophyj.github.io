# R code for chapter: 08-effect-measures.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch08
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
  library(epiR)
})


## -----------------------------------------------------------------------------
#| label: z-two-studies
two_studies <- data.frame(study = c("A", "B"), effect = c(0.50, 0.50), se = c(0.20, 0.05)) |>
  mutate(z = effect / se, p = 2 * pnorm(-abs(z)))
kable(two_studies, digits = c(0, 2, 2, 1, 4))


## -----------------------------------------------------------------------------
#| label: binary-hand
a <- 15; b <- 85; c <- 30; d <- 70
n1 <- a + b; n0 <- c + d
p1 <- a / n1; p0 <- c / n0

rd    <- p1 - p0
se_rd <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
log_rr    <- log(p1 / p0)
se_log_rr <- sqrt(1/a - 1/n1 + 1/c - 1/n0)
log_or    <- log((a * d) / (b * c))
se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)

binary_ex <- tibble::tibble(
  measure  = c("Risk difference", "Risk ratio", "Odds ratio"),
  estimate = c(rd, exp(log_rr), exp(log_or)),
  lower    = c(rd - 1.96 * se_rd, exp(log_rr - 1.96 * se_log_rr), exp(log_or - 1.96 * se_log_or)),
  upper    = c(rd + 1.96 * se_rd, exp(log_rr + 1.96 * se_log_rr), exp(log_or + 1.96 * se_log_or)))
kable(binary_ex, digits = 3)


## -----------------------------------------------------------------------------
#| label: binary-packages
tab <- as.table(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Treated", "Control"), c("Event", "No event"))))
epi_res <- epi.2by2(dat = tab, method = "cohort.count")
epi_res$massoc.summary[1:3, ]

sapply(c(RD = "RD", RR = "RR", OR = "OR"), function(m)
  unlist(escalc(measure = m, ai = a, bi = b, ci = c, di = d)[, c("yi", "vi")]))


## -----------------------------------------------------------------------------
#| label: noncollapsible
risk_from_or <- function(p0, or) { odds <- p0 / (1 - p0) * or; odds / (1 + odds) }
p0_strata <- c(low = 0.2, high = 0.6)
p1_strata <- risk_from_or(p0_strata, or = 2)
p0_all <- mean(p0_strata); p1_all <- mean(p1_strata)

collapse <- rbind(
  "Within each stratum" = c(OR = 2, RR_low = p1_strata[[1]] / p0_strata[[1]],
                            RR_high = p1_strata[[2]] / p0_strata[[2]]),
  "Whole trial"         = c(OR = (p1_all / (1 - p1_all)) / (p0_all / (1 - p0_all)),
                            RR_low = p1_all / p0_all, RR_high = p1_all / p0_all))
round(collapse, 2)


## -----------------------------------------------------------------------------
#| label: nnt-baseline
nnt_tab <- tibble::tibble(baseline_risk = c(0.30, 0.10, 0.02)) |>
  mutate(rd = baseline_risk * (0.5 - 1), nnt = ceiling(1 / abs(rd)))
kable(nnt_tab, digits = 3, col.names = c("Baseline risk", "Risk difference", "NNT"))


## -----------------------------------------------------------------------------
#| label: fig-baseline-risk
#| fig-cap: "Risk difference and odds ratio when the risk ratio is constant at 0.50, as baseline (control) risk varies."
#| fig-width: 8
#| fig-height: 4
tibble::tibble(p0 = seq(0.01, 0.6, by = 0.01)) |>
  mutate(`Risk ratio` = 0.5,
         `Odds ratio` = (0.5 * p0 / (1 - 0.5 * p0)) / (p0 / (1 - p0)),
         `Risk difference` = p0 * (0.5 - 1)) |>
  pivot_longer(-p0, names_to = "measure", values_to = "value") |>
  ggplot(aes(p0, value, colour = measure)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Risk ratio" = "#041E42", "Odds ratio" = "#F28E2B",
                                 "Risk difference" = "#C8102E"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Baseline (control) risk", y = "Value of the measure") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: md-variance
md_ex <- data.frame(m1i = 103, sd1i = 5.5, n1i = 40, m2i = 100, sd2i = 4.5, n2i = 60)
sp2   <- with(md_ex, ((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2) / (n1i + n2i - 2))
c(escalc_default  = escalc(measure = "MD", m1i = m1i, sd1i = sd1i, n1i = n1i,
                           m2i = m2i, sd2i = sd2i, n2i = n2i, data = md_ex)$vi,
  unequal_by_hand = with(md_ex, sd1i^2 / n1i + sd2i^2 / n2i),
  pooled_sd       = with(md_ex, sp2 * (1 / n1i + 1 / n2i))) |> round(3)


## -----------------------------------------------------------------------------
#| label: paired-r
paired <- tibble::tibble(r = c(0, 0.5, 0.9, 0.99)) |>
  mutate(sd_diff = sqrt(10^2 + 10^2 - 2 * r * 10 * 10), se_md = sd_diff / sqrt(50))
kable(paired, digits = 2, col.names = c("Correlation r", "SD of differences", "SE of MD"),
      caption = "Standard error of a paired mean difference with SDs of 10 and n = 50.")


## -----------------------------------------------------------------------------
#| label: smd-example
escalc(measure = "SMD", m1i = m1i, sd1i = sd1i, n1i = n1i,
       m2i = m2i, sd2i = sd2i, n2i = n2i, data = md_ex)


## -----------------------------------------------------------------------------
#| label: molloy
molloy <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat.molloy2014,
                 slab = paste(authors, year, sep = ", "))
res_molloy <- rma(yi, vi, data = molloy)
predict(res_molloy, transf = transf.ztor)


## -----------------------------------------------------------------------------
#| label: fig-molloy-forest
#| fig-cap: "Correlation between conscientiousness and medication adherence in 16 studies, pooled on Fisher's z scale and back-transformed [@molloy2014]."
#| fig-width: 8
#| fig-height: 6.5
forest(res_molloy, atransf = transf.ztor,
       at = transf.rtoz(c(-0.4, -0.2, 0, 0.2, 0.4, 0.6)),
       xlab = "Correlation (r)", header = "Study")


## -----------------------------------------------------------------------------
#| label: veteran-cox
vet <- survival::veteran |>
  mutate(trt = factor(trt, levels = 1:2, labels = c("Standard", "Test")))
fit_unadj <- coxph(Surv(time, status) ~ trt, data = vet)
fit_adj   <- coxph(Surv(time, status) ~ trt + karno + age, data = vet)
rbind(unadjusted = summary(fit_unadj)$conf.int["trtTest", c(1, 3, 4)],
      adjusted   = summary(fit_adj)$conf.int["trtTest", c(1, 3, 4)]) |> round(2)


## -----------------------------------------------------------------------------
#| label: fig-veteran-km
#| fig-cap: "Kaplan-Meier survival curves for standard and test chemotherapy in the veteran lung cancer trial. The curves cross."
#| fig-width: 7
#| fig-height: 4.5
km <- survfit(Surv(time, status) ~ trt, data = vet)
tibble::tibble(time = km$time, surv = km$surv,
               arm = rep(sub("trt=", "", names(km$strata)), km$strata)) |>
  ggplot(aes(time, surv, colour = arm)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(values = c(Standard = "#C8102E", Test = "#041E42"), name = NULL) +
  labs(x = "Days", y = "Survival probability") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: veteran-zph
zph <- cox.zph(fit_adj)
zph


## -----------------------------------------------------------------------------
#| label: fig-veteran-zph
#| fig-cap: "Scaled Schoenfeld residuals against time for each covariate in the adjusted Cox model, with a smoothed trend and ±2 SE band. A trend indicates a non-proportional effect."
#| fig-width: 10
#| fig-height: 3.5
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(zph)


## -----------------------------------------------------------------------------
#| label: veteran-rmst
rmst <- summary(km, rmean = 400)$table[, c("rmean", "se(rmean)")]
rmst


## -----------------------------------------------------------------------------
#| label: fig-whi
#| fig-cap: "Hazard ratios for coronary heart disease by year of follow-up in the Women's Health Initiative trial of estrogen plus progestin. Data from Table 2 of @manson2003whi; the dashed line is the overall HR of 1.24."
#| fig-width: 7
#| fig-height: 3.8
whi <- tibble::tibble(
  year = factor(c("1", "2", "3", "4", "5", "6+"), levels = c("1", "2", "3", "4", "5", "6+")),
  cases_ht = c(42, 38, 19, 32, 29, 28), rate_ht = c(0.50, 0.45, 0.23, 0.39, 0.41, 0.37),
  cases_pl = c(23, 28, 15, 25, 19, 37), rate_pl = c(0.29, 0.35, 0.19, 0.32, 0.28, 0.56),
  hr = c(1.81, 1.34, 1.27, 1.25, 1.45, 0.70),
  lo = c(1.09, 0.82, 0.64, 0.74, 0.81, 0.42),
  hi = c(3.01, 2.18, 2.50, 2.12, 2.59, 1.14))
ggplot(whi, aes(year, hr)) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 1.24, linetype = "dashed", colour = "#C8102E") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#041E42") +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(x = "Year of follow-up", y = "Hazard ratio (log scale)")

