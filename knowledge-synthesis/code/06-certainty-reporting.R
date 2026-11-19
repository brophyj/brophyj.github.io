# R code for chapter: 06-certainty-reporting.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch06
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
suppressPackageStartupMessages({
  library(metafor)
  library(robvis)
})


## -----------------------------------------------------------------------------
#| label: fig-emperor-rob
#| fig-cap: "RoB 2 assessment of the primary outcome of EMPEROR-Reduced."
#| fig-width: 5.5
#| fig-height: 2.3
#| out-width: "65%"
rob_traffic_light(data.frame(Study = "EMPEROR-Reduced", D1 = "Low", D2 = "Low",
                             D3 = "Low", D4 = "Low", D5 = "Low", Overall = "Low"),
                  tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: tbl-sglt2-two
#| tbl-cap: "Events and patients in four SGLT2 inhibitor trials."
sglt2 <- data.frame(
  study  = c("DAPA-HF", "EMPEROR-Reduced", "EMPEROR-Preserved", "DELIVER"),
  n_I    = c(2373, 1863, 2997, 3131),
  n_C    = c(2371, 1867, 2991, 3132),
  comp_I = c(386, 361, 415, 512),     # CV death or HF hospitalisation
  comp_C = c(502, 462, 511, 610),
  acm_I  = c(276, 249, 422, 497),     # death from any cause
  acm_C  = c(329, 266, 427, 526))
kable(sglt2, col.names = c("Trial", "N SGLT2i", "N placebo", "Composite SGLT2i",
                           "Composite placebo", "Deaths SGLT2i", "Deaths placebo"))


## -----------------------------------------------------------------------------
#| label: fit-both
comp <- escalc(measure = "RR", ai = comp_I, n1i = n_I, ci = comp_C, n2i = n_C,
               data = sglt2, slab = study)
acm  <- escalc(measure = "RR", ai = acm_I, n1i = n_I, ci = acm_C, n2i = n_C,
               data = sglt2, slab = study)
res_comp <- rma(yi, vi, data = comp, method = "REML")
res_acm  <- rma(yi, vi, data = acm,  method = "REML")


## -----------------------------------------------------------------------------
#| label: fig-fe-re
#| fig-cap: "The same five studies summarised by a fixed-effect model (left) and a random-effects model (right)."
#| fig-width: 12
#| out-width: "100%"
#| fig-height: 3.8
toy <- data.frame(study = paste("Study", LETTERS[1:5]),
                  yi = c(0.10, 0.85, 0.30, 1.05, 0.55),
                  vi = c(0.015, 0.020, 0.018, 0.022, 0.016))
fe_toy <- rma(yi, vi, data = toy, method = "FE",   slab = study)
re_toy <- rma(yi, vi, data = toy, method = "REML", slab = study)
par(mfrow = c(1, 2))
forest(fe_toy, header = "Study", xlab = "Effect size", mlab = "Fixed effect",
       xlim = c(-1.8, 2.7), alim = c(-0.5, 1.5))
forest(re_toy, header = "Study", xlab = "Effect size", mlab = "Random effects",
       xlim = c(-1.8, 2.7), alim = c(-0.5, 1.5))


## -----------------------------------------------------------------------------
#| label: het-sglt2
pi_comp <- predict(res_comp, transf = exp)
pi_acm  <- predict(res_acm,  transf = exp)
het <- data.frame(
  Outcome = c("CV death or HF hospitalisation", "Death from any cause"),
  RR      = c(pi_comp$pred, pi_acm$pred),
  CI      = sprintf("%.2f to %.2f", c(pi_comp$ci.lb, pi_acm$ci.lb), c(pi_comp$ci.ub, pi_acm$ci.ub)),
  PI      = sprintf("%.2f to %.2f", c(pi_comp$pi.lb, pi_acm$pi.lb), c(pi_comp$pi.ub, pi_acm$pi.ub)),
  Q_p     = c(res_comp$QEp, res_acm$QEp),
  I2      = c(res_comp$I2, res_acm$I2),
  tau2    = c(res_comp$tau2, res_acm$tau2))
kable(het, digits = c(0, 2, 0, 0, 2, 0, 4),
      col.names = c("Outcome", "Pooled RR", "95% CI", "95% PI", "P (Q)", "I² (%)", "τ²"))


## -----------------------------------------------------------------------------
#| label: fig-sglt2-pi
#| fig-cap: "Cardiovascular death or hospitalisation for heart failure. The line beneath the diamond is the 95% prediction interval."
#| fig-width: 8.5
#| fig-height: 3.6
forest(res_comp, atransf = exp, addpred = TRUE, header = "Trial",
       at = log(c(0.6, 0.8, 1, 1.25)), xlab = "Risk ratio (log scale)")


## -----------------------------------------------------------------------------
#| label: bcg-het
bcg <- rma(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
           data = dat.bcg, method = "REML")
pi_bcg <- predict(bcg, transf = exp)


## -----------------------------------------------------------------------------
#| label: ois
ois_binary <- function(p_C, RR, alpha = 0.05, power = 0.80) {
  p_I <- p_C * RR
  z   <- qnorm(1 - alpha / 2) + qnorm(power)
  ceiling(z^2 * (p_C * (1 - p_C) + p_I * (1 - p_I)) / (p_C - p_I)^2)   # per arm
}

p_C_comp <- sum(sglt2$comp_C) / sum(sglt2$n_C)    # pooled control-group risk
p_C_acm  <- sum(sglt2$acm_C)  / sum(sglt2$n_C)
n_total  <- sum(sglt2$n_I + sglt2$n_C)

ois <- data.frame(
  Outcome  = c("CV death or HF hospitalisation", "Death from any cause"),
  p_C      = c(p_C_comp, p_C_acm),
  RR       = c(pi_comp$pred, pi_acm$pred),
  OIS      = 2 * c(ois_binary(p_C_comp, pi_comp$pred), ois_binary(p_C_acm, pi_acm$pred)),
  N        = n_total)
kable(ois, digits = c(0, 3, 2, 0, 0),
      col.names = c("Outcome", "Control risk", "Pooled RR", "OIS (total)", "N in review"),
      format.args = list(big.mark = ","))


## -----------------------------------------------------------------------------
#| label: tbl-sof
#| tbl-cap: "Summary of findings: SGLT2 inhibitors compared with placebo in heart failure."
per1000 <- function(p_C, rr) round(1000 * p_C * (rr - 1))
sof <- data.frame(
  Outcome = c("CV death or HF hospitalisation", "Death from any cause"),
  Studies = sprintf("4 RCTs; %s patients", format(n_total, big.mark = ",")),
  Relative = sprintf("RR %.2f (%.2f to %.2f)",
                     c(pi_comp$pred, pi_acm$pred), c(pi_comp$ci.lb, pi_acm$ci.lb),
                     c(pi_comp$ci.ub, pi_acm$ci.ub)),
  Control = sprintf("%.0f per 1000", 1000 * c(p_C_comp, p_C_acm)),
  Absolute = sprintf("%d per 1000 (%d to %d)",
                     c(per1000(p_C_comp, pi_comp$pred), per1000(p_C_acm, pi_acm$pred)),
                     c(per1000(p_C_comp, pi_comp$ci.lb), per1000(p_C_acm, pi_acm$ci.lb)),
                     c(per1000(p_C_comp, pi_comp$ci.ub), per1000(p_C_acm, pi_acm$ci.ub))),
  Certainty = c("⊕⊕⊕⊕ High", "⊕⊕⊕◯ Moderate"),
  Reason = c("", "Rated down once for imprecision: CI includes effects too small to matter"))
kable(sof, col.names = c("Outcome", "Studies and participants", "Relative effect (95% CI)",
                         "Risk with placebo", "Difference with SGLT2i", "Certainty", "Comments"))


## -----------------------------------------------------------------------------
#| label: fig-sof-two
#| fig-cap: "Two answers for each outcome: the absolute effect per 1000 patients (point and 95% CI, left) and the certainty of the evidence (right). The idea of separating the two is captured in a cartoon by Leo Garel."
#| fig-width: 8.5
#| fig-height: 2.4
#| fig-alt: "For CV death or HF hospitalisation, a reduction of about 40 per 1000 with high certainty (four filled symbols); for death from any cause, a reduction of about 10 per 1000 whose interval reaches zero, with moderate certainty (three filled symbols)."
sof_plot <- data.frame(
  outcome = factor(c("CV death or HF hospitalisation", "Death from any cause"),
                   levels = rev(c("CV death or HF hospitalisation", "Death from any cause"))),
  est = c(per1000(p_C_comp, pi_comp$pred),  per1000(p_C_acm, pi_acm$pred)),
  lo  = c(per1000(p_C_comp, pi_comp$ci.lb), per1000(p_C_acm, pi_acm$ci.lb)),
  hi  = c(per1000(p_C_comp, pi_comp$ci.ub), per1000(p_C_acm, pi_acm$ci.ub)),
  certainty = c("\u2295\u2295\u2295\u2295  High", "\u2295\u2295\u2295\u25ef  Moderate"))
xmax <- max(sof_plot$hi) + 25
ggplot(sof_plot, aes(est, outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#041E42") +
  geom_text(aes(x = xmax, label = certainty), hjust = 0, size = 4) +
  scale_x_continuous(limits = c(min(sof_plot$lo) - 5, xmax + 40),
                     breaks = seq(-50, 10, 10)) +
  labs(x = "Difference with SGLT2 inhibitor, events per 1000 patients", y = NULL)

