# R code for chapter: 24-agentic-meta-analysis.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch24
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
navy <- "#041E42"; red <- "#C8102E"; green <- "#009E73"; blue <- "#56B4E9"; orange <- "#E69F00"


## -----------------------------------------------------------------------------
#| label: agentic-data
dcis <- read.csv("data/ch24/stuart2015_invasive_lr.csv")   # n_cases and 10-year LR (%) per arm
dcis <- dcis |>
  mutate(xi = round(n_cases * lr_rate_10yr_pct / 100),     # events, reconstructed from n and %
         treatment = factor(treatment, levels = unique(treatment)))
dcis |> group_by(treatment) |>
  summarise(arms = n(), patients = sum(n_cases), events = sum(xi), zero_event_arms = sum(xi == 0)) |>
  kable()


## -----------------------------------------------------------------------------
#| label: fig-agentic-extract
#| fig-cap: "Data extraction, 10 independent runs per model and configuration. ChatGPT results from Dobler et al. (2025); reas. = reasoning mode with web retrieval. Claude Opus 5.5 (2026): our rerun with Dobler's prompt verbatim (exact) and with a vaguer prompt that did not name the radiotherapy arm (vague)."
#| fig-width: 10
#| fig-height: 4.4
lvl <- c("Correct value", "22% (whole cohort) as answer,\n11.9% flagged for radiation arm",
         "22% (whole cohort)", "Other wrong value", "No rate / wrong study")
extr <- tribble(
  ~model,                 ~study,            ~outcome, ~n,
  "GPT-4o\nreas.",         "Vidali (9.6%)",   1, 10,
  "GPT-4o\nno reas.",      "Vidali (9.6%)",   1, 10,
  "o4-mini\nreas.",        "Vidali (9.6%)",   1,  9,
  "o4-mini\nreas.",        "Vidali (9.6%)",   4,  1,
  "o4-mini\nno reas.",     "Vidali (9.6%)",   1,  9,
  "o4-mini\nno reas.",     "Vidali (9.6%)",   4,  1,
  "Claude\nexact",         "Vidali (9.6%)",   1, 10,
  "Claude\nvague",         "Vidali (9.6%)",   1, 10,
  "GPT-4o\nreas.",         "Rudloff (11.9%)", 1,  8,
  "GPT-4o\nreas.",         "Rudloff (11.9%)", 3,  2,
  "GPT-4o\nno reas.",      "Rudloff (11.9%)", 1,  9,
  "GPT-4o\nno reas.",      "Rudloff (11.9%)", 3,  1,
  "o4-mini\nreas.",        "Rudloff (11.9%)", 1,  8,
  "o4-mini\nreas.",        "Rudloff (11.9%)", 3,  2,
  "o4-mini\nno reas.",     "Rudloff (11.9%)", 1,  2,
  "o4-mini\nno reas.",     "Rudloff (11.9%)", 3,  4,
  "o4-mini\nno reas.",     "Rudloff (11.9%)", 5,  4,
  "Claude\nexact",         "Rudloff (11.9%)", 1, 10,
  "Claude\nvague",         "Rudloff (11.9%)", 2, 10) |>
  mutate(model = factor(model, levels = unique(model)),
         study = factor(study, levels = c("Vidali (9.6%)", "Rudloff (11.9%)")),
         outcome = factor(lvl[outcome], levels = rev(lvl)))
ggplot(extr, aes(model, n, fill = outcome)) +
  geom_col(width = 0.7, colour = "white") +
  facet_wrap(~ study) +
  scale_fill_manual(values = setNames(c(green, "#7FB3D5", orange, "#E9967A", red), lvl), breaks = lvl) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  labs(x = NULL, y = "Runs (of 10)", fill = NULL) +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 8.5),
        strip.text = element_text(face = "bold"))


## -----------------------------------------------------------------------------
#| label: agentic-glmm
pool_glmm <- function(d) {
  r <- rma.glmm(measure = "PLO", xi = xi, ni = n_cases, data = d)   # random-effects logistic model
  p <- predict(r, transf = transf.ilogit)
  tibble(k = nrow(d), est = p$pred, lo = p$ci.lb, hi = p$ci.ub)
}
ours <- dcis |> group_by(treatment) |> group_modify(\(d, ...) pool_glmm(d)) |> ungroup() |>
  bind_rows(pool_glmm(dcis) |> mutate(treatment = "Overall"))
ours |> mutate(across(c(est, lo, hi), \(x) round(x, 3))) |> kable()


## -----------------------------------------------------------------------------
#| label: tbl-agentic-chatgpt
#| tbl-cap: "Pooled 10-year invasive LR rates in ten ChatGPT o4-mini chats (Dobler et al. 2025, Table 2). CC = continuity correction; Mixed = both models reported, choice by an I² threshold."
chatgpt <- tibble(
  Session = c(paste("Chat", 1:10), "Expert"),
  Mastect. = c(0.03, 0.02, 0.00, 0.02, 0.02, 0.00, 0.01, 0.00, 0.03, 0.02, 0.03),
  `S+R`    = c(0.07, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.07),
  `S-R`    = c(0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.21, 0.11, 0.11),
  Biopsy   = c(0.26, 0.26, 0.25, 0.26, 0.25, 0.25, 0.25, 0.25, 0.16, 0.26, 0.25),
  Overall  = c(0.09, 0.08, 0.00, 0.07, 0.15, 0.07, 0.07, 0.00, 0.11, 0.07, 0.08),
  Transformation = c("Logit", "CC all", "none", "CC zeros", "none", "CC zeros", "CC zeros", "none",
                     "none", "CC zeros", "Logit"),
  `FE/RE` = c("Mixed", "Mixed", "RE", "Mixed", "Mixed", "Mixed", "Mixed", "RE", "Mixed", "Mixed", "RE"))
kable(chatgpt, align = "lcccccll")


## -----------------------------------------------------------------------------
#| label: tbl-agentic-claude
#| tbl-cap: "Pooled 10-year invasive LR rates in ten Claude Opus 5.5 agent sessions (our rerun, November 2026). All used a logit GLMM with random effects chosen a priori."
claude <- tibble(
  Session = paste("Run", 1:10),
  Mastect. = 0.018, `S+R` = 0.067, `S-R` = 0.107, Biopsy = 0.232,
  Overall = c(rep(0.076, 6), 0.079, rep(0.076, 3)),
  CI = c("t(k-1)", "t(k-1)", "t(k-1)", "Wald z", "t(k-1)", "t(k-1)", "t(k-1)", "Wald z", "Wald z", "Wald z"),
  `Data decision` = c(rep("kept all 30 arms", 6), "dropped \"Lara\"", rep("kept all 30 arms", 3)))
kable(claude, align = "lccccccl")


## -----------------------------------------------------------------------------
#| label: fig-agentic-variability
#| fig-cap: "Pooled rates from each session. Blue: ChatGPT o4-mini (2025); orange: Claude Opus 5.5 (2026); black diamond: the expert benchmark of Dobler et al."
#| fig-width: 9
#| fig-height: 4
trt_lv <- c("Mastect.", "S+R", "S-R", "Biopsy", "Overall")
runs <- bind_rows(
  chatgpt |> filter(Session != "Expert") |> select(all_of(trt_lv)) |> mutate(source = "ChatGPT o4-mini (2025)"),
  claude |> select(all_of(trt_lv)) |> mutate(source = "Claude Opus 5.5 (2026)")) |>
  pivot_longer(all_of(trt_lv), names_to = "treatment", values_to = "estimate") |>
  mutate(treatment = factor(treatment, levels = trt_lv))
bench <- chatgpt |> filter(Session == "Expert") |>
  pivot_longer(all_of(trt_lv), names_to = "treatment", values_to = "estimate") |>
  mutate(treatment = factor(treatment, levels = trt_lv))
ggplot(runs, aes(treatment, estimate, colour = source)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1),
             alpha = 0.75, size = 2.4) +
  geom_point(data = bench, aes(treatment, estimate), inherit.aes = FALSE, shape = 18, size = 4.5) +
  scale_colour_manual(values = c(blue, orange), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Pooled 10-year rate") +
  theme(legend.position = "top", panel.grid.major.x = element_blank())


## -----------------------------------------------------------------------------
#| label: tbl-agentic-ci
#| tbl-cap: "The same logit GLMM with normal (z) and t(k-1) confidence intervals; pooled 10-year invasive LR, %."
ci_both <- function(d) {
  f <- function(test) predict(rma.glmm(measure = "PLO", xi = xi, ni = n_cases, data = d,
                                       test = test), transf = transf.ilogit)
  z <- f("z"); t <- f("t")
  tibble(k = nrow(d), estimate = 100 * z$pred,
         `95% CI, z` = sprintf("%.1f to %.1f", 100 * z$ci.lb, 100 * z$ci.ub),
         `95% CI, t` = sprintf("%.1f to %.1f", 100 * t$ci.lb, 100 * t$ci.ub))
}
dcis |> group_by(treatment) |> group_modify(\(d, ...) ci_both(d)) |> ungroup() |>
  mutate(estimate = round(estimate, 1)) |> kable()


## -----------------------------------------------------------------------------
#| label: fig-agentic-zeros
#| fig-cap: "The mastectomy group (7 arms, 3 with zero events) pooled three ways. Raw proportions with 0.5 added to zero cells give impossible negative limits; the logit with the same correction turns each zero into a 2 to 4% 'study' and inflates the pooled rate; the binomial GLMM needs no correction."
#| fig-width: 11
#| fig-height: 3.9
mx <- filter(dcis, treatment == "Mastectomy") |> transmute(study, xi, ni = n_cases)
lv <- c(mx$study, "Pooled")
tidy_forest <- function(est, lo, hi, pooled) bind_rows(
  tibble(study = mx$study, est = est, lo = lo, hi = hi),
  tibble(study = "Pooled", est = pooled[1], lo = pooled[2], hi = pooled[3])) |>
  mutate(study = factor(study, levels = rev(lv)), pooled = study == "Pooled")
d_raw <- escalc(measure = "PR",  xi = xi, ni = ni, data = mx, add = 0.5, to = "only0")
d_lo  <- escalc(measure = "PLO", xi = xi, ni = ni, data = mx, add = 0.5, to = "only0")
r_raw <- rma(yi, vi, data = d_raw); r_lo <- rma(yi, vi, data = d_lo)
r_gl  <- rma.glmm(measure = "PLO", xi = xi, ni = ni, data = mx)
se <- \(d) sqrt(d$vi)
f_raw <- tidy_forest(d_raw$yi, d_raw$yi - 1.96 * se(d_raw), d_raw$yi + 1.96 * se(d_raw),
                     c(coef(r_raw), r_raw$ci.lb, r_raw$ci.ub))
f_lo  <- tidy_forest(plogis(d_lo$yi), plogis(d_lo$yi - 1.96 * se(d_lo)), plogis(d_lo$yi + 1.96 * se(d_lo)),
                     plogis(c(coef(r_lo), r_lo$ci.lb, r_lo$ci.ub)))
cp <- map2(mx$xi, mx$ni, \(x, n) binom.test(x, n)$conf.int)          # exact study CIs
f_gl  <- tidy_forest(mx$xi / mx$ni, map_dbl(cp, 1), map_dbl(cp, 2),
                     plogis(c(coef(r_gl), r_gl$ci.lb, r_gl$ci.ub)))
panel <- function(df, ttl, shade = FALSE) {
  ggplot(df, aes(est, study, shape = pooled, colour = pooled)) +
    { if (shade) annotate("rect", xmin = -0.08, xmax = 0, ymin = -Inf, ymax = Inf, fill = red, alpha = 0.10) } +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 0.5, size = 0.35) +
    scale_shape_manual(values = c(15, 18), guide = "none") +
    scale_colour_manual(values = c("grey20", red), guide = "none") +
    scale_x_continuous(limits = c(-0.08, 0.46), labels = scales::percent_format(accuracy = 1)) +
    labs(title = ttl, subtitle = sprintf("pooled %.1f%%", 100 * df$est[df$pooled]), x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 11))
}
panel(f_raw, "Raw proportions + 0.5", shade = TRUE) + panel(f_lo, "Logit + 0.5") + panel(f_gl, "Logit GLMM (exact)")


## -----------------------------------------------------------------------------
#| label: agentic-samplesize
hr <- 12 / 18                                    # exponential survival: ratio of medians
p  <- 0.5                                        # 1:1 allocation
d  <- (qnorm(0.975) + qnorm(0.80))^2 / (p * (1 - p) * log(hr)^2)   # events needed
lam <- log(2) / c(standard = 12, new = 18)       # hazard per month
p_event <- mean(1 - exp(-lam * 24))              # chance of death within 24 months, averaged
n_eval  <- ceiling(d) / p_event                  # evaluable patients
n_total <- 2 * ceiling(n_eval / 2 / 0.85)        # inflate for 15% dropout, equal arms
c(events = ceiling(d), evaluable = ceiling(n_eval), recruit = n_total)

