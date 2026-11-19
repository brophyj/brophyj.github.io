# R code for chapter: 04-visualising-evidence.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch04
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
  library(tidyr)
  library(metafor)
})
set.seed(20260318)


## -----------------------------------------------------------------------------
#| label: sim-null
n   <- 400                       # number of comparisons made
se  <- 1                         # standard error of each log(HR)
yi_all <- rnorm(n, mean = 0, sd = se)   # true log(HR) = 0 for every comparison

# Show 8: the most extreme positive, the most extreme negative, and 6 others at random
idx_pos  <- which.max(yi_all)
idx_neg  <- which.min(yi_all)
idx_6    <- sample(setdiff(seq_len(n), c(idx_pos, idx_neg)), 6)
idx_pick <- c(idx_6[1:3], idx_pos, idx_neg, idx_6[4:6])

outcomes <- c("ADHD", "ASD composite", "Asperger syndrome", "Atypical autism",
              "Autistic disorder", "Developmental delay", "Language disorder",
              "Other PDD")
df <- tibble(
  outcome   = outcomes,
  yi        = yi_all[idx_pick],
  ci_lo     = yi - 1.96 * se,
  ci_hi     = yi + 1.96 * se,
  outcome_f = factor(outcome, levels = rev(outcomes)))


## -----------------------------------------------------------------------------
#| label: fig-test
#| fig-cap: "Association between exposure to chemical X and eight neurodevelopmental outcomes."
#| fig-width: 8
#| fig-height: 4.5
df |>
  mutate(label = sprintf("%.2f (%.2f, %.2f)", yi, ci_lo, ci_hi)) |>
  ggplot(aes(x = yi, y = outcome_f)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7) +
  geom_point(size = 3, colour = "#041E42") +
  geom_text(aes(x = max(ci_hi) + 0.4, label = label), hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(min(df$ci_lo) - 0.5, max(df$ci_hi) + 3.5)) +
  labs(x = "log(hazard ratio)", y = NULL)


## -----------------------------------------------------------------------------
#| label: bcg-data
bcg <- escalc(measure = "RR",
              ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg,
              slab = paste(author, year))
# Fitting from the raw counts (rather than yi and vi) keeps the arm-level
# data in the model object, which labbe() needs later
res_bcg <- rma(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
               data = dat.bcg, slab = paste(author, year), method = "REML")


## -----------------------------------------------------------------------------
#| label: fig-bcg-ggplot
#| fig-cap: "Forest plot of the BCG trials drawn with ggplot2. Red: confidence interval excludes no effect. Dotted line: 95% prediction interval."
#| fig-width: 8
#| fig-height: 5.5
pred <- predict(res_bcg)
k    <- nrow(bcg)
plot_df <- data.frame(study = attr(bcg$yi, "slab"),
                      yi = bcg$yi, se = sqrt(bcg$vi)) |>
  mutate(lo = yi - 1.96 * se, hi = yi + 1.96 * se,
         sig = hi < 0 | lo > 0,
         row = rev(seq_len(k)) + 1)          # leave rows 0-1 for the summary

ggplot(plot_df, aes(x = yi, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = sig)) +
  geom_point(aes(colour = sig, size = 1 / se^2), shape = 15) +
  annotate("polygon",
           x = c(pred$ci.lb, pred$pred, pred$ci.ub, pred$pred),
           y = c(0.5, 0.85, 0.5, 0.15), fill = "#041E42") +
  annotate("segment", x = pred$pi.lb, xend = pred$pi.ub, y = 0.5, yend = 0.5,
           linetype = "dotted") +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#C8102E"),
                      guide = "none") +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  scale_y_continuous(breaks = c(plot_df$row, 0.5),
                     labels = c(plot_df$study, "Random-effects pooled")) +
  scale_x_continuous(breaks = log(c(0.1, 0.25, 0.5, 1, 2)),
                     labels = c(0.1, 0.25, 0.5, 1, 2)) +
  labs(x = "Risk ratio (log scale)", y = NULL) +
  theme(panel.grid.minor = element_blank())


## -----------------------------------------------------------------------------
#| label: fig-bcg-metafor
#| fig-cap: "The same meta-analysis drawn by metafor."
#| fig-width: 8.5
#| fig-height: 5.5
forest(res_bcg, atransf = exp, at = log(c(0.05, 0.25, 1, 4)),
       showweights = TRUE, addpred = TRUE,
       xlab = "Risk ratio (log scale)", header = "Trial")


## -----------------------------------------------------------------------------
#| label: fig-bcg-meta
#| fig-cap: "The same meta-analysis drawn by meta."
#| fig-width: 9
#| fig-height: 5.5
library(meta)
m_bcg <- metabin(event.e = tpos, n.e = tpos + tneg,
                 event.c = cpos, n.c = cpos + cneg,
                 studlab = paste(author, year), data = dat.bcg,
                 sm = "RR", method.tau = "REML",
                 common = FALSE, random = TRUE, prediction = TRUE)
forest(m_bcg)


## -----------------------------------------------------------------------------
#| label: fig-reveal
#| fig-cap: "All 400 simulated log hazard ratios. Every true effect is zero; red points are the eight shown in the forest plot."
#| fig-width: 8
#| fig-height: 2.8
tibble(yi = yi_all, picked = seq_len(n) %in% idx_pick) |>
  mutate(jitter = runif(n, -0.1, 0.1)) |>
  ggplot(aes(x = yi, y = jitter, colour = picked)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C8102E"),
                      labels = c("Not shown", "Shown in the forest plot"),
                      name = NULL) +
  labs(x = "log(hazard ratio)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: blom
blom_max <- qnorm((n - 0.375) / (n + 0.25))
blom     <- qnorm((seq_len(n) - 0.375) / (n + 0.25))   # every order statistic
blom_max


## -----------------------------------------------------------------------------
#| label: monte-carlo
maxima <- replicate(5000, max(rnorm(n)))
mc <- c(mean = mean(maxima), quantile(maxima, c(0.025, 0.975)))
round(mc, 2)


## -----------------------------------------------------------------------------
#| label: fig-maxima
#| fig-cap: "Distribution of the largest of 400 standard normal values, in 5000 simulations. Red dashed line: Blom's approximation. Dotted lines: empirical 95% range."
#| fig-width: 7
#| fig-height: 3.8
ggplot(data.frame(maxima), aes(maxima)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "grey80", colour = "grey50") +
  geom_vline(xintercept = mc[2:3], linetype = "dotted") +
  geom_vline(xintercept = blom_max, linetype = "dashed",
             colour = "#C8102E", linewidth = 1) +
  labs(x = "Largest z value in each simulation", y = "Density")


## -----------------------------------------------------------------------------
#| label: fig-qq
#| fig-cap: "Normal Q–Q plot: the sorted simulated values against the Blom-expected order statistics."
#| fig-width: 5
#| fig-height: 4.5
#| out-width: "60%"
ggplot(tibble(expected = blom, observed = sort(yi_all)),
       aes(expected, observed)) +
  geom_point(size = 0.9, alpha = 0.6, colour = "#041E42") +
  geom_abline(linetype = "dashed", colour = "#C8102E") +
  labs(x = "Blom-expected normal order statistic", y = "Observed log(HR)")


## -----------------------------------------------------------------------------
#| label: tbl-fwer
#| tbl-cap: "Probability of at least one false-positive result among m independent tests of true null hypotheses."
m <- c(1, 5, 10, 20, 50, 400)
kable(data.frame(m = m, fwer = 1 - 0.95^m), digits = 2,
      col.names = c("Tests (m)", "P(at least one false positive)"))


## -----------------------------------------------------------------------------
#| label: corrections
pvals <- tibble(z = yi_all / se,
                p = 2 * pnorm(abs(z), lower.tail = FALSE)) |>
  mutate(Bonferroni = p.adjust(p, "bonferroni"),
         Holm       = p.adjust(p, "holm"),
         BH         = p.adjust(p, "BH"))

n_sig <- pvals |>
  summarise(Uncorrected = sum(p < 0.05),
            Bonferroni  = sum(Bonferroni < 0.05),
            Holm        = sum(Holm < 0.05),
            `Benjamini-Hochberg` = sum(BH < 0.05))
kable(pivot_longer(n_sig, everything(), names_to = "Method",
                   values_to = "Significant (of 400)"))


## -----------------------------------------------------------------------------
#| label: fig-corrections
#| fig-cap: "Adjusted P values against rank for the 400 null comparisons. With so many tests, Bonferroni and Holm adjusted values reach 1 for almost every rank and their lines coincide. Dashed line: 0.05."
#| fig-width: 7
#| fig-height: 4.5
pvals |>
  arrange(p) |>
  mutate(rank = row_number()) |>
  select(rank, Uncorrected = p, `Benjamini-Hochberg` = BH, Holm, Bonferroni) |>
  pivot_longer(-rank, names_to = "method", values_to = "adj_p") |>
  mutate(method = factor(method, c("Uncorrected", "Benjamini-Hochberg",
                                   "Holm", "Bonferroni"))) |>
  ggplot(aes(rank, adj_p, colour = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  scale_y_log10() +
  scale_colour_brewer(palette = "Set1", name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "22"), name = NULL) +
  labs(x = "Rank (smallest to largest P value)", y = "Adjusted P value (log scale)") +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: funnel-sim
set.seed(42)
n_studies <- 30
se_s  <- runif(n_studies, 0.05, 0.5)
eff_s <- rnorm(n_studies, mean = -0.25, sd = se_s)
published <- !(eff_s > 0 & se_s > 0.25)        # small studies showing harm vanish
df_funnel <- tibble(yi = eff_s[published], sei = se_s[published])

res_f <- rma(yi, sei = sei, data = df_funnel, method = "REML")
egger <- regtest(res_f, model = "lm")
begg  <- ranktest(res_f)


## -----------------------------------------------------------------------------
#| label: fig-funnel
#| fig-cap: "Funnel plot of the published studies. Small studies with unfavourable results are missing from the lower right."
#| fig-width: 6
#| fig-height: 4.5
#| out-width: "70%"
funnel(res_f, xlab = "log(risk ratio)", back = "grey95")


## -----------------------------------------------------------------------------
#| label: fig-funnel-contour
#| fig-cap: "Contour-enhanced funnel plot. Studies in the white region are not statistically significant (P > 0.10)."
#| fig-width: 6
#| fig-height: 4.5
#| out-width: "70%"
funnel(res_f, level = c(90, 95, 99), shade = c("white", "grey75", "grey55"),
       refline = 0, legend = TRUE, xlab = "log(risk ratio)")


## -----------------------------------------------------------------------------
#| label: trimfill
tf <- trimfill(res_f)


## -----------------------------------------------------------------------------
#| label: fig-trimfill
#| fig-cap: "Trim and fill. Open circles are imputed studies, not data."
#| fig-width: 6
#| fig-height: 4.5
#| out-width: "70%"
funnel(tf, xlab = "log(risk ratio)", back = "grey95")


## -----------------------------------------------------------------------------
#| label: fig-labbe
#| fig-cap: "L'Abbé plot of the BCG trials. Point area is proportional to trial size. Dashed line: no effect. Solid line: the pooled random-effects risk ratio."
#| fig-width: 6
#| fig-height: 6
#| out-width: "65%"
labbe(res_bcg, transf = exp, xlab = "Risk in control arm",
      ylab = "Risk in vaccine arm")


## -----------------------------------------------------------------------------
#| label: fig-radial
#| fig-cap: "Galbraith (radial) plot of the BCG trials."
#| fig-width: 6.5
#| fig-height: 5
#| out-width: "75%"
res_bcg_ce <- rma(yi, vi, data = bcg, method = "EE")   # common-effect model
radial(res_bcg_ce)


## -----------------------------------------------------------------------------
#| label: radial-outliers
z_resid <- rstandard(res_bcg_ce)$z
outside <- attr(bcg$yi, "slab")[abs(z_resid) > 1.96]


## -----------------------------------------------------------------------------
#| label: fig-rob-traffic
#| fig-cap: "Traffic-light plot of RoB 2 judgements for five example trials."
#| fig-width: 6
#| fig-height: 4
library(robvis)
rob_traffic_light(data_rob2[1:5, ], tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: fig-rob-summary
#| fig-cap: "Summary of RoB 2 judgements across the example trials."
#| fig-width: 8
#| fig-height: 3
rob_summary(data_rob2, tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: tbl-scurvy-raw
#| tbl-cap: "Lind's scurvy trial: the raw data."
scurvy <- medicaldata::scurvy
scurvy |>
  select(-dosing_regimen_for_scurvy) |>
  kable()


## -----------------------------------------------------------------------------
#| label: tbl-scurvy-summary
#| tbl-cap: "Lind's scurvy trial, summarised: mean symptom score at day 6 (four symptoms, each 0 to 3) and number fit for duty, by treatment."
grade <- function(x) as.integer(substr(as.character(x), 1, 1))
scurvy |>
  mutate(score = grade(gum_rot_d6) + grade(skin_sores_d6) +
                 grade(weakness_of_the_knees_d6) + grade(lassitude_d6),
         fit = grade(fit_for_duty_d6)) |>
  group_by(Treatment = treatment) |>
  summarise(Patients = n(),
            `Mean symptom score (0-12)` = mean(score),
            `Fit for duty` = sum(fit)) |>
  arrange(`Mean symptom score (0-12)`) |>
  kable(digits = 1)

