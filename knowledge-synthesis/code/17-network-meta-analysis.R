# R code for chapter: 17-network-meta-analysis.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch17
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
  library(readr)
  library(lme4)
})


## -----------------------------------------------------------------------------
#| label: plato-data
plato <- tibble::tibble(
  Region = c("US", "Eastern Europe", "Western Europe", "Asia", "Latin America", "Other"),
  tic_O = c(84, 299, 240, 90, 91, 60),   tic_N = c(707, 3820, 2725, 819, 621, 641),
  clo_O = c(67, 394, 281, 114, 104, 54), clo_N = c(706, 3825, 2704, 812, 615, 628)) |>
  mutate(yi  = log((tic_O / tic_N) / (clo_O / clo_N)),                # log RR, ticagrelor vs clopidogrel
         sei = sqrt(1 / tic_O - 1 / tic_N + 1 / clo_O - 1 / clo_N))
pool_mh <- rma.mh(measure = "RR", ai = tic_O, n1i = tic_N, ci = clo_O, n2i = clo_N, data = plato)
fmt3 <- \(x) sprintf("%.2f (%.2f to %.2f)", x[1], x[2], x[3])


## -----------------------------------------------------------------------------
#| label: plato-bayes
nn_grid <- function(yi, sei, mu_sd = 1, tau_sd = 0.5, grid = 300, n_draws = 40000) {
  g  <- expand.grid(mu = seq(-1.5, 1.5, length.out = grid), tau = seq(0, 1.5, length.out = grid))
  lp <- dnorm(g$mu, 0, mu_sd, log = TRUE) + dnorm(g$tau, 0, tau_sd, log = TRUE)
  for (i in seq_along(yi)) lp <- lp + dnorm(yi[i], g$mu, sqrt(sei[i]^2 + g$tau^2), log = TRUE)
  s <- g[sample(nrow(g), n_draws, replace = TRUE, prob = exp(lp - max(lp))), ]
  s$mu  <- s$mu + runif(n_draws, -1.5, 1.5) / (grid - 1)          # jitter within a grid cell
  s$tau <- abs(s$tau + runif(n_draws, -0.75, 0.75) / (grid - 1))
  # each region's true effect, given mu and tau (normal-normal conditional posterior)
  theta <- sapply(seq_along(yi), \(j) {
    w <- s$tau^2 / (s$tau^2 + sei[j]^2)
    rnorm(n_draws, w * yi[j] + (1 - w) * s$mu, sqrt(w) * sei[j])
  })
  list(mu = s$mu, tau = s$tau, new = rnorm(n_draws, s$mu, s$tau), theta = theta)
}
set.seed(629)
pp <- nn_grid(plato$yi, plato$sei)
q3 <- \(x) exp(quantile(x, c(0.5, 0.025, 0.975)))


## -----------------------------------------------------------------------------
#| label: plato-brms
#| eval: false
# library(brms)                               # the same model by MCMC
# priors <- c(prior(normal(0, 1), class = Intercept), prior(normal(0, 0.5), class = sd))
# fit_plato <- brm(yi | se(sei) ~ 1 + (1 | Region), data = plato, prior = priors,
#                  chains = 4, iter = 4000, seed = 629, control = list(adapt_delta = 0.99))


## -----------------------------------------------------------------------------
#| label: fig-plato
#| fig-cap: "PLATO by region: observed risk ratios (no pooling), Bayesian partially pooled estimates, and overall estimates under complete pooling, partial pooling, and for a new region."
#| fig-width: 9
#| fig-height: 5
fp <- bind_rows(
  plato |> transmute(Region, Analysis = "Observed (no pooling)",
                     est = exp(yi), lo = exp(yi - 1.96 * sei), hi = exp(yi + 1.96 * sei)),
  map_dfr(seq_len(nrow(plato)), \(j) {
    q <- q3(pp$theta[, j])
    tibble::tibble(Region = plato$Region[j], Analysis = "Partial pooling (Bayesian)",
                   est = q[1], lo = q[2], hi = q[3]) }),
  tibble::tibble(Region = "Overall",
                 Analysis = c("Complete pooling", "Partial pooling (Bayesian)", "Prediction, new region"),
                 est = c(exp(pool_mh$b), q3(pp$mu)[1], q3(pp$new)[1]),
                 lo  = c(exp(pool_mh$ci.lb), q3(pp$mu)[2], q3(pp$new)[2]),
                 hi  = c(exp(pool_mh$ci.ub), q3(pp$mu)[3], q3(pp$new)[3]))) |>
  mutate(Region = factor(Region, levels = rev(c(plato$Region, "Overall"))),
         Analysis = factor(Analysis, levels = c("Observed (no pooling)", "Partial pooling (Bayesian)",
                                                "Complete pooling", "Prediction, new region")))
ggplot(fp, aes(est, Region, colour = Analysis, shape = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  scale_x_log10(breaks = c(0.5, 0.7, 0.85, 1, 1.25, 1.5, 2)) +
  scale_colour_manual(values = c("#041E42", "#B22222", "grey45", "#E69F00"), name = NULL) +
  scale_shape_manual(values = c(16, 15, 17, 18), name = NULL) +
  labs(x = "Risk ratio, ticagrelor vs clopidogrel (log scale)", y = NULL) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: tbl-tc4
#| tbl-cap: "TC4: posterior risk ratio for the primary outcome (ticagrelor vs clopidogrel; events 50/450 vs 64/555) under each prior, with probabilities of meaningful benefit (RR < 0.9), practical equivalence (0.9 to 1.1) and meaningful harm (RR > 1.1), as published [@kutcher2025tc4]."
tibble::tribble(
  ~Prior,                    ~`RR (95% CrI)`,       ~`P(benefit)`, ~`P(equivalence)`, ~`P(harm)`,
  "Vague",                   "0.95 (0.67 to 1.35)", "37%", "42%", "21%",
  "Focused (primary)",       "1.12 (0.91 to 1.39)", "2%",  "41%", "57%",
  "Enthusiastic",            "0.85 (0.78 to 0.93)", "89%", "11%", "0%",
  "Modified enthusiastic",   "1.03 (0.73 to 1.44)", "22%", "35%", "43%",
  "Summary (earlier NMA)",   "0.99 (0.69 to 1.38)", "30%", "44%", "26%") |>
  kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: fig-dapt-direct
#| fig-cap: "Direct evidence only: Bayesian pairwise meta-analyses of the head-to-head trials (hazard ratio for MACE; T = ticagrelor, P = prasugrel, C = clopidogrel)."
#| fig-width: 7.5
#| fig-height: 3.2
dp <- read_csv("data/ch17/direct_pairwise.csv", show_col_types = FALSE)
dp |>
  mutate(label = factor(sprintf("%s (%d RCTs, %s patients)", comparison, trials, format(patients, big.mark = ",")),
                        levels = rev(sprintf("%s (%d RCTs, %s patients)", comparison, trials,
                                             format(patients, big.mark = ","))))) |>
  ggplot(aes(median_HR, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo95_HR, xmax = hi95_HR), colour = "#041E42") +
  geom_text(aes(label = sprintf("%.2f (%.2f to %.2f)", median_HR, lo95_HR, hi95_HR)), vjust = -1, size = 3.5) +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.2)) +
  labs(x = "Hazard ratio for MACE (log scale)", y = NULL)


## -----------------------------------------------------------------------------
#| label: bucher
se_ci <- \(lo, hi) (log(hi) - log(lo)) / (2 * 1.96)               # SE from a 95% interval
d_tc <- filter(dp, comparison == "T vs C"); d_pc <- filter(dp, comparison == "P vs C")
ind  <- log(d_pc$median_HR) - log(d_tc$median_HR)
se_ind <- sqrt(se_ci(d_pc$lo95_HR, d_pc$hi95_HR)^2 + se_ci(d_tc$lo95_HR, d_tc$hi95_HR)^2)
fmt3(exp(ind + c(0, -1.96, 1.96) * se_ind))


## -----------------------------------------------------------------------------
#| label: dapt-data
mace  <- read_csv("data/ch17/input_mace.csv", show_col_types = FALSE)
bleed <- read_csv("data/ch17/input_bleed.csv", show_col_types = FALSE)
mace |> filter(study %in% c("PLATO", "TRITON-TIMI 38", "ISAR-REACT 5"))   # one row per trial arm


## -----------------------------------------------------------------------------
#| label: fig-network
#| fig-cap: "The DAPT network: nodes sized by patients randomised, edges by number of trials."
#| fig-width: 7
#| fig-height: 4.5
drug_names <- c(C = "Clopidogrel", T = "Ticagrelor", P = "Prasugrel")
comps <- mace |>
  summarise(.by = study, comp = paste(sort(unique(treat), decreasing = TRUE), collapse = " vs "), N = sum(n)) |>
  mutate(comp = recode(comp, "T vs P" = "P vs T"))
nodes <- tibble::tibble(drug = c("Prasugrel", "Clopidogrel", "Ticagrelor"), x = c(0, 1, -1), y = c(1, -0.4, -0.4)) |>
  left_join(mace |> summarise(.by = treat, n = sum(n)) |> mutate(drug = drug_names[treat]), by = "drug")
edges <- tibble::tibble(from = c("Prasugrel", "Prasugrel", "Ticagrelor"),
                        to = c("Clopidogrel", "Ticagrelor", "Clopidogrel"), comp = c("P vs C", "P vs T", "T vs C")) |>
  left_join(comps |> summarise(.by = comp, trials = n(), N = sum(N)), by = "comp") |>
  left_join(select(nodes, from = drug, x1 = x, y1 = y), by = "from") |>
  left_join(select(nodes, to = drug, x2 = x, y2 = y), by = "to")
ggplot() +
  geom_segment(data = edges, aes(x1, y1, xend = x2, yend = y2, linewidth = trials), colour = "grey70") +
  geom_label(data = edges, aes((x1 + x2) / 2, (y1 + y2) / 2,
                               label = paste0(trials, " RCTs\n", format(N, big.mark = ","), " pts")),
             size = 3.5, label.size = 0, lineheight = 0.9) +
  geom_point(data = nodes, aes(x, y, size = n), shape = 21, fill = "#E69F00", colour = "#041E42") +
  geom_text(data = nodes, aes(x, y + ifelse(y > 0, 0.32, -0.32),
                              label = paste0(drug, "\nn = ", format(n, big.mark = ","))),
            fontface = "bold", size = 3.8, lineheight = 0.9) +
  scale_linewidth(range = c(2, 6), guide = "none") +
  scale_size(range = c(10, 20), guide = "none") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-0.95, 1.45)) +
  theme_void()


## -----------------------------------------------------------------------------
#| label: fig-nodesplit
#| fig-cap: "Node-splitting for prasugrel vs ticagrelor (MACE): direct, indirect and network estimates, and the ratio of direct to indirect. Shaded band: hazard ratios of 0.90 to 1.11."
#| fig-width: 7.5
#| fig-height: 3.3
ns <- read_csv("data/ch17/nodesplit_PT.csv", show_col_types = FALSE) |>
  mutate(estimate = factor(estimate, levels = rev(c("direct", "indirect", "network", "ratio_direct_indirect")),
                           labels = rev(c("Direct (5 P vs T trials)", "Indirect (15 C-controlled trials)",
                                          "Network (all 20 trials)", "Ratio direct / indirect"))))
ggplot(ns, aes(median_HR, estimate)) +
  annotate("rect", xmin = 0.90, xmax = 1.11, ymin = -Inf, ymax = Inf, fill = "#4393C3", alpha = 0.15) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = lo95_HR, xmax = hi95_HR, colour = grepl("Ratio", estimate)), show.legend = FALSE) +
  scale_colour_manual(values = c("#041E42", "#B22222")) +
  scale_x_log10(breaks = c(0.6, 0.8, 0.9, 1, 1.11, 1.3)) +
  labs(x = "Hazard ratio, prasugrel vs ticagrelor (log scale)", y = NULL)


## -----------------------------------------------------------------------------
#| label: brms-nma
#| eval: false
# library(brms)
# eff <- read_csv("data/ch17/switch_effective_counts.csv") |> filter(outcome == "MACE")
# dat <- mace |>                                     # SWITCH-SWEDEHEART enters as effective counts
#   mutate(r = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ eff$r0,
#                        study == "SWITCH-SWEDEHEART" & treat == "P" ~ eff$r1, TRUE ~ r),
#          n = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ eff$n0,
#                        study == "SWITCH-SWEDEHEART" & treat == "P" ~ eff$n1, TRUE ~ n),
#          treat = factor(treat, levels = c("C", "T", "P")), log_time = log(time_yr),
#          r_cc = ifelse(r == 0, 1, r), n_cc = ifelse(r == 0, n + 1, n))   # zero-event arms
# priors <- c(prior(normal(-2, 1), class = Intercept), prior(normal(0, 0.5), class = b),
#             prior(normal(0, 0.5), class = sd))
# fit_re <- brm(r_cc | trials(n_cc) ~ treat + offset(log_time) + (treat | study), data = dat,
#               family = binomial(link = "cloglog"), prior = priors,
#               chains = 4, iter = 4000, warmup = 2000, seed = 20260529)
# fit_fe <- update(fit_re, formula. = ~ treat + offset(log_time) + (1 | study))


## -----------------------------------------------------------------------------
#| label: draws
draws <- read_csv("data/ch17/posterior_draws.csv", show_col_types = FALSE) |>
  mutate(logHR_PT = logHR_PC - logHR_TC)
rope_lo <- log(0.90); rope_hi <- log(1.11)
pct <- \(p) ifelse(p < 0.005, "<1%", sprintf("%.0f%%", 100 * p))
summ <- \(x) tibble::tibble(HR = sprintf("%.2f (%.2f to %.2f)", exp(median(x)), exp(quantile(x, 0.025)),
                                         exp(quantile(x, 0.975))),
                            pben = mean(x < rope_lo), prope = mean(x >= rope_lo & x <= rope_hi),
                            pharm = mean(x > rope_hi))
long <- draws |>
  pivot_longer(starts_with("logHR"), names_to = "comp", values_to = "lh", names_prefix = "logHR_") |>
  mutate(comp = factor(comp, levels = c("TC", "PC", "PT"), labels = c("T vs C", "P vs C", "P vs T")))
res <- long |> group_by(model, comp) |> reframe(summ(lh))
tau_re <- draws |> filter(model == "MACE_RE")


## -----------------------------------------------------------------------------
#| label: tbl-fere
#| tbl-cap: "Fixed- and random-effects NMA of MACE: posterior hazard ratios and probabilities of a clinically meaningful benefit (HR < 0.90)."
res |> filter(model %in% c("MACE_FE", "MACE_RE")) |>
  select(model, comp, HR, pben) |>
  pivot_wider(names_from = model, values_from = c(HR, pben)) |>
  transmute(Comparison = comp, `Fixed effect: HR (95% CrI)` = HR_MACE_FE, `FE: P(HR < 0.9)` = pct(pben_MACE_FE),
            `Random effects: HR (95% CrI)` = HR_MACE_RE, `RE: P(HR < 0.9)` = pct(pben_MACE_RE)) |>
  kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: freq-check
eff  <- read_csv("data/ch17/switch_effective_counts.csv", show_col_types = FALSE)
dfq <- mace |>
  mutate(r = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ eff$r0[1],
                       study == "SWITCH-SWEDEHEART" & treat == "P" ~ eff$r1[1], TRUE ~ r),
         n = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ eff$n0[1],
                       study == "SWITCH-SWEDEHEART" & treat == "P" ~ eff$n1[1], TRUE ~ n),
         isT = as.numeric(treat == "T"), isP = as.numeric(treat == "P"), lt = log(time_yr)) |>
  filter(.by = study, sum(r) > 0)
f_fe <- glm(cbind(r, n - r) ~ 0 + factor(study) + isT + isP + offset(lt),
            family = binomial("cloglog"), data = dfq)
f_re <- glmer(cbind(r, n - r) ~ 0 + factor(study) + isT + isP + (0 + isT + isP || study) + offset(lt),
              family = binomial("cloglog"), data = dfq, control = glmerControl(optimizer = "bobyqa"))
hr_row <- \(b, se) sprintf("%.2f (%.2f to %.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
tibble::tibble(
  Model = c("Fixed effect (glm)", "Random effects (glmer, ML)"),
  `T vs C` = c(hr_row(coef(f_fe)["isT"], sqrt(vcov(f_fe)["isT", "isT"])),
               hr_row(fixef(f_re)["isT"], sqrt(vcov(f_re)["isT", "isT"]))),
  `P vs C` = c(hr_row(coef(f_fe)["isP"], sqrt(vcov(f_fe)["isP", "isP"])),
               hr_row(fixef(f_re)["isP"], sqrt(vcov(f_re)["isP", "isP"]))),
  `τ (T, P)` = c("0, 0", paste(sprintf("%.2f", attr(VarCorr(f_re)$study, "stddev")),
                               sprintf("%.2f", attr(VarCorr(f_re)$study.1, "stddev")), sep = ", "))) |>
  kable(align = "lccc")


## -----------------------------------------------------------------------------
#| label: effective-counts
adj_to_counts <- function(r0, n0, n1, or, lo, hi) {
  p0 <- r0 / n0                                       # observed control-arm risk, kept
  p1 <- p0 / (1 - p0) * or / (1 + p0 / (1 - p0) * or) # treated-arm risk implied by the adjusted OR
  se_adj <- (log(hi) - log(lo)) / (2 * qnorm(0.975))  # SE implied by the adjusted CI
  k <- (1 / (n0 * p0 * (1 - p0)) + 1 / (n1 * p1 * (1 - p1))) / se_adj^2
  N0 <- round(k * n0); N1 <- round(k * n1)
  tibble::tibble(r0 = round(N0 * p0), n0 = N0, r1 = round(N1 * p1), n1 = N1, design_effect = 1 / k)
}
sw <- bind_rows(MACE     = adj_to_counts(1113, 9444, 7651, 0.90, 0.77, 1.06),
                Bleeding = adj_to_counts(413, 9444, 7651, 0.80, 0.64, 0.99), .id = "Outcome")
sw
e <- sw[1, ]                                          # check: reproduces the adjusted estimate
predict(rma(measure = "OR", ai = e$r1, n1i = e$n1, ci = e$r0, n2i = e$n0, method = "EE"), transf = exp)


## -----------------------------------------------------------------------------
#| label: fig-mace-post
#| fig-cap: "Random-effects NMA: posterior distributions of the hazard ratio for MACE. Grey band: practical equivalence (0.90 to 1.11)."
#| fig-width: 8
#| fig-height: 4
comp_cols <- c("T vs C" = "#B22222", "P vs C" = "#1F78B4", "P vs T" = "#33A02C")
long |> filter(model == "MACE_RE") |>
  ggplot(aes(exp(lh), colour = comp, fill = comp)) +
  annotate("rect", xmin = 0.90, xmax = 1.11, ymin = 0, ymax = Inf, fill = "grey70", alpha = 0.3) +
  geom_density(alpha = 0.25, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = comp_cols, name = NULL) + scale_fill_manual(values = comp_cols, name = NULL) +
  scale_x_log10(limits = c(0.65, 1.45), breaks = c(0.7, 0.8, 0.9, 1, 1.11, 1.25, 1.4)) +
  labs(x = "Hazard ratio for MACE (log scale)", y = "Posterior density") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: tbl-nma-main
#| tbl-cap: "Random-effects NMA: hazard ratios and posterior probabilities of clinically meaningful benefit (HR < 0.90), practical equivalence (0.90 to 1.11) and meaningful harm (HR > 1.11). For bleeding, 'benefit' means less bleeding with the first-named drug."
res |> filter(model %in% c("MACE_RE", "Bleeding_RE")) |>
  mutate(Outcome = ifelse(model == "MACE_RE", "MACE", "Major bleeding"), across(c(pben, prope, pharm), pct)) |>
  arrange(desc(Outcome == "MACE")) |>
  transmute(Outcome, Comparison = comp, `HR (95% CrI)` = HR, `P(benefit)` = pben,
            `P(equivalence)` = prope, `P(harm)` = pharm) |>
  kable(align = "llcccc")


## -----------------------------------------------------------------------------
#| label: nma-probs
#| include: false
p_of <- \(m, cmp, col) res[[col]][res$model == m & res$comp == cmp]


## -----------------------------------------------------------------------------
#| label: absolute
abs_diff <- \(lh, p0) 100 * (p0 - (1 - (1 - p0)^exp(lh)))         # events prevented per 100
qs <- \(x) sprintf("%.1f (%.1f to %.1f)", median(x), quantile(x, 0.025), quantile(x, 0.975))
lhm <- draws$logHR_PC[draws$model == "MACE_RE"]; lhb <- draws$logHR_PC[draws$model == "Bleeding_RE"]
tibble::tibble(`Prasugrel vs clopidogrel, per 100 patients per year` =
                 c("MACE prevented (baseline 10%)", "Extra major bleeds (baseline 3%)",
                   "Extra major bleeds (baseline 8%, high bleeding risk)"),
               `Median (95% CrI)` = c(qs(abs_diff(lhm, 0.10)), qs(-abs_diff(lhb, 0.03)), qs(-abs_diff(lhb, 0.08)))) |>
  kable(align = "lc")


## -----------------------------------------------------------------------------
#| label: tbl-nma-sens
#| tbl-cap: "Pre-specified sensitivity analyses (MACE, random effects)."
read_csv("data/ch17/sensitivity_mace.csv", show_col_types = FALSE) |>
  mutate(hr = sprintf("%.2f (%.2f to %.2f)", median_HR, lo95_HR, hi95_HR),
         outcome = sub("Primary \\(cloglog RE\\)", "Primary analysis", outcome)) |>
  select(outcome, comparison, hr, pr_clin_benefit, pr_rope) |>
  pivot_wider(names_from = comparison, values_from = c(hr, pr_clin_benefit, pr_rope)) |>
  transmute(Analysis = outcome, `T vs C: P(equivalence)` = pct(`pr_rope_T vs C`),
            `P vs C: HR (95% CrI)` = `hr_P vs C`, `P vs C: P(benefit)` = pct(`pr_clin_benefit_P vs C`),
            `P vs T: P(benefit)` = pct(`pr_clin_benefit_P vs T`)) |>
  kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: tbl-recon
#| tbl-cap: "Trials in the 2026 frequentist NMA (JAMA Cardiology) and in the Bayesian NMA used here."
recon <- read_csv("data/ch17/reconciliation.csv", show_col_types = FALSE) |> mutate(Main = sub(" ★", "", Main))
recon |> transmute(Trial, Comparison, N = format(N, big.mark = ","), `Follow-up` = Followup,
                   `JAMA Cardiol 2026` = JAMA26, `Bayesian NMA` = Main, Note = tidyr::replace_na(Note, "")) |>
  kable()


## -----------------------------------------------------------------------------
#| label: tbl-sucra
#| tbl-cap: "Rankings for MACE from the random-effects posterior (rank 1 = lowest hazard)."
rk <- draws |> filter(model == "MACE_RE") |>
  transmute(draw, Clopidogrel = 0, Ticagrelor = logHR_TC, Prasugrel = logHR_PC) |>
  pivot_longer(-draw, names_to = "drug", values_to = "lh") |>
  mutate(.by = draw, rank = rank(lh))
rank_tab <- rk |> summarise(.by = drug, `P(best)` = mean(rank == 1), `Mean rank` = mean(rank),
                            SUCRA = (3 - mean(rank)) / 2) |> arrange(`Mean rank`)
kable(rank_tab, digits = 2, col.names = c("Treatment", "P(best)", "Mean rank", "SUCRA"))

