# R code for chapter: 03-question-to-data.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch03
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
  library(survival)
})


## -----------------------------------------------------------------------------
#| label: glp1-search
query <- '("glucagon-like peptide-1 receptor agonists"[MeSH Terms]
  OR liraglutide[Title/Abstract] OR semaglutide[Title/Abstract]
  OR dulaglutide[Title/Abstract] OR exenatide[Title/Abstract]
  OR lixisenatide[Title/Abstract] OR efpeglenatide[Title/Abstract])
  AND ("myocardial infarction"[MeSH Terms]
       OR "cardiovascular diseases"[MeSH Terms])
  AND "randomized controlled trial"[Publication Type]
  AND "humans"[MeSH Terms]'

glp1_file <- "data/glp1_mace_search.rds"
if (!file.exists(glp1_file)) {
  hits <- rentrez::entrez_search(db = "pubmed", term = query, retmax = 10)
  summ <- rentrez::entrez_summary(db = "pubmed", id = hits$ids)
  saveRDS(list(
    date  = Sys.Date(),
    query = query,
    count = hits$count,
    first_page = data.frame(
      PMID  = hits$ids,
      Year  = vapply(summ, function(x) substr(x$pubdate, 1, 4), ""),
      Title = vapply(summ, function(x) x$title, ""))),
    glp1_file)
}
glp1 <- readRDS(glp1_file)


## -----------------------------------------------------------------------------
#| label: tbl-glp1-first
#| tbl-cap: "The first page of results, unedited."
glp1$first_page |>
  mutate(Title = substr(Title, 1, 90)) |>
  kable()


## -----------------------------------------------------------------------------
#| label: itb-function
# One simulated cohort, analysed three ways. The drug has NO effect on survival.
simulate_itb <- function(n = 5000,           # patients per cohort
                         median_start = 18,  # median months from diagnosis to 1st Rx
                         p_user = 0.3,       # proportion who ever start the drug
                         median_surv = 60,   # median survival, months
                         max_fu = 120,       # administrative censoring, months
                         lm_time = 6) {      # landmark, months
  t_death <- rexp(n, log(2) / median_surv)
  user    <- runif(n) < p_user
  t_rx1   <- ifelse(user, rexp(n, log(2) / median_start), Inf)
  t_rx2   <- t_rx1 + runif(n, 1, 3)          # 2nd prescription 1-3 months later
  time    <- pmin(t_death, max_fu)
  status  <- as.integer(t_death <= max_fu)
  exposed <- t_rx2 < time                    # survived to fill a 2nd prescription
  d <- data.frame(id = seq_len(n), time, status, exposed, t_rx2)

  # 1. Paper's design: exposure fixed at baseline, clock starts at diagnosis
  hr_naive <- exp(coef(coxph(Surv(time, status) ~ exposed, data = d)))

  # 2. Landmark: exclude deaths before lm_time, restart the clock there
  dl <- subset(d, time > lm_time)
  hr_lm <- exp(coef(coxph(Surv(time - lm_time, status) ~ exposed, data = dl)))

  # 3. Aligned: unexposed until the 2nd prescription, exposed afterwards
  pre  <- data.frame(id = d$id, start = 0, stop = pmin(d$time, d$t_rx2),
                     event = ifelse(d$exposed, 0L, d$status), x = 0L)
  de   <- subset(d, exposed)
  post <- data.frame(id = de$id, start = de$t_rx2, stop = de$time,
                     event = de$status, x = 1L)
  hr_tv <- exp(coef(coxph(Surv(start, stop, event) ~ x,
                          data = rbind(pre, post), timefix = FALSE)))

  c(paper = unname(hr_naive), landmark = unname(hr_lm),
    aligned = unname(hr_tv))
}


## -----------------------------------------------------------------------------
#| label: itb-run
set.seed(2026)
scenarios <- c(6, 18, 36)
itb <- do.call(rbind, lapply(scenarios, function(m) {
  r <- t(replicate(200, simulate_itb(median_start = m)))
  data.frame(delay = m, r)
}))

itb_summary <- itb |>
  group_by(delay) |>
  summarise(across(c(paper, landmark, aligned), median))
kable(itb_summary, digits = 2,
      col.names = c("Median delay to 1st Rx (months)",
                    "As in the paper", "Landmark", "Time zero aligned"))


## -----------------------------------------------------------------------------
#| label: fig-itb
#| fig-cap: "Hazard ratios from 200 simulated cohorts per scenario in which the drug has no effect (true HR = 1, dashed line). The longer patients wait before starting treatment, the more immortal time is credited to the exposed group."
#| fig-width: 8
#| fig-height: 5
itb |>
  tidyr::pivot_longer(c(paper, landmark, aligned),
                      names_to = "analysis", values_to = "hr") |>
  mutate(analysis = factor(analysis, c("paper", "landmark", "aligned"),
                           c("As in the paper", "Landmark", "Time zero aligned")),
         delay = factor(paste("Median delay", delay, "months"),
                        levels = paste("Median delay", scenarios, "months"))) |>
  ggplot(aes(x = hr, fill = analysis)) +
  geom_density(alpha = 0.6, colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~ delay, ncol = 1) +
  scale_fill_manual(values = c("firebrick", "darkorange", "steelblue")) +
  labs(x = "Estimated hazard ratio", y = NULL, fill = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_blank())


## -----------------------------------------------------------------------------
#| label: inception-bayes
set.seed(123)
n_draws <- 200000
p_conv  <- rbeta(n_draws, 1 + 10, 1 + 62 - 10)   # conventional CPR: 10/62
p_ecpr  <- rbeta(n_draws, 1 + 14, 1 + 70 - 14)   # extracorporeal CPR: 14/70
or_post <- (p_ecpr / (1 - p_ecpr)) / (p_conv / (1 - p_conv))

inception <- c(
  median    = median(or_post),
  lower95   = unname(quantile(or_post, 0.025)),
  upper95   = unname(quantile(or_post, 0.975)),
  p_benefit = mean(or_post > 1),                     # any benefit
  p_clin    = mean(or_post > 1.1),                   # clinically meaningful benefit
  p_equiv   = mean(or_post > 0.9 & or_post < 1.1))   # practical equivalence
round(inception, 2)


## -----------------------------------------------------------------------------
#| label: fig-inception
#| fig-cap: "Posterior distribution of the odds ratio for a favourable outcome with extracorporeal versus conventional CPR in the INCEPTION trial (flat priors). The vertical lines bound a range of practical equivalence (0.9 to 1.1); the shaded areas show the probability of meaningful harm (left) and meaningful benefit (right)."
#| fig-width: 7
#| fig-height: 4
dens <- density(or_post, from = 0, to = 6, n = 2048)
dd   <- data.frame(x = dens$x, y = dens$y)
ggplot(dd, aes(x, y)) +
  geom_area(data = subset(dd, x <= 0.9), fill = "#f4a6b0") +
  geom_area(data = subset(dd, x >= 1.1), fill = "#3fbfbf") +
  geom_line(colour = "grey20") +
  geom_vline(xintercept = c(0.9, 1.1)) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(x = "Odds ratio (extracorporeal vs conventional CPR)", y = "Density")


## -----------------------------------------------------------------------------
#| label: inception-brms
#| eval: false
# library(brms)
# d_inc <- data.frame(arm = factor(c("conventional", "extracorporeal")),
#                     y = c(10, 14), n = c(62, 70))
# fit_inc <- brm(y | trials(n) ~ 0 + arm, data = d_inc,
#                family = binomial(link = "identity"),
#                prior  = prior(beta(1, 1), class = b, lb = 0, ub = 1),
#                seed = 123, refresh = 0)
# draws <- as_draws_df(fit_inc)
# p_c <- draws$b_armconventional; p_e <- draws$b_armextracorporeal
# quantile((p_e / (1 - p_e)) / (p_c / (1 - p_c)), c(0.025, 0.5, 0.975))


## -----------------------------------------------------------------------------
#| label: extraction-schema
# Fixed once, in the protocol
schema <- tribble(
  ~key,             ~Field,               ~Type,
  "study_id",       "Study ID",           "Text",
  "first_author",   "First author",       "Text",
  "year",           "Year",               "Integer",
  "design",         "Design",             "Categorical (RCT / cohort / ...)",
  "n_rand",         "N randomised",       "Integer",
  "mean_age",       "Mean age (SD)",      "Text",
  "pct_female",     "% female",           "Numeric",
  "intervention",   "Intervention",       "Text",
  "dose_duration",  "Dose / duration",    "Text",
  "comparator",     "Comparator",         "Text",
  "outcome",        "Primary outcome",    "Text",
  "followup_mo",    "Follow-up (months)", "Numeric",
  "events_I",       "Events, intervention", "Integer",
  "n_I",            "Total, intervention",  "Integer",
  "events_C",       "Events, control",      "Integer",
  "n_C",            "Total, control",       "Integer",
  "effect",         "Reported estimate",  "Numeric",
  "ci_lo",          "95% CI lower",       "Numeric",
  "ci_hi",          "95% CI upper",       "Numeric",
  "effect_measure", "Effect measure",     "Categorical (RR / OR / HR / RD)",
  "analysis_pop",   "Analysis population","Categorical (ITT / mITT / PP)",
  "registration",   "Registration",       "Text",
  "funding",        "Funding",            "Text",
  "notes",          "Notes / queries",    "Free text"
)

# One study, as a reviewer would fill it in
emperor_r <- list(
  study_id = "EMPEROR-Reduced", first_author = "Packer", year = 2020,
  design = "RCT", n_rand = 3730, mean_age = "66.9 (10.8)", pct_female = 24,
  intervention = "Empagliflozin", dose_duration = "10 mg daily",
  comparator = "Placebo", outcome = "CV death or HF hospitalisation",
  followup_mo = 16, events_I = 361, n_I = 1863, events_C = 462, n_C = 1867,
  effect = 0.75, ci_lo = 0.65, ci_hi = 0.86, effect_measure = "HR",
  analysis_pop = "ITT", registration = "NCT03057977",
  funding = "Boehringer Ingelheim; Eli Lilly", notes = "")


## -----------------------------------------------------------------------------
#| label: tbl-extraction
#| tbl-cap: "A minimal extraction form, filled in for one trial."
schema |>
  mutate(Value = vapply(key, function(k) as.character(emperor_r[[k]]), "")) |>
  select(Field, Type, Value) |>
  kable()


## -----------------------------------------------------------------------------
#| label: binary-first-principles
events_I <- 361; n_I <- 1863   # empagliflozin
events_C <- 462; n_C <- 1867   # placebo

p_I <- events_I / n_I
p_C <- events_C / n_C
RR  <- p_I / p_C
OR  <- (p_I / (1 - p_I)) / (p_C / (1 - p_C))
RD  <- p_I - p_C

# Standard errors on the log scale, used by meta-analysis
se_log_RR <- sqrt(1/events_I - 1/n_I + 1/events_C - 1/n_C)
se_log_OR <- sqrt(1/events_I + 1/(n_I - events_I) + 1/events_C + 1/(n_C - events_C))
se_RD     <- sqrt(p_I * (1 - p_I) / n_I + p_C * (1 - p_C) / n_C)

ci <- function(est, se, log = TRUE) {
  if (log) exp(log(est) + c(-1, 1) * 1.96 * se) else est + c(-1, 1) * 1.96 * se
}
rbind(RR = c(RR, ci(RR, se_log_RR)),
      OR = c(OR, ci(OR, se_log_OR)),
      RD = c(RD, ci(RD, se_RD, log = FALSE))) |>
  round(3)


## -----------------------------------------------------------------------------
#| label: epir
tab <- matrix(c(events_I, n_I - events_I,
                events_C, n_C - events_C),
              nrow = 2, byrow = TRUE,
              dimnames = list(Group = c("Empagliflozin", "Placebo"),
                              Outcome = c("Event", "No event")))
epiR::epi.2by2(as.table(tab), method = "cohort.count")


## -----------------------------------------------------------------------------
#| label: sd-backcalc
# (a) One arm reports its mean with a 95% CI: 24.1 (21.7 to 26.5), n = 60
n1 <- 60
se_arm <- (26.5 - 21.7) / (2 * qt(0.975, df = n1 - 1))
sd_arm <- se_arm * sqrt(n1)

# (b) A paper reports only the mean difference between arms, -5.0,
#     with 95% CI -8.5 to -1.5, and 60 participants per arm
n_I2 <- 60; n_C2 <- 60
df   <- n_I2 + n_C2 - 2
se_diff_ci <- (-1.5 - -8.5) / (2 * qt(0.975, df))
sd_pooled_ci <- se_diff_ci / sqrt(1/n_I2 + 1/n_C2)

# (c) ... or only the mean difference and its P value, P = 0.006
se_diff_p <- abs(-5.0) / qt(1 - 0.006/2, df)
sd_pooled_p <- se_diff_p / sqrt(1/n_I2 + 1/n_C2)

round(c(sd_one_arm = sd_arm,
        sd_pooled_from_ci = sd_pooled_ci,
        sd_pooled_from_p  = sd_pooled_p), 2)


## -----------------------------------------------------------------------------
#| label: escalc-smd
cont_dat <- data.frame(
  study = "Example trial",
  m1i = 24.1, sd1i = sd_arm, n1i = 60,   # intervention arm (SD recovered above)
  m2i = 28.6, sd2i = 9.2,    n2i = 60)   # control arm (SD reported)

escalc(measure = "SMD",
       m1i = m1i, sd1i = sd1i, n1i = n1i,
       m2i = m2i, sd2i = sd2i, n2i = n2i,
       data = cont_dat)


## -----------------------------------------------------------------------------
#| label: living-pipeline
#| eval: false
# # 1. Search for records added since the last update
# last_update <- as.Date(readLines("data/last_update.txt"))
# new_records <- fetch_rcts(mesh_terms = "heart failure",
#                           start = format(last_update, "%Y/%m/%d"),
#                           end   = format(Sys.Date(), "%Y/%m/%d"))
# 
# # 2. Screen new records against the inclusion criteria (in duplicate)
# # 3. Extract eligible studies (in duplicate) into new_extraction
# 
# # 4. Append to the existing extraction sheet, validate, and save
# extraction <- bind_rows(read.csv("data/extraction.csv"), new_extraction)
# validate_extraction(extraction)
# write.csv(extraction, "data/extraction.csv", row.names = FALSE)
# writeLines(format(Sys.Date()), "data/last_update.txt")
# 
# # 5. Rerun the analysis
# source("analysis/meta_analysis.R")


## -----------------------------------------------------------------------------
#| label: tbl-sglt2
#| tbl-cap: "Extracted data: cardiovascular death or hospitalisation for heart failure in four SGLT2 inhibitor trials."
sglt2 <- data.frame(
  study    = c("EMPEROR-Reduced 2020", "DAPA-HF 2019",
               "EMPEROR-Preserved 2021", "DELIVER 2022"),
  events_I = c(361, 386, 415, 512),
  n_I      = c(1863, 2373, 2997, 3131),
  events_C = c(462, 502, 511, 610),
  n_C      = c(1867, 2371, 2991, 3132))
kable(sglt2, col.names = c("Trial", "Events (SGLT2i)", "N (SGLT2i)",
                           "Events (placebo)", "N (placebo)"))


## -----------------------------------------------------------------------------
#| label: fig-sglt2-forest
#| fig-cap: "Random-effects meta-analysis of four SGLT2 inhibitor trials."
#| fig-width: 8
#| fig-height: 3.5
sglt2_es <- escalc(measure = "RR",
                   ai = events_I, n1i = n_I, ci = events_C, n2i = n_C,
                   data = sglt2, slab = study)
res_sglt2 <- rma(yi, vi, data = sglt2_es, method = "REML")
forest(res_sglt2, atransf = exp, at = log(c(0.6, 0.8, 1, 1.25)),
       xlab = "Risk ratio (log scale)", header = "Trial")


## -----------------------------------------------------------------------------
#| label: validate-function
validate_extraction <- function(df) {
  issues <- list()
  # Events cannot exceed totals, or be negative
  bad_I <- with(df, events_I > n_I | events_I < 0)
  bad_C <- with(df, events_C > n_C | events_C < 0)
  if (any(bad_I, na.rm = TRUE)) issues$impossible_counts_I <- df$study[which(bad_I)]
  if (any(bad_C, na.rm = TRUE)) issues$impossible_counts_C <- df$study[which(bad_C)]
  # Log risk ratio and its standard error
  log_rr <- with(df, log((events_I / n_I) / (events_C / n_C)))
  se     <- with(df, sqrt(1/events_I - 1/n_I + 1/events_C - 1/n_C))
  if (any(!is.finite(se) | se <= 0, na.rm = TRUE))
    issues$invalid_se <- df$study[which(!is.finite(se) | se <= 0)]
  # Implausibly large effects (|log RR| > 3, i.e. RR < 0.05 or > 20)
  if (any(abs(log_rr) > 3, na.rm = TRUE))
    issues$large_effect <- df$study[which(abs(log_rr) > 3)]
  # Duplicate studies
  if (any(duplicated(df$study)))
    issues$duplicates <- df$study[duplicated(df$study)]

  if (length(issues) == 0) message("All checks passed.") else print(issues)
  invisible(issues)
}

validate_extraction(sglt2)


## -----------------------------------------------------------------------------
#| label: validate-broken
broken <- rbind(sglt2, sglt2[2, ])        # DAPA-HF entered twice
broken$events_C[4] <- 6100                # 610 mistyped as 6100
validate_extraction(broken)

