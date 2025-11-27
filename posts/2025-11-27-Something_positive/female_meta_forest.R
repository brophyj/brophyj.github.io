# ---- Packages ----
library(tidyverse)
library(brms)
library(tidybayes)
library(cmdstanr)

options(mc.cores = parallel::detectCores())

# ---- Female data ----
df <- read_csv("female_meta.csv") %>%
  mutate(study = factor(study),
         group = factor(group, levels = c("1500-3000", ">3000")))

# ---- Priors ----
priors <- c(
  prior(normal(0, 50), class = "Intercept"),
  prior(student_t(3, 0, 10), class = "sd", group = "study")
)

# ---- Fit models ----
fit_all <- brm(yi | se((upper - lower)/(2*1.96)) ~ 1 + (1|study),
               data = df, family = gaussian(), prior = priors,
               backend = "cmdstanr", iter = 4000, chains = 4, seed = 123)

fit_1500 <- brm(yi | se((upper - lower)/(2*1.96)) ~ 1 + (1|study),
                data = filter(df, group == "1500-3000"), family = gaussian(),
                prior = priors, backend = "cmdstanr", iter = 4000, chains = 4, seed = 123)

fit_3000 <- brm(yi | se((upper - lower)/(2*1.96)) ~ 1 + (1|study),
                data = filter(df, group == ">3000"), family = gaussian(),
                prior = priors, backend = "cmdstanr", iter = 4000, chains = 4, seed = 123)

# ---- Summaries and plot ----
# (Same plotting logic as previous script: studies first, then subgroup pooled, then predicted, overall at bottom)
