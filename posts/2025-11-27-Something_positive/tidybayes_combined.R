
# ============================================================
# Combined Male + Female Bayesian Meta-Analysis (Full Script)
# ============================================================

library(tidyverse)
library(brms)
library(tidybayes)
library(cmdstanr)
library(patchwork)

options(mc.cores = parallel::detectCores())

# ---- Load Data ----
df <- read_csv("combined_meta.csv") %>%
  mutate(
    study = factor(study),
    group = factor(group, levels = c("1500-3000", ">3000")),
    sex = factor(sex, levels = c("Female", "Male")) # Female as reference
  )

# ---- Priors ----
priors <- c(
  prior(normal(0, 50), class = "Intercept"),
  prior(student_t(3, 0, 10), class = "sd", group = "study")
)

# ---- Fit Model ----
if (file.exists("fit_combined.rds")) {
  fit_all <- readRDS("fit_combined.rds")
} else {
  fit_all <- brm(
    yi | se((upper - lower)/(2*1.96)) ~ sex + (1|study),
    data = df, family = gaussian(),
    prior = priors, backend = "cmdstanr",
    iter = 4000, chains = 4, seed = 123,
    control = list(adapt_delta = 0.95)
  )
  saveRDS(fit_all, "fit_combined.rds")
}

# ---- Posterior Summaries ----
summary(fit_all)

# Extract pooled effects
draws <- fit_all %>%
  spread_draws(b_Intercept, b_sexMale, sd_study__Intercept)

# Compute predicted means for Female and Male
preds <- draws %>%
  mutate(
    Female = b_Intercept,
    Male = b_Intercept + b_sexMale
  ) %>%
  pivot_longer(cols = c(Female, Male), names_to = "Sex", values_to = "Effect")

# Summary table
summary_table <- preds %>%
  group_by(Sex) %>%
  mean_qi(Effect, .width = c(0.95, 0.80))
print(summary_table)

# ---- Forest Plot: Study-level effects ----
study_effects <- fit_all %>%
  spread_draws(b_Intercept, r_study[study,Intercept]) %>%
  mutate(study_effect = b_Intercept + r_study)

forest_plot <- study_effects %>%
  ggplot(aes(y = study, x = study_effect)) +
  stat_halfeye(.width = c(0.95, 0.80), fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Study-Level Effects", x = "Effect Size", y = "Study") +
  theme_minimal()

# ---- Pooled Effect Plot ----
pooled_plot <- preds %>%
  ggplot(aes(x = Effect, y = Sex, fill = Sex)) +
  stat_halfeye(.width = c(0.95, 0.80), alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Pooled Effects by Sex", x = "Effect Size", y = NULL) +
  theme_minimal()

# ---- Combine Plots ----
combined_plot <- pooled_plot / forest_plot + plot_annotation(title = "Combined Meta-Analysis")
print(combined_plot)

# ---- Posterior Predictive Check ----
