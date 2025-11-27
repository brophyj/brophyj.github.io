# ---- Packages ----
library(tidyverse)
library(brms)
library(tidybayes)
library(cmdstanr)

options(mc.cores = parallel::detectCores())

# ---- Data (extracted from your figure) ----

# ---- Data (COMBINED Male + Female) ----

# ---- Data (COMBINED Male + Female, grouped) ----

df <- tribble(
  ~group,        ~study,               ~yi,    ~lower,    ~upper,
  "1500-3000",   "Aengevaeren 2017",   -1.75,  -13.54,    10.04,
  "1500-3000",   "Pavlovic 2024",      13.25,   -2.58,    29.08,
  "1500-3000",   "Defina 2019",       -23.56,  -43.53,    -3.50,
  ">3000",       "Bosscher 2023",      17.40,    4.42,    30.38,
  ">3000",       "Mohlenkamp 2008",     7.04,  -14.44,    28.52
) %>%
  mutate(
    sei = (upper - lower) / (2 * 1.96),
    study = factor(study),
    group = factor(group, levels = c("1500-3000", ">3000"))
  )

    

# ---- Priors: vague but weakly informative ----
priors <- c(
  prior(normal(0, 50), class = "Intercept"),              # mu ~ N(0, 50^2)
  prior(student_t(3, 0, 10), class = "sd", group = "study")# tau ~ t(3,0,10)
)

# ---- Fit models ----
fit_all <- brm(yi | se(sei) ~ 1 + (1 | study),
               data = df, family = gaussian(),
               prior = priors, backend = "cmdstanr", 
               iter = 4000, chains = 4, seed = 123, refresh = 0)

# Save the model object to disk
saveRDS(fit_all, file = "model/combined_all.rds")


fit_1500 <- brm(yi | se(sei) ~ 1 + (1 | study),
                data = dplyr::filter(df, group == "1500-3000"),
                family = gaussian(), prior = priors,
                backend = "cmdstanr", iter = 4000, chains = 4, seed = 123, refresh = 0)
saveRDS(fit_all, file = "model/combined_1500.rds")

fit_3000 <- brm(yi | se(sei) ~ 1 + (1 | study),
                data = dplyr::filter(df, group == ">3000"),
                family = gaussian(), prior = priors,
                backend = "cmdstanr", iter = 4000, chains = 4, refresh = 0, seed = 123)
saveRDS(fit_all, file = "model/combined_3000.rds")

fit_all <- readRDS("model/combined_all.rds")
fit_1500 <- readRDS("model/combined_1500.rds")
fit_3000 <- readRDS("model/combined_3000.rds")

# ---- Posterior draws ----
draws_all <- fit_all %>% spread_draws(b_Intercept, sd_study__Intercept)
draws_1   <- fit_1500 %>% spread_draws(b_Intercept, sd_study__Intercept)
draws_2   <- fit_3000 %>% spread_draws(b_Intercept, sd_study__Intercept)

# ---- Summaries (μ and τ) ----
summarise_mu_tau <- function(dr) {
  tibble(
    mu_median  = median(dr$b_Intercept),
    mu_low     = quantile(dr$b_Intercept, 0.025),
    mu_high    = quantile(dr$b_Intercept, 0.975),
    tau_median = median(dr$sd_study__Intercept),
    tau_low    = quantile(dr$sd_study__Intercept, 0.025),
    tau_high   = quantile(dr$sd_study__Intercept, 0.975)
  )
}
summ_all <- summarise_mu_tau(draws_all)
summ_1   <- summarise_mu_tau(draws_1)
summ_2   <- summarise_mu_tau(draws_2)

# ---- Prediction intervals (new TRUE study effect) ----
summarise_pred <- function(dr) {
  theta_new <- dr$b_Intercept + rnorm(nrow(dr)) * dr$sd_study__Intercept
  tibble(
    pred_low  = quantile(theta_new, 0.025),
    pred_high = quantile(theta_new, 0.975)
  )
}
pi_all <- summarise_pred(draws_all)
pi_1   <- summarise_pred(draws_1)
pi_2   <- summarise_pred(draws_2)

# ---- Build plotting rows with exact order ----
# Desired order from TOP to BOTTOM:
# 1) Subgroup "1500-3000": its studies (top), then "1500-3000 Pooled", then "1500-3000 Predicted"
#    then a GAP
# 2) Subgroup ">3000":     its studies, then ">3000 Pooled", then ">3000 Predicted"
#    then a GAP
# 3) Overall:              "Overall Pooled", then "Overall Predicted" (very bottom)

# Helper to make study rows
make_study_rows <- function(dsub, label_group) {
  dsub %>%
    mutate(
      lab   = paste0(study, " (", label_group, ")"),
      type  = "study",
      xmin  = lower, x = yi, xmax = upper,
      col   = "#1f77b4", shape = NA_real_
    ) %>%
    select(lab, type, xmin, x, xmax, col, shape)
}

# Helper to make pooled + predicted rows (in that order)
make_summary_rows <- function(summ, pi, label_group) {
  pooled <- tibble(
    lab   = paste0(label_group, " Pooled"),
    type  = "pooled",
    xmin  = summ$mu_low, x = summ$mu_median, xmax = summ$mu_high,
    col   = "darkred", shape = 15
  )
  predicted <- tibble(
    lab   = paste0(label_group, " Predicted"),
    type  = "predicted",
    xmin  = pi$pred_low, x = summ$mu_median, xmax = pi$pred_high,
    col   = "orange", shape = 17
  )
  dplyr::bind_rows(pooled, predicted)
}

rows_top_to_bottom <- dplyr::bind_rows(
  # Subgroup 1500–3000 block
  make_study_rows(df %>% dplyr::filter(group == "1500-3000") %>% dplyr::arrange(dplyr::desc(study)), "1500-3000"),
  make_summary_rows(summ_1, pi_1, "1500-3000"),
  tibble(lab = "", type = "gap", xmin = NA, x = NA, xmax = NA, col = NA, shape = NA_real_),
  
  # Subgroup >3000 block
  make_study_rows(df %>% dplyr::filter(group == ">3000") %>% dplyr::arrange(dplyr::desc(study)), ">3000"),
  make_summary_rows(summ_2, pi_2, ">3000"),
  tibble(lab = "", type = "gap", xmin = NA, x = NA, xmax = NA, col = NA, shape = NA_real_),
  
  # Overall at very bottom: pooled then predicted
  make_summary_rows(summ_all, pi_all, "Overall")
)

rows_top_to_bottom <- rows_top_to_bottom %>%
  mutate(row_id = dplyr::row_number()) %>%
  mutate(y = rev(row_id))

library(ggplot2)

p <- ggplot(rows_top_to_bottom, aes(y = y, x = x)) +
  geom_errorbarh(data = dplyr::filter(rows_top_to_bottom, type == "study"),
                 aes(xmin = xmin, xmax = xmax), height = 0.2,
                 color = "#1f77b4") +
  geom_point(data = dplyr::filter(rows_top_to_bottom, type == "study"),
             color = "#1f77b4", size = 2) +
  geom_errorbarh(data = dplyr::filter(rows_top_to_bottom, type == "pooled"),
                 aes(xmin = xmin, xmax = xmax), height = 0.3,
                 color = "darkred", size = 1.1) +
  geom_point(data = dplyr::filter(rows_top_to_bottom, type == "pooled"),
             aes(shape = factor(shape)), color = "darkred", size = 3) +
  geom_errorbarh(data = dplyr::filter(rows_top_to_bottom, type == "predicted"),
                 aes(xmin = xmin, xmax = xmax), height = 0.3,
                 color = "orange", size = 1.1) +
  geom_point(data = dplyr::filter(rows_top_to_bottom, type == "predicted"),
             aes(shape = factor(shape)), color = "orange", size = 3) +
  geom_vline(xintercept = 0, color = "black") +
  geom_text(data = dplyr::filter(rows_top_to_bottom, type == "pooled"),
            aes(x = xmax + 3, y = y, label = "Pooled"),
            hjust = 0, size = 3.2) +
  geom_text(data = dplyr::filter(rows_top_to_bottom, type == "predicted"),
            aes(x = xmax + 3, y = y, label = "Predicted"),
            hjust = 0, size = 3.2) +
  scale_y_continuous(
    breaks = rows_top_to_bottom$y,
    labels = rows_top_to_bottom$lab
  ) +
  scale_shape_manual(values = c(`15` = 15, `17` = 17), guide = "none") +
  labs(
    x = "Mean difference", y = NULL,
    title = "Forest plot: CACs cores in male & female athletes engaging in different volumes of exercise",
    caption = "Points and horizontal lines represent posterior medians and 95% credible intervals") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold")
  )

# Save plot
ggsave("images/combined_forest_plot.png", p, width = 7.5, height = 8.5, dpi = 300)
ggsave("images/combined_forest_plot.pdf", p, width = 7.5, height = 8.5)
