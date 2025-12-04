# ================== Packages ==================

library(readr); library(dplyr); library(tidyr); library(ggplot2)
library(survival); library(survminer)
library(IPDfromKM)
library(brms)
library(cmdstanr)

# Use cmdstanr backend (faster, fewer divergences)
cmdstanr::set_cmdstan_path()
options(mc.cores = max(1, parallel::detectCores() - 1))

# ================== 1) Read & tidy your CSV ==================

df <- read_csv("input/decaf_incident_data.csv",
               col_names = c("coffee_time","coffee_inc","decaf_time","decaf_inc"),
               skip = 2, show_col_types = FALSE)

coffee_xy <- df %>% filter(!is.na(coffee_time) & !is.na(coffee_inc)) %>%
  transmute(time = coffee_time, inc = coffee_inc) %>%
  arrange(time) %>% distinct(time, .keep_all = TRUE)

decaf_xy  <- df %>% filter(!is.na(decaf_time)  & !is.na(decaf_inc)) %>%
  transmute(time = decaf_time,  inc = decaf_inc)  %>%
  arrange(time) %>% distinct(time, .keep_all = TRUE)

# Truncate raw digitized series at 180 days 
coffee_xy <- coffee_xy %>% filter(time <= 180)
decaf_xy  <- decaf_xy  %>% filter(time <= 180)

# ================== 2) Incidence(%) -> Survival(%) for reconstruction ==================
coffee_km <- coffee_xy %>% transmute(time, surv = pmax(0, 100 - inc))
decaf_km  <- decaf_xy  %>% transmute(time, surv = pmax(0, 100 - inc))
if (min(coffee_km$time) > 0) coffee_km <- bind_rows(tibble(time = 0, surv = 100), coffee_km)
if (min(decaf_km$time)  > 0) decaf_km  <- bind_rows(tibble(time = 0, surv = 100),  decaf_km)

# ================== 3) At-risk tables (yours) ==================
trisk <- c(0, 30, 60, 90, 120, 150, 180)
nrisk_coffee <- c(100, 79, 70, 61, 57, 56, 52)
nrisk_decaf  <- c(100, 64, 53, 47, 45, 38, 36)

# ================== 4) Reconstruct IPD (Guyot method via IPDfromKM) ==================

prep_coffee <- preprocess(as.matrix(coffee_km), trisk = trisk, nrisk = nrisk_coffee, maxy = 100)
prep_decaf  <- preprocess(as.matrix(decaf_km),  trisk = trisk, nrisk = nrisk_decaf,  maxy = 100)

ipd_coffee <- getIPD(prep_coffee, armID = 1)$IPD  # coffee arm
ipd_decaf  <- getIPD(prep_decaf,  armID = 0)$IPD  # decaf arm

ipd <- bind_rows(
  data.frame(time = ipd_coffee$time, status = ipd_coffee$status, arm_coffee = 1L),
  data.frame(time = ipd_decaf$time,  status = ipd_decaf$status,  arm_coffee = 0L)
)

# ================== 5) Truncate ANALYSIS at 180 days ==================
ipd <- ipd %>%
  mutate(time_orig = time, status_orig = status,
         time = pmin(time, 180),
         status = ifelse(time_orig > 180, 0L, status_orig)) %>%
  select(-time_orig, -status_orig)

# ================== 6) Frequentist Cox (authors’ coding: coffee=1, decaf=0) ==================
fit_cox <- coxph(Surv(time, status) ~ arm_coffee, data = ipd)
s_cox   <- summary(fit_cox)
hr      <- exp(coef(fit_cox)["arm_coffee"])
ci_95   <- exp(confint(fit_cox)["arm_coffee", ])
p_value <- s_cox$coefficients["arm_coffee", "Pr(>|z|)"]
cat(sprintf("\nFrequentist Cox (coffee vs decaf): HR = %.2f (95%% CI %.2f–%.2f); p = %.3f\n",
            hr, ci_95[1], ci_95[2], p_value))

# ================== 7) KM plot (truncated at 180 days) ==================
km_fit <- survfit(Surv(time, status) ~ arm_coffee, data = ipd)
p_km <- ggsurvplot(
  km_fit, data = ipd, conf.int = FALSE, risk.table = TRUE,
  xlim = c(0, 180),
  legend.labs = c("Decaf (control)", "Coffee (intervention)"),
  legend.title = "", palette = c("#D95F02","#1B9E77")
)

ggsave("output/km_plot_coffee_vs_decaf.png", p_km$plot, width = 7, height = 5, dpi = 300)

# ================== 8) Cumulative incidence overlay: digitized vs reconstructed ==================
# Reconstructed cumulative incidence from KM (in %)
km_df <- survminer::surv_summary(km_fit) %>%
  mutate(cuminc = 100 * (1 - surv),
         arm = ifelse(strata == "arm_coffee=0", "Decaf (KM)", "Coffee (KM)")) %>%
  select(time, cuminc, arm)

# Digitized cumulative incidence (already in %), labeled
dig_df <- bind_rows(
  coffee_xy %>% transmute(time, cuminc = inc, arm = "Coffee (digitized)"),
  decaf_xy  %>% transmute(time, cuminc = inc, arm = "Decaf (digitized)")
)

p_ci <- ggplot() +
  geom_step(data = km_df, aes(x = time, y = cuminc, color = arm), linewidth = 1) +
  geom_line(data = dig_df, aes(x = time, y = cuminc, color = arm), linetype = "dashed") +
  scale_color_manual(values = c("Coffee (KM)"="#1B9E77","Decaf (KM)"="#D95F02",
                                "Coffee (digitized)"="#1B9E77","Decaf (digitized)"="#D95F02")) +
  coord_cartesian(xlim = c(0,180), ylim = c(0,100)) +
  labs(title = "Cumulative incidence (reconstructed KM vs digitized points)",
       x = "Days", y = "Cumulative incidence (%)",
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(p_ci)
ggsave("output/cumulative_incidence_overlay.png", p_ci, width = 7, height = 5, dpi = 300)

# ================== 9) Bayesian Cox in brms (coffee vs decaf) ==================
# Prior: belief that decaf is better -> coffee harmful a priori.
# Encode as Normal(+0.871, 0.5) on the coffee coefficient (log-HR).
ipd_brm <- ipd %>% mutate(cens = 1L - status)  # brms uses 1=censored

prior_cox <- c(
  set_prior("normal(0.871, 0.5)", class = "b", coef = "arm_coffee")
)

fit_brm_cox <- brm(
  time | cens(cens) ~ arm_coffee,
  data    = ipd_brm,
  family  = cox(),
  prior   = prior_cox,
  backend = "cmdstanr",
  chains  = 4, iter = 2500, warmup = 1000, seed = 20251114,
  control = list(adapt_delta = 0.999, max_treedepth = 13)
)

print(fit_brm_cox)

# ---- Finish: summarize posterior HR and a simple probability statement ----
post_cox <- posterior_summary(fit_brm_cox, variable = "b_arm_coffee")
HR_post  <- exp(unlist(post_cox[1, c("Estimate","Q2.5","Q97.5")]))

cat(sprintf(
  "\nBayesian Cox (coffee vs decaf): HR = %.2f (95%% CrI %.2f–%.2f)\n",
  HR_post[1], HR_post[2], HR_post[3]
))

# Optional: posterior probability that coffee is beneficial (HR < 1)
draws_b <- as_draws_df(fit_brm_cox)[["b_arm_coffee"]]
cat(sprintf("Pr(HR < 1 | data, prior) = %.3f\n", mean(draws_b < 0)))