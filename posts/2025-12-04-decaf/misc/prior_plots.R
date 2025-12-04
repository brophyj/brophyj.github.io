# ======================================================
# Prior predictive for a two-arm Binomial trial
# and effect-measure plots with requested axis settings
# ======================================================

# Packages
library(ggplot2)
library(scales)

# -----------------------------
# 1) Prior specification (EDIT)
# -----------------------------
# Baseline (control) log-odds prior: theta0 ~ Normal(mu0, sd0^2)
mu0 <- 0.0       # logit(0.50) = 0
sd0 <- 1.5       # weakly-informative baseline prior

# Treatment effect prior on log-odds ratio: delta ~ Normal(mud, sdd^2)
mud <- -0.871    # ≈ logit(0.295) - logit(0.50), belief of ~29.5% at 50% baseline
sdd <- 0.5       # moderately informative

# Planned sample sizes (EDIT as needed)
nC <- 150
nT <- 150

# Number of prior draws (increase for smoother densities)
S <- 100000
set.seed(1)

# -----------------------------
# 2) Simulate from the priors
# -----------------------------
logistic <- function(x) 1 / (1 + exp(-x))

theta0 <- rnorm(S, mean = mu0, sd = sd0)   # baseline log-odds
delta  <- rnorm(S, mean = mud, sd = sdd)   # treatment log-OR

pC <- logistic(theta0)          # true control probability
pT <- logistic(theta0 + delta)  # true treatment probability

# -----------------------------
# 3) Prior predictive data
# -----------------------------
YC <- rbinom(S, size = nC, prob = pC)
YT <- rbinom(S, size = nT, prob = pT)

hat_pC <- YC / nC
hat_pT <- YT / nT

# -----------------------------
# 4) Effect measures
# -----------------------------
# Risk difference: RD = p_hat_T - p_hat_C (benefit shows as negative)
RD <- hat_pT - hat_pC

# Risk ratio
eps <- 1e-12
RR <- (hat_pT + eps) / (hat_pC + eps)

# Odds ratio
OR <- ((hat_pT + eps) / (1 - hat_pT + eps)) /
  ((hat_pC + eps) / (1 - hat_pC + eps))

logRR <- log(RR)
logOR <- log(OR)

# -----------------------------
# 5) Plots with requested axes
# -----------------------------

# (a) log(OR) with x in [-3, +1]
p_logOR <- ggplot(data.frame(x = logOR), aes(x = x)) +
  geom_density(fill = "#9ECAE1", color = "#3182BD", alpha = 0.3, adjust = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(xlim = c(-3, 1)) +
  labs(
    title = "Prior Predictive: log Odds ratio",
    x = "log(OR)", y = "Density"
  ) +
  theme_minimal(base_size = 12)


# (b) log(RR) with x in [-3, +1]
p_logRR <- ggplot(data.frame(x = logRR), aes(x = x)) +
  geom_density(fill = "#CFE8B9", color = "#31A354", alpha = 0.3, adjust = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(xlim = c(-3, 1)) +
  labs(
    title = "Prior Predictive: log Risk ratio",
    x = "log(RR)", y = "Density"
  ) +
  theme_minimal(base_size = 12)


# (c) Risk difference scaled per 100 patients
RD100 <- 100 * RD

# Probability that treatment is worse (RD > 0)
pr_RD_gt0 <- mean(RD > 0)

# Density and shading for RD > 0
den_RD <- density(RD100, adjust = 1.0, na.rm = TRUE)
df_RD  <- data.frame(x = den_RD$x, y = den_RD$y)

# Choose plotting limits (quantile-based for robustness)
xlim_rd <- quantile(RD100, c(0.005, 0.995), na.rm = TRUE)
# If you prefer fixed limits, e.g., -60 to +40 per 100, uncomment:
# xlim_rd <- c(-60, 40)

p_RD100 <- ggplot(df_RD, aes(x = x, y = y)) +
  geom_area(data = subset(df_RD, x > 0),
            fill = "#FDAE6B", alpha = 0.35) +
  geom_line(color = "#E6550D", linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(xlim = xlim_rd) +
  labs(
    title = "Prior Predictive: Absolute risk difference",
    subtitle = sprintf("Shaded area: P(RD > 0) = %.3f", pr_RD_gt0),
    x = "Absolute risk difference (per 100 patients)  [RD = p̂T − p̂C]",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

# Print
p_RD100

ggsave("output/prior_predictive_RD100.png", width = 6, height = 4, dpi = 300)
# -----------------------------
# 6) Quick numeric summaries
# -----------------------------
summ <- list(
  RD100_mean_median = c(mean = mean(RD100), median = median(RD100)),
  RD100_q025_q975   = quantile(RD100, c(0.025, 0.975)),
  P_RD_gt_0         = pr_RD_gt0
)
print(summ)
