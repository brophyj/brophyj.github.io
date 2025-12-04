# Required packages
library(ggplot2)
library(dplyr)
library(scales)

# Parameters of the prior (logit scale)
mu  <- -0.871   # mean on logit scale (logit(0.295) ≈ -0.871)
sdx <- 0.5      # standard deviation on logit scale

# ----- 1) Prior on the logit scale -----
# Choose a range wide enough to capture essentially all mass
x_min <- mu - 4*sdx
x_max <- mu + 4*sdx

df_logit <- data.frame(
  x = seq(x_min, x_max, length.out = 1000)
) |>
  mutate(
    density = dnorm(x, mean = mu, sd = sdx)
  )

subtitle_logit <- sprintf("X ~ Normal(%.3f, %.3f^2)", mu, sdx)

p1 <- ggplot(df_logit, aes(x = x, y = density)) +
  geom_line(color = "#2C7FB8", linewidth = 1) +
  geom_vline(xintercept = mu, linetype = "dashed", color = "grey40") +
  labs(
    title = "Baseline Prior on the Logit Scale",
    subtitle = subtitle_logit,
    x = "Log-odds (logit scale)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

# Print the logit-scale plot
p1

# ----- 2) Implied prior on the probability scale -----
# We map a fine grid on p in (0,1), avoiding the boundaries where the Jacobian blows up.
df_prob <- data.frame(
  p = seq(1e-6, 1 - 1e-6, length.out = 2000)
) |>
  mutate(
    x       = qlogis(p),                             # logit(p)
    fx      = dnorm(x, mean = mu, sd = sdx),         # density on logit scale
    jac     = 1 / (p * (1 - p)),                     # |d/dp logit(p)|
    density = fx * jac                               # transformed density on probability scale
  )

# For reference, the mean probability implied by mu
p_mean <- plogis(mu)

subtitle_prob <- sprintf(
  "P = logistic(X),  X ~ Normal(%.3f, %.3f^2)  (vertical line at p = %.3f)",
  mu, sdx, p_mean
)

p2 <- ggplot(df_prob, aes(x = p, y = density)) +
  geom_line(color = "#41AB5D", linewidth = 1) +
  geom_vline(xintercept = p_mean, linetype = "dashed", color = "grey40") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Implied Prior on the Probability Scale",
    subtitle = subtitle_prob,
    x = "Probability",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

# Print the probability-scale plot
p2