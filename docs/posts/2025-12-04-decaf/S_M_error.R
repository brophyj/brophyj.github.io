# ----------------------------------------
# Required packages
# ----------------------------------------
library(retrodesign)
library(ggplot2)
library(patchwork)

# ----------------------------------------
# Calculate SE from observed effect and CrI
# ----------------------------------------
observed_effect <- -0.076     # Observed RD (Caf - Decaf)
lower <- -0.195               # Lower bound of 95% CrI
upper <- 0.044                # Upper bound of 95% CrI

# Approximate SE from 95% CrI
se <- (upper - lower) / (2 * 1.96)

# z-score and p-value (frequentist approximation)
z <- observed_effect / se
p_value <- 2 * (1 - pnorm(abs(z)))

# ----------------------------------------
# Assumed effect from original design
# ----------------------------------------
assumed_effect <- 0.20  # DECAF assumed a large benefit for decaf

# ----------------------------------------
# Compute retrodesign metrics for observed and assumed effects
# ----------------------------------------
retro_observed <- retro_design_closed_form(observed_effect, se)
retro_assumed <- retro_design_closed_form(assumed_effect, se)

power_observed <- retro_observed$power
power_assumed <- retro_assumed$power
type_s_observed <- retro_observed$type_s
type_m_observed <- retro_observed$type_m

# ----------------------------------------
# Print numeric summary
# ----------------------------------------
cat("Observed Effect =", observed_effect, "\n")
cat("SE =", round(se, 4), "\n")
cat("Approx p-value =", round(p_value, 4), "\n")
cat("Power (Observed) =", round(power_observed, 3), "\n")
cat("Type S Error =", round(type_s_observed, 3), "\n")
cat("Type M Error =", round(type_m_observed, 3), "\n\n")

cat("Assumed Effect =", assumed_effect, "\n")
cat("Power (Assumed) =", round(power_assumed, 3), "\n\n")


# ----------------------------------------
# Range of true effects for plotting
# ----------------------------------------
true_effects <- seq(-0.20, 0.20, length.out = 200)
retro <- retro_design_closed_form(true_effects, se)

# Build data frame
df <- data.frame(
  TrueEffect = true_effects,
  Power = retro$power,
  TypeS = retro$type_s,
  Exaggeration = retro$type_m
)

# ----------------------------------------
# Custom theme for polished look
# ----------------------------------------
custom_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# ----------------------------------------
# Panel 1: Type S Error
# ----------------------------------------
p1 <- ggplot(df, aes(x = TrueEffect, y = TypeS)) +
  geom_line(color = "darkred", size = 1.5) +
  geom_vline(xintercept = assumed_effect, linetype = "dashed", color = "green", linewidth = 1) +
  geom_vline(xintercept = observed_effect, linetype = "dotted", color = "black", linewidth = 1) +
  annotate("text", x = .1, y = 0.42,
           label = paste("Assumed +0.20\nPower =", round(power_assumed, 2)),
           color = "green", hjust = -0.1, size = 5) +
  annotate("text", x = observed_effect, y = 0.42,
           label = paste("Observed =", observed_effect,
                         "\nPower =", round(power_observed, 2),
                         "\nType S =", round(type_s_observed, 2)),
           color = "black", hjust = 1.1, size = 5) +
  labs(title = "Type S Error vs True Effect Size",
       subtitle = paste("Observed =", observed_effect, "| SE =", round(se, 4)),
       x = "True Effect Size (Caf - Decaf)",
       y = "Type S Error Probability") +
  custom_theme

# ----------------------------------------
# Panel 2: Exaggeration Ratio
# ----------------------------------------
p2 <- ggplot(df, aes(x = TrueEffect, y = Exaggeration)) +
  geom_line(color = "steelblue", size = 1.5) +
  geom_vline(xintercept = assumed_effect, linetype = "dashed", color = "green", size = 1) +
  geom_vline(xintercept = observed_effect, linetype = "dotted", color = "black", size = 1) +
  annotate("text", x = .1, y = max(df$Exaggeration)*0.85,
           label = paste("Assumed +0.20\nPower =", round(power_assumed, 2)),
           color = "green", hjust = -0.1, size = 5) +
  annotate("text", x = 0, y = 30,
           label = paste("Observed =", observed_effect,
                         "\nPower =", round(power_observed, 2),
                         "\nType M =", round(type_m_observed, 2)),
           color = "black", hjust = 1.1, size = 5) +
  xlim(-0.10, 0.10) +
  ylim(0,50) +
  labs(title = "Exaggeration Ratio vs \nTrue Effect Size",
       x = "True Effect Size (Caf - Decaf)",
       y = "Exaggeration Ratio") +
  custom_theme

# ----------------------------------------
# Panel 3: Power
# ----------------------------------------
p3 <- ggplot(df, aes(x = TrueEffect, y = Power)) +
  geom_line(color = "darkgreen", size = 1.5) +
  geom_vline(xintercept = assumed_effect, linetype = "dashed", color = "green", size = 1) +
  geom_vline(xintercept = observed_effect, linetype = "dotted", color = "black", size = 1) +
  annotate("text", x = .05, y = 0.8,
           label = paste("Assumed +0.20\nPower =", round(power_assumed, 2)),
           color = "green", hjust = -0.1, size = 5) +
  annotate("text", x = 0, y = 0.8,
           label = paste("Observed =", observed_effect,
                         "\nPower =", round(power_observed, 2)),
           color = "black", hjust = 1.1, size = 5) +
  labs(title = "Power vs \nTrue Effect Size",
       x = "True Effect Size (Caf - Decaf)",
       y = "Power") +
  custom_theme

# ----------------------------------------
# Save combined plots
# ----------------------------------------

combined_plot <- p1 / p2 / p3
combined_plotMP <- p3 + p2
ggsave("output/decaf_retrodesign_panels_polished.png", combined_plot, width = 10, height = 14, dpi = 300)
ggsave("output/decaf_retrodesign_panels_polished_MP.png", combined_plotMP, width = 10, height = 7, dpi = 300)

