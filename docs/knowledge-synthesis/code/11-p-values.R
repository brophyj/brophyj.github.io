# R code for chapter: 11-p-values.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch11
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
suppressPackageStartupMessages(library(purrr))


## -----------------------------------------------------------------------------
#| label: fig-pval-tails
#| fig-cap: "The standard normal distribution of z under the null hypothesis. The shaded tails beyond ±1.96 hold 5% of the probability."
#| fig-width: 7
#| fig-height: 3
tibble::tibble(z = seq(-4, 4, by = 0.01), d = dnorm(z)) |>
  ggplot(aes(z, d)) +
  geom_area(data = \(x) filter(x, z <= -1.96), fill = "#56B4E9") +
  geom_area(data = \(x) filter(x, z >= 1.96), fill = "#56B4E9") +
  geom_line(linewidth = 1, colour = "#041E42") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey40") +
  annotate("text", x = c(-2.9, 2.9), y = 0.08, label = "2.5%", colour = "#0072B2") +
  labs(x = "z = estimate / SE", y = "Density under the null")


## -----------------------------------------------------------------------------
#| label: ten-pvalues
set.seed(123)
p_ten <- 2 * (1 - pnorm(abs(rnorm(10, mean = 2, sd = 1))))
round(p_ten, 4)


## -----------------------------------------------------------------------------
#| label: fig-pval-noise
#| fig-cap: "Left: the two-sided P value as a function of z. Right: P values from 2000 replications of a study whose true z is 2 (power about 50%), on a log scale."
#| fig-width: 10
#| fig-height: 3.8
#| out-width: "100%"
marks <- tibble::tibble(z = 0:3, p = 2 * pnorm(-z))
p_left <- tibble::tibble(z = seq(0, 4, 0.01), p = 2 * pnorm(-z)) |>
  ggplot(aes(z, p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#B22222") +
  geom_line(linewidth = 1.2, colour = "#041E42") +
  geom_point(data = marks, colour = "#B22222", size = 2.5) +
  geom_text(data = marks, aes(label = sprintf("%.3f", p)), nudge_x = 0.28, nudge_y = 0.04, size = 3.5) +
  labs(x = "z = estimate / SE", y = "Two-sided P value")
set.seed(2024)
z_rep <- rnorm(2000, 2, 1)
p_right <- tibble::tibble(p = 2 * pnorm(-abs(z_rep))) |>
  ggplot(aes(p, fill = p < 0.05)) +
  geom_histogram(bins = 40, boundary = 0.05) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "#B22222") +
  scale_x_log10(breaks = c(1e-4, 1e-3, 0.01, 0.05, 0.2, 1),
                labels = c("0.0001", "0.001", "0.01", "0.05", "0.2", "1")) +
  scale_fill_manual(values = c(`TRUE` = "#041E42", `FALSE` = "grey65"), guide = "none") +
  labs(x = "P value (log scale)", y = "Replications")
patchwork::wrap_plots(p_left, p_right)


## -----------------------------------------------------------------------------
#| label: fig-pval-uniform
#| fig-cap: "P values from 10 000 two-sample t tests (left) and 10 000 exact binomial tests (right), all with the null hypothesis true."
#| fig-width: 10
#| fig-height: 3.2
#| out-width: "100%"
set.seed(1234)
p_t   <- replicate(10000, t.test(rnorm(100), rnorm(100))$p.value)
p_bin <- replicate(10000, binom.test(rbinom(1, 100, 0.2), 100, p = 0.2)$p.value)
tibble::tibble(p = c(p_t, p_bin),
               test = rep(c("t test (continuous)", "Exact binomial test (discrete)"), each = 10000)) |>
  mutate(test = factor(test, levels = unique(test))) |>
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "#041E42", colour = "white") +
  geom_hline(yintercept = 500, linetype = "dashed", colour = "#B22222") +
  facet_wrap(~ test) +
  labs(x = "P value", y = "Count")


## -----------------------------------------------------------------------------
#| label: stich-svalues
s_curve <- function(hr, lo, hi, hr_test) {
  se <- (log(hi) - log(lo)) / (2 * 1.96)               # SE of log HR from the 95% CI
  p  <- 2 * pnorm(-abs(log(hr_test) - log(hr)) / se)
  tibble::tibble(hr_test, p, s = -log2(p))
}
stich <- s_curve(0.86, 0.72, 1.04, exp(seq(log(0.6), log(1.25), length.out = 400)))
s_at  <- s_curve(0.86, 0.72, 1.04, c(0.75, 1))
s_at


## -----------------------------------------------------------------------------
#| label: fig-stich
#| fig-cap: "Compatibility curves for STICH (HR 0.86, 95% CI 0.72 to 1.04). Left: P value for each hypothesised hazard ratio; the interval at P = 0.05 is the 95% confidence interval. Right: the same as S values (bits of information against each hazard ratio)."
#| fig-width: 10
#| fig-height: 3.8
#| out-width: "100%"
pts <- s_at |> mutate(lab = sprintf("HR %.2f\nP = %.2f, S = %.1f", hr_test, p, s))
p_curve <- ggplot(stich, aes(hr_test, p)) +
  geom_line(linewidth = 1.1, colour = "#041E42") +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey50") +
  geom_point(data = pts, colour = "#B22222", size = 2.5) +
  scale_x_log10(breaks = c(0.6, 0.75, 0.86, 1, 1.2)) +
  labs(x = "Hypothesised hazard ratio", y = "P value")
s_plot <- ggplot(stich, aes(hr_test, s)) +
  geom_line(linewidth = 1.1, colour = "#041E42") +
  geom_point(data = pts, colour = "#B22222", size = 2.5) +
  geom_label(data = pts, aes(y = 6.6, label = lab), size = 3.1, lineheight = 0.9, label.size = 0) +
  geom_segment(data = pts, aes(xend = hr_test, y = 5.9, yend = s + 0.3), colour = "grey50") +
  scale_x_log10(breaks = c(0.6, 0.75, 0.86, 1, 1.2)) +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = "Hypothesised hazard ratio", y = "S value (bits)")
patchwork::wrap_plots(p_curve, s_plot)


## -----------------------------------------------------------------------------
#| label: runs-sim
set.seed(629)
n <- 200
coin   <- rbinom(n, 1, 0.5)                                    # independent fair tosses
switch <- (rbinom(1, 1, 0.5) + cumsum(c(0, rbinom(n - 1, 1, 0.7)))) %% 2   # 70% switching
runs <- tibble::tibble(Sequence = c("Fair coin", "70% switching"),
                       switches = c(sum(diff(coin) != 0), sum(diff(switch) != 0)),
                       longest_run = c(max(rle(coin)$lengths), max(rle(switch)$lengths)))


## -----------------------------------------------------------------------------
#| label: fig-runs
#| fig-cap: "Two binary sequences of 200 draws, read left to right in rows of 50. Which is random?"
#| fig-width: 9
#| fig-height: 3.2
tibble::tibble(value = c(coin, switch),
               seq = rep(c("A", "B"), each = n), i = rep(seq_len(n) - 1, 2)) |>
  mutate(row = i %/% 50, col = i %% 50) |>
  ggplot(aes(col, row, fill = factor(value))) +
  geom_tile(colour = "white") +
  scale_y_reverse() +
  scale_fill_manual(values = c("grey85", "#041E42"), guide = "none") +
  facet_wrap(~ seq, ncol = 1, labeller = as_labeller(c(A = "Sequence A", B = "Sequence B"))) +
  theme_void() + theme(strip.text = element_text(face = "bold"))


## -----------------------------------------------------------------------------
#| label: tbl-simpson-data
#| tbl-cap: "Recovery with and without a drug, overall and by subgroup."
sp <- tibble::tribble(
  ~Subgroup,    ~drug_rec, ~drug_n, ~none_rec, ~none_n,
  "Stratum 1",        81,      87,      234,     270,
  "Stratum 2",       192,     263,       55,      80)
sp <- bind_rows(sp, summarise(sp, Subgroup = "All patients", across(-Subgroup, sum)))
sp |>
  transmute(Subgroup,
            Drug    = sprintf("%d/%d (%.0f%%)", drug_rec, drug_n, 100 * drug_rec / drug_n),
            `No drug` = sprintf("%d/%d (%.0f%%)", none_rec, none_n, 100 * none_rec / none_n)) |>
  kable(align = "lrr")


## -----------------------------------------------------------------------------
#| label: fig-simpson-dags
#| fig-cap: "Two causal structures consistent with the same 2×2 tables. Left: the stratifying variable is a confounder; adjust for it. Right: it is a mediator; adjusting for it blocks part of the effect."
#| fig-width: 9
#| fig-height: 2.8
dag_plot <- function(edges, nodes, title) {
  ggplot() +
    geom_curve(data = edges, aes(x, y, xend = xend, yend = yend), curvature = edges$curv[1],
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "#041E42") +
    geom_label(data = nodes, aes(x, y, label = lab), fill = "#F4F6F8", size = 3.8) +
    coord_cartesian(xlim = c(-0.4, 2.4), ylim = c(-0.4, 1.3)) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
n1 <- tibble::tibble(x = c(0, 2, 1), y = c(0, 0, 1), lab = c("Drug", "Recovery", "Sex"))
e1 <- tibble::tibble(x = c(0.25, 0.9, 1.1), y = c(0, 0.85, 0.85), xend = c(1.65, 0.15, 1.85),
                     yend = c(0, 0.12, 0.12), curv = 0)
n2 <- tibble::tibble(x = c(0, 2, 1), y = c(0, 0, 1), lab = c("Drug", "Recovery", "Blood\npressure"))
e2 <- tibble::tibble(x = c(0.25, 0.15, 1.1), y = c(0, 0.12, 0.8), xend = c(1.65, 0.8, 1.85),
                     yend = c(0, 0.85, 0.12), curv = 0)
patchwork::wrap_plots(dag_plot(e1, n1, "Confounder: stratify"),
                      dag_plot(e2, n2, "Mediator: do not stratify"))

