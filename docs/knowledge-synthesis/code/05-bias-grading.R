# R code for chapter: 05-bias-grading.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch05
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
  library(tidyr)
  library(robvis)
})


## -----------------------------------------------------------------------------
#| label: fig-phrases
#| fig-cap: "Probability assigned to 19 verbal phrases by more than 5000 respondents: median (circle), interquartile range (coloured bar, darker red for wider ranges) and 10th to 90th percentiles (grey). Adapted from J. Lopes' TidyTuesday visualisation 'The Scale of Chaos' [@lopes2026chaos]; data from the CAPphrase dataset [@kucharski2026capphrase]."
#| fig-width: 8
#| fig-height: 6
#| fig-alt: "Horizontal interval chart of 19 probability phrases ordered by median. Phrases such as Realistic Possibility, Could Happen, Might Happen and May Happen have very wide ranges; Almost No Chance and Will Happen are narrow."
phrases <- readr::read_csv("data/ch05/absolute_judgements.csv", show_col_types = FALSE)

phrase_summary <- phrases |>
  group_by(term) |>
  summarise(q10 = quantile(probability, 0.10), q25 = quantile(probability, 0.25),
            median = median(probability),
            q75 = quantile(probability, 0.75), q90 = quantile(probability, 0.90)) |>
  mutate(iqr = q75 - q25) |>
  arrange(median, q25) |>
  mutate(term = factor(term, levels = term))

ggplot(phrase_summary, aes(y = term)) +
  geom_linerange(aes(xmin = q10, xmax = q90), linewidth = 3, colour = "grey85") +
  geom_linerange(aes(xmin = q25, xmax = q75, colour = iqr), linewidth = 3) +
  geom_point(aes(x = median), shape = 21, fill = "white", size = 2.2) +
  scale_colour_gradient(low = "#2C4A7E", high = "#C8102E",
                        name = "Interquartile range (percentage points)") +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  labs(x = "Probability assigned", y = NULL) +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: fig-bias-stages
#| fig-cap: "Where bias enters a trial, and the main tools that detect it."
#| fig-width: 9
#| fig-height: 2.2
#| fig-alt: "Four boxes in sequence: design, conduct, analysis and reporting, with examples of bias and the relevant tools under each."
stages <- tibble::tibble(
  x     = 1:4,
  stage = c("Design", "Conduct", "Analysis", "Reporting"),
  what  = c("sequence generation\nallocation concealment",
            "blinding, co-interventions\nlosses to follow-up",
            "not as randomised\nchoice among analyses",
            "selective outcomes\nnon-publication"),
  tool  = c("RoB 2 domain 1", "RoB 2 domains 2-4",
            "RoB 2 domains 2 and 5", "RoB 2 domain 5; ROB-ME;\nregisters, funnel plots"))
ggplot(stages) +
  geom_tile(aes(x, 2, fill = stage), width = 0.92, height = 0.9, show.legend = FALSE) +
  geom_text(aes(x, 2.25, label = stage), colour = "white", fontface = "bold", size = 5) +
  geom_text(aes(x, 1.85, label = what), colour = "white", size = 3.3, lineheight = 0.9) +
  geom_text(aes(x, 1.38, label = tool), size = 3.3, lineheight = 0.9) +
  annotate("segment", x = 1:3 + 0.47, xend = 1:3 + 0.53, y = 2, yend = 2,
           arrow = arrow(length = unit(2, "mm"))) +
  scale_fill_manual(values = c("#041E42", "#2C4A7E", "#5B7DB1", "#C8102E")) +
  coord_cartesian(ylim = c(1.2, 2.5)) +
  theme_void()


## -----------------------------------------------------------------------------
#| label: fig-meta-epi
#| fig-cap: "Ratios of odds ratios comparing trials with and without each flaw, from three meta-epidemiological studies. Values below 1 mean that flawed trials report more benefit. Data from @schulz1995 (adapted from its Table 2), @wood2008 and @savovic2012 (95% credible intervals)."
#| fig-width: 8.5
#| fig-height: 5
#| fig-alt: "Forest plot of ratios of odds ratios. Schulz 1995 inadequate concealment 0.59; Wood 2008 subjective outcomes 0.69 and 0.75, objective outcomes 0.91 and 1.01; Savovic 2012 0.89, 0.93 and 0.87."
meta_epi <- tibble::tribble(
  ~source,                ~flaw,                                        ~ror, ~lo,  ~hi,
  "Schulz 1995",          "Allocation concealment inadequate",          0.59, 0.48, 0.73,
  "Schulz 1995",          "Allocation concealment unclear",             0.67, 0.60, 0.75,
  "Schulz 1995",          "Not double blind",                           0.83, 0.71, 0.96,
  "Wood 2008",            "Concealment inadequate/unclear: subjective", 0.69, 0.59, 0.82,
  "Wood 2008",            "Concealment inadequate/unclear: objective",  0.91, 0.80, 1.03,
  "Wood 2008",            "No blinding: subjective outcomes",           0.75, 0.61, 0.93,
  "Wood 2008",            "No blinding: objective outcomes",            1.01, 0.92, 1.10,
  "Savovi\u0107 2012",    "Sequence generation inadequate/unclear",     0.89, 0.82, 0.96,
  "Savovi\u0107 2012",    "Concealment inadequate/unclear",             0.93, 0.87, 0.99,
  "Savovi\u0107 2012",    "Lack of double blinding",                    0.87, 0.79, 0.96) |>
  mutate(source = factor(source, unique(source)),
         flaw = factor(flaw, rev(unique(flaw))))

ggplot(meta_epi, aes(ror, flaw)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#041E42") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)) +
  labs(x = "Ratio of odds ratios (log scale): < 1 = flawed trials show more benefit", y = NULL) +
  theme(strip.text.y = element_text(angle = 0))


## -----------------------------------------------------------------------------
#| label: fig-sponsor
#| fig-cap: "Association between industry sponsorship and favourable results. Data from @lexchin2003 (odds ratio) and @lundh2017 (risk ratios); the two measures are not directly comparable."
#| fig-width: 8
#| fig-height: 2.8
#| fig-alt: "Forest plot: Lexchin 2003 odds ratio 4.05 for outcomes favouring the sponsor; Lundh 2017 risk ratios 1.27 for favourable efficacy results, 1.34 for favourable conclusions and 1.37 with a wide interval for harms."
sponsor <- tibble::tribble(
  ~label,                                                ~est, ~lo,  ~hi,
  "Lexchin 2003: outcome favours sponsor (OR)",          4.05, 2.98, 5.51,
  "Lundh 2017: favourable efficacy results (RR)",        1.27, 1.17, 1.37,
  "Lundh 2017: favourable conclusions (RR)",             1.34, 1.19, 1.51,
  "Lundh 2017: favourable harms results (RR)",           1.37, 0.64, 2.93) |>
  mutate(label = factor(label, rev(label)))
ggplot(sponsor, aes(est, label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#C8102E") +
  geom_text(aes(label = sprintf("%.2f (%.2f to %.2f)", est, lo, hi)), x = log10(7),
            hjust = 0, size = 3.3) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4), limits = c(0.5, 16)) +
  labs(x = "Industry vs other sponsorship (log scale)", y = NULL)


## -----------------------------------------------------------------------------
#| label: fig-early-stop
#| fig-cap: "Simulated trials with a true relative risk of 0.80, analysed after every 50 events and stopped for benefit if Z < -2.96 (two-sided P < 0.003), with a maximum of 1000 events. Each box shows the ratio of the estimated to the true relative risk. The shaded band is the pooled ratio (95% CI) reported by @bassler2010 for truncated trials."
#| fig-width: 8
#| fig-height: 4
#| fig-alt: "Box plots of the ratio of estimated to true relative risk by number of events when the trial stopped: trials stopped with few events overestimate the benefit most; trials run to completion are close to unbiased."
# Normal approximation on the log RR scale: with D events, var(log RR) is about 4 / D
one_sequential_trial <- function(true_rr = 0.8, looks = seq(50, 1000, by = 50),
                                 z_stop = -2.96) {
  info  <- looks / 4                                   # statistical information
  inc   <- rnorm(length(looks), log(true_rr) * diff(c(0, info)), sqrt(diff(c(0, info))))
  score <- cumsum(inc)
  z     <- score / sqrt(info)
  k     <- which(z < z_stop)[1]                       # first look crossing the boundary
  stopped <- !is.na(k)
  if (!stopped) k <- length(looks)
  c(events = looks[k], stopped = stopped, rr_hat = exp(score[k] / info[k]))
}
set.seed(2010)
seq_trials <- as.data.frame(t(replicate(5000, one_sequential_trial()))) |>
  mutate(ratio = rr_hat / 0.8,
         group = case_when(!stopped ~ "Ran to 1000 events",
                           events < 200 ~ "Stopped, < 200 events",
                           events < 500 ~ "Stopped, 200-499 events",
                           TRUE ~ "Stopped, 500+ events"),
         group = factor(group, c("Stopped, < 200 events", "Stopped, 200-499 events",
                                 "Stopped, 500+ events", "Ran to 1000 events")))

ggplot(seq_trials, aes(group, ratio)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.65, ymax = 0.77,
           fill = "#C8102E", alpha = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.2, fill = "#B9C6DE") +
  scale_y_log10() +
  labs(x = NULL, y = "Estimated / true relative risk (log scale)")


## -----------------------------------------------------------------------------
#| label: early-stop-summary
early_gm <- seq_trials |> group_by(group) |>
  summarise(n = n(), geo_mean_ratio = exp(mean(log(ratio))))
early_gm |> kable(digits = 2, col.names = c("", "Trials", "Geometric mean ratio"))


## -----------------------------------------------------------------------------
#| label: fig-rob2-traffic
#| fig-cap: "RoB 2 judgements for four example trials (robvis example data)."
#| fig-width: 6
#| fig-height: 4
#| out-width: "70%"
rob_traffic_light(data_rob2[1:4, ], tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: fig-rob2-summary
#| fig-cap: "Summary of RoB 2 judgements across the example trials."
#| fig-width: 8
#| fig-height: 3
rob_summary(data_rob2, tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: fig-robins
#| fig-cap: "ROBINS-I judgements for three example non-randomised studies."
#| fig-width: 6.5
#| fig-height: 3.6
#| out-width: "75%"
# The example data are called data_robins in robvis 0.3.0 (CRAN), data_robins_i in later versions
robins_example <- if (exists("data_robins_i")) data_robins_i else data_robins
rob_traffic_light(robins_example[1:3, ], tool = "ROBINS-I")


## -----------------------------------------------------------------------------
#| label: tbl-stride
#| tbl-cap: "STRIDE primary endpoint. Adapted from Table 2 of @bonaca2025stride. Medians (IQR) are observed; the estimated treatment ratio (ETR) is from the trial's model."
stride_tab <- data.frame(
  ` ` = "Maximum walking distance at week 52, ratio to baseline",
  Semaglutide = "n = 338; 1.21 (0.95 to 1.55)",
  Placebo     = "n = 345; 1.08 (0.86 to 1.36)",
  `ETR (95% CI)` = "1.13 (1.06 to 1.21)",
  P = "0.0004", check.names = FALSE)
kable(stride_tab)


## -----------------------------------------------------------------------------
#| label: rtm-sim
one_trial <- function(n = 396, effect = 10) {
  true_bl  <- rnorm(2 * n, mean = 185, sd = 40)       # true walking capacity
  obs_bl   <- true_bl + rnorm(2 * n, sd = 30)         # baseline measured with error
  treat    <- rep(0:1, each = n)
  followup <- true_bl + effect * treat + rnorm(2 * n, sd = 30)
  change   <- followup - obs_bl
  c(followup  = unname(diff(tapply(followup, treat, mean))),
    change    = unname(diff(tapply(change, treat, mean))),
    ancova    = unname(coef(lm(followup ~ treat + obs_bl))["treat"]),
    imbalance = unname(diff(tapply(obs_bl, treat, mean))))   # treated - control baseline
}
set.seed(2027)
sims <- as.data.frame(t(replicate(2000, one_trial())))

sims |>
  summarise(across(c(followup, change, ancova), list(mean = mean, sd = sd))) |>
  round(1)


## -----------------------------------------------------------------------------
#| label: fig-rtm
#| fig-cap: "Estimated treatment effect against the chance difference in observed baseline means (treated minus control), in 2000 simulated trials with a true effect of 10 m (dashed line). Follow-up-only estimates rise and change-score estimates fall as the treated arm's baseline rises; ANCOVA estimates do not depend on it."
#| fig-width: 9
#| fig-height: 3.8
sims |>
  pivot_longer(c(followup, change, ancova), names_to = "analysis", values_to = "estimate") |>
  mutate(analysis = factor(analysis, c("followup", "change", "ancova"),
                           c("Follow-up only", "Change from baseline", "ANCOVA"))) |>
  ggplot(aes(imbalance, estimate)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "#041E42") +
  geom_smooth(method = "lm", se = FALSE, colour = "#C8102E") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  facet_wrap(~ analysis) +
  labs(x = "Chance baseline imbalance (m)", y = "Estimated effect (m)")


## -----------------------------------------------------------------------------
#| label: rtm-slope
slopes <- sapply(c(followup = "followup", change = "change", ancova = "ancova"),
                 function(a) unname(coef(lm(sims[[a]] ~ sims$imbalance))[2]))
reliability <- 40^2 / (40^2 + 30^2)   # var(true baseline) / var(observed baseline)
round(c(slopes, reliability = reliability), 2)


## -----------------------------------------------------------------------------
#| label: stride-identity
set.seed(2025)
n <- 396
trial_rtb <- tibble::tibble(
  treat    = rep(0:1, each = n),
  true_bl  = rlnorm(2 * n, log(185), 0.4),
  baseline = true_bl * exp(rnorm(2 * n, sd = 0.2))) |>
  mutate(week52 = true_bl * exp(log(1.13) * treat + rnorm(2 * n, sd = 0.2)))

fit_ratio <- lm(log(week52 / baseline) ~ treat + log(baseline), data = trial_rtb)
fit_level <- lm(log(week52)            ~ treat + log(baseline), data = trial_rtb)
rbind(ratio_outcome = coef(fit_ratio), level_outcome = coef(fit_level)) |> round(4)


## -----------------------------------------------------------------------------
#| label: missing-sim
set.seed(2028)
n <- 396
trial <- tibble(
  treat    = rep(0:1, each = n),
  baseline = rnorm(2 * n, 185, 40) + rnorm(2 * n, sd = 30)) |>
  mutate(followup = baseline + 10 * treat + rnorm(2 * n, sd = 35))
fit_full <- coef(lm(followup ~ treat + baseline, data = trial))["treat"]

# Extreme scenario: the worst 15% of treated and the best 15% of controls go missing
trial_obs <- trial |>
  mutate(gain = followup - baseline) |>
  group_by(treat) |>
  arrange(if_else(treat == 1, -gain, gain), .by_group = TRUE) |>
  slice_head(prop = 0.85) |>
  ungroup()
fit_miss <- coef(lm(followup ~ treat + baseline, data = trial_obs))["treat"]
round(c(complete_data = fit_full, after_selective_loss = fit_miss), 1)


## -----------------------------------------------------------------------------
#| label: fig-stride-rob
#| fig-cap: "An illustrative RoB 2 assessment of STRIDE's primary outcome."
#| fig-width: 5.5
#| fig-height: 2.3
#| out-width: "65%"
stride_rob <- data.frame(
  Study = "STRIDE (2025)",
  D1 = "Low", D2 = "Low", D3 = "High", D4 = "Low", D5 = "Some concerns",
  Overall = "High")
rob_traffic_light(stride_rob, tool = "ROB2")


## -----------------------------------------------------------------------------
#| label: fig-pyramids
#| fig-cap: "Left: a traditional hierarchy of evidence by study design. Right: an alternative pyramid, adapted from a post by Darren Dahly [@dahly2017pyramid]. GRADE sits between the two: design sets the starting point, and the quality of the specific studies decides the rating."
#| fig-width: 10
#| fig-height: 4.2
#| fig-alt: "Two triangles. Left: seven horizontal bands from expert opinion at the base to systematic reviews and meta-analyses at the apex. Right: a small blue apex labelled thoughtful, well-conducted studies of any design, above a large red base labelled the other garbage."
levels_trad <- c("Expert opinion",
                 "Case reports and case series", "Case-control studies",
                 "Cohort studies", "Non-randomised controlled trials",
                 "Randomised controlled trials",
                 "Systematic reviews and meta-analyses")
band <- function(i, n) {               # trapezium i (1 = bottom) of a pyramid of n bands
  y0 <- (i - 1) / n; y1 <- i / n
  data.frame(x = c(-(1 - y0), 1 - y0, 1 - y1, -(1 - y1)) / 2, y = c(y0, y0, y1, y1))
}
trad <- do.call(rbind, lapply(seq_along(levels_trad), function(i)
  cbind(band(i, length(levels_trad)), id = i, panel = "Traditional hierarchy")))
labs_trad <- data.frame(y = (seq_along(levels_trad) - 0.5) / length(levels_trad),
                        label = levels_trad, panel = "Traditional hierarchy") |>
  mutate(x = (1 - y) / 2 + 0.03)      # just right of each band
alt <- rbind(
  data.frame(x = c(-0.5, 0.5, 0.15, -0.15), y = c(0, 0, 0.7, 0.7), id = 1),
  data.frame(x = c(-0.15, 0.15, 0, 0), y = c(0.7, 0.7, 1, 1), id = 2)) |>
  mutate(panel = "An alternative view")
labs_alt <- data.frame(x = 0, y = c(0.3, 0.78),
                       label = c("The other garbage",
                                 "Thoughtful,\nwell-conducted\nstudies of\nany design"),
                       panel = "An alternative view")
pyr <- rbind(trad, alt) |>
  mutate(panel = factor(panel, c("Traditional hierarchy", "An alternative view")),
         fill = case_when(panel == "An alternative view" & id == 1 ~ "#E05A55",
                          panel == "An alternative view" ~ "#87C8F5",
                          TRUE ~ colorRampPalette(c("#C8102E", "#F2C14E", "#5B9B5B",
                                                    "#2C4A7E", "#6A3D9A"))(7)[id]))
ggplot(pyr, aes(x, y)) +
  geom_polygon(aes(group = id, fill = fill), colour = "white") +
  geom_text(data = rbind(labs_trad, labs_alt) |>
              mutate(panel = factor(panel, levels(pyr$panel))),
            aes(label = label, hjust = ifelse(panel == "Traditional hierarchy", 0, 0.5)),
            size = 3.2, lineheight = 0.85, fontface = "bold", colour = "black") +
  scale_fill_identity() +
  facet_wrap(~ panel, scales = "free_x") +
  scale_x_continuous(expand = expansion(add = c(0.02, 0.02))) +
  geom_blank(data = data.frame(x = c(-0.5, 1.3), y = 0, panel = "Traditional hierarchy") |>
               mutate(panel = factor(panel, levels(pyr$panel)))) +
  theme_void() +
  theme(strip.text = element_text(size = 12, face = "bold"))

