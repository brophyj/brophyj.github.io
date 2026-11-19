# R code for chapter: 07-appraising-reviews.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch07
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


## -----------------------------------------------------------------------------
#| label: collider-sim
set.seed(123)
n <- 5000
population <- tibble::tibble(
  income = rnorm(n),                       # standardised; independent of BP
  bp     = rnorm(n)) |>
  mutate(visit = rbinom(n, 1, plogis(-2 + 2 * income + 2 * bp)))   # the collider
sampled <- filter(population, visit == 1)

rbind(
  population = coef(summary(lm(bp ~ income, data = population)))["income", c(1, 4)],
  visited    = coef(summary(lm(bp ~ income, data = sampled)))["income", c(1, 4)]) |>
  round(3)


## -----------------------------------------------------------------------------
#| label: fig-collider
#| fig-cap: "Income and blood pressure are independent in the population (grey line), but negatively associated among the people who had a medical visit (red line), because a visit is a common effect of both."
#| fig-width: 7
#| fig-height: 4.5
ggplot(population, aes(income, bp)) +
  geom_point(aes(colour = factor(visit)), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30") +
  geom_smooth(data = sampled, method = "lm", se = FALSE, colour = "#C8102E") +
  scale_colour_manual(values = c("0" = "grey75", "1" = "#C8102E"),
                      labels = c("No visit", "Visit"), name = NULL) +
  labs(x = "Income (standardised)", y = "Blood pressure (standardised)")


## -----------------------------------------------------------------------------
#| label: fig-myopia-dag
#| fig-cap: "Two causal models for fluoroquinolone (FQ) use, myopia and retinal detachment (RD). (A) Myopia as a confounder. (B) Myopia as a collider of low socioeconomic status (SES) and unknown causes of detachment (U): an M-structure. The FQ to RD arrow is the effect under study."
#| fig-width: 9
#| fig-height: 3.4
nodes <- tibble::tribble(
  ~panel, ~node,    ~x,  ~y,
  "A",    "FQ",     0,   1,
  "A",    "RD",     2,   1,
  "A",    "Myopia", 1,   0,
  "B",    "SES",    0,   2,
  "B",    "U",      2,   2,
  "B",    "Myopia", 1,   1.2,
  "B",    "FQ",     0,   0,
  "B",    "RD",     2,   0)
edges <- tibble::tribble(
  ~panel, ~from,    ~to,
  "A",    "FQ",     "RD",
  "A",    "Myopia", "FQ",
  "A",    "Myopia", "RD",
  "B",    "SES",    "FQ",
  "B",    "SES",    "Myopia",
  "B",    "U",      "Myopia",
  "B",    "U",      "RD",
  "B",    "FQ",     "RD") |>
  left_join(nodes, by = c("panel", from = "node")) |>
  left_join(nodes, by = c("panel", to = "node"), suffix = c("", "end")) |>
  mutate(shrink = 0.12,
         dx = xend - x, dy = yend - y, len = sqrt(dx^2 + dy^2),
         x = x + shrink * dx / len, y = y + shrink * dy / len,
         xend = xend - shrink * dx / len, yend = yend - shrink * dy / len)
ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
               colour = "grey35") +
  geom_label(data = nodes, aes(x, y, label = node,
                               fill = node == "Myopia"),
             size = 3.6, label.size = 0, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#F2C14E", "FALSE" = "#DCE3EE")) +
  facet_wrap(~ panel, labeller = as_labeller(c(A = "(A) Myopia as confounder",
                                               B = "(B) Myopia as collider"))) +
  coord_cartesian(xlim = c(-0.3, 2.3), ylim = c(-0.3, 2.3)) +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 11))


## -----------------------------------------------------------------------------
#| label: mbias-sim
one_mbias <- function(n = 20000, strength = 1.5) {
  ses    <- rnorm(n)                                    # low SES (standardised)
  u      <- rnorm(n)                                    # unknown causes of RD
  myopia <- rbinom(n, 1, plogis(-1 + strength * ses + strength * u))
  fq     <- rbinom(n, 1, plogis(-1.5 + strength * ses))
  rd     <- rbinom(n, 1, plogis(-3 + strength * u))     # no effect of fq
  c(crude    = exp(coef(glm(rd ~ fq, family = binomial))[["fq"]]),
    adjusted = exp(coef(glm(rd ~ fq + myopia, family = binomial))[["fq"]]))
}
set.seed(2012)
mbias <- rbind(
  data.frame(arrows = "Strong (log OR 1.5)", t(replicate(200, one_mbias(strength = 1.5)))),
  data.frame(arrows = "Modest (log OR 0.5)", t(replicate(200, one_mbias(strength = 0.5)))))
mbias_summary <- mbias |>
  group_by(arrows) |>
  summarise(across(c(crude, adjusted), median))
kable(mbias_summary, digits = 3,
      col.names = c("Arrows into myopia", "Crude OR (median)", "Adjusted for myopia (median)"))


## -----------------------------------------------------------------------------
#| label: amstar-function
critical_items <- c(2, 4, 7, 9, 11, 13, 15)

amstar2_rating <- function(ratings) {
  # ratings: character vector of 16 responses ("Yes", "Partial Yes", "No", "No MA")
  flaw <- ratings == "No"                 # "Partial Yes" and "No MA" are not flaws
  n_critical     <- sum(flaw[critical_items])
  n_noncritical  <- sum(flaw[-critical_items])
  if (n_critical > 1) "Critically low"
  else if (n_critical == 1) "Low"
  else if (n_noncritical > 1) "Moderate"
  else "High"
}


## -----------------------------------------------------------------------------
#| label: amstar-profiles
review_a <- rep("Yes", 16); review_a[c(2, 4)] <- "No"       # no protocol; narrow search
review_b <- rep("Yes", 16); review_b[c(11, 13)] <- "No"     # wrong model; RoB ignored
c(A = amstar2_rating(review_a), B = amstar2_rating(review_b))


## -----------------------------------------------------------------------------
#| label: amstar-data
amstar_reviews <- readr::read_csv("data/ch07/amstar2vis_reviews.csv",
                                  show_col_types = FALSE)
item_cols <- paste0("Item", 1:16)
amstar_reviews$overall <- apply(amstar_reviews[, item_cols], 1, amstar2_rating)
overall_tab <- table(factor(amstar_reviews$overall,
                            levels = c("High", "Moderate", "Low", "Critically low")))
overall_tab
pct_no <- round(100 * colMeans(amstar_reviews[, item_cols] == "No"))   # % rated "No", by item


## -----------------------------------------------------------------------------
#| label: fig-amstar-items
#| fig-cap: "AMSTAR 2 item ratings across 495 systematic reviews; critical items are marked with an asterisk. Data distributed with the amstar2Vis package [@bougioukas2024amstar2vis]."
#| fig-width: 9
#| fig-height: 5.5
amstar_reviews |>
  pivot_longer(all_of(item_cols), names_to = "item", values_to = "rating") |>
  mutate(item_no = as.integer(sub("Item", "", item)),
         item = paste0(item_no, ifelse(item_no %in% critical_items, "*", "")),
         item = factor(item, levels = rev(unique(item[order(item_no)]))),
         rating = factor(rating, levels = c("No", "Partial Yes", "Yes", "No MA"))) |>
  ggplot(aes(y = item, fill = rating)) +
  geom_bar(position = "fill") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#C8102E", "Partial Yes" = "#F2C14E",
                               "Yes" = "#2E7D32", "No MA" = "grey70"), name = NULL) +
  labs(x = "Percentage of reviews", y = "AMSTAR 2 item") +
  theme(legend.position = "bottom")

