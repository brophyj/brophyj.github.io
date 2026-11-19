# R code for chapter: 01-introduction.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch01
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


## -----------------------------------------------------------------------------
#| label: fig-konstam
#| fig-cap: "Main results of the 2001 pooled analysis of rofecoxib trials: relative risk of the combined end point of cardiovascular, haemorrhagic and unknown deaths, non-fatal myocardial infarction and non-fatal stroke. Squares are sized by patient-years. Data from @konstam2001."
#| fig-alt: "Forest plot of three comparisons: rofecoxib vs placebo 0.84, vs non-naproxen NSAIDs 0.79, vs naproxen 1.69."
#| fig-width: 8.5
#| fig-height: 3.2
#| out-width: "90%"
#| code-fold: true
#| code-summary: "Show the code"
konstam <- tibble::tribble(
  ~comparison,                        ~rr,  ~lo,  ~hi,  ~pyears,
  "Rofecoxib vs placebo",              0.84, 0.51, 1.38, 3867,
  "Rofecoxib vs non-naproxen NSAIDs",  0.79, 0.40, 1.55, 2918,
  "Rofecoxib vs naproxen",             1.69, 1.07, 2.69, 8364) |>
  mutate(comparison = factor(comparison, levels = rev(comparison)),
         label = sprintf("%.2f (%.2f to %.2f)", rr, lo, hi))
ggplot(konstam, aes(x = rr, y = comparison)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.8) +
  geom_point(aes(size = pyears), shape = 15, colour = "#041E42") +
  geom_text(aes(x = 4.2, label = label), hjust = 0, size = 3.8) +
  annotate("text", x = c(0.93, 1.07), y = 3.55, hjust = c(1, 0),
           label = c("\u2190 Lower risk with rofecoxib", "Higher risk with rofecoxib \u2192"),
           size = 3.5, colour = "grey30") +
  annotate("text", x = 4.2, y = 3.55, hjust = 0, label = "RR (95% CI)", size = 3.8, fontface = "bold") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4"), limits = c(0.25, 9)) +
  scale_size_area(max_size = 7, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Relative risk (log scale)", y = NULL) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())


## -----------------------------------------------------------------------------
#| label: tbl-ross
#| tbl-cap: "Placebo-controlled rofecoxib trials of 4 weeks or longer conducted by the manufacturer, as listed by @ross2009 (Table 1). Trials 098/103 and 120/121 were reported together, so the 28 rows cover 30 trials. LPO: last patient out."
#| code-fold: true
#| code-summary: "Show the data"
rofecoxib_trials <- tribble(
  ~source,                                ~trial,       ~indication,                      ~dose_mg,          ~weeks, ~lpo,
  "Ehrich et al, 1999",                   "010",        "Osteoarthritis",                 "25 and 125",       6,     "1996-02-08",
  "Ehrich et al, 2001",                   "029",        "Osteoarthritis",                 "12.5, 25 and 50",  6,     "1997-02-05",
  "Saag et al, 2000",                     "033",        "Osteoarthritis",                 "12.5 and 25",      6,     "1997-11-18",
  "Day et al, 2000",                      "040",        "Osteoarthritis",                 "12.5 and 25",      6,     "1998-01-01",
  "Laine et al, 1999",                    "044",        "Osteoarthritis",                 "25 and 50",       24,     "1998-02-18",
  "Hawkey et al, 2000",                   "045",        "Osteoarthritis",                 "25 and 50",       24,     "1998-02-18",
  "Truitt et al, 2001",                   "058",        "Osteoarthritis",                 "12.5 and 25",      6,     "1998-04-01",
  "Unpublished",                          "083",        "Osteoarthritis",                 "25",              64,     "2000-02-09",
  "Kivitz et al, 2004",                   "085",        "Osteoarthritis",                 "12.5",             6,     "1999-03-03",
  "Weaver et al, 2006",                   "090",        "Osteoarthritis",                 "12.5",             6,     "1999-05-17",
  "Smugar et al, 2006",                   "112",        "Osteoarthritis",                 "12.5 and 25",      6,     "2000-09-08",
  "Smugar et al, 2006",                   "116",        "Osteoarthritis",                 "25",               6,     "2000-06-22",
  "Laine et al, 2004",                    "136",        "Osteoarthritis",                 "25",              12,     "2002-02-05",
  "Birbara et al, 2006",                  "219",        "Osteoarthritis",                 "12.5",             6,     "2003-11-28",
  "Birbara et al, 2006",                  "220",        "Osteoarthritis",                 "12.5",             6,     "2003-11-24",
  "Unpublished",                          "017",        "Rheumatoid arthritis",           "125 and 175",      6,     "1997-05-21",
  "Schnitzer et al, 1999",                "068",        "Rheumatoid arthritis",           "25 and 50",        8,     "1998-09-10",
  "Truitt et al, 2001 (abstract only)",   "096",        "Rheumatoid arthritis",           "12.5 and 25",     12,     "2000-07-21",
  "Geusens et al, 2002",                  "097",        "Rheumatoid arthritis",           "25 and 50",       12,     "2000-06-06",
  "Hawkey et al, 2003",                   "098 and 103","Rheumatoid arthritis",           "50",              12,     "2000-07-06",
  "Thal et al, 2005",                     "078",        "Alzheimer disease",              "25",             208,     "2003-04-23",
  "Reines et al, 2004",                   "091",        "Alzheimer disease",              "25",              52,     "2000-11-30",
  "Unpublished",                          "126",        "Alzheimer disease",              "25",              52,     "2001-05-30",
  "Nickel et al, 2003",                   "118",        "Chronic nonbacterial prostatitis","25 and 50",       6,     "2000-07-26",
  "Katz et al, 2003",                     "120 and 121","Low back pain",                  "25 and 50",        4,     "2000-06-27",
  "Bresalier et al, 2005; Baron et al, 2008", "122",    "Colorectal adenomas",            "25",             156,     "2004-09-30",
  "Unpublished",                          "125",        "Migraine prophylaxis",           "25",              12,     "2001-06-29",
  "Unpublished",                          "129",        "Familial adenomatous polyposis", "25",              24,     "2002-05-14"
) |>
  mutate(lpo = as.Date(lpo))

rofecoxib_trials |>
  arrange(lpo) |>
  mutate(lpo = format(lpo, "%d %b %Y")) |>
  kable(col.names = c("Source", "Trial no.", "Indication", "Rofecoxib dose (mg)",
                      "Planned duration (wk)", "LPO"))


## -----------------------------------------------------------------------------
#| label: fig-ross
#| fig-cap: "Cumulative pooled analysis of cardiovascular thrombotic events or death in the 30 placebo-controlled trials of rofecoxib, by the date each trial's data became available. Filled points and lines: rate ratio and 95% CI published by @ross2009 (red: lower limit above 1; arrows: upper limit beyond 10). Open circles: our reconstruction from the summary data."
#| fig-alt: "Cumulative forest plot from February 1998 to September 2004: estimates above 1 throughout, with intervals narrowing and excluding 1 from 2002."
#| fig-width: 8.5
#| fig-height: 7.5
#| out-width: "90%"
#| code-fold: true
#| code-summary: "Show the code"
library(metafor)
ross_cum <- tibble::tribble(
  ~trial, ~available, ~ev_r, ~pyr_r, ~ev_p, ~pyr_p, ~rr,  ~lo,  ~hi,
  "010", "1996-05",  2,   16,  0,    7,   NA,   NA,    NA,
  "029", "1997-05",  3,   46,  0,   16,   NA,   NA,    NA,
  "017", "1997-08",  1,    8,  0,    7,   NA,   NA,    NA,
  "033", "1998-02",  1,   66,  1,    9, 2.58, 0.38,  60.8,
  "040", "1998-04",  2,   72,  0,   11, 2.92, 0.45,  67.0,
  "044", "1998-05",  3,  154,  0,   52, 3.96, 0.66,  87.2,
  "045", "1998-05",  3,  157,  3,   61, 1.18, 0.41,  4.13,
  "058", "1998-06",  1,   21,  0,    6, 1.35, 0.47,  4.77,
  "068", "1998-12",  1,   49,  0,   24, 1.47, 0.51,  5.15,
  "085", "1999-06",  1,   61,  0,   28, 1.58, 0.56,  5.50,
  "090", "1999-08",  5,   56,  0,   27, 2.16, 0.79,  7.37,
  "083", "2000-05",  0,   21,  1,   21, 1.80, 0.71,  5.40,
  "097", "2000-09",  0,  137,  0,   62, 1.80, 0.71,  5.40,
  "121", "2000-09",  0,   23,  0,   11, 1.80, 0.71,  5.40,
  "116", "2000-09",  1,   54,  0,   15, 1.86, 0.74,  5.56,
  "120", "2000-09",  1,   28,  0,   14, 1.91, 0.77,  5.62,
  "098", "2000-10",  0,   11,  1,   12, 1.70, 0.71,  4.61,
  "103", "2000-10",  0,   44,  0,   45, 1.70, 0.71,  4.61,
  "096", "2000-10",  5,   97,  0,   58, 2.18, 0.93,  5.81,
  "118", "2000-10",  0,   15,  0,    8, 2.18, 0.93,  5.81,
  "112", "2000-12",  0,  104,  0,   15, 2.18, 0.93,  5.81,
  "091", "2001-02", 23,  273, 20,  189, 1.43, 0.88,  2.31,
  "078", "2001-06", 57, 1257, 46, 1478, 1.35, 1.00,  1.96,
  "126", "2001-08", 12,  193, 11,  100, 1.32, 1.00,  1.75,
  "125", "2001-09",  0,   23,  0,   22, 1.32, 1.00,  1.75,
  "078", "2002-04", 78, 1148, 59, 1374, 1.39, 1.07,  1.80,
  "136", "2002-05",  1,   95,  1,  201, 1.39, 1.07,  1.81,
  "129", "2002-08",  0,    3,  0,    4, 1.39, 1.07,  1.81,
  "078", "2003-07", 91, 1047, 70, 1075, 1.38, 1.08,  1.77,
  "220", "2004-02",  0,   18,  0,    8, 1.38, 1.08,  1.77,
  "219", "2004-02",  0,   18,  1,    8, 1.37, 1.07,  1.75,
  "122", "2004-09", 64, 3429, 42, 3598, 1.43, 1.16,  1.76) |>
  mutate(available = as.Date(paste0(available, "-01")), step = row_number())

# Our reconstruction: at each step, a Mantel-Haenszel rate ratio over the trials
# available so far, using the most recent data for each trial (078 was updated twice)
mh_step <- function(k) {
  d <- ross_cum[1:k, ] |> group_by(trial) |> slice_tail(n = 1) |> ungroup()
  if (sum(d$ev_p) == 0) return(NA_real_)
  fit <- rma.mh(x1i = ev_r, t1i = pyr_r, x2i = ev_p, t2i = pyr_p, data = d, measure = "IRR")
  exp(fit$beta[1])
}
ross_cum$rr_mh <- sapply(ross_cum$step, mh_step)

plot_df <- ross_cum |>
  filter(!is.na(rr)) |>
  mutate(label = paste0(format(available, "%b %Y"), " - trial ", trial),
         label = factor(label, levels = rev(unique(label))),
         hi_plot = pmin(hi, 10), clipped = hi > 10,
         sig = lo > 1)
ggplot(plot_df, aes(y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = lo, xend = hi_plot, yend = label, colour = sig), linewidth = 0.7,
               arrow = NULL) +
  geom_segment(data = filter(plot_df, clipped), aes(x = 8, xend = 10, yend = label),
               arrow = arrow(length = unit(0.12, "cm")), colour = "grey40", linewidth = 0.7) +
  geom_point(aes(x = rr, colour = sig, shape = "Ross et al. (subject-level Cox models)"), size = 2.4) +
  geom_point(aes(x = rr_mh, shape = "Reconstructed from summary data (Mantel-Haenszel)"),
             size = 2.2, colour = "grey45") +
  scale_shape_manual(values = c(1, 16), name = NULL) +
  scale_colour_manual(values = c(`FALSE` = "#041E42", `TRUE` = "#C8102E"), guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10), labels = c("0.5", "1", "2", "5", "10")) +
  coord_cartesian(xlim = c(0.35, 10)) +
  labs(x = "Cumulative rate ratio, rofecoxib vs placebo (log scale)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())


## -----------------------------------------------------------------------------
#| label: fig-steps
#| fig-cap: "The 24 steps of a systematic review and meta-analysis, grouped into six phases. Adapted from the 24-step guide of @muka2020steps."
#| fig-alt: "Six columns of numbered boxes: Plan (steps 1-6), Search (7-9), Select (10-15), Extract and appraise (16-17), Analyse (18-22) and Report (23-24)."
#| fig-width: 12
#| fig-height: 6.2
#| out-width: "100%"
#| code-fold: true
#| code-summary: "Show the code"
steps <- tibble::tibble(
  step = 1:24,
  text = c("Define the research question", "Establish the team",
           "Define the search strategy", "Define the selection criteria",
           "Design the data collection form", "Write the protocol and register the review",
           "Run the search in multiple databases", "Collect all references in a single file",
           "Remove duplicates", "Two reviewers screen titles and abstracts",
           "Compare decisions and select for retrieval", "Retrieve full texts and apply the criteria",
           "Contact experts", "Search for additional references",
           "Final selection and flow diagram", "Extract data in pairs",
           "Assess quality and risk of bias", "Prepare the database for analysis",
           "Descriptive synthesis", "Meta-analyse or not?",
           "Explore heterogeneity", "Check for reporting bias",
           "Rate the certainty of the evidence", "Update, report and submit"),
  phase = rep(c("Plan", "Search", "Select", "Extract & appraise", "Analyse", "Report"),
              times = c(6, 3, 6, 2, 5, 2))) |>
  mutate(phase = factor(phase, levels = unique(phase)),
         col = as.integer(phase)) |>
  group_by(phase) |> mutate(row = row_number()) |> ungroup()

phase_cols <- c("#041E42", "#2E6DA4", "#3B8C88", "#8C6D1F", "#C8102E", "#6B4C9A")
ggplot(steps) +
  # phase headers
  geom_label(data = distinct(steps, phase, col),
             aes(x = col, y = 0, label = phase, fill = phase),
             colour = "white", fontface = "bold", size = 4, label.size = 0,
             label.padding = unit(0.35, "lines")) +
  # step boxes
  geom_tile(aes(x = col, y = row, colour = phase), fill = "white",
            width = 0.92, height = 0.82, linewidth = 0.6) +
  geom_point(aes(x = col - 0.36, y = row, colour = phase), size = 7.5) +
  geom_text(aes(x = col - 0.36, y = row, label = step), colour = "white",
            fontface = "bold", size = 3.4) +
  geom_text(aes(x = col - 0.27, y = row, label = stringr::str_wrap(text, 20)),
            hjust = 0, size = 3, lineheight = 0.9) +
  scale_y_reverse() +
  scale_fill_manual(values = phase_cols, guide = "none") +
  scale_colour_manual(values = phase_cols, guide = "none") +
  coord_cartesian(xlim = c(0.52, 6.48), clip = "off") +
  labs(subtitle = "Before step 1: has the question already been answered? Search PROSPERO and the main databases for existing reviews.") +
  theme_void() +
  theme(plot.subtitle = element_text(size = 11, colour = "grey30", margin = margin(b = 8)),
        plot.margin = margin(5, 10, 5, 10), plot.title.position = "plot")

