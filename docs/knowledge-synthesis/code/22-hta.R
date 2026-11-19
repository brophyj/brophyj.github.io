# R code for chapter: 22-hta.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch22
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
suppressPackageStartupMessages({
  library(metafor)
  library(purrr)
  library(mgcv)      # generalised additive models, for EVPPI
})
navy <- "#041E42"; red <- "#B22222"; green <- "#2E7D32"; orange <- "#E69F00"


## -----------------------------------------------------------------------------
#| label: fig-mcgregor
#| fig-cap: "Extra cost per coronary bypass operation by annual volume, reproducing the cost model of McGregor and Pelletier (1978; 1976 Canadian dollars, physicians' fees excluded)."
#| fig-width: 7
#| fig-height: 3.6
ops <- 5:1000
cost_per_op <- 18993 * ceiling(ops / 400) / ops + 733
ggplot(tibble::tibble(ops, cost_per_op), aes(ops, cost_per_op)) +
  geom_line(linewidth = 1, colour = navy) +
  geom_vline(xintercept = c(50, 200), linetype = 2, colour = red) +
  annotate("text", x = c(55, 205), y = 3800, hjust = 0, size = 3.3, colour = red,
           label = c("~50: cost levels off", "~200: suggested\nfor good outcomes")) +
  labs(x = "Operations per year", y = "Extra cost per operation ($)")


## -----------------------------------------------------------------------------
#| label: tbl-absolute
#| tbl-cap: "The same relative risk, 0.75, at four baseline risks."
tibble::tibble(`Baseline risk per year` = c(0.08, 0.15, 0.25, 0.35)) |>
  mutate(`Absolute risk reduction` = `Baseline risk per year` * 0.25,
         NNT = round(1 / `Absolute risk reduction`),
         across(1:2, \(x) sprintf("%.2f%%", 100 * x))) |>
  kable(align = "ccc")


## -----------------------------------------------------------------------------
#| label: bucher
trials <- tibble::tribble(
  ~trial,            ~treat,          ~hr,  ~lo,  ~hi,
  "EMPHASIS-HF",     "Eplerenone",    0.63, 0.54, 0.74,
  "PARADIGM-HF",     "Sac/val",       0.80, 0.73, 0.87,
  "DAPA-HF",         "Dapagliflozin", 0.74, 0.65, 0.85,
  "EMPEROR-Reduced", "Empagliflozin", 0.75, 0.65, 0.86) |>
  mutate(loghr = log(hr), se = (log(hi) - log(lo)) / (2 * qnorm(0.975)))
nma <- rma(yi = loghr, sei = se, mods = ~ 0 + treat, data = trials, method = "EE")
names(coef(nma))                                    # Dapa, Empa, Eplerenone, Sac/val
bd <- predict(nma, newmods = c(1, 0, 0, -1), transf = exp)      # dapagliflozin vs sac/val
round(c(HR = bd$pred, lower = bd$ci.lb, upper = bd$ci.ub), 2)


## -----------------------------------------------------------------------------
#| label: fig-combination
#| fig-cap: "Each drug against its own trial's control (navy), and the combination of sacubitril/valsartan, a mineralocorticoid antagonist and an SGLT2 inhibitor estimated by multiplying hazard ratios (red). Outcome: cardiovascular death or HF hospitalisation."
#| fig-width: 8
#| fig-height: 3.4
combo_l  <- log(0.80) + log(0.63) + log(0.75)
combo_se <- sqrt(trials$se[2]^2 + trials$se[1]^2 + ((log(0.85) - log(0.65)) / 3.92)^2)
bind_rows(mutate(trials, type = "Each drug vs its own control"),
          tibble::tibble(treat = "ARNI + MRA + SGLT2i vs ACEi", hr = exp(combo_l),
                         lo = exp(combo_l - 1.96 * combo_se), hi = exp(combo_l + 1.96 * combo_se),
                         type = "Combination")) |>
  mutate(treat = forcats::fct_rev(forcats::fct_inorder(treat))) |>
  ggplot(aes(hr, treat, colour = type)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi)) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1)) +
  scale_colour_manual(values = c(Combination = red, `Each drug vs its own control` = navy), name = NULL) +
  labs(x = "Hazard ratio (log scale)", y = NULL) + theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: markov
p2r <- function(p) -log(1 - p); r2p <- function(r) 1 - exp(-r)
hf_model <- function(hr_hosp = 1, hr_cvd = 1, p_hosp = .10, p_cvd = .06, p_cvd_post = .15,
                     p_oth = .03, u_stable = .78, u_post = .65, dis_hosp = .10,
                     c_stable = 1500, c_post = 2500, c_hosp = 4000, c_drug = 0,
                     cycles = 30, disc = .035) {
  ph  <- r2p(p2r(p_hosp) * hr_hosp)              # HRs act on rates
  pc  <- r2p(p2r(p_cvd) * hr_cvd)
  pcp <- r2p(p2r(p_cvd_post) * hr_cvd)
  P <- matrix(c(1 - ph - pc - p_oth, ph,              pc + p_oth,     # from Stable
                0,                   1 - pcp - p_oth, pcp + p_oth,    # from Post-hosp
                0,                   0,               1), 3, byrow = TRUE)
  s <- c(1, 0, 0); cost <- qaly <- 0
  for (t in 1:cycles) {
    hosp  <- s[1] * ph                           # hospitalised this cycle
    s_new <- as.vector(s %*% P)
    d     <- 1 / (1 + disc)^(t - 1)              # discounting
    alive <- (s[1:2] + s_new[1:2]) / 2           # half-cycle correction
    qaly  <- qaly + d * (sum(alive * c(u_stable, u_post)) - hosp * dis_hosp)
    cost  <- cost + d * (sum(alive * (c(c_stable, c_post) + c_drug)) + hosp * c_hosp)
    s <- s_new
  }
  c(cost = cost, qaly = qaly)
}
sc <- hf_model(); tx <- hf_model(hr_hosp = .70, hr_cvd = .82, c_drug = 477)
rbind(`Standard care` = sc, Dapagliflozin = tx, Difference = tx - sc) |> round(2)
icer_det <- unname((tx - sc)[1] / (tx - sc)[2])


## -----------------------------------------------------------------------------
#| label: psa
set.seed(2026); n_psa <- 5000
beta_ms  <- function(n, m, s) { a <- m * (m * (1 - m) / s^2 - 1); rbeta(n, a, a * (1 - m) / m) }
gamma_ms <- function(n, m, s) rgamma(n, shape = (m / s)^2, rate = m / s^2)
pars <- tibble::tibble(
  lhr_hosp = rnorm(n_psa, log(.70), (log(.83) - log(.59)) / 3.92),   # DAPA-HF, worsening HF
  lhr_cvd  = rnorm(n_psa, log(.82), (log(.98) - log(.69)) / 3.92),   # DAPA-HF, CV death
  p_hosp = rbeta(n_psa, 50, 450), p_cvd = rbeta(n_psa, 30, 470), p_cvd_post = rbeta(n_psa, 30, 170),
  u_stable = beta_ms(n_psa, .78, .03), u_post = beta_ms(n_psa, .65, .05),
  c_hosp = gamma_ms(n_psa, 4000, 800), c_post = gamma_ms(n_psa, 2500, 500),
  c_stable = gamma_ms(n_psa, 1500, 300))
psa <- pmap_dfr(pars, function(lhr_hosp, lhr_cvd, ...) {
  a <- hf_model(...)
  b <- hf_model(hr_hosp = exp(lhr_hosp), hr_cvd = exp(lhr_cvd), c_drug = 477, ...)
  tibble::tibble(dC = b[["cost"]] - a[["cost"]], dQ = b[["qaly"]] - a[["qaly"]])
})
icer_psa <- mean(psa$dC) / mean(psa$dQ)
psa |> summarise(across(c(dC, dQ), list(mean = mean, lo = \(x) quantile(x, .025),
                                         hi = \(x) quantile(x, .975)))) |> round(2)


## -----------------------------------------------------------------------------
#| label: fig-ceac
#| fig-cap: "Left: 5000 PSA draws on the cost-effectiveness plane, with the £20,000/QALY line and the mean (diamond). Right: the cost-effectiveness acceptability curve, the probability that dapagliflozin is cost-effective at each threshold."
#| fig-width: 10
#| fig-height: 4
wtp  <- seq(0, 40000, 500)
ceac <- tibble::tibble(wtp, p = map_dbl(wtp, \(w) mean(w * psa$dQ - psa$dC > 0)))
p20 <- ceac$p[ceac$wtp == 20000]; p5 <- ceac$p[ceac$wtp == 5000]
g1 <- ggplot(psa, aes(dQ, dC)) +
  geom_point(alpha = .15, size = .6, colour = navy) +
  geom_hline(yintercept = 0, colour = "grey50") + geom_vline(xintercept = 0, colour = "grey50") +
  geom_abline(slope = 20000, intercept = 0, colour = red) +
  annotate("point", x = mean(psa$dQ), y = mean(psa$dC), colour = red, size = 3, shape = 18) +
  labs(x = "Incremental QALYs", y = "Incremental cost (£)")
g2 <- ggplot(ceac, aes(wtp / 1000, p)) +
  geom_line(linewidth = 1.2, colour = navy) +
  geom_vline(xintercept = c(5, 20), linetype = c(3, 2), colour = red) +
  annotate("text", x = c(5, 20) + 0.8, y = c(p5, p20) - 0.06, hjust = 0, size = 3.3,
           label = sprintf("%.0f%% at £%s", 100 * c(p5, p20), c("5,000", "20,000"))) +
  labs(x = "Willingness to pay (£000 per QALY)", y = "P(cost-effective)") + ylim(0, 1)
patchwork::wrap_plots(g1, g2)


## -----------------------------------------------------------------------------
#| label: magnesium
mag <- tibble::tribble(~study, ~year, ~Ee, ~Ne, ~Ec, ~Nc,
  "Morton",1984,1,40,2,36,       "Rasmussen",1986,9,135,23,135, "Smith",1986,2,200,7,200,
  "Abraham",1987,1,48,1,46,      "Feldstedt",1988,10,150,8,148, "Shechter",1989,1,59,9,56,
  "Ceremuzynski",1989,1,25,3,23, "Bertschat",1989,0,22,1,21,    "Singh",1990,6,76,11,75,
  "Pereira",1990,1,27,7,27,      "Shechter",1991,2,89,12,80,    "Golf",1991,5,23,13,33,
  "Thogersen",1991,4,130,8,122,  "LIMIT-2",1992,90,1159,118,1157, "Shechter",1995,4,107,17,108) |>
  mutate(study = paste(study, year))
mag_es <- escalc(measure = "OR", ai = Ee, n1i = Ne, ci = Ec, n2i = Nc, data = mag) |>
  mutate(sei = sqrt(vi))
re_mag <- rma(yi, vi, data = mag_es, method = "DL")
lim <- mag$study == "LIMIT-2 1992"
re_small <- rma(yi, vi, data = mag_es[!lim, ], method = "DL")
or_limit <- exp(mag_es$yi[lim] + c(0, -1.96, 1.96) * mag_es$sei[lim])
pi_mag <- predict(re_mag, transf = exp)
round(c(OR = pi_mag$pred, lower = pi_mag$ci.lb, upper = pi_mag$ci.ub,
        PI_lower = pi_mag$pi.lb, PI_upper = pi_mag$pi.ub, I2 = re_mag$I2), 2)


## -----------------------------------------------------------------------------
#| label: magnesium-brms
#| eval: false
# library(brms)
# fit_mag <- brm(yi | se(sei) ~ 1 + (1 | study), data = mag_es,
#                prior = c(prior(normal(0, 1), class = Intercept), prior(normal(0, 0.5), class = sd)),
#                control = list(adapt_delta = 0.95), seed = 629)
# dr <- as_draws_df(fit_mag)[, c("b_Intercept", "sd_study__Intercept")]


## -----------------------------------------------------------------------------
#| label: fig-magnesium
#| fig-cap: "Magnesium after myocardial infarction: posterior for the mean effect and posterior predictive distribution for a new trial (random-effects model of the 15 trials), and the result of ISIS-4 (red band, 95% CI)."
#| fig-width: 7.5
#| fig-height: 3.8
dr <- read.csv("data/ch22/magnesium_draws.csv")          # 4000 posterior draws from fit_mag
set.seed(629)
theta_new <- rnorm(nrow(dr), dr$b_Intercept, dr$sd_study)  # predictive: a new trial
isis <- escalc(measure = "OR", ai = 2216, n1i = 29011, ci = 2103, n2i = 29039)
isis_or <- exp(isis$yi + c(0, -1, 1) * 1.96 * sqrt(isis$vi))
tibble::tibble(`Mean effect (posterior)` = exp(dr$b_Intercept), `New trial (predictive)` = exp(theta_new)) |>
  pivot_longer(everything()) |>
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = .5, colour = NA) +
  annotate("rect", xmin = isis_or[2], xmax = isis_or[3], ymin = 0, ymax = Inf, fill = red, alpha = .35) +
  annotate("text", x = isis_or[1], y = 3.5, label = "ISIS-4", colour = red, hjust = -0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10(limits = c(0.1, 2.5), breaks = c(0.1, 0.25, 0.5, 1, 2)) +
  scale_fill_manual(values = c(navy, orange), name = NULL) +
  labs(x = "Odds ratio for death (log scale)", y = "Density") + theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: evpi
evpi_at <- function(lam) {
  inb <- lam * psa$dQ - psa$dC                      # incremental NB of dapagliflozin
  mean(pmax(inb, 0)) - max(mean(inb), 0)            # perfect information - decide now
}
evpi_curve <- tibble::tibble(wtp = seq(0, 30000, 250)) |> mutate(evpi = map_dbl(wtp, evpi_at))
pop <- 75000 * sum(1 / 1.035^(0:9))   # illustrative: 75,000 new patients/yr for 10 years, discounted
ev6 <- evpi_at(6000); ev20 <- evpi_at(20000)
round(c(EVPI_6000 = ev6, EVPI_20000 = ev20, patients = pop))


## -----------------------------------------------------------------------------
#| label: fig-evpi
#| fig-cap: "Per-patient EVPI for the dapagliflozin model at each willingness to pay. EVPI peaks near the ICER (dashed), where the two options are closest."
#| fig-width: 7
#| fig-height: 3.6
ggplot(evpi_curve, aes(wtp / 1000, evpi)) +
  geom_line(linewidth = 1.2, colour = navy) +
  geom_vline(xintercept = icer_psa / 1000, linetype = 2, colour = red) +
  annotate("point", x = c(6, 20), y = c(ev6, ev20), size = 3, colour = red) +
  annotate("text", x = c(6, 20) + 0.5, y = c(ev6, ev20) + 12, hjust = 0, size = 3.4,
           label = sprintf("£%s: £%.0f", c("6,000", "20,000"), c(ev6, ev20))) +
  labs(x = "Willingness to pay (£000 per QALY)", y = "EVPI per patient (£)")


## -----------------------------------------------------------------------------
#| label: tbl-evppi
#| tbl-cap: "EVPPI by parameter group at £6,000 per QALY, estimated by GAM regression on the PSA draws."
evppi <- function(formula, lam = 6000) {
  dat <- bind_cols(pars, inb = lam * psa$dQ - psa$dC)
  g <- fitted(gam(update(formula, inb ~ .), data = dat))     # E[INB | phi]
  mean(pmax(g, 0)) - max(mean(g), 0)
}
groups <- list(`Treatment effects (HRs)` = ~ s(lhr_hosp) + s(lhr_cvd),
               `Baseline event risks`    = ~ s(p_hosp) + s(p_cvd) + s(p_cvd_post),
               `Utilities`               = ~ s(u_stable) + s(u_post),
               `Costs`                   = ~ s(c_hosp) + s(c_post) + s(c_stable))
ev <- tibble::tibble(Group = names(groups), evppi = map_dbl(groups, evppi))
ev |> arrange(desc(evppi)) |>
  mutate(`£ per patient` = round(evppi), `£ million` = round(evppi * pop / 1e6, 1),
         `% of EVPI` = round(100 * evppi / ev6),
         `Research it points to` = c(`Treatment effects (HRs)` = "Trial or updated meta-analysis",
                                     Costs = "Costing study", Utilities = "Preference study",
                                     `Baseline event risks` = "Cohort study")[Group]) |>
  select(-evppi) |> kable(align = "lcccl")


## -----------------------------------------------------------------------------
#| label: fig-xigris
#| fig-cap: "Cumulative meta-analysis of the three trials of drotrecogin alfa (28-day mortality, random effects after the first trial)."
#| fig-width: 7
#| fig-height: 2.8
xig <- tibble::tibble(study = c("PROWESS 2001", "ADDRESS 2005", "PROWESS-SHOCK 2012"),
                      ai = c(210, 247, 223), n1 = c(850, 1333, 846),
                      ci = c(259, 222, 202), n2 = c(840, 1307, 834))
cum <- map_dfr(1:3, function(k) {
  r <- rma(measure = "RR", ai = ai, n1i = n1, ci = ci, n2i = n2, data = xig[1:k, ],
           method = if (k == 1) "EE" else "DL")
  tibble::tibble(label = paste0("+ ", xig$study[k]), est = exp(r$b[1]), lo = exp(r$ci.lb), hi = exp(r$ci.ub))
}) |> mutate(label = forcats::fct_rev(forcats::fct_inorder(label)))
ggplot(cum, aes(est, label)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = navy) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.5)) +
  labs(x = "Relative risk of death (log scale)", y = NULL)


## -----------------------------------------------------------------------------
#| label: fig-ipi
#| fig-cap: "Illustrative survival curves calibrated to the published medians and 2-year survival of ipilimumab and control, extrapolated beyond 24 months under two assumptions about late mortality."
#| fig-width: 8
#| fig-height: 4
fit_cure <- function(med, s24) {
  f <- function(p) (p[1] + (1 - p[1]) * exp(-p[2] * med) - 0.5)^2 +
                   (p[1] + (1 - p[1]) * exp(-p[2] * 24) - s24)^2
  o <- optim(c(0.1, 0.1), f, method = "L-BFGS-B", lower = c(0, 0.001), upper = c(0.5, 1))$par
  list(pi = o[1], lam = o[2])
}
ctl <- fit_cure(6.4, 0.137); trt <- fit_cure(10.0, 0.216)
S <- function(t, p, h_late) {
  s24 <- p$pi + (1 - p$pi) * exp(-p$lam * 24)
  ifelse(t <= 24, p$pi + (1 - p$pi) * exp(-p$lam * t), s24 * exp(-h_late * (t - 24)))
}
tt <- 0:120
curves <- bind_rows(
  tibble::tibble(t = tt, arm = "Ipilimumab", A = S(tt, trt, 0.005), B = S(tt, trt, 0.06)),
  tibble::tibble(t = tt, arm = "Control",    A = S(tt, ctl, 0.005), B = S(tt, ctl, 0.06)))
gain <- function(col, upto) sum(curves[[col]][curves$arm == "Ipilimumab" & curves$t <= upto] -
                                curves[[col]][curves$arm == "Control" & curves$t <= upto])
g24 <- gain("A", 24); gA <- gain("A", 120); gB <- gain("B", 120)
curves |> pivot_longer(c(A, B), names_to = "extrap") |>
  mutate(extrap = recode(extrap, A = "2-year survivors cured (hazard 0.5%/month)",
                         B = "Disease continues (hazard 6%/month)")) |>
  ggplot(aes(t, value, colour = arm, linetype = extrap)) +
  geom_line(linewidth = 0.9) + geom_vline(xintercept = 24, colour = "grey50") +
  annotate("text", x = 26, y = 0.9, label = "end of\nobserved data", hjust = 0, size = 3) +
  scale_colour_manual(values = c(Control = "grey40", Ipilimumab = red), name = NULL) +
  labs(x = "Months", y = "Survival", linetype = NULL) +
  guides(linetype = guide_legend(nrow = 2)) + theme(legend.position = "top")

