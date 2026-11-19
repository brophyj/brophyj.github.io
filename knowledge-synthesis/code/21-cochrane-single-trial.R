# R code for chapter: 21-cochrane-single-trial.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch21
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
})
clr_prior <- "#3182bd"; clr_lik <- "#e6550d"; clr_post <- "#31a354"; clr_sig <- "#756bb1"


## -----------------------------------------------------------------------------
#| label: fig-living
#| fig-cap: "A simulated living review: 12 trials of a treatment with a true OR of 0.75, pooled (common effect) after each new trial. Whether the 95% CI excludes 1 changes from one update to the next."
#| fig-width: 7
#| fig-height: 4.2
set.seed(42)
true_or  <- 0.75; n_trials <- 12
trial_n  <- sample(80:400, n_trials)
trial_se <- sqrt(4 / trial_n)                       # approximate SE of a log OR
trial_y  <- rnorm(n_trials, log(true_or), 0.25)     # observed log ORs
cumulative <- map_dfr(1:n_trials, function(k) {
  m <- rma(trial_y[1:k], sei = trial_se[1:k], method = "EE")
  tibble::tibble(update = k, pooled = exp(coef(m)), lo = exp(m$ci.lb), hi = exp(m$ci.ub),
                 sig = m$ci.ub < 0)
})
ggplot(cumulative, aes(update, pooled, colour = sig)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = true_or, linetype = 3, colour = "#041E42") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.9) +
  geom_point(size = 2.5) +
  annotate("text", x = 12, y = true_or - 0.04, label = "true OR = 0.75", hjust = 1, colour = "#041E42") +
  scale_colour_manual(values = c(`FALSE` = "#de2d26", `TRUE` = clr_post),
                      labels = c("CI includes 1", "CI excludes 1"), name = NULL) +
  scale_x_continuous(breaks = 1:12) +
  labs(x = "Update (trials included so far)", y = "Pooled OR (95% CI)") +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: decaf-inputs
HR <- 0.61; CI_lo <- 0.42; CI_hi <- 0.89
y_obs <- log(HR)                                        # observed log HR, coffee vs abstinence
se    <- (log(CI_hi) - log(CI_lo)) / (2 * qnorm(0.975))  # its standard error
z_obs <- y_obs / se
round(c(log_HR = y_obs, se = se, se_from_events = sqrt(1 / 47 + 1 / 64), z = z_obs), 3)


## -----------------------------------------------------------------------------
#| label: fig-sigfilter
#| fig-cap: "Risk differences in 10,000 simulated trials with 100 patients per arm and a true difference of −7.5 percentage points (green line). Trials reaching p < 0.05 are in purple; DECAF's observed difference is the dashed line."
#| fig-width: 7
#| fig-height: 4
set.seed(629)
p0 <- 0.50; p1 <- 0.425; n <- 100
sims <- tibble::tibble(e0 = rbinom(1e4, n, p0), e1 = rbinom(1e4, n, p1)) |>
  mutate(rd  = (e1 - e0) / n,
         se  = sqrt(e1 / n * (1 - e1 / n) / n + e0 / n * (1 - e0 / n) / n),
         sig = abs(rd / se) > 1.96)
ggplot(sims, aes(rd, fill = sig)) +
  geom_histogram(binwidth = 0.01, colour = "white", boundary = 0) +
  geom_vline(xintercept = p1 - p0, linewidth = 1, colour = clr_post) +
  geom_vline(xintercept = -0.17, linetype = 2, colour = clr_lik, linewidth = 1) +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = clr_sig),
                    labels = c("p > 0.05", "p < 0.05"), name = NULL) +
  labs(x = "Observed risk difference", y = "Simulated trials") + theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: retrodesign
retrodesign <- function(A, s, alpha = 0.05, n.sims = 1e5) {
  z     <- qnorm(1 - alpha / 2)
  p.hi  <- 1 - pnorm(z - A / s)              # significant, correct sign
  p.lo  <- pnorm(-z - A / s)                 # significant, wrong sign
  power <- p.hi + p.lo
  est   <- A + s * rnorm(n.sims)             # replicate estimates
  sig   <- abs(est) > s * z
  c(power = power, type_s = p.lo / power, type_m = mean(abs(est[sig])) / A)
}
s_design <- sqrt(0.5 * 0.5 / 100 + 0.425 * 0.575 / 100)
set.seed(1)
retro <- rbind(`Design: RD 20.5%`    = retrodesign(0.205, s_design),
               `Realistic: RD 7.5%`  = retrodesign(0.075, s_design))
round(retro, 2)


## -----------------------------------------------------------------------------
#| label: conjugate
m0 <- 0.871; s0 <- 0.5
w0 <- 1 / s0^2; w <- 1 / se^2
m1 <- (w0 * m0 + w * y_obs) / (w0 + w); s1 <- sqrt(1 / (w0 + w))
ma <- rma(yi = c(m0, y_obs), sei = c(s0, se), method = "EE")     # the same, as a meta-analysis
rbind(closed_form = exp(c(m1, m1 - 1.96 * s1, m1 + 1.96 * s1)),
      rma         = exp(c(ma$b, ma$ci.lb, ma$ci.ub))) |> round(2) |> `colnames<-`(c("HR", "lower", "upper"))
setNames(round(pnorm(log(c(1, 0.9, 0.8)), m1, s1), 2), c("P(HR<1)", "P(HR<0.9)", "P(HR<0.8)"))


## -----------------------------------------------------------------------------
#| label: fig-decaf-hr
#| fig-cap: "Design prior (centred on harm from coffee), DECAF likelihood, and posterior for the HR of coffee versus abstinence. Shaded: P(HR < 0.9)."
#| fig-width: 8
#| fig-height: 3.8
th <- seq(-1.3, 2.4, length.out = 800)
dd <- bind_rows(
  tibble::tibble(th, d = dnorm(th, m0, s0),    what = "Prior (design: abstinence better)"),
  tibble::tibble(th, d = dnorm(th, y_obs, se), what = "Likelihood (DECAF)"),
  tibble::tibble(th, d = dnorm(th, m1, s1),    what = "Posterior")) |>
  mutate(what = factor(what, levels = unique(what)))
ggplot(dd, aes(th, d, colour = what)) +
  geom_area(data = filter(dd, what == "Posterior", th < log(0.9)), aes(th, d),
            fill = clr_post, alpha = 0.25, colour = NA, inherit.aes = FALSE) +
  geom_line(linewidth = 1.1) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(clr_prior, clr_lik, clr_post), name = NULL) +
  scale_x_continuous(breaks = log(c(0.3, 0.5, 0.75, 1, 1.5, 2.5, 5, 10)),
                     labels = c(0.3, 0.5, 0.75, 1, 1.5, 2.5, 5, 10)) +
  labs(x = "HR, coffee vs abstinence (log scale)", y = "Density") + theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: decaf-brms
#| eval: false
# library(brms)
# decaf <- tibble::tibble(abst = c(0, 1), events = c(47, 64), n = c(100, 100))
# pri_rd <- c(prior(normal(0, 1.5), class = b, coef = Intercept),
#             prior(normal(-0.871, 0.5), class = b, coef = abst))
# fit_rd <- brm(events | trials(n) ~ 0 + Intercept + abst, data = decaf, family = binomial,
#               prior = pri_rd, chains = 4, iter = 6000, warmup = 1000, seed = 2026)
# draws <- as_draws_df(fit_rd)[, c("b_Intercept", "b_abst")]


## -----------------------------------------------------------------------------
#| label: decaf-rd
draws <- read.csv("data/ch21/decaf_rd_draws.csv")
rd <- plogis(draws$b_Intercept) - plogis(draws$b_Intercept + draws$b_abst)
round(100 * c(mean = mean(rd), quantile(rd, c(0.025, 0.975))), 1)
round(c(`P(RD < 0)` = mean(rd < 0), `P(RD < -2%)` = mean(rd < -0.02)), 2)


## -----------------------------------------------------------------------------
#| label: tbl-decaf-sens
#| tbl-cap: "DECAF combined with five priors for the HR of coffee versus abstinence (normal approximation)."
post_norm <- function(m0, s0) {
  w0 <- 1 / s0^2; m <- (w0 * m0 + w * y_obs) / (w0 + w); s <- sqrt(1 / (w0 + w))
  tibble::tibble(post_hr = exp(m), post_lo = exp(m - 1.96 * s), post_hi = exp(m + 1.96 * s),
                 p1 = pnorm(0, m, s), p08 = pnorm(log(0.8), m, s))
}
hr_design <- log(1 - 0.295) / log(1 - 0.50)        # abstinence vs coffee, proportional hazards
sens <- tibble::tibble(
  Prior = c("1. Vague: N(0, 10²)", "2. Sceptical: N(0, 0.5²)",
            sprintf("3. Design as an HR: N(%.3f, 0.5²)", log(1 / hr_design)),
            "4. Design as an OR: N(0.871, 0.5²)", "5. Design, held firmly: N(0.871, 0.25²)"),
  m0 = c(0, 0, log(1 / hr_design), 0.871, 0.871), s0 = c(10, 0.5, 0.5, 0.5, 0.25)) |>
  mutate(prior_p1 = pnorm(0, m0, s0), post = map2(m0, s0, post_norm)) |> unnest(post)
sens |>
  transmute(Prior, `Prior P(HR < 1)` = sprintf("%.2f", prior_p1),
            `Posterior HR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", post_hr, post_lo, post_hi),
            `P(HR < 1)` = sprintf("%.2f", p1), `P(HR < 0.8)` = sprintf("%.2f", p08)) |>
  kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: eb
sd_z  <- c(1.17, 1.74, 2.38, 5.73); p_mix <- c(0.32, 0.31, 0.30, 0.07)
sd_snr <- sqrt(sd_z^2 - 1)                         # remove the unit noise
post_w <- p_mix * dnorm(z_obs, 0, sd_z); post_w <- post_w / sum(post_w)
Bk <- sd_snr^2 / (sd_snr^2 + 1)
set.seed(2026); N <- 2e5
comp   <- sample.int(4, N, TRUE, post_w)
snr_EB <- rnorm(N, Bk[comp] * z_obs, sqrt(Bk[comp]))
hr_EB  <- exp(snr_EB * se)                         # back to the HR scale
rbind(prior_weight = p_mix, sd_true_snr = sd_snr, posterior_weight = post_w, shrinkage_B = Bk) |>
  round(2) |> `colnames<-`(paste0("k=", 1:4))
round(c(HR = median(hr_EB), quantile(hr_EB, c(0.025, 0.975)),
        `P(HR<1)` = mean(hr_EB < 1), `P(HR<0.8)` = mean(hr_EB < 0.8)), 2)


## -----------------------------------------------------------------------------
#| label: tbl-four-priors
#| tbl-cap: "DECAF under four priors on the signal-to-noise scale. Type S: posterior probability that the true effect is on the other side of 1. Type M: median factor by which |z| overstates the true |SNR|. Replication: probability that a same-size repeat is significant in the same direction."
normal_post_snr <- function(mu0, tau) {
  B <- tau^2 / (tau^2 + 1); rnorm(N, B * z_obs + (1 - B) * mu0, sqrt(B))
}
set.seed(2026)
snr <- list(`Empirical Bayes (Cochrane)` = snr_EB,
            `Vague N(0, 10²)`            = normal_post_snr(0, 10),
            `Sceptical N(0, 1²)`         = normal_post_snr(0, 1),
            `Design prior`               = normal_post_snr(0.871 / se, 0.5 / se))
imap_dfr(snr, function(s, nm) {
  z2 <- s + rnorm(length(s))
  tibble::tibble(`Prior on SNR` = nm,
                 `Posterior HR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", median(exp(s * se)),
                                                    quantile(exp(s * se), .025), quantile(exp(s * se), .975)),
                 `Type S` = round(mean(sign(s) != sign(z_obs)), 3),
                 `Type M` = round(median(abs(z_obs / s)), 2),
                 Replication = round(mean(z2 < -1.96), 2))
}) |> kable(align = "lcccc")

