# R code for chapter: 18-nma-in-practice.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch18
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
  library(readr)
})
fmt_ci <- \(m, l, h, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), m, l, h)
fmt_p  <- \(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))


## -----------------------------------------------------------------------------
#| label: published
published <- tibble::tribble(                         # Zhou et al., Figure 4A (SMD vs pill placebo)
  ~treat, ~smd, ~lo, ~hi,
  "Fluoxetine+CBT", -0.73, -1.39, -0.07, "Fluoxetine", -0.51, -0.84, -0.18,
  "Nortriptyline", 1.14, 0.46, 1.81, "WL", 0.99, 0.18, 1.79)


## -----------------------------------------------------------------------------
#| label: nma-tools
arm_data <- function(arms, ref) {
  arms |>
    mutate(.by = study, sd_pool = sqrt(sum((n - 1) * sd^2) / (sum(n) - n()))) |>
    mutate(y = mean / sd_pool, v = 1 / n, se = sqrt(v),
           trt = gsub("[^A-Za-z0-9]", "", treat),
           trt = relevel(factor(trt), gsub("[^A-Za-z0-9]", "", ref)),
           study = factor(study), arm = paste(study, trt))
}
nma_freq <- \(ad) rma.mv(y, v, mods = ~ 0 + study + trt, random = ~ 1 | arm, data = ad)
nma_tau  <- \(fit) sqrt(2 * fit$sigma2)
nma_effects <- function(fit, ad) {                # each treatment vs the reference, with 95% PI
  b <- coef(fit); keep <- startsWith(names(b), "trt"); tau2 <- 2 * fit$sigma2
  tibble::tibble(trt = sub("^trt", "", names(b)[keep]), est = b[keep],
                 se = sqrt(diag(vcov(fit)))[keep]) |>
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se,
           pi_lo = est - 1.96 * sqrt(se^2 + tau2), pi_hi = est + 1.96 * sqrt(se^2 + tau2)) |>
    left_join(distinct(ad, trt = as.character(trt), treat), by = "trt")
}
nma_contrast <- function(fit, A, B) {             # network estimate of A vs B
  b <- coef(fit); cv <- setNames(rep(0, length(b)), names(b))
  if (paste0("trt", A) %in% names(b)) cv[paste0("trt", A)] <- 1
  if (paste0("trt", B) %in% names(b)) cv[paste0("trt", B)] <- -1
  c(est = sum(cv * b), var = as.numeric(t(cv) %*% vcov(fit) %*% cv))
}
node_split <- function(ad, fit, A, B) {           # direct vs back-calculated indirect evidence
  A <- gsub("[^A-Za-z0-9]", "", A); B <- gsub("[^A-Za-z0-9]", "", B)
  dd <- ad |> filter(.by = study, all(c(A, B) %in% trt)) |> filter(trt %in% c(A, B)) |>
    mutate(trt = relevel(droplevels(trt), B), study = droplevels(study))
  form <- if (nlevels(dd$study) > 1) ~ 0 + study + trt else ~ trt
  fd   <- rma.mv(y, v, mods = form, random = ~ 1 | arm, data = dd, sigma2 = fit$sigma2)
  dir  <- coef(fd)[[paste0("trt", A)]]; vdir <- vcov(fd)[paste0("trt", A), paste0("trt", A)]
  net  <- nma_contrast(fit, A, B)
  vind <- 1 / (1 / net[["var"]] - 1 / vdir)
  ind  <- (net[["est"]] / net[["var"]] - dir / vdir) * vind
  tibble::tibble(comparison = paste(A, "vs", B), k_direct = n_distinct(dd$study),
                 direct = dir, se_dir = sqrt(vdir), indirect = ind, se_ind = sqrt(vind),
                 network = net[["est"]], p = 2 * pnorm(-abs((dir - ind) / sqrt(vdir + vind))))
}


## -----------------------------------------------------------------------------
#| label: repro
arms_posted <- read_csv("data/ch18/zhou_arms_posted.csv", show_col_types = FALSE)   # as posted
arms        <- read_csv("data/ch18/zhou_arms.csv", show_col_types = FALSE)          # corrected
ad_posted <- arm_data(arms_posted, "Pill-PBO"); f_posted <- nma_freq(ad_posted)
ad1 <- arm_data(arms, "Pill-PBO");              f1 <- nma_freq(ad1)
fx <- \(f, ad) with(filter(nma_effects(f, ad), treat == "Fluoxetine"), fmt_ci(est, lo, hi))
audit <- full_join(select(arms_posted, study, treat, m_post = mean, n_post = n),
                   select(arms, study, treat, m_corr = mean, n_corr = n), by = c("study", "treat")) |>
  summarise(.by = study, bad = any(abs(m_post - m_corr) > 0.01 | n_post != n_corr, na.rm = TRUE),
            pbo = any(treat == "Pill-PBO"))


## -----------------------------------------------------------------------------
#| label: tbl-repro
#| tbl-cap: "Reproducing Zhou et al.: published estimates, the frequentist arm-level NMA (metafor), and the Bayesian version (brms posterior draws). SMD vs pill placebo; the prediction interval is for a new trial."
draws <- read_csv("data/ch18/zhou_draws.csv", show_col_types = FALSE)
d_full <- filter(draws, fit == "zhou_full")
bq <- function(d, t, pred = FALSE) {
  x <- d[[paste0("b_trt", gsub("[^A-Za-z0-9]", "", t))]]
  if (pred) { set.seed(629); x <- x + rnorm(length(x), 0, sqrt(2) * d$sd_arm__Intercept) }
  fmt_ci(median(x), quantile(x, 0.025), quantile(x, 0.975))
}
e1 <- nma_effects(f1, ad1)
published |>
  transmute(Treatment = treat, Published = fmt_ci(smd, lo, hi),
            `Frequentist (metafor)` = map_chr(treat, \(t) with(filter(e1, treat == t), fmt_ci(est, lo, hi))),
            `Bayesian (brms)` = map_chr(treat, \(t) bq(d_full, t)),
            `95% prediction interval (brms)` = map_chr(treat, \(t) bq(d_full, t, pred = TRUE))) |>
  kable(align = "lcccc")


## -----------------------------------------------------------------------------
#| label: nma-brms
#| eval: false
# library(brms)
# nma_priors <- c(prior(normal(0, 10), class = b), prior(normal(0, 1), class = sd))
# b1 <- brm(y | se(se) ~ 0 + study + trt + (1 | arm), data = ad1, prior = nma_priors,
#           chains = 4, iter = 3000, warmup = 1000, seed = 629, control = list(adapt_delta = 0.95))


## -----------------------------------------------------------------------------
#| label: contrasts
make_contrasts <- function(arms) {             # each arm vs its trial's first arm
  arms |>
    mutate(.by = study, ref = first(treat), mref = first(mean), sref = first(sd),
           nref = first(n), N = sum(n), k = row_number()) |>
    filter(k > 1) |>
    escalc(measure = "SMD", m1i = mean, sd1i = sd, n1i = n, m2i = mref, sd2i = sref, n2i = nref, data = _)
}
contr <- make_contrasts(arms)
att_c <- filter(contr, study == "Attari 2006")
tau2_a <- 2 * f1$sigma2
c(SMD = att_c$yi, `variance, arm-level` = 1 / 20 + 1 / 20, `variance, contrast-level` = att_c$vi,
  `weight ratio, arm vs contrast` = (1 / (0.10 + tau2_a)) / (1 / (att_c$vi + tau2_a))) |> round(2)


## -----------------------------------------------------------------------------
#| label: fig-flx-pairwise
#| fig-cap: "Pairwise meta-analysis of the nine fluoxetine-placebo trials (direct evidence only; HKSJ confidence interval). Grey band: ±0.35 SD. Red: 95% prediction interval."
#| fig-width: 8
#| fig-height: 4.2
pw <- arms |>
  filter(.by = study, all(c("Fluoxetine", "Pill-PBO") %in% treat)) |>
  reframe(.by = study,
          m1 = mean[treat == "Fluoxetine"], s1 = sd[treat == "Fluoxetine"], n1 = n[treat == "Fluoxetine"],
          m2 = mean[treat == "Pill-PBO"], s2 = sd[treat == "Pill-PBO"], n2 = n[treat == "Pill-PBO"]) |>
  escalc(measure = "SMD", m1i = m1, sd1i = s1, n1i = n1, m2i = m2, sd2i = s2, n2i = n2, data = _)
pw_fit <- rma(yi, vi, data = pw, method = "REML", test = "knha"); pw_pr <- predict(pw_fit)
bind_rows(
  as.data.frame(pw) |> transmute(label = study, est = yi, lo = yi - 1.96 * sqrt(vi), hi = yi + 1.96 * sqrt(vi), type = "Trial"),
  tibble::tibble(label = c("Pooled (95% CI)", "95% prediction interval"), est = pw_pr$pred,
                 lo = c(pw_pr$ci.lb, pw_pr$pi.lb), hi = c(pw_pr$ci.ub, pw_pr$pi.ub), type = c("Pooled", "PI"))) |>
  mutate(label = factor(label, levels = rev(unique(label)))) |>
  ggplot(aes(est, label, colour = type)) +
  annotate("rect", xmin = -0.35, xmax = 0.35, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi)) +
  scale_colour_manual(values = c(Trial = "#041E42", Pooled = "#041E42", PI = "#B22222"), guide = "none") +
  labs(x = "SMD, fluoxetine vs placebo (< 0 favours fluoxetine)", y = NULL)


## -----------------------------------------------------------------------------
#| label: tbl-nodesplit
#| tbl-cap: "Node-splitting in the fluoxetine-nortriptyline-placebo loop (frequentist back-calculation). SMD (SE); < 0 favours the first-named treatment."
loop <- list(c("Fluoxetine", "Pill-PBO"), c("Nortriptyline", "Pill-PBO"), c("Nortriptyline", "Fluoxetine"))
map_dfr(loop, \(x) node_split(ad1, f1, x[1], x[2])) |>
  transmute(Comparison = comparison, `Direct trials` = k_direct,
            Direct = sprintf("%.2f (%.2f)", direct, se_dir), Indirect = sprintf("%.2f (%.2f)", indirect, se_ind),
            Network = sprintf("%.2f", network), `P (direct = indirect)` = fmt_p(p)) |>
  kable(align = "lccccc")


## -----------------------------------------------------------------------------
#| label: fig-bayes-split
#| fig-cap: "Bayesian node-split, fluoxetine vs placebo: posterior distributions of the direct and indirect SMDs, with all trials and without Attari 2006."
#| fig-width: 8
#| fig-height: 3.8
sp <- read_csv("data/ch18/zhou_split_draws.csv", show_col_types = FALSE) |>
  mutate(Analysis = recode(fit, zhou_split = "All trials", zhou_split_noatt = "Attari 2006 excluded"),
         Indirect = b_trtFluoxetine, Direct = b_trtFluoxetine + b_split)
sp |> pivot_longer(c(Direct, Indirect), names_to = "Evidence", values_to = "smd") |>
  ggplot(aes(smd, fill = Evidence, colour = Evidence)) +
  geom_density(alpha = 0.3) + geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ Analysis) +
  scale_fill_manual(values = c(Direct = "#041E42", Indirect = "#B22222")) +
  scale_colour_manual(values = c(Direct = "#041E42", Indirect = "#B22222")) +
  coord_cartesian(xlim = c(-3, 1)) +
  labs(x = "SMD, fluoxetine vs placebo", y = "Posterior density") + theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: fig-outlier
#| fig-cap: "Every within-trial comparison in the network (each arm vs its trial's first arm), by SMD and standard error. The sign reflects only which arm is listed first."
#| fig-width: 7
#| fig-height: 4
as.data.frame(contr) |>
  mutate(se = sqrt(vi), attari = study == "Attari 2006") |>
  ggplot(aes(yi, se, colour = attari)) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 2.3, alpha = 0.8) +
  annotate("text", x = 4.1, y = 0.5, hjust = 1, colour = "#B22222", label = "Attari 2006\n(n = 40)") +
  scale_y_reverse() + scale_colour_manual(values = c(`FALSE` = "#041E42", `TRUE` = "#B22222"), guide = "none") +
  labs(x = "Trial-level SMD", y = "Standard error")


## -----------------------------------------------------------------------------
#| label: loo
loo <- map_dfr(levels(ad1$study), \(s) {
  d <- filter(ad1, study != s) |> mutate(study = droplevels(study), trt = droplevels(trt))
  if (!all(c("Fluoxetine", "PillPBO") %in% d$trt)) return(NULL)
  f <- nma_freq(d); e <- nma_contrast(f, "Fluoxetine", "PillPBO")
  tibble::tibble(omitted = s, est = e[["est"]], tau = nma_tau(f),
                 p_incons = node_split(d, f, "Fluoxetine", "Pill-PBO")$p)
})
others <- filter(loo, omitted != "Attari 2006"); att <- filter(loo, omitted == "Attari 2006")


## -----------------------------------------------------------------------------
#| label: fig-attari-r
#| fig-cap: "Attari 2006: the change-score SMD implied by the reported baseline and week-8 means and SDs, for each assumed pre-post correlation r. Only r ≈ 0.97 reproduces the SMD that entered both Lancet NMAs."
#| fig-width: 6.5
#| fig-height: 3.8
att <- list(bF = 28.9, sbF = 8.4, eF = 17.9, seF = 8.1,      # CDI, mean (SD), 20 per arm
            bN = 28.4, sbN = 8.7, eN = 25.8, seN = 8.6)
sd_change <- \(sb, se, r) sqrt(sb^2 + se^2 - 2 * r * sb * se)
tibble::tibble(r = seq(0, 0.99, 0.01)) |>
  mutate(smd = ((att$bF - att$eF) - (att$bN - att$eN)) /
               sqrt((sd_change(att$sbF, att$seF, r)^2 + sd_change(att$sbN, att$seN, r)^2) / 2)) |>
  ggplot(aes(r, smd)) +
  geom_line(linewidth = 1, colour = "#041E42") +
  geom_hline(yintercept = 4.24, colour = "#B22222", linetype = "dashed") +
  annotate("text", x = 0.02, y = 4.45, hjust = 0, colour = "#B22222", label = "SMD used in the NMAs") +
  labs(x = "Assumed pre-post correlation r", y = "Implied change-score SMD")


## -----------------------------------------------------------------------------
#| label: ancova
anc_attari <- function(r) {                         # nortriptyline vs fluoxetine, ANCOVA-type SMD
  sd_e <- sqrt((att$seF^2 + att$seN^2) / 2); sd_b <- sqrt((att$sbF^2 + att$sbN^2) / 2)
  d <- ((att$eN - att$eF) - r * sd_e / sd_b * (att$bN - att$bF)) / sd_e
  c(d = d, v = (2 / 20) * (1 - r^2) + d^2 / 80)
}
sapply(c(0.3, 0.5, 0.7), anc_attari) |> round(2)


## -----------------------------------------------------------------------------
#| label: fragility
fragility_index <- function(e1, n1, e2, n2, alpha = 0.05) {
  p <- \(a, c) fisher.test(matrix(c(a, n1 - a, c, n2 - c), 2))$p.value
  if (p(e1, e2) >= alpha) return(NA_integer_)
  k <- 0
  while (p(e1, e2) < alpha) { if (e1 / n1 < e2 / n2) e1 <- e1 + 1 else e2 <- e2 + 1; k <- k + 1 }
  k
}
fragility_index(10, 20, 2, 20)


## -----------------------------------------------------------------------------
#| label: tbl-nma-sens
#| tbl-cap: "Zhou et al. NMA, with Attari 2006 as published, corrected and excluded (frequentist arm-level NMA). SMD vs pill placebo."
ad2 <- arm_data(filter(arms, study != "Attari 2006"), "Pill-PBO"); f2 <- nma_freq(ad2)
ad3 <- ad1; a <- anc_attari(0.5)
ad3$y[ad3$arm == "Attari 2006 Fluoxetine"] <- 0
ad3$y[ad3$arm == "Attari 2006 Nortriptyline"] <- a[["d"]]
ad3$v[ad3$study == "Attari 2006"] <- a[["v"]] / 2; ad3$se <- sqrt(ad3$v)
f3 <- nma_freq(ad3)
sens <- list(`All trials (as published)` = list(f = f1, ad = ad1), `Attari corrected (ANCOVA)` = list(f = f3, ad = ad3),
             `Attari excluded` = list(f = f2, ad = ad2))
imap_dfr(sens, \(s, nm) {
  e <- nma_effects(s$f, s$ad); g <- \(t) with(filter(e, treat == t), fmt_ci(est, lo, hi))
  ef <- filter(e, treat == "Fluoxetine")
  tibble::tibble(Analysis = nm, Fluoxetine = g("Fluoxetine"), `Fluoxetine 95% PI` = sprintf("%.2f to %.2f", ef$pi_lo, ef$pi_hi),
                 Nortriptyline = g("Nortriptyline"), τ = sprintf("%.2f", nma_tau(s$f)),
                 `Inconsistency P` = fmt_p(node_split(s$ad, s$f, "Fluoxetine", "Pill-PBO")$p))
}) |> kable(align = "lccccc")


## -----------------------------------------------------------------------------
#| label: tbl-nma-probs
#| tbl-cap: "Bayesian NMA: fluoxetine vs pill placebo, with the posterior probabilities that the SMD is below 0 (any benefit) and below −0.35 (a clinically important benefit)."
lab <- c(zhou_full = "All trials (as published)", zhou_anc = "Attari corrected (ANCOVA)", zhou_noatt = "Attari excluded")
draws |> group_by(fit) |>
  summarise(`SMD (95% CrI)` = fmt_ci(median(b_trtFluoxetine), quantile(b_trtFluoxetine, 0.025), quantile(b_trtFluoxetine, 0.975)),
            `P(SMD < 0)` = sprintf("%.3f", mean(b_trtFluoxetine < 0)),
            `P(SMD < −0.35)` = sprintf("%.2f", mean(b_trtFluoxetine < -0.35))) |>
  mutate(fit = lab[fit]) |> rename(Analysis = fit) |>
  arrange(match(Analysis, lab)) |> kable(align = "lccc")


## -----------------------------------------------------------------------------
#| label: tbl-zhou-sucra
#| tbl-cap: "SUCRA and 95% rank intervals (among 32 treatments, 1 = best) for antidepressants in the Bayesian NMA."
sucra <- function(d) {
  m <- as.matrix(select(d, starts_with("b_trt"))); colnames(m) <- sub("^b_trt", "", colnames(m))
  m <- cbind(PillPBO = 0, m); rk <- t(apply(m, 1, rank))
  tibble::tibble(trt = colnames(m), SUCRA = (ncol(m) - colMeans(rk)) / (ncol(m) - 1),
                 rank_ci = apply(rk, 2, \(r) sprintf("%d to %d", as.integer(quantile(r, 0.025)), as.integer(quantile(r, 0.975)))))
}
drugs <- c("Fluoxetine", "Duloxetine", "Venlafaxine", "Mirtazapine", "Sertraline", "Paroxetine", "Escitalopram", "Citalopram")
left_join(sucra(d_full) |> rename(`SUCRA, all trials` = SUCRA, `Rank, 95% interval` = rank_ci),
          sucra(filter(draws, fit == "zhou_noatt")) |> select(trt, `SUCRA, Attari excluded` = SUCRA), by = "trt") |>
  filter(trt %in% drugs) |> arrange(desc(`SUCRA, all trials`)) |> rename(Antidepressant = trt) |>
  kable(digits = 2)


## -----------------------------------------------------------------------------
#| label: tbl-lyus
#| tbl-cap: "Fluoxetine vs placebo in three NMAs, with and without Attari 2006, re-analysed with the arm-level NMA of this chapter (SMD, 95% CI)."
cip <- read_csv("data/ch18/cipriani_arms.csv", show_col_types = FALSE)
het <- read_csv("data/ch18/hetrick_arms.csv", show_col_types = FALSE)
cip_arms <- filter(cip, !is.na(mean)) |> select(study, treat, mean, sd, n)
adc1 <- arm_data(cip_arms, "Placebo"); fc1 <- nma_freq(adc1)
adc2 <- arm_data(filter(cip_arms, study != "Attari 2006"), "Placebo"); fc2 <- nma_freq(adc2)
adh  <- arm_data(het, "Placebo"); fh <- nma_freq(adh)
tibble::tribble(~NMA, ~Analysis, ~fit, ~ad,
  "Cipriani 2016", "All trials", fc1, adc1, "Cipriani 2016", "Attari excluded", fc2, adc2,
  "Zhou 2020", "All trials", f1, ad1, "Zhou 2020", "Attari excluded", f2, ad2,
  "Hetrick 2021", "No Attari (by design)", fh, adh) |>
  mutate(`Fluoxetine vs placebo` = map2_chr(fit, ad, fx),
         `Inconsistency P` = map2_chr(fit, ad, \(f, a) {
           ref <- levels(a$trt)[1]; fmt_p(node_split(a, f, "Fluoxetine", ref)$p) })) |>
  select(-fit, -ad) |> kable(align = "llcc")


## -----------------------------------------------------------------------------
#| label: response
resp <- cip |>
  filter(.by = study, all(c("Fluoxetine", "Placebo") %in% treat)) |>
  filter(treat %in% c("Fluoxetine", "Placebo"), !is.na(responders)) |>
  reframe(.by = study, ai = responders[treat == "Fluoxetine"], n1i = responders_n[treat == "Fluoxetine"],
          ci = responders[treat == "Placebo"], n2i = responders_n[treat == "Placebo"])
rr <- rma(measure = "RR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = resp, method = "REML", test = "knha")
predict(rr, transf = exp)

