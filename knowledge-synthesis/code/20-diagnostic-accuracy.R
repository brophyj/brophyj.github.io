# R code for chapter: 20-diagnostic-accuracy.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch20
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
  library(mada)      # bivariate model for diagnostic accuracy: reitsma()
})


## -----------------------------------------------------------------------------
#| label: fig-threshold
#| fig-cap: "Simulated troponin concentrations in 1000 patients with MI and 1000 without (left), and the ROC curve traced by moving the threshold (right). Each marked point is one threshold."
#| fig-width: 10
#| fig-height: 4
#| out-width: "100%"
set.seed(42)
trop_d  <- rlnorm(1000, 3.0, 0.6)          # MI
trop_nd <- rlnorm(1000, 1.8, 0.5)          # no MI
ops <- tibble::tibble(t = c(5, 14, 40)) |>
  mutate(sens = map_dbl(t, \(x) mean(trop_d >= x)), fpr = map_dbl(t, \(x) mean(trop_nd >= x)),
         lab = sprintf("%g ng/L\nsens %.2f, spec %.2f", t, sens, 1 - fpr))
roc <- tibble::tibble(t = sort(unique(c(0, trop_d, trop_nd, Inf)))) |>
  mutate(sens = map_dbl(t, \(x) mean(trop_d >= x)), fpr = map_dbl(t, \(x) mean(trop_nd >= x)))
p1 <- tibble::tibble(trop = c(trop_d, trop_nd), group = rep(c("MI", "No MI"), each = 1000)) |>
  ggplot(aes(trop, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = ops$t, linetype = c(2, 1, 3)) +
  scale_x_continuous(limits = c(0, 120)) +
  scale_fill_manual(values = c(MI = "#B22222", `No MI` = "#041E42"), name = NULL) +
  labs(x = "Troponin (ng/L)", y = "Patients") + theme(legend.position = "top")
p2 <- ggplot(roc, aes(fpr, sens)) +
  geom_abline(linetype = 2, colour = "grey60") +
  geom_path(colour = "#041E42", linewidth = 1) +
  geom_point(data = ops, colour = c("#2E7D32", "#E69F00", "#B22222"), size = 3) +
  geom_text(data = ops, aes(label = lab), hjust = -0.1, vjust = 1, size = 3, lineheight = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "1 − specificity (false-positive rate)", y = "Sensitivity")
patchwork::wrap_plots(p1, p2)


## -----------------------------------------------------------------------------
#| label: tbl-posttest
#| tbl-cap: "Post-test probabilities for a test with sensitivity 0.98 and specificity 0.75, at three pre-test probabilities."
post <- \(p, lr) { o <- p / (1 - p) * lr; o / (1 + o) }
lr_pos <- 0.98 / 0.25; lr_neg <- 0.02 / 0.75
tibble::tibble(`Pre-test probability` = c(0.05, 0.15, 0.40)) |>
  mutate(`After a positive test` = post(`Pre-test probability`, lr_pos),
         `After a negative test` = post(`Pre-test probability`, lr_neg),
         across(everything(), \(x) sprintf("%.1f%%", 100 * x))) |>
  kable(align = "ccc")


## -----------------------------------------------------------------------------
#| label: tbl-verification
#| tbl-cap: "One hypothetical study of 1000 patients (150 with MI; true sensitivity 0.90, specificity 0.80) under three verification schemes. Only 25% of test-negatives receive the full reference standard in the last two."
tp <- 135; fn <- 15; fp <- 170; tn <- 680; v <- 0.25
tibble::tibble(Verification = c("Full", "Partial: unverified patients excluded",
                                "Differential: unverified negatives assumed MI-free"),
               TP = tp, FN = c(fn, v * fn, v * fn), FP = fp, TN = c(tn, v * tn, tn + (1 - v) * fn)) |>
  mutate(Sensitivity = TP / (TP + FN), Specificity = TN / (TN + FP)) |>
  kable(digits = 2, align = "lcccccc")


## -----------------------------------------------------------------------------
#| label: dta-sim
set.seed(2026)
k <- 14
cut  <- sample(c(3, 5, 7, 9), k, replace = TRUE)
n    <- round(runif(k, 400, 2500)); prev <- runif(k, 0.08, 0.20)
u    <- MASS::mvrnorm(k, c(0, 0), matrix(c(0.30, 0.05, 0.05, 0.12), 2))
lsens <- qlogis(0.97) - 0.30 * (cut - 5) + u[, 1]
lfpr  <- qlogis(0.21) - 0.18 * (cut - 5) + u[, 2]
dis <- rbinom(k, n, prev)
TP <- rbinom(k, dis, plogis(lsens)); FN <- dis - TP
FP <- rbinom(k, n - dis, plogis(lfpr)); TN <- n - dis - FP
dta <- data.frame(study = paste0("S", 1:k), cut, TP, FN, FP, TN)
fit <- reitsma(dta)                                    # bivariate model, REML
sm  <- summary(fit)$coefficients
tibble::tibble(Parameter = c("Sensitivity", "Specificity", "SD logit sensitivity", "SD logit FPR", "Correlation ρ"),
               Estimate = c(sprintf("%.3f (%.3f to %.3f)", sm[3, 1], sm[3, 5], sm[3, 6]),
                            sprintf("%.3f (%.3f to %.3f)", 1 - sm[4, 1], 1 - sm[4, 6], 1 - sm[4, 5]),
                            sprintf("%.2f", sqrt(fit$Psi[1, 1])), sprintf("%.2f", sqrt(fit$Psi[2, 2])),
                            sprintf("%.2f", cov2cor(fit$Psi)[1, 2]))) |>
  kable(align = "lc")


## -----------------------------------------------------------------------------
#| label: fig-sroc
#| fig-cap: "The bivariate model for 14 simulated hs-cTn studies: study estimates (points, coloured by threshold), summary point (diamond), 95% confidence region (shaded), 95% prediction region (dashed) and summary ROC curve (red)."
#| fig-width: 7
#| fig-height: 5.5
mu  <- coef(fit)[1, ]; V <- vcov(fit); Psi <- fit$Psi
sw  <- \(m) m[2:1, 2:1]                                # order (FPR, sens)
conf <- as.data.frame(plogis(ellipse::ellipse(sw(V), centre = mu[2:1], level = 0.95)))
pred <- as.data.frame(plogis(ellipse::ellipse(sw(Psi + V), centre = mu[2:1], level = 0.95)))
names(conf) <- names(pred) <- c("fpr", "sens")
srocc <- as.data.frame(sroc(fit, fpr = seq(0.01, 0.6, length.out = 300))); names(srocc) <- c("fpr", "sens")
ggplot() +
  geom_polygon(data = pred, aes(fpr, sens), fill = NA, colour = "#041E42", linetype = 2) +
  geom_polygon(data = conf, aes(fpr, sens), fill = "#041E42", alpha = 0.25) +
  geom_line(data = srocc, aes(fpr, sens), colour = "#B22222", linewidth = 1.1) +
  geom_point(data = dta, aes(FP / (FP + TN), TP / (TP + FN), colour = factor(cut)), size = 2.5) +
  annotate("point", x = plogis(mu[2]), y = plogis(mu[1]), shape = 18, size = 5, colour = "#B22222") +
  scale_colour_manual(values = c("#0072B2", "#009E73", "#E69F00", "#B22222"), name = "Threshold (ng/L)") +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0.5, 1)) +
  labs(x = "1 − specificity", y = "Sensitivity")


## -----------------------------------------------------------------------------
#| label: hsroc
sA <- sqrt(Psi[1, 1]); sB <- sqrt(Psi[2, 2]); sAB <- -Psi[1, 2]   # A = sensitivity, B = specificity
muA <- mu[[1]]; muB <- -mu[[2]]; r <- sB / sA
hs <- c(beta = log(r), Lambda = sqrt(r) * muA + sqrt(1 / r) * muB,
        Theta = (sqrt(r) * muA - sqrt(1 / r) * muB) / 2,
        sd_accuracy = sqrt(2 * (sA * sB + sAB)), sd_threshold = sqrt((sA * sB - sAB) / 2))
round(hs, 2)


## -----------------------------------------------------------------------------
#| label: dta-clinical
s_sum <- plogis(mu[[1]]); f_sum <- plogis(mu[[2]])
lrn <- (1 - s_sum) / (1 - f_sum)


## -----------------------------------------------------------------------------
#| label: dta-metareg
fit_cut <- reitsma(dta, formula = cbind(tsens, tfpr) ~ cut)
b <- coef(fit_cut)
tibble::tibble(`Threshold (ng/L)` = c(3, 5, 7, 9)) |>
  mutate(Sensitivity = plogis(b[1, 1] + b[2, 1] * `Threshold (ng/L)`),
         Specificity = 1 - plogis(b[1, 2] + b[2, 2] * `Threshold (ng/L)`)) |>
  kable(digits = 3, align = "ccc")


## -----------------------------------------------------------------------------
#| label: fig-cstat
#| fig-cap: "Nine simulated external validations of one risk model: C-statistics pooled on the logit scale (REML, Hartung-Knapp), with the 95% prediction interval. Colour: setting."
#| fig-width: 8.5
#| fig-height: 4.8
val <- tibble::tribble(~study, ~auc, ~se, ~setting,
  "Cohort A", 0.77, 0.02, "Hospital", "Cohort B", 0.70, 0.03, "Community", "Cohort C", 0.80, 0.02, "Hospital",
  "Cohort D", 0.70, 0.03, "Primary care", "Cohort E", 0.67, 0.03, "Community", "Cohort F", 0.75, 0.02, "Hospital",
  "Cohort G", 0.73, 0.03, "Community", "Cohort H", 0.68, 0.02, "Primary care", "Cohort I", 0.79, 0.02, "Hospital") |>
  mutate(yi = qlogis(auc), sei = se / (auc * (1 - auc)))      # logit C and delta-method SE
res <- rma(yi, sei = sei, data = val, method = "REML", test = "knha", slab = study)
pr  <- predict(res, transf = plogis)
set_col <- c(Hospital = "#041E42", Community = "#2E7D32", `Primary care` = "#B22222")
metafor::forest(res, transf = plogis, refline = NA, addpred = TRUE, at = seq(0.5, 0.9, 0.1),
                colout = set_col[val$setting], ilab = val$setting, ilab.lab = "Setting",
                header = c("Validation", "C [95% CI]"), mlab = "Random effects (REML, HKSJ)", xlab = "C-statistic")


## -----------------------------------------------------------------------------
#| label: cstat-bayes
#| eval: false
# library(brms)
# fit_b <- brm(yi | se(sei) ~ 1 + (1 | study), data = val, seed = 629,
#              prior = c(prior(normal(0, 2), class = Intercept), prior(student_t(3, 0, 0.5), class = sd)))
# d <- as_draws_df(fit_b)
# c_new <- plogis(d$b_Intercept + rnorm(nrow(d), 0, d$sd_study__Intercept))
# quantile(c_new, c(0.025, 0.975)); mean(c_new > 0.70)

