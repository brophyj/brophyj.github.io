import math, numpy as np, pandas as pd

# ----- Table 2 (transcribed) -----
published = pd.DataFrame([
    ["RECOVERY",   0.09, 0.01, 0.67],
    ["EVOLVED",    0.79, 0.44, 1.43],
    ["AVATAR",     0.46, 0.23, 0.90],
    ["EARLY_TAVR", 0.50, 0.40, 0.63],
], columns=["study","HR","lower","upper"])

ipd = pd.DataFrame([
    ["RECOVERY",   0.05, 0.01, 0.40],
    ["EVOLVED",    0.71, 0.43, 1.18],
    ["AVATAR",     0.51, 0.27, 0.94],
    ["EARLY_TAVR", 0.57, 0.46, 0.70],
], columns=["study","HR","lower","upper"])

# ----- Convert to log-HR and SE from 95% CI -----
Z = 1.96
def add_log_effects(df):
    df = df.copy()
    df["logHR"]    = np.log(df["HR"])
    df["logLower"] = np.log(df["lower"])
    df["logUpper"] = np.log(df["upper"])
    df["SE"]       = (df["logUpper"] - df["logLower"]) / (2*Z)
    return df

pub = add_log_effects(published)
ipd2 = add_log_effects(ipd)

# ----- Bayesian RE meta-analysis with mu ~ N(0,1.5^2), tau ~ HalfNormal(0.5) -----
mu0, sigma0 = 0.0, 1.5
scale_tau   = 0.5

def mu_posterior(y, se, tau):
    V = se**2 + tau**2
    w = 1.0 / V
    S = w.sum()
    num = (mu0 / (sigma0**2)) + (w * y).sum()
    den = (1.0 / (sigma0**2)) + S
    mean = num / den
    var  = 1.0 / den
    return mean, var

def log_marginal_lik(y, se, tau):
    V = se**2 + tau**2
    w = 1.0 / V
    S = w.sum()
    sum_wy  = (w * y).sum()
    sum_wy2 = (w * (y**2)).sum()
    prec0   = 1.0 / (sigma0**2)
    den     = prec0 + S
    log_norm = 0.5 * (np.log(prec0) - np.log(den))
    quad     = -0.5 * (sum_wy2 + prec0*(mu0**2) - ((sum_wy + prec0*mu0)**2)/den)
    return log_norm + quad

def log_prior_tau(tau):
    if tau < 0: return -np.inf
    return np.log(np.sqrt(2/(np.pi*(scale_tau**2)))) - tau**2/(2*(scale_tau**2))

def sample_posterior(y, se, n_draws=10000):
    tau_grid = np.linspace(0, 2.0, 2001)
    lp = np.array([log_prior_tau(t) + log_marginal_lik(y, se, t) for t in tau_grid])
    lp -= lp.max()
    w  = np.exp(lp); w /= w.sum()
    tau_samples = np.random.choice(tau_grid, size=n_draws, p=w)
    mu_means, mu_sds = np.zeros(n_draws), np.zeros(n_draws)
    for i, t in enumerate(tau_samples):
        m, v = mu_posterior(y, se, t)
        mu_means[i] = m; mu_sds[i] = np.sqrt(v)
    mu_samples   = np.random.normal(mu_means, mu_sds)         # pooled mean, log scale
    pred_samples = np.random.normal(mu_samples, tau_samples)   # new study, log scale
    return mu_samples, pred_samples, tau_samples

pub_y, pub_se = pub["logHR"].values,  pub["SE"].values
ipd_y, ipd_se = ipd2["logHR"].values, ipd2["SE"].values

np.random.seed(1234)
mu_pub,  pred_pub,  tau_pub  = sample_posterior(pub_y,  pub_se)
mu_ipd,  pred_ipd,  tau_ipd  = sample_posterior(ipd_y,  ipd_se)

def summarize_log_to_HR(draws_log):
    draws = np.exp(draws_log)
    return np.mean(draws), np.quantile(draws, [0.025, 0.975]), np.median(draws)

# Pooled mean (Average) and Prediction, both on HR scale
avg_pub_mean,  avg_pub_ci,  avg_pub_median  = summarize_log_to_HR(mu_pub)
avg_ipd_mean,  avg_ipd_ci,  avg_ipd_median  = summarize_log_to_HR(mu_ipd)
pred_pub_mean, pred_pub_ci, pred_pub_median = summarize_log_to_HR(pred_pub)
pred_ipd_mean, pred_ipd_ci, pred_ipd_median = summarize_log_to_HR(pred_ipd)

# Posterior probabilities of benefit
p_benefit_pub  = float(np.mean(mu_pub < 0.0))              # HR < 1
p_benefit_ipd  = float(np.mean(mu_ipd < 0.0))              # HR < 1
p_benefit_pub_09 = float(np.mean(mu_pub < np.log(0.9)))    # HR < 0.9
p_benefit_ipd_09 = float(np.mean(mu_ipd < np.log(0.9)))    # HR < 0.9
