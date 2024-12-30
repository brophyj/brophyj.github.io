
# Data for the model
data_list <- list(
  n2 = 2366,
  y2 = as.integer(2366 * 0.055),
  n1 = 2379,
  y1 = as.integer(2379 * 0.071)
)

# Compile and fit the model
mod <- cmdstan_model("binom_2.stan")
fit <- mod$sample(data = data_list, chains = 4, parallel_chains = 4,        refresh = 0, seed = 123)

# Extracting posterior samples
posterior_samples <- fit$draws()

# Plotting
color_scheme_set("blue")
mcmc_trace(posterior_samples, pars = c("p1", "p2", "rr"), nrow = 3)

# Plot relative risk distributions with no formatting
# mcmc_areas(posterior_samples, pars = "rr", prob = 0.95)

# print summary
fit$summary()

# Correct extraction of relative risk samples
rr_samples <- fit$draws(variables = "rr")
rr_vector <- as.vector(rr_samples)  # Convert to a simple vector for easier handling

# Calculating probabilities
prob_rr_less_09 <- mean(rr_vector < 0.9)
prob_rr_between_09_11 <- mean(rr_vector >= 0.9 & rr_vector <= 1.1)
prob_rr_greater_11 <- mean(rr_vector > 1.1)
prob_rr_less_08 <- mean(rr_vector < 0.8)

# Print the probabilities
cat("Probability RR < 0.8: ", prob_rr_less_08, 
    "Probability RR < 0.9: ", prob_rr_less_09,
    "\nProbability RR 0.9 to 1.1: ", prob_rr_between_09_11, 
    "\nProbability RR > 1.1: ", prob_rr_greater_11, "\n")
