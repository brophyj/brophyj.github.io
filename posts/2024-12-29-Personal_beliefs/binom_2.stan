
data {
  int<lower=0> n1;  // Total number of trials in the placebo arm
  int<lower=0> y1;  // Number of successes in the placebo arm
  int<lower=0> n2;  // Total number of trials in the intervention arm
  int<lower=0> y2;  // Number of successes in the intervention arm
}

parameters {
  real<lower=0, upper=1> p1; // Probability of success in the placebo arm
  real<lower=0, upper=1> p2; // Probability of success in the intervention arm
}

model {
  p1 ~ beta(1, 1); // Default prior for p1
  p2 ~ beta(1, 1); // Default prior for p2

  y1 ~ binomial(n1, p1); // Likelihood for the placebo arm
  y2 ~ binomial(n2, p2); // Likelihood for the intervention arm
}

generated quantities {
  real rr = p2 / p1;  // Relative Risk of intervention arm over placebo arm
}

