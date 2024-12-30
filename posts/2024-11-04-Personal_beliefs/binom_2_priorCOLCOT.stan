
data {
  int<lower=0> n1;  // Total number of trials in the placebo arm
  int<lower=0> y1;  // Number of successes in the placebo arm
  int<lower=0> n2;  // Total number of trials in the intervention arm
  int<lower=0> y2;  // Number of successes in the intervention arm
  real<lower=0> alpha1; // Alpha for beta prior on p1
  real<lower=0> beta1;  // Beta for beta prior on p1
  real<lower=0> alpha2; // Alpha for beta prior on p2
  real<lower=0> beta2;  // Beta for beta prior on p2
}

parameters {
  real<lower=0, upper=1> p1; // Probability of success in the placebo arm
  real<lower=0, upper=1> p2; // Probability of success in the intervention arm
}

model {
  p1 ~ beta(alpha1, beta1); // Informative prior for p1 based on COLCOT
  p2 ~ beta(alpha2, beta2); // Informative prior for p2 based on COLCOT

  y1 ~ binomial(n1, p1); // Likelihood for the placebo arm
  y2 ~ binomial(n2, p2); // Likelihood for the intervention arm
}

generated quantities {
  real rr = p2 / p1;  // Relative Risk of intervention arm over placebo arm
}
