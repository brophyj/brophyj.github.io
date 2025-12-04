
data {
  int<lower=1> S;                      // number of studies
  array[S] int<lower=0> cross;         // number of crossovers per study
  array[S] int<lower=1> N;             // total subjects per study
}
parameters {
  real mu;                             // global average on logit scale
  real<lower=0> tau;                   // between-study heterogeneity
  vector[S] alpha;                     // random effects on logit scale
}
transformed parameters {
  vector[S] p;
  for (s in 1:S)
    p[s] = inv_logit(mu + tau * alpha[s]); // p[s] is the crossover probability for study s
}
model {
  // Priors
  mu ~ normal(0, 1.5);
  tau ~ cauchy(0, 1);
  alpha ~ normal(0, 1);
  
  // Likelihood
  for (s in 1:S)
    cross[s] ~ binomial(N[s], p[s]);
}

