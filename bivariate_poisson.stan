functions {
  real bivariate_poisson_lpmf(array[] int y, vector lambda) {
    int m = y[1];
    int n = y[2];
    real lp = -sum(lambda);  // log(e^{-sum lambda})
    real lsum = negative_infinity();
    
    for (k in 0:min(m, n)) {
      real term = (
        (m - k) * log(lambda[1])
        + (n - k) * log(lambda[2])
        + k * log(lambda[3])
        - lgamma(m - k + 1)
        - lgamma(n - k + 1)
        - lgamma(k + 1)
      );
      lsum = log_sum_exp(lsum, term);
    }
    
    return lp + lsum;
  }
}
data {
  int<lower=1> N;
  array[N] int<lower=0> measles;
  array[N] int<lower=0> rubella;
  vector[N] log_total_population;
}
parameters {
  real alpha1;
  real alpha2;
  real alpha3;
  
  real beta1;
  real beta2;
  real beta3;
}
model {
  for (n in 1:N) {
    vector[3] lambda;
    lambda[1] = exp(alpha1 + beta1 * log_total_population[n]);
    lambda[2] = exp(alpha2 + beta2 * log_total_population[n]);
    lambda[3] = exp(alpha3 + beta3 * log_total_population[n]);
    target += bivariate_poisson_lpmf({measles[n], rubella[n]} | lambda);
  }

  // priors
  alpha1 ~ normal(0, 5);
  alpha2 ~ normal(0, 5);
  alpha3 ~ normal(0, 5);
  beta1 ~ normal(0, 2);
  beta2 ~ normal(0, 2);
  beta3 ~ normal(0, 2);
}
