functions {
  matrix L_cov_exp_quad_ARD(vector[] x,
                            real alpha,
                            vector theta,
                            real delta,
                            real sigma) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    real sq_sigma = square(sigma);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + sq_sigma + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-1 * dot_self((x[i] - x[j]) .* theta));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + sq_sigma + delta;
    return cholesky_decompose(K);
  }
}
data {
  int<lower=1> N;
  int<lower=1> D;
  vector[D] x[N];//array[N] vector[D] x;
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  vector<lower=0>[D] theta;
  real<lower=0> sigma;
  real<lower=0> alpha;
  vector[D] beta;
}
model {
  vector[N] mu;  
  matrix[N, N] L_K = L_cov_exp_quad_ARD(x, alpha, theta, delta, sigma);
  for (n in 1:N) {
  mu[n] = dot_product(x[n], beta);
  }

  theta ~ exponential(1);
  sigma ~ inv_gamma(50,5);
  alpha ~ inv_gamma(5,5);
  beta ~ normal(0,1);

  y ~ multi_normal_cholesky(mu, L_K);
}
