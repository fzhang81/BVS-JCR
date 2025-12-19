# Bayesian Variable Selection via Joint Credible Regions

This repository contains the R implementation for a Bayesian variable selection method applied to computer experiments (specifically demonstrated on the Borehole function).

The algorithm identifies active variables by constructing a **Joint Credible Region (JCR)** based on a generated "inert" (dummy) variable. Variables whose posterior parameters fall outside this reference region are classified as active.

## Key Features

* **Inert Variable Construction:** Generates an orthogonal inert variable to serve as a baseline for noise.
* **Bayesian Inference:** Uses `rstan` to sample from the posterior distribution of the model parameters (Gaussian Process with ARD kernel).
* **Selection Mechanism:** Constructs a 95% Joint Credible Region (ellipse) using the posterior samples of the inert variable's parameters ($\beta$ and transformed $\theta$).
* **Visualization:** Visually maps the selection boundary and variable positions.

## Dependencies

* **R** (>= 4.0.0)
* **R Packages:**
  * `rstan` (for MCMC sampling)
  * `MASS` (for statistical functions)
  * `car` (for data ellipses)

## Usage

1. Ensure your data files (`X.csv`, `y.csv`) and the Stan model file (`L_cov_exp_quad_ARD5_samplesigma.stan`) are in the working directory.
2. Run the main R script. The code will:
  * Perform MCMC sampling.
  * Compute the median $\beta$ and length-scale parameters.
  * Generate a selection plot.
  * Output the selection status (Selected/Not Selected) for each variable.
