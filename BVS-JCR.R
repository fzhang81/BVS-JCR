################################################################################
# Project: Bayesian Variable Selection via Joint Credible Regions
#
# Description:
#   This script implements a Bayesian variable selection method for computer 
#   experiments (demonstrated here on the Borehole function). It identifies 
#   active variables by constructing a Joint Credible Region (JCR) based on 
#   a generated orthogonal "inert" (dummy) variable. Variables whose posterior 
#   parameters fall outside this reference region are classified as active.
#
# Methodology:
#   1. Generates an inert variable orthogonal to the original design matrix.
#   2. Fits a Gaussian Process model with an ARD kernel using Stan (MCMC).
#   3. Constructs a 95% Joint Credible Region (ellipse) using the posterior 
#      samples of the inert variable's parameters.
#   4. Visualizes the selection boundary and classifies variables.
#
# Dependencies:
#   - R (>= 4.0.0)
#   - Libraries: rstan, MASS, car
#   - External Files:
#       - Stan model: 'L_cov_exp_quad_ARD5_samplesigma.stan'
#       - Data: Design matrix (X) and Response (y) CSV files
#
# Outputs:
#   - Console: Progress of MCMC replication and final list of selected variables.
#   - Plot: A visualization showing the elliptical selection boundary (JCR) 
#     and the position of each variable relative to the inert baseline.
#
# Date: December 2025
################################################################################

# Bayesian Variable Selection using Joint Credible Regions
library(rstan)
library(MASS)
library(car)

# Configuration
set.seed(202504510)
rep <- 100
X_path <- 'newdesign_borehole_1_200/20250331borehole2_X_n50p6_2.csv'
y_path <- 'newdesign_borehole_1_200/20250331borehole2_y_n50p6_2.csv'
stan_file <- 'L_cov_exp_quad_ARD5_samplesigma.stan'
output_name <- 'borehole2_n50_2_ARD5_noise10'

# Read and preprocess data
X <- as.matrix(read.csv(X_path, header = TRUE))
y <- as.matrix(read.csv(y_path, header = TRUE))
y_center <- scale(y, center = TRUE, scale = TRUE)
X_scale_ori <- scale(X, center = TRUE, scale = TRUE)
basis2 <- Null(X_scale_ori)

N <- nrow(X)
P <- ncol(X) + 1

# Initialize storage
median_b <- matrix(nrow = P, ncol = rep)
median_t <- matrix(nrow = P, ncol = rep)
median_sigma <- numeric(rep)
median_alpha <- numeric(rep)
all_b <- vector('list', P)
all_t <- vector('list', P)
all_sigma <- numeric(0)
all_alpha <- numeric(0)

for (i in 1:P) {
  all_b[[i]] <- numeric(0)
  all_t[[i]] <- numeric(0)
}

# Main simulation loop
for (i in 1:rep) {
  cat("Running replication", i, "of", rep, "\n")
  
  # Generate inert variable
  c <- matrix(runif(N - P + 1, min = 0, max = 1), ncol = 1)
  x_inert <- basis2 %*% c
  X_scale <- scale(cbind(X_scale_ori, x_inert), center = TRUE, scale = TRUE)
  
  # Fit Stan model
  data_list <- list(N = N, D = P, x = X_scale, y = as.vector(y_center))
  fit <- stan(file = stan_file, data = data_list, iter = 2000, chains = 1, 
              control = list(adapt_delta = 0.95))
  
  # Extract samples (theta, sigma, alpha, beta)
  samples <- fit@sim$samples[[1]]
  burn_in <- 1001:2000
  
  # Extract parameters
  for (j in 1:P) {
    theta_samples <- exp(-as.numeric(samples[[j]])[burn_in])
    median_t[j, i] <- median(theta_samples)
    all_t[[j]] <- c(all_t[[j]], theta_samples)
  }
  
  sigma_samples <- as.numeric(samples[[P + 1]])[burn_in]
  alpha_samples <- as.numeric(samples[[P + 2]])[burn_in]
  median_sigma[i] <- median(sigma_samples)
  median_alpha[i] <- median(alpha_samples)
  all_sigma <- c(all_sigma, sigma_samples)
  all_alpha <- c(all_alpha, alpha_samples)
  
  for (j in 1:P) {
    beta_samples <- as.numeric(samples[[P + 2 + j]])[burn_in]
    median_b[j, i] <- median(beta_samples)
    all_b[[j]] <- c(all_b[[j]], beta_samples)
  }
}

# Variable selection function
is_point_in_region <- function(point, ell_l, x_min, x_max) {
  x <- point[1]
  y <- point[2]
  if (x < x_min || x > x_max) return(1)
  ell_y <- approx(ell_l[, 1], ell_l[, 2], xout = x)$y
  if (y < ell_y) return(1)
  return(0)
}

# Create selection plot and analysis
inert_idx <- P
x_data <- as.numeric(median_b[inert_idx, ])
y_data <- as.numeric(exp(-(-log(median_t[inert_idx, ]))^2))

# Create ellipse
ell <- dataEllipse(x_data, y_data, 
                   main = paste('Variable Selection:', output_name),
                   levels = 0.95,
                   xlim = range(x_data) + c(-0.02, 0.02),
                   ylim = range(y_data) + c(-0.02, 0.02),
                   xlab = paste('Median of beta', inert_idx),
                   ylab = paste('Median of exp(-theta', inert_idx, ')'))

# Process ellipse for region definition  
x_min_index <- which(ell[, 1] == min(ell[, 1]))
x_max_index <- which(ell[, 1] == max(ell[, 1]))

if (length(x_min_index) == 1 && length(x_max_index) == 1) {
  ell <- rbind(ell[(x_max_index + 1):nrow(ell), ], ell[1:x_max_index, ])
  ell_l <- ell[1:(nrow(ell) / 2), ]
  ell_u <- ell[(nrow(ell) / 2 + 1):nrow(ell), ]
} else {
  ell_l <- ell[(nrow(ell) / 2 + 1):nrow(ell), ]
  ell_u <- ell[1:(nrow(ell) / 2), ]
}

x_min <- min(ell[, 1])
x_max <- max(ell[, 1])

# Create detailed plot
plot(ell_l[, 1], ell_l[, 2], type = 'n',
     xlim = range(x_data) + c(-0.005, 0.005),
     ylim = range(y_data) + c(-0.001, 0.001),
     xlab = paste('Median of beta', inert_idx),
     ylab = paste('Median of exp(-theta', inert_idx, ')'),
     main = paste('Variable Selection:', output_name))

# Fill selection region
x_fill <- if (ell_l[1, 1] > ell_l[nrow(ell_l), 1]) {
  c(x_max, ell_l[, 1], x_min, x_min, x_max)
} else {
  c(x_min, ell_l[, 1], x_max, x_max, x_min)
}
y_fill <- c(1, ell_l[, 2], ell_l[which(ell_l[, 1] == max(ell_l[, 1])), 2], 1, 1)

polygon(x_fill, y_fill, col = rgb(0.5, 0.5, 1, alpha = 0.3), border = NA)
lines(ell_l[, 1], ell_l[, 2], col = 'blue', lwd = 2)
lines(ell_u[, 1], ell_u[, 2], col = 'blue', lty = 2, lwd = 2)
points(x_data, y_data, pch = 20, col = 'blue')

# Add boundary lines
segments(x_min, ell_l[which(ell_l[, 1] == min(ell_l[, 1])), 2], x_min, 1, 
         col = 'blue', lwd = 2)
segments(x_max, ell_l[which(ell_l[, 1] == max(ell_l[, 1])), 2], x_max, 1, 
         col = 'blue', lwd = 2)
abline(h = 1, col = 'blue', lty = 2)

# Add labels for original variables
for (j in 1:(P - 1)) {
  text(x = median(all_b[[j]]), 
       y = median(exp(-(-log(all_t[[j]]))^2)),
       labels = paste0('x', j), cex = 1.5, pos = 4, col = 'red')
}

# Perform variable selection and display results
selection_results <- matrix(0, nrow = 1, ncol = P - 1)
for (j in 1:(P - 1)) {
  point <- c(median(all_b[[j]]), median(exp(-(-log(all_t[[j]]))^2)))
  selection_results[1, j] <- is_point_in_region(point, ell_l, x_min, x_max)
}

# Display results
cat("\nVariable Selection Results:\n")
for (i in 1:(P - 1)) {
  status <- ifelse(selection_results[1, i] == 0, "Selected", "Not Selected")
  cat("Variable", i, ":", status, "\n")
}
