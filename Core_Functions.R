library(tidyverse)
library(qad)

# density_to_kernel creates a linear interpolation of the ecdf based on the checkerboard copula density estimator from qad:ECBC()
# Input:
# col = A column of a density matrix created by qad::EBCB()

# Output:
# A vector of numerical values (k_A(x, [0, y_1]), ..., k_A(x, [0, y_N])) that correspond to the ecdf at the centerpoints of the checkerboard fields.

density_to_kernel <- function(col){
  col_2 <- col / 2
  col <- c(0, col)
  col <- cumsum(col)
  total <- col[length(col)]
  col <- col[1:(length(col) - 1)]
  col <- col + col_2
  col <- col / total
  return(col)
}

# pairwise_diff_phi computes the pairwise differences between ecdfs at one interpolation point

# Input:
# col = Vector of values of the ecdfs at given interpolation point

# Output:
# Numerical vector of pariwise differences phi(k_A(x_1, [0, y]) - k_A(x_2, [0, y])) between ecdfs at given interpolation point
pairwise_diff_phi <- function(col, phi){
  temp_grid <- expand.grid(x1 = col, x2 = col)
  return(phi(temp_grid[, 1] - temp_grid[, 2]))
}


# Lambda_phi estimates the functional Lambda_phi

# Input:
# ecbc = A matrix that is the output of qad::ECBC(),
# phi = The function phi to be used

# Output:
# The estimate of Lambda_phi (numeric)
Lambda_phi <- function(ecbc, phi){
  N <- ncol(ecbc)
  alpha <- 2/3 * phi(0) + 1/6 * phi(1) + 1/6 * phi(-1)
  # Compute k_A(x_i, [0, y_j]) for all grid points of the checkerboard copula
  kernel_mat <- t(apply(ecbc, 1, FUN = density_to_kernel))
  # Compute all the N * (N - 1) * N / 2 L^p distances: |k_A(x_i, [0, y_j]) - k_A(x_i', [0, y_j])|^p
  lambda_temp <- map(1:N, .f = function(i) pairwise_diff_phi(col = kernel_mat[, i], phi = phi))
  # Sum them all up for the finite sample version of the integral
  lambda_temp <- sum(unlist(lambda_temp))
  # Divide by the number of elements that are summed up
  lambda_temp <- lambda_temp / (N * (N + 1) * (N - 1))
  lambda_temp <- lambda_temp / alpha
  return(lambda_temp)
}