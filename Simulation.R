library(tidyverse)
library(qad)
library(copula)
library(cubature)
source("./Core_Functions.R")

# mo_kernel returns the kernel k_A(x, [0, y]) of a Marshall-Olkin copula with parameters alpha and beta for the interval [0, y] at point x

# Input:
# alpha = First Parameter
# beta = Second Parameter
# x = x-value
# y = y-value

# Output:
# k_A(x, [0, y]) (numeric)
mo_kernel <- function(alpha, beta, x, y){
  if(x == 0 | y == 0){
    return(0)
  }
  if(beta * log(y) <= alpha * log(x)){
    return((1 - alpha) * x^(-alpha) * y)
  }
  if(beta * log(y) > alpha * log(x)){
    return(y^(1 - beta))
  }
}

# kernel_diff takes the difference phi(k_A(x_1, [0,y]) - k_A(x_2, [0, y])) between the conditional cdfs (i.e. kernels) of a MO-Copula

# Input:
# x = vector with components (x_1, x_2, y), i.e. the two-values to condition on and the upper bound of interval over which the mass is calculated
# alpha = First Parameter of MO Copula
# beta = Second Parameter of MO Copula
# phi = function phi to be applied to the difference of kernels
# const = normalizng constant alpha. Usually (phi(-1) + phi(1)) / 6

# Output:
# phi(k_A(x[1], [0, x[3]]) - k_A(x[2], [0, x[3]])) (numeric)
kernel_diff <- function(x, alpha, beta, phi, const){
  return((1/const) * phi(mo_kernel(alpha = alpha, beta = beta, x = x[1], y = x[3]) - mo_kernel(alpha = alpha, beta = beta, x = x[2], y = x[3])))
}


# mo_int_fun calculates the true value of Lambda_phi for a MO-Copula with parameters alpha and beta and function phi

# Input:
# alpha = First Parameter of MO Copula
# beta = Second Parameter of MO Copula
# phi_char = function phi to be applied to the difference of kernels, given as character

# Output
# Lambda_phi (numeric)
mo_int_fun <- function(alpha, beta, phi_char, ...){
  phi <- eval(parse(text = phi_char))
  const <- 2/3 * phi(0) + 1/6 * phi(1) + 1/6 * phi(-1)
  int_object <- pcubature(f = kernel_diff, lowerLimit = c(0, 0, 0), upperLimit = c(1, 1, 1), tol = 0.01, alpha = alpha, beta = beta, phi = phi, const = const)
  limit <- int_object$integral
  error <- int_object$error
  fun_eval <- int_object$functionEvaluations
  return_code <- int_object$returnCode
  output_frame <- data.frame(alpha = alpha, beta = beta, phi = phi_char, limit = limit, error = error, fun_eval = fun_eval, return_code = return_code, ...)
  print(paste0("alpha = ", alpha, ", beta = ", beta, ", phi = ", phi_char))
  return(output_frame)
}

# mo_sim_fun estimates Lambda_phi for a MO-Copula in the specified setting

# Input:
# alpha = First Parameter of MO Copula
# beta = Second Parameter of MO Copula
# phi_char = function phi to be applied to the difference of kernels, given as character
# n = Sample Size
# run = Seed to be used

# Output:
# A data frame with one row containing the values of 
# alpha = First Parameter of MO Copula
# beta = Second Parameter of MO Copula
# phi_char = function phi to be applied to the difference of kernels, given as character
# n = Sample Size
# run = Seed to be used
# lambda = estimate of Lambda_phi

mo_sim_fun <- function(alpha, beta, phi_char, n, run, ...){
  phi <- eval(parse(text = phi_char))
  temp_copula <- moCopula(param = c(alpha, beta))
  set.seed(run)
  temp_sample <- rCopula(n = n, copula = temp_copula)
  ecbc <- ECBC(X = temp_sample[,1], Y = temp_sample[,2])
  lambda <- Lambda_phi(ecbc, phi = phi)
  output_frame <- data.frame(alpha = alpha, beta = beta, n = n, run = run, phi = phi_char, lambda = lambda, ...)
  print(paste0("alpha = ", alpha, ", beta = ", beta, ", n = ", n, ", phi = ", phi_char, ", run = ", run))
  return(output_frame)
}

# Creating a data.frame sim_frame containing all combinations of simulation parameters and a data.frame
# limit_frame which will contain the true values of Lambda_phi
par_frame <- data.frame(alpha = c(1, 1, 0.2, 0.3, 1, 0.5), beta = c(0, 1, 0.7, 1, 0.7, 0.5))

setup_frame <- expand.grid(run = 1:1000, n = c(10, 50, 100, 500, 1000, 5000, 10000))

phi_frame <- data.frame(phi_char = c("function(x) abs(x)", "function(x) abs(x)^2", "function(x) abs(x)^3",
                                     "function(x) exp(1/5 * x) - 1", "function(x) exp(x) - 1", "function(x) exp(5 * x) - 1",
                                     "function(x) exp(abs(1/5 * x)) - 1", "function(x) exp(abs(x)) - 1", "function(x) exp(abs(5 * x)) - 1"),
                        type = rep(c("abs(x)^p", "exp(c * x) - 1", "exp(abs(c * x)) - 1"), each = 3),
                        p_c = c(c("1", "2", "3"), rep(c("1/5", "1", "5"), times = 2)))

sim_frame <- crossing(par_frame, setup_frame, phi_frame)

sim_frame <- crossing(par_frame, setup_frame, phi_frame)

limit_frame <- crossing(par_frame, phi_frame)

# In special cases the true value is known
limit_num_frame <- limit_frame %>%
  filter(!(alpha == 1 & beta == 0)) %>%
  filter(!(alpha == 1 & beta == 1))

limit_man_frame <- limit_frame %>%
  filter((alpha == 1 & beta == 0) | (alpha == 1 & beta == 1)) %>%
  mutate(limit = case_when(alpha == 1 & beta == 0 ~ 0,
                           alpha == 1 & beta == 1 ~ 1)) %>%
  rename(phi = phi_char)

# In the other cases it is calculated by numerical integration
limits <- pmap_dfr(limit_num_frame, mo_int_fun)

limits <- bind_rows(limits, limit_man_frame)

# True Values are saved
save(limits, file = "./Simulation_Results/True_Values.Rda")

# Estimates for each simulation scenario are created and saved
results <- pmap_dfr(.l = sim_frame, .f = mo_sim_fun)
save(results, file = "./Simulation_Results/Simulation_MO_1000.Rda")