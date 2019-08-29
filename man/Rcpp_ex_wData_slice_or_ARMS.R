## Load R packages; make sure to inclde Rcpp packages
require(Rcpp)
require(RcppArmadillo)
require(RcppNumerical)
library(MASS)
library(HI)
library(tidyverse)

### Write or read-in C++ code
sourceCpp("ARMS_slice_LOD_test.cpp", verbose = TRUE)

### Univariate normal distribution
## Generate data point from slice under LOD rate of 10%
# Like function
sigma_2_y = 1
beta = c(1, 2, 2)
mean_x_preds = c(0,0)
sigma_x_preds = matrix(c(1, 0, 0, 1), byrow = TRUE, nrow=2)

x_data <- c(1, 0, NA)
y_data <- 0
LOD_rate = 0.1
set.seed(2)

sample <- slice_sample_normal_ex(y_data, x_data, mean_x_preds, beta, sigma_2_y, sigma_x_preds,
                                 LOD_lower = -10, LOD_upper = qnorm(LOD_rate), sample_size = 1000000)

hist(sample, breaks=100)

### Multivariate regression, normal distribution
# Example 1: No LOD
n <- 100
x_data <- cbind(rep(1, n),mvrnorm(n,c(0,0),matrix(c(1,0,0,1),2)))
y_data <- x_data%*%t(t(c(1,1,1)))+rnorm(n,0,1)

LOD_rate=NA; LOD_mat = rbind(c(NA, NA),
                  c(NA, NA),
                  c(-100, qnorm(LOD_rate)),
                  c(-100, qnorm(LOD_rate)))

x_data_obs <- x_data
data_obs <- data.frame(y_data, x_data_obs)
fit <- lm(y_data~X2+X3, data=data_obs)

set.seed(1)

test <- LOD_fit_multiple(y_data=y_data, x_data=x_data_obs, 
                         mean_x_preds=colMeans(x_data_obs[,-1], na.rm = TRUE), 
                         beta=fit$coefficients, 
                         sigma_2_y = sigma(fit)^2, 
                         sigma_x_preds = cov(x_data_obs[,-1], use = "complete.obs"),
                         no_of_samples=250, threshold = 0.001, max_iterations = 100,
                         LOD_u_l = LOD_mat,
                         sampler = 0)

fit$coefficients # Same as Betas from LOD fn

test_bt <- bootstrap_multi_test(num_of_boots = 25,
                                y_data=y_data, x_data=x_data_obs,
                                no_of_samples=250, threshold = 0.001, max_iterations = 100,
                                LOD_u_l = LOD_mat,
                                sampler = 0)

## Example 2: 10% LOD
n <- 100
x_data <- cbind(rep(1, n),mvrnorm(n,c(0,0),matrix(c(1,0,0,1),2)))
y_data <- x_data%*%t(t(c(1,1,1)))+rnorm(n,0,1)

LOD_rate <- 0.1; LOD_mat = rbind(c(NA, NA),
                                 c(NA, NA),
                                 c(-100, qnorm(LOD_rate)),
                                 c(-100, qnorm(LOD_rate)))
x_data_obs <- x_data
x_data_obs[x_data_obs[,3] < qnorm(LOD_rate),3] = NA
data_obs <- data.frame(y_data, x_data_obs)
fit <- lm(y_data~X2+X3, data=data_obs)

set.seed(1)

test <- LOD_fit_multiple(y_data=y_data, x_data=x_data_obs, 
                         mean_x_preds=colMeans(x_data_obs[,-1], na.rm = TRUE), 
                         beta=fit$coefficients, 
                         sigma_2_y = sigma(fit)^2, 
                         sigma_x_preds = cov(x_data_obs[,-1], use = "complete.obs"),
                         no_of_samples=250, threshold = 0.001, max_iterations = 100,
                         LOD_u_l = LOD_mat,
                         sampler = 0)

test_bt <- bootstrap_multi_test(num_of_boots = 25,
                                y_data=y_data, x_data=x_data_obs,
                                no_of_samples=250, threshold = 0.001, max_iterations = 100,
                                LOD_u_l = LOD_mat,
                                sampler = 0)

## Example 3: 50% LOD
n <- 100
x_data <- cbind(rep(1, n),mvrnorm(n,c(0,0,0),matrix(c(1,0,0,0,1,0,0,0,1),3)))
y_data <- x_data%*%t(t(c(1,1,1,1)))+rnorm(n,0,1)

LOD_rate <- 0.5; LOD_mat = rbind(c(NA, NA),
                                 c(NA, NA),
                                 c(-100, qnorm(LOD_rate)),
                                 c(-100, qnorm(LOD_rate)))
x_data_obs <- x_data
x_data_obs[x_data_obs[,4] < qnorm(LOD_rate),4] = NA
data_obs <- data.frame(y_data, x_data_obs)
fit <- lm(y_data~X2+X3+X4, data=data_obs)

set.seed(1)

test <- LOD_fit_multiple(y_data=y_data, x_data=x_data_obs, 
                         mean_x_preds=colMeans(x_data_obs[,-1], na.rm = TRUE), 
                         beta=fit$coefficients, 
                         sigma_2_y = sigma(fit)^2, 
                         sigma_x_preds = cov(x_data_obs[,-1], use = "complete.obs"),
                         no_of_samples=250, threshold = 0.001, max_iterations = 100,
                         LOD_u_l = LOD_mat,
                         sampler = 0)

start <- Sys.time()
test_bt <- bootstrap_multi_test(num_of_boots = 25,
                                y_data=y_data, x_data=x_data_obs,
                                no_of_samples=250, threshold = 0.001, max_iterations = 100,
                                LOD_u_l = LOD_mat,
                                sampler = 0)
end <- Sys.time() - start

### Multivariate regression, normal distribution, many LOD covariates
## Example 1: 50% LOD
n <- 100
x_data <- cbind(rep(1, n),mvrnorm(n,c(0,0,0),matrix(c(1,0,0,0,1,0,0,0,1),3)))
y_data <- x_data%*%t(t(c(1,1,1,1)))+rnorm(n,0,1)

LOD_rate <- 0.50; LOD_mat = rbind(c(NA, NA),
                                 c(NA, NA),
                                 c(-100, qnorm(LOD_rate)),
                                 c(-100, qnorm(LOD_rate)))
x_data_obs <- x_data
x_data_obs[x_data_obs[,dim(x_data)[2]-1] < qnorm(LOD_rate),dim(x_data)[2]-1] = NA
x_data_obs[x_data_obs[,dim(x_data)[2]] < qnorm(LOD_rate),dim(x_data)[2]] = NA

data_obs <- data.frame(y_data, x_data_obs)
fit <- lm(y_data~X2+X3+X4, data=data_obs)

test <- LOD_fit_multiple(y_data=y_data, x_data=x_data_obs, 
                         mean_x_preds=colMeans(x_data_obs[,-1], na.rm = TRUE), 
                         beta=fit$coefficients, 
                         sigma_2_y = sigma(fit)^2, 
                         sigma_x_preds = cov(x_data_obs[,-1], use = "complete.obs"),
                         no_of_samples=250, threshold = 0.001, max_iterations = 100,
                         LOD_u_l = LOD_mat,
                         sampler = 0)

start <- Sys.time()
test_bt <- bootstrap_multi_test(num_of_boots = 25,
                                y_data=y_data, x_data=x_data_obs,
                                no_of_samples=250, threshold = 0.001, max_iterations = 100,
                                LOD_u_l = LOD_mat,
                                sampler = 0)
end <- Sys.time() - start

boot_results <- do.call("rbind", test_bt)

summary(fit)$coefficients[,"Std. Error"]
apply(boot_results, 2, sd)