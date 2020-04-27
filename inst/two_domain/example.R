# devtools::install_github("jatotterdell/varapproxr")

library(ascotsims)
library(parallel)
library(tidyverse)

cores <- parallel::detectCores() - 2

# Design matrix for domain
D <- expand.grid(HCQ = c(0, 1), LR = c(0, 1), CP = c(0, 1))
X <- model.matrix( ~ HCQ * LR + CP, data = D)
colnames(X)[1] <- "SoC"
rownames(X) <- sapply(1:nrow(X), function(i) paste(colnames(X)[1:4][X[i, 1:4] == 1], collapse = " + "))
XA <- X[, c(1,2,3,5)]
XB <- X[, 4, drop = F]
M0 <- rep(0, ncol(X))
S0 <- diag(2, ncol(X))
X <- cbind(XA, XB)

# Data generation parameters
p <- rep(0.25, nrow(X))
b <- solve(crossprod(X)) %*% t(X) %*% qlogis(p)

# Generate data
full_dat <- mclapply(1:100, function(i) gen_potential_outcomes(2000, p), mc.cores = cores)

# Run trial on data
res <- mclapply(full_dat, function(i) lr_brar_trial2(seq(500, 2000, 500), i, XA, XB, M0, S0), mc.cores = cores)

# Summarise trial
rest <- tibble_trial_data(res, mc.cores = cores)
resp <- tibble_model_data(res, mc.cores = cores)
