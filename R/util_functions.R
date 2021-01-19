#' Basic pseudoinverse
#'
#' @param X A matrix
#' @return Pseudoinverse of X
#' @export
psolve <- function(X) {
  solve(crossprod(X), t(X))
}

#' @export
findfirst <- function(x, v = NA) {
  j <- which(x)
  if(length(j)) min(j) else v
}

#' Multiply every column combination between two matrices
#'
#' Function is used for generating interactions between domain matrices.
#' Note that terms on matrix A change first, i.e. `out[, 1] = A[, 1]B[, 1]`, `out[, 2] = A[, 2]B[,1]` etc.
#'
#' @param A matrix 1
#' @param B matrix 2
#' @param sep Separation in column names.
#' @return Matrix of interactions between A and B
colwise_mult <- function(A, B, sep = "") {
  out <- matrix(apply(A, 2, FUN = function(x) x * B), nrow(A), ncol(A) * ncol(B))
  rownames(out) <- rownames(A)
  colnames(out) <- paste(colnames(A), rep(colnames(B), each = ncol(A)), sep = sep)
  return(out)
}

#' Mass-weighted urn randomisation
#'
#' @param target_alloc The target allocation ratios
#' @param sample_size The number of allocations to generate
#' @param alpha Parameter to control imbalance between arms
#' @return A list detailing the mass-weighted-urn process.
#' @export
mass_weighted_urn_design <- function(
  target_alloc,
  sample_size,
  alpha = 4
) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)

  imbalance_cap <- sqrt(sum(((alpha - 1)*(1 - prob_alloc) + (arms - 1))^2))

  for(i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    p[i - 1, ] <- pmax(alpha * prob_alloc - n[i - 1, ] + (i - 1)*prob_alloc, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i-1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i-1]] <- n[i, trt[i-1]] + 1
    # Update urn masses
    x[i, trt[i-1]] <- x[i - 1, trt[i-1]] - 1 + prob_alloc[trt[i-1]]
    x[i, -trt[i-1]] <- x[i - 1, -trt[i-1]] + prob_alloc[-trt[i-1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1)*prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p))
}

#' Variance of a beta random variable with parameters a and b
#'
#' @param a Par 1
#' @param b Par 2
#' @examples
#' beta_var(5, 5)
#' @export
beta_var <- function(a, b) {
  return( a * b / ( (a + b)^2 * (a + b + 1)) )
}


#' Variance of a added beta random variables
#'
#' Always taken in reference to control, e.g. `(a[1], b[1])`.
#'
#' @param a Collection of shape 1 par
#' @param b Collection of shape 2 par
#' @examples
#' diff_beta_var(c(5, 5, 4), c(5, 5, 4))
#' @export
diff_beta_var <- function(a, b) {
  return(beta_var(a[1], b[1]) + beta_var(a[-1], b[-1]))
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller if reverse) to another Beta RV by integration.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @param ... Other arguments to `pracma::quadgk`
#' @examples
#' beta_prob_supr_numeric(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr_numeric(c(5, 10, 15), c(15, 10, 5), reverse = TRUE)
#' @importFrom pracma quadgk
#' @importFrom stats dbeta pbeta
#' @export
beta_prob_supr_numeric <- function(a, b, reverse = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k-1)
  for(i in 2:k) {
    f <- function(z)  dbeta(z, a[1], b[1])*(1 - pbeta(z, a[i], b[i]))
    ans[i-1] <- tryCatch(
      pracma::quadgk(f, 0, 1, ...),
      error = function(err) NA)
  }
  if(reverse) ans <- 1 - ans
  return(ans)
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller if reverse) to another Beta RV by normal approximation.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @examples
#' beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5), reverse = TRUE)
#' @importFrom stats pnorm
#' @export
beta_prob_supr_approx <- function(a, b, reverse = FALSE) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  m <- a / (a + b)
  v <- a*b / ( (a + b)^2 * (a + b + 1))
  z <- (m[-1] - m[1]) / sqrt(v[-1] + v[1])
  ans <- pnorm(z)
  if(reverse) ans <- 1 - ans
  return(ans)
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller) to another Beta RV by normal approximation.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @param approx Use normal approximation instead of numerical integration
#' @param ... Other arguments to `pracma::quadgk`
#' @examples
#' beta_prob_supr(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr(c(5, 10, 15), c(15, 10, 5), approx = TRUE)
#' @export
beta_prob_supr <- function(a, b, reverse = FALSE, approx = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k-1)
  if(!approx) {
    return(beta_prob_supr_numeric(a, b, reverse, ...))
  } else {
    return(beta_prob_supr_approx(a, b, reverse))
  }
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param ... Other arguments to `pracma::quadgk`
#' @return A vector of `length(a)` giving Pr(a[i] is max)
#' @examples
#' beta_prob_best_numeric(1 + seq(1, 9, 2), 1 + rep(10, 5))
#' @export
beta_prob_best_numeric <- function(a, b, minimum = FALSE, ...) {
  k <- length(a)
  ans <- numeric(k)
  for(i in 1:k) {
    if(!minimum) {
      f <- function(z)
        sapply(z, function(x) dbeta(x, a[i], b[i]) * prod(pbeta(x, a[-i], b[-i])))
    } else {
      f <- function(z)
        sapply(z, function(x) dbeta(x, a[i], b[i]) * prod(1 - pbeta(x, a[-i], b[-i])))
    }
    ans[i] =   tryCatch(
      pracma::quadgk(f, 0, 1, ...),
      error = function(err) NA)
  }
  return(ans)
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param nsim Number of `rbeta` draws to use in assessment.
#' @return A vector of `length(a)` giving Pr(a[i] is max)
#' @examples
#' beta_prob_best_approx(1 + seq(1, 9, 2), 1 + rep(10, 5))
#' @importFrom stats rbeta
#' @export
beta_prob_best_approx <- function(a, b, minimum = FALSE, nsim = 2e4) {
  k <- length(a)
  ans <- numeric(k)
  r <- matrix(rbeta(k*nsim, a, b), nsim, k, byrow = T)
  if(!minimum) {
    ans <- prop.table(table(factor(max.col(r), levels = 1:k, labels = 2:(k+1))))
  } else {
    ans <- prop.table(table(factor(max.col(1 - r), levels = 1:k, labels = 2:(k+1))))
  }
  return(ans)
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param approx Use approximation instead of integration
#' @param ... Other arguments to `pracma::quadgk`
#' @return A vector of `length(a)` giving Pr(a[i] is max)
beta_prob_best <- function(a, b, minimum = FALSE, approx = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k)
  if(!approx) {
    return(beta_prob_best_numeric(a, b, minimum, ...))
  } else {
    return(beta_prob_best_approx(a, b, minimum))
  }
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param ... Other arguments to `quadgk`.
#' @importFrom pracma quadgk
ent_best_arm_numeric <- function(a, b, minimum = FALSE, ...) {
  k <- length(a)
  if(!minimum) {
    f <- function(x) {
      ss <- sum(sapply(1:k, function(i) {
        dbeta(x, a[i], b[i]) * prod(pbeta(x, a[-i], b[-i]))
      }))
      ss <- ss * log(ss)
      if(is.na(ss))
        ss <- 0
      return(-ss)
    }
  } else {
    f <- function(x) {
      ss <- sum(sapply(1:k, function(i) {
        dbeta(x, a[i], b[i]) * prod(1 - pbeta(x, a[-i], b[-i]))
      }))
      ss <- ss * log(ss)
      if(is.na(ss))
        ss <- 0
      return(-ss)
    }
  }
  ans <- pracma::quadgk(Vectorize(f), 0, 1, ...)
  return(ans)
}


#' Evaluate joint density of maximum of Beta RVs
#'
#' @param x Value to evaluate density at
#' @param a Vector of shape1 parameters
#' @param b Vector of shape1=2 parameters
#' @param minimum Use minimum instead of maximum.
beta_best_dens <- function(x, a, b, minimum = FALSE) {
  k <- length(a)
  if(!minimum) {
    sapply(x, function(z) sum(sapply(1:k, function(i) dbeta(z, a[i], b[i]) * prod(pbeta(z, a[-i], b[-i])))))
  } else {
    sapply(x, function(z) sum(sapply(1:k, function(i) dbeta(z, a[i], b[i]) * prod(1 - pbeta(z, a[-i], b[-i])))))
  }
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param nn Grid points.
#' @importFrom pracma quadgk
ent_best_arm_approx <- function(a, b, minimum = FALSE, nn = 100) {
  k <- length(a)
  del <- 1/nn
  x <- seq(0 + del, 1 - del, length.out = nn - 1)
  f <- beta_best_dens(x, a, b, minimum)
  f[f == 0] <- .Machine$double.eps^12
  return(-sum(f*log(f) * del))
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param approx Use approximation instead of numerical.
#' @param ... Other arguments to `quadgk`.
ent_best_arm <- function(a, b, minimum = F, approx = F, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  if(!approx) {
    return(ent_best_arm_numeric(a, b, minimum, ...))
  } else {
    return(ent_best_arm_approx(a, b, minimum))
  }
}


#' Calculate change in utility using variance
#'
#' @param n Sample size in each arm
#' @param y Response in each arm
#' @param a Prior shape
#' @param b Prior shape
utility_var <- function(n, y, a = 1, b = 1) {
  K <- length(n)
  # Prior variance
  priorV <- diff_beta_var(a, b)
  priorU <- sum(priorV)
  # Current utility
  V <- diff_beta_var(y + a, n - y + b)
  U <- sum(V)
  # Possible allocations
  new_n <- t(n + diag(1, K))
  # If response to observed from allocation
  new_y <- t(y + diag(1, K))
  # For each arm, what is the expected gain in utility
  Eu <- numeric(K)
  for(j in 1:K) {
    # predicted probability of event
    m <- (a + y[j]) / (a + b + n[j])
    # utility if event occurs
    H1 <- diff_beta_var(new_y[j, ] + a, new_n[j, ] - new_y[j, ] + b)
    # utility if no event occurs
    H0 <- diff_beta_var(y + a, new_n[j, ] - y + b)
    # Expected utility
    Eu[j] <- m * sum(H1) + (1 - m) * sum(H0)
  }
  # Change in utility
  return(U - Eu)
}


#' Calculate change in utility using entropy
#'
#' @param n Sample size in each arm
#' @param y Response in each arm
#' @param a Prior shape
#' @param b Prior shape
#' @param ... Other function arguments
#' @examples
#' ascotsims:::utility_ent(rep(0, 4), rep(0, 4))
#' ascotsims:::utility_ent(c(20, 19, 18, 17), c(8,9,10,11))
utility_ent <- function(n, y, a = 1, b = 1, ...) {
  K <- length(n)
  # Current utility
  H <- ent_best_arm(y + a, n - y + b, ...)
  # Possible allocations
  new_n <- t(n + diag(1, K))
  # If response to observed from allocation
  new_y <- t(y + diag(1, K))
  # For each arm, what is the expected gain in utility
  Eu <- numeric(K)
  for(j in 1:K) {
    # predicted probability of event
    m <- (a + y[j]) / (a + b + n[j])
    # utility if event occurs
    H1 <- ent_best_arm(new_y[j, ] + a, new_n[j, ] - new_y[j, ] + b, ...)
    # utility if no event occurs
    H0 <- ent_best_arm(y + a, new_n[j, ] - y + b, ...)
    # Expected utility
    Eu[j] <- m * H1 + (1 - m) * H0
  }
  # Change in utility
  return(H - Eu)
}



#' Generate potential binary outcomes
#'
#' @param n Sample size
#' @param p Vector of probabilities of outcome
#' @param accrual_rate A function target integer argument which
#' generates that number of arrival times.
#' @param delay The time until endpoint is observed.
#' @param ... Optional additional arguments to the accrual_rate function.
#' @return A matrix of potential outcomes for each arm
#' @importFrom stats rbinom
#' @export
gen_potential_outcomes <- function(
  n,
  p,
  accrual_rate = function(n) (1:n)/2,
  delay = 14,
  ...) {

  k <- length(p)
  tim <- sort(accrual_rate(n, ...))
  obs <- tim + delay
  cbind(
    t_acc = tim,
    t_obs = obs,
    matrix(
      rbinom(n*k, 1, p), n, k, byrow = T,
      dimnames = list('id' = 1:n, 'o' = paste0('y_', 1:k))
    ))
}


#' Normalise vector to sum to 1
#'
#' @param x A numeric vector
#' @return The vector `x` normalised to sum to 1.
normalise <- function(x) x / sum(x)


#' Probability each column is maximum (or minimum)
#'
#' @param mat A matrix of draws.
#' @param minimum If `TRUE` determine best based on minimum instead of maximum value.
#' @export
prob_best <- function(mat, minimum = F) {
  if(!minimum) {
    as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
  } else {
    as.numeric(prop.table(table(factor(max.col(-mat), levels = 1:ncol(mat)))))
  }
}


#' BRAR with fixed ratio to control
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param h Scaling factor
fix_ctrl_brar <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!any(active[-1])) {
    w[1] <- 1
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- 1 / sum(active)
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}


#' BRAR with fixed ratio to control
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param h Scaling factor
fix_ctrl_brar2 <- function(quantity, active, h, fixq) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!any(active[-1])) {
    w[1] <- 1
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- fixq
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}


#' BRAR with matching control on biggest sample size (Trippa)
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param n Sample size in each arm
#' @param h Scaling factor
#' @param bpar Control matching parameter
match_ctrl_brar <- function(quantity, active, n, h, bpar = 1) {
  arms <- length(active) + 1
  w <- numeric(arms)
  if(!any(active)) {
    w[1] <- 1
    return(w)
  } else {
    w[-1] <- quantity ^ h
    w[-1][!active] <- 0
    C <- sum(w[-1]) / sum(active)
    w[1] <- C * exp(bpar * (max(n[-1]) - n[1]))
    w <- normalise(w)
    return(w)
  }
}


#' BRAR with constant proportion to control
#'
#' @param quantity The quantity to be used for allocation ratios.
#' @param active Flag indicating active arms.
#' @param h Scaling factor value
#' @return Allocation probabilities to each arm.
#' @examples
#' ascotsims:::const_ctrl_brar(runif(5), sample(c(TRUE, FALSE), 5, rep = TRUE), 1)
#' @export
const_ctrl_brar <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- 1 / arms
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}

#' BRAR on all arms.
#'
#' @param quantity The quantity to be used for allocation ratios.
#' @param active Flag indicating active arms.
#' @param h Scaling factor value
#' @return Allocation probabilities to each arm.
#' @examples
#' ascotsims:::brar_all(runif(5), sample(c(TRUE, FALSE), 5, rep = TRUE), 1)
brar_all <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else {
    w <- quantity ^ h
    w[!active] <- 0
    w <- normalise(w)
  }
  return(w)
}



brar <- function(quantity, active, h, min_rar, miter = 10, tol = 1e-8) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
    w[1] <- 1 # If nothing active, then reactivate SoC
  } else if(active[1] & !any(active[-1])) {
    w <- rep(0, arms)
    w[1] <- 1
  } else  {
    w <- quantity ^ h
    w[!active] <- 0
    w <- w / sum(w)
  }
  iter <- 1
  # Fixed point updates to achieve min_rar
  while(iter < miter & abs(min(w)-min_rar) > tol) {
    below_min <- w < min_rar & active
    w[below_min] <- min_rar
    w[!below_min] <- w[!below_min]*(1 - sum(below_min)*min_rar)/sum(w[!below_min])
  }
  return(w)
}


#' @export
brar_fix <- function(quantity, active, h, min_rar, miter = 10, tol = 1e-8) {
  arms <- length(active)
  active_arms <- sum(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
    w[1] <- 1 # If nothing active, then reactivate SoC
  } else if(active[1] & !any(active[-1])) {
    w <- rep(0, arms)
    w[1] <- 1
  } else {
    w <- quantity^h
    w[1] <- 1 / active_arms
    w[!active] <- 0
    w[-1] <- w[-1] / sum(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  iter <- 1
  # Fixed point updates to achieve min_rar
  while(iter < miter & abs(min(w)-min_rar) > tol) {
    below_min <- c(FALSE, w[-1] < min_rar & active[-1])
    w[below_min] <- min_rar
    w[!below_min] <- w[!below_min]*(1 - sum(below_min)*min_rar)/sum(w[!below_min])
    iter <- iter + 1
  }
  return(w)
}


#' Variational approximation to logistic regression with Normal prior.
#'
#' @param X Design matrix
#' @param y Response events
#' @param n Response trials
#' @param M0 Prior normal mean on regression parameter
#' @param S0 Prior normal variance on regression parameter
#' @param ... Additional arguments to `varapproxr::vb_logistic_n`
#' @return A list giving results of ELBO optimised variational approximation.
#' @importFrom varapproxr vb_logistic_n
#' @export
vb_mod <- function(
  X = diag(1, 4),
  y,
  n,
  M0 = rep(0, ncol(X)),
  S0 = diag(2^2, ncol(X)),
  ...
  ) {
  vb_logistic_n(
    X, y, n,
    mu0 = M0, Sigma0 = S0,
    mu_init = M0, Sigma_init = diag(1, ncol(X)),
    maxiter_jj = 100, alg = "sj", ...)
}

#' Correct RAR
#'
#' @param p The current regimen allocations
#' @param min_p The targeted allocation for SoC
#' @param X The indicator design for regimens
#' @param tol The tolerance level
#' @param maxiter The maximum iterations
correct_rar <- function(p, min_p, X, tol = 1e-4, maxiter = 20) {
  pnew <- p
  pagg <- drop(p %*% XI)
  paggnew <- pagg
  soc <- pagg[grepl("0", colnames(XI))]
  iter <- 1
  id_a <- grepl("a[1-9]", names(pnew))
  id_b <- grepl("b[1-9]", names(pnew))
  id_c <- grepl("c[1-9]", names(pnew))
  while(iter <= maxiter & any(abs(soc - min_p) > tol)) {
    pnew[id_a] <- (1 - min_p[1])*pnew[id_a]/sum(pnew[id_a])
    pnew[!id_a] <- min_p[1]*pnew[!id_a]/sum(pnew[!id_a])
    pnew[id_b] <- (1 - min_p[2])*pnew[id_b]/sum(pnew[id_b])
    pnew[!id_b] <- min_p[2]*pnew[!id_b]/sum(pnew[!id_b])
    pnew[id_c] <- (1-min_p)[3]*pnew[id_c]/sum(pnew[id_c])
    pnew[!id_c] <- min_p[3]*pnew[!id_c]/sum(pnew[!id_c])
    paggnew <- drop(pnew %*% XI)
    soc <- paggnew[grepl("0", colnames(XI))]
    iter <- iter + 1
  }
  return(pnew)
}

#' Correct RAR to achieve minimum on SoC
#'
#' @param p The current regimen allocations
#' @param min_p The minimum allocation to SoC
#' @param X The indicator design for regimens
#' @param tol The tolerance level
#' @param maxiter The maximum iterations
correct_min_rar_soc <- function(p, min_p, act, X, tol = 1e-4, maxiter = 20) {
  pnew <- p
  pagg <- drop(p %*% X)
  paggnew <- pagg
  soc <- pagg[grepl("0", colnames(X))]
  correct <- (soc < min_p) * act
  iter <- 1
  id_a <- grepl("a[1-9]", names(pnew))
  id_b <- grepl("b[1-9]", names(pnew))
  id_c <- grepl("c[1-9]", names(pnew))
  while(iter <= maxiter & any(abs(soc - min_p)*correct > tol)) {
    if(correct[1]) {
      pnew[id_a] <- (1 - min_p[1])*pnew[id_a]/sum(pnew[id_a])
      pnew[!id_a] <- min_p[1]*pnew[!id_a]/sum(pnew[!id_a])
    }
    if(correct[2]) {
      pnew[id_b] <- (1 - min_p[2])*pnew[id_b]/sum(pnew[id_b])
      pnew[!id_b] <- min_p[2]*pnew[!id_b]/sum(pnew[!id_b])
    }
    if(correct[3]) {
      pnew[id_c] <- (1-min_p)[3]*pnew[id_c]/sum(pnew[id_c])
      pnew[!id_c] <- min_p[3]*pnew[!id_c]/sum(pnew[!id_c])
    }
    paggnew <- drop(pnew %*% X)
    soc <- paggnew[grepl("0", colnames(X))]
    iter <- iter + 1
  }
  return(pnew)
}

#' Correct RAR to achieve minimum across all interventions
#'
#' @param p The current regimen allocations
#' @param min_p The minimum allocation to any active arm in the domain
#' @param act Active intervention status
#' @param X The indicator design for regimens
#' @param tol The tolerance level
#' @param maxiter The maximum iterations
#' @export
correct_min_rar_any <- function(p, min_p, act, X, tol = 1e-4, maxiter = 20) {
  pnew <- p
  pagg <- drop(p %*% X)
  paggnew <- pagg

  labs <- c("a", "b", "c")
  min_pl <- rep(min_p, sapply(labs, function(x) sum(grepl(x, names(pagg)))))

  correct <- (pagg < min_pl) * act
  iter <- 1
  while(iter <= maxiter & any(abs(paggnew - min_pl)*correct > tol) & sum(correct) >= 1) {
    id_correct <- which(correct == 1)
    id_reg <- sapply(names(id_correct), function(x) grepl(x, names(p)))
    for(i in 1:sum(correct)) {
      pnew[id_reg[, i]] <- min_pl[id_correct[i]]*pnew[id_reg[, i]]/sum(pnew[id_reg[, i]])
      pnew[!id_reg[, i]] <- (1 - min_pl[id_correct[i]])*pnew[!id_reg[, i]]/sum(pnew[!id_reg[, i]])
    }
    paggnew <- drop(pnew %*% X)
    iter <- iter + 1
    correct <- (paggnew < min_pl) * act
  }
  return(pnew)
}
