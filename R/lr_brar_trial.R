#' Run a logistic regression BRAR trial.
#'
#' @param n_seq Sequence of interim analyses
#' @param full_dat   Potential outcomes data from `gen_potential_outcomes`
#' @param X Design matrix
#' @param M0 Prior mean
#' @param S0 Prior variance
#' @param impute Impute missing outcomes using early event information
#' @param impute_sim Number of imputation draws
#' @param early_t Time to early outcome (should be less than value of `delay` in `gen_potential_outcomes`)
#' @param trans_prob Probability of transitioning from no early event to a final event
#' @param perm_drop Permanently drop arms
#' @param rar Turn on RAR
#' @param rar_control Include control arm in RAR or fix
#' @param hpar1 Parameter to control amount of RAR
#' @param hpar2 Parameter to control scaling of RAR
#' @param h Function for RAR
#' @param fpar1 Parameter to control lack of effect boundary (futility)
#' @param fpar2 Parameter to control lack of effect boundary (futility)
#' @param f Function for utility boundary
#' @param gpar1 Parameter to control most effective boundary (success)
#' @param gpar2 Parameter to control most effective boundary (success)
#' @param g Function for most effective boundary
#' @param supr_ref If dropping control, what threshold must an intervention arm achieve in terms of superiority
#' @param mc_draw Number of Monte Carlo draws form approximating posterior
#' @param ... Other function arguments
#' @return A list of the trial quantities
#'
#' @importFrom HDInterval hdi
#' @import matrixStats
#' @importFrom mvnfast rmvn
#' @importFrom stats plogis
#'
#' @export
lr_brar_trial <- function(
  n_seq,
  full_dat,
  X,
  M0,
  S0,
  impute = F,
  impute_sim = 100,
  early_t = 7,
  trans_prob = 0.75,
  perm_drop = F,
  rar = F,
  rar_control = F,
  hpar1 = 0.5,
  hpar2 = 0,
  h = function(m) hpar1 * (m / length(n_seq)) ^ hpar2,
  fpar1 = 0,
  fpar2 = 0,
  f = function(m) fpar1 * (m / length(n_seq)) ^ fpar2,
  gpar1 = 1 - 0.95^3,
  gpar2 = 0,
  g = function(m) 1 - gpar1 * (m / length(n_seq)) ^ gpar2,
  supr_ref = 0.99,
  mc_draw = 2e4,
  ...
) {

  if(max(n_seq) != nrow(full_dat)) stop("dimension mismatch, max(n_seq) should equal nrow(dat)")

  tdat <- full_dat[, 1:2]
  dat <- full_dat[, -(1:2)]

  if(nrow(X) != ncol(dat)) stop("Design matrix inconsistent with full_dat")

  P <- ncol(X)
  Xt <- t(X)
  lr_mod <- function(y, n) {
    vb_mod(X = X, y, n, S0 = S0)
  }


  # Interim schedule
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_seq_aug <- c(0, n_seq)
  n_new <- diff(n_seq_aug)
  A <- ncol(dat)
  x <- numeric(n_max)
  y <- rep(NA_real_, n_max)
  z <- rep(NA_real_, n_max)
  n_enrolled <- 0

  arm_names <- rownames(X)
  par_names <- colnames(X)
  dm <- list("interim" = 1:n_int, "arm" = arm_names)
  dm2 <- list("interim" = 1:n_int, "arm" = par_names)

  n_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  y_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  z_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  n_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  y_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  z_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  n_agg_ear <- matrix(0, n_int, A, dimnames = dm)
  z_agg_ear <- matrix(0, n_int, A, dimnames = dm)
  p_trans <- matrix(0, n_int, A, dimnames = dm)

  p_supr <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  is_supr <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  p_best_trt <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  p_best <- matrix(0, n_int, A, dimnames = dm)
  p_best_active <- matrix(0, n_int, A, dimnames = dm)
  p_rand <- matrix(1/A, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = arm_names))
  a_active <- matrix(TRUE, n_int, A, dimnames = dm)
  is_best <- matrix(TRUE, n_int, A, dimnames = dm)
  c_best <- matrix(TRUE, n_int, A, dimnames = dm)
  a_best_by <- matrix(TRUE, n_int, A, dimnames = dm)
  ever_drop_by <- matrix(TRUE, n_int, A, dimnames = dm)

  beta_mu <- matrix(0, n_int, P, dimnames = dm2)
  beta_lb <- matrix(0, n_int, P, dimnames = dm2)
  beta_ub <- matrix(0, n_int, P, dimnames = dm2)
  beta_eff <- matrix(0, n_int, P, dimnames = dm2)
  arm_mu_lo <- matrix(0, n_int, A, dimnames = dm)
  arm_mu_p <- matrix(0, n_int, A, dimnames = dm)
  arm_var_lo <- matrix(0, n_int, A, dimnames = dm)
  arm_var_p <- matrix(0, n_int, A, dimnames = dm)
  eff_mu_p <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  eff_var_p <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  arm_mu_lb <- matrix(0, n_int, A, dimnames = dm)
  arm_mu_ub <- matrix(0, n_int, A, dimnames = dm)

  # Loop interim analyses
  for(i in 1:n_int) {

    t_ref <- tdat[, 2][n_seq[i]]
    n_enrol <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
    id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
    n_enrolled <- n_enrolled + n_enrol

    new_x <- factor(sample.int(A, n_enrol, prob = p_rand[i, ], replace = T), levels = 1:A)
    new_y <- dat[cbind(id_enrolled, new_x)]
    new_z <- (new_y == 1) * rbinom(n_enrol, 1, trans_prob)

    new_y_agg <- aggregate(new_y, list(new_x), sum, drop = F)[, 2]
    new_n_agg <- table(new_x, dnn = NULL)

    x[id_enrolled] <- new_x
    y[id_enrolled] <- new_y
    z[id_enrolled] <- new_z
    id_analysis <- which(tdat[, 2] <= t_ref)

    # What is imputed?
    id_impute_early <- which(tdat[, 2] - t_ref > 0 & tdat[, 2] - t_ref <= early_t)
    id_impute_late  <- which(tdat[, 1] <= t_ref & tdat[, 2] - t_ref > early_t)

    n_agg_obs[i, ] <- table(x[id_analysis], dnn = NULL)
    n_agg_enr[i, ] <- table(x[1:n_enrolled], dnn = NULL)
    n_agg_ear[i, ] <- table(x[c(id_analysis, id_impute_early)], dnn = NULL)
    y_agg_obs[i, ] <- aggregate(y[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    y_agg_enr[i, ] <- aggregate(y[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
    z_agg_ear[i, ] <- aggregate(z[c(id_analysis, id_impute_early)], by = list(x[c(id_analysis, id_impute_early)]), FUN = sum, drop = F)[, 2]
    z_agg_obs[i, ] <- aggregate(z[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    z_agg_enr[i, ] <- aggregate(z[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]

    # Transitions
    w01 <- aggregate(y[id_analysis] == 1 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    w00 <- aggregate(y[id_analysis] == 0 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    p_trans[i, ] <- (1 + w01) / (2 + w01 + w00)


    # Do the imputation here
    if(impute) {
      ptran <- matrix(rbeta(impute_sim, 1 + w01, 1 + w00), impute_sim, A, byrow = T)
      presp <- matrix(0.25, impute_sim + 1, A)
      y_impute <- y[1:n_enrolled]
      # Cycle the imputation
      mult <- mc_draw / impute_sim
      b_draws <- matrix(0, mult*impute_sim, P)
      for(j in 1:impute_sim) {
        y_impute[id_impute_early] <- rbinom(length(id_impute_early), 1, ptran[j, x[id_impute_early]])
        y_impute[id_impute_early][z[id_impute_early] == 1] <- 1
        y_impute[id_impute_late] <- rbinom(length(id_impute_late), 1, presp[j, x[id_impute_late]])
        y_agg <- aggregate(y_impute, by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
        mod <- lr_mod(y = y_agg, n = n_agg_enr[i, ])
        b_draws[(mult*(j-1) + 1):(mult*j), ] <- mvnfast::rmvn(mult, mod$mu, mod$Sigma)
        presp[j + 1, ] <- plogis(b_draws[mult*j, , drop = F] %*% Xt)
      }
    } else {
      mod <- lr_mod(y = y_agg_obs[i, ], n = n_agg_obs[i, ])
      b_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    }

    m_draws <- b_draws %*% Xt
    p_draws <- plogis(m_draws)
    eff_draws <- sweep(p_draws[, -1], 1, p_draws[, 1])

    p_supr[i, ] <- matrixStats::colMeans2(eff_draws < 0)
    is_supr[i, ] <- p_supr[i, ] > supr_ref
    p_best[i, ] <- prob_best(m_draws, minimum = TRUE)
    p_best_trt[i, ] <- prob_best(m_draws[, -1], minimum = TRUE)

    is_best[i, ] <- p_best[i, ] > g(i)
    a_best_by[i, ] <- apply(is_best[1:i, , drop = F], 2, any)
    c_best[i, ] <- p_best[i, ] == max(p_best[i, ])

    # Include control group in RAR or not?
    if(rar_control) {

      a_active[i, ] <-  p_best[i, ] >= f(i)
      if(perm_drop & i > 1) {
        a_active[i, ] <- a_active[i, ] & a_active[i - 1, ]
      }
      if(rar) p_rand[i + 1, ] <- brar_all(p_best[i, ], a_active[i, ], h(i))

    } else {

      # If control not in RAR, still drop if something superior
      a_active[i, -1] <-  p_best[i, -1] >= f(i)
      a_active[i, 1] <- !any(p_supr[i, ] > supr_ref)
      if(perm_drop & i > 1) {
        a_active[i, ] <- a_active[i, ] & a_active[i - 1, ]
      }
      if(rar) p_rand[i + 1, ] <- const_ctrl_brar(p_best[i, ], a_active[i, ], h(i))

    }

    p_best_active[i, a_active[i, ]] <- prob_best(m_draws[, a_active[i, ], drop = F], minimum = TRUE)
    ever_drop_by[i, ] <- apply(!a_active[1:i, , drop = F], 2, any)

    # Parameter summaries
    hdival_p <- HDInterval::hdi(p_draws)
    hdival_b <- HDInterval::hdi(b_draws)
    beta_mu[i, ] <- matrixStats::colMeans2(b_draws)
    beta_lb[i, ] <- hdival_b[1, ]
    beta_ub[i, ] <- hdival_b[2, ]
    beta_eff[i, ] <- matrixStats::colMeans2(b_draws < 0)
    arm_mu_p[i, ] <- matrixStats::colMeans2(p_draws)
    arm_var_p[i, ] <- matrixStats::colVars(p_draws)
    eff_mu_p[i, ] <- matrixStats::colMeans2(eff_draws)
    eff_var_p[i, ] <- matrixStats::colVars(eff_draws)
    arm_mu_lb[i, ] <- hdival_p[1, ]
    arm_mu_ub[i, ] <- hdival_p[2, ]
  }
  p_rand <- p_rand[-1, , drop = F]

  # Interim summaries
  # n_enr <- matrix(apply(n_agg_enr, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # n_obs <- matrix(apply(n_agg_obs, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # y_enr <- matrix(apply(y_agg_enr, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # y_obs <- matrix(apply(y_agg_obs, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # best_arm <- matrix(apply(is_best, 1, function(x) ifelse(any(x), which(x), 0)), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))

  # # Arm summaries
  # announce_supr <- matrix(suppressWarnings(apply(is_supr, 2, function(x) {min(which(x == TRUE))})),
  #                         1, A - 1, dimnames = list("interim" = as.character(i), "arm" = arm_names[-1]))
  # announce_best <- matrix(suppressWarnings(apply(is_best, 2, function(x) {min(which(x == TRUE))})),
  #                         1, A, dimnames = list("interim" = as.character(i), "arm" = arm_names))
  # drop_at <- matrix(suppressWarnings(apply(a_active, 2, function(x) {min(which(x == FALSE))})),
  #                   1, A, dimnames = list("interim" = as.character(i), "arm" = arm_names))

  # is_best <- cbind(is_best, "Total" = apply(is_best, 1, any))

  trial_quant <- nlist(
    n_agg_enr, y_agg_enr, z_agg_enr,
    n_agg_obs, y_agg_obs, z_agg_obs,
    n_agg_ear, z_agg_ear,
    p_trans,
    p_supr, is_supr,
    p_best, is_best,
    p_best_trt, p_best_active,
    p_rand, a_active, a_best_by, c_best, ever_drop_by,
    beta_mu, beta_lb, beta_ub, beta_eff,
    arm_mu_p, eff_mu_p,
    arm_mu_lb, arm_mu_ub,
    arm_var_p, eff_var_p
  )

  model_quant <- nlist(

  )

  return(nlist(
    n_agg_enr, y_agg_enr, z_agg_enr,
    n_agg_obs, y_agg_obs, z_agg_obs,
    n_agg_ear, z_agg_ear,
    p_trans,
    p_supr, is_supr,
    p_best, is_best,
    p_best_trt, p_best_active,
    p_rand, a_active, a_best_by, c_best, ever_drop_by,
    beta_mu, beta_lb, beta_ub, beta_eff,
    arm_mu_p, eff_mu_p,
    arm_mu_lb, arm_mu_ub,
    arm_var_p, eff_var_p
    # n_enr, n_obs,
    # y_enr, y_obs,
    # best_arm,
    # announce_best, announce_supr, drop_at
  ))
}

#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_trial_data <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    lapply(1:length(i),
           function(x) {
             tidyr::gather(tidyr::as_tibble(i[[x]], rownames = "interim"), "arm", !!names(i)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "arm")) %>%
      dplyr::mutate(arm = forcats::fct_inorder(arm)) %>%
      dplyr::arrange(interim, arm)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial))
}
#' Run a logistic regression BRAR trial.
#'
#' @param n_seq Sequence of interim analyses
#' @param full_dat   Potential outcomes data from `gen_potential_outcomes`
#' @param X Design matrix
#' @param M0 Prior mean
#' @param S0 Prior variance
#' @param impute Impute missing outcomes using early event information
#' @param impute_sim Number of imputation draws
#' @param early_t Time to early outcome (should be less than value of `delay` in `gen_potential_outcomes`)
#' @param trans_prob Probability of transitioning from no early event to a final event
#' @param perm_drop Permanently drop arms
#' @param rar Turn on RAR
#' @param rar_control Include control arm in RAR or fix
#' @param hpar1 Parameter to control amount of RAR
#' @param hpar2 Parameter to control scaling of RAR
#' @param h Function for RAR
#' @param fpar1 Parameter to control lack of effect boundary (futility)
#' @param fpar2 Parameter to control lack of effect boundary (futility)
#' @param f Function for utility boundary
#' @param gpar1 Parameter to control most effective boundary (success)
#' @param gpar2 Parameter to control most effective boundary (success)
#' @param g Function for most effective boundary
#' @param supr_ref If dropping control, what threshold must an intervention arm achieve in terms of superiority
#' @param mc_draw Number of Monte Carlo draws form approximating posterior
#' @param ... Other function arguments
#' @return A list of the trial quantities
#'
#' @importFrom HDInterval hdi
#' @import matrixStats
#' @importFrom mvnfast rmvn
#' @importFrom stats plogis
#'
#' @export
lr_brar_trial <- function(
  n_seq,
  full_dat,
  X,
  M0,
  S0,
  impute = F,
  impute_sim = 100,
  early_t = 7,
  trans_prob = 0.75,
  perm_drop = F,
  rar = F,
  rar_control = F,
  hpar1 = 0.5,
  hpar2 = 0,
  h = function(m) hpar1 * (m / length(n_seq)) ^ hpar2,
  fpar1 = 0,
  fpar2 = 0,
  f = function(m) fpar1 * (m / length(n_seq)) ^ fpar2,
  gpar1 = 1 - 0.95^3,
  gpar2 = 0,
  g = function(m) 1 - gpar1 * (m / length(n_seq)) ^ gpar2,
  supr_ref = 1,
  mc_draw = 2e4,
  ...
) {

  if(max(n_seq) != nrow(full_dat)) stop("dimension mismatch, max(n_seq) should equal nrow(dat)")

  tdat <- full_dat[, 1:2]
  dat <- full_dat[, -(1:2)]

  if(nrow(X) != ncol(dat)) stop("Design matrix inconsistent with full_dat")

  P <- ncol(X)
  Xt <- t(X)
  lr_mod <- function(y, n) {
    vb_mod(X = X, y, n, S0 = S0)
  }


  # Interim schedule
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_seq_aug <- c(0, n_seq)
  n_new <- diff(n_seq_aug)
  A <- ncol(dat)
  x <- numeric(n_max)
  y <- rep(NA_real_, n_max)
  z <- rep(NA_real_, n_max)
  n_enrolled <- 0

  arm_names <- rownames(X)
  par_names <- colnames(X)
  dm <- list("interim" = 1:n_int, "arm" = arm_names)
  dm2 <- list("interim" = 1:n_int, "arm" = par_names)

  n_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  y_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  z_agg_enr <- matrix(0, n_int, A, dimnames = dm)
  n_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  y_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  z_agg_obs <- matrix(0, n_int, A, dimnames = dm)
  n_agg_ear <- matrix(0, n_int, A, dimnames = dm)
  z_agg_ear <- matrix(0, n_int, A, dimnames = dm)
  p_trans <- matrix(0, n_int, A, dimnames = dm)

  p_supr <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  is_supr <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  p_best_trt <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  p_best <- matrix(0, n_int, A, dimnames = dm)
  p_best_active <- matrix(0, n_int, A, dimnames = dm)
  p_rand <- matrix(1/A, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = arm_names))
  a_active <- matrix(TRUE, n_int, A, dimnames = dm)
  is_best <- matrix(TRUE, n_int, A, dimnames = dm)
  c_best <- matrix(TRUE, n_int, A, dimnames = dm)
  a_best_by <- matrix(TRUE, n_int, A, dimnames = dm)
  ever_drop_by <- matrix(TRUE, n_int, A, dimnames = dm)

  beta_mu <- matrix(0, n_int, P, dimnames = dm2)
  beta_lb <- matrix(0, n_int, P, dimnames = dm2)
  beta_ub <- matrix(0, n_int, P, dimnames = dm2)
  beta_eff <- matrix(0, n_int, P, dimnames = dm2)
  arm_mu_lo <- matrix(0, n_int, A, dimnames = dm)
  arm_mu_p <- matrix(0, n_int, A, dimnames = dm)
  arm_var_lo <- matrix(0, n_int, A, dimnames = dm)
  arm_var_p <- matrix(0, n_int, A, dimnames = dm)
  eff_mu_p <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  eff_var_p <- matrix(0, n_int, A - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  arm_mu_lb <- matrix(0, n_int, A, dimnames = dm)
  arm_mu_ub <- matrix(0, n_int, A, dimnames = dm)

  # Loop interim analyses
  for(i in 1:n_int) {

    t_ref <- tdat[, 2][n_seq[i]]
    n_enrol <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
    id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
    n_enrolled <- n_enrolled + n_enrol

    new_x <- factor(sample.int(A, n_enrol, prob = p_rand[i, ], replace = T), levels = 1:A)
    new_y <- dat[cbind(id_enrolled, new_x)]
    new_z <- (new_y == 1) * rbinom(n_enrol, 1, trans_prob)

    new_y_agg <- aggregate(new_y, list(new_x), sum, drop = F)[, 2]
    new_n_agg <- table(new_x, dnn = NULL)

    x[id_enrolled] <- new_x
    y[id_enrolled] <- new_y
    z[id_enrolled] <- new_z
    id_analysis <- which(tdat[, 2] <= t_ref)

    # What is imputed?
    id_impute_early <- which(tdat[, 2] - t_ref > 0 & tdat[, 2] - t_ref <= early_t)
    id_impute_late  <- which(tdat[, 1] <= t_ref & tdat[, 2] - t_ref > early_t)

    n_agg_obs[i, ] <- table(x[id_analysis], dnn = NULL)
    n_agg_enr[i, ] <- table(x[1:n_enrolled], dnn = NULL)
    n_agg_ear[i, ] <- table(x[c(id_analysis, id_impute_early)], dnn = NULL)
    y_agg_obs[i, ] <- aggregate(y[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    y_agg_enr[i, ] <- aggregate(y[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
    z_agg_ear[i, ] <- aggregate(z[c(id_analysis, id_impute_early)], by = list(x[c(id_analysis, id_impute_early)]), FUN = sum, drop = F)[, 2]
    z_agg_obs[i, ] <- aggregate(z[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    z_agg_enr[i, ] <- aggregate(z[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]

    # Transitions
    w01 <- aggregate(y[id_analysis] == 1 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    w00 <- aggregate(y[id_analysis] == 0 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    p_trans[i, ] <- (1 + w01) / (2 + w01 + w00)


    # Do the imputation here
    if(impute) {
      ptran <- matrix(rbeta(impute_sim, 1 + w01, 1 + w00), impute_sim, A, byrow = T)
      presp <- matrix(0.25, impute_sim + 1, A)
      y_impute <- y[1:n_enrolled]
      # Cycle the imputation
      mult <- mc_draw / impute_sim
      b_draws <- matrix(0, mult*impute_sim, P)
      for(j in 1:impute_sim) {
        y_impute[id_impute_early] <- rbinom(length(id_impute_early), 1, ptran[j, x[id_impute_early]])
        y_impute[id_impute_early][z[id_impute_early] == 1] <- 1
        y_impute[id_impute_late] <- rbinom(length(id_impute_late), 1, presp[j, x[id_impute_late]])
        y_agg <- aggregate(y_impute, by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
        mod <- lr_mod(y = y_agg, n = n_agg_enr[i, ])
        b_draws[(mult*(j-1) + 1):(mult*j), ] <- mvnfast::rmvn(mult, mod$mu, mod$Sigma)
        presp[j + 1, ] <- plogis(b_draws[mult*j, , drop = F] %*% Xt)
      }
    } else {
      mod <- lr_mod(y = y_agg_obs[i, ], n = n_agg_obs[i, ])
      b_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    }

    m_draws <- b_draws %*% Xt
    p_draws <- plogis(m_draws)
    eff_draws <- sweep(p_draws[, -1], 1, p_draws[, 1])

    p_supr[i, ] <- matrixStats::colMeans2(eff_draws < 0)
    is_supr[i, ] <- p_supr[i, ] > supr_ref
    p_best[i, ] <- prob_best(m_draws, minimum = TRUE)
    p_best_trt[i, ] <- prob_best(m_draws[, -1], minimum = TRUE)

    is_best[i, ] <- p_best[i, ] > g(i)
    a_best_by[i, ] <- apply(is_best[1:i, , drop = F], 2, any)
    c_best[i, ] <- p_best[i, ] == max(p_best[i, ])

    # Include control group in RAR or not?
    if(rar_control) {

      a_active[i, ] <-  p_best[i, ] >= f(i)
      if(perm_drop & i > 1) {
        a_active[i, ] <- a_active[i, ] & a_active[i - 1, ]
      }
      if(rar) p_rand[i + 1, ] <- brar_all(p_best[i, ], a_active[i, ], h(i))

    } else {

      # If control not in RAR, still drop if something superior
      a_active[i, -1] <-  p_best[i, -1] >= f(i)
      a_active[i, 1] <- !any(p_supr[i, ] > supr_ref)
      if(perm_drop & i > 1) {
        a_active[i, ] <- a_active[i, ] & a_active[i - 1, ]
      }
      if(rar) p_rand[i + 1, ] <- const_ctrl_brar(p_best[i, ], a_active[i, ], h(i))

    }

    p_best_active[i, a_active[i, ]] <- prob_best(m_draws[, a_active[i, ], drop = F], minimum = TRUE)
    ever_drop_by[i, ] <- apply(!a_active[1:i, , drop = F], 2, any)

    # Parameter summaries
    hdival_p <- HDInterval::hdi(p_draws)
    hdival_b <- HDInterval::hdi(b_draws)
    beta_mu[i, ] <- matrixStats::colMeans2(b_draws)
    beta_lb[i, ] <- hdival_b[1, ]
    beta_ub[i, ] <- hdival_b[2, ]
    beta_eff[i, ] <- matrixStats::colMeans2(b_draws < 0)
    arm_mu_p[i, ] <- matrixStats::colMeans2(p_draws)
    arm_var_p[i, ] <- matrixStats::colVars(p_draws)
    eff_mu_p[i, ] <- matrixStats::colMeans2(eff_draws)
    eff_var_p[i, ] <- matrixStats::colVars(eff_draws)
    arm_mu_lb[i, ] <- hdival_p[1, ]
    arm_mu_ub[i, ] <- hdival_p[2, ]
  }
  p_rand <- p_rand[-1, , drop = F]

  # Interim summaries
  # n_enr <- matrix(apply(n_agg_enr, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # n_obs <- matrix(apply(n_agg_obs, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # y_enr <- matrix(apply(y_agg_enr, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # y_obs <- matrix(apply(y_agg_obs, 1, sum), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))
  # best_arm <- matrix(apply(is_best, 1, function(x) ifelse(any(x), which(x), 0)), n_int, 1, dimnames = list("interim" = as.character(1:n_int), "arm" = "Total"))

  # # Arm summaries
  # announce_supr <- matrix(suppressWarnings(apply(is_supr, 2, function(x) {min(which(x == TRUE))})),
  #                         1, A - 1, dimnames = list("interim" = as.character(i), "arm" = arm_names[-1]))
  # announce_best <- matrix(suppressWarnings(apply(is_best, 2, function(x) {min(which(x == TRUE))})),
  #                         1, A, dimnames = list("interim" = as.character(i), "arm" = arm_names))
  # drop_at <- matrix(suppressWarnings(apply(a_active, 2, function(x) {min(which(x == FALSE))})),
  #                   1, A, dimnames = list("interim" = as.character(i), "arm" = arm_names))

  # is_best <- cbind(is_best, "Total" = apply(is_best, 1, any))

  trial_quant <- nlist(
    n_agg_enr, y_agg_enr, z_agg_enr,
    n_agg_obs, y_agg_obs, z_agg_obs,
    n_agg_ear, z_agg_ear,
    p_trans,
    p_supr, is_supr,
    p_best, is_best,
    p_best_trt, p_best_active,
    p_rand, a_active, a_best_by, c_best, ever_drop_by,
    arm_mu_p, eff_mu_p,
    arm_mu_lb, arm_mu_ub,
    arm_var_p, eff_var_p
  )

  model_quant <- nlist(
    beta_mu, beta_lb, beta_ub, beta_eff
  )

  return(
    nlist(trial_quant, model_quant)
  )
}

#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_trial_data <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["trial_quant"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "arm", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "arm")) %>%
      dplyr::mutate(arm = forcats::fct_inorder(arm)) %>%
      dplyr::arrange(interim, arm)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial))
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_model_data <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["model_quant"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "parameter", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "parameter")) %>%
      dplyr::mutate(parameter = forcats::fct_inorder(parameter)) %>%
      dplyr::arrange(interim, parameter)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial))
}
