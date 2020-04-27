#' Run a logistic regression BRAR trial.
#'
#' @param n_seq Sequence of interim analyses
#' @param full_dat   Potential outcomes data from `gen_potential_outcomes`
#' @param activate_at Integer giving the interim at which each arm is activated
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
lr_brar_add_arms_trial <- function(
  n_seq,
  full_dat,
  activate_at,
  X,
  M0,
  S0,
  impute = F,
  impute_sim = 100,
  early_t = 7,
  trans_prob = 0.75,
  perm_drop = T,
  rar_control = F,
  drop_best = F,
  hpar1 = 0.5,
  hpar2 = 0,
  h = function(m) hpar1 * (m / length(n_seq)) ^ hpar2,
  fpar1 = 0,
  fpar2 = 0,
  f = function(m) fpar1 * (m / length(n_seq)) ^ fpar2,
  gpar1 = 0.05,
  gpar2 = 0,
  g = function(m) 1 - gpar1 * (m / length(n_seq)) ^ gpar2,
  supr_thres = 0.975,
  mc_draw = 2e4,
  ...
) {

  if(max(n_seq) != nrow(full_dat)) stop("dimension mismatch, max(n_seq) should equal nrow(dat)")

  tdat <- full_dat[, 1:2]
  dat <- full_dat[, -(1:2)]

  if(nrow(X) != ncol(dat)) stop("Design matrix inconsistent with full_dat")

  # Interim schedule
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_seq_aug <- c(0, n_seq)
  n_new <- diff(n_seq_aug)
  Apar <- ncol(X)
  Aarm <- nrow(X)

  x <- factor(numeric(n_max), levels = 1:Aarm)
  y <- rep(NA_real_, n_max)
  z <- rep(NA_real_, n_max)
  n_enrolled <- 0

  arm_names <- rownames(X)
  par_names <- colnames(X)
  dm_arm <- list("interim" = 1:n_int, "arm" = arm_names)
  dm_par <- list("interim" = 1:n_int, "parameter" = par_names)

  n_agg_enr <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  y_agg_enr <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  z_agg_enr <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  n_agg_obs <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  y_agg_obs <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  z_agg_obs <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  n_agg_ear <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  z_agg_ear <- matrix(0, n_int, Aarm, dimnames = dm_arm)
  p_trans <- matrix(0, n_int, Aarm, dimnames = dm_arm)

  # Values at interim analysis
  p_supr <- matrix(NA_real_, n_int, Aarm - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  p_best <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)
  p_best_active <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)

  # State at interim analysis
  is_supr <- matrix(FALSE, n_int, Aarm - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  is_best <- matrix(FALSE, n_int, Aarm, dimnames = dm_arm)
  is_current_best <- matrix(FALSE, n_int, Aarm, dimnames = dm_arm)
  is_active <- matrix(1:n_int, n_int, Aarm, dimnames = dm_arm)
  is_active <- sweep(is_active, 2, activate_at, ">=")
  p_rand <- matrix(0, n_int + 1, Aarm, dimnames = list("interim" = 0:n_int, "arm" = arm_names))
  p_rand[1, is_active[1, ]] <- 1 / sum(is_active[1, ])

  beta_mu <- matrix(NA_real_, n_int, Apar, dimnames = dm_par)
  beta_lb <- matrix(NA_real_, n_int, Apar, dimnames = dm_par)
  beta_ub <- matrix(NA_real_, n_int, Apar, dimnames = dm_par)
  beta_eff <- matrix(NA_real_, n_int, Apar, dimnames = dm_par)

  arm_mu_p <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)
  arm_var_p <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)
  arm_mu_lb <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)
  arm_mu_ub <- matrix(NA_real_, n_int, Aarm, dimnames = dm_arm)

  eff_mu_p <- matrix(NA_real_, n_int, Aarm - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))
  eff_var_p <- matrix(NA_real_, n_int, Aarm - 1, dimnames = list("interim" = 1:n_int, "arm" = arm_names[-1]))

  # Loop interim analyses
  for(i in 1:n_int) {

    if(any(activate_at == i) & i > 1) {
      activated <- which(activate_at == i)
      n_activated <- sum(activate_at == i)
      n_active <- sum(is_active[i-1, ])
      n_total <- n_active + n_activated
      p_rand[i, is_active[i - 1, ]] <- p_rand[i, is_active[i-1, ]]*n_active/n_total
      p_rand[i, activated] <- (1 - sum(p_rand[i, ]))/n_activated
      is_active[i - 1, activated] <- TRUE
    }
    activated_arms <- which(activate_at <= i)
    active_arms <- which(is_active[i, ])
    Xactivated <- X[activated_arms, ]
    activated_pars <- colSums(Xactivated) > 0
    Xactivated <- Xactivated[activated_arms, activated_pars]

    t_ref <- tdat[, 2][n_seq[i]]
    n_enrol <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
    id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
    n_enrolled <- n_enrolled + n_enrol

    new_x <- factor(sample.int(Aarm, n_enrol, prob = p_rand[i, ], replace = T), levels = 1:Aarm)
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
    y_agg_obs[i, ] <- ifelse(is.na(y_agg_obs[i, ]), 0, y_agg_obs[i, ])
    y_agg_enr[i, ] <- aggregate(y[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
    y_agg_enr[i, ] <- ifelse(is.na(y_agg_enr[i, ]), 0, y_agg_enr[i, ])
    z_agg_ear[i, ] <- aggregate(z[c(id_analysis, id_impute_early)], by = list(x[c(id_analysis, id_impute_early)]), FUN = sum, drop = F)[, 2]
    z_agg_ear[i, ] <- ifelse(is.na(z_agg_ear[i, ]), 0, z_agg_ear[i, ])
    z_agg_obs[i, ] <- aggregate(z[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    z_agg_obs[i, ] <- ifelse(is.na(z_agg_obs[i, ]), 0, z_agg_obs[i, ])
    z_agg_enr[i, ] <- aggregate(z[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
    z_agg_enr[i, ] <- ifelse(is.na(z_agg_enr[i, ]), 0, z_agg_enr[i, ])

    # Transitions
    w01 <- aggregate(y[id_analysis] == 1 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    w01 <- ifelse(is.na(w01), 0, w01)
    w00 <- aggregate(y[id_analysis] == 0 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    w00 <- ifelse(is.na(w00), 0, w00)
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
      mod <- vb_mod(Xactivated, y_agg_obs[i, activated_arms], n = n_agg_obs[i, activated_arms])
      b_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    }

    m_draws <- b_draws %*% t(Xactivated)
    p_draws <- plogis(m_draws)
    eff_draws <- sweep(p_draws[, -1], 1, p_draws[, 1])

    p_supr[i, activated_arms[-1] - 1] <- matrixStats::colMeans2(eff_draws < 0)
    is_supr[i, activated_arms[-1] - 1] <- p_supr[i, activated_arms[-1] - 1] > supr_thres
    p_best[i, activated_arms] <- prob_best(m_draws, minimum = TRUE)
    is_best[i, activated_arms] <- p_best[i, activated_arms] > g(i)
    is_current_best[i, ] <- p_best[i, ] == max(p_best[i, ])

    if(drop_best) {
      is_active[i, activated_arms] <- p_best[i, activated_arms] >= f(i)
    } else {
      is_active[i, activated_arms][-1] <- p_supr[i, activated_arms[-1] - 1] >= f(i)
      is_active[i, 1] <- !any(is_supr[i, ])
    }
    if(perm_drop & i > 1) {
      is_active[i, activated_arms] <- is_active[i, activated_arms] & is_active[i - 1, activated_arms]
    }

    # Include control group in RAR or not?
    if(rar_control) {
      p_rand[i + 1, ] <- brar_all(p_best[i, ], is_active[i, ], h(i))
    } else {
      p_rand[i + 1, ] <- fix_ctrl_brar(p_best[i, ], is_active[i, ], h(i))
    }

    p_best_active[i, is_active[i, ]] <- prob_best(m_draws[, is_active[i, activated_arms], drop = F], minimum = TRUE)

    # Parameter summaries
    hdival_p <- HDInterval::hdi(p_draws)
    hdival_b <- HDInterval::hdi(b_draws)
    beta_mu[i, activated_pars] <- matrixStats::colMeans2(b_draws)
    beta_lb[i, activated_pars] <- hdival_b[1, ]
    beta_ub[i, activated_pars] <- hdival_b[2, ]
    beta_eff[i, activated_pars] <- matrixStats::colMeans2(b_draws < 0)
    arm_mu_p[i, activated_arms] <- matrixStats::colMeans2(p_draws)
    arm_var_p[i, activated_arms] <- matrixStats::colVars(p_draws)
    arm_mu_lb[i, activated_arms] <- hdival_p[1, ]
    arm_mu_ub[i, activated_arms] <- hdival_p[2, ]
    eff_mu_p[i, activated_arms[-1] - 1] <- matrixStats::colMeans2(eff_draws)
    eff_var_p[i, activated_arms[-1] - 1] <- matrixStats::colVars(eff_draws)

  }
  p_rand <- p_rand[-1, , drop = F]

  interim_quant <- loo::nlist(
    n_agg_enr, y_agg_enr, z_agg_enr,
    n_agg_obs, y_agg_obs, z_agg_obs,
    n_agg_ear, z_agg_ear,
    p_trans,
    p_supr, is_supr,
    p_best, is_best,
    p_best_active,
    p_rand, is_active, is_current_best,
    arm_mu_p, eff_mu_p,
    arm_mu_lb, arm_mu_ub,
    arm_var_p, eff_var_p
  )

  model_quant <- loo::nlist(
    beta_mu, beta_lb, beta_ub, beta_eff
  )

  trial_quant <- list(

  )

  return(
    loo::nlist(trial_quant, interim_quant, model_quant)
  )
}
