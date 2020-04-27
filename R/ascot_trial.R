#' Run ASCOT trial simulation
#'
#'
#' @export
ascot_trial <- function(
  b,
  n_seq,
  M0 = rep(0, 5),
  S0 = diag(c(2^2, 2^2, 2^2, 0.5^2, 2^2), 5),
  delta = 0.2,
  effective_thres = 0.975,
  equivalent_thres = 0.9,
  ineffective_thres = 0.1,
  mod_effective_thres = 0.05,
  best_thres = 0.95,
  rar_control = FALSE,
  rar_best = TRUE,
  perm_drop = TRUE,
  use_mwud = FALSE,
  use_optimal = FALSE,
  early_t = 7,
  h = 0.5,
  trans_prob = 0.75,
  mc_draw = 2e4,
  return_all = FALSE,
  ...
) {

  # Design matrix for domains
  D <- expand.grid(HCQ = c(0, 1), LR = c(0, 1), CP = c(0, 1))
  X <- model.matrix( ~ HCQ * LR + CP, data = D)
  colnames(X)[1] <- "SoC"
  rownames(X) <- sapply(1:nrow(X), function(i) paste(colnames(X)[1:4][X[i, 1:4] == 1], collapse = " + "))
  XA <- X[, c(1,2,3,5)]
  XB <- X[, 4, drop = F]
  X <- cbind(XA, XB)
  p <- plogis(X %*% b)[, 1]
  n_arms <- nrow(X)
  n_pars <- ncol(X)
  Xarm <- rbind(cbind(1, rbind(0, diag(1, 3)), 0), cbind(1, rbind(0, diag(1, 3)), 1))
  colnames(Xarm) <- colnames(X)
  rownames(Xarm) <- rownames(X)
  if(use_optimal) {
    ratio <- c(sqrt(n_arms - 1), rep(1, n_arms - 1))
    ratio <- ratio / sum(ratio)
  } else {
    ratio <- rep(1 / n_arms, n_arms)
  }


  # Data
  full_dat <- gen_potential_outcomes(max(n_seq), p, ...)
  tdat <- full_dat[, 1:2]
  dat <- full_dat[, -(1:2)]

  # Interim schedule
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_seq_aug <- c(0, n_seq)
  n_new <- diff(n_seq_aug)
  x <- factor(numeric(n_max), levels = 1:n_arms)
  y <- rep(NA_real_, n_max)
  z <- rep(NA_real_, n_max)
  n_enrolled <- 0

  # Output labels
  arm_names <- rownames(X)
  par_names <- colnames(X)
  arm_dm <- list("interim" = 1:n_int, "arm" = arm_names)
  arm_dm2 <- list("interim" = 1:n_int, "arm" = arm_names[-1])
  par_dm <- list("interim" = 1:n_int, "treatment" = par_names)
  par_dm2 <- list("interim" = 1:n_int, "treatment" = par_names[-1])

  # Data aggregation storage
  n_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  # z_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  n_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  # z_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  # n_agg_ear <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  # z_agg_ear <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  # p_trans <- matrix(0, n_int, n_arms, dimnames = arm_dm)

  # Model parameter storage
  is_par_active <- matrix(TRUE, n_int, n_pars, dimnames = par_dm)
  is_par_effective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_par_equivalent <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_par_ineffective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  beta_mu <- matrix(0, n_int, n_pars, dimnames = par_dm)
  beta_var <- matrix(0, n_int, n_pars, dimnames = par_dm)
  beta_lb <- matrix(0, n_int, n_pars, dimnames = par_dm)
  beta_ub <- matrix(0, n_int, n_pars, dimnames = par_dm)
  beta_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  beta_mod_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  beta_equ <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)

  is_trt_active <- matrix(TRUE, n_int, n_pars, dimnames = par_dm)
  is_trt_effective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_trt_equivalent <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_trt_ineffective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_trt_not_mod_effective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  trt_mu <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_lb <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_ub <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_mod_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_equ <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)

  p_best_antiviral <- matrix(0, n_int, 3, dimnames = list("interim" = 1:n_int, "arm" = par_names[2:4]))

  # Arm storage
  p_superior <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  p_equivalent <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  p_best <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  p_best_active <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  p_best_treat <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  p_rand <- matrix(ratio, n_int + 1, n_arms, dimnames = list("interim" = 0:n_int, "arm" = arm_names))

  is_arm_superior <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2)
  is_arm_inferior <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2)
  is_arm_equivalent <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2)
  is_arm_active <- matrix(TRUE, n_int, n_arms, dimnames = arm_dm)
  is_arm_best <- matrix(FALSE, n_int, n_arms, dimnames = arm_dm)

  arm_mu <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_var <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_lb <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_ub <- matrix(0, n_int, n_arms, dimnames = arm_dm)

  # futile <- FALSE
  stopped <- FALSE
  final <- FALSE

  # Loop interims
  for(i in 1:n_int) {

    # What gets analysed
    if(!stopped) {
      t_ref <- tdat[, 2][n_seq[i]]
      n_enrol <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
      id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
      n_enrolled <- n_enrolled + n_enrol

      if(use_mwud) {
        new_x <- factor(mass_weighted_urn_design(p_rand[i, ], n_enrol, alpha = 5)$trt, levels = 1:n_arms)
      } else {
        new_x <- factor(sample.int(n_arms, n_enrol, prob = p_rand[i, ], replace = T), levels = 1:n_arms)
      }

      new_y <- dat[cbind(id_enrolled, new_x)]
      new_z <- (new_y == 1) * rbinom(n_enrol, 1, trans_prob)
      new_y_agg <- aggregate(new_y, list(new_x), sum, drop = F)[, 2]
      new_n_agg <- table(new_x, dnn = NULL)

      x[id_enrolled] <- new_x
      y[id_enrolled] <- new_y
      z[id_enrolled] <- new_z
      id_analysis <- which(tdat[, 2] <= t_ref)
    } else {
      final <- TRUE
      id_analysis <- 1:n_enrolled
    }

    # What is imputed?
    # id_impute_early <- which(tdat[, 2] - t_ref > 0 & tdat[, 2] - t_ref <= early_t)
    # id_impute_late  <- which(tdat[, 1] <= t_ref & tdat[, 2] - t_ref > early_t)

    # Do the aggregations
    n_agg_obs[i, ] <- table(x[id_analysis], dnn = NULL)
    n_agg_enr[i, ] <- table(x[1:n_enrolled], dnn = NULL)
    # n_agg_ear[i, ] <- table(x[c(id_analysis, id_impute_early)], dnn = NULL)
    y_agg_obs[i, ] <- aggregate(y[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    y_agg_enr[i, ] <- aggregate(y[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]
    # z_agg_ear[i, ] <- aggregate(z[c(id_analysis, id_impute_early)], by = list(x[c(id_analysis, id_impute_early)]), FUN = sum, drop = F)[, 2]
    # z_agg_obs[i, ] <- aggregate(z[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    # z_agg_enr[i, ] <- aggregate(z[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]

    # Do the transitions
    # w01 <- aggregate(y[id_analysis] == 1 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    # w00 <- aggregate(y[id_analysis] == 0 & z[id_analysis] == 0, by = list(x[id_analysis]), sum, drop = F)[, 2]
    # p_trans[i, ] <- (1 + w01) / (2 + w01 + w00)

    # Fit the model
    mod <- vb_mod(X, y = y_agg_obs[i, ], n = n_agg_obs[i, ], M0 = M0, S0 = S0)
    b_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    m_draws <- b_draws %*% t(X)
    p_draws <- plogis(m_draws)
    eff_draws <- b_draws[, -1] %*% t(X[-1, -1] )

    # Parameter summaries
    hdival <- HDInterval::hdi(b_draws)
    beta_mu[i, ] <- matrixStats::colMeans2(b_draws)
    beta_var[i, ] <- diag(mod$Sigma)
    beta_lb[i, ] <- hdival[1, ]
    beta_ub[i, ] <- hdival[2, ]
    beta_eff[i, ] <- matrixStats::colMeans2(b_draws[, -1] < 0)
    beta_mod_eff[i, ] <- matrixStats::colMeans2(b_draws[, -1] < -delta)
    beta_equ[i, ] <- matrixStats::colMeans2(abs(b_draws[, -1]) < delta)

    hdival <- HDInterval::hdi(eff_draws[, 1:4])
    trt_mu[i, ] <- matrixStats::colMeans2(eff_draws[, 1:4])
    trt_lb[i, ] <- hdival[1, ]
    trt_ub[i, ] <- hdival[2, ]
    trt_eff[i, ] <- matrixStats::colMeans2(eff_draws[, 1:4] < 0)
    trt_mod_eff[i, ] <- matrixStats::colMeans2(eff_draws[, 1:4] < -delta)
    trt_equ[i, ] <- matrixStats::colMeans2(abs(eff_draws[, 1:4]) < delta)

    p_best_antiviral[i, ] <- prob_best(eff_draws[, 1:3], minimum = T)

    hdival <- HDInterval::hdi(p_draws)
    arm_mu[i, ] <- matrixStats::colMeans2(p_draws)
    arm_var[i, ] <- matrixStats::colVars(p_draws)
    arm_lb[i, ] <- hdival[1, ]
    arm_ub[i, ] <- hdival[2, ]

    is_par_effective[i, ] <- beta_eff[i, ] > effective_thres
    is_par_equivalent[i, ] <- beta_equ[i, ] > equivalent_thres
    is_par_ineffective[i, ] <- beta_eff[i, ] < ineffective_thres
    is_par_active[i, 1] <- !any(is_par_effective[i, ])
    is_par_active[i, -1] <- !(is_par_ineffective)[i, ]

    is_trt_effective [i, ] <- trt_eff[i, ] > effective_thres
    is_trt_equivalent[i, ] <- trt_equ[i, ] > equivalent_thres
    is_trt_ineffective[i, ] <- trt_eff[i, ] < ineffective_thres
    is_trt_not_mod_effective[i, ] <- trt_mod_eff[i, ] < mod_effective_thres
    is_trt_active[i, 1] <- !any(is_trt_effective[i, ])
    is_trt_active[i, -1] <- !(is_trt_ineffective)[i, ]

    # If drop permanently, was the arm already inactive?
    if(perm_drop & i > 1) {
      is_trt_active[i, ] <- is_trt_active[i, ] & is_trt_active[i - 1, ]
      is_par_active[i, ] <- is_par_active[i, ] & is_par_active[i - 1, ]
    }

    # Arm summaries
    p_superior[i, ] <- matrixStats::colMeans2(eff_draws < 0)
    p_equivalent[i, ] <- matrixStats::colMeans2(abs(eff_draws) < delta)
    p_best[i, ] <- prob_best(m_draws, minimum = TRUE)
    p_best_treat[i, ] <- prob_best(m_draws[, -1], minimum = TRUE)
    is_arm_superior[i, ] <- p_superior[i, ] > effective_thres
    is_arm_equivalent[i, ] <- p_equivalent[i, ] > equivalent_thres
    is_arm_inferior[i, ] <- p_superior[i, ] < ineffective_thres
    is_arm_best[i, ] <- p_best[i, ] > best_thres
    is_arm_active[i, 1] <- is_trt_active[i, 1]
    is_arm_active[i, -1] <- !(rowSums(Xarm[-1, -1][, !is_trt_active[i, -1], drop = F]) > 0)

    # Include control group in RAR or not?
    if(rar_control) {
      p_rand[i + 1, ] <- brar_all(p_best[i, ], is_arm_active[i, ], h)
    } else {
      if(rar_best) {
        p_rand[i + 1, ] <- fix_ctrl_brar(p_best[i, ], is_arm_active[i, ], h)
      } else {
        if(use_optimal) {
          ratio <- c(sqrt(sum(is_arm_active[i, -1])), rep(1, sum(is_arm_active[i, -1])))
          ratio <- ratio / sum(ratio)
        } else {
          ratio <- 1 / sum(is_arm_active[i, ])
        }
        p_rand[i + 1, ] <- fix_ctrl_brar2(c(0, p_superior[i, ]), is_arm_active[i, ], h, ratio[1])
      }
    }

    # If everything or everything but SoC has been dropped, just stop
    if(!any(is_trt_active[i, -1])) stopped <- TRUE
    # If only one active and it's superior then stop
    # if(!any(is_trt_active[i, -1] & !is_trt_effective[i, ])) stopped <- TRUE

    # if(final) break
    if(stopped) break

  }
  p_rand <- p_rand[-1, ]

  drop_trt_at <- apply(is_trt_active, 2, function(x) which(x == FALSE)[1])

  trigger_effective_at <- apply(is_trt_effective, 2, function(x) which(x)[1])
  trigger_equivalent_at <- apply(is_trt_equivalent, 2, function(x) which(x)[1])
  trigger_ineffective_at <- apply(is_trt_ineffective, 2, function(x) which(x)[1])
  trigger_ineffective <- ifelse(is.na(trigger_ineffective_at), FALSE, ifelse(trigger_ineffective_at < n_int, TRUE, FALSE))
  trigger_equivalent <- ifelse(is.na(trigger_equivalent_at), FALSE, ifelse(trigger_equivalent_at < n_int, TRUE, FALSE))
  trigger_effective <- ifelse(is.na(trigger_effective_at), FALSE, ifelse(trigger_effective_at < n_int, TRUE, FALSE))
  trigger_futility <- trigger_equivalent | trigger_ineffective
  final_effective <- is_trt_effective[i, ]
  final_ineffective <- is_trt_ineffective[i, ]
  final_equivalent <- is_trt_equivalent[i, ]
  final_not_mod_effective <- is_trt_not_mod_effective[i, ]
  stop_early <- i < n_int
  stop_at <- i
  drop_soc_at <- findfirst(!is_trt_active[, 1])
  futility <- all(trigger_futility)
  success <-any(trigger_effective)
  effective <- any(final_effective)
  ineffective <- all(final_ineffective)


  if(return_all) indx <- 1:i else indx <- i

  n_agg_enr <- n_agg_enr[indx, , drop = F]
  y_agg_enr <- y_agg_enr[indx, , drop = F]
  n_agg_obs <- n_agg_obs[indx, , drop = F]
  y_agg_obs <- y_agg_obs[indx, , drop = F]
  p_rand <- p_rand[indx, , drop = F]
  p_superior <- p_superior[indx, , drop = F]
  p_equivalent <-  p_equivalent[indx, , drop = F]
  p_best <- p_best[indx, , drop = F]
  p_best_treat <- p_best_treat[indx, , drop = F]
  arm_mu <- arm_mu[indx, , drop = F]
  arm_var <- arm_var[indx, , drop = F]
  arm_lb <- arm_lb[indx, , drop = F]
  arm_ub <- arm_ub[indx, , drop = F]
  is_arm_superior <- is_arm_superior[indx, , drop = F]
  is_arm_equivalent <- is_arm_equivalent[indx, , drop = F]
  is_arm_inferior <- is_arm_inferior[indx, , drop = F]
  is_arm_best <- is_arm_best[indx, , drop = F]
  is_arm_active <- is_arm_active[indx, , drop = F]

  beta_mu <- beta_mu[indx, , drop = F]
  beta_var <- beta_var[indx, , drop = F]
  beta_lb <- beta_lb[indx, , drop = F]
  beta_ub <- beta_ub[indx, , drop = F]
  beta_eff <- beta_eff[indx, , drop = F]
  beta_equ <- beta_equ[indx, , drop = F]
  beta_mod_eff <- beta_mod_eff[indx, , drop = F]
  trt_mu <- trt_mu[indx, , drop = F]
  trt_lb <- trt_lb[indx, , drop = F]
  trt_ub <- trt_ub[indx, , drop = F]
  trt_eff <- trt_eff[indx, , drop = F]
  trt_equ <- trt_equ[indx, , drop = F]
  trt_mod_eff <- trt_mod_eff[indx, , drop = F]
  is_trt_effective <-is_trt_effective[indx, , drop = F]
  is_trt_equivalent <- is_trt_equivalent[indx, , drop = F]
  is_trt_ineffective <- is_trt_ineffective[indx, , drop = F]
  is_trt_not_mod_effective <- is_trt_not_mod_effective[indx, , drop = F]
  is_trt_active <- is_trt_active[indx, , drop = F]
  is_par_active <- is_par_active[indx, , drop = F]
  is_par_effective <- is_par_effective[indx, , drop = F]
  is_par_equivalent <- is_par_equivalent[indx, , drop = F]
  is_par_ineffective <- is_par_ineffective[indx, , drop = F]
  p_best_antiviral <- p_best_antiviral[indx, , drop = F]

  interim_quantities <- loo::nlist(
    n_agg_enr, y_agg_enr, #z_agg_enr,
    n_agg_obs, y_agg_obs, #z_agg_obs,
    # n_agg_ear, z_agg_ear, p_trans,
    p_rand, p_superior, p_equivalent, p_best, p_best_treat,
    arm_mu, arm_var, arm_lb, arm_ub,
    is_arm_superior, is_arm_equivalent,
    is_arm_inferior, is_arm_best, is_arm_active
  )

  model_quantities <- loo::nlist(
    beta_mu, beta_var, beta_lb, beta_ub, beta_eff, beta_equ, beta_mod_eff,
    trt_mu, trt_lb, trt_ub, trt_eff, trt_equ, trt_mod_eff,
    is_trt_effective, is_trt_equivalent, is_trt_ineffective, is_trt_active, is_trt_not_mod_effective,
    is_par_active, is_par_effective, is_par_equivalent, is_par_ineffective,
    p_best_antiviral
  )

  trial_quantities <- loo::nlist(
    drop_trt_at,
    trigger_effective_at, trigger_equivalent_at, trigger_ineffective_at,
    trigger_effective, trigger_equivalent, trigger_ineffective, trigger_futility,
    final_effective, final_ineffective, final_equivalent, final_not_mod_effective
  )

  result_quantities <- c(
    stop_early = stop_early,
    stop_at = stop_at,
    drop_soc_at = unname(drop_soc_at),
    futility = futility,
    success = success,
    effective = effective,
    ineffective = ineffective)

  return(loo::nlist(
    interim_quantities,
    model_quantities,
    trial_quantities,
    result_quantities
  ))
}

#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_interim_quantities <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["interim_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "arm", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "arm")) %>%
      dplyr::mutate(arm = forcats::fct_inorder(arm)) %>%
      dplyr::arrange(interim, arm)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_model_quantities <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["model_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "treatment", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "treatment")) %>%
      dplyr::mutate(treatment = forcats::fct_inorder(treatment)) %>%
      dplyr::arrange(interim, treatment)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_trial_quantities <- function(dat, ...) {
  dplyr::bind_rows(lapply(dat, function(i) {
    tq <- i[["trial_quantities"]]
    tibble::enframe(tq) %>%
      tidyr::unnest_wider(value) %>%
      tidyr::gather(treatment, value, -name) %>%
      tidyr::spread(name, value)
  }), .id = "trial")
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `lr_brar_trial` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_result_quantities <- function(dat, ...) {
  enframe(lapply(dat, function(x) x[["result_quantities"]]), name = "trial") %>%
    unnest_wider(value)
}

