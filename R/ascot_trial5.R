#' Run updated ASCOT trial simulation
#'
#'@param b The true log-odds parameter
#'@param n_seq The sequence of interim analysis in terms of sample size
#'@param S0 The prior covariance matrix
#'@param delta The reference log-odds value for futility
#'@param effective_thres The threshold for declaring a treatment effective
#'@param ineffective_thres The threshold for declaring harm
#'@param inferior_thres The threshold for declaring inferiority, scaled by number of arms in domain
#'@param equivalent_thres The threshold for declaring a treatment equivalent to SoC
#'@param futility_thres The threshold for declaring a treatment futile
#'@param superior_thres The threshold for declaring superiority, scaled by number of arms in domain
#'@param scale_sup_thres Scale the superiorirty threshold by number of arms in domain?
#'@param rar_control Indicate if control should be included in RAR or fixed to 1/(num active arms)
#'@param rar_best Indicate if RAR should be proportional to Pr(best) or Pr(effective)
#'@param rar_scale Scaling factor for RAR
#'@param rar_use_ss Use sample size in RAR
#'@param use_optimal Use optimal allocation ratios
#'@param early_t Timing of early event
#'@param return_all Return all interim results
#'@param min_rar The minimum allocation to an arm in each domain (e.g. if Pr(alloc=1)<min_rar then set Pr(alloc=1)=min_rar)
#'@param perm_drop Permanently drop arms or only suspend them
#'@param use_mwud Use Mass-Weighted Urn for randomisation or just simple randomisation
#'@param mc_draw The number of Monte Carlo draws to use for computing posterior quantities
#'@param ... Other arguments
#'
#' @export
ascot_trial5 <- function(
  b,
  n_seq,
  M0 = rep(0, 7),
  S0 = diag(c(100, rep(2.5^2, 6))),
  delta = log(1.1),
  effective_thres = 0.99,
  futility_thres = 0.95,
  inferior_thres = 0.01,
  superior_thres = 0.99,
  min_any = TRUE,
  rar_scale = 0.5,
  rar_use_ss = TRUE,
  min_rar = rep(1/5,3),
  perm_drop = TRUE,
  use_mwud = TRUE,
  early_t = 0,
  mc_draw = 2e4,
  return_all = FALSE,
  ...
) {

  if(length(b) != 7) stop("wrong number of parameters in b. Should be 7")

  # DESIGN #
  #--------#

  # Arm indicator designs, note that "x0" indicates nothing (or SoC) from domain x
  # I.e. does the regimen involve this treatment.
  XIa <- diag(1, 4)
  colnames(XIa) <- rownames(XIa) <- c("a0", "a1", "a2", "a3")
  XIb <- diag(1, 3)
  colnames(XIb) <- rownames(XIb) <- c("b0", "b1", "b2")
  XIc <- diag(1, 2)
  colnames(XIc) <- rownames(XIc) <- c("c0", "c1")
  XI <- tidyr::expand_grid(as.data.frame(XIa), as.data.frame(XIb), as.data.frame(XIc))
  XI <- as.matrix(XI)
  rownames(XI) <- sapply(1:nrow(XI), function(i)
    paste(gsub(":", "", colnames(XI)[XI[i, ] == 1]), collapse = ""))


  # Design matrix for domains
  Xa <- eigen(diag(1,4) - 1/4)$vectors[, 1:3]
  colnames(Xa) <- c("a1", "a2", "a3")
  rownames(Xa) <- c("a0", "a1", "a2", "a3")
  Xb <- eigen(diag(1,3) - 1/3)$vectors[, 1:2]
  colnames(Xb) <- c("b1", "b2")
  rownames(Xb) <- c("b0", "b1", "b2")
  Xc <- eigen(diag(1,2) - 1/2)$vectors[, 1, drop = F]
  colnames(Xc) <- c("c1")
  rownames(Xc) <- c("c0", "c1")
  X <- tidyr::expand_grid(as.data.frame(Xa), as.data.frame(Xb), as.data.frame(Xc))
  X <- cbind("a0b0c0" = 1, as.matrix(X))
  rownames(X) <- rownames(XI)

  # Arm to domain map
  XM <- tidyr::expand_grid(a = 1:4, b = 1:3, c = 1:2)

  # Initialise allocation ratios
  ratio_dom <- sapply(list(A = Xa, B = Xb, C = Xc), function(j) rep(1 / (dim(j)[1]), dim(j)[1]))
  ratio <- ratio_dom$A[XM$a]*ratio_dom$B[XM$b]*ratio_dom$C[XM$c]


  #-----------------------------------

  # STORAGE #
  #---------#

  # Data
  n_arms <- nrow(X)
  n_pars <- ncol(X)

  p <- plogis(X %*% b)[, 1]
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
  n_enrolled <- 0

  # Output labels
  arm_names <- rownames(XI)
  trt_names <- colnames(XI)
  par_names <- colnames(X)
  act_names <- trt_names[!grepl("0", trt_names)]

  arm_dm <- list("interim" = 1:n_int, "arm" = arm_names)
  arm_dm2 <- list("interim" = 1:n_int, "arm" = arm_names[-1])
  trt_dm <- list("interim" = 1:n_int, "treatment" = trt_names)
  trt_dm2 <- list("interim" = 1:n_int, "treatment" = trt_names[-1])
  act_dm <- list("interim" = 1:n_int, "treatment" = act_names)
  par_dm <- list("interim" = 1:n_int, "parameter" = par_names)
  par_dm2 <- list("interim" = 1:n_int, "parameter" = par_names[-1])

  # Data aggregation storage
  n_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  n_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  n_agg_enr_trt <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)
  y_agg_enr_trt <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)
  n_agg_obs_trt <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)
  y_agg_obs_trt <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)

  # Allocation probabilities by domain and overall
  p_rand_a <- matrix(ratio_dom$A, n_int + 1, nrow(Xa), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xa)))
  p_rand_b <- matrix(ratio_dom$B, n_int + 1, nrow(Xb), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xb)))
  p_rand_c <- matrix(ratio_dom$C, n_int + 1, nrow(Xc), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xc)))
  p_rand_trt <- matrix(c(ratio_dom$A, ratio_dom$B, ratio_dom$C),
                       n_int + 1, ncol(XI), byrow = T,
                       dimnames = list("interim" = 0:n_int,
                                       "treatment" = colnames(XI)))
  p_rand_trt_raw <- p_rand_trt
  p_rand <- matrix(ratio, n_int + 1, n_arms, dimnames = list("interim" = 0:n_int, "arm" = arm_names))

  # Model parameter storage
  par_mu  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_var <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_lb  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_ub  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2) # Pr(b<0)
  par_fut <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2) # Pr(b<d)

  # Treatment storage
  trt_mu  <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_var <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_lb  <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_ub  <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_eff <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_fut <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_bes <- matrix(0, n_int, n_pars - 1, dimnames = act_dm)
  trt_bes_all <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)
  trt_bes_act <- trt_bes_all
  trt_in_best <- matrix(0, n_int, ncol(XI), dimnames = trt_dm)
  trt_in_best_act <- trt_in_best

  is_trt_active <- matrix(TRUE, n_int, ncol(XI), dimnames = trt_dm)
  is_trt_eff    <- matrix(FALSE, n_int, n_pars - 1, dimnames = act_dm)
  is_trt_fut    <- matrix(FALSE, n_int, n_pars - 1, dimnames = act_dm) # Pr(-d<sum(b)<d)
  is_trt_sup    <- matrix(FALSE, n_int, ncol(XI), dimnames = trt_dm)
  is_trt_inf    <- matrix(FALSE, n_int, ncol(XI), dimnames = trt_dm)

  # Regimen storage
  # reg_mu  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_var <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_lb  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_ub  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_eff <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_fut <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  # reg_equ <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_bes <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  reg_bes_act <- reg_bes
  is_reg_active <- matrix(TRUE, n_int, n_arms, dimnames = arm_dm)

  arm_mu <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_var <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_lb <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  arm_ub <- matrix(0, n_int, n_arms, dimnames = arm_dm)

  # futile <- FALSE
  stopped <- FALSE
  final <- FALSE

  #--------------------------------


  # INTERIMS #
  #----------#

  for(i in 1:n_int) {

    if(!stopped) {

      # What gets analysed
      t_ref       <- tdat[, 2][n_seq[i]]
      n_enrol     <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
      id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
      n_enrolled  <- n_enrolled + n_enrol

      # If not stopped generate new data
      if(use_mwud) {
        new_x <- factor(mass_weighted_urn_design(p_rand[i, ], n_enrol, alpha = 5)$trt, levels = 1:n_arms)
      } else {
        new_x <- factor(sample.int(n_arms, n_enrol, prob = p_rand[i, ], replace = T), levels = 1:n_arms)
      }

      new_y <- dat[cbind(id_enrolled, new_x)]
      new_y_agg <- aggregate(new_y, list(new_x), sum, drop = F)[, 2]
      new_n_agg <- table(new_x, dnn = NULL)

      x[id_enrolled] <- new_x
      y[id_enrolled] <- new_y
      id_analysis <- which(tdat[, 2] <= t_ref)

    } else {

      final <- TRUE
      id_analysis <- 1:n_enrolled

    }

    # Do the aggregations
    n_agg_obs[i, ] <- table(x[id_analysis], dnn = NULL)
    n_agg_enr[i, ] <- table(x[1:n_enrolled], dnn = NULL)
    y_agg_obs[i, ] <- aggregate(y[id_analysis], by = list(x[id_analysis]), FUN = sum, drop = F)[, 2]
    y_agg_enr[i, ] <- aggregate(y[1:n_enrolled], by = list(x[1:n_enrolled]), FUN = sum, drop = F)[, 2]

    n_agg_obs_trt[i, ] <- matrixStats::colSums2(XI[as.integer(x[id_analysis]), ])
    n_agg_enr_trt[i, ] <- matrixStats::colSums2(XI[as.integer(x[1:n_enrolled]), ])
    y_agg_obs_trt[i, 1:4] <- aggregate(y[id_analysis], by = XM[as.integer(x[id_analysis]), 1], FUN = sum, drop = F)[, 2]
    y_agg_obs_trt[i, 5:7] <- aggregate(y[id_analysis], by = XM[as.integer(x[id_analysis]), 2], FUN = sum, drop = F)[, 2]
    y_agg_obs_trt[i, 8:9] <- aggregate(y[id_analysis], by = XM[as.integer(x[id_analysis]), 3], FUN = sum, drop = F)[, 2]
    y_agg_enr_trt[i, 1:4] <- aggregate(y[1:n_enrolled], by = XM[as.integer(x[1:n_enrolled]), 1], FUN = sum, drop = F)[, 2]
    y_agg_enr_trt[i, 5:7] <- aggregate(y[1:n_enrolled], by = XM[as.integer(x[1:n_enrolled]), 2], FUN = sum, drop = F)[, 2]
    y_agg_enr_trt[i, 8:9] <- aggregate(y[1:n_enrolled], by = XM[as.integer(x[1:n_enrolled]), 3], FUN = sum, drop = F)[, 2]

    #
    if(i == 1) {
      last_active_all <- is_trt_active[i, ]
      last_active_trt <- is_trt_active[i, !grepl("0", colnames(is_trt_active))]
    } else {
      last_active_all <- is_trt_active[i-1, ]
      last_active_trt <- is_trt_active[i-1, !grepl("0", colnames(is_trt_active))]
    }

    # Fit the model
    mod <- vb_mod(X, y = y_agg_obs[i, ], n = n_agg_obs[i, ], M0 = M0, S0 = S0)
    par_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    colnames(par_draws) <- colnames(X)

    # Linear predictor
    eta_draws <- par_draws %*% t(X)
    eta_draws_a <- par_draws[, grepl("a", colnames(par_draws))] %*% t(cbind(1, Xa))
    eta_draws_b <- par_draws[, grepl("b", colnames(par_draws))] %*% t(cbind(1, Xb))
    eta_draws_c <- par_draws[, grepl("c", colnames(par_draws)), drop = F] %*% t(cbind(1, Xc))

    # Probability of response
    p_draws <- plogis(eta_draws)

    # Treatment Effects
    trt_draws_a <- par_draws[, grepl("a[1-9]", colnames(par_draws))] %*% t(Xa[-1, ])
    trt_draws_b <- par_draws[, grepl("b[1-9]", colnames(par_draws))] %*% t(Xb[-1, ])
    trt_draws_c <- par_draws[, grepl("c[1-9]", colnames(par_draws)), drop = F] %*% t(Xc[-1, , drop = F])
    trt_draws   <- cbind(trt_draws_a, trt_draws_b, trt_draws_c)

    # Regimen Effects
    # reg_draws <- par_draws[, -1] %*% t(X[, -1][-1, ])

    # Parameter summaries
    hdival <- HDInterval::hdi(par_draws)
    par_mu[i, ]  <- matrixStats::colMeans2(par_draws)
    par_var[i, ] <- diag(mod$Sigma)
    par_lb[i, ]  <- hdival[1, ]
    par_ub[i, ]  <- hdival[2, ]

    # Arm summaries
    hdival       <- HDInterval::hdi(p_draws)
    arm_mu[i, ]  <- matrixStats::colMeans2(p_draws)
    arm_var[i, ] <- matrixStats::colVars(p_draws)
    arm_lb[i, ]  <- hdival[1, ]
    arm_ub[i, ]  <- hdival[2, ]

    # Regimen summaries
    reg_bes[i, ] <- prob_best(eta_draws, minimum = T)

    # Treatment summaries
    hdival       <- HDInterval::hdi(trt_draws)
    trt_mu[i, ]  <- matrixStats::colMeans2(trt_draws)
    trt_lb[i, ]  <- hdival[1, ]
    trt_ub[i, ]  <- hdival[2, ]
    trt_eff[i, ] <- matrixStats::colMeans2(trt_draws < 0)
    trt_fut[i, ] <- matrixStats::colMeans2(trt_draws > -delta)
    trt_bes[i, ] <- c(
      prob_best(trt_draws[, grepl("a", colnames(trt_draws))], minimum = T),
      prob_best(trt_draws[, grepl("b", colnames(trt_draws))], minimum = T),
      prob_best(trt_draws[, grepl("c", colnames(trt_draws)), drop = F], minimum = T)
    )
    trt_bes_all[i, ] <- c(
      prob_best(eta_draws_a, minimum = T),
      prob_best(eta_draws_b, minimum = T),
      prob_best(eta_draws_c, minimum = T)
    )
    # Best among active
    if(i == 1) {
      trt_bes_act[i, ] <- trt_bes_all[i, ]
      reg_bes_act[i, ] <- reg_bes[i, ]
    } else {
      trt_bes_act[i, is_trt_active[i-1, ]] <- c(
        prob_best(eta_draws_a[, is_trt_active[i-1, grepl("a", colnames(is_trt_active))], drop = F], minimum = T),
        prob_best(eta_draws_b[, is_trt_active[i-1, grepl("b", colnames(is_trt_active))], drop = F], minimum = T),
        prob_best(eta_draws_c[, is_trt_active[i-1, grepl("c", colnames(is_trt_active))], drop = F], minimum = T)
      )
      reg_bes_act[i,  is_reg_active[i-1, ]] <- prob_best(eta_draws[,  is_reg_active[i-1, ], drop = F], minimum = T)
    }

    # Calc from prob(regimen is best)
    trt_in_best[i, ]     <- reg_bes[i, ] %*% XI
    trt_in_best_act[i, ] <- reg_bes_act[i, ] %*% XI

    # Treatment triggers
    # - decision can only be made while treatment is active
    # - once declared effective, can't declare futile
    # - once declared futile, can't declare effective
    is_trt_eff[i, ] <- (trt_eff[i, ] * last_active_trt) > effective_thres   # Better than SoC (only if still active)
    is_trt_fut[i, ] <- (trt_fut[i, ] * last_active_trt) > futility_thres    # Insufficiently better than SoC (only if still active)

    # Active treatments in domain
    # Refer to model document for the decision triggers, but basically
    # - If something superior to nothing, drop nothing
    # - If something equivalent to nothing drop the something
    # - If something superior drop everything else
    # - If something inferior drop the something
    # - If something futile (inadequate) drop the something
    # - If SoC (nothing) has been dropped, but subsequently every treatment also dropped, then reactivate SoC
    for(dom in c("a", "b", "c")) {
      idx1 <- grep(dom, colnames(is_trt_active)) # Interventions plus SoC
      idx2 <- grep(dom, colnames(is_trt_eff))    # Interventions excl SoC
      if(i == 1) {
        idx3 <- idx1
      } else {
        idx3 <- idx1[last_active_all[idx1]]
      }
      nact <- length(idx3)
      narm <- length(idx1)
      sup_ref <- superior_thres
      is_trt_sup[i, idx1] <- (trt_in_best_act[i, idx1] > sup_ref)  # Superior (best in domain)
      if(nact > 1) {
        is_trt_inf[i, idx3] <- (trt_in_best_act[i, idx3]) < inferior_thres / (nact - 1)  # Inferior (not best in domain)
      }
      is_trt_active[i, idx1[-1]] <- !(is_trt_fut[i, idx2] | is_trt_inf[i, idx1[-1]])
      is_trt_active[i, idx1[1]] <- !any(is_trt_eff[i, idx2] | is_trt_inf[i, idx1[1]])
      if(any(is_trt_sup[i, idx1])) {
        bidx <- which(is_trt_sup[i, idx1])
        is_trt_active[i, idx1][bidx] <- TRUE
        is_trt_active[i, idx1][-bidx] <- FALSE
      }
    }

    # If drop permanently, was the arm already inactive?
    if(perm_drop & i > 1) {
      is_trt_active[i, ] <- is_trt_active[i, ] & is_trt_active[i - 1, ]
    }
    is_reg_active[i, ] <- is_trt_active[i, ] %*% t(XI) == 3

    # Include control group in RAR or not?
    if(min_any) {
      if(rar_use_ss) {
        rar_quant <- reg_bes_act[i, ] / (n_agg_obs[i, ] + 1)
      } else {
        rar_quant <- reg_bes_act[i, ]
      }
      p_rand[i + 1, ] <- brar_fix(rar_quant, is_reg_active[i, ], rar_scale, 0)
      p_rand_trt[i + 1, ] <- p_rand[i + 1, ] %*% XI
      p_rand_trt_raw[i + 1, ] <- p_rand_trt[i + 1, ]
      p_rand[i + 1, ] <- correct_min_rar_any(p_rand[i + 1, ], min_rar, is_trt_active[i, ],  XI)
      p_rand_trt[i + 1, ] <- p_rand[i + 1, ] %*% XI
    } else {
      if(rar_use_ss) {
        rar_quant <- reg_bes_act[i, ] / (n_agg_obs[i, ] + 1)
      } else {
        rar_quant <- reg_bes_act[i, ]
      }
      p_rand[i + 1, ] <- brar_fix(rar_quant, is_reg_active[i, ], rar_scale, 0)
      p_rand_trt[i + 1, ] <- p_rand[i + 1, ] %*% XI
      min_rar <- pmin(
        1/3,
        c(1 / sum(is_trt_active[i, grepl("a", colnames(is_trt_active))]),
          1 / sum(is_trt_active[i, grepl("b", colnames(is_trt_active))]),
          1 / sum(is_trt_active[i, grepl("c", colnames(is_trt_active))])))
      p_rand_trt_raw[i + 1, ] <- p_rand_trt[i + 1, ]
      p_rand[i + 1, ] <- correct_min_rar_soc(p_rand[i + 1, ], min_rar, is_trt_active[i, grepl("0", colnames(is_trt_active))],  XI)
      p_rand_trt[i + 1, ] <- p_rand[i + 1, ] %*% XI
    }

    if(!all.equal(sum(p_rand_trt[i + 1, ]), 3)) {
      return(p_rand_trt[i + 1, ])
    }

    # Need to check here if SoC has been reactivated within domain
    # (e.g. because futility of all active treatment)
    # for(dom in c("a", "b", "c")) {
    #   idx <- grep(dom, colnames(is_trt_active))
    #   if(sum(is_trt_active[i, idx]) == 0) {
    #     is_trt_active[i, idx][1] <- 1
    #     p_rand_trt[i + 1, idx][1] <- 1
    #   }
    # }
    # is_reg_active[i, ] <- is_trt_active[i, ] %*% t(XI) == 3
  }
  p_rand <- p_rand[-1, ]
  p_rand_trt <- p_rand_trt[-1, ]
  p_rand_trt_raw <- p_rand_trt_raw[-1, ]

  drop_trt_at <- apply(is_trt_active, 2, function(x) which(x == FALSE)[1])

  # Summarise when triggers first occurred and whether they occurred at all
  trig_eff_at   <- apply(is_trt_eff, 2, findfirst)
  trig_fut_at   <- apply(is_trt_fut, 2, findfirst)
  trig_sup_at   <- apply(is_trt_sup, 2, findfirst)
  trig_inf_at   <- apply(is_trt_inf, 2, findfirst)
  trig_drop_at  <- apply(!is_trt_active, 2, findfirst)

  trig_eff   <- !is.na(trig_eff_at)
  trig_fut   <- !is.na(trig_fut_at)
  trig_sup   <- !is.na(trig_sup_at)
  trig_inf   <- !is.na(trig_inf_at)
  trig_drop  <- !is.na(trig_drop_at)

  final_eff    <- is_trt_eff[i, ]
  final_fut    <- is_trt_fut[i, ]
  final_sup    <- is_trt_sup[i, ]
  final_inf    <- is_trt_inf[i, ]
  final_active <- is_trt_active[i, ]

  stop_early <- stopped
  stop_at <- ifelse(stopped, i - 1, NA)

  final_results <- list(

    arm_quantities = list(
      n_agg_enr = n_agg_enr[i, , drop = F],
      y_agg_enr = y_agg_enr[i, , drop = F],
      n_agg_obs = n_agg_obs[i, , drop = F],
      y_agg_obs = y_agg_obs[i, , drop = F],
      p_rand = p_rand[i, , drop = F],
      reg_bes = reg_bes[i, , drop = F],
      reg_bes_act = reg_bes_act[i, , drop = F],
      arm_mu = arm_mu[i, , drop = F],
      arm_var = arm_var[i, , drop = F],
      arm_lb = arm_lb[i, , drop = F],
      arm_ub = arm_ub[i, , drop = F]
    ),

    par_quantities = list(
      par_mu = par_mu[i, , drop = F],
      par_var = par_var[i, , drop = F],
      par_lb = par_lb[i, , drop = F],
      par_ub = par_ub[i, , drop = F]
    ),

    trt_quantities = list(
      p_rand_trt = p_rand_trt[i, , drop = F],
      p_rand_trt_raw = p_rand_trt_raw[i, , drop = F],
      n_agg_enr_trt = n_agg_enr_trt[i, , drop = F],
      y_agg_enr_trt = y_agg_enr_trt[i, , drop = F],
      n_agg_obs_trt = n_agg_enr_trt[i, , drop = F],
      y_agg_obs_trt = y_agg_enr_trt[i, , drop = F],
      trt_mu = trt_mu[i, , drop = F],
      trt_lb = trt_lb[i, , drop = F],
      trt_ub = trt_ub[i, , drop = F],
      trt_eff = trt_eff[i, , drop = F],
      trt_fut = trt_fut[i, , drop = F],
      trt_bes = trt_bes[i, , drop = F],
      trt_bes_all = trt_bes_all[i, , drop = F],
      trt_bes_act = trt_bes_act[i, , drop = F],
      trt_in_best = trt_in_best[i, , drop = F],
      trt_in_best_act = trt_in_best_act[i, , drop = F],
      is_trt_eff = is_trt_eff[i, , drop = F],
      is_trt_fut = is_trt_fut[i, , drop = F],
      is_trt_active = is_trt_active[i, , drop = F],
      is_trt_sup = is_trt_sup[i, , drop = F],
      is_trt_inf = is_trt_inf[i, , drop = F]
    ),

    trial_quantities = loo::nlist(
      trig_drop_at,
      trig_eff_at,
      trig_fut_at,
      trig_sup_at,
      trig_inf_at,
      trig_drop,
      trig_eff,
      trig_fut,
      trig_sup,
      trig_inf,
      final_active,
      final_eff,
      final_fut,
      final_sup,
      final_inf
    ),

    result_quantities = c(
      stop_early = stop_early,
      stop_at = stop_at)
  )

  if(return_all) indx <- 1:(i-1) else indx <- i-1

  n_agg_enr <- n_agg_enr[indx, , drop = F]
  y_agg_enr <- y_agg_enr[indx, , drop = F]
  n_agg_obs <- n_agg_obs[indx, , drop = F]
  y_agg_obs <- y_agg_obs[indx, , drop = F]
  n_agg_enr_trt <- n_agg_enr_trt[indx, , drop = F]
  y_agg_enr_trt <- y_agg_enr_trt[indx, , drop = F]
  n_agg_obs_trt <- n_agg_obs_trt[indx, , drop = F]
  y_agg_obs_trt <- y_agg_obs_trt[indx, , drop = F]
  p_rand <- p_rand[indx, , drop = F]
  p_rand_trt_raw <- p_rand_trt_raw[indx, , drop = F]
  reg_bes <- reg_bes[indx, , drop = F]
  reg_bes_act <- reg_bes_act[indx, , drop = F]
  arm_mu <- arm_mu[indx, , drop = F]
  arm_var <- arm_var[indx, , drop = F]
  arm_lb <- arm_lb[indx, , drop = F]
  arm_ub <- arm_ub[indx, , drop = F]

  par_mu <- par_mu[indx, , drop = F]
  par_var <- par_var[indx, , drop = F]
  par_lb <- par_lb[indx, , drop = F]
  par_ub <- par_ub[indx, , drop = F]

  trt_mu <- trt_mu[indx, , drop = F]
  trt_lb <- trt_lb[indx, , drop = F]
  trt_ub <- trt_ub[indx, , drop = F]
  trt_eff <- trt_eff[indx, , drop = F]
  trt_fut <- trt_fut[indx, , drop = F]
  trt_bes <- trt_bes[indx, , drop = F]
  trt_bes_all <- trt_bes_all[indx, , drop = F]
  trt_bes_act <- trt_bes_act[indx, , drop = F]
  trt_in_best <- trt_in_best[indx, , drop = F]
  trt_in_best_act <- trt_in_best_act[indx, , drop = F]
  is_trt_eff <-is_trt_eff[indx, , drop = F]
  is_trt_fut <- is_trt_fut[indx, , drop = F]
  is_trt_active <- is_trt_active[indx, , drop = F]
  is_trt_sup = is_trt_sup[indx, , drop = F]
  is_trt_inf = is_trt_inf[indx, , drop = F]
  p_rand_trt <- p_rand_trt[indx, , drop = F]

  interim_results <- list(
    arm_quantities = loo::nlist(
      n_agg_enr,
      y_agg_enr,
      n_agg_obs,
      y_agg_obs,
      p_rand,
      reg_bes,
      arm_mu,
      arm_var,
      arm_lb,
      arm_ub
    ),

    par_quantities = loo::nlist(
      par_mu,
      par_var,
      par_lb,
      par_ub
    ),

    trt_quantities = loo::nlist(
      n_agg_enr_trt,
      y_agg_enr_trt,
      n_agg_obs_trt,
      y_agg_obs_trt,
      p_rand_trt,
      p_rand_trt_raw,
      trt_mu,
      trt_lb,
      trt_ub,
      trt_eff,
      trt_fut,
      trt_bes,
      trt_bes_all,
      trt_bes_act,
      trt_in_best,
      trt_in_best_act,
      is_trt_eff,
      is_trt_active,
      is_trt_fut,
      is_trt_sup,
      is_trt_inf
    )
  )

  return(loo::nlist(
    interim_results,
    final_results
  ))
}




#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_arm_quantities <- function(dat, final = TRUE, ...) {
  quant <- ifelse(final, "final_results", "interim_results")
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[[quant]][["arm_quantities"]]
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
#' @param dat The results of `ascot_trial2` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_par_quantities <- function(dat, final = TRUE, ...) {
  quant <- ifelse(final, "final_results", "interim_results")
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[[quant]][["par_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "parameter", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "parameter")) %>%
      dplyr::mutate(parameter = forcats::fct_inorder(parameter)) %>%
      dplyr::arrange(interim, parameter)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}



#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_trt_quantities <- function(dat, final = TRUE, ...) {
  quant <- ifelse(final, "final_results", "interim_results")
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[[quant]][["trt_quantities"]]
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
#' @param dat The results of `ascot_trial2` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_ctr_quantities <- function(dat, final = TRUE, ...) {
  quant <- ifelse(final, "final_results", "interim_results")
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[[quant]][["ctr_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "contrast", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "contrast")) %>%
      dplyr::mutate(contrast = forcats::fct_inorder(contrast)) %>%
      dplyr::arrange(interim, contrast)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_trial_quantities <- function(dat, ...) {
  dplyr::bind_rows(lapply(dat, function(i) {
    tq <- i[["final_results"]][["trial_quantities"]]
    tibble::enframe(tq) %>%
      tidyr::unnest_wider(value) %>%
      tidyr::gather(treatment, value, -name) %>%
      tidyr::spread(name, value)
  }), .id = "trial")
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_result_quantities <- function(dat, ...) {
  tibble::enframe(lapply(dat, function(x) x[["final_results"]][["result_quantities"]]), name = "trial") %>%
    tidyr::unnest_wider(value)
}

