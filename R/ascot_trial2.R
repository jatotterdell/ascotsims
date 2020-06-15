#' Run ASCOT trial simulation
#'
#'
#' @export
ascot_trial2 <- function(
  b,
  n_seq,
  delta = log(1.2),
  effective_thres = 0.99,
  equivalent_thres = 0.9,
  futility_thres = 0.01,
  ineffective_thres = 0.1,
  best_thres = 0.95,
  av_less_eff_thres = 0.01,
  rar_control = FALSE,
  rar_best = FALSE,
  perm_drop = TRUE,
  use_mwud = TRUE,
  use_optimal = FALSE,
  early_t = 7,
  h = 0.5,
  trans_prob = 0.75,
  mc_draw = 2e4,
  return_all = FALSE,
  ...
) {

  if(length(b) != 10) stop("wrong number of parameters in b. Should be 10")

  # Priors
  M0 <- rep(0, 10)
  S0 <- diag(100, 10)

  # Design matrix for domains
  Xa <- rbind(0, cbind(diag(1, 4), 0), c(1, 1, 0, 0, 1))
  colnames(Xa) <- c("a1", "a2", "a3", "a4", "a1:a2")
  rownames(Xa) <- c("SoC", "a1", "a2", "a3", "a4", "a1+a2")
  Xb <- rbind(0, diag(1, 3))
  colnames(Xb) <- c("b1", "b2", "b3")
  rownames(Xb) <- c("SoC", "b1", "b2", "b3")
  Xc <- rbind(0, diag(1))
  colnames(Xc) <- c("c1")
  rownames(Xc) <- c("SoC", "c1")
  X <- expand_grid(as.data.frame(Xa), as.data.frame(Xb), as.data.frame(Xc))
  X <- cbind("SoC" = 1, as.matrix(X))
  rownames(X) <- sapply(1:nrow(X), function(i)
    paste(colnames(X)[!grepl("(SoC|:)", colnames(X))][X[i, !grepl("(SoC|:)", colnames(X))] == 1], collapse = "+"))
  rownames(X)[1] <- "SoC"

  # Arm indicator designs
  XIa <- diag(1, 6)
  colnames(XIa) <- rownames(XIa) <- c("a0", "a1", "a2", "a3", "a4", "a1+a2")
  XIb <- diag(1, 4)
  colnames(XIb) <- rownames(XIb) <- c("b0", "b1", "b2", "b3")
  XIc <- diag(1, 2)
  colnames(XIc) <- rownames(XIc) <- c("c0", "c1")
  XI <- expand_grid(as.data.frame(XIa), as.data.frame(XIb), as.data.frame(XIc))
  XI <- as.matrix(XI)
  rownames(XI) <- sapply(1:nrow(XI), function(i)
    paste(colnames(XI)[!grepl("(0|:)", colnames(XI))][XI[i, !grepl("(0|:)", colnames(XI))] == 1], collapse = "+"))
  rownames(XI)[1] <- "SoC"

  # Arm to domain map
  XM <- expand_grid(a = 1:6, b = 1:4, c = 1:2)

  # Initialise allocation ratios
  if(use_optimal) {
    ratio_dom <- sapply(list(A = Xa, B = Xb, C = Xc),
                    function(j) {
                      J <- dim(j)[2]
                      r <- c(sqrt(J), rep(1, J))
                      r / sum(r)
                    })
  } else {
    ratio_dom <- sapply(list(A = Xa, B = Xb, C = Xc),
                    function(j) rep(1 / (dim(j)[1]), dim(j)[1]))
  }
  ratio <- ratio_dom$A[XM$a]*ratio_dom$B[XM$b]*ratio_dom$C[XM$c]

  # Data
  p <- plogis(X %*% b)[, 1]
  n_arms <- nrow(X)
  n_pars <- ncol(X)
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
  arm_names <- rownames(X)
  par_names <- colnames(X)
  arm_dm <- list("interim" = 1:n_int, "arm" = arm_names)
  arm_dm2 <- list("interim" = 1:n_int, "arm" = arm_names[-1])
  par_dm <- list("interim" = 1:n_int, "treatment" = par_names)
  par_dm2 <- list("interim" = 1:n_int, "treatment" = par_names[-1])

  # Data aggregation storage
  n_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_enr <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  n_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  y_agg_obs <- matrix(0, n_int, n_arms, dimnames = arm_dm)
  n_agg_enr_trt <- matrix(0, n_int, ncol(XI), dimnames = list(interim = 1:n_int, treatment = colnames(XI)))
  y_agg_enr_trt <- matrix(0, n_int, ncol(XI), dimnames = list(interim = 1:n_int, treatment = colnames(XI)))
  n_agg_obs_trt <- matrix(0, n_int, ncol(XI), dimnames = list(interim = 1:n_int, treatment = colnames(XI)))
  y_agg_obs_trt <- matrix(0, n_int, ncol(XI), dimnames = list(interim = 1:n_int, treatment = colnames(XI)))

  # Allocation probabilities by domain and overall
  p_rand_a <- matrix(ratio_dom$A, n_int + 1, nrow(Xa), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xa)))
  p_rand_b <- matrix(ratio_dom$B, n_int + 1, nrow(Xb), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xb)))
  p_rand_c <- matrix(ratio_dom$C, n_int + 1, nrow(Xc), byrow = T,
                     dimnames = list("interim" = 0:n_int, "arm" = rownames(Xc)))
  p_rand <- matrix(ratio, n_int + 1, n_arms, dimnames = list("interim" = 0:n_int, "arm" = arm_names))

  # Model parameter storage
  par_mu  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_var <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_lb  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_ub  <- matrix(0, n_int, n_pars, dimnames = par_dm)
  par_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2) # Pr(b<0)
  par_fut <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2) # Pr(b<d)
  par_equ <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2) # Pr(-d<b<d)

  is_par_active      <- matrix(TRUE, n_int, n_pars, dimnames = par_dm)
  is_par_effective   <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_par_equivalent  <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_par_ineffective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_par_futile      <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)

  # Treatment storage
  trt_mu  <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_var <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_lb  <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_ub  <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_eff <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_fut <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_equ <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)
  trt_bes <- matrix(0, n_int, n_pars - 1, dimnames = par_dm2)

  is_trt_active      <- matrix(TRUE, n_int, n_pars, dimnames = par_dm)
  is_trt_effective   <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2)
  is_trt_equivalent  <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2) # Pr(sum(b)<0)
  is_trt_ineffective <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2) # Pr(sum(b)<d)
  is_trt_futile      <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2) # Pr(-d<sum(b)<d)
  is_trt_best        <- matrix(FALSE, n_int, n_pars - 1, dimnames = par_dm2) # Pr(sum(b)=min(sum(b)))

  # Regimen storage
  reg_mu  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_var <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_lb  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_ub  <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_eff <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_fut <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)
  reg_equ <- matrix(0, n_int, n_arms - 1, dimnames = arm_dm2)

  is_reg_active      <- matrix(TRUE, n_int, n_arms, dimnames = arm_dm)
  is_reg_effective   <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2)
  is_reg_equivalent  <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2) # Pr(sum(b)<0)
  is_reg_ineffective <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2) # Pr(sum(b)<d)
  is_reg_futile      <- matrix(FALSE, n_int, n_arms - 1, dimnames = arm_dm2) # Pr(-d<sum(b)<d)

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
      t_ref       <- tdat[, 2][n_seq[i]]
      n_enrol     <- sum(tdat[(n_enrolled + 1):n_max, 1] <= t_ref)
      id_enrolled <- (n_enrolled + 1):(n_enrolled + n_enrol)
      n_enrolled  <- n_enrolled + n_enrol

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
    # y_agg_obs_trt[i, 1:6] <- aggregate(y, by = XM[as.integer(x[id_analysis]), 1], FUN = sum, drop = F)[, 2]
    # y_agg_obs_trt[i, 7:10] <- aggregate(y, by = XM[as.integer(x[id_analysis]), 2], FUN = sum, drop = F)[, 2]
    # y_agg_obs_trt[i, 11:12] <- aggregate(y, by = XM[as.integer(x[id_analysis]), 3], FUN = sum, drop = F)[, 2]
    # y_agg_enr_trt[i, 1:6] <- aggregate(y, by = XM[as.integer(x[1:n_enrolled]), 1], FUN = sum, drop = F)[, 2]
    # y_agg_enr_trt[i, 7:10] <- aggregate(y, by = XM[as.integer(x[1:n_enrolled]), 2], FUN = sum, drop = F)[, 2]
    # y_agg_enr_trt[i, 11:12] <- aggregate(y, by = XM[as.integer(x[1:n_enrolled]), 3], FUN = sum, drop = F)[, 2]

    # Fit the model
    mod <- vb_mod(X, y = y_agg_obs[i, ], n = n_agg_obs[i, ], M0 = M0, S0 = S0)
    par_draws <- mvnfast::rmvn(mc_draw, mod$mu, mod$Sigma)
    colnames(par_draws) <- colnames(X)

    # Linear predictor
    eta_draws <- par_draws %*% t(X)

    # Probability of response
    p_draws <- plogis(eta_draws)

    # Effects
    reg_draws <- par_draws[, -1] %*% t(X[, -1][-1, ])
    trt_draws_a <- par_draws[, grepl("a", colnames(par_draws))] %*% t(Xa[-1, ])
    trt_draws_b <- par_draws[, grepl("b", colnames(par_draws))] %*% t(Xb[-1, ])
    trt_draws_c <- par_draws[, grepl("c", colnames(par_draws)), drop = F] %*% t(Xc[-1, , drop = F])
    trt_draws   <- cbind(trt_draws_a, trt_draws_b, trt_draws_c)

    # Parameter summaries
    hdival <- HDInterval::hdi(par_draws)
    par_mu[i, ]  <- matrixStats::colMeans2(par_draws)
    par_var[i, ] <- diag(mod$Sigma)
    par_lb[i, ]  <- hdival[1, ]
    par_ub[i, ]  <- hdival[2, ]
    par_eff[i, ] <- matrixStats::colMeans2(par_draws[, -1] < 0)
    par_fut[i, ] <- matrixStats::colMeans2(par_draws[, -1] < -delta)
    par_equ[i, ] <- matrixStats::colMeans2(abs(par_draws[, -1]) < delta)

    is_par_effective[i, ] <- par_eff[i, ] > effective_thres
    is_par_equivalent[i, ] <- par_equ[i, ] > equivalent_thres
    is_par_ineffective[i, ] <- par_eff[i, ] < ineffective_thres
    is_par_futile[i, ] <- par_fut[i, ] < futility_thres
    is_par_active[i, 1] <- !any(is_par_effective[i, ])
    is_par_active[i, -1] <- !(is_par_ineffective)[i, ]

    # Treatment summaries
    hdival <- HDInterval::hdi(trt_draws)
    trt_mu[i, ]  <- matrixStats::colMeans2(trt_draws)
    trt_lb[i, ]  <- hdival[1, ]
    trt_ub[i, ]  <- hdival[2, ]
    trt_eff[i, ] <- matrixStats::colMeans2(trt_draws < 0)
    trt_fut[i, ] <- matrixStats::colMeans2(trt_draws < -delta)
    trt_equ[i, ] <- matrixStats::colMeans2(abs(trt_draws) < delta)
    trt_bes[i, ] <- c(
      prob_best(trt_draws[, grepl("a", colnames(trt_draws))], minimum = T),
      prob_best(trt_draws[, grepl("b", colnames(trt_draws))], minimum = T),
      prob_best(trt_draws[, grepl("c", colnames(trt_draws)), drop = F], minimum = T))

    is_trt_effective[i, ] <- trt_eff[i, ] > effective_thres
    is_trt_equivalent[i, ] <- trt_equ[i, ] > equivalent_thres
    is_trt_ineffective[i, ] <- trt_eff[i, ] < ineffective_thres
    is_trt_futile[i, ] <- trt_fut[i, ] < futility_thres
    is_trt_active[i, 1] <- !any(is_trt_effective[i, ])
    is_trt_active[i, -1] <- !(is_trt_ineffective)[i, ]

    # Regimen summaries
    hdival <- HDInterval::hdi(reg_draws)
    reg_mu[i, ] <- matrixStats::colMeans2(reg_draws)
    reg_lb[i, ] <- hdival[1, ]
    reg_ub[i, ] <- hdival[2, ]
    reg_eff[i, ] <- matrixStats::colMeans2(reg_draws < 0)
    reg_fut[i, ] <- matrixStats::colMeans2(reg_draws < -delta)
    reg_equ[i, ] <- matrixStats::colMeans2(abs(reg_draws) < delta)

    is_reg_effective[i, ] <- reg_eff[i, ] > effective_thres
    is_reg_equivalent[i, ] <- reg_equ[i, ] > equivalent_thres
    is_reg_ineffective[i, ] <- reg_eff[i, ] < ineffective_thres
    is_reg_futile[i, ] <- reg_fut[i, ] < futility_thres
    is_reg_active[i, 1] <- !any(is_reg_effective[i, ])
    is_reg_active[i, -1] <- !(is_reg_ineffective)[i, ]

    # Arm summaries
    hdival <- HDInterval::hdi(p_draws)
    arm_mu[i, ] <- matrixStats::colMeans2(p_draws)
    arm_var[i, ] <- matrixStats::colVars(p_draws)
    arm_lb[i, ] <- hdival[1, ]
    arm_ub[i, ] <- hdival[2, ]


    # If drop permanently, was the arm already inactive?
    # if(perm_drop & i > 1) {
    #   is_trt_active[i, ] <- is_trt_active[i, ] & is_trt_active[i - 1, ]
    #   is_par_active[i, ] <- is_par_active[i, ] & is_par_active[i - 1, ]
    # }


    # Include control group in RAR or not?
    # if(rar_control) {
    #   p_rand[i + 1, ] <- brar_all(p_best[i, ], is_arm_active[i, ], h)
    # } else {
    #   if(rar_best) {
    #     p_rand[i + 1, ] <- fix_ctrl_brar(p_best[i, ], is_arm_active[i, ], h)
    #   } else {
    #     if(use_optimal) {
    #       ratio <- c(sqrt(sum(is_arm_active[i, -1])), rep(1, sum(is_arm_active[i, -1])))
    #       ratio <- ratio / sum(ratio)
    #     } else {
    #       ratio <- 1 / sum(is_arm_active[i, ])
    #     }
    #     p_rand[i + 1, ] <- fix_ctrl_brar2(c(0, p_superior[i, ]), is_arm_active[i, ], h, ratio[1])
    #   }
    # }

    # If everything or everything but SoC has been dropped, just stop
    # if(!any(is_trt_active[i, -1])) stopped <- TRUE
    # if(final) break
    # if(!any(is_trt_active[i, -1])) stopped <- TRUE
    # If only one active and it's superior then stop
    # if(!any(is_trt_active[i, -1] & !is_trt_effective[i, ])) stopped <- TRUE

    # if(final) break
    # if(stopped) break

  }
  p_rand <- p_rand[-1, ]

  drop_trt_at <- apply(is_trt_active, 2, function(x) which(x == FALSE)[1])

  trigger_effective_at <- apply(is_trt_effective, 2, function(x) which(x)[1])
  trigger_equivalent_at <- apply(is_trt_equivalent, 2, function(x) which(x)[1])
  trigger_ineffective_at <- apply(is_trt_ineffective, 2, function(x) which(x)[1])
  trigger_futile_at <- apply(is_trt_futile, 2, function(x) which(x)[1])

  trigger_ineffective <- ifelse(is.na(trigger_ineffective_at), FALSE, ifelse(trigger_ineffective_at < n_int, TRUE, FALSE))
  trigger_equivalent <- ifelse(is.na(trigger_equivalent_at), FALSE, ifelse(trigger_equivalent_at < n_int, TRUE, FALSE))
  trigger_effective <- ifelse(is.na(trigger_effective_at), FALSE, ifelse(trigger_effective_at < n_int, TRUE, FALSE))
  trigger_futility <- ifelse(is.na(trigger_futile_at), FALSE, ifelse(trigger_futile_at < n_int, TRUE, FALSE))

  final_effective <- is_trt_effective[i, ]
  final_ineffective <- is_trt_ineffective[i, ]
  final_equivalent <- is_trt_equivalent[i, ]
  final_futile <- is_trt_futile[i, ]

  stop_early <- stopped
  stop_at <- ifelse(stopped, i - 1, NA)
  drop_soc_at <- findfirst(!is_trt_active[, 1])

  futility <- all(trigger_futility)
  success <-any(trigger_effective)
  effective <- any(final_effective)
  ineffective <- all(final_ineffective)

  final_results <- list(
    arm_quantities = list(
      n_agg_obs = n_agg_obs[i, , drop = F],
      y_agg_obs = y_agg_obs[i, , drop = F],
      p_rand = p_rand[i, , drop = F],
      arm_mu = arm_mu[i, , drop = F],
      arm_var = arm_var[i, , drop = F],
      arm_lb = arm_lb[i, , drop = F],
      arm_ub = arm_ub[i, , drop = F]
    ),

    trt_quantities = list(
      par_mu = par_mu[i, , drop = F],
      par_var = par_var[i, , drop = F],
      par_lb = par_lb[i, , drop = F],
      par_ub = par_ub[i, , drop = F],
      par_eff = par_eff[i, , drop = F],
      par_equ = par_equ[i, , drop = F],
      par_fut = par_fut[i, , drop = F],
      trt_mu = trt_mu[i, , drop = F],
      trt_lb = trt_lb[i, , drop = F],
      trt_ub = trt_ub[i, , drop = F],
      trt_eff = trt_eff[i, , drop = F],
      trt_equ = trt_equ[i, , drop = F],
      trt_fut = trt_fut[i, , drop = F],
      is_trt_effective =is_trt_effective[i, , drop = F],
      is_trt_equivalent = is_trt_equivalent[i, , drop = F],
      is_trt_ineffective = is_trt_ineffective[i, , drop = F],
      is_trt_futile = is_trt_futile[i, , drop = F],
      is_trt_active = is_trt_active[i, , drop = F],
      is_par_active = is_par_active[i, , drop = F],
      is_par_effective = is_par_effective[i, , drop = F],
      is_par_equivalent = is_par_equivalent[i, , drop = F],
      is_par_ineffective = is_par_ineffective[i, , drop = F],
      is_par_futile = is_par_futile[i, , drop = F]
    ),

    trial_quantities = loo::nlist(
      drop_trt_at,
      trigger_effective_at, trigger_equivalent_at, trigger_ineffective_at,
      trigger_effective, trigger_equivalent, trigger_ineffective, trigger_futility,
      final_effective, final_ineffective, final_equivalent, final_futile
    ),

    result_quantities = c(
      stop_early = stop_early,
      stop_at = stop_at,
      drop_soc_at = unname(drop_soc_at),
      futility = futility,
      success = success,
      effective = effective,
      ineffective = ineffective)
  )


  if(return_all) indx <- 1:(i-1) else indx <- i-1

  n_agg_enr <- n_agg_enr[indx, , drop = F]
  y_agg_enr <- y_agg_enr[indx, , drop = F]
  n_agg_obs <- n_agg_obs[indx, , drop = F]
  y_agg_obs <- y_agg_obs[indx, , drop = F]
  p_rand <- p_rand[indx, , drop = F]
  arm_mu <- arm_mu[indx, , drop = F]
  arm_var <- arm_var[indx, , drop = F]
  arm_lb <- arm_lb[indx, , drop = F]
  arm_ub <- arm_ub[indx, , drop = F]

  par_mu <- par_mu[indx, , drop = F]
  par_var <- par_var[indx, , drop = F]
  par_lb <- par_lb[indx, , drop = F]
  par_ub <- par_ub[indx, , drop = F]
  par_eff <- par_eff[indx, , drop = F]
  par_equ <- par_equ[indx, , drop = F]
  par_fut <- par_fut[indx, , drop = F]
  trt_mu <- trt_mu[indx, , drop = F]
  trt_lb <- trt_lb[indx, , drop = F]
  trt_ub <- trt_ub[indx, , drop = F]
  trt_eff <- trt_eff[indx, , drop = F]
  trt_equ <- trt_equ[indx, , drop = F]
  trt_fut <- trt_fut[indx, , drop = F]
  is_trt_effective <-is_trt_effective[indx, , drop = F]
  is_trt_equivalent <- is_trt_equivalent[indx, , drop = F]
  is_trt_ineffective <- is_trt_ineffective[indx, , drop = F]
  is_trt_futile <- is_trt_futile[indx, , drop = F]
  is_trt_active <- is_trt_active[indx, , drop = F]
  is_par_active <- is_par_active[indx, , drop = F]
  is_par_effective <- is_par_effective[indx, , drop = F]
  is_par_equivalent <- is_par_equivalent[indx, , drop = F]
  is_par_ineffective <- is_par_ineffective[indx, , drop = F]
  is_par_futile <- is_par_futile[indx, , drop = F]

  interim_results <- list(
    arm_quantities = loo::nlist(
      n_agg_enr, y_agg_enr,
      n_agg_obs, y_agg_obs,
      p_rand, arm_mu, arm_var, arm_lb, arm_ub
    ),

    trt_quantities = loo::nlist(
      par_mu, par_var, par_lb, par_ub, par_eff, par_equ, par_fut,
      trt_mu, trt_lb, trt_ub, trt_eff, trt_equ, trt_fut,
      is_trt_effective, is_trt_equivalent, is_trt_ineffective, is_trt_active, is_trt_futile,
      is_par_active, is_par_effective, is_par_equivalent, is_par_ineffective, is_par_futile
    )
  )

  return(loo::nlist(
    interim_results,
    final_results
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
tibble_arm_quantities <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["arm_quantities"]]
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
tibble_trt_quantities <- function(dat, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["trt_quantities"]]
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

