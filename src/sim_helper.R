library(parallel)
library(dplyr)
mc.cores <- parallel::detectCores() - 1

checkTruePartition <- function(part, p) {
  homogeneous_groups <- all(vapply(part, function(x) {
    length(unique(p[x])) == 1
  }, logical(1)))
  if (!homogeneous_groups) {
    FALSE
  } else {
    block_elements <- vapply(part, function(x) unique(p[x]), numeric(1))
    length(unique(block_elements)) == length(block_elements)
  }
}


calibrateGlobalNullQ <- function(n_b, p0, FWER_limit, n_sim = 1000,
                                 efficacyFunction, y_sim = NULL) {
  if (is.null(y_sim)) {
    y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                           prob = rep(p0, times = n_sim * length(n_b))),
                    nrow = n_sim, byrow = TRUE)
  }

  sim_Pb <- vapply(seq_len(n_sim), function(i) efficacyFunction(n_b, y_sim[i,], p0),
                   numeric(length(n_b)))
  max_Pbs <- apply(sim_Pb, 2, max)
  sorted_Pbs <- sort(max_Pbs)
  FWER_est <- vapply(sorted_Pbs, function(x) mean(sorted_Pbs < x), numeric(1))
  q <- sorted_Pbs[which.min(FWER_est <= (1 - FWER_limit))]
  ECD <-  sum(sim_Pb < q) / n_sim
  list(threshold = q,
       estimated_FWER = mean(max_Pbs >= q),
       ECD = ECD)
}

computeGlobalPower <- function(n_b, p0, p1, q, n_sim = 1000,
                               efficacyFunction, y_sim = NULL) {
  if(is.null(y_sim)) {
    y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                           prob = rep(p1, times = n_sim * length(n_b))),
                    nrow = n_sim, byrow = TRUE)
  }
  sim_Pb <- vapply(seq_len(n_sim), function(i) efficacyFunction(n_b, y_sim[i,], p0),
                   numeric(length(n_b)))
  ECD <- sum(sim_Pb >= q) / n_sim
  basket_power <- apply(sim_Pb, 1, function(x) mean(x >= q))
  trial_power <- sum(basket_power * n_b / sum(n_b))
  list(basket_power = basket_power,
       trial_power = trial_power,
       n_b = n_b,
       ECD = ECD)
}

computeMixedPower <- function(n_b, p0, p, promising, q, n_sim = 1000,
                              efficacyFunction, y_sim = NULL) {
  stopifnot(length(n_b) == length(p))
  if(is.null(y_sim)) {
    y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                           prob = rep(p, times = n_sim)), nrow = n_sim, byrow = TRUE)
  }
  sim_Pb <- vapply(seq_len(n_sim), function(i) efficacyFunction(n_b, y_sim[i,], p0),
                   numeric(length(n_b)))
  sim_FWER <- vapply(seq_len(n_sim), function(i) any(sim_Pb[!promising,i] >= q),
                     logical(1))
  sim_FWER <- mean(sim_FWER)
  correct <- vapply(seq_len(n_sim), function(i) {
    true_positive <- sum(sim_Pb[promising, i] >= q)
    true_negative <- sum(sim_Pb[!promising, i] < q)
    sum(true_positive + true_negative)
  }, numeric(1))
  basket_positive <- apply(sim_Pb, 1, function(x) mean(x >= q))
  basket_power <- rep(NA, length(n_b))
  basket_type1_error <- rep(NA, length(n_b))
  basket_power[promising] <- basket_positive[promising]
  basket_type1_error[!promising] <- basket_positive[!promising]
  trial_power <- sum(basket_power[promising] * n_b[promising] / sum(n_b[promising]))
  list(sim_Pb = sim_Pb,
       basket_power = basket_power,
       trial_power = trial_power,
       basket_type1_error = basket_type1_error,
       FWER = sim_FWER,
       ECD = mean(correct),
       n_b = n_b,
       promising = promising)
}


simPartitionRecovery <- function(sim_p, n_b, n_sim = 1000,
                                 partition_function, y_sim = NULL) {
  stopifnot(length(sim_p) == length(n_b), n_sim >= 1)
  if(is.null(y_sim)) {
    y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                           prob = rep(sim_p, times = n_sim)), nrow = n_sim, byrow = TRUE)
  }
  correct <- vapply(seq_len(n_sim), function(i) {
    part <- partition_function(n_b, y_sim[i,])
    checkTruePartition(part, sim_p)
  }, logical(1))
  mean(correct)
}

simPartitions <- function(sim_p, n_b, n_sim = 1000, partition_function) {
  stopifnot(length(sim_p) == length(n_b), n_sim >= 1)
  y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                         prob = rep(sim_p, times = n_sim)), nrow = n_sim, byrow = TRUE)
  lapply(seq_len(n_sim), function(i) {
    partition_function(n_b, y_sim[i,])
  })
}

evaluateScenarios <- function(n_b, scenarios, p0, efficacyFunction, q = NULL,
                              FWER_limit = 0.05, n_sim = 1000, y_sim_list = NULL) {
  n_s <- length(scenarios)
  k <- length(n_b)
  if (is.null(names(scenarios))) {
    scenario_names <- seq_along(scenarios)
  } else {
    scenario_names <- names(scenarios)
  }
  if(is.null(y_sim_list)) {
    y_sim_list <- lapply(scenarios, function(scenario) {
      tp <- as.vector(scenario)
      y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
                             prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
    })
  }
  if (is.null(q)) {
    calibrate <- calibrateGlobalNullQ(n_b = n_b, p0 = p0, FWER_limit = FWER_limit,
                                      n_sim = n_sim, efficacyFunction = efficacyFunction)
    q <- calibrate$threshold
  }
  evaluateScenario <- function(tp, promising, y_sim) {
    computeMixedPower(n_b = n_b, p0 = p0, p = tp, promising, q = q,
                      n_sim = n_sim, efficacyFunction = efficacyFunction,
                      y_sim = y_sim)
  }
  res <- mclapply(seq_len(n_s), function(i) {
    tp <- as.vector(scenarios[,i])
    promising <- tp > p0
    evaluateScenario(tp, promising, y_sim_list[[i]])
  }, mc.cores = mc.cores)
  type1_errors <- t(vapply(res, function(x) x$basket_type1_error, numeric(k)))
  basket_power <- t(vapply(res, function(x) x$basket_power, numeric(k)))
  ecds <- vapply(res, function(x) x$ECD, numeric(1))
  FWERs <- vapply(res, function(x) x$FWER, numeric(1))
  list(scenarios = t(scenarios),
       type1_errors = type1_errors,
       basket_power = basket_power,
       ecds = ecds,
       mean_ecds = mean(ecds),
       ef = efficacyFunction,
       FWERs = FWERs,
       FWER_limit = FWER_limit,
       n_sim = n_sim)
}


evalInterim <- function(n_i, y, p0, q) {
  pbinom(sum(n_i), sum(y), p0, lower.tail = FALSE) > q
}

computeExactBinomSS <- function(p0, p1, alpha1, beta1, B) {
  n <- B
  stopped <- FALSE
  while(!stopped) {
    r <- qbinom(alpha1, size = n, prob = p0, lower.tail = FALSE)
    beta <- pbinom(r, size = n, prob = p1)
    if (beta <= beta1) {
      stopped <- TRUE
    } else {
      n <- n + B
    }
  }
  return(n)
}

calibrateTwoStage <- function(B, p0, p1, alpha1, beta1, alpha2, beta2, n_sim = 1000,
                              efficacyFunction, y_sim = NULL, increment_start = 1) {
  total_ni <- computeExactBinomSS(p0, p1, alpha1, beta1, B)
  # pooled interim analysis: reject if total responses is less than or equal to r
  interim_threshold <- qbinom(1 - alpha1, size = total_ni, prob = p0, lower.tail = TRUE)
  p_passinterim <- pbinom(interim_threshold, size = total_ni, prob = p0, lower.tail = FALSE)
  HA_passinterim <- pbinom(interim_threshold, size = total_ni, prob = p1, lower.tail = FALSE)
  if (HA_passinterim < (1 - beta2)) {
    stop("Interim power loss is not recoverable under alternate, increase alpha1 or decrease beta1")
  }
  twostage_ef <- function(n_b, yi, y, p0) {
    if (sum(yi) <= interim_threshold) return(rep(0, length(n_b)))
    efficacyFunction(n_b, y, p0)
  }
  n_i <- rep(total_ni / B, B)
  n_b <- rep(total_ni / B + increment_start, B)
  yi_H0 <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                         prob = rep(p0, times = n_sim)), nrow = n_sim, byrow = TRUE)
  yi_HA <-  matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                          prob = rep(p1, times = n_sim)), nrow = n_sim, byrow = TRUE)
  y_H0 <-  yi_H0 + rbinom(length(n_i) * n_sim, n_b - n_i, p0)
  y_HA <- yi_HA + rbinom(length(n_i) * n_sim, n_b - n_i, p1)
  accumulating_ss <- TRUE
  while(accumulating_ss) {
    cat("Interim Basket sample size:", n_i[1],
        ", Stage 2 Basket sample size:", n_b[1], "\n")
    sim_Pb <- vapply(seq_len(n_sim), function(i) twostage_ef(n_b, yi_H0[i,], y_H0[i,], p0),
                     numeric(length(n_b)))
    sorted_Pbs <- sort(unique(as.vector(sim_Pb)), decreasing = FALSE)
    possible_alpha2s <- vapply(sorted_Pbs, function(x) mean(sim_Pb >= x), numeric(1))
    # obtain threshold that controls type I error
    # assumes equal weighting between baskets and equal sample sizes
    pp_threshold <- sorted_Pbs[which.max(possible_alpha2s <= alpha2)]
    alpha2_estimate <- mean(sim_Pb >= pp_threshold)
    # estimate trial power
    sim_Pb <- vapply(seq_len(n_sim), function(i) twostage_ef(n_b, yi_HA[i,], y_HA[i,], p0),
                     numeric(length(n_b)))
    power <- mean(sim_Pb >= pp_threshold)
    cat("Power: ", power, "\n")
    if (power < (1 - beta2)) {
      n_b <- n_b + 1
      yi_H0 <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                             prob = rep(p0, times = n_sim)), nrow = n_sim, byrow = TRUE)
      yi_HA <-  matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                              prob = rep(p1, times = n_sim)), nrow = n_sim, byrow = TRUE)
      y_H0 <-  yi_H0 + rbinom(length(n_i) * n_sim, n_b - n_i, p0)
      y_HA <- yi_HA + rbinom(length(n_i) * n_sim, n_b - n_i, p1)
    } else {
      accumulating_ss <- FALSE
    }
  }
  ESS <- sum(n_i) + sum(n_b - n_i) * p_passinterim
  list(n_i = n_i,
       n_b = n_b,
       interim_threshold = interim_threshold,
       H0_passinterim = p_passinterim,
       alpha2_estimate = alpha2_estimate,
       pp_threshold = pp_threshold,
       ESS = ESS,
       HA_power = power,
       ef = efficacyFunction,
       alpha1 = alpha1,
       beta1 = beta1,
       alpha2 = alpha2,
       beta2 = beta2)
}

evaluateTwoStageScenario <- function(tp, promising, n_i, n_b, p0,
                                     interim_threshold, pp_threshold,
                                     efficacyFunction,
                                     yi_sim, y_sim) {
  stopifnot(identical(dim(yi_sim), dim(y_sim)),
            length(n_i) == length(n_b),
            length(promising) == length(tp))
  n_sim <- nrow(yi_sim)
  B <- ncol(yi_sim)
  interim_responses <- rowSums(yi_sim)
  interim_accept <- matrix(interim_responses > interim_threshold,
                           nrow = n_sim, ncol = B)
  sim_Pb <- vapply(seq_len(n_sim), function(i) efficacyFunction(n_b, y_sim[i,], p0),
                   numeric(B))
  final_accept <- sim_Pb >= pp_threshold
  accept <- interim_accept & t(final_accept)
  sim_FWER <- apply(accept, MARGIN = 1, function(x) any(x[!promising]))
  sim_FWER <- mean(sim_FWER)
  correct <- vapply(seq_len(n_sim), function(i) {
    true_positive <- sum(accept[i, promising])
    true_negative <- sum(!accept[i, !promising])
    sum(true_positive + true_negative)
  }, numeric(1))
  accept_means <- apply(accept, 2, mean)
  basket_power <- rep(NA, length(n_b))
  basket_type1_error <- rep(NA, length(n_b))
  basket_power[promising] <- accept_means[promising]
  basket_type1_error[!promising] <- accept_means[!promising]
  trial_power <- sum(basket_power[promising] * n_b[promising] / sum(n_b[promising]))
  ESS <- sum(n_i) + mean(interim_accept) * sum(n_b - n_i)
  list(basket_power = basket_power,
       trial_power = trial_power,
       basket_type1_error = basket_type1_error,
       FWER = sim_FWER,
       ECD = mean(correct),
       ESS = ESS,
       n_b = n_b,
       interim_acceptance = mean(interim_accept),
       promising = promising)
}

evaluateTwoStageScenarios <- function(n_i, n_b, scenarios, p0,
                                      interim_threshold,
                                      pp_threshold,
                                      efficacyFunction,
                                      n_sim = 1000,
                                      yi_sim_list = NULL,
                                      y_sim_list = NULL) {
  n_s <- length(scenarios)
  k <- length(n_b)
  if (is.null(names(scenarios))) {
    scenario_names <- seq_along(scenarios)
  } else {
    scenario_names <- names(scenarios)
  }
  if(is.null(yi_sim_list)) {
    yi_sim_list <- lapply(scenarios, function(scenario) {
      tp <- as.vector(scenario)
      yi_sim <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                             prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
    })
  }
  if(is.null(y_sim_list)) {
    y_sim_list <- lapply(seq_along(scenarios), function(i) {
      scenario <- scenarios[[i]]
      tp <- as.vector(scenario)
      y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b - n_i, times = n_sim),
                             prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
      y_sim <- yi_sim_list[[i]] + y_sim
    })
  }
  res <- mclapply(seq_len(n_s), function(i) {
    tp <- as.vector(scenarios[,i])
    promising <- tp > p0
    evaluateTwoStageScenario(tp, promising, n_i, n_b, p0,
                             interim_threshold, pp_threshold,
                             efficacyFunction,
                             yi_sim_list[[i]], y_sim_list[[i]])
  }, mc.cores = mc.cores)
  type1_errors <- t(vapply(res, function(x) x$basket_type1_error, numeric(k)))
  basket_power <- t(vapply(res, function(x) x$basket_power, numeric(k)))
  ecds <- vapply(res, function(x) x$ECD, numeric(1))
  FWERs <- vapply(res, function(x) x$FWER, numeric(1))
  ESSs <- vapply(res, function(x) x$ESS, numeric(1))
  list(scenarios = t(scenarios),
       type1_errors = type1_errors,
       basket_power = basket_power,
       ecds = ecds,
       mean_ecds = mean(ecds),
       ef = efficacyFunction,
       FWERs = FWERs,
       ESS = ESSs,
       n_b = n_b,
       n_i = n_i,
       interim_threshold = interim_threshold,
       pp_threshold = pp_threshold,
       n_sim = n_sim)
}

simResponses <- function(n_i, n, pb, n_sim, sort = TRUE) {
  yi <- rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
               prob = rep(pb, times = n_sim)) %>%
    matrix(nrow = n_sim, byrow = TRUE)
  y <- rbinom(length(n) * n_sim, rep(n - n_i, times = n_sim),
              prob = rep(pb, times = n_sim)) %>%
    matrix(nrow = n_sim, byrow = TRUE)
  y <- yi + y
  if (sort) {
    y <- t(apply(yi + y, 1, sort))
    yi <- t(apply(yi, 1, sort))
  }
  out <- cbind(yi, y) %>%
    as_tibble() %>%
    group_by_all() %>%
    summarize(count = n()) %>%
    ungroup()
  colnames(out) <- c(paste0("ib", seq_along(n_i)), paste0("b", seq_along(n_b)), "n")
  out
}

calibrateTwoStageEF <- function(p0 = 0.05, p1 = 0.20, n_i, n_b,
                            efficacyFunction, alpha2, n_sim = 10000) {
  k <- length(n_i)
  H0_resp <- simResponses(n_i, n_b, p0, n_sim)
  yi_H0 <- as.matrix(H0_resp[,1:k])
  y_H0 <- as.matrix(H0_resp[,(k+1):(2*k)])
  H0_w <- H0_resp$n / sum(H0_resp$n)
  HA_resp <- simResponses(n_i, n_b, p1, n_sim)
  yi_HA <- as.matrix(HA_resp[,1:k])
  y_HA <- as.matrix(HA_resp[,(k+1):(2*k)])
  HA_w <- HA_resp$n / sum(HA_resp$n)
  Pbi_H0 <- vapply(seq_len(nrow(yi_H0)), function(i) efficacyFunction(n_i, yi_H0[i,], p0), numeric(k))
  Pbi_HA <- vapply(seq_len(nrow(yi_HA)), function(i) efficacyFunction(n_i, yi_HA[i,], p0), numeric(k))
  sorted_Pbis <- sort(unique(as.vector(Pbi_H0)), decreasing = FALSE)
  possible_thresholds <- data.table(it = sorted_Pbis)
  res <- lapply(seq_len(nrow(possible_thresholds)), function(i) {
    it = possible_thresholds[i,]$it
    cat("Testing interim threshold:", it, "\n")
    interim_accept <- apply(Pbi_H0, 2, function(x) x > it)
    Pb_H0 <- lapply(seq_len(nrow(y_H0)), function(j) {
      indx <- interim_accept[,j]
      if (any(indx)) {
        efficacyFunction(n_b[indx], y_H0[j,indx], p0)
      } else {
        numeric(0)
      }
    })
    w <- rep(H0_w, times = vapply(Pb_H0, length, numeric(1)))
    #indx <- order(unlist(Pb_H0), decreasing = FALSE)
    Pb_H0 <- unlist(Pb_H0)
    Pb_H0_unique <- sort(unique(Pb_H0), decreasing = TRUE)
    w_unique <- vapply(Pb_H0_unique, function(p) {
      indx <- Pb_H0 == p
      sum(w[indx])
    }, numeric(1))
    w_unique <- w_unique / sum(w_unique)
    cw <- cumsum(w_unique)
    # declare efficacy if > final threshold
    # examine cumulative prob mass for ft and below
    ft <- Pb_H0_unique[max(which(cw <= alpha2))]
    #sum(w_unique[Pb_H0_unique > ft])
    # find power
    interim_accept <- apply(Pbi_HA, 2, function(x) x > it)
    HA_decisions <- lapply(seq_len(nrow(yi_HA)), function(j) {
      i_accept <- interim_accept[,j]
      if (all(i_accept == FALSE)) {
        out <- i_accept
      } else {
        Pb <- efficacyFunction(n_b[i_accept], y_HA[j, i_accept], p0)
        out <- rep(FALSE, k)
        out[i_accept] <- Pb > ft
      }
      out
    })
    power <- vapply(seq_along(HA_decisions), function(i) {
      mean(HA_decisions[[i]] * HA_w[i])
    }, numeric(1))
    list(ft = ft, power = sum(power))
  }, numeric(1))
  possible_thresholds$power <- vapply(res, function(x) x$power, numeric(1))
  possible_thresholds$ft <- vapply(res, function(x) x$ft, numeric(1))
  possible_thresholds
}

