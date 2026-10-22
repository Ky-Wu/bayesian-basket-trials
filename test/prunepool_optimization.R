library(poisbinom)
library(parallel)
mc.cores <- detectCores() - 1

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  #"Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  #"Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  #"Half" = c(0.2, 0.20, 0.05, 0.05)
)
p0 <- 0.05
p1 <- 0.20
B <- 4
n_i <- c(9, 9, 9, 10)
n_b <- c(21, 21, 21, 21)
n_sim <- 20000
p0 <- 0.05
n_s <- length(scenarios)
k <- length(n_b)
constructPPRule <- function(R1, r, alpha2) {
  function(n_i, n_b, yi, y, p0) {
    stopifnot(length(yi) == length(y))
    K <- length(y)
    futility <- sum(yi) < R1
    prune <- y < r
    if (futility) {
      SS <- sum(n_i)
      decisions <- rep(FALSE, K)
    } else {
      y_notpruned <- y[!prune]
      SS <- sum(n_b)
      decisions <- rep(FALSE, K)
      if (!all(prune)) {
        R2 <- qpoisbinom(1 - alpha2, rep(p0, sum(n_b[!prune])))
        #p_value <- pbinom(sum(y_notpruned), sum(n_b[!prune]), p0, lower.tail = FALSE)
        decisions[!prune] <- (sum(y_notpruned) >= R2)
        #decisions[!prune] <- (p_value < alpha2)
      }
    }
    list(decisions = decisions,
         prune = prune,
         SS = SS)
  }
}

grid <- expand.grid(n_i = 5:15, R1 = 1:5, n_b = 15:30,
                    ESS = NA, typeI_error = NA,
                    power = NA, r = NA, alpha2 = NA)
grid <- grid[grid$n_b > grid$n_i,]

evaluateESS <- function(n_i, R1, n_b, B = 4, p0 = 0.05, p1 = 0.20) {
  interim_accept <- 1 - pbinom(R1, B * n_i, p0)
  B * (n_i + interim_accept * (n_b - n_i))
}
evaluateESS(grid$n_i[2], grid$R1[2], grid$n_b[2])
grid$ESS <- mapply(evaluateESS,
                   grid$n_i, grid$R1, grid$n_b)
grid <- grid[order(grid$ESS),]
grid$acceptable <- NA

typeI_error_requirement <- 0.05
power_requirement <- 0.80
found_optimal_design <- FALSE
i <- 1
for (i in seq_len(nrow(grid))) {
  n_i <- grid$n_i[i]
  R1 <- grid$R1[i]
  n_b <- grid$n_b[i]
  print(paste("Evaluating: n_i =", n_i, "| R1 = ", R1, "| n_b =", n_b))
  HA_interim_accept <- 1 - pbinom(R1, B * n_i, p1)
  if (HA_interim_accept <= power_requirement) {
    grid$acceptable[i] <- FALSE
    next
  }
  ### H0 sims ###
  yi_H0 <- matrix(rbinom(n_sim * B, rep(n_i, times = n_sim * B),
                         prob = rep(p0, times = n_sim * B)), nrow = n_sim, byrow = TRUE)
  interim_accept <- rowSums(yi_H0) >= R1
  interim_accept_p <- mean(interim_accept)
  n2_sim <- sum(interim_accept)
  y_H0 <- matrix(rbinom(n2_sim * B, rep(n_b - n_i, times = n2_sim * B),
                        prob = rep(p0, times = n2_sim * B)), nrow = n2_sim, ncol = B,
                 byrow = TRUE)
  y_H0 <- yi_H0[interim_accept,] + y_H0

  ### HA sims ###

  yi_HA <- matrix(rbinom(n_sim * B, rep(n_i, times = n_sim * B),
                         prob = rep(p1, times = n_sim * B)), nrow = n_sim, byrow = TRUE)
  interim_accept_HA <- rowSums(yi_HA) >= R1
  interim_accept_HA_p <- mean(interim_accept_HA)
  n2_sim_HA <- sum(interim_accept_HA)
  y_HA <- matrix(rbinom(n2_sim_HA * B, rep(n_b - n_i, times = n2_sim_HA * B),
                        prob = rep(p1, times = n2_sim_HA * B)), nrow = n2_sim_HA, ncol = B,
                 byrow = TRUE)
  y_HA <- yi_HA[interim_accept_HA,] + y_HA

  subgrid <- expand.grid(r = 0:7,
                         p_cutoff = NA, typeI_error = NA, power = 0)
  for (j in seq_len(nrow(subgrid))) {
    # initialize settings #
    r <- subgrid$r[j]

    # compute p-value cutoff using H0 sims #
    H0_prune <- y_H0 < r
    H0_p_values <- apply(y_H0, MARGIN = 1, function(y) {
      1 - pbinom(sum(y[y >= r]) - 1, n_b * sum(y >= r), p0)
    })
    possible_cutoffs <- unique(c(H0_p_values, 1.0))
    cutoff_typeI_errors <- vapply(possible_cutoffs, function(p) {
      # decision rule: reject if p-value < cutoff
      mean(H0_p_values < p) * interim_accept_p
    }, numeric(1))
    indx <- cutoff_typeI_errors <= typeI_error_requirement
    if (any(indx)) {
      cutoff <- max(possible_cutoffs[indx])
      typeI_error <- mean(H0_p_values < cutoff) * interim_accept_p

      ### evaluate power ###
      HA_prune <- y_HA < r
      HA_p_values <- apply(y_HA, MARGIN = 1, function(y) {
        1 - pbinom(sum(y[y >= r]) - 1, n_b * sum(y >= r), p0)
      })
      HA_efficacy <- !HA_prune & (HA_p_values < cutoff)
      power <- mean(HA_efficacy) * interim_accept_HA_p

      subgrid$p_cutoff[j] <- cutoff
      subgrid$typeI_error[j] <- typeI_error
      subgrid$power[j] <- power
    }
  }
  stopifnot(any(subgrid$power > 0))
  optim_indx <- subgrid$power == max(subgrid$power)
  r <- subgrid$r[optim_indx]
  alpha2 <- subgrid$p_cutoff[optim_indx]
  grid$typeI_error[i] <- subgrid$typeI_error[optim_indx]
  grid$power[i] <- subgrid$power[optim_indx]
  grid$r[i] <- r
  grid$alpha2[i] <- alpha2
  if (subgrid$typeI_error[optim_indx] <= typeI_error_requirement &
      subgrid$power[optim_indx] >= power_requirement) {
    grid$acceptable[i] <- TRUE
    print(paste("Found optimal design: n_i =", n_i, "| R1 = ", R1, "| n_b =", n_b,
                "| r =", r, "| alpha2 =", alpha2))
  }
  grid$acceptable[i] <- FALSE
}

grid[!is.na(grid$acceptable) & grid$acceptable,]

# for fixed sample size comparison

grid[grid$n_i == 9 & grid$n_b == 21,]

PPRule <- constructPPRule(3, 2, 0.018)

yi_sim_list <- lapply(scenarios, function(scenario) {
  tp <- as.vector(scenario)
  yi_sim <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                          prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
  yi_sim
})
y_sim_list <- lapply(seq_along(scenarios), function(i) {
  scenario <- scenarios[[i]]
  tp <- as.vector(scenario)
  y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b - n_i, times = n_sim),
                         prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
  y_sim <- yi_sim_list[[i]] + y_sim
  y_sim
})



evaluatePPScenario <- function(PPRule, yi_sim, y_sim, promising) {
  res <- lapply(seq_len(nrow(yi_sim)),
                function(i) PPRule(n_i, n_b, yi_sim[i,], y_sim[i,], p0 = p0))
  SS <- vapply(res, function(x) x$SS, numeric(1))
  reject <- vapply(res, function(x) x$decisions, logical(ncol(yi_sim)))
  FWER <- vapply(res, function(x) any(x$decisions[!promising]), logical(1))
  list(ESS = mean(SS),
       reject = rowMeans(reject),
       FWER = mean(FWER))
}
