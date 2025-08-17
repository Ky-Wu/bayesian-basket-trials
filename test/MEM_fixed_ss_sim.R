library(ggplot2)
library(data.table)
rm(list = ls())

# source(file.path(getwd(), "src", "LMEM.R"))
# source(file.path(getwd(), "src", "uniform.R"))
source(file.path(getwd(), "src", "sim_helper.R"))
# source(file.path(getwd(), "src", "pooled.R"))
#source(file.path(getwd(), "src", "twostage_prunepool.R"))
source(file.path(getwd(), "src", "MEM.R"))

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.20, 0.05, 0.05)
)

n_i <- c(9, 9, 9, 9)
n_b <- c(21, 21, 21, 21)
p0 <- 0.05
p1 <- 0.20
alpha2 <- 0.07
beta2 <- 0.07

computeESSFixedSS <- function(p0 = 0.05, p1 = 0.20, n_i, n_b, r,
                              efficacyFunction, alpha2, n_sim = 10000) {
  alpha1 <- 1 - pbinom(r, size = sum(n_i), prob = p0, lower.tail = TRUE)
  ESS <- sum(n_i) + alpha1 * (sum(n_b) - sum(n_i))
  twostage_ef <- function(n_b, yi, y, p0) {
    if (sum(yi) <= r) return(rep(0, length(n_b)))
    efficacyFunction(n_b, y, p0)
  }
  yi_H0 <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                         prob = rep(p0, times = n_sim)), nrow = n_sim, byrow = TRUE)
  yi_HA <-  matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                          prob = rep(p1, times = n_sim)), nrow = n_sim, byrow = TRUE)
  y_H0 <-  yi_H0 + rbinom(length(n_i) * n_sim, n_b - n_i, p0)
  y_HA <- yi_HA + rbinom(length(n_i) * n_sim, n_b - n_i, p1)
  interim_accept <- matrix(rowSums(yi_H0) > r,
                           nrow = n_sim, ncol = length(n_i))
  Pb_H0 <- vapply(seq_len(n_sim), function(i) twostage_ef(n_b, yi_H0[i,], y_H0[i,], p0),
                  numeric(length(n_b)))
  sorted_Pbs <- sort(unique(as.vector(Pb_H0)), decreasing = FALSE)
  possible_alpha2s <- vapply(sorted_Pbs, function(x) mean(t(interim_accept) & Pb_H0 >= x), numeric(1))
  # obtain threshold that controls type I error
  # assumes equal weighting between baskets and equal sample sizes
  if (any(possible_alpha2s <= alpha2)) {
    pp_threshold <- sorted_Pbs[which.max(possible_alpha2s <= alpha2)]
  } else {
    pp_threshold <- 1
  }
  alpha2_estimate <- mean(t(interim_accept) & Pb_H0 >= pp_threshold)
  # estimate trial power
  Pb_HA <- vapply(seq_len(n_sim), function(i) twostage_ef(n_b, yi_HA[i,], y_HA[i,], p0),
                  numeric(length(n_b)))
  HA_interim_accept <- matrix(rowSums(yi_HA) > r,
                              nrow = n_sim, ncol = length(n_i))
  power <- mean(t(HA_interim_accept) & Pb_HA >= pp_threshold)
  list(n_i = n_i,
       n_b = n_b,
       interim_threshold = r,
       H0_passinterim = alpha1,
       alpha2_estimate = alpha2_estimate,
       pp_threshold = pp_threshold,
       ESS = ESS,
       HA_power = power,
       ef = efficacyFunction,
       alpha1 = alpha1,
       alpha2 = alpha2)
}

param_grid <- data.table(ESS = 0,
                         power = 0,
                         type1_error = 0,
                         pp_threshold = 0,
                         r = 2)

### MEM ###

for (i in 1:nrow(param_grid)) {
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           MEMEfficacy, alpha2, n_sim = 400)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
#LMEM_setting <- param_grid[ESS == min(ESS[power >= 1 - beta2]),]
MEM_setting <- param_grid
MEM_res <- evaluateTwoStageScenarios(n_i,
                                     n_b,
                                     scenarios, p0,
                                     MEM_setting$r,
                                     MEM_setting$pp_threshold,
                                     MEMEfficacy, n_sim = 400)

dir.create(file.path(getwd(), "output", "fixedss_sim"), showWarnings = FALSE)
saveRDS(MEM_res, file.path(getwd(), "output", "fixedss_sim", "MEM_res.rds"))
