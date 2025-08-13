library(ggplot2)
library(data.table)
rm(list = ls())
source(file.path(getwd(), "src", "maxent_MEM.R"))
source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "sim_helper.R"))
scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.20, 0.05, 0.05)
)

LMEM_ef <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1/2, b = 1/2,
                                                   d1 = 0, d2 = 3,
                                                   WW_method = FALSE)

n_i <- c(9, 9, 9, 9)
n_b <- c(21, 21, 21, 21)
p0 <- 0.05
p1 <- 0.20
alpha2 <- 0.07
n_sim <- 5000
k <- length(n_i)
H0_resp <- simResponses(n_i, n_b, p0, n_sim)
yi_H0 <- as.matrix(H0_resp[,1:k])
y_H0 <- as.matrix(H0_resp[,(k+1):(2*k)])
H0_w <- H0_resp$n / sum(H0_resp$n)
HA_resp <- simResponses(n_i, n_b, p1, n_sim)
yi_HA <- as.matrix(HA_resp[,1:k])
y_HA <- as.matrix(HA_resp[,(k+1):(2*k)])
HA_w <- HA_resp$n / sum(HA_resp$n)
Pbi_H0 <- vapply(seq_len(nrow(yi_H0)), function(i) LMEM_ef(n_i, yi_H0[i,], p0), numeric(k))
Pbi_HA <- vapply(seq_len(nrow(yi_HA)), function(i) LMEM_ef(n_i, yi_HA[i,], p0), numeric(k))
sorted_Pbis <- sort(unique(as.vector(Pbi_H0)), decreasing = FALSE)
possible_thresholds <- data.table(it = sorted_Pbis)
res <- lapply(seq_len(nrow(possible_thresholds)), function(i) {
  it = possible_thresholds[i,]$it
  cat("Testing interim threshold:", it, "\n")
  interim_accept <- apply(Pbi_H0, 2, function(x) x > it)
  ESS <- vapply(seq_len(ncol(interim_accept)), function(i) {
    acc <- interim_accept[,i]
    ess <- sum(n_b[acc]) + sum(n_i[!acc])
    ess * H0_w[i]
  }, numeric(1)) %>%
    sum()
  Pb_H0 <- vapply(seq_len(nrow(y_H0)), function(j) {
    indx <- interim_accept[,j]
    out <- numeric(k)
    if (any(indx)) {
      out[indx] <- LMEM_ef(n_b[indx], y_H0[j,indx], p0)
    }
    out
  }, numeric(k))
  w <- matrix(H0_w, nrow = nrow(Pb_H0), ncol = ncol(Pb_H0), byrow = TRUE)
  #indx <- order(unlist(Pb_H0), decreasing = FALSE)
  Pb_H0_unique <- sort(unique(as.vector(Pb_H0)), decreasing = TRUE)
  alpha_estimates <- vapply(Pb_H0_unique, function(p) {
    H0_decisions <- vapply(seq_len(ncol(Pb_H0)), function(i) {
      mean((Pb_H0[,i] > p) * H0_w[i])
    }, numeric(1))
    sum(H0_decisions)
  }, numeric(1))
  # declare efficacy if final probability > final threshold
  # examine cumulative prob mass for ft and below
  ft <- Pb_H0_unique[max(which(alpha_estimates <= alpha2))]
  alpha2_estimate <- alpha_estimates[max(which(alpha_estimates <= alpha2))]
  # find power
  interim_accept <- apply(Pbi_HA, 2, function(x) x > it)
  HA_decisions <- lapply(seq_len(nrow(yi_HA)), function(j) {
    i_accept <- interim_accept[,j]
    if (all(i_accept == FALSE)) {
      out <- i_accept
    } else {
      Pb <- LMEM_ef(n_b[i_accept], y_HA[j, i_accept], p0)
      out <- rep(FALSE, k)
      out[i_accept] <- Pb > ft
    }
    out
  })
  power <- vapply(seq_along(HA_decisions), function(i) {
    mean(HA_decisions[[i]] * HA_w[i])
  }, numeric(1))
  list(ft = ft,
       power = sum(power),
       ESS = ESS,
       alpha2_estimate = alpha2_estimate)
})
possible_thresholds$power <- vapply(res, function(x) x$power, numeric(1))
possible_thresholds$ft <- vapply(res, function(x) x$ft, numeric(1))
possible_thresholds$ESS <- vapply(res, function(x) x$ESS, numeric(1))
possible_thresholds$alpha2_estimate <- vapply(res, function(x) x$alpha2_estimate, numeric(1))
indx <- which.max(possible_thresholds$ESS <= 50.2)
possible_thresholds[indx,]
it <- possible_thresholds[indx,]$it
ft <- possible_thresholds[indx,]$ft


scenarios_res <- lapply(seq_len(ncol(scenarios)), function(i) {
  pb <- scenarios[[i]]
  promising <- pb > p0
  resp <- simResponses(n_i, n_b, pb, n_sim, sort = FALSE)
  yi <- as.matrix(resp[,1:k])
  y <- as.matrix(resp[,(k+1):(2*k)])
  w <- resp$n / sum(resp$n)
  Pbi <- vapply(seq_len(nrow(yi)), function(i) LMEM_ef(n_i, yi[i,], p0), numeric(k))
  interim_accept <- apply(Pbi, 2, function(x) x > it)
  ESS <- vapply(seq_len(ncol(interim_accept)), function(i) {
    acc <- interim_accept[,i]
    ess <- sum(n_b[acc]) + sum(n_i[!acc])
    ess * w[i]
  }, numeric(1)) %>%
    sum()
  decisions <- vapply(seq_len(nrow(yi)), function(j) {
    i_accept <- interim_accept[,j]
    if (all(i_accept == FALSE)) {
      out <- i_accept
    } else {
      Pb <- LMEM_ef(n_b[i_accept], y[j, i_accept], p0)
      out <- rep(FALSE, k)
      out[i_accept] <- Pb > ft
    }
    out
  }, logical(k))
  accept_means <- apply(decisions, 1, function(x) sum(x * w))
  basket_power <- rep(NA, length(n_b))
  basket_type1_error <- rep(NA, length(n_b))
  basket_power[promising] <- accept_means[promising]
  basket_type1_error[!promising] <- accept_means[!promising]
  list(basket_power = basket_power,
       basket_type1_error = basket_type1_error,
       n_i = n_i,
       n_b = n_b,
       ESS = ESS)
})
type1_errors <- t(vapply(scenarios_res, function(x) x$basket_type1_error, numeric(k)))
basket_power <- t(vapply(scenarios_res, function(x) x$basket_power, numeric(k)))
ESS <- t(vapply(scenarios_res, function(x) x$ESS, numeric(1)))
LMEM_twostage <- list(scenarios = t(scenarios),
                      type1_errors = type1_errors,
                      basket_power = basket_power,
                      ESS = ESS,
                      it = it,
                      ft = ft,
                      n_i = n_i,
                      n_b = n_b
)
LMEM_twostage

save(LMEM_twostage, file = file.path(getwd(), "output", "LMEM_twostage.RData"))
