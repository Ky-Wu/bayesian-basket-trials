library(poisbinom)
library(parallel)
mc.cores <- detectCores() - 1

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.20, 0.05, 0.05)
)
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
        decisions[!prune] <- (sum(y_notpruned) >= R2)
      }
    }
    list(decisions = decisions,
         prune = prune,
         SS = SS)
  }
}

PPRule <- constructPPRule(3, 2, 0.018)

yi_sim_list <- lapply(scenarios, function(scenario) {
  tp <- as.vector(scenario)
  yi_sim <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                          prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
})
y_sim_list <- lapply(seq_along(scenarios), function(i) {
  scenario <- scenarios[[i]]
  tp <- as.vector(scenario)
  y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b - n_i, times = n_sim),
                         prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
  y_sim <- yi_sim_list[[i]] + y_sim
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

res <- mclapply(seq_len(n_s), function(i) {
  tp <- as.vector(scenarios[,i])
  promising <- tp > p0
  c(evaluatePPScenario(PPRule, yi_sim_list[[i]], y_sim_list[[i]], promising),
    list(promising = promising))
}, mc.cores = mc.cores)
type1_errors <- t(vapply(res, function(x) {
  out <- rep(NA, k)
  out[!x$promising] <- x$reject[!x$promising]
  out
}, numeric(k)))
basket_power <- t(vapply(res, function(x) {
  out <- rep(NA, k)
  out[x$promising] <- x$reject[x$promising]
  out
}, numeric(k)))
ESS <- vapply(res, function(x) x$ESS, numeric(1))
FWERs <- vapply(res, function(x) x$FWER, numeric(1))
prunepool_res <- list(scenarios = t(scenarios),
     basket_power = basket_power,
     type1_errors = type1_errors,
     ESS = ESS,
     FWER = FWERs,
     n_i = n_i,
     n_b = n_b,
     R1 = 3,
     n_sim = n_sim)

saveRDS(prunepool_res, file.path(getwd(), "output", "twostage_comparison", "prunepool_res.rds"))
