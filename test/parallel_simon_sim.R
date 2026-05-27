# Evaluate Parallel-Simon's design (PS) for optimal design comparison

# Load in libraries
library(clinfun)
library(parallel)
mc.cores <- detectCores() - 1

# scenarios to evaluate
scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.20, 0.05, 0.05)
)
# Find optimal design
ph2simon(0.05, 0.20, 0.05, 0.20)
# Set optimal design parameters
B <- 4
n_i <- rep(10, B)
n_b <- rep(29, B)
R1 <- rep(0, B)
r <- rep(3, B)
# number of simulated trials to use in each scenario
n_sim <- 10000
# historical/control response rate
p0 <- 0.05
# number of scenarios
n_s <- length(scenarios)
# number of baskets (equal to B in manuscript)
k <- length(n_b)

# Function that outputs decisions from observed data under PS design
PSRule <- function(n_i, n_b, yi, y, p0) {
  futile <- yi <= R1
  accept <- y > r
  SS <- sum(n_i[futile]) + sum(n_b[!futile])
  decisions <- !futile & accept
  list(decisions = decisions,
       futile = futile,
       SS = SS)
}

# pre-simulate data for each scenario

# interim data
yi_sim_list <- lapply(scenarios, function(scenario) {
  tp <- as.vector(scenario)
  yi_sim <- matrix(rbinom(length(n_i) * n_sim, rep(n_i, times = n_sim),
                          prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
  yi_sim
})
# second-stage data
y_sim_list <- lapply(seq_along(scenarios), function(i) {
  scenario <- scenarios[[i]]
  tp <- as.vector(scenario)
  y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b - n_i, times = n_sim),
                         prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
  y_sim <- yi_sim_list[[i]] + y_sim
  y_sim
})

# evaluation function
evaluatePPScenario <- function(PPRule, yi_sim, y_sim, promising) {
  res <- lapply(seq_len(nrow(yi_sim)),
                function(i) PSRule(n_i, n_b, yi_sim[i,], y_sim[i,], p0 = p0))
  SS <- vapply(res, function(x) x$SS, numeric(1))
  reject <- vapply(res, function(x) x$decisions, logical(ncol(yi_sim)))
  FWER <- vapply(res, function(x) any(x$decisions[!promising]), logical(1))
  list(ESS = mean(SS),
       FWER = mean(FWER),
       reject = rowMeans(reject))
}

# evaluate PS design under each scenario
res <- mclapply(seq_len(n_s), function(i) {
  tp <- as.vector(scenarios[,i])
  promising <- tp > p0
  c(evaluatePPScenario(PSRule, yi_sim_list[[i]], y_sim_list[[i]], promising),
    list(promising = promising))
}, mc.cores = mc.cores)
# extract basket-wise Type I error rates
type1_errors <- t(vapply(res, function(x) {
  out <- rep(NA, k)
  out[!x$promising] <- x$reject[!x$promising]
  out
}, numeric(k)))
# extract basket-wise power
basket_power <- t(vapply(res, function(x) {
  out <- rep(NA, k)
  out[x$promising] <- x$reject[x$promising]
  out
}, numeric(k)))
# extract effective sample size and family-wise error rates
ESS <- vapply(res, function(x) x$ESS, numeric(1))
FWERs <- vapply(res, function(x) x$FWER, numeric(1))
# save in list
PS_res <- list(scenarios = t(scenarios),
               basket_power = basket_power,
               type1_errors = type1_errors,
               ESS = ESS,
               FWERs = FWERs,
               n_i = n_i,
               n_b = n_b,
               R1 = R1,
               r = r,
               n_sim = n_sim)
# save to file
saveRDS(PS_res, file.path(getwd(), "output", "twostage_comparison", "PS_res.rds"))
