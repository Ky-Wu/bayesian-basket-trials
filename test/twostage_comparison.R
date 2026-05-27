# Evaluation of optimal two-stage Bayesian and frequentist basket trial designs
# (optimal meaning that they minimize expected sample size under the global null)

# load libraries
library(data.table)
library(partitions)
rm(list = ls())
set.seed(1130)

# load helper functions
source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "LMEM2.R"))
source(file.path(getwd(), "src", "uniform.R"))
source(file.path(getwd(), "src", "pooled.R"))
source(file.path(getwd(), "src", "MEM.R"))

# Due to computational complexity of MEM method, evaluation of MEM
# designs is handled in a separate file, MEM_simulation.R

# state scenarios with different response rates for comparison
scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.2, 0.05, 0.05)
)
# set seed for
set.seed(1130)
# number of simulations to use to evaluate each design
n_sim <- 2000
# historical/control response rate
p0 <- 0.05
# minimal clinically meaningful response rate
p1 <- 0.20
# second-stage basket-wise Type I error (under global null) constraint
alpha2 <- 0.05
# second-stage basket-wise Type II error (under global alternative) constraint
beta2 <- 0.20
# Optimization grid: possible first-stage Type I/II error constraints
param_grid <- expand.grid(alpha1 = seq(0.10, 0.45, by = 0.05),
                          beta1 = seq(0.02, 0.15, by = 0.01))
param_grid$ESS <- NA

# Reference method: basket-wise (independent) evaluation
# Evaluation of each method follows this general structure:
# 1. Copy optimization grid
uparam_grid <- param_grid
# ESS: expected sample size under global null: p_i = p_0 for all i
uparam_grid$ESS <- NA
cat("Calibrating basket-wise two-stage design...")
# 2. Evaluate grid of designs
for(i in seq(1, nrow(param_grid))) {
  alpha1 <- uparam_grid[i,]$alpha1
  beta1 <- uparam_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 UniformBasketEfficacy)
  uparam_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
# 3. Extract optimal design parameters
alpha1 <- uparam_grid[which.min(uparam_grid$ESS),]$alpha1
beta1 <- uparam_grid[which.min(uparam_grid$ESS),]$beta1
# 4. Re-calibrate with optimal design parameters
u_calibrate2 <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 UniformBasketEfficacy)
# 5. Evaluate optimal design under all designs
system.time({
  u_res2 <- evaluateTwoStageScenarios(u_calibrate2$n_i,
                                      u_calibrate2$n_b,
                                      scenarios, p0 = p0,
                                      u_calibrate2$interim_threshold,
                                      u_calibrate2$pp_threshold,
                                      UniformBasketEfficacy, n_sim = n_sim)
})
u_res2

# Reference Method 2: Pooled evaluation

pparam_grid <- param_grid
pparam_grid$ESS <- NA
cat("Calibrating pooled two-stage design...")
for(i in seq(1, nrow(param_grid))) {
  alpha1 <- uparam_grid[i,]$alpha1
  beta1 <- uparam_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 PooledBasketEfficacy)
  pparam_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- pparam_grid[which.min(pparam_grid$ESS),]$alpha1
beta1 <- pparam_grid[which.min(pparam_grid$ESS),]$beta1
p_calibrate2 <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                  alpha1 = alpha1, beta1 = beta1,
                                  alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                  PooledBasketEfficacy)
system.time({
  p_res2 <- evaluateTwoStageScenarios(p_calibrate2$n_i,
                                      p_calibrate2$n_b,
                                      scenarios, p0 = p0,
                                      p_calibrate2$interim_threshold,
                                      p_calibrate2$pp_threshold,
                                      PooledBasketEfficacy,
                                      n_sim = n_sim)
})
p_res2

# LMEM method with delta = 0
# Prior: p_i ~ Beta(a, b) for i = 1, ..., B
LMEM_ef <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                   d1 = 0, d2 = 0)
LMEM_param_grid <- param_grid
LMEM_param_grid$ESS <- 0
cat("Calibrating LMEM(0) two-stage design...")
for(i in seq(1, nrow(LMEM_param_grid))) {
  alpha1 <- LMEM_param_grid[i,]$alpha1
  beta1 <- LMEM_param_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0,
                                 p1 = p1,
                                 alpha1 = alpha1,
                                 beta1 = beta1,
                                 alpha2 = alpha2,
                                 beta2 = beta2,
                                 n_sim = n_sim,
                                 LMEM_ef)
  LMEM_param_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- LMEM_param_grid[which.min(LMEM_param_grid$ESS),]$alpha1
beta1 <- LMEM_param_grid[which.min(LMEM_param_grid$ESS),]$beta1
LMEM_calibrate <- calibrateTwoStage(B = 4, p0, p1, alpha1 = alpha1, beta1 = beta1,
                                    alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                    LMEM_ef)
system.time({
  LMEM_res <- evaluateTwoStageScenarios(LMEM_calibrate$n_i,
                                        LMEM_calibrate$n_b,
                                        scenarios, p0,
                                        LMEM_calibrate$interim_threshold,
                                        LMEM_calibrate$pp_threshold,
                                        LMEM_ef, n_sim = n_sim)
})
LMEM_res


# LMEM2: substantial borrowing prior with delta1 = -4, delta2 = -2
LMEM2_ef <- function(n_b, y, p0) LMEM2BasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                    d1 = -4, d2 = -2)
LMEM2_param_grid <- param_grid
LMEM2_param_grid$ESS <- 0
cat("Calibrating LMEM(-4, -2) two-stage design...")
for(i in seq(1, nrow(LMEM2_param_grid))) {
  alpha1 <- LMEM2_param_grid[i,]$alpha1
  beta1 <- LMEM2_param_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 LMEM2_ef)
  LMEM2_param_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- LMEM2_param_grid[which.min(LMEM_param_grid$ESS),]$alpha1 # 0.35
beta1 <- LMEM2_param_grid[which.min(LMEM_param_grid$ESS),]$beta1 # 0.04
LMEM2_calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1, alpha1 = alpha1, beta1 = beta1,
                                    alpha2 = alpha2,
                                    beta2 = beta2,
                                    n_sim = n_sim,
                                    LMEM2_ef)
system.time({
  LMEM2_res <- evaluateTwoStageScenarios(LMEM2_calibrate$n_i,
                                         LMEM2_calibrate$n_b,
                                         scenarios, p0 = p0,
                                         LMEM2_calibrate$interim_threshold,
                                         LMEM2_calibrate$pp_threshold,
                                         LMEM2_ef, n_sim = n_sim)
})
LMEM2_res



# LMEM2: conservative borrowing prior with delta1 = 0, delta2 = 0.5
LMEM2_ef2 <- function(n_b, y, p0) LMEM2BasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                      d1 = 0, d2 = 0.5)
LMEM2_param_grid2 <- param_grid
LMEM2_param_grid2$ESS <- 0
cat("Calibrating LMEM2(0, 0.5) two-stage design...")
for(i in seq(1, nrow(LMEM2_param_grid2))) {
  alpha1 <- LMEM2_param_grid2[i,]$alpha1
  beta1 <- LMEM2_param_grid2[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 LMEM2_ef2)
  LMEM2_param_grid2[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- LMEM2_param_grid2[which.min(LMEM2_param_grid2$ESS),]$alpha1 # 0.35
beta1 <- LMEM2_param_grid2[which.min(LMEM2_param_grid2$ESS),]$beta1 # 0.04
LMEM2_calibrate2 <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1, alpha1 = alpha1, beta1 = beta1,
                                     alpha2 = alpha2,
                                     beta2 = beta2,
                                     n_sim = n_sim,
                                     LMEM2_ef2)
system.time({
  LMEM2_res2 <- evaluateTwoStageScenarios(LMEM2_calibrate2$n_i,
                                         LMEM2_calibrate2$n_b,
                                         scenarios, p0 = p0,
                                         LMEM2_calibrate2$interim_threshold,
                                         LMEM2_calibrate2$pp_threshold,
                                         LMEM2_ef2, n_sim = n_sim)
})
LMEM2_res2

# Save all results
cat("All done! Saving results...")
dir.create(file.path(getwd(), "output", "twostage_comparison"), showWarnings = FALSE)
saveRDS(u_res2, file.path(getwd(), "output", "twostage_comparison", "basketwise.rds"))
saveRDS(p_res2, file.path(getwd(), "output", "twostage_comparison", "pooled.rds"))
saveRDS(LMEM_res, file.path(getwd(), "output", "twostage_comparison", "LMEM_0.rds"))
saveRDS(LMEM2_res, file.path(getwd(), "output", "twostage_comparison", "LMEM2_neg4_neg2.rds"))
saveRDS(LMEM2_res2, file.path(getwd(), "output", "twostage_comparison", "LMEM2_0_half.rds"))
