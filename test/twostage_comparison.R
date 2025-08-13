library(ggplot2)
library(data.table)
rm(list = ls())
set.seed(1130)

source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "uniform.R"))
source(file.path(getwd(), "src", "pooled.R"))
source(file.path(getwd(), "src", "MEM.R"))
# source(file.path(getwd(), "src", "parallel_simon.R"))

# p <- c(0.05, 0.05, 0.2, 0.2)
# n_b <- c(20, 20, 20, 20)
# promising <- p > 0.05

# compare methods
scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.2, 0.05, 0.05)
)

set.seed(1130)
n_calibrate <- 10000
n_sim <- 5000
p0 <- 0.05
p1 <- 0.20
# single stage evaluation
# y_calibrate <- matrix(rbinom(length(n_b) * n_calibrate, rep(n_b, times = n_calibrate),
#                              prob = rep(p0, times = n_calibrate)), nrow = n_calibrate, byrow = TRUE)
# y_sim_list <- lapply(scenarios, function(scenario) {
#   tp <- as.vector(scenario)
#   y_sim <- matrix(rbinom(length(n_b) * n_sim, rep(n_b, times = n_sim),
#                          prob = rep(tp, times = n_sim)), nrow = n_sim, byrow = TRUE)
# })

alpha2 <- 0.05
beta2 <- 0.20
param_grid <- expand.grid(alpha1 = seq(0.10, 0.45, by = 0.05),
                          beta1 = seq(0.02, 0.15, by = 0.01))
param_grid$ESS <- NA

# Reference method: basket-wise (independent) evaluation
uparam_grid <- param_grid
uparam_grid$ESS <- NA
cat("Calibrating basket-wise two-stage design...")
for(i in seq(1, nrow(param_grid))) {
  alpha1 <- uparam_grid[i,]$alpha1
  beta1 <- uparam_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 UniformBasketEfficacy)
  uparam_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- uparam_grid[which.min(uparam_grid$ESS),]$alpha1
beta1 <- uparam_grid[which.min(uparam_grid$ESS),]$beta1
u_calibrate2 <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = 0.05, beta2 = 0.20, n_sim = n_sim,
                                 UniformBasketEfficacy)
system.time({
  u_res2 <- evaluateTwoStageScenarios(u_calibrate2$n_i,
                                      u_calibrate2$n_b,
                                      scenarios, p0 = 0.05,
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
  calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 PooledBasketEfficacy)
  pparam_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- pparam_grid[which.min(pparam_grid$ESS),]$alpha1
beta1 <- pparam_grid[which.min(pparam_grid$ESS),]$beta1
p_calibrate2 <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                  alpha1 = alpha1, beta1 = beta1,
                                  alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                  PooledBasketEfficacy)
system.time({
  p_res2 <- evaluateTwoStageScenarios(p_calibrate2$n_i,
                                      p_calibrate2$n_b,
                                      scenarios, p0 = 0.05,
                                      p_calibrate2$interim_threshold,
                                      p_calibrate2$pp_threshold,
                                      PooledBasketEfficacy,
                                      n_sim = n_sim)
})
p_res2

# LMEM method
LMEM_ef <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                   d1 = 0, d2 = 3,
                                                   WW_method = FALSE)
LMEM_param_grid <- param_grid
LMEM_param_grid$ESS <- 0
cat("Calibrating LMEM(0, 3) two-stage design...")
for(i in seq(1, nrow(LMEM_param_grid))) {
  alpha1 <- LMEM_param_grid[i,]$alpha1
  beta1 <- LMEM_param_grid[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1,
                                 beta1 = beta1,
                                 alpha2 = alpha2,
                                 beta2 = beta2,
                                 n_sim = n_sim,
                                 LMEM_ef)
  LMEM_param_grid[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- LMEM_param_grid[which.min(LMEM_param_grid$ESS),]$alpha1 # 0.35
beta1 <- LMEM_param_grid[which.min(LMEM_param_grid$ESS),]$beta1 # 0.04
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


# LMEM method 2: more conservative borrowing prior
LMEM_ef2 <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1, b = 1, d1 = 0, d2 = 5,
                                                   WW_method = FALSE)
LMEM_param_grid2 <- param_grid
LMEM_param_grid2$ESS <- 0
cat("Calibrating LMEM(0, 5) two-stage design...")
for(i in seq(1, nrow(LMEM_param_grid2))) {
  alpha1 <- LMEM_param_grid2[i,]$alpha1
  beta1 <- LMEM_param_grid2[i,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 LMEM_ef2)
  LMEM_param_grid2[i,]$ESS <- calibrate$ESS
  #cat("ESS: ", calibrate$ESS, "\n")
}
alpha1 <- LMEM_param_grid2[which.min(LMEM_param_grid2$ESS),]$alpha1 # 0.35
beta1 <- LMEM_param_grid2[which.min(LMEM_param_grid2$ESS),]$beta1 # 0.04
LMEM_calibrate2 <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20, alpha1 = alpha1, beta1 = beta1,
                                    alpha2 = alpha2,
                                    beta2 = beta2,
                                    n_sim = n_sim,
                                    LMEM_ef2)
system.time({
  LMEM_res2 <- evaluateTwoStageScenarios(LMEM_calibrate2$n_i,
                                         LMEM_calibrate2$n_b,
                                         scenarios, p0 = 0.05,
                                         LMEM_calibrate2$interim_threshold,
                                         LMEM_calibrate2$pp_threshold,
                                         LMEM_ef2, n_sim = n_sim)
})
LMEM_res2

cat("All done! Saving results...")
dir.create(file.path(getwd(), "output", "twostage_comparison"), showWarnings = FALSE)
saveRDS(u_res2, file.path(getwd(), "output", "twostage_comparison", "basketwise.rds"))
saveRDS(p_res2, file.path(getwd(), "output", "twostage_comparison", "pooled.rds"))
saveRDS(LMEM_res, file.path(getwd(), "output", "twostage_comparison", "LMEM_0_3.rds"))
saveRDS(LMEM_res2, file.path(getwd(), "output", "twostage_comparison", "LMEM_0_5.rds"))

