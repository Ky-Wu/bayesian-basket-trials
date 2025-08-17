library(data.table)
set.seed(1130)
source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "MEM.R"))

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.2, 0.05, 0.05)
)
files <- dir(file.path(getwd(), "output", "MEM_simulation"))
param_grid <- data.frame()
for (file in files) {
  param_grid <- rbind(param_grid,
                      fread(file.path(getwd(), "output", "MEM_simulation", file)))
}
alpha1 <- param_grid[which.min(param_grid$ESS),]$alpha1
beta1 <- param_grid[which.min(param_grid$ESS),]$beta1
p0 <- 0.05
p1 <- 0.20
alpha2 <- 0.05
beta2 <- 0.20
n_sim <- 200
cal_time <- system.time({
  MEM_calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 MEMEfficacy, increment_start = 13)
})

system.time({
  MEM_res <- evaluateTwoStageScenarios(MEM_calibrate$n_i,
                                      MEM_calibrate$n_b,
                                      scenarios, p0 = 0.05,
                                      MEM_calibrate$interim_threshold,
                                      MEM_calibrate$pp_threshold,
                                      MEMEfficacy, n_sim = n_sim)
})

saveRDS(MEM_res, file.path(getwd(), "output", "twostage_comparison", "MEM_0.1.rds"))
