library(data.table)
rm(list = ls())
set.seed(1130)

source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "MEM.R"))

args <- commandArgs(trailingOnly = TRUE)
sim_i <- as.integer(args[1])

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.2, 0.05, 0.05)
)

n_sim <- 200
p0 <- 0.05
p1 <- 0.20
alpha2 <- 0.05
beta2 <- 0.20
param_grid <- expand.grid(alpha1 = seq(0.10, 0.45, by = 0.05),
                          beta1 = seq(0.02, 0.15, by = 0.01))
# MEM method
MEM_param_grid <- param_grid
MEM_param_grid$ESS <- 0
cat("Calibrating MEM two-stage design...")
alpha1 <- MEM_param_grid[sim_i,]$alpha1
beta1 <- MEM_param_grid[sim_i,]$beta1
#cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
cal_time <- system.time({
  calibrate <- calibrateTwoStage(B = 4, p0 = 0.05, p1 = 0.20,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 MEMEfficacy)
})

res <- data.frame(alpha1 = alpha1, beta1 = beta1,
                  ESS = calibrate$ESS, cal_time = cal_time[[3]])

name <- paste0("MEM_sim", sim_i, ".csv")
dir.create(file.path(getwd(), "output", "MEM_simulation"), showWarnings = FALSE)
fwrite(res, file.path(getwd(), "output", "MEM_simulation", name))
