# load libraries
library(data.table)
library(partitions)
rm(list = ls())
set.seed(1130)

# load command line args:
# sim_i: index of param grid to evaluate
# This is meant for a script that loops through the indices (1:25)
# of the param grid for parallel evaluation.
args <- commandArgs(trailingOnly = TRUE)
sim_i <- as.integer(args[1])

# load helper functions
source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "LMEM2.R"))

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

# evaluate different possible values for d1 and d2
d_df <- expand.grid(d1 = seq(-4, 4, by = 2),
                    d2 = seq(-2, 2, by = 1))
format_model_name <- function(d1, d2) {
  d1_part <- ifelse(d1 < 0, paste0("neg", abs(d1)), as.character(d1))
  d2_part <- ifelse(d2 < 0, paste0("neg", abs(d2)), as.character(d2))
  paste0("LMEM2_", d1_part, "_", d2_part, ".rds")
}

d1 <- d_df[sim_i,]$d1
d2 <- d_df[sim_i,]$d2
print(paste0("Calibrating LMEM2(", d1, ", ", d2, "):"))
LMEM2_ef <- function(n_b, y, p0) LMEM2BasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                     d1 = d1, d2 = d2)
LMEM2_param_grid <- param_grid
LMEM2_param_grid$ESS <- 0
for(j in seq(1, nrow(LMEM2_param_grid))) {
  alpha1 <- LMEM2_param_grid[j,]$alpha1
  beta1 <- LMEM2_param_grid[j,]$beta1
  #cat("alpha1: ", alpha1, "| beta1: ", beta1, "\n")
  calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                 alpha1 = alpha1, beta1 = beta1,
                                 alpha2 = alpha2, beta2 = beta2, n_sim = n_sim,
                                 LMEM2_ef, verbose = FALSE)
  LMEM2_param_grid[j,]$ESS <- calibrate$ESS
}
alpha1 <- LMEM2_param_grid[which.min(LMEM2_param_grid$ESS),]$alpha1
beta1 <- LMEM2_param_grid[which.min(LMEM2_param_grid$ESS),]$beta1
LMEM2_calibrate <- calibrateTwoStage(B = 4, p0 = p0, p1 = p1,
                                     alpha1 = alpha1,
                                     beta1 = beta1,
                                     alpha2 = alpha2,
                                     beta2 = beta2,
                                     n_sim = n_sim,
                                     LMEM2_ef)
LMEM2_res <- evaluateTwoStageScenarios(LMEM2_calibrate$n_i,
                                       LMEM2_calibrate$n_b,
                                       scenarios, p0 = p0,
                                       LMEM2_calibrate$interim_threshold,
                                       LMEM2_calibrate$pp_threshold,
                                       LMEM2_ef, n_sim = n_sim)
model_fp <- format_model_name(d1, d2)
saveRDS(LMEM2_res, file.path(getwd(), "output", "twostage_comparison", model_fp))
