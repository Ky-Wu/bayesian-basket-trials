library(ggplot2)
library(data.table)
rm(list = ls())

source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "uniform.R"))
source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "pooled.R"))
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
LMEM_ef <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1/2, b = 1/2,
                                                   d1 = 0, d2 = 3,
                                                   WW_method = FALSE)

for (i in 1:nrow(param_grid)) {
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           LMEM_ef, alpha2, n_sim = 4000)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
#LMEM_setting <- param_grid[ESS == min(ESS[power >= 1 - beta2]),]
LMEM_setting <- param_grid
LMEM_res <- evaluateTwoStageScenarios(n_i,
                                      n_b,
                                      scenarios, p0,
                                      LMEM_setting$r,
                                      LMEM_setting$pp_threshold,
                                      LMEM_ef, n_sim = 4000)

for (i in 1:nrow(param_grid)) {
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           PooledBasketEfficacy, alpha2, n_sim = 20000)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
pooled_setting <- param_grid[ESS == min(ESS[power >= 1 - beta2]),]
pooled_res <- evaluateTwoStageScenarios(n_i, n_b,
                                        scenarios, p0,
                                        pooled_setting$r,
                                        pooled_setting$pp_threshold,
                                        PooledBasketEfficacy, n_sim = 20000)
for (i in 1:nrow(param_grid)) {
  print(i)
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           UniformBasketEfficacy, alpha2 = 0.09,
                           n_sim = 20000)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
uniform_setting <- param_grid[ESS == min(ESS[power == max(power)]),]

uniform_res <- evaluateTwoStageScenarios(n_i, n_b,
                                         scenarios, p0,
                                         uniform_setting$r,
                                         uniform_setting$pp_threshold,
                                         UniformBasketEfficacy,
                                         n_sim = 20000)

LMEM_ef2 <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1/2, b = 1/2,
                                                    d1 = 0, d2 = 5,
                                                    WW_method = FALSE)
for (i in 1:nrow(param_grid)) {
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           LMEM_ef2, alpha2, n_sim = 4000)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
#LMEM_setting <- param_grid[ESS == min(ESS[power >= 1 - beta2]),]
LMEM_setting2 <- param_grid
LMEM_res2 <- evaluateTwoStageScenarios(n_i,
                                       n_b,
                                       scenarios, p0,
                                       LMEM_setting2$r,
                                       LMEM_setting2$pp_threshold,
                                       LMEM_ef2, n_sim = 4000)

LMEM_ef3 <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1/2, b = 1/2,
                                                    d1 = 5, d2 = 5,
                                                    WW_method = FALSE)
for (i in 1:nrow(param_grid)) {
  res <- computeESSFixedSS(p0, p1, n_i, n_b, param_grid[i,]$r,
                           LMEM_ef3, alpha2, n_sim = 4000)
  param_grid[i,]$ESS <- res$ESS
  param_grid[i,]$power <- res$HA_power
  param_grid[i,]$type1_error <- res$alpha2_estimate
  param_grid[i,]$pp_threshold <- res$pp_threshold
}
#LMEM_setting <- param_grid[ESS == min(ESS[power >= 1 - beta2]),]
LMEM_setting3 <- param_grid
LMEM_res3 <- evaluateTwoStageScenarios(n_i,
                                       n_b,
                                       scenarios, p0,
                                       LMEM_setting3$r,
                                       LMEM_setting3$pp_threshold,
                                       LMEM_ef3, n_sim = 4000)
### MEM ###

MEM_res <- readRDS(file.path(getwd(), "output", "fixedss_sim", "MEM_res.rds"))


getResLongDT <- function(method_name, power_tab, type1_errors, sc, p0 = 0.05, p1 = 0.2) {
  A <- power_tab
  A[is.na(A)] <- type1_errors[is.na(A)]
  data <- data.table()
  for (i in 1:nrow(A)) {
    promising <- sc[i,] > p0
    for (j in 1:ncol(A)) {
      newdat <- data.table(scenario = i, basket = j, promising = promising[j],
                           n_promising = sum(promising),
                           accept_prob = A[i, j],
                           method = method_name)
      data <- rbind(data, newdat)
    }
  }
  data
}

source(file.path(getwd(), "test", "prunepool_sim.R"))

indx <- c(1,2,5,6,7)
LMEM_dt <- getResLongDT("LMEM(0, 3)", LMEM_res$basket_power, LMEM_res$type1_error,
                        t(scenarios), 0.05, 0.2)
LMEM_dt2 <- getResLongDT("LMEM(0, 5)", LMEM_res2$basket_power, LMEM_res2$type1_error,
                         t(scenarios), 0.05, 0.2)
LMEM_dt3 <- getResLongDT("LMEM(5, 5)", LMEM_res3$basket_power, LMEM_res3$type1_error,
                         t(scenarios), 0.05, 0.2)
PP_dt <- getResLongDT("Prune-pool", prunepool_res$basket_power, prunepool_res$type1_error,
                      t(scenarios), 0.05, 0.2)
u_dt <- getResLongDT("Uniform", uniform_res$basket_power, uniform_res$type1_error,
                     t(scenarios), 0.05, 0.2)
p_dt <- getResLongDT("Pooled", pooled_res$basket_power, pooled_res$type1_error,
                     t(scenarios), 0.05, 0.2)
MEM_dt <- getResLongDT("MEM(0.1)", MEM_res$basket_power, MEM_res$type1_error,
                       t(scenarios), 0.05, 0.2)
all_dt <- rbindlist(list(LMEM_dt, LMEM_dt2, LMEM_dt3, PP_dt, u_dt, p_dt, MEM_dt))
#all_dt <- all_dt[scenario %in% indx]
all_dt[, scenario_label := as.character(n_promising)]
all_dt[scenario == 3, scenario_label := "One in the Middle"]
all_dt[scenario == 4, scenario_label := "Linear"]
all_dt[, active := ifelse(promising, "Active", "Inactive")]
comparison_plot <- ggplot(data = all_dt[scenario %in% indx,]) +
  geom_jitter(aes(x = scenario_label, y = accept_prob, shape = method, color = method),
              height = 0, width = 0.2, alpha = 0.8, size = 3) +
  theme_bw() +
  facet_wrap(~active) +
  labs(x = "Scenario", y = "Acceptance Probability") +
  scale_shape_manual(values = c(0, 2, 4, 7, 8, 16, 17))
comparison_plot

ggsave(file.path(getwd(), "output", "eqss_comparison_plot.png"),
       comparison_plot,
       dpi = 500, width = 8, height = 6)


pp_rules <- data.table(
  method = c("LMEM(0, 3)", "LMEM(0, 5)", "LMEM(5, 5)", "MEM(0.1)", "Uniform", "Pooled"),
  interim_ss = sum(n_i),
  interim_threshold = c(2, 2, 2, 2, 2, 2),
  #interim_pp_threshold = c("-", "-", "-", signif(LMEM_twostage$it, 3), "-", "-"),
  total_ss = sum(n_b),
  pp_threshold = c(signif(LMEM_res$pp_threshold, 5),
                   signif(LMEM_res2$pp_threshold, 5),
                   signif(LMEM_res3$pp_threshold, 5),
                   signif(MEM_res$pp_threshold, 5),
                   signif(uniform_res$pp_threshold, 5),
                   signif(pooled_res$pp_threshold, 5))
)
pp_rules$interim_ss <- as.integer(pp_rules$interim_ss)
pp_rules$interim_threshold <- as.integer(pp_rules$interim_threshold)
pp_rules$total_ss <- as.integer(pp_rules$total_ss)
colnames(pp_rules) <- c("Method", "Interim SS", "Interim Threshold", "Total Trial Size", "PP Threshold")
print(xtable::xtable(pp_rules, caption = "Minimum posterior probabilities to declare efficacy in a basket under
                     each design with fixed interim sample size (9 per basket), interim threshold (2 total responses), and
                     second-stage cumulative sample size (21 per basket). The Type I error rate is controlled at 7\\% under the global null.",
                     label = "tab:eqss_design_pps", align = rep("c", ncol(pp_rules) + 1), digits = 3),
      type = "latex", include.rownames = FALSE,
      file.path(getwd(), "output", "eqss_design_pps.tex"))
