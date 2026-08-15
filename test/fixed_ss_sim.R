# Fixed sample size evaluation of Bayesian and frequentist two-stage basket trial designs

# load in libraries
library(ggplot2)
library(data.table)
rm(list = ls())
set.seed(1130)

# load helper functions
source(file.path(getwd(), "src", "LMEM.R"))
source(file.path(getwd(), "src", "LMEM2.R"))
source(file.path(getwd(), "src", "uniform.R"))
source(file.path(getwd(), "src", "sim_helper.R"))
source(file.path(getwd(), "src", "pooled.R"))
#source(file.path(getwd(), "src", "twostage_prunepool.R"))
source(file.path(getwd(), "src", "MEM.R"))

# state scenarios with different response rates for comparison
scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.20, 0.05, 0.05)
)

# Simulation settings same as MEM optimal design
# number of patients per basket in interim stage
n_i <- c(4, 4, 4, 4)
# number of patients per basket total in second stage
# (including first-stage accrual)
n_b <- c(19, 19, 19, 19)
# historical/control response rate
p0 <- 0.05
# minimal clinically meaningful response rate
p1 <- 0.20
# second-stage basket-wise type I error constraint under global null
alpha2 <- 0.05
# interim threshold (reject/stop if total Stage I < R1, or <= r)
R1 <- 2
r <- R1 - 1

#' Fixed Sample Size Two-Stage Design Effective Sample Size (ESS) Evaluation
#'
#' @param p0 Historical/control response rate.
#' @param p1 Minimal clinically meaningful response rate.
#' @param n_i Vector of interim stage sample sizes.
#' @param n_b Vector of second-stage cumulative sample sizes.
#' @param r Interim stopping rule. Stop for futility if <= r responses are
#' observed in total across all functions.
#' @param efficacyFunction Function with args n_b, y, p_0 that outputs posterior efficacy probabilities.
#' @param alpha2 Second-stage basket-wise type I error constraint under global null.
#' @param n_sim Number of simulated trials for evaluation of each scenario and design combination.
#'
#' @return List containing metrics, most notable posterior predictive probability threshold $f^{\star}$.
#' @export
#'
#' @examples
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
  # obtain threshold that controls type I error while maximizing power
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


# LMEM2(-4, -2) efficacy function
LMEM2_ef <- function(n_b, y, p0) LMEM2BasketEfficacy(n_b, y, p0,
                                                     a = 1, b = 1,
                                                     d1 = -4, d2 = -2)
# obtain posterior probability threshold for efficacy decision
LMEM2_setting <- computeESSFixedSS(p0, p1, n_i, n_b, r,
                                   LMEM2_ef, alpha2, n_sim = 4000)
# use obtained threshold for evaluation of each scenario
LMEM2_res <- evaluateTwoStageScenarios(n_i, n_b,
                                      scenarios, p0, r,
                                      LMEM2_setting$pp_threshold,
                                      LMEM2_ef, n_sim = 4000)

# LMEM2(0, 0.5)
LMEM2_ef2 <- function(n_b, y, p0) LMEM2BasketEfficacy(n_b, y, p0,
                                                      a = 1, b = 1,
                                                      d1 = 0, d2 = 0.5)
LMEM2_setting2 <- computeESSFixedSS(p0, p1, n_i, n_b, r,
                                   LMEM2_ef2, alpha2, n_sim = 4000)
LMEM2_res2 <- evaluateTwoStageScenarios(n_i, n_b,
                                       scenarios, p0, r,
                                       LMEM2_setting2$pp_threshold,
                                       LMEM2_ef2, n_sim = 4000)

# LMEM(0)
LMEM_ef <- function(n_b, y, p0) LMEMBasketEfficacy(n_b, y, p0, a = 1, b = 1,
                                                   d1 = 0, d2 = 0)

LMEM_setting <- computeESSFixedSS(p0, p1, n_i, n_b, r,
                                  LMEM_ef, alpha2, n_sim = 4000)
LMEM_res <- evaluateTwoStageScenarios(n_i,
                                      n_b,
                                      scenarios, p0, r,
                                      LMEM_setting$pp_threshold,
                                      LMEM_ef, n_sim = 4000)

# Pooled analysis
pooled_setting <- computeESSFixedSS(p0, p1, n_i, n_b, r,
                                    PooledBasketEfficacy, alpha2, n_sim = 20000)
pooled_res <- evaluateTwoStageScenarios(n_i, n_b,
                                        scenarios, p0, r,
                                        pooled_setting$pp_threshold,
                                        PooledBasketEfficacy, n_sim = 20000)

# Uniform analysis: aggregated futility analysis and independent final analysis
uniform_setting <- computeESSFixedSS(p0, p1, n_i, n_b, r,
                                     UniformBasketEfficacy, alpha2 = alpha2,
                                     n_sim = 20000)

uniform_res <- evaluateTwoStageScenarios(n_i, n_b,
                                         scenarios, p0, r,
                                         uniform_setting$pp_threshold,
                                         UniformBasketEfficacy,
                                         n_sim = 20000)

# MEM analysis: resuses design from optimal design simulation
# load in saved results
# prunepool design found in src/prunepool/prunepool_optimize.R
MEM_res <- readRDS(file.path(getwd(), "output", "twostage_comparison", "MEM_0.1.rds"))
prunepool_res <- readRDS(file.path(getwd(), "output", "fixedss_sim", "prunepool_res.rds"))


# formatting function
getResLongDT <- function(method_name, power_tab, type1_errors,
                         sc, p0 = 0.05, p1 = 0.2) {
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


indx <- c(1,2,5,6,7)
# format all results into data tables
LMEM_dt <- getResLongDT("LMEM(0)", LMEM_res$basket_power, LMEM_res$type1_error,
                        t(scenarios), 0.05, 0.2)
PP_dt <- getResLongDT("Prune-pool", prunepool_res$basket_power, prunepool_res$type1_error,
                      t(scenarios), 0.05, 0.2)
u_dt <- getResLongDT("Uniform", uniform_res$basket_power, uniform_res$type1_error,
                     t(scenarios), 0.05, 0.2)
p_dt <- getResLongDT("Pooled", pooled_res$basket_power, pooled_res$type1_error,
                     t(scenarios), 0.05, 0.2)
MEM_dt <- getResLongDT("MEM(0.1)", MEM_res$basket_power, MEM_res$type1_error,
                       t(scenarios), 0.05, 0.2)
LMEM2_dt <- getResLongDT("LMEM2(-4, -2)", LMEM2_res$basket_power, LMEM2_res$type1_error,
                         t(scenarios), 0.05, 0.2)
LMEM2_dt2 <- getResLongDT("LMEM2(0, 0.5)", LMEM2_res2$basket_power, LMEM2_res2$type1_error,
                         t(scenarios), 0.05, 0.2)
all_dt <- rbindlist(list(LMEM2_dt, LMEM2_dt2, LMEM_dt, MEM_dt, PP_dt, u_dt, p_dt))
#all_dt <- all_dt[scenario %in% indx]
all_dt[, scenario_label := as.character(n_promising)]
all_dt[scenario == 3, scenario_label := "One in the Middle"]
all_dt[scenario == 4, scenario_label := "Linear"]
all_dt[, active := ifelse(promising, "Active", "Inactive")]
all_dt <- all_dt[, .(accept_prob = mean(accept_prob)),
                 by = .(scenario, active, method, scenario_label)]
# plot basket Type I and power metrics across 0-5 active scenarios
comparison_plot <- ggplot(data = all_dt[scenario %in% indx,]) +
  geom_jitter(aes(x = scenario_label, y = accept_prob, shape = method, color = method),
              height = 0, width = 0.2, alpha = 0.8, size = 3) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, data = all_dt[active == "Inactive",],
             linewidth = 0.35) +
  theme_bw() +
  facet_wrap(~active) +
  labs(x = "Number of Active Baskets", y = "Acceptance Probability") +
  scale_shape_manual(values = c(0, 2, 4, 7, 16, 17, 18)) +
  scale_color_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67",
                                "#00A9FF", "#C77CFF", "#FF61CC"))

comparison_plot <- ggplot(data = all_dt[scenario %in% indx,]) +
  geom_point(aes(x = as.numeric(scenario_label), y = accept_prob, shape = method, color = method),
             alpha = 0.8) +
  geom_line(aes(x = as.numeric(scenario_label), y = accept_prob, linetype = method, color = method),
            alpha = 0.8) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, data = all_dt[active == "Inactive",],
             linewidth = 0.35) +
  theme_bw() +
  labs(x = "Number of Active Baskets", y = "Acceptance Probability") +
  facet_wrap(~active) +
  scale_shape_manual(values = c(0, 2, 4, 7, 8, 16, 17, 18)) +
  scale_color_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67",
                                "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"))
comparison_plot
# save plot
ggsave(file.path(getwd(), "output", "eqss_comparison_plot.png"),
       comparison_plot,
       dpi = 500, width = 8, height = 6)

# save posterior probability thresholds in data table
pp_rules <- data.table(
  method = c("LMEM(0)", "LMEM2(-4, -2)", "LMEM2(0, 0.5)", "MEM(0.1)", "Uniform", "Pooled"),
  #interim_ss = sum(n_i),
  #total_ss = sum(n_b),
  pp_threshold = c(signif(LMEM_res$pp_threshold, 5),
                   signif(LMEM2_res$pp_threshold, 5),
                   signif(LMEM2_res2$pp_threshold, 5),
                   signif(MEM_res$pp_threshold, 5),
                   signif(uniform_res$pp_threshold, 5),
                   signif(pooled_res$pp_threshold, 5))
)
#pp_rules$interim_ss <- as.integer(pp_rules$interim_ss)
#pp_rules$interim_threshold <- as.integer(pp_rules$interim_threshold)
#pp_rules$total_ss <- as.integer(pp_rules$total_ss)
colnames(pp_rules) <- c("Method", #"Interim SS", "Total Trial Size",
                        "Posterior Probability Threshold")
print(xtable::xtable(pp_rules, caption = "Minimum posterior probabilities to declare efficacy in a basket under
                     each design in the fixed sample size design. The minimum posterior probability for each framework was found by maximizing basket-wise power while controlling
                     the basket-wise Type I error rate at 5\\% under the global null.
                     To align with the optimal MEM design from the optimal design simulation study,
                     each design employs a first-stage sample size of 4 patients per basket,
                     a total sample size of 19 patients per basket, and stops for futility if there are one or zero total responses in the first stage.",
                     label = "tab:eqss_design_pps", align = rep("c", ncol(pp_rules) + 1), digits = 3),
      type = "latex", include.rownames = FALSE,
      file.path(getwd(), "output", "eqss_design_pps.tex"))
