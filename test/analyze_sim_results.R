getResLongDT <- function(method_name, power_tab, type1_errors, ESS, FWERs, sc, p0 = 0.05, p1 = 0.2) {
  A <- power_tab
  A[is.na(A)] <- type1_errors[is.na(A)]
  data <- data.table()
  for (i in 1:nrow(A)) {
    promising <- sc[i,] > p0
    for (j in 1:ncol(A)) {
      newdat <- data.table(scenario = i, basket = j, promising = promising[j],
                           n_promising = sum(promising),
                           accept_prob = A[i, j],
                           ESS = ESS[i],
                           FWER = FWERs[i],
                           method = method_name)
      data <- rbind(data, newdat)
    }
  }
  data
}

LMEM_dt <- getResLongDT("LMEM(0, 3)", LMEM_res$basket_power, LMEM_res$type1_error,
                        LMEM_res$ESS, LMEM_res$FWERs,
                        t(scenarios), 0.05, 0.2)
LMEM_dt2 <- getResLongDT("LMEM(0, 5)", LMEM_res2$basket_power, LMEM_res2$type1_error,
                         LMEM_res2$ESS, LMEM_res2$FWERs,
                         t(scenarios), 0.05, 0.2)
PP_dt <- getResLongDT("Prune-pool", prunepool_res$basket_power, prunepool_res$type1_error,
                      prunepool_res$ESS, prunepool_res$FWER,
                      t(scenarios), 0.05, 0.2)
u_dt <- getResLongDT("Uniform", u_res2$basket_power, u_res2$type1_error,
                     u_res2$ESS, u_res2$FWERs,
                     t(scenarios),0.05, 0.2)
p_dt <- getResLongDT("Pooled", p_res2$basket_power, p_res2$type1_error,
                     p_res2$ESS, p_res2$FWERs,
                     t(scenarios),0.05, 0.2)
ps_dt <- getResLongDT("Parallel Simon", PS_res$basket_power, PS_res$type1_error,
                      PS_res$ESS, PS_res$FWERs,
                      t(scenarios),0.05, 0.2)
all_dt <- rbindlist(list(LMEM_dt, LMEM_dt2, PP_dt, u_dt, p_dt, ps_dt))
#all_dt <- all_dt[scenario %in% indx]
all_dt[, scenario_label := paste(as.character(n_promising), "Active")]
all_dt[scenario == 3, scenario_label := "One in the Middle"]
all_dt[scenario == 4, scenario_label := "Linear"]
comparison_plot <- ggplot(data = all_dt) +
  geom_jitter(aes(x = scenario_label, y = accept_prob, shape = promising, color = method),
              height = 0, width = 0.15, alpha = 0.8, size = 2) +
  theme_bw() +
  labs(x = "Scenario", y = "Acceptance Probability") +
  scale_shape_manual(values = c(0, 17))

ess_plot <- ggplot(data = all_dt[basket == 1,]) +
  geom_col(aes(x = scenario_label, y = ESS, fill = method), color = "black", alpha = 0.8,
           position = "dodge") +
  theme_bw() +
  coord_cartesian(ylim = c(0, max(all_dt$ESS)))+
  labs(x = "Scenario", y = "Expected Total Sample Size")
fwer_plot <- ggplot(data = all_dt[basket == 1 & scenario %in% c(1,4,5,6,7),]) +
  geom_col(aes(x = scenario_label, y = FWER, fill = method), color = "black", alpha = 0.8,
           position = "dodge") +
  theme_bw() +
  coord_cartesian(ylim = c(0, max(all_dt$FWER))) +
  labs(x = "Scenario", y = "Family-wise Type I Error Rate")
#test6.Rdata: alpha2 = 0.04, beta2 = 0.8
# LMEM has similar power to uniform analysis and lower ESS but higher type I errors in mixed
# max FWER: 28% in 2 active
# prune-pool has most power but similar ESS to uniform and significantly inflated type I error in mixed scenarios
# prune-pool has FWER 43% in 1 or 2 active

ggsave(file.path(getwd(), "output", "min_ess_designs.png"),
       comparison_plot, dpi = 500, width = 8, height = 6)
ggsave(file.path(getwd(), "output", "ess_plot.png"), ess_plot,
       dpi = 500, width = 8, height = 6)
ggsave(file.path(getwd(), "output", "fwer_plot.png"), fwer_plot,
       dpi = 500, width = 8, height = 6)


sc <- as.data.table(t(scenarios[c(1, 6, 7, 5, 2, 4, 3)]))
sc <- cbind(data.table(Scenario = c("0 Active (Global Null)", "1 Active", "2 Active",
                                    "3 Active", "4 Active (Global Alternative)",
                                    "Linear", "One in the Middle")), sc)
colnames(sc) <- c("Scenario", "Basket 1", "Basket 2", "Basket 3", "Basket 4")

library(xtable)

print(xtable::xtable(sc, caption = "Response rates of simulated responses under global and mixed simulation scenarios.",
                     label = "tab:sim_scenarios"),
      type = "latex", include.rownames = FALSE,
      file.path(getwd(), "output", "sim_scenarios.tex"))

rules <- data.table(
  method = c("Uniform", "Pooled", "Parallel Simon", "LMEM(0, 3)", "LMEM(0, 5)", "Prune-pool"),
  interim_ss = as.integer(c(sum(u_res2$n_i), sum(p_res2$n_i), sum(PS_res$n_i),
                            sum(LMEM_res$n_i), sum(LMEM_res2$n_i), sum(prunepool_res$n_i))),
  aggregated_futility = c("Yes", "Yes", "No", "Yes", "Yes", "Yes"),
  interim_threshold = as.integer(c(u_res2$interim_threshold, p_res2$interim_threshold,
                                   PS_res$R1[1], LMEM_res$interim_threshold, LMEM_res2$interim_threshold,
                                   prunepool_res$R1)),
  total_ss = as.integer(c(sum(u_res2$n_b), sum(p_res2$n_b), sum(PS_res$n_b),
                          sum(LMEM_res$n_b), sum(LMEM_res2$n_b), sum(prunepool_res$n_b))),
  pp_threshold = c(signif(u_res2$pp_threshold, 3),
                   signif(p_res2$pp_threshold, 3),
                   "-",
                   signif(LMEM_res$pp_threshold, 3),
                   signif(LMEM_res2$pp_threshold, 3),
                   "-")
  # prune_threshold = c("-", "-", 5, "-", "-", 1),
  #  pooled_threshold = c("-", "-", "-", "-", "-", "Q(M; 21, 0.018)")
)
rules
colnames(rules) <- c("Method", "Interim SS", "Aggregated Futility",
                     "Interim Threshold", "Total Sample Size",
                     "PP Threshold")
print(xtable::xtable(rules, caption = "Sample sizes and posterior probability thresholds of optimal two-stage designs controlling global Type II error rate at 0.20 and global Type I error rate at 0.04.
                     Interim threshold denotes the number of responses to be exceeded in order for the trial to proceed to the second stage.",
                     label = "tab:optimal_designs", align = rep("c", ncol(rules) + 1)),
      type = "latex", include.rownames = FALSE,
      file.path(getwd(), "output", "optimal_designs.tex"))
