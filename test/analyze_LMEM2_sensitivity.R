library(data.table)
library(ggplot2)
library(ggpubr)

d_df <- expand.grid(d1 = seq(-4, 4, by = 2),
                    d2 = seq(-2, 2, by = 1))
format_model_name <- function(d1, d2) {
  d1_part <- ifelse(d1 < 0, paste0("neg", abs(d1)), as.character(d1))
  d2_part <- ifelse(d2 < 0, paste0("neg", abs(d2)), as.character(d2))
  paste0("LMEM2_", d1_part, "_", d2_part, ".rds")
}
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

scenarios <- data.frame(
  "Global_Null" = c(0.05, 0.05, 0.05, 0.05),
  "Global_Alternative" = c(0.2, 0.2, 0.2, 0.2),
  "One_in_the_Middle" = c(0.2, 0.2, 0.1, 0.3),
  "Linear" = c(0.05, 0.15, 0.25, 0.35),
  "Good_Nugget" = c(0.2, 0.2, 0.2, 0.05),
  "Bad_Nugget" = c(0.2, 0.05, 0.05, 0.05),
  "Half" = c(0.2, 0.2, 0.05, 0.05)
)

all_results <- lapply(seq_len(nrow(d_df)), function(i) {
  d1 <- d_df[i,]$d1
  d2 <- d_df[i,]$d2
  model_fp <- format_model_name(d1, d2)
  model_name <- paste0("LMEM2(", d1, ", ", d2, ")")
  res <- readRDS(file.path(getwd(), "output", "twostage_comparison", model_fp))
  res_dt <- getResLongDT(model_name, res$basket_power, res$type1_error,
                         res$ESS, res$FWERs, t(scenarios), 0.05, 0.20)
  res_dt[, d1 := d1]
  res_dt[, d2 := d2]
  res_dt[, n_i := res$n_i[1]]
  res_dt[, n_b := res$n_b[1]]
})
all_dt <- rbindlist(all_results)


all_dt[, scenario_label := paste(as.character(n_promising), "Active")]
all_dt[scenario == 3, scenario_label := "One in the Middle"]
all_dt[scenario == 4, scenario_label := "Linear"]
indx <- c(1,2,5,6,7)
all_dt <- all_dt[scenario %in% indx,]
all_dt[, scenario_label := as.character(n_promising)]
all_dt[, active := ifelse(promising, "Active", "Inactive")]
graph_data <- all_dt[, .(accept_prob = mean(accept_prob)),
       by = .(scenario_label, active, method, d1, d2)]
comparison_plot <- ggplot(data = graph_data) +
  geom_line(aes(x = as.numeric(scenario_label),
                y = accept_prob, linetype = factor(d1), color = factor(d2))) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, data = all_dt[active == "Inactive",],
             linewidth = 0.35) +
  theme_bw() +
  labs(x = "Number of Active Baskets", y = "Acceptance Probability") +
  scale_color_discrete(name = "delta_2") +
  scale_linetype_discrete(name = "delta_1") +
  facet_wrap(~active)
comparison_plot

ess_data <- all_dt[, .(ess = mean(ESS)),
                     by = .(scenario_label, method, d1, d2)]
ess_plot <- ggplot(data = ess_data) +
  geom_line(aes(x = as.numeric(scenario_label),
                y = ess, linetype = factor(d1), color = factor(d2))) +
  scale_color_discrete(name = "delta_2") +
  scale_linetype_discrete(name = "delta_1") +
  theme_bw() +
  labs(x = "Number of Active Baskets", y = "Expected Sample Size")
ess_plot

trial_size_data <- melt(all_dt, measure.vars = c("n_i", "n_b"),
                        variable.name = "sample_size_type",
                        value.name = "sample_size")
trial_size_data <- unique(trial_size_data[, .(d1, d2, sample_size_type, sample_size)])
trial_size_plot <- ggplot(data = trial_size_data) +
  geom_col(aes(x = d2, y = sample_size),
           position = "dodge", fill = "dodgerblue") +
  facet_grid(
    rows = vars(sample_size_type),
    cols = vars(d1),
    labeller = labeller(
      d1 = function(x) paste0("delta_1: ", x),
      sample_size_type = as_labeller(c(
        "n_i"  = "First Stage",
        "n_b" = "Second Stage"
      ))
    )
  ) +
  labs(x = "delta_2", y = "Sample Size Per Basket") +
  theme_bw()
trial_size_plot

ggsave(file.path(getwd(), "output", "LMEM2_sensitivity_analysis", "optimal_oc_plot.png"),
       comparison_plot,
       dpi = 500, width = 8, height = 6)
ggsave(file.path(getwd(), "output", "LMEM2_sensitivity_analysis", "ess_plot.png"),
       ess_plot,
       dpi = 500, width = 8, height = 6)

p <- ggarrange(comparison_plot, ess_plot, nrow = 1, common.legend = TRUE, legend = "right")
p

ggsave(file.path(getwd(), "output", "LMEM2_sensitivity_analysis", "LMEM2_optimal_comparison_plot.png"),
       p,
       dpi = 500, width = 8, height = 6)
### Pre-specified stage size comparison ###

all_results <- lapply(seq_len(nrow(d_df)), function(i) {
  d1 <- d_df[i,]$d1
  d2 <- d_df[i,]$d2
  model_fp <- format_model_name(d1, d2)
  model_name <- paste0("LMEM2(", d1, ", ", d2, ")")
  res <- readRDS(file.path(getwd(), "output", "fixedss_sim", model_fp))
  res_dt <- getResLongDT(model_name, res$basket_power, res$type1_error,
                         res$ESS, res$FWERs, t(scenarios), 0.05, 0.20)
  res_dt[, d1 := d1]
  res_dt[, d2 := d2]
})
all_dt <- rbindlist(all_results)


all_dt[, scenario_label := paste(as.character(n_promising), "Active")]
all_dt[scenario == 3, scenario_label := "One in the Middle"]
all_dt[scenario == 4, scenario_label := "Linear"]
indx <- c(1,2,5,6,7)
all_dt <- all_dt[scenario %in% indx,]
all_dt[, scenario_label := as.character(n_promising)]
all_dt[, active := ifelse(promising, "Active", "Inactive")]
graph_data <- all_dt[, .(accept_prob = mean(accept_prob)),
                     by = .(scenario_label, active, method, d1, d2)]
comparison_plot <- ggplot(data = graph_data) +
  geom_line(aes(x = as.numeric(scenario_label),
                y = accept_prob, linetype = factor(d1), color = factor(d2))) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, data = all_dt[active == "Inactive",],
             linewidth = 0.35) +
  theme_bw() +
  labs(x = "Number of Active Baskets", y = "Acceptance Probability") +
  scale_color_discrete(name = "delta_2") +
  scale_linetype_discrete(name = "delta_1") +
  facet_wrap(~active)
comparison_plot

ggsave(file.path(getwd(), "output", "LMEM2_sensitivity_analysis", "LMEM2_eqss_oc.png"),
       comparison_plot,
       dpi = 500, width = 8, height = 6)
