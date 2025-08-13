library(bhmbasket)

EXNEXEfficacy <- function(n_b, y, p0) {
  trial <- bhmbasket::createTrial(n_b, y)
  analysis <- performAnalyses(trial, target_rates = rep(p0, length(n_b)))
  decisions <- bhmbasket::getGoDecisions(analysis,
                            evidence_levels = rep(0.5, length(n_b)),
                            cohort_names = paste("p", seq_along(n_b), sep = "_"),
                            boundary_rules = quote(c(x[1] > 0.10, x[2] > 0.10,
                                                     x[3] > 0.10, x[4] > 0.10)),
                            overall_min_gos = 1)
}
