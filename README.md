# Optimal Local Exchangeability Two-Stage Designs for Basket Trials with Aggregated Futility Analysis

This repository contains code implementation for the simulation comparison in "Optimal Local Exchangeability Two-Stage Designs for Basket Trials with Aggregated Futility Analysis."

Simulation files:

- Optimal design comparison: `test/twostage_comparison.R`
   - MEM simulations: `test/MEM_simulation.R`
- Pre-specified stage size design comparison: `test/fixed_ss_sim.R`
   - MEM simulations: `test/MEM_fixed_ss_sim.R`
   - For manuscript results, ran in parallel on Hoffman2 Cluster with command arg 1 to 112 for each parameter combination
- LMEM2 Sensitivity Analysis: `test/LMEM2_sensitivity_analysis.R` (optimal sample size design) and `test/LMEM2_sensitivity_analysis_fixedss.R`
  - Scripts take command arg for row of hyperparameter grid (defined in variable named `d_df`, by default 25 rows)
  - Example run for first hyperparameter value: `Rscript test/LMEM2_sensitivity_analysis.R 1`
