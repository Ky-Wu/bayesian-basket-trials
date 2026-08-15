library(data.table)

#' Efficiently Compute Basketwise Acceptance Rates for a Homogeneous Basket Trial
#'
#' @param K Number of tumor indications (baskets)
#' @param n1 Stage I sample size per basket
#' @param N Total sample size per basket (Stage I + Stage II)
#' @param p0 Null response rate for all baskets
#' @param p True response rate for all baskets
#' @param r Threshold of total responses for a basket to be pooled (not pruned), note: this is separate from r as defined in general framework
#' @param R1 Overall Stage I efficacy boundary (reject/stop if total Stage I < R1)
#' @param alpha2_star Significance level used to determine Stage II boundaries R2
#' @return The exact basket-wise acceptance rate
compute_basket_acceptance <- function(K, n1, N, p0, p, r, R1, alpha2_star) {

  # Stage I and ESS calculation
  prob_pass_stage1 <- 1 - pbinom(R1 - 1, size = K * n1, prob = p)
  expected_N <- (K * n1) + (prob_pass_stage1 * K * (N - n1))

  # Stage II boundaries
  R2_vec <- numeric(K)
  for (M in 1:K) {
    R2_vec[M] <- qbinom(1 - alpha2_star, size = M * N, prob = p0)
  }

  # Step 1: Base Joint PMF
  dt_base <- as.data.table(expand.grid(x1 = 0:n1, x2 = 0:(N - n1)))
  dt_base[, prob := dbinom(x1, size = n1, prob = p) * dbinom(x2, size = N - n1, prob = p)]

  dt_base[, z := x1]
  dt_base[, v := as.integer((x1 + x2) >= r)]
  dt_base[, w := (x1 + x2) * v]

  # Group and initialize log_prob (matching your original logic)
  dt_total <- dt_base[, .(prob = sum(prob)), by = .(z, v, w)]
  dt_total[, log_prob := log(prob)]

  # EXPLICIT ISOLATION: Copy and rename columns to avoid any inner-join collisions
  dt_single <- copy(dt_total)
  setnames(dt_single,
           old = c("z", "v", "w", "prob", "log_prob"),
           new = c("z_new", "v_new", "w_new", "prob_new", "log_prob_new"))

  # Add dummy key for forced Cartesian join
  dt_total[, dummy := 1L]
  dt_single[, dummy := 1L]

  # Step 2: Convolutions
  if (K > 1) {
    for (k in 2:K) {

      # Forced cross-join on the single shared dummy column
      merged <- merge(dt_total, dt_single, by = "dummy", allow.cartesian = TRUE)

      # Compute new states explicitly into temporary columns
      merged[, z_sum := z + z_new]
      merged[, v_sum := v + v_new]
      merged[, w_sum := w + w_new]

      # Use log_prob addition to prevent floating point underflow
      merged[, log_prob_sum := log_prob + log_prob_new]
      merged[, prob_prod := exp(log_prob_sum)]

      # Aggregate and overwrite dt_total with standard names
      dt_total <- merged[, .(prob = sum(prob_prod)), by = .(z = z_sum, v = v_sum, w = w_sum)]

      # Re-initialize variables for the next iteration
      dt_total[, log_prob := log(prob)]
      dt_total[, dummy := 1L]
    }
  }

  # Step 3: Evaluate rejection rules
  # Filter only successful trials (passed R1 and at least 1 basket pooled)
  dt_success <- dt_total[z >= R1 & v > 0]

  if (nrow(dt_success) > 0) {
    # Map the correct R2 boundary for each row's 'v' (number of pooled baskets)
    dt_success[, R2_thresh := R2_vec[v]]

    # Keep only outcomes that pass the Stage II boundary
    dt_success <- dt_success[w >= R2_thresh]

    global_acceptance <- sum(dt_success$prob)
    basket_acceptance <- sum(dt_success$prob * dt_success$v) / K
  } else {
    global_acceptance <- 0
    basket_acceptance <- 0
  }

  return(list(
    Global_acceptance = global_acceptance,
    Basket_acceptance = basket_acceptance,
    Expected_Sample_Size = expected_N
  ))
}

# --- Example Usage ---

alpha <- compute_basket_acceptance(
  K = 4, n1 = 9, N = 21, p0 = 0.05, p = 0.05, r = 2, R1 = 3, alpha2_star = 0.020
)

power <- compute_basket_acceptance(
  K = 4, n1 = 9, N = 21, p0 = 0.05, p = 0.20, r = 2, R1 = 3, alpha2_star = 0.020
)
print(paste("Computed Global Type I Error:", alpha$Global_acceptance))
print(paste("Computed Basket-wise Type I Error:", alpha$Basket_acceptance))

#' Find Boundaries of alpha2_star
#'
#' Calculates the exact values of alpha2_star at which the critical
#' boundaries R2_vec will change for any number of pooled baskets M.
#'
#' @param K Maximum number of pooled baskets
#' @param N Total sample size per basket
#' @param p0 Null response rate
#' @param max_alpha Optional limit to only return boundaries below a practical threshold (e.g., 0.20)
#' @return A sorted numeric vector of alpha2_star boundary values
find_alpha2_star_boundaries <- function(K, N, p0, max_alpha = 1.0) {

  all_boundaries <- c()

  # Iterate over all possible pooling scenarios (M)
  for (M in 1:K) {
    # Generate all possible numbers of successes for M * N patients
    x <- 0:(M * N)

    # Calculate the upper tail probability (1 - CDF)
    # These are the exact points where qbinom steps to the next integer
    alpha_bounds_M <- 1 - pbinom(x, size = M * N, prob = p0)

    # Collect boundaries
    all_boundaries <- c(all_boundaries, alpha_bounds_M)
  }

  # Due to floating point arithmetic, very close numbers might not be treated
  # as duplicates. Rounding to a high precision (e.g., 10 decimal places)
  # ensures we get true mathematical unique boundaries.
  all_boundaries <- round(all_boundaries, digits = 10)
  unique_boundaries <- sort(unique(all_boundaries))

  # Filter to only return realistic alpha levels (alpha > 0 and alpha <= max_alpha)
  unique_boundaries <- unique_boundaries[unique_boundaries > 0 & unique_boundaries <= max_alpha]

  return(unique_boundaries)
}

# --- Example Usage ---
# K=4 baskets, N=21 total patients per basket, p0=0.05
boundaries <- find_alpha2_star_boundaries(K = 4, N = 21, p0 = 0.05, max_alpha = 1)

print("Boundary values of alpha2_star where R2_vec changes:")
print(boundaries)

#' Optimize Trial Design to Minimize ESS under the Null Hypothesis
#'
#' @param K Number of tumor indications (baskets)
#' @param p0 Null response rate (used for Type I error and ESS null)
#' @param p1 Alternative response rate (used for Power)
#' @param min_power Lower bound on global power
#' @param max_type1 Upper bound on global Type I error
#' @param N_min Minimum total sample size per basket to search
#' @param N_max Maximum total sample size per basket to search
#' @param max_alpha Upper bound for alpha2_star search space
#' @return A data frame with the optimal configuration, or NULL if no design is feasible
optimize_basket_design <- function(K, p0, p1, min_power, max_type1,
                                   N_min = 10, N_max = 30, max_alpha = 0.20,
                                   n1_min = 5, n1_max = NULL, R1_max = 5) {

  cat("Generating search space...\n")

  # 1. Build the grid of all possible (N, n1, r, R1) combinations
  designs <- list()
  for (N in N_min:N_max) {
    actual_n1_min <- max(1, n1_min)
    actual_n1_max <- if (is.null(n1_max)) (N - 1) else min(n1_max, N - 1)
    if (actual_n1_min > N - 1) next

    for (n1 in actual_n1_min:actual_n1_max) {
      actual_R1_max <- K * n1
      if (!is.null(R1_max)) {
        actual_R1_max <- min(actual_R1_max, R1_max)
      }
      for (R1 in 1:actual_R1_max) {

        # Compute ESS under the null hypothesis (p0)
        prob_pass_null <- 1 - pbinom(R1 - 1, size = K * n1, prob = p0)
        ess_null <- (K * n1) + (prob_pass_null * K * (N - n1))

        for (r in 1:N) {
          designs[[length(designs) + 1]] <- c(N = N, n1 = n1, R1 = R1, r = r, ESS = ess_null)
        }
      }
    }
  }

  search_space <- as.data.frame(do.call(rbind, designs))

  # 2. Sort the search space by ESS (ascending)
  search_space <- search_space[order(search_space$ESS), ]
  total_configs <- nrow(search_space)
  cat(sprintf("Evaluating %d combinations ordered by minimum ESS...\n", total_configs))

  # 3. Iterate through candidate designs
  for (i in 1:total_configs) {

    cfg <- search_space[i, ]

    # Optional: Print progress every 500 configurations evaluated
    if (i %% 2 == 0) cat(sprintf("Evaluated %d / %d...\n", i, total_configs))

    # Get unique boundaries for alpha2_star for this specific N
    boundaries <- find_alpha2_star_boundaries(K, cfg$N, p0, max_alpha)
    if (length(boundaries) == 0) next

    # Generate midpoints to test each distinct R2 step region
    intervals <- c(0, boundaries)
    test_alphas <- (intervals[1:(length(intervals)-1)] + intervals[2:length(intervals)]) / 2

    # Test each alpha2_star boundary state
    for (a_star in test_alphas) {

      # Evaluate Type I Error first (often fails constraint, saving computation time)
      null_res <- compute_basket_acceptance(
        K = K, n1 = cfg$n1, N = cfg$N, p0 = p0, p = p0,
        r = cfg$r, R1 = cfg$R1, alpha2_star = a_star
      )

      if (null_res$Global_acceptance > max_type1) next # Fails Type I Error bound

      # If Type I error passes, evaluate Power
      alt_res <- compute_basket_acceptance(
        K = K, n1 = cfg$n1, N = cfg$N, p0 = p0, p = p1,
        r = cfg$r, R1 = cfg$R1, alpha2_star = a_star
      )

      if (alt_res$Global_acceptance >= min_power) {
        # Both constraints met! Because we sorted by ESS, this is mathematically the optimal design.
        cat("\n--- Optimal Design Found! ---\n")

        result <- data.frame(
          K = K,
          N = cfg$N,
          n1 = cfg$n1,
          r = cfg$r,
          R1 = cfg$R1,
          alpha2_star = a_star,
          ESS_Null = cfg$ESS,
          Global_Type1 = null_res$Global_acceptance,
          Global_Power = alt_res$Global_acceptance,
          Basket_Type1 = null_res$Basket_acceptance,
          Basket_Power = alt_res$Basket_acceptance
        )
        return(result)
      }
    }
  }

  cat("\nNo feasible design found within the specified N limits and constraints.\n")
  return(NULL)
}

# --- Optimal design simulation study ---
# Seeking a design for K=4 with at least 80% power and at most 5% Type I error
# N limited between 15 and 25 to restrict search space time for the example
optimal_design <- optimize_basket_design(
  K = 4,
  p0 = 0.05,
  p1 = 0.20,
  min_power = 0.80,
  max_type1 = 0.05,
  N_min = 21,
  N_max = 30,
  n1_min = 5,
  R1_max = 5
)

# K = 4, N = 21, n1 = 7, r = 1, R1 = 4, alpha2_star = 0.0003560104
print(optimal_design)
#K  N n1 r R1  alpha2_star ESS_Null Global_Type1 Global_Power Basket_Type1 Basket_Power
#4  21 7 1  4 0.0003560104 30.74814  0.001432363    0.8130395  0.001168068    0.8075914

# --- Fixed sample size design ---
# Seeking a design for K=4 with at most 5% basket-wise Type I error under global null
# Criterion: maximize power


#' Maximize Power for a Fixed (K, n1, N, R1) Design
#'
#' Sweeps over pruning threshold r and alpha2_star, holding K, n1, N, R1 fixed,
#' and returns the combination with the highest Global_Power subject to
#' Global_Type1 <= max_type1.
#'
#' @param K Number of baskets
#' @param n1 Fixed Stage I sample size per basket
#' @param N Fixed total sample size per basket
#' @param R1 Fixed Stage I efficacy boundary
#' @param p0 Null response rate
#' @param p1 Alternative response rate
#' @param max_type1 Upper bound on global Type I error
#' @param max_alpha Upper bound for alpha2_star search space
#' @return A data frame of the best design found (max power), or NULL if none feasible
maximize_power_fixed_design <- function(K, r_max, n1, N, R1, p0, p1, max_type1, max_alpha = 0.20) {

  boundaries <- find_alpha2_star_boundaries(K, N, p0, max_alpha)
  if (length(boundaries) == 0) {
    cat("No alpha2_star boundaries found in the given range.\n")
    return(NULL)
  }
  intervals <- c(0, boundaries)
  test_alphas <- (intervals[1:(length(intervals) - 1)] + intervals[2:length(intervals)]) / 2

  best_result <- NULL
  best_power <- -Inf

  total_alpha <- length(test_alphas)
  cat(sprintf("Scanning %d values of r x %d values of alpha2_star = %d combinations...\n",
              r_max, total_alpha, r_max * total_alpha))

  for (r in 1:r_max) {
    for (a_star in test_alphas) {

      null_res <- compute_basket_acceptance(
        K = K, n1 = n1, N = N, p0 = p0, p = p0,
        r = r, R1 = R1, alpha2_star = a_star
      )
      if (null_res$Global_acceptance > max_type1) next  # Type I error constraint violated

      alt_res <- compute_basket_acceptance(
        K = K, n1 = n1, N = N, p0 = p0, p = p1,
        r = r, R1 = R1, alpha2_star = a_star
      )

      if (alt_res$Global_acceptance > best_power) {
        best_power <- alt_res$Global_acceptance
        best_result <- data.frame(
          K = K,
          N = N,
          n1 = n1,
          r = r,
          R1 = R1,
          alpha2_star = a_star,
          Global_Type1 = null_res$Global_acceptance,
          Global_Power = alt_res$Global_acceptance,
          Basket_Type1 = null_res$Basket_acceptance,
          Basket_Power = alt_res$Basket_acceptance,
          Expected_N = alt_res$Expected_Sample_Size
        )
      }
    }
  }

  if (is.null(best_result)) {
    cat("\nNo feasible design found (Type I error constraint never satisfied).\n")
  } else {
    cat("\n--- Best Power-Maximizing Design Found! ---\n")
  }

  return(best_result)
}

# --- Run for fixed design: K = 4, n1 = 10, N = 29, R1 = 2 ---
best_fixedss_design <- maximize_power_fixed_design(
  K = 4,
  r_max = 5,
  n1 = 4,
  N = 19,
  R1 = 2,
  p0 = 0.05,
  p1 = 0.20,
  max_type1 = 0.05
)

print(best_fixedss_design)
#K  N n1 r R1 alpha2_star Global_Type1 Global_Power Basket_Type1 Basket_Power Expected_N
#4 19  4 1  2  0.01829471   0.04173792    0.8567951   0.03161046    0.8465016   67.55575
