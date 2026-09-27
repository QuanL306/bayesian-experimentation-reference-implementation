#!/usr/bin/env Rscript
# =============================================================================
# Multi-Armed Bandit Simulation Capstone Project
# Comprehensive comparison of bandit algorithms across diverse scenarios
# =============================================================================

set.seed(42)

library(ggplot2)
library(gridExtra)
library(grid)

cat("============================================================\n")
cat("  Multi-Armed Bandit Simulation Study\n")
cat("  Algorithms: EG, UCB1, UCB-Tuned, TS, EXP3, + non-stationary variants\n")
cat("  Scenarios: 5 (standard, hard, many-arms, non-stationary, normal)\n")
cat("============================================================\n\n")

# =============================================================================
# SECTION 1: ALGORITHM IMPLEMENTATIONS
# =============================================================================

# --- Epsilon-Greedy ---
epsilon_greedy <- function(arms, T, epsilon = 0.1, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  counts <- integer(K)
  sums <- numeric(K)

  for (t in seq_len(T)) {
    if (t <= K) {
      a <- t
    } else if (runif(1) < epsilon) {
      a <- sample.int(K, 1)
    } else {
      a <- which.max(sums / pmax(counts, 1))
    }
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    counts[a] <- counts[a] + 1L
    sums[a] <- sums[a] + r
    arm_means_history[t, ] <- sums / pmax(counts, 1)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- UCB1 ---
ucb1 <- function(arms, T, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  counts <- integer(K)
  sums <- numeric(K)

  for (t in seq_len(T)) {
    if (t <= K) {
      a <- t
    } else {
      ucb <- sums / pmax(counts, 1) + sqrt(2 * log(t) / pmax(counts, 1))
      a <- which.max(ucb)
    }
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    counts[a] <- counts[a] + 1L
    sums[a] <- sums[a] + r
    arm_means_history[t, ] <- sums / pmax(counts, 1)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- UCB-Tuned (variance-adaptive) ---
ucb_tuned <- function(arms, T, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  counts <- integer(K)
  sums <- numeric(K)
  sum_sq <- numeric(K)

  for (t in seq_len(T)) {
    if (t <= K) {
      a <- t
    } else {
      means <- sums / pmax(counts, 1)
      variances <- pmax(
        sum_sq / pmax(counts, 1) - means^2 + sqrt(2 * log(t) / pmax(counts, 1)),
        0
      )
      ucb <- means + sqrt(log(t) / pmax(counts, 1) * pmin(variances, 1/4))
      a <- which.max(ucb)
    }
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    counts[a] <- counts[a] + 1L
    sums[a] <- sums[a] + r
    sum_sq[a] <- sum_sq[a] + r^2
    arm_means_history[t, ] <- sums / pmax(counts, 1)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- Thompson Sampling (Beta-Binomial) ---
thompson_sampling <- function(arms, T, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  alpha <- rep(1, K)
  beta <- rep(1, K)

  for (t in seq_len(T)) {
    theta <- rbeta(K, alpha, beta)
    a <- which.max(theta)
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    alpha[a] <- alpha[a] + r
    beta[a] <- beta[a] + (1 - r)
    arm_means_history[t, ] <- alpha / (alpha + beta)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- Thompson Sampling Normal-Normal conjugate ---
thompson_sampling_normal <- function(arms_means, arms_sds, T) {
  K <- length(arms_means)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  # Prior: mu ~ N(0, 100), known variance approximated from first few samples
  # Use conjugate update: posterior mean and variance
  prior_mu <- rep(0, K)
  prior_var <- rep(100, K)
  # Track sufficient stats
  post_mu <- prior_mu
  post_var <- prior_var
  # We'll estimate variance online
  count <- integer(K)
  sum_x <- numeric(K)
  sum_x2 <- numeric(K)

  for (t in seq_len(T)) {
    # Estimate variance per arm (use sample variance, fallback to 1)
    est_var <- ifelse(count >= 2,
                      (sum_x2 - sum_x^2 / pmax(count, 1)) / pmax(count - 1, 1),
                      1)
    est_var <- pmax(est_var, 0.01)
    # Sample from posterior
    theta <- rnorm(K, post_mu, sqrt(post_var))
    a <- which.max(theta)
    r <- rnorm(1, arms_means[a], arms_sds[a])
    actions[t] <- a
    rewards[t] <- r
    count[a] <- count[a] + 1L
    sum_x[a] <- sum_x[a] + r
    sum_x2[a] <- sum_x2[a] + r^2
    # Conjugate update
    obs_var <- est_var[a]
    post_var[a] <- 1 / (1/prior_var[a] + count[a]/obs_var)
    post_mu[a] <- post_var[a] * (prior_mu[a]/prior_var[a] + sum_x[a]/obs_var)
    arm_means_history[t, ] <- post_mu
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- EXP3 (adversarial bandit) ---
exp3 <- function(arms, T, gamma = 0.1, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  weights <- rep(1, K)
  cum_weights <- numeric(K)  # for tracking

  for (t in seq_len(T)) {
    probs <- (1 - gamma) * weights / sum(weights) + gamma / K
    a <- sample.int(K, 1, prob = probs)
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    # Importance-weighted reward estimate
    est_r <- r / probs[a]
    weights[a] <- weights[a] * exp(gamma * est_r / K)
    # Normalize to prevent overflow
    weights <- weights / sum(weights)
    cum_weights <- cum_weights + weights
    arm_means_history[t, ] <- cum_weights / t
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- Sliding Window Thompson Sampling ---
sw_thompson_sampling <- function(arms, T, window = 500, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  # Store recent rewards per arm
  recent_rewards <- vector("list", K)
  for (k in seq_len(K)) recent_rewards[[k]] <- numeric(0)
  recent_actions <- integer(0)
  recent_reward_vals <- numeric(0)

  for (t in seq_len(T)) {
    # Compute alpha/beta from recent window
    alpha <- rep(1, K)
    beta_param <- rep(1, K)
    start_idx <- max(1, t - window + 1)
    if (t > 1 && start_idx < t) {
      for (idx in start_idx:(t-1)) {
        a_idx <- actions[idx]
        r_idx <- rewards[idx]
        alpha[a_idx] <- alpha[a_idx] + r_idx
        beta_param[a_idx] <- beta_param[a_idx] + (1 - r_idx)
      }
    }
    theta <- rbeta(K, alpha, beta_param)
    a <- which.max(theta)
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    arm_means_history[t, ] <- alpha / (alpha + beta_param)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# --- Discounted UCB ---
discounted_ucb <- function(arms, T, discount = 0.98, reward_fn = NULL) {
  K <- length(arms)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  disc_counts <- numeric(K)
  disc_sums <- numeric(K)

  for (t in seq_len(T)) {
    # Discount all counts
    disc_counts <- disc_counts * discount
    disc_sums <- disc_sums * discount

    if (t <= K) {
      a <- t
    } else {
      n_total <- sum(disc_counts)
      ucb <- disc_sums / pmax(disc_counts, 0.01) +
        sqrt(2 * log(n_total) / pmax(disc_counts, 0.01))
      a <- which.max(ucb)
    }
    r <- if (is.null(reward_fn)) rbinom(1, 1, arms[a]) else reward_fn(a, t)
    actions[t] <- a
    rewards[t] <- r
    disc_counts[a] <- disc_counts[a] + 1
    disc_sums[a] <- disc_sums[a] + r
    arm_means_history[t, ] <- disc_sums / pmax(disc_counts, 0.01)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}

# UCB1 for Normal rewards (no binary assumption)
ucb1_normal <- function(arms_means, arms_sds, T) {
  K <- length(arms_means)
  actions <- integer(T)
  rewards <- numeric(T)
  arm_means_history <- matrix(0, nrow = T, ncol = K)
  counts <- integer(K)
  sums <- numeric(K)

  for (t in seq_len(T)) {
    if (t <= K) {
      a <- t
    } else {
      ucb <- sums / pmax(counts, 1) + sqrt(2 * log(t) / pmax(counts, 1))
      a <- which.max(ucb)
    }
    r <- rnorm(1, arms_means[a], arms_sds[a])
    actions[t] <- a
    rewards[t] <- r
    counts[a] <- counts[a] + 1L
    sums[a] <- sums[a] + r
    arm_means_history[t, ] <- sums / pmax(counts, 1)
  }
  list(actions = actions, rewards = rewards, arm_means_history = arm_means_history)
}


# =============================================================================
# SECTION 2: SIMULATION RUNNER
# =============================================================================

run_simulation <- function(algo_fn, arms, T, n_sims, algo_name, best_arm_rate, ...) {
  cum_regret_all <- matrix(0, nrow = T, ncol = n_sims)
  simple_regret_all <- numeric(n_sims)
  best_arm_pct_all <- numeric(n_sims)
  switching_all <- numeric(n_sims)

  for (i in seq_len(n_sims)) {
    if (i %% 50 == 0) cat(sprintf("  %s sim %d/%d\n", algo_name, i, n_sims))
    result <- algo_fn(arms, T, ...)

    # Cumulative regret
    cum_optimal <- cumsum(rep(best_arm_rate, T))
    cum_actual <- cumsum(result$rewards)
    cum_regret_all[, i] <- cum_optimal - cum_actual

    # Simple regret at end
    simple_regret_all[i] <- best_arm_rate - max(result$arm_means_history[T, ])

    # % best arm pulls
    best_arm <- which.max(arms)
    best_arm_pct_all[i] <- mean(result$actions == best_arm) * 100

    # Switching rate
    switches <- sum(diff(result$actions) != 0)
    switching_all[i] <- switches / (T - 1) * 100
  }

  list(
    algo_name = algo_name,
    cum_regret = cum_regret_all,
    simple_regret = simple_regret_all,
    best_arm_pct = best_arm_pct_all,
    switching_rate = switching_all
  )
}

# Runner for Normal-reward scenarios (different regret calc)
run_simulation_normal <- function(algo_fn, arms_means, arms_sds, T, n_sims, algo_name, ...) {
  best_rate <- max(arms_means)
  cum_regret_all <- matrix(0, nrow = T, ncol = n_sims)
  simple_regret_all <- numeric(n_sims)
  best_arm_pct_all <- numeric(n_sims)
  switching_all <- numeric(n_sims)

  for (i in seq_len(n_sims)) {
    if (i %% 50 == 0) cat(sprintf("  %s sim %d/%d\n", algo_name, i, n_sims))
    result <- algo_fn(arms_means, arms_sds, T, ...)

    cum_optimal <- cumsum(rep(best_rate, T))
    cum_actual <- cumsum(result$rewards)
    cum_regret_all[, i] <- cum_optimal - cum_actual

    simple_regret_all[i] <- best_rate - max(result$arm_means_history[T, ])
    best_arm <- which.max(arms_means)
    best_arm_pct_all[i] <- mean(result$actions == best_arm) * 100

    switches <- sum(diff(result$actions) != 0)
    switching_all[i] <- switches / (T - 1) * 100
  }

  list(
    algo_name = algo_name,
    cum_regret = cum_regret_all,
    simple_regret = simple_regret_all,
    best_arm_pct = best_arm_pct_all,
    switching_rate = switching_all
  )
}

# Runner for non-stationary scenario (regret computed per-segment)
run_simulation_nonstat <- function(algo_fn, T, n_sims, algo_name, rate_fn, ...) {
  cum_regret_all <- matrix(0, nrow = T, ncol = n_sims)
  switching_all <- numeric(n_sims)
  # Track arm pulls over time for adaptation plot
  arm_pulls_matrix <- matrix(0, nrow = T, ncol = 5)  # 5 arms

  for (i in seq_len(n_sims)) {
    if (i %% 50 == 0) cat(sprintf("  %s sim %d/%d\n", algo_name, i, n_sims))
    # Create a reward function that uses time-varying rates
    reward_fn_local <- function(a, t) {
      rates <- rate_fn(t)
      rbinom(1, 1, rates[a])
    }
    arms_placeholder <- rate_fn(1)  # just for K
    result <- algo_fn(arms_placeholder, T, reward_fn = reward_fn_local, ...)

    # Compute regret against oracle that knows current best
    for (t in seq_len(T)) {
      current_rates <- rate_fn(t)
      if (t == 1) {
        cum_regret_all[t, i] <- max(current_rates) - result$rewards[t]
      } else {
        cum_regret_all[t, i] <- cum_regret_all[t-1, i] + max(current_rates) - result$rewards[t]
      }
    }

    switches <- sum(diff(result$actions) != 0)
    switching_all[i] <- switches / (T - 1) * 100

    # Aggregate arm pulls
    for (a in 1:5) {
      arm_pulls_matrix[, a] <- arm_pulls_matrix[, a] +
        cumsum(result$actions == a)
    }
  }

  list(
    algo_name = algo_name,
    cum_regret = cum_regret_all,
    switching_rate = switching_all,
    arm_pulls = arm_pulls_matrix / n_sims
  )
}


# =============================================================================
# SECTION 3: HELPER FOR SUMMARY AND PLOTTING
# =============================================================================

summary_stats <- function(results) {
  data.frame(
    Algorithm = sapply(results, function(x) x$algo_name),
    Mean_Cum_Regret = sapply(results, function(x) mean(x$cum_regret[nrow(x$cum_regret), ])),
    SD_Cum_Regret = sapply(results, function(x) sd(x$cum_regret[nrow(x$cum_regret), ])),
    Mean_Simple_Regret = sapply(results, function(x) mean(x$simple_regret)),
    Mean_BestArm_Pct = sapply(results, function(x) mean(x$best_arm_pct)),
    Mean_Switching = sapply(results, function(x) mean(x$switching_rate)),
    stringsAsFactors = FALSE
  )
}

make_regret_df <- function(results, T) {
  dfs <- lapply(results, function(res) {
    means <- rowMeans(res$cum_regret)
    sds <- apply(res$cum_regret, 1, sd)
    n <- ncol(res$cum_regret)
    se <- sds / sqrt(n)
    data.frame(
      t = seq_len(T),
      mean_regret = means,
      lower = means - 1.96 * se,
      upper = means + 1.96 * se,
      Algorithm = res$algo_name,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, dfs)
}

make_allocation_df <- function(results, T, K, sample_every = 100) {
  # We need arm_means_history — rebuild from stored data
  # Actually, we stored arm_means_history per sim but not across sims
  # Instead, track arm pulls over time from actions
  # For this we need to re-extract from results — let me store actions
  # Actually, we can use arm_means_history averaged. But we didn't save it across sims.
  # Let's just use the regret curves for visualization and skip allocation trajectory
  # in the main runner. We'll add a special run for allocation tracking.
  NULL
}


# =============================================================================
# SECTION 4: SCENARIO 1 — STANDARD BERNOULLI (K=5)
# =============================================================================

cat("\n============================================================\n")
cat("SCENARIO 1: Standard Bernoulli (K=5 arms)\n")
cat("True rates: (0.1, 0.3, 0.5, 0.7, 0.9)\n")
cat("============================================================\n\n")

arms1 <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T1 <- 10000
n_sims1 <- 200

cat("Running Epsilon-Greedy (eps=0.05)...\n")
res1_eg05 <- run_simulation(epsilon_greedy, arms1, T1, n_sims1, "EG-0.05", 0.9, epsilon = 0.05)

cat("Running Epsilon-Greedy (eps=0.10)...\n")
res1_eg10 <- run_simulation(epsilon_greedy, arms1, T1, n_sims1, "EG-0.10", 0.9, epsilon = 0.10)

cat("Running Epsilon-Greedy (eps=0.20)...\n")
res1_eg20 <- run_simulation(epsilon_greedy, arms1, T1, n_sims1, "EG-0.20", 0.9, epsilon = 0.20)

cat("Running UCB1...\n")
res1_ucb <- run_simulation(ucb1, arms1, T1, n_sims1, "UCB1", 0.9)

cat("Running UCB-Tuned...\n")
res1_ucbt <- run_simulation(ucb_tuned, arms1, T1, n_sims1, "UCB-Tuned", 0.9)

cat("Running Thompson Sampling...\n")
res1_ts <- run_simulation(thompson_sampling, arms1, T1, n_sims1, "TS", 0.9)

cat("Running EXP3...\n")
res1_exp3 <- run_simulation(exp3, arms1, T1, n_sims1, "EXP3", 0.9, gamma = 0.1)

scenario1_results <- list(res1_eg05, res1_eg10, res1_eg20, res1_ucb, res1_ucbt, res1_ts, res1_exp3)
scenario1_summary <- summary_stats(scenario1_results)
cat("\n--- Scenario 1 Summary ---\n")
print(scenario1_summary)


# =============================================================================
# SECTION 5: SCENARIO 2 — HARD DISCRIMINATION (K=5, close rates)
# =============================================================================

cat("\n============================================================\n")
cat("SCENARIO 2: Hard Discrimination (K=5, close rates)\n")
cat("True rates: (0.40, 0.42, 0.45, 0.48, 0.50)\n")
cat("============================================================\n\n")

arms2 <- c(0.40, 0.42, 0.45, 0.48, 0.50)
T2 <- 10000
n_sims2 <- 200

cat("Running Epsilon-Greedy (eps=0.10)...\n")
res2_eg <- run_simulation(epsilon_greedy, arms2, T2, n_sims2, "EG-0.10", 0.50, epsilon = 0.10)

cat("Running UCB1...\n")
res2_ucb <- run_simulation(ucb1, arms2, T2, n_sims2, "UCB1", 0.50)

cat("Running UCB-Tuned...\n")
res2_ucbt <- run_simulation(ucb_tuned, arms2, T2, n_sims2, "UCB-Tuned", 0.50)

cat("Running Thompson Sampling...\n")
res2_ts <- run_simulation(thompson_sampling, arms2, T2, n_sims2, "TS", 0.50)

cat("Running EXP3...\n")
res2_exp3 <- run_simulation(exp3, arms2, T2, n_sims2, "EXP3", 0.50, gamma = 0.1)

scenario2_results <- list(res2_eg, res2_ucb, res2_ucbt, res2_ts, res2_exp3)
scenario2_summary <- summary_stats(scenario2_results)
cat("\n--- Scenario 2 Summary ---\n")
print(scenario2_summary)


# =============================================================================
# SECTION 6: SCENARIO 3 — MANY ARMS (K=20)
# =============================================================================

cat("\n============================================================\n")
cat("SCENARIO 3: Many Arms (K=20)\n")
cat("Rates drawn from Beta(2,5), best arm identified per simulation\n")
cat("============================================================\n\n")

T3 <- 20000
n_sims3 <- 200

# For many arms, we need a modified runner that generates arms per sim
run_many_arms_sim <- function(algo_fn, T, n_sims, algo_name, ...) {
  cum_regret_all <- matrix(0, nrow = T, ncol = n_sims)
  best_arm_pct_all <- numeric(n_sims)
  switching_all <- numeric(n_sims)

  for (i in seq_len(n_sims)) {
    if (i %% 50 == 0) cat(sprintf("  %s sim %d/%d\n", algo_name, i, n_sims))
    set.seed(42 + i)  # reproducible per sim
    arms_k <- rbeta(20, 2, 5)
    best_rate <- max(arms_k)
    result <- algo_fn(arms_k, T, ...)

    cum_optimal <- cumsum(rep(best_rate, T))
    cum_actual <- cumsum(result$rewards)
    cum_regret_all[, i] <- cum_optimal - cum_actual

    best_arm <- which.max(arms_k)
    best_arm_pct_all[i] <- mean(result$actions == best_arm) * 100

    switches <- sum(diff(result$actions) != 0)
    switching_all[i] <- switches / (T - 1) * 100
  }

  list(
    algo_name = algo_name,
    cum_regret = cum_regret_all,
    simple_regret = numeric(0),  # not meaningful here
    best_arm_pct = best_arm_pct_all,
    switching_rate = switching_all
  )
}

cat("Running UCB1...\n")
res3_ucb <- run_many_arms_sim(ucb1, T3, n_sims3, "UCB1")

cat("Running UCB-Tuned...\n")
res3_ucbt <- run_many_arms_sim(ucb_tuned, T3, n_sims3, "UCB-Tuned")

cat("Running Thompson Sampling...\n")
res3_ts <- run_many_arms_sim(thompson_sampling, T3, n_sims3, "TS")

cat("Running Epsilon-Greedy (eps=0.10)...\n")
res3_eg <- run_many_arms_sim(epsilon_greedy, T3, n_sims3, "EG-0.10", epsilon = 0.10)

cat("Running EXP3...\n")
res3_exp3 <- run_many_arms_sim(exp3, T3, n_sims3, "EXP3", gamma = 0.1)

scenario3_results <- list(res3_eg, res3_ucb, res3_ucbt, res3_ts, res3_exp3)
cat("\n--- Scenario 3 Summary ---\n")
cat(sprintf("  UCB1:         Mean final regret = %.1f, Best arm %% = %.1f%%\n",
            mean(res3_ucb$cum_regret[T3,]), mean(res3_ucb$best_arm_pct)))
cat(sprintf("  UCB-Tuned:    Mean final regret = %.1f, Best arm %% = %.1f%%\n",
            mean(res3_ucbt$cum_regret[T3,]), mean(res3_ucbt$best_arm_pct)))
cat(sprintf("  TS:           Mean final regret = %.1f, Best arm %% = %.1f%%\n",
            mean(res3_ts$cum_regret[T3,]), mean(res3_ts$best_arm_pct)))
cat(sprintf("  EG-0.10:      Mean final regret = %.1f, Best arm %% = %.1f%%\n",
            mean(res3_eg$cum_regret[T3,]), mean(res3_eg$best_arm_pct)))
cat(sprintf("  EXP3:         Mean final regret = %.1f, Best arm %% = %.1f%%\n",
            mean(res3_exp3$cum_regret[T3,]), mean(res3_exp3$best_arm_pct)))


# =============================================================================
# SECTION 7: SCENARIO 4 — NON-STATIONARY (K=5, drifting rates)
# =============================================================================

cat("\n============================================================\n")
cat("SCENARIO 4: Non-Stationary (K=5, sudden shifts)\n")
cat("Rate changes at t=3000 and t=7000\n")
cat("============================================================\n\n")

T4 <- 10000
n_sims4 <- 200

# Define rate function with shifts
# Phase 1 (t<=3000):    (0.1, 0.3, 0.5, 0.7, 0.9)  — arm 5 best
# Phase 2 (3001-7000):  (0.9, 0.7, 0.5, 0.3, 0.1)  — arm 1 best
# Phase 3 (t>7000):     (0.3, 0.9, 0.5, 0.7, 0.1)  — arm 2 best
rate_phases <- function(t) {
  if (t <= 3000) {
    c(0.1, 0.3, 0.5, 0.7, 0.9)
  } else if (t <= 7000) {
    c(0.9, 0.7, 0.5, 0.3, 0.1)
  } else {
    c(0.3, 0.9, 0.5, 0.7, 0.1)
  }
}

# Vectorized version for plotting
rate_at_time <- Vectorize(function(t) {
  if (t <= 3000) return(5)
  if (t <= 7000) return(1)
  return(2)
})

cat("Running UCB1 (non-stationary)...\n")
res4_ucb <- run_simulation_nonstat(ucb1, T4, n_sims4, "UCB1", rate_phases)

cat("Running Thompson Sampling (non-stationary)...\n")
res4_ts <- run_simulation_nonstat(thompson_sampling, T4, n_sims4, "TS", rate_phases)

cat("Running Sliding-Window TS (non-stationary)...\n")
res4_swts <- run_simulation_nonstat(sw_thompson_sampling, T4, n_sims4, "SW-TS", rate_phases, window = 500)

cat("Running Discounted UCB (non-stationary)...\n")
res4_ducb <- run_simulation_nonstat(discounted_ucb, T4, n_sims4, "D-UCB", rate_phases, discount = 0.98)

cat("Running EG-0.10 (non-stationary)...\n")
res4_eg <- run_simulation_nonstat(epsilon_greedy, T4, n_sims4, "EG-0.10", rate_phases, epsilon = 0.10)

cat("\n--- Scenario 4 Summary ---\n")
cat(sprintf("  UCB1:     Mean final regret = %.1f\n", mean(res4_ucb$cum_regret[T4,])))
cat(sprintf("  TS:       Mean final regret = %.1f\n", mean(res4_ts$cum_regret[T4,])))
cat(sprintf("  SW-TS:    Mean final regret = %.1f\n", mean(res4_swts$cum_regret[T4,])))
cat(sprintf("  D-UCB:    Mean final regret = %.1f\n", mean(res4_ducb$cum_regret[T4,])))
cat(sprintf("  EG-0.10:  Mean final regret = %.1f\n", mean(res4_eg$cum_regret[T4,])))


# =============================================================================
# SECTION 8: SCENARIO 5 — NORMAL REWARDS
# =============================================================================

cat("\n============================================================\n")
cat("SCENARIO 5: Normal Rewards (K=5)\n")
cat("Means: (0.0, 0.5, 1.0, 1.5, 2.0), SDs: (1.0, 0.8, 1.2, 0.5, 1.5)\n")
cat("============================================================\n\n")

arms5_means <- c(0.0, 0.5, 1.0, 1.5, 2.0)
arms5_sds <- c(1.0, 0.8, 1.2, 0.5, 1.5)
T5 <- 10000
n_sims5 <- 200

cat("Running UCB1 (Normal)...\n")
res5_ucb <- run_simulation_normal(ucb1_normal, arms5_means, arms5_sds, T5, n_sims5, "UCB1")

cat("Running Thompson Sampling Normal...\n")
res5_ts <- run_simulation_normal(thompson_sampling_normal, arms5_means, arms5_sds, T5, n_sims5, "TS-Normal")

scenario5_results <- list(res5_ucb, res5_ts)
scenario5_summary <- summary_stats(scenario5_results)
cat("\n--- Scenario 5 Summary ---\n")
print(scenario5_summary)


# =============================================================================
# SECTION 9: VISUALIZATION
# =============================================================================

cat("\n============================================================\n")
cat("Generating figures...\n")
cat("============================================================\n\n")

# Set ggplot2 theme
theme_set(theme_minimal(base_size = 12))

# --- Figure 1: Regret curves for Scenario 1 ---
df1 <- make_regret_df(scenario1_results, T1)
p1 <- ggplot(df1, aes(x = t, y = mean_regret, color = factor(Algorithm))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(Algorithm)),
              alpha = 0.15, color = NA) +
  labs(title = "Scenario 1: Standard Bernoulli — Cumulative Regret",
       subtitle = "K=5 arms, rates=(0.1, 0.3, 0.5, 0.7, 0.9), T=10000, 200 sims",
       x = "Time step", y = "Cumulative Regret",
       color = "Algorithm", fill = "Algorithm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

ggsave("fig1_regret_standard.png",
       p1, width = 10, height = 6, dpi = 300)
cat("Saved fig1_regret_standard.png\n")

# --- Figure 2: Regret curves for Scenario 2 ---
df2 <- make_regret_df(scenario2_results, T2)
p2 <- ggplot(df2, aes(x = t, y = mean_regret, color = factor(Algorithm))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(Algorithm)),
              alpha = 0.15, color = NA) +
  labs(title = "Scenario 2: Hard Discrimination — Cumulative Regret",
       subtitle = "K=5 arms, rates=(0.40, 0.42, 0.45, 0.48, 0.50), T=10000, 200 sims",
       x = "Time step", y = "Cumulative Regret",
       color = "Algorithm", fill = "Algorithm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

ggsave("fig2_regret_hard.png",
       p2, width = 10, height = 6, dpi = 300)
cat("Saved fig2_regret_hard.png\n")

# --- Figure 3: Regret for Many Arms (K=20) ---
df3 <- make_regret_df(scenario3_results, T3)
p3 <- ggplot(df3, aes(x = t, y = mean_regret, color = factor(Algorithm))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(Algorithm)),
              alpha = 0.15, color = NA) +
  labs(title = "Scenario 3: Many Arms (K=20) — Cumulative Regret",
       subtitle = "Rates ~ Beta(2,5), T=20000, 200 sims",
       x = "Time step", y = "Cumulative Regret",
       color = "Algorithm", fill = "Algorithm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

ggsave("fig3_many_arms.png",
       p3, width = 10, height = 6, dpi = 300)
cat("Saved fig3_many_arms.png\n")

# --- Figure 4: Non-stationary adaptation ---
ns_results <- list(res4_ucb, res4_ts, res4_swts, res4_ducb, res4_eg)
df4 <- do.call(rbind, lapply(ns_results, function(res) {
  T_ns <- nrow(res$cum_regret)
  means <- rowMeans(res$cum_regret)
  sds <- apply(res$cum_regret, 1, sd)
  n <- ncol(res$cum_regret)
  se <- sds / sqrt(n)
  data.frame(
    t = seq_len(T_ns),
    mean_regret = means,
    lower = means - 1.96 * se,
    upper = means + 1.96 * se,
    Algorithm = res$algo_name,
    stringsAsFactors = FALSE
  )
}))

p4 <- ggplot(df4, aes(x = t, y = mean_regret, color = factor(Algorithm))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(Algorithm)),
              alpha = 0.12, color = NA) +
  geom_vline(xintercept = c(3000, 7000), linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 3000, y = Inf, label = "Shift 1", vjust = 2, hjust = -0.2, color = "red", size = 3) +
  annotate("text", x = 7000, y = Inf, label = "Shift 2", vjust = 2, hjust = -0.2, color = "red", size = 3) +
  labs(title = "Scenario 4: Non-Stationary — Cumulative Regret Under Distribution Shifts",
       subtitle = "K=5, sudden shifts at t=3000 and t=7000, T=10000, 200 sims",
       x = "Time step", y = "Cumulative Regret",
       color = "Algorithm", fill = "Algorithm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

ggsave("fig4_nonstationary.png",
       p4, width = 10, height = 6, dpi = 300)
cat("Saved fig4_nonstationary.png\n")

# --- Figure 5: Arm Allocation Trajectories ---
# Run a single simulation per algorithm to get arm pull trajectories
cat("Running single-sim arm allocation tracking for Figure 5...\n")
set.seed(42)
alloc_data <- list()
algo_labels <- c("UCB1", "TS", "EG-0.10", "UCB-Tuned", "EXP3")
algo_fns <- list(
  function(a, T) ucb1(a, T),
  function(a, T) thompson_sampling(a, T),
  function(a, T) epsilon_greedy(a, T, epsilon = 0.10),
  function(a, T) ucb_tuned(a, T),
  function(a, T) exp3(a, T, gamma = 0.1)
)

T_alloc <- T1
sample_every <- 50
time_points <- seq(1, T_alloc, by = sample_every)
n_tp <- length(time_points)

for (j in seq_along(algo_labels)) {
  result <- algo_fns[[j]](arms1, T_alloc)
  # Compute cumulative proportion per arm at each time point
  alloc_matrix <- matrix(0, nrow = n_tp, ncol = 5)
  for (a in 1:5) {
    cum_pulls <- cumsum(result$actions == a)
    alloc_matrix[, a] <- cum_pulls[time_points] / time_points
  }
  for (a in 1:5) {
    alloc_data[[length(alloc_data) + 1]] <- data.frame(
      t = time_points,
      proportion = alloc_matrix[, a],
      Arm = paste0("Arm ", a, " (p=", arms1[a], ")"),
      Algorithm = algo_labels[j],
      stringsAsFactors = FALSE
    )
  }
}
df5 <- do.call(rbind, alloc_data)

p5 <- ggplot(df5, aes(x = t, y = proportion, color = factor(Arm))) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ factor(Algorithm), ncol = 3) +
  labs(title = "Arm Allocation Trajectories Over Time",
       subtitle = "Scenario 1: K=5, rates=(0.1, 0.3, 0.5, 0.7, 0.9)",
       x = "Time step", y = "Cumulative pull proportion",
       color = "Arm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d()

ggsave("fig5_arm_allocation.png",
       p5, width = 12, height = 8, dpi = 300)
cat("Saved fig5_arm_allocation.png\n")

# --- Figure 6: Performance Matrix Heatmap ---
# Build a matrix: rows = algorithms, columns = scenarios
# Metrics: normalized final regret (lower = better)
perf_algos <- c("EG-0.10", "UCB1", "UCB-Tuned", "TS", "EXP3")
perf_scenarios <- c("S1: Standard", "S2: Hard", "S3: K=20", "S4: NonStat", "S5: Normal")

# Extract final regrets
s1_vals <- c(
  mean(res1_eg10$cum_regret[T1,]),
  mean(res1_ucb$cum_regret[T1,]),
  mean(res1_ucbt$cum_regret[T1,]),
  mean(res1_ts$cum_regret[T1,]),
  mean(res1_exp3$cum_regret[T1,])
)
s2_vals <- c(
  mean(res2_eg$cum_regret[T2,]),
  mean(res2_ucb$cum_regret[T2,]),
  mean(res2_ucbt$cum_regret[T2,]),
  mean(res2_ts$cum_regret[T2,]),
  mean(res2_exp3$cum_regret[T2,])
)
s3_vals <- c(
  mean(res3_eg$cum_regret[T3,]),
  mean(res3_ucb$cum_regret[T3,]),
  mean(res3_ucbt$cum_regret[T3,]),
  mean(res3_ts$cum_regret[T3,]),
  mean(res3_exp3$cum_regret[T3,])
)
s4_vals <- c(
  mean(res4_eg$cum_regret[T4,]),
  mean(res4_ucb$cum_regret[T4,]),
  NA,  # no ucb-tuned for nonstat
  mean(res4_ts$cum_regret[T4,]),
  NA   # no exp3 for nonstat in this form
)
s5_vals <- c(
  NA,  # no EG for normal
  mean(res5_ucb$cum_regret[T5,]),
  NA,
  mean(res5_ts$cum_regret[T5,]),
  NA
)

# For non-stationary, fill in D-UCB and SW-TS
# Use the algorithms we actually ran
s4_vals_full <- c(
  mean(res4_eg$cum_regret[T4,]),
  mean(res4_ucb$cum_regret[T4,]),
  mean(res4_ducb$cum_regret[T4,]),
  mean(res4_ts$cum_regret[T4,]),
  mean(res4_swts$cum_regret[T4,])
)
# For S5, we only have UCB1 and TS-Normal
s5_vals_full <- c(
  NA,
  mean(res5_ucb$cum_regret[T5,]),
  NA,
  mean(res5_ts$cum_regret[T5,]),
  NA
)

# Build heatmap data: use the 5 algorithms for S1-S3, special for S4/S5
heatmap_df <- expand.grid(
  Algorithm = perf_algos,
  Scenario = perf_scenarios,
  stringsAsFactors = FALSE
)

# Fill values
heatmap_df$value <- NA
# S1
for (i in 1:5) heatmap_df$value[heatmap_df$Algorithm == perf_algos[i] & heatmap_df$Scenario == perf_scenarios[1]] <- s1_vals[i]
# S2
for (i in 1:5) heatmap_df$value[heatmap_df$Algorithm == perf_algos[i] & heatmap_df$Scenario == perf_scenarios[2]] <- s2_vals[i]
# S3
for (i in 1:5) heatmap_df$value[heatmap_df$Algorithm == perf_algos[i] & heatmap_df$Scenario == perf_scenarios[3]] <- s3_vals[i]
# S4
for (i in 1:5) heatmap_df$value[heatmap_df$Algorithm == perf_algos[i] & heatmap_df$Scenario == perf_scenarios[4]] <- s4_vals_full[i]
# S5
for (i in 1:5) heatmap_df$value[heatmap_df$Algorithm == perf_algos[i] & heatmap_df$Scenario == perf_scenarios[5]] <- s5_vals_full[i]

# Normalize within each scenario (0-1 scale)
heatmap_df$norm_value <- NA
for (sc in perf_scenarios) {
  idx <- heatmap_df$Scenario == sc
  vals <- heatmap_df$value[idx]
  if (any(!is.na(vals))) {
    rng <- range(vals, na.rm = TRUE)
    if (rng[2] > rng[1]) {
      heatmap_df$norm_value[idx] <- (vals - rng[1]) / (rng[2] - rng[1])
    } else {
      heatmap_df$norm_value[idx] <- 0.5
    }
  }
}

# Add label with actual value
heatmap_df$label <- ifelse(is.na(heatmap_df$value), "N/A",
                           sprintf("%.0f", heatmap_df$value))

p6 <- ggplot(heatmap_df, aes(x = factor(Scenario, levels = perf_scenarios),
                              y = factor(Algorithm, levels = rev(perf_algos)))) +
  geom_tile(aes(fill = norm_value), color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#2ecc71", high = "#e74c3c", na.value = "grey90",
                      name = "Normalized\nRegret\n(low=good)") +
  labs(title = "Algorithm x Scenario Performance Matrix",
       subtitle = "Cell values = mean final cumulative regret; color = normalized (green=low, red=high)",
       x = "Scenario", y = "Algorithm") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave("fig6_performance_matrix.png",
       p6, width = 10, height = 6, dpi = 300)
cat("Saved fig6_performance_matrix.png\n")

# --- Figure 7: Epsilon Sensitivity ---
eps_values <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50)
eps_regrets <- numeric(length(eps_values))
eps_bestarm <- numeric(length(eps_values))

cat("Running epsilon sensitivity analysis...\n")
for (j in seq_along(eps_values)) {
  res_tmp <- run_simulation(epsilon_greedy, arms1, T1, 100, paste0("EG-", eps_values[j]), 0.9, epsilon = eps_values[j])
  eps_regrets[j] <- mean(res_tmp$cum_regret[T1, ])
  eps_bestarm[j] <- mean(res_tmp$best_arm_pct)
}

eps_df <- data.frame(
  epsilon = eps_values,
  regret = eps_regrets,
  best_arm_pct = eps_bestarm
)

p7a <- ggplot(eps_df, aes(x = epsilon, y = regret)) +
  geom_line(linewidth = 1, color = "#3498db") +
  geom_point(size = 3, color = "#3498db") +
  scale_x_continuous(breaks = eps_values) +
  labs(title = "Epsilon Sensitivity: Cumulative Regret",
       x = expression(epsilon), y = "Mean Final Cumulative Regret")

p7b <- ggplot(eps_df, aes(x = epsilon, y = best_arm_pct)) +
  geom_line(linewidth = 1, color = "#e74c3c") +
  geom_point(size = 3, color = "#e74c3c") +
  scale_x_continuous(breaks = eps_values) +
  labs(title = "Epsilon Sensitivity: Best Arm Identification",
       x = expression(epsilon), y = "% Pulls of Best Arm")

# Combine using a grid approach
p7_combined <- arrangeGrob(p7a, p7b, ncol = 2,
                           top = textGrob("Epsilon-Greedy Sensitivity Analysis (Scenario 1)",
                                          gp = gpar(fontsize = 14, fontface = "bold")))

ggsave("fig7_epsilon_sensitivity.png",
       p7_combined, width = 12, height = 5, dpi = 300)
cat("Saved fig7_epsilon_sensitivity.png\n")

# --- Figure 8: Normal Rewards ---
df8 <- make_regret_df(scenario5_results, T5)
p8 <- ggplot(df8, aes(x = t, y = mean_regret, color = factor(Algorithm))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(Algorithm)),
              alpha = 0.15, color = NA) +
  labs(title = "Scenario 5: Normal Rewards — Cumulative Regret",
       subtitle = "K=5 arms, means=(0.0, 0.5, 1.0, 1.5, 2.0), heterogeneous variances, T=10000, 200 sims",
       x = "Time step", y = "Cumulative Regret",
       color = "Algorithm", fill = "Algorithm") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

ggsave("fig8_normal_rewards.png",
       p8, width = 10, height = 6, dpi = 300)
cat("Saved fig8_normal_rewards.png\n")


# =============================================================================
# SECTION 10: ADDITIONAL VISUALIZATIONS
# =============================================================================

# --- Figure 9: Bar chart of best arm identification across scenarios ---
bar_df <- data.frame(
  Algorithm = c(
    rep(c("EG-0.10", "UCB1", "UCB-Tuned", "TS", "EXP3"), 2),
    rep(c("EG-0.10", "UCB1", "UCB-Tuned", "TS", "EXP3"), 1)
  ),
  Scenario = c(
    rep("S1: Standard", 5),
    rep("S2: Hard", 5),
    rep("S3: K=20", 5)
  ),
  BestArmPct = c(
    sapply(scenario1_results[c(2,4,5,6,7)], function(x) mean(x$best_arm_pct)),  # eg10, ucb, ucbt, ts, exp3
    sapply(scenario2_results, function(x) mean(x$best_arm_pct)),
    sapply(scenario3_results, function(x) mean(x$best_arm_pct))
  ),
  stringsAsFactors = FALSE
)

p9 <- ggplot(bar_df, aes(x = factor(Algorithm, levels = c("TS", "UCB-Tuned", "UCB1", "EG-0.10", "EXP3")),
                          y = BestArmPct, fill = factor(Scenario))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Best Arm Identification Rate Across Scenarios",
       x = "Algorithm", y = "% Pulls of Best Arm",
       fill = "Scenario") +
  theme(legend.position = "bottom") +
  scale_fill_viridis_d()

ggsave("fig9_bestarm_barchart.png",
       p9, width = 10, height = 6, dpi = 300)
cat("Saved fig9_bestarm_barchart.png\n")


# --- Figure 10: Switching rate comparison ---
switch_df <- data.frame(
  Algorithm = factor(c("EG-0.05", "EG-0.10", "EG-0.20", "UCB1", "UCB-Tuned", "TS", "EXP3")),
  SwitchingRate = sapply(scenario1_results, function(x) mean(x$switching_rate))
)

p10 <- ggplot(switch_df, aes(x = Algorithm, y = SwitchingRate, fill = Algorithm)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Switching Rate by Algorithm (Scenario 1)",
       subtitle = "Lower switching = more stable exploitation",
       x = "Algorithm", y = "% Steps with Arm Switch") +
  theme(legend.position = "none") +
  scale_fill_viridis_d()

ggsave("fig10_switching_rate.png",
       p10, width = 10, height = 5, dpi = 300)
cat("Saved fig10_switching_rate.png\n")


# =============================================================================
# SECTION 11: FINAL SUMMARY TABLES
# =============================================================================

cat("\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n\n")

cat("--- Scenario 1: Standard Bernoulli ---\n")
print(format(scenario1_summary, digits = 3))

cat("\n--- Scenario 2: Hard Discrimination ---\n")
print(format(scenario2_summary, digits = 3))

cat("\n--- Scenario 3: Many Arms (K=20) ---\n")
cat(sprintf("%-12s %12s %12s\n", "Algorithm", "Final Regret", "Best Arm %"))
cat(sprintf("%-12s %12.1f %12.1f\n", "EG-0.10", mean(res3_eg$cum_regret[T3,]), mean(res3_eg$best_arm_pct)))
cat(sprintf("%-12s %12.1f %12.1f\n", "UCB1", mean(res3_ucb$cum_regret[T3,]), mean(res3_ucb$best_arm_pct)))
cat(sprintf("%-12s %12.1f %12.1f\n", "UCB-Tuned", mean(res3_ucbt$cum_regret[T3,]), mean(res3_ucbt$best_arm_pct)))
cat(sprintf("%-12s %12.1f %12.1f\n", "TS", mean(res3_ts$cum_regret[T3,]), mean(res3_ts$best_arm_pct)))
cat(sprintf("%-12s %12.1f %12.1f\n", "EXP3", mean(res3_exp3$cum_regret[T3,]), mean(res3_exp3$best_arm_pct)))

cat("\n--- Scenario 4: Non-Stationary ---\n")
cat(sprintf("%-12s %12s %12s\n", "Algorithm", "Final Regret", "Switching %"))
cat(sprintf("%-12s %12.1f %12.1f\n", "UCB1", mean(res4_ucb$cum_regret[T4,]), mean(res4_ucb$switching_rate)))
cat(sprintf("%-12s %12.1f %12.1f\n", "TS", mean(res4_ts$cum_regret[T4,]), mean(res4_ts$switching_rate)))
cat(sprintf("%-12s %12.1f %12.1f\n", "SW-TS", mean(res4_swts$cum_regret[T4,]), mean(res4_swts$switching_rate)))
cat(sprintf("%-12s %12.1f %12.1f\n", "D-UCB", mean(res4_ducb$cum_regret[T4,]), mean(res4_ducb$switching_rate)))
cat(sprintf("%-12s %12.1f %12.1f\n", "EG-0.10", mean(res4_eg$cum_regret[T4,]), mean(res4_eg$switching_rate)))

cat("\n--- Scenario 5: Normal Rewards ---\n")
print(format(scenario5_summary, digits = 3))

cat("\n============================================================\n")
cat("ALL FIGURES SAVED SUCCESSFULLY\n")
cat("Total figures: 10\n")
cat("============================================================\n")
