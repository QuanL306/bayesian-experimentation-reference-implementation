# Week 20 Capstone: Comprehensive Experiment Design
# Three scenarios: A/B test, multi-variant, switchback
# Power analysis, multiple comparisons, and design trade-offs

library(ggplot2)

set.seed(2024)

cat("=============================================================\n")
cat("=== Week 20 Capstone: Experiment Design Portfolio ===\n")
cat("=============================================================\n\n")

# Output directory for figures
fig_dir <- ""

# ============================================================================
# SCENARIO 1: E-commerce Checkout A/B Test
# ============================================================================

cat("=== SCENARIO 1: E-commerce Checkout A/B Test ===\n\n")

# --- Parameters ---
p_A <- 0.08       # Control conversion rate
mde <- 0.01       # Minimum Detectable Effect
p_B <- p_A + mde  # Treatment conversion rate (under alternative)
alpha <- 0.05     # Type I error rate
target_power <- 0.80

# --- 1a. Frequentist Analytical Power ---
cat("--- 1a. Frequentist Analytical Power ---\n")

# Sample size formula for two-proportion z-test
# n = [(z_{1-alpha/2} + z_{1-beta})^2 * (p_A*(1-p_A) + p_B*(1-p_B))] / (p_B - p_A)^2
z_alpha <- qnorm(1 - alpha/2)
pooled_var <- p_A * (1 - p_A) + p_B * (1 - p_B)

# Analytical power for given n
freq_power <- function(n, p_A, p_B, alpha) {
  se <- sqrt(p_A*(1-p_A)/n + p_B*(1-p_B)/n)
  delta <- p_B - p_A
  z_crit <- qnorm(1 - alpha/2)
  # Power = P(reject H0 | H1 true)
  power <- pnorm(delta/se - z_crit) + pnorm(-delta/se - z_crit)
  return(power)
}

# Required sample size for target power
z_beta <- qnorm(target_power)
n_required <- ceiling((z_alpha + z_beta)^2 * pooled_var / mde^2)
cat(sprintf("Required n per group (analytical, power=%.0f%%): %d\n",
            target_power*100, n_required))
cat(sprintf("Total sample size: %d\n\n", n_required * 2))

# Power at various sample sizes
sample_sizes <- c(500, 1000, 2000, 5000, 10000)
freq_powers <- sapply(sample_sizes, function(n) freq_power(n, p_A, p_B, alpha))

cat("Frequentist power by sample size:\n")
for (i in seq_along(sample_sizes)) {
  cat(sprintf("  n = %5d per group: power = %.3f\n", sample_sizes[i], freq_powers[i]))
}
cat("\n")

# --- 1b. Bayesian Power Analysis (ROPE-based simulation) ---
cat("--- 1b. Bayesian Power Analysis (ROPE simulation, 500 experiments each) ---\n")

n_sims <- 500
rope_bounds <- c(-0.005, 0.005)  # ROPE: effects within +/- 0.5pp considered negligible

bayesian_power_sim <- function(n, p_A, p_B, n_sims, rope_bounds, n_posterior = 5000) {
  detected <- 0
  for (i in 1:n_sims) {
    # Simulate data
    y_A <- rbinom(1, n, p_A)
    y_B <- rbinom(1, n, p_B)

    # Posterior via Beta-Binomial conjugacy
    # Prior: Beta(1, 1) — uniform
    post_A <- rbeta(n_posterior, y_A + 1, n - y_A + 1)
    post_B <- rbeta(n_posterior, y_B + 1, n - y_B + 1)
    delta_post <- post_B - post_A

    # ROPE decision
    prop_in_rope <- mean(delta_post >= rope_bounds[1] & delta_post <= rope_bounds[2])
    prop_above_rope <- mean(delta_post > rope_bounds[2])

    # "Detected" if > 95% posterior mass outside ROPE and above
    if (prop_above_rope > 0.95) {
      detected <- detected + 1
    }
  }
  return(detected / n_sims)
}

bayes_powers <- numeric(length(sample_sizes))
for (i in seq_along(sample_sizes)) {
  bayes_powers[i] <- bayesian_power_sim(sample_sizes[i], p_A, p_B, n_sims, rope_bounds)
  cat(sprintf("  n = %5d per group: Bayesian power = %.3f (of %d sims)\n",
              sample_sizes[i], bayes_powers[i], n_sims))
}
cat("\n")

# --- 1c. Comparison Table ---
cat("--- 1c. Scenario 1: Power Comparison ---\n")
power_df_s1 <- data.frame(
  n = sample_sizes,
  Freq_Power = freq_powers,
  Bayes_Power = bayes_powers
)
cat(sprintf("  %-10s  %-15s  %-15s\n", "n/group", "Freq Power", "Bayes Power"))
for (i in 1:nrow(power_df_s1)) {
  cat(sprintf("  %-10d  %-15.3f  %-15.3f\n",
              power_df_s1$n[i], power_df_s1$Freq_Power[i], power_df_s1$Bayes_Power[i]))
}
cat("\n")


# ============================================================================
# SCENARIO 2: SaaS Pricing Page (Multi-Variant Test)
# ============================================================================

cat("=== SCENARIO 2: SaaS Pricing Page (Multi-Variant) ===\n\n")

# --- Parameters ---
k_variants <- 4                    # control + 3 treatments
p_control <- 0.15                  # control trial-to-paid conversion
mde_s2 <- 0.03                     # MDE for the best variant
# True conversion rates for each variant
p_variants <- c(p_control, p_control + 0.01, p_control + 0.02, p_control + mde_s2)
names(p_variants) <- c("Control", "T1", "T2", "T3")

alpha_s2 <- 0.05
n_comparisons <- k_variants - 1  # 3 comparisons vs control

cat("Variant true conversion rates:\n")
for (i in seq_along(p_variants)) {
  cat(sprintf("  %s: %.3f\n", names(p_variants)[i], p_variants[i]))
}
cat(sprintf("\nNumber of comparisons vs control: %d\n", n_comparisons))

# --- 2a. Multiple Comparison Corrections ---
cat("\n--- 2a. Multiple Comparison Corrections ---\n")

# Bonferroni-corrected alpha
alpha_bonf <- alpha_s2 / n_comparisons
cat(sprintf("Bonferroni-corrected alpha: %.4f (= %.2f / %d)\n",
            alpha_bonf, alpha_s2, n_comparisons))

# Holm-Bonferroni: step-down procedure (less conservative)
# For k=3 comparisons, smallest p gets alpha/3, next alpha/2, largest alpha/1
cat("Holm-Bonferroni: step-down procedure (alpha/(k+1-i) for i-th smallest p)\n")

# Sample size per variant (Bonferroni-corrected)
z_alpha_bonf <- qnorm(1 - alpha_bonf/2)
z_beta_s2 <- qnorm(target_power)
p_B_s2 <- p_control + mde_s2
pooled_var_s2 <- p_control*(1-p_control) + p_B_s2*(1-p_B_s2)
n_required_s2_bonf <- ceiling((z_alpha_bonf + z_beta_s2)^2 * pooled_var_s2 / mde_s2^2)
cat(sprintf("\nRequired n per variant (Bonferroni, 80%% power): %d\n", n_required_s2_bonf))

# Sample size without correction (for comparison)
z_alpha_unc <- qnorm(1 - alpha_s2/2)
n_required_s2_unc <- ceiling((z_alpha_unc + z_beta_s2)^2 * pooled_var_s2 / mde_s2^2)
cat(sprintf("Required n per variant (uncorrected, 80%% power): %d\n", n_required_s2_unc))
cat(sprintf("Inflation factor due to Bonferroni: %.2fx\n\n",
            n_required_s2_bonf / n_required_s2_unc))

# --- 2b. Simulate Family-Wise Error Rate ---
cat("--- 2b. Simulating Family-Wise Error Rate (1000 experiments) ---\n")
cat("    (All variants actually have same rate p=0.15 — null is true)\n\n")

n_sims_s2 <- 1000
n_per_variant <- 2000  # Use a moderate sample size

simulate_multi_variant <- function(n, p_true, k, alpha, method = "none") {
  # Simulate k variants all from same p_true (null hypothesis)
  y <- rbinom(k, n, p_true)
  p_hat <- y / n

  # P-values for each treatment vs control (two-proportion z-test)
  p_control_hat <- p_hat[1]
  p_treat_hats <- p_hat[-1]

  p_values <- numeric(k - 1)
  for (j in 1:(k-1)) {
    p_pool <- (y[1] + y[j+1]) / (2 * n)
    se <- sqrt(p_pool * (1 - p_pool) * (2/n))
    z_stat <- abs(p_treat_hats[j] - p_control_hat) / se
    p_values[j] <- 2 * pnorm(-z_stat)
  }

  # Apply correction
  any_significant <- FALSE
  if (method == "none") {
    any_significant <- any(p_values < alpha)
  } else if (method == "bonferroni") {
    any_significant <- any(p_values < alpha / (k - 1))
  } else if (method == "holm") {
    sorted_p <- sort(p_values)
    m <- length(sorted_p)
    for (i in 1:m) {
      threshold <- alpha / (m + 1 - i)
      if (sorted_p[i] < threshold) {
        any_significant <- TRUE
      } else {
        break
      }
    }
  }
  return(any_significant)
}

# Run simulations under null for each method
methods <- c("none", "bonferroni", "holm")
fwer_results <- numeric(length(methods))
names(fwer_results) <- methods

for (m in methods) {
  false_positives <- 0
  for (s in 1:n_sims_s2) {
    if (simulate_multi_variant(n_per_variant, p_control, k_variants, alpha_s2, m)) {
      false_positives <- false_positives + 1
    }
  }
  fwer_results[m] <- false_positives / n_sims_s2
}

cat("Family-Wise Error Rate (target alpha = 0.05):\n")
for (m in methods) {
  cat(sprintf("  %-15s: FWER = %.4f\n", m, fwer_results[m]))
}
cat("\n")

# --- 2c. Bayesian Hierarchical Model Simulation ---
cat("--- 2c. Bayesian Hierarchical Model ---\n")
cat("    Simulating hierarchical shrinkage effect on false positive rate\n\n")

# In a Bayesian hierarchical model, variant effects share a common prior:
# theta_j ~ Normal(mu, tau)
# This partial pooling shrinks estimates toward the group mean,
# naturally controlling for multiple comparisons

n_sims_hier <- 1000
bayes_fp <- 0

for (s in 1:n_sims_hier) {
  # Simulate data under null (all variants same rate)
  y <- rbinom(k_variants, n_per_variant, p_control)

  # Empirical Bayes: estimate group mean and variance
  p_hat <- y / n_per_variant
  mu_hat <- mean(p_hat)
  tau2_hat <- max(0, var(p_hat) - mean(p_hat * (1 - p_hat) / n_per_variant))

  # Shrink estimates toward group mean
  # Posterior mean = weighted average of data and prior
  if (tau2_hat > 0) {
    shrinkage <- (p_hat * (1/n_per_variant) / (p_hat*(1-p_hat)/n_per_variant) +
                    mu_hat / tau2_hat) /
                 (1/n_per_variant / (p_hat*(1-p_hat)/n_per_variant) + 1/tau2_hat)
    # Simplified shrinkage: use James-Stein-like estimator
    B <- p_hat * (1 - p_hat) / n_per_variant / (tau2_hat + p_hat * (1 - p_hat) / n_per_variant)
    shrink_est <- (1 - B) * p_hat + B * mu_hat
  } else {
    shrink_est <- rep(mu_hat, k_variants)
  }

  # Check if any treatment credibly different from control
  # Using simple posterior interval approach
  for (j in 2:k_variants) {
    # Posterior variance for difference
    se_diff <- sqrt(shrink_est[1]*(1-shrink_est[1])/n_per_variant +
                      shrink_est[j]*(1-shrink_est[j])/n_per_variant)
    diff_est <- shrink_est[j] - shrink_est[1]

    # 95% credible interval excludes zero?
    if (abs(diff_est) > 1.96 * se_diff) {
      bayes_fp <- bayes_fp + 1
      break
    }
  }
}

fwer_bayes <- bayes_fp / n_sims_hier
cat(sprintf("Bayesian hierarchical FWER: %.4f\n", fwer_bayes))
cat(sprintf("  (Shrinkage naturally reduces false positives vs uncorrected)\n\n"))

# FWER comparison data
fwer_df <- data.frame(
  Method = factor(c("Uncorrected", "Bonferroni", "Holm", "Bayes Hierarchical"),
                  levels = c("Uncorrected", "Bonferroni", "Holm", "Bayes Hierarchical")),
  FWER = c(fwer_results["none"], fwer_results["bonferroni"],
           fwer_results["holm"], fwer_bayes)
)


# ============================================================================
# SCENARIO 3: Marketplace Driver Incentive (Switchback Experiment)
# ============================================================================

cat("=== SCENARIO 3: Marketplace Driver Incentive (Switchback) ===\n\n")

# --- Parameters ---
n_days <- 30
hours_per_day <- 24
total_hours <- n_days * hours_per_day  # 720 hours
block_size <- 2  # 2-hour blocks
n_blocks_per_day <- hours_per_day / block_size  # 12 blocks per day
total_blocks <- n_days * n_blocks_per_day  # 360 blocks

true_effect <- 2.00   # $2.00 per driver per hour
n_drivers <- 500      # average active drivers per hour
base_earnings <- 25.00 # base hourly earnings per driver
sigma_hourly <- 3.00  # hourly noise SD
rho_ar <- 0.6         # temporal autocorrelation (AR(1))

cat(sprintf("Design parameters:\n"))
cat(sprintf("  Duration: %d days (%d hours, %d blocks of %d hours)\n",
            n_days, total_hours, total_blocks, block_size))
cat(sprintf("  True treatment effect: $%.2f per driver per hour\n", true_effect))
cat(sprintf("  Average active drivers: %d per hour\n", n_drivers))
cat(sprintf("  Temporal autocorrelation (AR(1) rho): %.2f\n\n", rho_ar))

# --- 3a. Simulate Hourly Data with Switchback Design ---
cat("--- 3a. Simulating 30 days of hourly data ---\n")

# Generate switchback assignment: alternating blocks
block_assignment <- rep(0:1, length.out = total_blocks)  # 0=control, 1=treatment

# Assign treatment to hours based on blocks
hourly_treatment <- rep(block_assignment, each = block_size)

# Generate AR(1) noise
ar_noise <- numeric(total_hours)
ar_noise[1] <- rnorm(1, 0, sigma_hourly)
for (t in 2:total_hours) {
  ar_noise[t] <- rho_ar * ar_noise[t-1] + rnorm(1, 0, sigma_hourly * sqrt(1 - rho_ar^2))
}

# Generate outcome: hourly average earnings per driver
hourly_earnings <- base_earnings +
  hourly_treatment * true_effect +
  ar_noise +
  # Add time-of-day effect (higher during rush hours)
  2.0 * sin(2 * pi * rep(0:23, n_days) / 24 - pi/2) +
  # Weekend effect
  ifelse(rep(c(rep(0, 5), rep(1, 2)), length.out = n_days)[
    rep(1:n_days, each = 24)] == 1, -1.5, 0)

# Create data frame
hour_data <- data.frame(
  hour = 1:total_hours,
  day = rep(1:n_days, each = 24),
  hour_of_day = rep(1:24, n_days),
  treatment = hourly_treatment,
  earnings = hourly_earnings,
  block_id = rep(1:total_blocks, each = block_size)
)

cat(sprintf("  Data generated: %d hourly observations\n", nrow(hour_data)))
cat(sprintf("  Treatment hours: %d (%.1f%%)\n",
            sum(hour_data$treatment == 1),
            100 * mean(hour_data$treatment)))
cat(sprintf("  Mean earnings (control):   $%.2f\n",
            mean(hour_data$earnings[hour_data$treatment == 0])))
cat(sprintf("  Mean earnings (treatment): $%.2f\n",
            mean(hour_data$earnings[hour_data$treatment == 1])))
cat(sprintf("  Observed difference: $%.2f (true: $%.2f)\n\n",
            mean(hour_data$earnings[hour_data$treatment == 1]) -
            mean(hour_data$earnings[hour_data$treatment == 0]),
            true_effect))

# --- 3b. Three Estimators ---
cat("--- 3b. Comparing Estimators ---\n\n")

# Estimator 1: Naive t-test (ignores switchback, assumes iid)
# Compute treatment effect as treatment mean - control mean
naive_test <- t.test(hour_data$earnings[hour_data$treatment == 1],
                     hour_data$earnings[hour_data$treatment == 0])
naive_estimate <- naive_test$estimate[1] - naive_test$estimate[2]
naive_se <- naive_test$stderr
cat(sprintf("1. Naive t-test (ignores switchback):\n"))
cat(sprintf("   Estimate: $%.3f\n", naive_estimate))
cat(sprintf("   SE:       $%.3f\n", naive_se))
cat(sprintf("   95%% CI:   [$%.3f, $%.3f]\n",
            naive_test$conf.int[1], naive_test$conf.int[2]))
cat(sprintf("   p-value:  %.6f\n", naive_test$p.value))
cat(sprintf("   WARNING: SE underestimated due to ignored autocorrelation!\n\n"))

# Estimator 2: Block-aggregated estimator
# Average within each block, then compare treatment vs control blocks
block_means <- tapply(hour_data$earnings, hour_data$block_id, mean)
block_treatment <- tapply(hour_data$treatment, hour_data$block_id, function(x) x[1])

block_data <- data.frame(
  block_id = as.numeric(names(block_means)),
  mean_earnings = as.numeric(block_means),
  treatment = as.numeric(block_treatment)
)

block_data$treatment_f <- factor(block_data$treatment, levels = c(0, 1))
block_test <- t.test(block_data$mean_earnings[block_data$treatment == 1],
                     block_data$mean_earnings[block_data$treatment == 0])
block_estimate <- block_test$estimate[1] - block_test$estimate[2]
block_se <- block_test$stderr
cat(sprintf("2. Block-aggregated estimator:\n"))
cat(sprintf("   Estimate: $%.3f\n", block_estimate))
cat(sprintf("   SE:       $%.3f\n", block_se))
cat(sprintf("   95%% CI:   [$%.3f, $%.3f]\n",
            block_test$conf.int[1], block_test$conf.int[2]))
cat(sprintf("   p-value:  %.6f\n", block_test$p.value))
cat(sprintf("   N blocks: %d (control=%d, treatment=%d)\n\n",
            nrow(block_data),
            sum(block_data$treatment == 0),
            sum(block_data$treatment == 1)))

# Estimator 3: HAC-robust standard errors (Newey-West)
# Manual implementation of Newey-West HAC estimator
cat("3. HAC-robust (Newey-West) estimator:\n")

# OLS: earnings = beta_0 + beta_1 * treatment + epsilon
X <- cbind(1, hour_data$treatment)
y <- hour_data$earnings
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
residuals <- y - X %*% beta_hat
n_obs <- length(y)
k_params <- 2

# Standard OLS SE (for comparison)
sigma2_ols <- sum(residuals^2) / (n_obs - k_params)
var_ols <- sigma2_ols * solve(t(X) %*% X)
se_ols <- sqrt(diag(var_ols))

# Newey-West HAC with bandwidth = floor(4*(n/100)^(2/9))
bandwidth <- floor(4 * (n_obs/100)^(2/9))
cat(sprintf("   Bandwidth (Newey-West): %d lags\n", bandwidth))

# Meat of the sandwich estimator
meat <- matrix(0, k_params, k_params)
# X'e e'X component (contemporaneous)
for (t in 1:n_obs) {
  xt <- matrix(X[t,], ncol = 1)
  meat <- meat + residuals[t]^2 * (xt %*% t(xt))
}

# Add autocovariance terms with Bartlett kernel
for (j in 1:bandwidth) {
  weight <- 1 - j / (bandwidth + 1)  # Bartlett kernel
  gamma_j <- matrix(0, k_params, k_params)
  for (t in (j+1):n_obs) {
    xt <- matrix(X[t,], ncol = 1)
    xt_j <- matrix(X[t-j,], ncol = 1)
    gamma_j <- gamma_j + residuals[t] * residuals[t-j] * (xt %*% t(xt_j))
  }
  meat <- meat + weight * (gamma_j + t(gamma_j))
}

bread <- solve(t(X) %*% X)
var_hac <- bread %*% meat %*% bread
se_hac <- sqrt(diag(var_hac))

hac_estimate <- beta_hat[2]
hac_se <- se_hac[2]
hac_ci_lower <- hac_estimate - 1.96 * hac_se
hac_ci_upper <- hac_estimate + 1.96 * hac_se
hac_z <- hac_estimate / hac_se
hac_pvalue <- 2 * pnorm(-abs(hac_z))

cat(sprintf("   Estimate: $%.3f\n", hac_estimate))
cat(sprintf("   SE:       $%.3f (OLS SE was $%.3f)\n", hac_se, se_ols[2]))
cat(sprintf("   SE inflation factor: %.2fx\n", hac_se / se_ols[2]))
cat(sprintf("   95%% CI:   [$%.3f, $%.3f]\n", hac_ci_lower, hac_ci_upper))
cat(sprintf("   p-value:  %.6f\n\n", hac_pvalue))

# Summary comparison
cat("--- Estimator Comparison Summary ---\n")
cat(sprintf("  %-25s  %-10s  %-10s  %-20s\n", "Method", "Estimate", "SE", "95% CI"))
cat(sprintf("  %-25s  $%-9.3f  $%-9.3f  [$%.3f, $%.3f]\n",
            "Naive t-test", naive_estimate, naive_se,
            naive_test$conf.int[1], naive_test$conf.int[2]))
cat(sprintf("  %-25s  $%-9.3f  $%-9.3f  [$%.3f, $%.3f]\n",
            "Block-aggregated", block_estimate, block_se,
            block_test$conf.int[1], block_test$conf.int[2]))
cat(sprintf("  %-25s  $%-9.3f  $%-9.3f  [$%.3f, $%.3f]\n",
            "HAC (Newey-West)", hac_estimate, hac_se,
            hac_ci_lower, hac_ci_upper))
cat(sprintf("\n  True effect: $%.2f\n\n", true_effect))


# ============================================================================
# CROSS-SCENARIO COMPARISON
# ============================================================================

cat("=== CROSS-SCENARIO COMPARISON ===\n\n")

# Comparison table
cat("--- Design Comparison Table ---\n")
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Feature", "S1: E-comm A/B", "S2: SaaS Multi", "S3: Marketplace"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "--------------------", "--------------------",
            "--------------------", "--------------------"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Randomization", "User-level", "User-level", "Time-block"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Variants", "2 (A/B)", "4 (1 ctrl + 3 trt)", "2 (on/off)"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Primary Metric", "Conversion Rate", "Trial-to-Paid", "$/Driver/Hour"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "MDE", "0.01 (1pp)", "0.03 (3pp)", "$2.00"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "n per group", sprintf("%d", n_required), sprintf("%d", n_required_s2_bonf),
            sprintf("%d blocks", total_blocks)))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Key Challenge", "Seasonality", "Multiple Comparisons", "Spillover/SUTVA"))
cat(sprintf("  %-20s  %-20s  %-20s  %-20s\n",
            "Analysis", "Beta-Binomial", "Hierarchical", "Block + HAC"))
cat("\n")


# ============================================================================
# FIGURE 1: Power Curves for All Scenarios
# ============================================================================

cat("Generating Figure 1: Power curves...\n")

# Generate finer grid for smooth curves
fine_sizes <- seq(200, 12000, by = 200)
freq_curve <- sapply(fine_sizes, function(n) freq_power(n, p_A, p_B, alpha))

# For S2, generate power curve with Bonferroni correction
freq_curve_s2 <- sapply(fine_sizes, function(n) {
  freq_power(n, p_control, p_control + mde_s2, alpha_bonf)
})

# For S3, approximate power based on block-level analysis
# With block_size=2, n_blocks = n_hours/2, SE scales as sigma/sqrt(n_blocks)
# Approximate: use block-level t-test power
power_s3 <- sapply(fine_sizes, function(n_hours) {
  n_blocks_est <- n_hours / block_size
  se_block <- sigma_hourly * sqrt(2 / n_blocks_est)  # block-level SE
  delta <- true_effect
  z_crit <- qnorm(1 - alpha/2)
  power <- pnorm(delta/se_block - z_crit) + pnorm(-delta/se_block - z_crit)
  min(power, 1.0)
})

power_long <- data.frame(
  n = rep(fine_sizes, 3),
  power = c(freq_curve, freq_curve_s2, power_s3),
  Scenario = factor(rep(c("S1: E-commerce A/B", "S2: SaaS Multi-variant",
                           "S3: Marketplace Switchback"),
                         each = length(fine_sizes)),
                     levels = c("S1: E-commerce A/B", "S2: SaaS Multi-variant",
                                "S3: Marketplace Switchback"))
)

p1 <- ggplot(power_long, aes(x = n, y = power, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  # Add simulation points for S1
  geom_point(data = data.frame(n = sample_sizes, power = freq_powers,
                                Scenario = factor("S1: E-commerce A/B",
                                                  levels = levels(power_long$Scenario))),
             aes(x = n, y = power), size = 2.5, shape = 16) +
  geom_point(data = data.frame(n = sample_sizes, power = bayes_powers,
                                Scenario = factor("S1: E-commerce A/B",
                                                  levels = levels(power_long$Scenario))),
             aes(x = n, y = power), size = 2.5, shape = 17) +
  annotate("text", x = 9000, y = 0.82, label = "80% power", color = "gray40", size = 3) +
  annotate("text", x = 9000, y = 0.07, label = "alpha = 0.05", color = "gray60", size = 3) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(x = "Sample Size per Group (or Hours for S3)",
       y = "Power",
       title = "Figure 1: Power Curves Across Experiment Designs",
       subtitle = sprintf("S1: p_A=%.2f, MDE=%.2f | S2: p=%.2f, MDE=%.2f (Bonferroni) | S3: effect=$%.0f",
                           p_A, mde, p_control, mde_s2, true_effect),
       color = "Scenario", linetype = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, color = "gray40"))

ggsave(paste0(fig_dir, "fig1_power_curves.png"), p1, width = 10, height = 6, dpi = 300)
cat("  Saved fig1_power_curves.png\n\n")


# ============================================================================
# FIGURE 2: Multiple Comparison FWER
# ============================================================================

cat("Generating Figure 2: Multiple comparison FWER...\n")

# Also simulate FWER at different numbers of variants
n_variants_seq <- c(2, 3, 4, 5, 6, 8, 10)
fwer_matrix <- matrix(NA, nrow = length(n_variants_seq), ncol = 3,
                      dimnames = list(NULL, c("Uncorrected", "Bonferroni", "Holm")))

n_sims_fwer <- 500
for (vi in seq_along(n_variants_seq)) {
  k <- n_variants_seq[vi]
  for (m in c("none", "bonferroni", "holm")) {
    fp <- 0
    for (s in 1:n_sims_fwer) {
      if (simulate_multi_variant(1500, 0.15, k, 0.05, m)) {
        fp <- fp + 1
      }
    }
    mname <- switch(m, "none" = "Uncorrected", "bonferroni" = "Bonferroni", "holm" = "Holm")
    fwer_matrix[vi, mname] <- fp / n_sims_fwer
  }
}

fwer_viz <- data.frame(
  n_variants = rep(n_variants_seq, 3),
  FWER = as.vector(fwer_matrix),
  Method = factor(rep(c("Uncorrected", "Bonferroni", "Holm"),
                       each = length(n_variants_seq)),
                   levels = c("Uncorrected", "Bonferroni", "Holm"))
)

p2 <- ggplot(fwer_viz, aes(x = factor(n_variants), y = FWER,
                             fill = Method, group = Method)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = 8.5, y = 0.06, label = "Nominal alpha = 0.05",
           color = "red", size = 3, hjust = 1) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_y_continuous(limits = c(0, 0.55)) +
  labs(x = "Number of Variants (including control)",
       y = "Family-Wise Error Rate (FWER)",
       title = "Figure 2: FWER Control Methods Comparison",
       subtitle = sprintf("Simulated with %d experiments per condition (null hypothesis true)", n_sims_fwer),
       fill = "Correction Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 10))

ggsave(paste0(fig_dir, "fig2_multiple_comparison.png"), p2, width = 10, height = 6, dpi = 300)
cat("  Saved fig2_multiple_comparison.png\n\n")


# ============================================================================
# FIGURE 3: Switchback Design Timeline and Estimates
# ============================================================================

cat("Generating Figure 3: Switchback design visualization...\n")

# Timeline plot (first 7 days for clarity)
hours_show <- 7 * 24
timeline_data <- hour_data[1:hours_show, ]
timeline_data$treatment_f <- factor(timeline_data$treatment,
                                     levels = c(0, 1),
                                     labels = c("Control", "Treatment"))

p3a <- ggplot(timeline_data, aes(x = hour, y = earnings, color = treatment_f)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  geom_rect(data = data.frame(xmin = seq(1, hours_show, by = block_size * 2),
                               xmax = seq(1, hours_show, by = block_size * 2) + block_size - 1,
                               ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.15, color = NA, inherit.aes = FALSE) +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_x_continuous(breaks = seq(0, hours_show, by = 24),
                     labels = paste0("Day ", 0:7)) +
  labs(x = "Time", y = "Hourly Avg Earnings ($)",
       color = "Assignment",
       title = "Switchback Timeline (First 7 Days)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# Estimator comparison plot
est_df <- data.frame(
  Method = factor(c("Naive t-test", "Block-aggregated", "HAC (Newey-West)"),
                  levels = c("Naive t-test", "Block-aggregated", "HAC (Newey-West)")),
  Estimate = c(naive_estimate, block_estimate, hac_estimate),
  Lower = c(naive_test$conf.int[1], block_test$conf.int[1], hac_ci_lower),
  Upper = c(naive_test$conf.int[2], block_test$conf.int[2], hac_ci_upper)
)

p3b <- ggplot(est_df, aes(x = Estimate, y = Method, color = Method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.2, linewidth = 0.8,
                orientation = "y") +
  geom_vline(xintercept = true_effect, linetype = "dashed", color = "darkred", linewidth = 0.7) +
  annotate("text", x = true_effect + 0.15, y = 3.3,
           label = sprintf("True effect = $%.2f", true_effect),
           color = "darkred", size = 3.5, hjust = 0) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(x = "Estimated Treatment Effect ($)",
       y = "",
       title = "Estimator Comparison",
       subtitle = "Dashed line = true effect; bars = 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Combine using grid arrangement (manually)
# Save as a two-panel figure
png(paste0(fig_dir, "fig3_switchback_design.png"), width = 12, height = 10,
    units = "in", res = 300)

# Manual 2-panel layout
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# Panel A: Timeline
plot(timeline_data$hour, timeline_data$earnings,
     type = "l", col = ifelse(timeline_data$treatment == 0, "#2166AC", "#B2182B"),
     xlab = "Time (hours)", ylab = "Hourly Avg Earnings ($)",
     main = "A) Switchback Timeline (First 7 Days)",
     lwd = 0.5, cex.main = 1.2)
# Add treatment blocks shading
for (b_start in seq(1, hours_show, by = block_size * 2)) {
  rect(b_start, par("usr")[3], b_start + block_size - 1, par("usr")[4],
       col = rgb(0.7, 0.85, 1, 0.2), border = NA)
}
# Re-draw lines on top
lines(timeline_data$hour[timeline_data$treatment == 0],
      timeline_data$earnings[timeline_data$treatment == 0],
      col = "#2166AC", lwd = 0.5)
lines(timeline_data$hour[timeline_data$treatment == 1],
      timeline_data$earnings[timeline_data$treatment == 1],
      col = "#B2182B", lwd = 0.5)
legend("topright", legend = c("Control", "Treatment"),
       col = c("#2166AC", "#B2182B"), lwd = 2, cex = 0.9)
# Day markers
abline(v = seq(0, hours_show, by = 24), col = "gray70", lty = 2)
axis(1, at = seq(0, hours_show, by = 24),
     labels = paste0("Day ", 0:7), tick = FALSE, line = -0.5, cex.axis = 0.8)

# Panel B: Estimator comparison
plot(1, 1, type = "n",
     xlim = c(min(est_df$Lower) - 0.5, max(est_df$Upper) + 0.5),
     ylim = c(0.5, 3.5),
     xlab = "Estimated Treatment Effect ($)", ylab = "",
     main = "B) Estimator Comparison",
     yaxt = "n", cex.main = 1.2)
axis(2, at = 1:3, labels = as.character(est_df$Method), las = 1, cex.axis = 0.9)
abline(v = true_effect, lty = 2, col = "darkred", lwd = 2)
colors_est <- c("#E41A1C", "#377EB8", "#4DAF4A")
for (i in 1:3) {
  points(est_df$Estimate[i], i, pch = 16, cex = 2, col = colors_est[i])
  segments(est_df$Lower[i], i, est_df$Upper[i], i, lwd = 3, col = colors_est[i])
}
mtext(sprintf("True effect = $%.2f", true_effect), side = 3, line = -1,
      at = true_effect + 0.8, cex = 0.8, col = "darkred")

dev.off()
cat("  Saved fig3_switchback_design.png\n\n")


# ============================================================================
# FIGURE 4: Sample Size Requirements Across Scenarios
# ============================================================================

cat("Generating Figure 4: Sample size comparison...\n")

# Required sample sizes for 80% power at various MDEs
mde_grid_s1 <- seq(0.005, 0.03, by = 0.001)
mde_grid_s2 <- seq(0.01, 0.06, by = 0.002)

n_req_s1 <- sapply(mde_grid_s1, function(d) {
  p_b <- p_A + d
  pv <- p_A*(1-p_A) + p_b*(1-p_b)
  ceiling((z_alpha + z_beta)^2 * pv / d^2)
})

n_req_s2 <- sapply(mde_grid_s2, function(d) {
  p_b <- p_control + d
  pv <- p_control*(1-p_control) + p_b*(1-p_b)
  ceiling((z_alpha_bonf + z_beta)^2 * pv / d^2)
})

# For S3: required hours for 80% power at various effect sizes
effect_grid_s3 <- seq(0.5, 5.0, by = 0.25)
n_req_s3 <- sapply(effect_grid_s3, function(d) {
  # Block-level: need n_blocks such that d / (sigma_block / sqrt(n_blocks/2)) > z_alpha + z_beta
  # sigma_block ~ sigma_hourly * sqrt(block_size) for block average
  sigma_block <- sigma_hourly * sqrt(block_size)
  n_blocks_needed <- 2 * ((z_alpha + z_beta) * sigma_block / d)^2
  n_hours_needed <- n_blocks_needed * block_size
  ceiling(n_hours_needed)
})

# Build the visualization
ss_df <- data.frame(
  MDE = c(mde_grid_s1, mde_grid_s2, effect_grid_s3),
  n_required = c(n_req_s1, n_req_s2, n_req_s3),
  Scenario = factor(c(rep("S1: E-commerce A/B", length(mde_grid_s1)),
                       rep("S2: SaaS Multi-variant", length(mde_grid_s2)),
                       rep("S3: Marketplace Switchback", length(effect_grid_s3))),
                     levels = c("S1: E-commerce A/B", "S2: SaaS Multi-variant",
                                "S3: Marketplace Switchback"))
)

p4 <- ggplot(ss_df, aes(x = MDE, y = n_required, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(data = data.frame(
    MDE = c(mde, mde_s2, true_effect),
    n_required = c(n_required, n_required_s2_bonf, total_hours),
    Scenario = factor(c("S1: E-commerce A/B", "S2: SaaS Multi-variant",
                         "S3: Marketplace Switchback"),
                       levels = levels(ss_df$Scenario))
  ), aes(x = MDE, y = n_required), size = 4, shape = 18) +
  scale_y_log10() +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(x = "Minimum Detectable Effect (proportion for S1/S2, $ for S3)",
       y = "Required Sample Size per Group (log scale)",
       title = "Figure 4: Sample Size Requirements vs Effect Size",
       subtitle = "Diamonds show the design MDE for each scenario (80% power)",
       color = "Scenario", linetype = "Scenario") +
  annotation_logticks(sides = "l") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))

ggsave(paste0(fig_dir, "fig4_sample_size_comparison.png"), p4, width = 10, height = 6, dpi = 300)
cat("  Saved fig4_sample_size_comparison.png\n\n")


# ============================================================================
# FIGURE 5: Design Trade-offs Visualization
# ============================================================================

cat("Generating Figure 5: Design trade-offs...\n")

# Create a trade-off comparison across multiple dimensions
tradeoff_df <- data.frame(
  Dimension = factor(rep(c("Statistical Power", "Implementation Cost",
                             "Internal Validity", "External Validity",
                             "Speed to Result", "Scalability"), 3),
                      levels = c("Statistical Power", "Implementation Cost",
                                 "Internal Validity", "External Validity",
                                 "Speed to Result", "Scalability")),
  Score = c(
    # S1: E-commerce A/B
    8, 3, 9, 6, 7, 8,
    # S2: SaaS Multi-variant
    6, 5, 7, 6, 5, 5,
    # S3: Marketplace Switchback
    5, 7, 6, 8, 4, 4
  ),
  Scenario = factor(rep(c("S1: E-commerce A/B",
                            "S2: SaaS Multi-variant",
                            "S3: Marketplace Switchback"),
                          each = 6),
                     levels = c("S1: E-commerce A/B",
                                "S2: SaaS Multi-variant",
                                "S3: Marketplace Switchback"))
)

p5 <- ggplot(tradeoff_df, aes(x = Dimension, y = Score, fill = Scenario)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  labs(x = "",
       y = "Score (1-10, higher = better)",
       title = "Figure 5: Experiment Design Trade-offs",
       subtitle = "Comparing three designs across six practical dimensions",
       fill = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(paste0(fig_dir, "fig5_design_tradeoffs.png"), p5, width = 10, height = 6, dpi = 300)
cat("  Saved fig5_design_tradeoffs.png\n\n")


# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("=============================================================\n")
cat("=== FINAL SUMMARY ===\n")
cat("=============================================================\n\n")

cat("Scenario 1: E-commerce Checkout A/B Test\n")
cat(sprintf("  - Design: User-level randomization, 2 variants\n"))
cat(sprintf("  - Required n: %d per group (freq), Bayesian sim agrees\n", n_required))
cat(sprintf("  - Bayesian power at n=5000: %.1f%% (ROPE-based)\n", bayes_powers[4]*100))
cat(sprintf("  - Freq power at n=5000: %.1f%%\n\n", freq_powers[4]*100))

cat("Scenario 2: SaaS Pricing Multi-Variant Test\n")
cat(sprintf("  - Design: 4 variants, Bonferroni correction\n"))
cat(sprintf("  - Required n: %d per variant (Bonferroni-adjusted)\n", n_required_s2_bonf))
cat(sprintf("  - FWER: uncorrected=%.3f, Bonferroni=%.3f, Holm=%.3f\n",
            fwer_results["none"], fwer_results["bonferroni"], fwer_results["holm"]))
cat(sprintf("  - Bayesian hierarchical FWER: %.3f\n\n", fwer_bayes))

cat("Scenario 3: Marketplace Driver Incentive (Switchback)\n")
cat(sprintf("  - Design: %d 2-hour blocks over %d days\n", total_blocks, n_days))
cat(sprintf("  - True effect: $%.2f\n", true_effect))
cat(sprintf("  - Naive estimate: $%.3f (SE=$%.3f) -- anti-conservative\n",
            naive_estimate, naive_se))
cat(sprintf("  - Block-aggregated: $%.3f (SE=$%.3f) -- appropriate\n",
            block_estimate, block_se))
cat(sprintf("  - HAC-robust: $%.3f (SE=$%.3f) -- robust to autocorrelation\n\n",
            hac_estimate, hac_se))

cat("Figures generated:\n")
cat("  1. fig1_power_curves.png\n")
cat("  2. fig2_multiple_comparison.png\n")
cat("  3. fig3_switchback_design.png\n")
cat("  4. fig4_sample_size_comparison.png\n")
cat("  5. fig5_design_tradeoffs.png\n\n")

cat("=== Complete ===\n")
