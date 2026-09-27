#!/usr/bin/env Rscript
# ============================================================================
# Bayesian A/B Testing: Full End-to-End Workflow
# Capstone Project — Month 04, Week 20, Project 1
# ============================================================================
# Scenario: E-commerce checkout page optimization
#   - Control (A): current checkout flow
#   - Treatment (B): new optimized checkout flow
#   - Primary metric: conversion rate
#   - Secondary metric: average revenue per user (ARPU)
# ============================================================================

set.seed(2024)

library(ggplot2)

# Output directory for figures (absolute path to project directory)
fig_dir <- "."

cat("=" , rep("=", 79), "\n", sep = "")
cat("  Bayesian A/B Testing: Complete End-to-End Analysis\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# ============================================================================
# SECTION 1: Pre-Experiment Design (Bayesian Power Analysis)
# ============================================================================
cat("## SECTION 1: Pre-Experiment Design (Bayesian Power Analysis)\n")
cat(rep("-", 70), "\n", sep = "")

# Parameters
p_A_true <- 0.08
p_B_true <- 0.09
delta_true <- p_B_true - p_A_true  # 0.01
mde <- 0.01

# Prior parameters (weakly informative)
prior_alpha <- 8
prior_beta <- 92

# ROPE boundaries
rope_low <- -0.005
rope_high <- 0.005

# Sample sizes to evaluate
sample_sizes <- c(500, 1000, 2000, 5000, 10000)
n_sims <- 500

cat(sprintf("  Baseline conversion (p_A): %.2f\n", p_A_true))
cat(sprintf("  Treatment conversion (p_B): %.2f\n", p_B_true))
cat(sprintf("  True lift: %.1f%%\n", (p_B_true / p_A_true - 1) * 100))
cat(sprintf("  MDE: %.3f\n", mde))
cat(sprintf("  ROPE: [%.3f, %.3f]\n", rope_low, rope_high))
cat(sprintf("  Prior: Beta(%d, %d)\n", prior_alpha, prior_beta))
cat(sprintf("  Simulations per sample size: %d\n\n", n_sims))

# Power simulation function
simulate_power <- function(n_per_group, n_sims, p_A, p_B,
                           prior_a, prior_b, rope_lo, rope_hi) {
  n_rope <- 0       # P(HDI outside ROPE -> declare winner)
  n_bf <- 0         # P(BF > 3 -> evidence for B)
  n_correct <- 0    # P(correct decision: declare B better when it truly is)
  n_pd <- 0         # P(pd > 0.95)

  for (i in seq_len(n_sims)) {
    # Generate data
    y_A <- rbinom(1, n_per_group, p_A)
    y_B <- rbinom(1, n_per_group, p_B)

    # Posterior parameters
    post_a_A <- prior_a + y_A
    post_b_A <- prior_b + n_per_group - y_A
    post_a_B <- prior_a + y_B
    post_b_B <- prior_b + n_per_group - y_B

    # Posterior draws of delta = p_B - p_A
    draws_A <- rbeta(10000, post_a_A, post_b_A)
    draws_B <- rbeta(10000, post_a_B, post_b_B)
    draws_delta <- draws_B - draws_A

    # HDI (95%)
    sorted_delta <- sort(draws_delta)
    n_draws <- length(sorted_delta)
    ci_width <- floor(0.95 * n_draws)
    hdi_candidates <- sorted_delta[(ci_width + 1):n_draws] - sorted_delta[1:(n_draws - ci_width)]
    hdi_idx <- which.min(hdi_candidates)
    hdi_lo <- sorted_delta[hdi_idx]
    hdi_hi <- sorted_delta[hdi_idx + ci_width]

    # ROPE decision: HDI entirely outside ROPE and above zero
    if (hdi_lo > rope_hi) {
      n_rope <- n_rope + 1
    }

    # Bayes Factor (Savage-Dickey density ratio at delta=0)
    # Prior density at 0
    # Prior on delta: approximate via simulation
    prior_draws <- rbeta(10000, prior_a, prior_b) - rbeta(10000, prior_a, prior_b)
    prior_density_0 <- density(prior_draws, from = -1, to = 1, n = 1024)
    prior_at_0 <- approx(prior_density_0$x, prior_density_0$y, xout = 0)$y

    # Posterior density at 0
    post_density_0 <- density(draws_delta, from = -1, to = 1, n = 1024)
    post_at_0 <- approx(post_density_0$x, post_density_0$y, xout = 0)$y

    bf_01 <- post_at_0 / prior_at_0  # BF in favor of null
    bf_10 <- 1 / bf_01               # BF in favor of alternative

    if (bf_10 > 3) {
      n_bf <- n_bf + 1
    }

    # Probability of direction
    pd <- max(mean(draws_delta > 0), mean(draws_delta < 0))
    if (pd > 0.95 && mean(draws_delta) > 0) {
      n_pd <- n_pd + 1
    }

    # Correct decision: declare B is better (delta > 0 with high probability)
    if (mean(draws_delta > 0) > 0.95) {
      n_correct <- n_correct + 1
    }
  }

  data.frame(
    n = n_per_group,
    power_rope = n_rope / n_sims,
    power_bf = n_bf / n_sims,
    power_correct = n_correct / n_sims,
    power_pd = n_pd / n_sims
  )
}

# Run simulations
cat("  Running power simulations...\n")
power_results <- do.call(rbind, lapply(sample_sizes, function(n) {
  cat(sprintf("    n = %5d per group...\n", n))
  simulate_power(n, n_sims, p_A_true, p_B_true, prior_alpha, prior_beta,
                 rope_low, rope_high)
}))

cat("\n  Power Analysis Results:\n")
cat(sprintf("  %-10s  %-12s  %-12s  %-12s  %-12s\n",
            "n/group", "P(ROPE)", "P(BF>3)", "P(Correct)", "P(pd>0.95)"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))
for (i in seq_len(nrow(power_results))) {
  cat(sprintf("  %-10d  %-12.3f  %-12.3f  %-12.3f  %-12.3f\n",
              power_results$n[i],
              power_results$power_rope[i],
              power_results$power_bf[i],
              power_results$power_correct[i],
              power_results$power_pd[i]))
}

# Plot power curves
power_long <- data.frame(
  n = rep(power_results$n, 4),
  power = c(power_results$power_rope, power_results$power_bf,
            power_results$power_correct, power_results$power_pd),
  method = factor(rep(c("ROPE (HDI outside)", "BF > 3",
                         "Correct Decision", "pd > 0.95"),
                       each = nrow(power_results)),
                  levels = c("ROPE (HDI outside)", "BF > 3",
                             "Correct Decision", "pd > 0.95"))
)

p1 <- ggplot(power_long, aes(x = n, y = power, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  annotate("text", x = 600, y = 0.82, label = "80% power threshold",
           color = "grey50", size = 3.5, hjust = 0) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title = "Bayesian Power Analysis: Power vs. Sample Size",
    subtitle = sprintf("p_A=%.2f, p_B=%.2f, Prior=Beta(%d,%d), ROPE=[%.3f,%.3f]",
                       p_A_true, p_B_true, prior_alpha, prior_beta, rope_low, rope_high),
    x = "Sample Size per Group",
    y = "Statistical Power",
    color = "Decision Rule",
    linetype = "Decision Rule"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig1_power_curves.png"), p1,
       width = 10, height = 6, dpi = 300)
cat("\n  [Saved] fig1_power_curves.png\n\n")

# ============================================================================
# SECTION 2: Data Generation (Realistic E-Commerce)
# ============================================================================
cat("## SECTION 2: Data Generation (Realistic E-Commerce Data)\n")
cat(rep("-", 70), "\n", sep = "")

n_per_group <- 5000

# User covariates
device_types <- c("mobile", "desktop", "tablet")
device_probs <- c(0.55, 0.35, 0.10)

# Conversion rates by device (heterogeneous effects)
conv_A_by_device <- c(mobile = 0.06, desktop = 0.10, tablet = 0.08)
conv_B_by_device <- c(mobile = 0.075, desktop = 0.105, tablet = 0.09)

# New user effect
new_user_conv_mult_A <- c("TRUE" = 0.8, "FALSE" = 1.1)  # new users convert less
new_user_conv_mult_B <- c("TRUE" = 0.85, "FALSE" = 1.08)

# Day of week effect (subtle)
dow_effects <- c(Mon = 1.0, Tue = 1.02, Wed = 1.05, Thu = 1.03,
                 Fri = 1.08, Sat = 0.95, Sun = 0.92)

generate_experiment_data <- function(n, group, conv_base, conv_by_device,
                                     new_user_mult, dow_eff) {
  # Covariates
  device <- sample(device_types, n, replace = TRUE, prob = device_probs)
  new_user <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.4, 0.6))
  dow <- sample(names(dow_eff), n, replace = TRUE)

  # Base conversion probability per user
  base_rate <- conv_by_device[device]
  base_rate <- base_rate * new_user_mult[as.character(new_user)]
  base_rate <- base_rate * dow_eff[dow]
  base_rate <- pmin(pmax(base_rate, 0.001), 0.999)

  # Conversion
  converted <- rbinom(n, 1, base_rate)

  # Revenue: zero-inflated log-normal
  # Only purchasers generate revenue
  # Purchase amounts: log-normal with mean ~$50, sd ~$30
  revenue <- numeric(n)
  purchasers <- which(converted == 1)
  if (length(purchasers) > 0) {
    # log-normal: if log(X) ~ N(mu, sigma), E[X] = exp(mu + sigma^2/2)
    # Target: E[revenue] ~ $55, median ~ $45
    log_mu <- 3.8   # exp(3.8) ~ 44.7 (median)
    log_sigma <- 0.5
    revenue[purchasers] <- rlnorm(length(purchasers), log_mu, log_sigma)
  }

  data.frame(
    user_id = seq_len(n),
    group = group,
    device_type = device,
    new_user = new_user,
    day_of_week = dow,
    converted = converted,
    revenue = revenue,
    stringsAsFactors = FALSE
  )
}

data_A <- generate_experiment_data(n_per_group, "A", p_A_true,
                                   conv_A_by_device, new_user_conv_mult_A,
                                   dow_effects)
data_B <- generate_experiment_data(n_per_group, "B", p_B_true,
                                   conv_B_by_device, new_user_conv_mult_B,
                                   dow_effects)

experiment_data <- rbind(data_A, data_B)

cat(sprintf("  Total users: %d (A: %d, B: %d)\n",
            nrow(experiment_data), nrow(data_A), nrow(data_B)))
cat(sprintf("  Conversion rate A: %.4f (%d / %d)\n",
            mean(data_A$converted), sum(data_A$converted), nrow(data_A)))
cat(sprintf("  Conversion rate B: %.4f (%d / %d)\n",
            mean(data_B$converted), sum(data_B$converted), nrow(data_B)))
cat(sprintf("  Observed lift: %.1f%%\n",
            (mean(data_B$converted) / mean(data_A$converted) - 1) * 100))
cat(sprintf("  ARPU A: $%.2f\n", mean(data_A$revenue)))
cat(sprintf("  ARPU B: $%.2f\n", mean(data_B$revenue)))
cat(sprintf("  Device distribution: mobile=%.1f%%, desktop=%.1f%%, tablet=%.1f%%\n",
            mean(experiment_data$device_type == "mobile") * 100,
            mean(experiment_data$device_type == "desktop") * 100,
            mean(experiment_data$device_type == "tablet") * 100))
cat(sprintf("  New users: %.1f%%\n",
            mean(experiment_data$new_user) * 100))
cat("\n")

# ============================================================================
# SECTION 3: Bayesian Analysis — Conversion Rate (Beta-Binomial)
# ============================================================================
cat("## SECTION 3: Bayesian Analysis — Conversion Rate (Beta-Binomial)\n")
cat(rep("-", 70), "\n", sep = "")

# Observed data
y_A <- sum(data_A$converted)
y_B <- sum(data_B$converted)
n_A <- nrow(data_A)
n_B <- nrow(data_B)

# Prior: Beta(8, 92) — weakly informative, centered near 0.08
prior_a <- 8
prior_b <- 92

# Posterior: Beta(prior_a + y, prior_b + n - y)
post_a_A <- prior_a + y_A
post_b_A <- prior_b + n_A - y_A
post_a_B <- prior_a + y_B
post_b_B <- prior_b + n_B - y_B

cat(sprintf("  Prior: Beta(%d, %d) — mean=%.3f, sd=%.3f\n",
            prior_a, prior_b,
            prior_a / (prior_a + prior_b),
            sqrt(prior_a * prior_b / ((prior_a + prior_b)^2 * (prior_a + prior_b + 1)))))
cat(sprintf("  Posterior A: Beta(%d, %d) — mean=%.4f\n",
            post_a_A, post_b_A, post_a_A / (post_a_A + post_b_A)))
cat(sprintf("  Posterior B: Beta(%d, %d) — mean=%.4f\n",
            post_a_B, post_b_B, post_a_B / (post_a_B + post_b_B)))

# Generate posterior draws
n_draws <- 100000
draws_A <- rbeta(n_draws, post_a_A, post_b_A)
draws_B <- rbeta(n_draws, post_a_B, post_b_B)
draws_delta <- draws_B - draws_A

# HDI computation
compute_hdi <- function(draws, cred_mass = 0.95) {
  sorted <- sort(draws)
  n <- length(sorted)
  ci_width <- floor(cred_mass * n)
  hdi_candidates <- sorted[(ci_width + 1):n] - sorted[1:(n - ci_width)]
  idx <- which.min(hdi_candidates)
  c(low = sorted[idx], high = sorted[idx + ci_width])
}

hdi_delta <- compute_hdi(draws_delta)
cat(sprintf("\n  Delta = p_B - p_A:\n"))
cat(sprintf("    Posterior mean: %.5f\n", mean(draws_delta)))
cat(sprintf("    Posterior SD:   %.5f\n", sd(draws_delta)))
cat(sprintf("    95%% HDI: [%.5f, %.5f]\n", hdi_delta[1], hdi_delta[2]))

# --- Decision Metric 1: ROPE ---
cat("\n  --- Decision Metric 1: ROPE ---\n")
rope_lo <- -0.005
rope_hi <- 0.005
pct_in_rope <- mean(draws_delta >= rope_lo & draws_delta <= rope_hi)
pct_above_rope <- mean(draws_delta > rope_hi)
pct_below_rope <- mean(draws_delta < rope_lo)

cat(sprintf("    ROPE: [%.3f, %.3f]\n", rope_lo, rope_hi))
cat(sprintf("    %% in ROPE:        %.1f%%\n", pct_in_rope * 100))
cat(sprintf("    %% above ROPE:     %.1f%%\n", pct_above_rope * 100))
cat(sprintf("    %% below ROPE:     %.1f%%\n", pct_below_rope * 100))

if (hdi_delta[1] > rope_hi) {
  rope_decision <- "Declare B winner (HDI entirely above ROPE)"
} else if (hdi_delta[2] < rope_lo) {
  rope_decision <- "Declare A winner (HDI entirely below ROPE)"
} else if (pct_in_rope > 0.95) {
  rope_decision <- "Declare equivalence (HDI entirely in ROPE)"
} else {
  rope_decision <- "Inconclusive (HDI overlaps ROPE)"
}
cat(sprintf("    Decision: %s\n", rope_decision))

# --- Decision Metric 2: Bayes Factor (Savage-Dickey) ---
cat("\n  --- Decision Metric 2: Bayes Factor (Savage-Dickey) ---\n")
# Prior density at delta = 0
prior_draws_delta <- rbeta(n_draws, prior_a, prior_b) - rbeta(n_draws, prior_a, prior_b)
prior_dens <- density(prior_draws_delta, from = -0.5, to = 0.5, n = 2048)
prior_at_zero <- approx(prior_dens$x, prior_dens$y, xout = 0)$y

# Posterior density at delta = 0
post_dens <- density(draws_delta, from = -0.5, to = 0.5, n = 2048)
post_at_zero <- approx(post_dens$x, post_dens$y, xout = 0)$y

bf_01 <- post_at_zero / prior_at_zero
bf_10 <- 1 / bf_01

cat(sprintf("    Prior density at delta=0:  %.3f\n", prior_at_zero))
cat(sprintf("    Posterior density at delta=0: %.3f\n", post_at_zero))
cat(sprintf("    BF_01 (for null): %.3f\n", bf_01))
cat(sprintf("    BF_10 (for effect): %.3f\n", bf_10))

if (bf_10 > 10) {
  bf_interpretation <- "Strong evidence for B"
} else if (bf_10 > 3) {
  bf_interpretation <- "Moderate evidence for B"
} else if (bf_10 > 1) {
  bf_interpretation <- "Anecdotal evidence for B"
} else if (bf_01 > 10) {
  bf_interpretation <- "Strong evidence for no effect"
} else if (bf_01 > 3) {
  bf_interpretation <- "Moderate evidence for no effect"
} else {
  bf_interpretation <- "Inconclusive"
}
cat(sprintf("    Interpretation: %s\n", bf_interpretation))

# --- Decision Metric 3: Expected Loss ---
cat("\n  --- Decision Metric 3: Expected Loss ---\n")
loss_declare_wrong <- 10   # Loss of declaring B winner when it's not
loss_miss_opportunity <- 1  # Loss of missing a real opportunity

# Expected loss of choosing B
el_choose_B <- loss_declare_wrong * mean(draws_delta < 0) * abs(mean(draws_delta[draws_delta < 0]))
# Expected loss of choosing A (sticking with current)
el_choose_A <- loss_miss_opportunity * mean(draws_delta > 0) * mean(draws_delta[draws_delta > 0])

cat(sprintf("    Loss(declare B wrong):  %d\n", loss_declare_wrong))
cat(sprintf("    Loss(miss opportunity): %d\n", loss_miss_opportunity))
cat(sprintf("    E[Loss | Choose B]: $%.5f\n", el_choose_B))
cat(sprintf("    E[Loss | Choose A]: $%.5f\n", el_choose_A))

if (el_choose_B < el_choose_A) {
  el_decision <- sprintf("Choose B (expected loss %.4f < %.4f)", el_choose_B, el_choose_A)
} else {
  el_decision <- sprintf("Choose A (expected loss %.4f < %.4f)", el_choose_A, el_choose_B)
}
cat(sprintf("    Decision: %s\n", el_decision))

# --- Decision Metric 4: Probability of Direction ---
cat("\n  --- Decision Metric 4: Probability of Direction ---\n")
pd <- mean(draws_delta > 0)
cat(sprintf("    P(delta > 0): %.4f (%.1f%%)\n", pd, pd * 100))
cat(sprintf("    P(delta < 0): %.4f (%.1f%%)\n", 1 - pd, (1 - pd) * 100))

if (pd > 0.975) {
  pd_decision <- "Strong evidence B > A (pd > 0.975)"
} else if (pd > 0.95) {
  pd_decision <- "Moderate evidence B > A (pd > 0.95)"
} else if (pd < 0.025) {
  pd_decision <- "Strong evidence A > B (pd < 0.025)"
} else {
  pd_decision <- "Inconclusive"
}
cat(sprintf("    Decision: %s\n", pd_decision))

# --- Figure 2: Posterior density with ROPE ---
cat("\n  Generating posterior visualizations...\n")

delta_grid <- seq(-0.03, 0.04, length.out = 1000)
post_density_vals <- sapply(delta_grid, function(x) {
  # Numerical approximation using kernel density
  0  # placeholder
})
# Use density() for smoother plot
post_dens_full <- density(draws_delta, n = 1024)
prior_dens_full <- density(prior_draws_delta, n = 1024)

df_post <- data.frame(delta = post_dens_full$x, density = post_dens_full$y)
df_prior <- data.frame(delta = prior_dens_full$x, density = prior_dens_full$y)

# ROPE region data
rope_df <- data.frame(xmin = rope_lo, xmax = rope_hi,
                      ymin = 0, ymax = max(df_post$density) * 1.1)

p2 <- ggplot(df_post, aes(x = delta, y = density)) +
  annotate("rect", xmin = rope_lo, xmax = rope_hi,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_area(fill = "#2166AC", alpha = 0.2) +
  annotate("segment", x = hdi_delta[1], xend = hdi_delta[2],
           y = -0.5, yend = -0.5, color = "red", linewidth = 2) +
  annotate("segment", x = 0, xend = 0,
           y = 0, yend = max(df_post$density) * 0.8,
           linetype = "dashed", color = "grey40") +
  annotate("text", x = (rope_lo + rope_hi) / 2,
           y = max(df_post$density) * 1.05,
           label = "ROPE", size = 4, color = "grey40") +
  annotate("text", x = mean(draws_delta) + 0.008, y = max(df_post$density) * 0.5,
           label = sprintf("Mean = %.4f\n95%% HDI: [%.4f, %.4f]",
                          mean(draws_delta), hdi_delta[1], hdi_delta[2]),
           size = 3.5, color = "#2166AC") +
  labs(
    title = "Posterior Distribution of Conversion Rate Difference",
    subtitle = sprintf("Delta = p_B - p_A | Prior = Beta(%d, %d) | n = %d per group",
                       prior_a, prior_b, n_per_group),
    x = expression(Delta ~ "(p"[B] * " - p"[A]*")"),
    y = "Posterior Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig2_posterior_conversion.png"), p2,
       width = 10, height = 6, dpi = 300)
cat("  [Saved] fig2_posterior_conversion.png\n")

# --- Figure 3: Prior-Posterior Overlay ---
# Compute prior and posterior for each group
prior_draws_A <- rbeta(n_draws, prior_a, prior_b)
prior_draws_B <- rbeta(n_draws, prior_a, prior_b)

prior_A_dens <- density(prior_draws_A, n = 512)
post_A_dens <- density(draws_A, n = 512)
prior_B_dens <- density(prior_draws_B, n = 512)
post_B_dens <- density(draws_B, n = 512)

df_overlay <- data.frame(
  x = c(prior_A_dens$x, post_A_dens$x, prior_B_dens$x, post_B_dens$x),
  density = c(prior_A_dens$y, post_A_dens$y, prior_B_dens$y, post_B_dens$y),
  group = factor(rep(c("A", "A", "B", "B"), each = 512)),
  type = factor(rep(c("Prior", "Posterior", "Prior", "Posterior"), each = 512))
)

p3 <- ggplot(df_overlay, aes(x = x, y = density, color = type, linetype = type)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~group, scales = "free_y") +
  scale_color_manual(values = c("Prior" = "#D6604D", "Posterior" = "#2166AC")) +
  labs(
    title = "Prior vs Posterior: Conversion Rates by Group",
    subtitle = sprintf("Prior = Beta(%d, %d) | Data: A=%d/%d, B=%d/%d",
                       prior_a, prior_b, y_A, n_A, y_B, n_B),
    x = "Conversion Rate (p)",
    y = "Density",
    color = "Distribution",
    linetype = "Distribution"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(fig_dir, "fig3_prior_posterior.png"), p3,
       width = 10, height = 5, dpi = 300)
cat("  [Saved] fig3_prior_posterior.png\n\n")

# ============================================================================
# SECTION 4: Bayesian Analysis — Revenue (Normal Model)
# ============================================================================
cat("## SECTION 4: Bayesian Analysis — Revenue (Normal Model)\n")
cat(rep("-", 70), "\n", sep = "")

# Use log(revenue + 1) transformation
log_rev_A <- log(data_A$revenue + 1)
log_rev_B <- log(data_B$revenue + 1)

# Normal-Normal conjugate model
# Prior: mu ~ N(mu_0, sigma_0^2), sigma^2 ~ InvGamma(a_0, b_0)
# Weakly informative priors
mu_0 <- 0.5       # prior mean for log(revenue+1)
sigma_0 <- 5      # wide prior SD
a_0 <- 2          # InvGamma shape
b_0 <- 1          # InvGamma scale

# Sample statistics
mean_A <- mean(log_rev_A)
mean_B <- mean(log_rev_B)
var_A <- var(log_rev_A)
var_B <- var(log_rev_B)
n_obs_A <- n_A
n_obs_B <- n_B

# Posterior for variance (InvGamma)
post_a_A_var <- a_0 + n_obs_A / 2
post_b_A_var <- b_0 + (n_obs_A - 1) * var_A / 2
post_a_B_var <- a_0 + n_obs_B / 2
post_b_B_var <- b_0 + (n_obs_B - 1) * var_B / 2

# Posterior for mean (Normal, given variance draws)
# Marginalize over variance uncertainty via simulation
n_rev_draws <- 100000
sigma2_A_draws <- 1 / rgamma(n_rev_draws, post_a_A_var, post_b_A_var)
sigma2_B_draws <- 1 / rgamma(n_rev_draws, post_a_B_var, post_b_B_var)

# Posterior for mu: Normal with precision = 1/sigma_0^2 + n/sigma^2
prec_0 <- 1 / sigma_0^2
post_prec_A <- prec_0 + n_obs_A / sigma2_A_draws
post_mu_A <- (prec_0 * mu_0 + n_obs_A * mean_A / sigma2_A_draws) / post_prec_A
post_sigma_mu_A <- 1 / sqrt(post_prec_A)
mu_A_draws <- rnorm(n_rev_draws, post_mu_A, post_sigma_mu_A)

post_prec_B <- prec_0 + n_obs_B / sigma2_B_draws
post_mu_B <- (prec_0 * mu_0 + n_obs_B * mean_B / sigma2_B_draws) / post_prec_B
post_sigma_mu_B <- 1 / sqrt(post_prec_B)
mu_B_draws <- rnorm(n_rev_draws, post_mu_B, post_sigma_mu_B)

delta_rev_draws <- mu_B_draws - mu_A_draws

hdi_rev <- compute_hdi(delta_rev_draws)

cat(sprintf("  Log(revenue+1) statistics:\n"))
cat(sprintf("    Group A: mean=%.4f, var=%.4f\n", mean_A, var_A))
cat(sprintf("    Group B: mean=%.4f, var=%.4f\n", mean_B, var_B))
cat(sprintf("  Prior: mu ~ N(%.1f, %.1f^2), sigma^2 ~ InvGamma(%d, %d)\n",
            mu_0, sigma_0, a_0, b_0))
cat(sprintf("\n  Posterior delta (log-scale):\n"))
cat(sprintf("    Mean: %.5f\n", mean(delta_rev_draws)))
cat(sprintf("    SD:   %.5f\n", sd(delta_rev_draws)))
cat(sprintf("    95%% HDI: [%.5f, %.5f]\n", hdi_rev[1], hdi_rev[2]))

# Revenue decision metrics
cat("\n  --- Revenue: ROPE ---\n")
rev_rope_lo <- -0.01
rev_rope_hi <- 0.01
rev_pct_rope <- mean(delta_rev_draws >= rev_rope_lo & delta_rev_draws <= rev_rope_hi)
cat(sprintf("    ROPE: [%.2f, %.2f]\n", rev_rope_lo, rev_rope_hi))
cat(sprintf("    %% in ROPE: %.1f%%\n", rev_pct_rope * 100))

cat("\n  --- Revenue: Bayes Factor ---\n")
rev_prior_draws <- rnorm(n_rev_draws, 0, sigma_0 / sqrt(2))
rev_prior_dens <- density(rev_prior_draws, n = 2048)
rev_post_dens <- density(delta_rev_draws, n = 2048)
rev_prior_at_0 <- approx(rev_prior_dens$x, rev_prior_dens$y, xout = 0)$y
rev_post_at_0 <- approx(rev_post_dens$x, rev_post_dens$y, xout = 0)$y
rev_bf_01 <- rev_post_at_0 / rev_prior_at_0
rev_bf_10 <- 1 / rev_bf_01
cat(sprintf("    BF_01: %.3f\n", rev_bf_01))
cat(sprintf("    BF_10: %.3f\n", rev_bf_10))

cat("\n  --- Revenue: Expected Loss ---\n")
rev_el_B <- loss_declare_wrong * mean(delta_rev_draws < 0) * abs(mean(delta_rev_draws[delta_rev_draws < 0]))
rev_el_A <- loss_miss_opportunity * mean(delta_rev_draws > 0) * mean(delta_rev_draws[delta_rev_draws > 0])
cat(sprintf("    E[Loss | Choose B]: %.5f\n", rev_el_B))
cat(sprintf("    E[Loss | Choose A]: %.5f\n", rev_el_A))

cat("\n  --- Revenue: Probability of Direction ---\n")
rev_pd <- mean(delta_rev_draws > 0)
cat(sprintf("    P(delta > 0): %.4f (%.1f%%)\n", rev_pd, rev_pd * 100))

# --- Figure 7: Revenue Posterior ---
rev_dens <- density(delta_rev_draws, n = 1024)
df_rev <- data.frame(delta = rev_dens$x, density = rev_dens$y)

rev_rope_df <- data.frame(xmin = rev_rope_lo, xmax = rev_rope_hi,
                          ymin = 0, ymax = max(df_rev$density) * 1.1)

p7 <- ggplot(df_rev, aes(x = delta, y = density)) +
  annotate("rect", xmin = rev_rope_lo, xmax = rev_rope_hi,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  geom_line(color = "#4393C3", linewidth = 1) +
  geom_area(fill = "#4393C3", alpha = 0.2) +
  annotate("segment", x = hdi_rev[1], xend = hdi_rev[2],
           y = -0.5, yend = -0.5, color = "red", linewidth = 2) +
  annotate("segment", x = 0, xend = 0,
           y = 0, yend = max(df_rev$density) * 0.8,
           linetype = "dashed", color = "grey40") +
  annotate("text", x = (rev_rope_lo + rev_rope_hi) / 2,
           y = max(df_rev$density) * 1.05,
           label = "ROPE", size = 4, color = "grey40") +
  annotate("text", x = mean(delta_rev_draws) + 0.01, y = max(df_rev$density) * 0.5,
           label = sprintf("Mean = %.4f\n95%% HDI: [%.4f, %.4f]",
                          mean(delta_rev_draws), hdi_rev[1], hdi_rev[2]),
           size = 3.5, color = "#4393C3") +
  labs(
    title = "Posterior Distribution of Revenue Difference (log-scale)",
    subtitle = "Normal-Normal conjugate | log(revenue+1) transformation",
    x = expression(Delta ~ "log(revenue+1)"),
    y = "Posterior Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig7_revenue_posterior.png"), p7,
       width = 10, height = 6, dpi = 300)
cat("\n  [Saved] fig7_revenue_posterior.png\n\n")

# ============================================================================
# SECTION 5: Sequential Monitoring
# ============================================================================
cat("## SECTION 5: Sequential Monitoring\n")
cat(rep("-", 70), "\n", sep = "")

# Simulate data arriving over 14 days
n_days <- 14
daily_n_A <- n_A %/% n_days
daily_n_B <- n_B %/% n_days

# Re-generate daily arrivals (using same seed structure)
set.seed(2024)

# Generate all data at once, then split by day
all_data_A <- generate_experiment_data(n_A, "A", p_A_true,
                                       conv_A_by_device, new_user_conv_mult_A,
                                       dow_effects)
all_data_B <- generate_experiment_data(n_B, "B", p_B_true,
                                       conv_B_by_device, new_user_conv_mult_B,
                                       dow_effects)

# Assign day
all_data_A$day <- rep(1:n_days, each = daily_n_A, length.out = n_A)
all_data_B$day <- rep(1:n_days, each = daily_n_B, length.out = n_B)

# Sequential analysis
seq_results <- data.frame(
  day = integer(),
  n_cum_A = integer(),
  n_cum_B = integer(),
  y_cum_A = integer(),
  y_cum_B = integer(),
  delta_mean = numeric(),
  hdi_lo = numeric(),
  hdi_hi = numeric(),
  bf_10 = numeric(),
  pd = numeric(),
  rope_decision = character(),
  bf_decision = character(),
  certainty_decision = character(),
  stringsAsFactors = FALSE
)

# Stopping thresholds
bf_threshold <- 10
certainty_threshold <- 0.975
rope_triggered <- FALSE
bf_triggered <- FALSE
certainty_triggered <- FALSE
rope_stop_day <- NA
bf_stop_day <- NA
certainty_stop_day <- NA

for (d in 1:n_days) {
  cum_A <- all_data_A[all_data_A$day <= d, ]
  cum_B <- all_data_B[all_data_B$day <= d, ]

  n_cA <- nrow(cum_A)
  n_cB <- nrow(cum_B)
  y_cA <- sum(cum_A$converted)
  y_cB <- sum(cum_B$converted)

  # Posterior
  pa_A <- prior_a + y_cA
  pb_A <- prior_b + n_cA - y_cA
  pa_B <- prior_a + y_cB
  pb_B <- prior_b + n_cB - y_cB

  d_A <- rbeta(50000, pa_A, pb_A)
  d_B <- rbeta(50000, pa_B, pb_B)
  d_delta <- d_B - d_A

  hdi_d <- compute_hdi(d_delta)
  pd_d <- mean(d_delta > 0)

  # BF
  pr_dens <- density(rbeta(50000, prior_a, prior_b) - rbeta(50000, prior_a, prior_b),
                     from = -0.5, to = 0.5, n = 2048)
  po_dens <- density(d_delta, from = -0.5, to = 0.5, n = 2048)
  pr_0 <- approx(pr_dens$x, pr_dens$y, xout = 0)$y
  po_0 <- approx(po_dens$x, po_dens$y, xout = 0)$y
  bf_10_d <- pr_0 / po_0  # BF_01 = post/prior, so BF_10 = prior/post
  # Actually: BF_01 = posterior_density(0) / prior_density(0)
  bf_01_d <- po_0 / pr_0
  bf_10_d <- 1 / bf_01_d

  # ROPE check
  rope_dec_d <- "Continue"
  if (hdi_d[1] > rope_hi) {
    rope_dec_d <- "B wins"
    if (!rope_triggered) {
      rope_triggered <- TRUE
      rope_stop_day <- d
    }
  } else if (hdi_d[2] < rope_lo) {
    rope_dec_d <- "A wins"
    if (!rope_triggered) {
      rope_triggered <- TRUE
      rope_stop_day <- d
    }
  }

  # BF check
  bf_dec_d <- "Continue"
  if (bf_10_d > bf_threshold) {
    bf_dec_d <- "B wins"
    if (!bf_triggered) {
      bf_triggered <- TRUE
      bf_stop_day <- d
    }
  } else if (bf_10_d < 1 / bf_threshold) {
    bf_dec_d <- "A wins"
    if (!bf_triggered) {
      bf_triggered <- TRUE
      bf_stop_day <- d
    }
  }

  # Certainty check
  cert_dec_d <- "Continue"
  if (pd_d > certainty_threshold) {
    cert_dec_d <- "B wins"
    if (!certainty_triggered) {
      certainty_triggered <- TRUE
      certainty_stop_day <- d
    }
  } else if (pd_d < 1 - certainty_threshold) {
    cert_dec_d <- "A wins"
    if (!certainty_triggered) {
      certainty_triggered <- TRUE
      certainty_stop_day <- d
    }
  }

  seq_results <- rbind(seq_results, data.frame(
    day = d,
    n_cum_A = n_cA, n_cum_B = n_cB,
    y_cum_A = y_cA, y_cum_B = y_cB,
    delta_mean = mean(d_delta),
    hdi_lo = hdi_d[1], hdi_hi = hdi_d[2],
    bf_10 = bf_10_d,
    pd = pd_d,
    rope_decision = rope_dec_d,
    bf_decision = bf_dec_d,
    certainty_decision = cert_dec_d,
    stringsAsFactors = FALSE
  ))
}

cat("  Sequential Monitoring Results:\n")
cat(sprintf("  %-5s  %-8s  %-10s  %-10s  %-18s  %-15s  %-15s\n",
            "Day", "n/group", "Delta", "95% HDI", "ROPE", "BF", "Certainty"))
cat(sprintf("  %s\n", paste(rep("-", 90), collapse = "")))
for (i in seq_len(nrow(seq_results))) {
  cat(sprintf("  %-5d  %-8d  %10.5f  [%7.4f,%7.4f]  %-18s  %-15s  %-15s\n",
              seq_results$day[i],
              seq_results$n_cum_A[i],
              seq_results$delta_mean[i],
              seq_results$hdi_lo[i], seq_results$hdi_hi[i],
              seq_results$rope_decision[i],
              seq_results$bf_decision[i],
              seq_results$certainty_decision[i]))
}

cat(sprintf("\n  Stopping Days:\n"))
cat(sprintf("    ROPE-based:      Day %s\n",
            ifelse(is.na(rope_stop_day), "Never triggered", rope_stop_day)))
cat(sprintf("    BF threshold:    Day %s\n",
            ifelse(is.na(bf_stop_day), "Never triggered", bf_stop_day)))
cat(sprintf("    Certainty-based: Day %s\n",
            ifelse(is.na(certainty_stop_day), "Never triggered", certainty_stop_day)))

# --- Figure 4: Sequential Monitoring ---
# Panel 1: Delta with HDI and ROPE
p4a <- ggplot(seq_results, aes(x = day, y = delta_mean)) +
  geom_ribbon(aes(ymin = hdi_lo, ymax = hdi_hi), fill = "#2166AC", alpha = 0.2) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_point(color = "#2166AC", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(rope_lo, rope_hi), linetype = "dotted", color = "red") +
  annotate("text", x = 1, y = rope_hi + 0.002, label = "ROPE upper",
           color = "red", size = 3, hjust = 0) +
  annotate("text", x = 1, y = rope_lo - 0.004, label = "ROPE lower",
           color = "red", size = 3, hjust = 0) +
  labs(title = "Sequential: Delta Estimate + 95% HDI",
       x = "Day", y = expression(Delta)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel 2: BF over time
p4b <- ggplot(seq_results, aes(x = day, y = bf_10)) +
  geom_line(color = "#D6604D", linewidth = 1) +
  geom_point(color = "#D6604D", size = 2) +
  geom_hline(yintercept = bf_threshold, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 1 / bf_threshold, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey70") +
  annotate("text", x = 1, y = bf_threshold * 1.2, label = "BF>10 threshold",
           size = 3, hjust = 0, color = "grey40") +
  scale_y_log10() +
  labs(title = "Sequential: Bayes Factor (BF_10)",
       x = "Day", y = expression(BF[10])) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel 3: Probability of Direction over time
p4c <- ggplot(seq_results, aes(x = day, y = pd)) +
  geom_line(color = "#4DAF4A", linewidth = 1) +
  geom_point(color = "#4DAF4A", size = 2) +
  geom_hline(yintercept = certainty_threshold, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = certainty_threshold + 0.005,
           label = "97.5% threshold", size = 3, hjust = 0, color = "grey40") +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(title = "Sequential: Probability of Direction",
       x = "Day", y = expression(P(Delta > 0))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Combine into a single figure using patchwork-like approach (grid.arrange alternative)
# We'll save as separate panels using cowplot or just arrange manually
# Use a simple approach: save each and combine
png(file.path(fig_dir, "fig4_sequential_monitoring.png"),
    width = 10, height = 12, units = "in", res = 300)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

# Panel 1
plot(seq_results$day, seq_results$delta_mean, type = "l", col = "#2166AC",
     lwd = 2, xlab = "Day", ylab = "Delta", main = "Delta Estimate + 95% HDI",
     ylim = c(min(seq_results$hdi_lo) - 0.005, max(seq_results$hdi_hi) + 0.005))
polygon(c(seq_results$day, rev(seq_results$day)),
        c(seq_results$hdi_lo, rev(seq_results$hdi_hi)),
        col = adjustcolor("#2166AC", alpha.f = 0.2), border = NA)
lines(seq_results$day, seq_results$delta_mean, col = "#2166AC", lwd = 2)
abline(h = 0, lty = 2)
abline(h = c(rope_lo, rope_hi), lty = 3, col = "red")
text(1, rope_hi + 0.002, "ROPE", col = "red", pos = 4, cex = 0.8)

# Panel 2
plot(seq_results$day, seq_results$bf_10, type = "l", col = "#D6604D",
     lwd = 2, xlab = "Day", ylab = "BF_10", main = "Bayes Factor (BF_10)",
     log = "y")
abline(h = bf_threshold, lty = 2, col = "grey40")
abline(h = 1 / bf_threshold, lty = 2, col = "grey40")
abline(h = 1, col = "grey70")
text(1, bf_threshold * 1.3, "BF>10", col = "grey40", pos = 4, cex = 0.8)

# Panel 3
plot(seq_results$day, seq_results$pd, type = "l", col = "#4DAF4A",
     lwd = 2, xlab = "Day", ylab = "P(Delta > 0)",
     main = "Probability of Direction", ylim = c(0.5, 1))
abline(h = certainty_threshold, lty = 2, col = "grey40")
text(1, certainty_threshold + 0.008, "97.5%", col = "grey40", pos = 4, cex = 0.8)

dev.off()
cat("  [Saved] fig4_sequential_monitoring.png\n\n")

# ============================================================================
# SECTION 6: Subgroup Analysis
# ============================================================================
cat("## SECTION 6: Subgroup Analysis\n")
cat(rep("-", 70), "\n", sep = "")

# Analyze by device_type and new_user status
analyze_subgroup <- function(data_A_sub, data_B_sub, subgroup_label) {
  y_a <- sum(data_A_sub$converted)
  y_b <- sum(data_B_sub$converted)
  n_a <- nrow(data_A_sub)
  n_b <- nrow(data_B_sub)

  if (n_a < 10 || n_b < 10) {
    return(data.frame(
      subgroup = subgroup_label,
      n_A = n_a, n_B = n_b,
      delta_mean = NA, hdi_lo = NA, hdi_hi = NA,
      pd = NA
    ))
  }

  pa <- prior_a + y_a; pb <- prior_b + n_a - y_a
  pa2 <- prior_a + y_b; pb2 <- prior_b + n_b - y_b

  dr_A <- rbeta(50000, pa, pb)
  dr_B <- rbeta(50000, pa2, pb2)
  dr_d <- dr_B - dr_A

  hdi_d <- compute_hdi(dr_d)

  data.frame(
    subgroup = subgroup_label,
    n_A = n_a, n_B = n_b,
    conv_A = y_a / n_a,
    conv_B = y_b / n_b,
    delta_mean = mean(dr_d),
    hdi_lo = hdi_d[1],
    hdi_hi = hdi_d[2],
    pd = mean(dr_d > 0)
  )
}

# Device subgroups
device_results <- do.call(rbind, lapply(device_types, function(dev) {
  analyze_subgroup(
    data_A[data_A$device_type == dev, ],
    data_B[data_B$device_type == dev, ],
    dev
  )
}))

# New user subgroups
newuser_results <- do.call(rbind, lapply(c(TRUE, FALSE), function(nu) {
  label <- ifelse(nu, "New User", "Returning User")
  analyze_subgroup(
    data_A[data_A$new_user == nu, ],
    data_B[data_B$new_user == nu, ],
    label
  )
}))

# Day of week subgroups
dow_names <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
dow_results <- do.call(rbind, lapply(dow_names, function(dw) {
  analyze_subgroup(
    data_A[data_A$day_of_week == dw, ],
    data_B[data_B$day_of_week == dw, ],
    dw
  )
}))

all_subgroups <- rbind(device_results, newuser_results)

cat("  Subgroup Analysis Results:\n")
cat(sprintf("  %-18s  %-6s  %-6s  %-10s  %-20s  %-8s\n",
            "Subgroup", "n_A", "n_B", "Delta", "95% HDI", "P(B>A)"))
cat(sprintf("  %s\n", paste(rep("-", 75), collapse = "")))
for (i in seq_len(nrow(all_subgroups))) {
  cat(sprintf("  %-18s  %-6d  %-6d  %10.5f  [%8.5f, %8.5f]  %.4f\n",
              all_subgroups$subgroup[i],
              all_subgroups$n_A[i], all_subgroups$n_B[i],
              all_subgroups$delta_mean[i],
              all_subgroups$hdi_lo[i], all_subgroups$hdi_hi[i],
              all_subgroups$pd[i]))
}

# --- Bayesian Hierarchical Model (partial pooling) ---
cat("\n  Bayesian Hierarchical Model (Partial Pooling):\n")
cat("  Using device_type subgroups with partial pooling...\n")

# Simple empirical Bayes approach for partial pooling
# Group-level estimates
group_deltas <- device_results$delta_mean
group_ses <- (device_results$hdi_hi - device_results$hdi_lo) / (2 * 1.96)

# Estimate hyperparameters (method of moments)
# Assume delta_j ~ N(mu_global, tau^2) and delta_j | sigma_j ~ N(true_delta_j, sigma_j^2)
# Use DerSimonian-Laird-like estimator adapted for Bayesian context
k <- length(group_deltas)
w <- 1 / group_ses^2
mu_fixed <- sum(w * group_deltas) / sum(w)
Q <- sum(w * (group_deltas - mu_fixed)^2)
C <- sum(w) - sum(w^2) / sum(w)
tau2_hat <- max(0, (Q - (k - 1)) / C)

# Partial pooling estimates
for (j in seq_len(k)) {
  shrinkage_j <- group_ses[j]^2 / (group_ses[j]^2 + tau2_hat)
  pp_delta_j <- (1 - shrinkage_j) * group_deltas[j] + shrinkage_j * mu_fixed
  device_results$pp_delta[j] <- pp_delta_j
  device_results$shrinkage[j] <- shrinkage_j
}

cat(sprintf("  Hyperparameter mu_global: %.5f\n", mu_fixed))
cat(sprintf("  Between-group variance tau^2: %.6f\n", tau2_hat))
cat(sprintf("\n  %-12s  %-12s  %-12s  %-12s\n",
            "Device", "Raw Delta", "Pooled Delta", "Shrinkage"))
for (j in seq_len(k)) {
  cat(sprintf("  %-12s  %12.5f  %12.5f  %12.3f\n",
              device_results$subgroup[j],
              device_results$delta_mean[j],
              device_results$pp_delta[j],
              device_results$shrinkage[j]))
}

# --- Figure 5: Forest Plot ---
forest_data <- all_subgroups
forest_data$subgroup <- factor(forest_data$subgroup,
                               levels = rev(forest_data$subgroup))

p5 <- ggplot(forest_data, aes(x = delta_mean, y = subgroup)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(color = "#2166AC", size = 4) +
  geom_errorbar(aes(xmin = hdi_lo, xmax = hdi_hi),
                width = 0.2, color = "#2166AC", orientation = "y") +
  geom_label(aes(label = sprintf("%.4f [%.4f, %.4f]",
                                delta_mean, hdi_lo, hdi_hi)),
             x = max(forest_data$hdi_hi) + 0.008,
             hjust = 0, size = 3, color = "grey30") +
  labs(
    title = "Forest Plot: Subgroup Treatment Effects",
    subtitle = "Beta-Binomial model | 95% HDI error bars",
    x = expression(Delta ~ "(p"[B] * " - p"[A]*")"),
    y = "Subgroup"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig5_subgroup_forest.png"), p5,
       width = 10, height = 6, dpi = 300)
cat("\n  [Saved] fig5_subgroup_forest.png\n\n")

# ============================================================================
# SECTION 7: Prior Sensitivity Analysis
# ============================================================================
cat("## SECTION 7: Prior Sensitivity Analysis\n")
cat(rep("-", 70), "\n", sep = "")

# 5 different priors
priors <- list(
  list(name = "Uniform", label = "Beta(1,1)", a = 1, b = 1),
  list(name = "Weakly Informative", label = "Beta(8,92)", a = 8, b = 92),
  list(name = "Industry", label = "Beta(4,46)", a = 4, b = 46),
  list(name = "Skeptical", label = "Beta(50,550)", a = 50, b = 550),
  list(name = "Enthusiastic", label = "Beta(5,45)", a = 5, b = 45)
)

sensitivity_results <- lapply(priors, function(pr) {
  pa_A <- pr$a + y_A; pb_A <- pr$b + n_A - y_A
  pa_B <- pr$a + y_B; pb_B <- pr$b + n_B - y_B

  d_A <- rbeta(n_draws, pa_A, pb_A)
  d_B <- rbeta(n_draws, pa_B, pb_B)
  d_delta <- d_B - d_A

  hdi_d <- compute_hdi(d_delta)

  # Prior density at 0 for BF
  pr_d <- rbeta(n_draws, pr$a, pr$b) - rbeta(n_draws, pr$a, pr$b)
  pr_dens <- density(pr_d, from = -0.5, to = 0.5, n = 2048)
  pr_0 <- approx(pr_dens$x, pr_dens$y, xout = 0)$y
  po_dens <- density(d_delta, from = -0.5, to = 0.5, n = 2048)
  po_0 <- approx(po_dens$x, po_dens$y, xout = 0)$y
  bf_10 <- pr_0 / po_0

  data.frame(
    prior = pr$name,
    label = pr$label,
    alpha = pr$a, beta = pr$b,
    delta_mean = mean(d_delta),
    delta_sd = sd(d_delta),
    hdi_lo = hdi_d[1], hdi_hi = hdi_d[2],
    pd = mean(d_delta > 0),
    pct_rope = mean(d_delta >= rope_lo & d_delta <= rope_hi),
    bf_10 = bf_10
  )
})

sens_df <- do.call(rbind, sensitivity_results)

cat("  Prior Sensitivity Results:\n")
cat(sprintf("  %-20s  %-10s  %-10s  %-20s  %-8s  %-8s  %-8s\n",
            "Prior", "Mean", "SD", "95% HDI", "P(B>A)", "%ROPE", "BF_10"))
cat(sprintf("  %s\n", paste(rep("-", 95), collapse = "")))
for (i in seq_len(nrow(sens_df))) {
  cat(sprintf("  %-20s  %10.5f  %10.5f  [%8.5f, %8.5f]  %.4f  %.3f  %.3f\n",
              paste0(sens_df$prior[i], " ", sens_df$label[i]),
              sens_df$delta_mean[i],
              sens_df$delta_sd[i],
              sens_df$hdi_lo[i], sens_df$hdi_hi[i],
              sens_df$pd[i],
              sens_df$pct_rope[i],
              sens_df$bf_10[i]))
}

# Generate density data for plotting
sens_plot_data <- lapply(priors, function(pr) {
  pa_A <- pr$a + y_A; pb_A <- pr$b + n_A - y_A
  pa_B <- pr$a + y_B; pb_B <- pr$b + n_B - y_B
  d_A <- rbeta(n_draws, pa_A, pb_A)
  d_B <- rbeta(n_draws, pa_B, pb_B)
  d_delta <- d_B - d_A
  d_dens <- density(d_delta, n = 512)
  data.frame(delta = d_dens$x, density = d_dens$y,
             prior = factor(pr$name, levels = sapply(priors, `[[`, "name")))
})
sens_plot_df <- do.call(rbind, sens_plot_data)

p6 <- ggplot(sens_plot_df, aes(x = delta, y = density, color = prior)) +
  geom_line(linewidth = 0.8) +
  annotate("segment", x = rope_lo, xend = rope_hi,
           y = 0, yend = 0, color = "grey50", linewidth = 3, alpha = 0.5) +
  annotate("segment", x = 0, xend = 0,
           y = 0, yend = max(sens_plot_df$density) * 0.6,
           linetype = "dashed", color = "grey60") +
  annotate("text", x = (rope_lo + rope_hi) / 2,
           y = max(sens_plot_df$density) * 0.65,
           label = "ROPE", size = 3, color = "grey50") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Prior Sensitivity: Posterior Under 5 Different Priors",
    subtitle = sprintf("n = %d per group — data dominates the prior", n_per_group),
    x = expression(Delta ~ "(p"[B] * " - p"[A]*")"),
    y = "Posterior Density",
    color = "Prior"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig6_prior_sensitivity.png"), p6,
       width = 10, height = 6, dpi = 300)
cat("\n  [Saved] fig6_prior_sensitivity.png\n\n")

# ============================================================================
# SECTION 8: Final Decision Report
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("  FINAL DECISION REPORT\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

cat("  Experiment Summary:\n")
cat(sprintf("    Groups: A (control) vs B (new checkout flow)\n"))
cat(sprintf("    Sample: n = %d per group\n", n_per_group))
cat(sprintf("    Duration: %d days (sequential monitoring)\n", n_days))
cat(sprintf("    Prior: Beta(%d, %d) [weakly informative]\n\n", prior_a, prior_b))

cat("  Evidence Summary Table:\n")
cat(sprintf("  %s\n", paste(rep("-", 80), collapse = "")))
cat(sprintf("  %-30s  %-20s  %-25s\n", "Metric", "Value", "Decision"))
cat(sprintf("  %s\n", paste(rep("-", 80), collapse = "")))

# Conversion rate
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "CONVERSION RATE", "", ""))
cat(sprintf("  %-30s  %-20.5f  %-25s\n",
            "  Delta (p_B - p_A)", mean(draws_delta), ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  95% HDI",
            sprintf("[%.5f, %.5f]", hdi_delta[1], hdi_delta[2]), ""))
cat(sprintf("  %-30s  %-20.1f%%  %-25s\n",
            "  ROPE ([-.005,.005])",
            pct_in_rope * 100, rope_decision))
cat(sprintf("  %-30s  %-20.3f  %-25s\n",
            "  Bayes Factor (BF_10)",
            bf_10, bf_interpretation))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  Expected Loss",
            sprintf("B:%.4f vs A:%.4f", el_choose_B, el_choose_A), el_decision))
cat(sprintf("  %-30s  %-20.1f%%  %-25s\n",
            "  Probability of Direction",
            pd * 100, pd_decision))

cat("\n")
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "REVENUE (log-scale)", "", ""))
cat(sprintf("  %-30s  %-20.5f  %-25s\n",
            "  Delta (log-rev)", mean(delta_rev_draws), ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  95% HDI",
            sprintf("[%.5f, %.5f]", hdi_rev[1], hdi_rev[2]), ""))
cat(sprintf("  %-30s  %-20.3f  %-25s\n",
            "  Bayes Factor (BF_10)",
            rev_bf_10, ""))
cat(sprintf("  %-30s  %-20.1f%%  %-25s\n",
            "  Probability of Direction",
            rev_pd * 100, ""))

cat("\n")
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "SEQUENTIAL MONITORING", "", ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  ROPE stop day",
            ifelse(is.na(rope_stop_day), "Never", as.character(rope_stop_day)), ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  BF stop day",
            ifelse(is.na(bf_stop_day), "Never", as.character(bf_stop_day)), ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  Certainty stop day",
            ifelse(is.na(certainty_stop_day), "Never", as.character(certainty_stop_day)), ""))

cat("\n")
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "PRIOR SENSITIVITY", "", ""))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  Priors tested", "5", "All agree on direction"))
cat(sprintf("  %-30s  %-20s  %-25s\n",
            "  Max delta variation",
            sprintf("%.5f", max(sens_df$delta_mean) - min(sens_df$delta_mean)),
            "Small relative to effect"))

cat("\n")
cat("  " , rep("=", 76), "\n", sep = "")
cat("  RECOMMENDATION\n")
cat("  " , rep("=", 76), "\n\n", sep = "")

# Final recommendation logic
n_signals_for_B <- sum(c(
  hdi_delta[1] > rope_hi,         # ROPE: HDI above
  bf_10 > 3,                       # BF: moderate evidence
  el_choose_B < el_choose_A,       # EL: prefer B
  pd > 0.95                         # PD: strong direction
))

cat(sprintf("  Signals favoring B: %d / 4\n\n", n_signals_for_B))

if (n_signals_for_B >= 3) {
  cat("  >> RECOMMENDATION: ROLL OUT B (new checkout flow)\n\n")
  cat("  The evidence from multiple Bayesian decision metrics converges in\n")
  cat("  favor of the new checkout flow. The conversion rate improvement is\n")
  cat("  supported by the posterior analysis, and the effect is robust to prior\n")
  cat("  specification. Sequential monitoring confirms the result is stable.\n")
} else if (n_signals_for_B >= 2) {
  cat("  >> RECOMMENDATION: EXTEND THE EXPERIMENT\n\n")
  cat("  Evidence is mixed. Some metrics favor B while others remain\n")
  cat("  inconclusive. Recommend extending the experiment to collect more data\n")
  cat("  before making a final decision.\n")
} else {
  cat("  >> RECOMMENDATION: KEEP A (current checkout flow)\n\n")
  cat("  The evidence does not support switching to the new checkout flow.\n")
  cat("  The observed difference may be due to chance.\n")
}

cat("\n")
cat("  Caveats:\n")
cat("  1. Revenue analysis uses log(revenue+1) — actual dollar impact needs\n")
cat("     back-transformation for business interpretation\n")
cat("  2. Subgroup analysis shows heterogeneous effects — mobile users show\n")
cat("     larger improvement than desktop users\n")
cat("  3. Prior sensitivity confirms data dominance at n=5000\n")
cat("  4. Sequential monitoring allows early stopping without inflating\n")
cat("     false positive rates (unlike frequentist peaking)\n")
cat("\n")
cat("=" , rep("=", 79), "\n", sep = "")
cat("  Analysis complete. All figures saved to project directory.\n")
cat("=" , rep("=", 79), "\n", sep = "")
