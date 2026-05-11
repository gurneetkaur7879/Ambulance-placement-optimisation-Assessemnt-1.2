# Simulated Annealing for the ambulance p-median model

library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(scales)

source("02_Shared_helpers.R")

p <- 10
coverage_threshold <- 5000
min_coverage <- 0.80
coverage_penalty <- 1e11
seeds <- 1:30

max_iter <- 3000
initial_temp <- 1e9
cooling_rate <- 0.998

mode_cols <- c(deterministic = "#0072B2", stochastic = "#E69F00")
extra_cols <- c("#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#D55E00")

if (p > ncol(distance_matrix)) stop("p cannot be larger than the number of candidate sites.")
stopifnot(coverage_penalty == 1e11)

# Select best result using strict 80 percent rule
select_best_strict <- function(results_df) {
  results_df %>%
    group_by(mode) %>%
    arrange(desc(meets_coverage), mean_distance, desc(coverage), .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup()
}

# Run one SA seed
run_sa <- function(mode, seed) {
  set.seed(seed)
  
  # Create initial solution
  current_metrics <- evaluate_full(
    generate_random_solution(p = p),
    mode = mode,
    p = p,
    threshold = coverage_threshold,
    min_coverage = min_coverage,
    penalty = coverage_penalty
  )
  
  current_solution <- current_metrics$solution
  current_objective <- current_metrics$penalised_objective
  best_solution <- current_solution
  best_objective <- current_objective
  
  # Store search history
  best_history <- numeric(max_iter)
  mean_history <- numeric(max_iter)
  coverage_history <- numeric(max_iter)
  temp_history <- numeric(max_iter)
  accepted_history <- logical(max_iter)
  
  temp <- initial_temp
  
  # Main SA loop
  for (iter in seq_len(max_iter)) {
    
    # Generate neighbour
    candidate_solution <- generate_neighbor(current_solution, p)
    
    # Evaluate neighbour
    candidate_objective <- evaluate_penalised(
      candidate_solution,
      mode = mode,
      p = p,
      threshold = coverage_threshold,
      min_coverage = min_coverage,
      penalty = coverage_penalty
    )
    
    # Calculate objective change
    delta <- candidate_objective - current_objective
    accepted <- FALSE
    
    # Accept better or probabilistic worse move
    if (delta < 0 || runif(1) < exp(-delta / temp)) {
      current_solution <- candidate_solution
      current_objective <- candidate_objective
      accepted <- TRUE
    }
    
    # Update best solution
    if (current_objective < best_objective) {
      best_solution <- current_solution
      best_objective <- current_objective
    }
    
    # Store iteration metrics
    best_history[iter] <- best_objective
    mean_history[iter] <- current_objective
    coverage_history[iter] <- compute_coverage(best_solution, coverage_threshold, mode = mode, p = p)
    temp_history[iter] <- temp
    accepted_history[iter] <- accepted
    
    # Cool temperature
    temp <- temp * cooling_rate
  }
  
  # Evaluate final best solution
  best_metrics <- evaluate_full(
    best_solution,
    mode = mode,
    p = p,
    threshold = coverage_threshold,
    min_coverage = min_coverage,
    penalty = coverage_penalty
  )
  
  # Return run result
  list(
    method = "SA",
    mode = mode,
    seed = seed,
    best_solution = best_metrics$solution,
    selected_sites = which(best_metrics$solution == 1L),
    weighted_distance = best_metrics$weighted_distance,
    mean_distance = best_metrics$mean_distance,
    penalised_objective = best_metrics$penalised_objective,
    coverage = best_metrics$coverage,
    best_history = best_history,
    mean_history = mean_history,
    coverage_history = coverage_history,
    temp_history = temp_history,
    accepted_history = accepted_history
  )
}

# Prepare storage
runs <- list()
rows <- list()
k <- 1

cat("\nRunning SA experiments...\n\n")

# Run deterministic and stochastic experiments
for (mode in c("deterministic", "stochastic")) {
  for (seed in seeds) {
    
    # Time each run
    start_time <- Sys.time()
    res <- run_sa(mode, seed)
    runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Check strict coverage
    meets_coverage <- res$coverage >= min_coverage
    
    # Print run result
    cat(
      toupper(mode), "Seed", seed,
      "| mean distance =", round(res$mean_distance / 1000, 2), "km",
      "| coverage =", sprintf("%.4f", res$coverage),
      "| target = 0.8000",
      "| feasible =", meets_coverage,
      "| runtime =", round(runtime, 2), "s\n"
    )
    
    # Store full object
    runs[[k]] <- res
    
    # Store table row
    rows[[k]] <- data.frame(
      method = "SA",
      mode = mode,
      seed = seed,
      weighted_distance = res$weighted_distance,
      mean_distance = res$mean_distance,
      penalised_objective = res$penalised_objective,
      coverage = res$coverage,
      meets_coverage = meets_coverage,
      runtime_seconds = runtime,
      selected_sites = paste(res$selected_sites, collapse = ", ")
    )
    
    k <- k + 1
  }
}

# Combine results
sa_results <- bind_rows(rows)

# Select best strict solution
best_sa_by_mode <- select_best_strict(sa_results)

# Save main outputs
saveRDS(runs, "sa_30_run_objects.rds")
write_csv(sa_results, "sa_30_run_results.csv")
#write_csv(best_sa_by_mode, "sa_best_solution_by_mode.csv")

# Summarise runs
sa_summary <- sa_results %>%
  group_by(method, mode) %>%
  summarise(
    mean_weighted_distance = mean(weighted_distance),
    sd_weighted_distance = sd(weighted_distance),
    mean_distance_km = mean(mean_distance) / 1000,
    sd_distance_km = sd(mean_distance) / 1000,
    best_weighted_distance = min(weighted_distance),
    worst_weighted_distance = max(weighted_distance),
    mean_coverage = mean(coverage),
    sd_coverage = sd(coverage),
    best_coverage = max(coverage),
    coverage_success_rate = mean(meets_coverage),
    mean_runtime = mean(runtime_seconds),
    sd_runtime = sd(runtime_seconds),
    .groups = "drop"
  )

# Save summary
write_csv(sa_summary, "sa_summary_table.csv")

# Save report-friendly summary
sa_summary_display <- sa_summary %>%
  select(
    method,
    mode,
    mean_distance_km,
    sd_distance_km,
    best_weighted_distance,
    worst_weighted_distance,
    mean_coverage,
    best_coverage,
    coverage_success_rate,
    mean_runtime
  )

#write_csv(sa_summary_display, "sa_summary_display.csv")

# Build convergence history
sa_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$best_history),
    best = x$best_history,
    mean = x$mean_history
  )
}))

# Average convergence history
sa_history_avg <- sa_history %>%
  group_by(mode, iteration) %>%
  summarise(
    mean_best = mean(best),
    sd_best = sd(best),
    ci_best = 1.96 * sd_best / sqrt(n()),
    .groups = "drop"
  )

#write_csv(sa_history, "sa_convergence_history_all_runs.csv")
#write_csv(sa_history_avg, "sa_convergence_summary.csv")

# Build coverage history
sa_coverage_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$coverage_history),
    coverage = x$coverage_history
  )
}))

# Average coverage history
sa_coverage_avg <- sa_coverage_history %>%
  group_by(mode, iteration) %>%
  summarise(
    mean_coverage = mean(coverage),
    sd_coverage = sd(coverage),
    ci_coverage = 1.96 * sd_coverage / sqrt(n()),
    .groups = "drop"
  )

#write_csv(sa_coverage_history, "sa_coverage_history_all_runs.csv")
#write_csv(sa_coverage_avg, "sa_coverage_summary.csv")

# Build acceptance history
sa_accept_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$accepted_history),
    accepted = as.integer(x$accepted_history)
  )
}))

# Average acceptance history
sa_accept_avg <- sa_accept_history %>%
  group_by(mode, iteration) %>%
  summarise(mean_acceptance = mean(accepted), .groups = "drop")

#write_csv(sa_accept_history, "sa_acceptance_history_all_runs.csv")
#write_csv(sa_accept_avg, "sa_acceptance_summary.csv")

# Build temperature history
sa_temp <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$temp_history),
    temperature = x$temp_history
  )
}))

# Average temperature history
sa_temp_avg <- sa_temp %>%
  group_by(iteration) %>%
  summarise(mean_temperature = mean(temperature), .groups = "drop")

#write_csv(sa_temp, "sa_temperature_history_all_runs.csv")
#write_csv(sa_temp_avg, "sa_temperature_summary.csv")

# Plot SA objective convergence
p_sa_convergence <- ggplot(sa_history_avg, aes(iteration, mean_best, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_best - ci_best, ymax = mean_best + ci_best), alpha = 0.20, colour = NA) +
  scale_colour_manual(values = mode_cols) +
  scale_fill_manual(values = mode_cols) +
  labs(
    title = "SA Objective Convergence",
    x = "Iteration",
    y = "Mean best penalised objective",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_sa_convergence)

# Plot SA coverage progress
p_sa_coverage <- ggplot(sa_coverage_avg, aes(iteration, mean_coverage, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_coverage - ci_coverage, ymax = mean_coverage + ci_coverage), alpha = 0.18, colour = NA) +
  geom_hline(yintercept = min_coverage, linetype = "dotted", linewidth = 1.1, colour = extra_cols[5]) +
  scale_colour_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  scale_fill_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  labs(
    title = "SA Coverage Progression",
    x = "Iteration",
    y = "Mean coverage",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_sa_coverage)

# Plot SA cooling schedule
p_sa_temp <- ggplot(sa_temp_avg, aes(iteration, mean_temperature)) +
  geom_line(linewidth = 1.2, colour = extra_cols[2]) +
  labs(
    title = "SA Cooling Schedule",
    x = "Iteration",
    y = "Temperature"
  ) +
  theme_minimal(base_size = 14)

print(p_sa_temp)

# Plot SA acceptance trend
p_sa_accept <- ggplot(sa_accept_avg, aes(iteration, mean_acceptance, colour = mode)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.12) +
  scale_colour_manual(values = c(deterministic = extra_cols[1], stochastic = extra_cols[3])) +
  labs(
    title = "SA Acceptance Trend",
    x = "Iteration",
    y = "Mean acceptance rate",
    colour = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_sa_accept)

# Print main results
cat("\n==============================\n")
cat("SA MAIN RESULTS\n")
cat("==============================\n\n")
cat("Best SA solutions by model:\n\n")

# Print best solution per model
for (i in seq_len(nrow(best_sa_by_mode))) {
  cat("Model:", best_sa_by_mode$mode[i], "\n")
  cat("Seed:", best_sa_by_mode$seed[i], "\n")
  cat("Mean distance:", round(best_sa_by_mode$mean_distance[i] / 1000, 2), "km\n")
  cat("Coverage:", sprintf("%.4f", best_sa_by_mode$coverage[i]), "\n")
  cat("Target coverage: 0.8000\n")
  cat("Meets 80% coverage:", best_sa_by_mode$meets_coverage[i], "\n")
  cat("Selected sites:", best_sa_by_mode$selected_sites[i], "\n\n")
}

cat("Summary across 30 runs:\n\n")

# Print clean summary table
sa_summary_clean <- sa_summary %>%
  transmute(
    Model = mode,
    `Mean distance (km)` = round(mean_distance_km, 2),
    `Std. dev. distance` = round(sd_distance_km, 4),
    `Mean coverage` = sprintf("%.4f", mean_coverage),
    `Best coverage` = sprintf("%.4f", best_coverage),
    `Coverage success rate` = round(coverage_success_rate, 3),
    `Mean runtime (s)` = round(mean_runtime, 2)
  )

print(as.data.frame(sa_summary_clean), row.names = FALSE)

cat("\nSA experiment complete.\n")
cat("Saved: sa_30_run_results.csv and sa_summary_table.csv. \n")