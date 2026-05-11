# 05_binary_pso.R
# Binary PSO for the ambulance p-median model

library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(scales)

source("02_Shared_helpers.R")

p <- 10
coverage_threshold <- 5000
min_coverage <- 0.80
near_feasible_threshold <- min_coverage - 0.005
coverage_penalty <- 1e11
seeds <- 1:30

swarm_size <- 24
n_iter <- 100
w_inertia <- 0.7
c1 <- 1.6
c2 <- 1.6
vmax <- 4

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

# Convert velocity to probability
sigmoid <- function(x) {
  x <- pmax(pmin(x, 30), -30)
  1 / (1 + exp(-x))
}

# Fast PSO repair with global-best preservation
fast_repair_pso <- function(solution, p, gbest = NULL) {
  solution <- as.integer(solution > 0)
  open_sites <- which(solution == 1L)
  
  # Drop excess sites
  if (length(open_sites) > p) {
    n_drop <- length(open_sites) - p
    
    if (!is.null(gbest)) {
      not_in_gbest <- open_sites[gbest[open_sites] == 0L]
      
      if (length(not_in_gbest) >= n_drop) {
        drop_sites <- sample(not_in_gbest, n_drop)
      } else {
        remaining_needed <- n_drop - length(not_in_gbest)
        remaining_pool <- setdiff(open_sites, not_in_gbest)
        drop_sites <- c(not_in_gbest, sample(remaining_pool, remaining_needed))
      }
    } else {
      drop_sites <- sample(open_sites, n_drop)
    }
    
    solution[drop_sites] <- 0L
  }
  
  # Add missing sites
  if (sum(solution) < p) {
    closed_sites <- which(solution == 0L)
    n_add <- p - sum(solution)
    
    if (!is.null(gbest)) {
      gbest_closed <- closed_sites[gbest[closed_sites] == 1L]
      
      if (length(gbest_closed) >= n_add) {
        add_sites <- sample(gbest_closed, n_add)
      } else {
        remaining_needed <- n_add - length(gbest_closed)
        remaining_pool <- setdiff(closed_sites, gbest_closed)
        add_sites <- c(gbest_closed, sample(remaining_pool, remaining_needed))
      }
    } else {
      add_sites <- sample(closed_sites, n_add)
    }
    
    solution[add_sites] <- 1L
  }
  
  solution
}

# Evaluate one PSO solution
evaluate_pso_solution <- function(solution, mode, gbest = NULL) {
  solution <- fast_repair_pso(solution, p, gbest)
  open_sites <- which(solution == 1L)
  nearest_dist <- apply(distance_matrix[, open_sites, drop = FALSE], 1, min)
  
  weighted_distance <- weighted_distance_from_nearest(nearest_dist, mode = mode)
  coverage <- coverage_from_nearest(nearest_dist, threshold = coverage_threshold, mode = mode)
  
  penalised_objective <- if (coverage >= min_coverage) {
    weighted_distance
  } else {
    weighted_distance + coverage_penalty * (min_coverage - coverage)^2
  }
  
  denom <- if (mode == "deterministic") deterministic_denom else stochastic_denom
  
  list(
    solution = solution,
    weighted_distance = weighted_distance,
    mean_distance = weighted_distance / denom,
    coverage = coverage,
    penalised_objective = penalised_objective,
    feasible = coverage >= min_coverage
  )
}

# Compare two solutions
is_better <- function(new, old) {
  if (new$feasible && !old$feasible) return(TRUE)
  if (!new$feasible && old$feasible) return(FALSE)
  if (new$feasible && old$feasible) return(new$weighted_distance < old$weighted_distance)
  if (abs(new$coverage - old$coverage) > 1e-10) return(new$coverage > old$coverage)
  new$penalised_objective < old$penalised_objective
}

# Create initial swarm
make_swarm <- function() {
  swarm <- matrix(0L, nrow = swarm_size, ncol = ncol(distance_matrix))
  for (i in seq_len(swarm_size)) {
    swarm[i, sample(seq_len(ncol(distance_matrix)), p)] <- 1L
  }
  swarm
}

# Run one PSO seed
run_pso <- function(mode, seed) {
  set.seed(seed)
  
  n_sites <- ncol(distance_matrix)
  
  # Initialise swarm and velocity
  swarm <- make_swarm()
  velocity <- matrix(
    runif(swarm_size * n_sites, -1, 1),
    nrow = swarm_size,
    ncol = n_sites
  )
  
  # Initialise personal bests
  pbest <- swarm
  pbest_metrics <- lapply(seq_len(swarm_size), function(i) {
    evaluate_pso_solution(swarm[i, ], mode)
  })
  
  # Initialise global best
  gbest_id <- which.min(sapply(pbest_metrics, function(x) x$penalised_objective))
  gbest_metrics <- pbest_metrics[[gbest_id]]
  gbest <- gbest_metrics$solution
  
  # Store search history
  best_history <- numeric(n_iter)
  mean_history <- numeric(n_iter)
  coverage_history <- numeric(n_iter)
  swarm_coverage_history <- numeric(n_iter)
  
  # Main PSO loop
  for (iter in seq_len(n_iter)) {
    current_objectives <- numeric(swarm_size)
    current_coverages <- numeric(swarm_size)
    
    # Update each particle
    for (i in seq_len(swarm_size)) {
      r1 <- runif(n_sites)
      r2 <- runif(n_sites)
      
      # Update velocity
      velocity[i, ] <- w_inertia * velocity[i, ] +
        c1 * r1 * (pbest[i, ] - swarm[i, ]) +
        c2 * r2 * (gbest - swarm[i, ])
      
      # Limit velocity
      velocity[i, ] <- pmax(pmin(velocity[i, ], vmax), -vmax)
      
      # Create binary position
      move_probability <- sigmoid(velocity[i, ])
      new_position <- as.integer(runif(n_sites) < move_probability)
      
      # Repair to exactly p sites
      new_position <- fast_repair_pso(new_position, p, gbest)
      
      # Evaluate particle
      this_metrics <- evaluate_pso_solution(new_position, mode, gbest)
      
      # Try local swap every 10 iterations
      if (iter %% 10 == 0) {
        trial <- this_metrics$solution
        open_sites <- which(trial == 1L)
        closed_sites <- which(trial == 0L)
        
        if (length(open_sites) > 0 && length(closed_sites) > 0) {
          trial[sample(open_sites, 1)] <- 0L
          trial[sample(closed_sites, 1)] <- 1L
          trial <- fast_repair_pso(trial, p, gbest)
          
          trial_metrics <- evaluate_pso_solution(trial, mode, gbest)
          
          if (is_better(trial_metrics, this_metrics)) {
            this_metrics <- trial_metrics
          }
        }
      }
      
      # Store particle state
      swarm[i, ] <- this_metrics$solution
      current_objectives[i] <- this_metrics$penalised_objective
      current_coverages[i] <- this_metrics$coverage
      
      # Update personal best
      if (is_better(this_metrics, pbest_metrics[[i]])) {
        pbest[i, ] <- this_metrics$solution
        pbest_metrics[[i]] <- this_metrics
      }
      
      # Update global best
      if (is_better(this_metrics, gbest_metrics)) {
        gbest <- this_metrics$solution
        gbest_metrics <- this_metrics
      }
    }
    
    # Store iteration metrics
    best_history[iter] <- gbest_metrics$penalised_objective
    mean_history[iter] <- mean(current_objectives)
    coverage_history[iter] <- gbest_metrics$coverage
    swarm_coverage_history[iter] <- mean(current_coverages)
  }
  
  # Return run result
  list(
    method = "PSO",
    mode = mode,
    seed = seed,
    best_solution = gbest_metrics$solution,
    selected_sites = which(gbest_metrics$solution == 1L),
    weighted_distance = gbest_metrics$weighted_distance,
    mean_distance = gbest_metrics$mean_distance,
    penalised_objective = gbest_metrics$penalised_objective,
    coverage = gbest_metrics$coverage,
    best_history = best_history,
    mean_history = mean_history,
    coverage_history = coverage_history,
    swarm_coverage_history = swarm_coverage_history
  )
}

# Prepare storage
runs <- list()
rows <- list()
k <- 1

cat("\nRunning PSO experiments...\n\n")

# Run deterministic and stochastic experiments
for (mode in c("deterministic", "stochastic")) {
  for (seed in seeds) {
    
    # Time each run
    start_time <- Sys.time()
    res <- run_pso(mode, seed)
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
      method = "PSO",
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
pso_results <- bind_rows(rows)

# Select best strict solution
best_pso_by_mode <- select_best_strict(pso_results)

# Save main outputs
saveRDS(runs, "pso_30_run_objects.rds")
write_csv(pso_results, "pso_30_run_results.csv")
#write_csv(best_pso_by_mode, "pso_best_solution_by_mode.csv")

# Summarise runs
pso_summary <- pso_results %>%
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
    near_feasible_rate = mean(coverage >= near_feasible_threshold),
    mean_runtime = mean(runtime_seconds),
    sd_runtime = sd(runtime_seconds),
    .groups = "drop"
  )

# Save summary
write_csv(pso_summary, "pso_summary_table.csv")

# Save report-friendly summary
pso_summary_display <- pso_summary %>%
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
    near_feasible_rate,
    mean_runtime
  )

#write_csv(pso_summary_display, "pso_summary_display.csv")

# Build convergence history
pso_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$best_history),
    best = x$best_history,
    mean = x$mean_history
  )
}))

# Average convergence history
pso_history_avg <- pso_history %>%
  group_by(mode, iteration) %>%
  summarise(
    mean_best = mean(best),
    sd_best = sd(best),
    ci_best = 1.96 * sd_best / sqrt(n()),
    .groups = "drop"
  )

#write_csv(pso_history, "pso_convergence_history_all_runs.csv")
#write_csv(pso_history_avg, "pso_convergence_summary.csv")

# Build coverage history
pso_coverage_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$coverage_history),
    coverage = x$coverage_history
  )
}))

# Average coverage history
pso_coverage_avg <- pso_coverage_history %>%
  group_by(mode, iteration) %>%
  summarise(
    mean_coverage = mean(coverage),
    sd_coverage = sd(coverage),
    ci_coverage = 1.96 * sd_coverage / sqrt(n()),
    .groups = "drop"
  )

#write_csv(pso_coverage_history, "pso_coverage_history_all_runs.csv")
#write_csv(pso_coverage_avg, "pso_coverage_summary.csv")

# Build swarm coverage history
pso_swarm_coverage_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    iteration = seq_along(x$swarm_coverage_history),
    swarm_coverage = x$swarm_coverage_history
  )
}))

# Average swarm coverage history
pso_swarm_coverage_avg <- pso_swarm_coverage_history %>%
  group_by(mode, iteration) %>%
  summarise(mean_swarm_coverage = mean(swarm_coverage), .groups = "drop")

#write_csv(pso_swarm_coverage_history, "pso_swarm_coverage_history_all_runs.csv")
#write_csv(pso_swarm_coverage_avg, "pso_swarm_coverage_summary.csv")

# Plot PSO objective convergence
p_pso_gbest <- ggplot(pso_history_avg, aes(iteration, mean_best, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_best - ci_best, ymax = mean_best + ci_best), alpha = 0.20, colour = NA) +
  scale_colour_manual(values = mode_cols) +
  scale_fill_manual(values = mode_cols) +
  labs(
    title = "PSO Global-Best Convergence",
    x = "Iteration",
    y = "Mean global-best penalised objective",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_pso_gbest)

# Plot PSO global-best coverage
p_pso_coverage <- ggplot(pso_coverage_avg, aes(iteration, mean_coverage, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_coverage - ci_coverage, ymax = mean_coverage + ci_coverage), alpha = 0.18, colour = NA) +
  geom_hline(yintercept = min_coverage, linetype = "dotted", linewidth = 1.1, colour = extra_cols[5]) +
  scale_colour_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  scale_fill_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  labs(
    title = "PSO Coverage Progression",
    x = "Iteration",
    y = "Mean global-best coverage",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_pso_coverage)

# Plot PSO swarm coverage
p_pso_swarm_cov <- ggplot(pso_swarm_coverage_avg, aes(iteration, mean_swarm_coverage, colour = mode)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = min_coverage, linetype = "dotted", linewidth = 1.1, colour = extra_cols[5]) +
  scale_colour_manual(values = c(deterministic = extra_cols[1], stochastic = extra_cols[2])) +
  labs(
    title = "PSO Mean Swarm Coverage",
    x = "Iteration",
    y = "Mean swarm coverage",
    colour = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_pso_swarm_cov)

# Print main results
cat("\n==============================\n")
cat("PSO MAIN RESULTS\n")
cat("==============================\n\n")
cat("Best PSO solutions by model:\n\n")

# Print best solution per model
for (i in seq_len(nrow(best_pso_by_mode))) {
  cat("Model:", best_pso_by_mode$mode[i], "\n")
  cat("Seed:", best_pso_by_mode$seed[i], "\n")
  cat("Mean distance:", round(best_pso_by_mode$mean_distance[i] / 1000, 2), "km\n")
  cat("Coverage:", sprintf("%.4f", best_pso_by_mode$coverage[i]), "\n")
  cat("Target coverage: 0.8000\n")
  cat("Meets 80% coverage:", best_pso_by_mode$meets_coverage[i], "\n")
  cat("Selected sites:", best_pso_by_mode$selected_sites[i], "\n\n")
}

cat("Summary across 30 runs:\n\n")

# Print clean summary table
pso_summary_clean <- pso_summary %>%
  transmute(
    Model = mode,
    `Mean distance (km)` = round(mean_distance_km, 2),
    `Std. dev. distance` = round(sd_distance_km, 4),
    `Mean coverage` = sprintf("%.4f", mean_coverage),
    `Best coverage` = sprintf("%.4f", best_coverage),
    `Coverage success rate` = round(coverage_success_rate, 3),
    `Near-feasible rate (>=0.795)` = round(near_feasible_rate, 3),
    `Mean runtime (s)` = round(mean_runtime, 2)
  )

print(as.data.frame(pso_summary_clean), row.names = FALSE)

cat(sprintf("\nNote: PSO uses 30 deterministic and 30 stochastic runs, with swarm size %d and %d iterations.\n", swarm_size, n_iter))


cat("\nPSO experiment complete.\n")
cat("Saved: pso_30_run_results.csv and pso_summary_table.csv \n")