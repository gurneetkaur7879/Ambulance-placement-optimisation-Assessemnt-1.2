# Genetic Algorithm for the ambulance p-median model

# Load packages
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(scales)

# Load shared helper functions and data
source("02_Shared_helpers.R")

# Model settings
p <- 10
coverage_threshold <- 5000
min_coverage <- 0.80
near_feasible_threshold <- min_coverage - 0.005
coverage_penalty <- 1e11
seeds <- 1:30

# GA settings
pop_size <- 60
generations <- 140
crossover_rate <- 0.8
mutation_rate <- 0.08
elite_size <- 1
tournament_k <- 3

# Plot colours
mode_cols <- c(deterministic = "#0072B2", stochastic = "#E69F00")
extra_cols <- c("#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#D55E00")

# Check model size
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

# Create random initial population
make_population <- function() {
  pop <- matrix(0L, nrow = pop_size, ncol = ncol(distance_matrix))
  for (i in seq_len(pop_size)) {
    pop[i, sample(seq_len(ncol(distance_matrix)), p)] <- 1L
  }
  pop
}

# Select one parent by tournament
select_parent <- function(pop, fit) {
  idx <- sample(seq_len(nrow(pop)), tournament_k)
  pop[idx[which.min(fit[idx])], ]
}

# Create two children from two parents
make_children <- function(parent1, parent2) {
  if (runif(1) > crossover_rate) {
    return(list(repair_solution(parent1, p), repair_solution(parent2, p)))
  }
  
  mask <- sample(c(TRUE, FALSE), length(parent1), replace = TRUE)
  
  list(
    repair_solution(ifelse(mask, parent1, parent2), p),
    repair_solution(ifelse(mask, parent2, parent1), p)
  )
}

# Mutate by swapping one open and one closed site
mutate_solution <- function(x) {
  x <- repair_solution(x, p)
  
  if (runif(1) < mutation_rate) {
    open_sites <- which(x == 1L)
    closed_sites <- which(x == 0L)
    x[sample(open_sites, 1)] <- 0L
    x[sample(closed_sites, 1)] <- 1L
  }
  
  repair_solution(x, p)
}

# Run one GA seed
run_ga <- function(mode, seed) {
  set.seed(seed)
  
  # Build starting population
  pop <- make_population()
  
  # Evaluate initial population
  fit <- apply(
    pop,
    1,
    evaluate_penalised,
    mode = mode,
    p = p,
    threshold = coverage_threshold,
    min_coverage = min_coverage,
    penalty = coverage_penalty
  )
  
  # Store search history
  best_history <- numeric(generations)
  mean_history <- numeric(generations)
  coverage_history <- numeric(generations)
  
  # Main GA loop
  for (gen in seq_len(generations)) {
    
    # Keep elite solution
    elite_order <- order(fit)[seq_len(elite_size)]
    new_pop <- pop[elite_order, , drop = FALSE]
    new_fit <- fit[elite_order]
    
    # Fill new generation
    while (nrow(new_pop) < pop_size) {
      parent1 <- select_parent(pop, fit)
      parent2 <- select_parent(pop, fit)
      children <- make_children(parent1, parent2)
      
      # Add first child
      child1 <- mutate_solution(children[[1]])
      new_pop <- rbind(new_pop, child1)
      new_fit <- c(
        new_fit,
        evaluate_penalised(
          child1,
          mode = mode,
          p = p,
          threshold = coverage_threshold,
          min_coverage = min_coverage,
          penalty = coverage_penalty
        )
      )
      
      # Add second child if needed
      if (nrow(new_pop) < pop_size) {
        child2 <- mutate_solution(children[[2]])
        new_pop <- rbind(new_pop, child2)
        new_fit <- c(
          new_fit,
          evaluate_penalised(
            child2,
            mode = mode,
            p = p,
            threshold = coverage_threshold,
            min_coverage = min_coverage,
            penalty = coverage_penalty
          )
        )
      }
    }
    
    # Update population
    pop <- new_pop[seq_len(pop_size), , drop = FALSE]
    fit <- new_fit[seq_len(pop_size)]
    
    # Save generation metrics
    best_id_gen <- which.min(fit)
    best_history[gen] <- min(fit)
    mean_history[gen] <- mean(fit)
    coverage_history[gen] <- compute_coverage(
      pop[best_id_gen, ],
      coverage_threshold,
      mode = mode,
      p = p
    )
  }
  
  # Evaluate final best solution
  best_id <- which.min(fit)
  best_metrics <- evaluate_full(
    pop[best_id, ],
    mode = mode,
    p = p,
    threshold = coverage_threshold,
    min_coverage = min_coverage,
    penalty = coverage_penalty
  )
  
  # Return run result
  list(
    method = "GA",
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
    coverage_history = coverage_history
  )
}

# Prepare storage
runs <- list()
rows <- list()
k <- 1

cat("\nRunning GA experiments...\n\n")

# Run deterministic and stochastic experiments
for (mode in c("deterministic", "stochastic")) {
  for (seed in seeds) {
    
    # Time each run
    start_time <- Sys.time()
    res <- run_ga(mode, seed)
    runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Check strict coverage
    meets_coverage <- res$coverage >= min_coverage
    coverage_gap <- max(0, min_coverage - res$coverage)
    
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
      method = "GA",
      mode = mode,
      seed = seed,
      weighted_distance = res$weighted_distance,
      mean_distance = res$mean_distance,
      penalised_objective = res$penalised_objective,
      coverage = res$coverage,
      coverage_gap = coverage_gap,
      meets_coverage = meets_coverage,
      runtime_seconds = runtime,
      selected_sites = paste(res$selected_sites, collapse = ", ")
    )
    
    k <- k + 1
  }
}

# Combine results
ga_results <- bind_rows(rows)

# Select best strict solution
best_ga_by_mode <- select_best_strict(ga_results)

# Save main outputs
saveRDS(runs, "ga_30_run_objects.rds")
write_csv(ga_results, "ga_30_run_results.csv")

# Summarise runs
ga_summary <- ga_results %>%
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
    mean_coverage_gap = mean(coverage_gap),
    best_coverage = max(coverage),
    coverage_success_rate = mean(meets_coverage),
    near_feasible_rate = mean(coverage >= near_feasible_threshold),
    mean_runtime = mean(runtime_seconds),
    sd_runtime = sd(runtime_seconds),
    .groups = "drop"
  )

# Save summary
write_csv(ga_summary, "ga_summary_table.csv")

# Save summary
ga_summary_display <- ga_summary %>%
  select(
    method,
    mode,
    mean_distance_km,
    sd_distance_km,
    best_weighted_distance,
    worst_weighted_distance,
    mean_coverage,
    mean_coverage_gap,
    best_coverage,
    coverage_success_rate,
    near_feasible_rate,
    mean_runtime
  )


# Build convergence history
ga_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    generation = seq_along(x$best_history),
    best = x$best_history,
    mean = x$mean_history
  )
}))

# Average convergence history
ga_history_avg <- ga_history %>%
  group_by(mode, generation) %>%
  summarise(
    mean_best = mean(best),
    sd_best = sd(best),
    ci_best = 1.96 * sd_best / sqrt(n()),
    .groups = "drop"
  )


write_csv(ga_history_avg, "ga_convergence_summary.csv")

# Build coverage history
ga_coverage_history <- bind_rows(lapply(runs, function(x) {
  data.frame(
    mode = x$mode,
    seed = x$seed,
    generation = seq_along(x$coverage_history),
    coverage = x$coverage_history
  )
}))

# Average coverage history
ga_coverage_avg <- ga_coverage_history %>%
  group_by(mode, generation) %>%
  summarise(
    mean_coverage = mean(coverage),
    sd_coverage = sd(coverage),
    ci_coverage = 1.96 * sd_coverage / sqrt(n()),
    .groups = "drop"
  )


write_csv(ga_coverage_avg, "ga_coverage_summary.csv")

# Plot GA objective convergence
p_ga_convergence <- ggplot(ga_history_avg, aes(generation, mean_best, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_best - ci_best, ymax = mean_best + ci_best), alpha = 0.20, colour = NA) +
  scale_colour_manual(values = mode_cols) +
  scale_fill_manual(values = mode_cols) +
  labs(
    title = "GA Convergence Behaviour",
    x = "Generation",
    y = "Mean best penalised objective",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_ga_convergence)

# Plot GA improvement rate
p_ga_improve <- ga_history_avg %>%
  group_by(mode) %>%
  mutate(improvement = abs(mean_best - lag(mean_best))) %>%
  filter(!is.na(improvement)) %>%
  ggplot(aes(generation, improvement, colour = mode)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c(deterministic = extra_cols[1], stochastic = extra_cols[2])) +
  labs(
    title = "GA Improvement Rate Over Generations",
    x = "Generation",
    y = "Absolute improvement in objective",
    colour = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_ga_improve)

# Plot GA coverage progress
p_ga_coverage <- ggplot(ga_coverage_avg, aes(generation, mean_coverage, colour = mode, fill = mode)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_coverage - ci_coverage, ymax = mean_coverage + ci_coverage), alpha = 0.18, colour = NA) +
  geom_hline(yintercept = min_coverage, linetype = "dotted", linewidth = 1.1, colour = extra_cols[5]) +
  scale_colour_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  scale_fill_manual(values = c(deterministic = extra_cols[4], stochastic = extra_cols[5])) +
  labs(
    title = "GA Coverage Progression",
    x = "Generation",
    y = "Mean coverage",
    colour = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14)

print(p_ga_coverage)

# Print main results
cat("\n==============================\n")
cat("GA MAIN RESULTS\n")
cat("==============================\n\n")
cat("Best GA solutions by model:\n\n")

# Print best solution per model
for (i in seq_len(nrow(best_ga_by_mode))) {
  cat("Model:", best_ga_by_mode$mode[i], "\n")
  cat("Seed:", best_ga_by_mode$seed[i], "\n")
  cat("Mean distance:", round(best_ga_by_mode$mean_distance[i] / 1000, 2), "km\n")
  cat("Coverage:", sprintf("%.4f", best_ga_by_mode$coverage[i]), "\n")
  cat("Target coverage: 0.8000\n")
  cat("Meets 80% coverage:", best_ga_by_mode$meets_coverage[i], "\n")
  cat("Selected sites:", best_ga_by_mode$selected_sites[i], "\n\n")
}

cat("Summary across 30 runs:\n\n")

# Print clean summary table
ga_summary_clean <- ga_summary %>%
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

print(as.data.frame(ga_summary_clean), row.names = FALSE)

cat("\nGA experiment complete.\n")
cat("Saved: ga_30_run_results.csv and ga_summary_table.csv.\n")


