# Common functions used by GA, SA and PSO

# Load package
library(dplyr)

# Load prepared input files
distance_matrix <- readRDS("distance_matrix.rds")
demand_points_det <- readRDS("demand_points_det.rds")
candidate_sites <- readRDS("candidate_sites.rds")
demand_centroids_det <- readRDS("demand_centroids_det.rds")
demand_weights_det <- readRDS("demand_weights_det.rds")
scenario_matrix <- readRDS("scenario_matrix.rds")
scenario_weights_list <- readRDS("scenario_weights_list.rds")
scenario_probs <- readRDS("scenario_probs.rds")

# Check input dimensions
stopifnot(nrow(distance_matrix) == length(demand_weights_det))
stopifnot(nrow(distance_matrix) == nrow(scenario_matrix))
stopifnot(ncol(distance_matrix) == nrow(candidate_sites))
stopifnot(abs(sum(scenario_probs) - 1) < 1e-8)

# Calculate demand denominators
deterministic_denom <- sum(demand_weights_det)
scenario_totals <- sapply(scenario_weights_list, sum)
stochastic_denom <- sum(scenario_probs * scenario_totals)

# Repair solution to exactly p sites
repair_solution <- function(solution, p, D = distance_matrix, w = demand_weights_det) {
  if (p > length(solution)) stop("p cannot be larger than the number of candidate sites.")
  solution <- as.integer(round(solution))
  solution[solution < 0] <- 0L
  solution[solution > 1] <- 1L
  
  open_sites <- which(solution == 1L)
  
  # Remove excess sites
  while (length(open_sites) > p) {
    selected_dist <- D[, open_sites, drop = FALSE]
    nearest_pos <- apply(selected_dist, 1, which.min)
    contribution <- tapply(w, nearest_pos, sum)
    contribution_full <- rep(0, length(open_sites))
    contribution_full[as.integer(names(contribution))] <- contribution
    drop_site <- open_sites[which.min(contribution_full)]
    solution[drop_site] <- 0L
    open_sites <- which(solution == 1L)
  }
  
  # Add missing sites
  while (length(open_sites) < p) {
    closed_sites <- which(solution == 0L)
    
    if (length(open_sites) == 0) {
      candidate_cost <- colSums(D[, closed_sites, drop = FALSE] * w)
      add_site <- closed_sites[which.min(candidate_cost)]
    } else {
      current_dist <- apply(D[, open_sites, drop = FALSE], 1, min)
      improvement <- sapply(closed_sites, function(j) sum(pmax(current_dist - D[, j], 0) * w))
      add_site <- closed_sites[which.max(improvement)]
    }
    
    solution[add_site] <- 1L
    open_sites <- which(solution == 1L)
  }
  
  solution
}

# Get nearest open-site distance
nearest_open_distance <- function(solution, D = distance_matrix, p) {
  solution <- repair_solution(solution, p, D = D)
  open_sites <- which(solution == 1L)
  if (length(open_sites) == 0) return(rep(Inf, nrow(D)))
  apply(D[, open_sites, drop = FALSE], 1, min)
}

# Calculate weighted distance
weighted_distance_from_nearest <- function(nearest_dist, mode = "deterministic", w = demand_weights_det, scenarios = scenario_weights_list, probs = scenario_probs) {
  if (mode == "deterministic") return(sum(nearest_dist * w))
  
  if (mode == "stochastic") {
    total <- 0
    for (s in seq_along(scenarios)) total <- total + probs[s] * sum(nearest_dist * scenarios[[s]])
    return(total)
  }
  
  stop("mode must be either 'deterministic' or 'stochastic'")
}

# Calculate coverage
coverage_from_nearest <- function(nearest_dist, threshold, mode = "deterministic", w = demand_weights_det, scenarios = scenario_weights_list, probs = scenario_probs) {
  covered <- nearest_dist <= threshold
  
  if (mode == "deterministic") return(sum(w[covered]) / sum(w))
  
  if (mode == "stochastic") {
    cov_value <- 0
    for (s in seq_along(scenarios)) {
      scenario_weight <- scenarios[[s]]
      scenario_total <- sum(scenario_weight)
      if (scenario_total > 0) cov_value <- cov_value + probs[s] * sum(scenario_weight[covered]) / scenario_total
    }
    return(cov_value)
  }
  
  stop("mode must be either 'deterministic' or 'stochastic'")
}

# Evaluate objective value
evaluate_solution <- function(solution, mode = "deterministic", D = distance_matrix, p) {
  nearest_dist <- nearest_open_distance(solution, D = D, p = p)
  weighted_distance_from_nearest(nearest_dist, mode = mode)
}

# Compute solution coverage
compute_coverage <- function(solution, threshold, mode = "deterministic", D = distance_matrix, p) {
  nearest_dist <- nearest_open_distance(solution, D = D, p = p)
  coverage_from_nearest(nearest_dist, threshold = threshold, mode = mode)
}

# Evaluate penalised objective
evaluate_penalised <- function(solution, mode = "deterministic", p, threshold = 5000, min_coverage = 0.80, penalty = 1e11, D = distance_matrix) {
  solution <- repair_solution(solution, p, D = D)
  open_sites <- which(solution == 1L)
  nearest_dist <- apply(D[, open_sites, drop = FALSE], 1, min)
  weighted_distance <- weighted_distance_from_nearest(nearest_dist, mode = mode)
  coverage <- coverage_from_nearest(nearest_dist, threshold = threshold, mode = mode)
  if (coverage >= min_coverage) return(weighted_distance)
  weighted_distance + penalty * (min_coverage - coverage)^2
}

# Return full solution metrics
evaluate_full <- function(solution, mode = "deterministic", p, threshold = 5000, min_coverage = 0.80, penalty = 1e11, D = distance_matrix) {
  solution <- repair_solution(solution, p, D = D)
  open_sites <- which(solution == 1L)
  nearest_dist <- apply(D[, open_sites, drop = FALSE], 1, min)
  weighted_distance <- weighted_distance_from_nearest(nearest_dist, mode = mode)
  coverage <- coverage_from_nearest(nearest_dist, threshold = threshold, mode = mode)
  penalised_objective <- if (coverage >= min_coverage) weighted_distance else weighted_distance + penalty * (min_coverage - coverage)^2
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

# Generate random p-site solution
generate_random_solution <- function(n_sites = ncol(distance_matrix), p) {
  if (p > n_sites) stop("p cannot be larger than n_sites.")
  solution <- rep(0L, n_sites)
  solution[sample(seq_len(n_sites), p)] <- 1L
  solution
}

# Generate neighbour solution
generate_neighbor <- function(solution, p) {
  solution <- repair_solution(solution, p)
  open_idx <- which(solution == 1L)
  closed_idx <- which(solution == 0L)
  if (length(open_idx) == 0 || length(closed_idx) == 0) return(solution)
  neighbor <- solution
  neighbor[sample(open_idx, 1)] <- 0L
  neighbor[sample(closed_idx, 1)] <- 1L
  repair_solution(neighbor, p)
}

# Assign demand to nearest site
get_assignments <- function(solution, D = distance_matrix, p) {
  solution <- repair_solution(solution, p, D = D)
  open_sites <- which(solution == 1L)
  selected_distances <- D[, open_sites, drop = FALSE]
  nearest_pos <- apply(selected_distances, 1, which.min)
  nearest_dist <- apply(selected_distances, 1, min)
  
  list(
    nearest_site_index = open_sites[nearest_pos],
    nearest_distance = nearest_dist
  )
}

# Select best solution by model
select_best_by_mode <- function(results_df) {
  results_df %>%
    group_by(mode) %>%
    arrange(desc(meets_coverage), mean_distance, desc(coverage), .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup()
}

# Calculate 95 percent confidence interval
ci95 <- function(x) {
  if (length(x) <= 1 || sd(x) == 0) return(0)
  1.96 * sd(x) / sqrt(length(x))
}

# Print input summary
cat("Shared helpers loaded.\n")
cat("Demand points:", nrow(distance_matrix), "\n")
cat("Candidate sites:", ncol(distance_matrix), "\n")
cat("Scenarios:", length(scenario_weights_list), "\n")