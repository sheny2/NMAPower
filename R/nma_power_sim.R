#' Generalized Power Simulation for Arbitrary Node Network Meta-Analysis
#'
#' Evaluates statistical power, ranking accuracy, and bias for a Bayesian NMA
#' using a unified edge-list dataframe. Supports any network topology.
#'
#' @param S Number of simulation iterations.
#' @param target_contrast A character vector of length 2 defining the comparison
#'   of interest (e.g., c("B", "C")). Evaluates the effect of the 2nd vs the 1st.
#' @param network_design A data.frame defining all edges in the network. Must contain:
#'   \itemize{
#'     \item \code{t1}: Name of treatment 1 (e.g., "A").
#'     \item \code{t2}: Name of treatment 2 (e.g., "B").
#'     \item \code{k}: Number of studies on this direct comparison (can be 0).
#'     \item \code{OR}: True Odds Ratio of t2 vs t1. If NA, the function will attempt
#'           to infer it assuming statistical consistency from other connected edges.
#'   }
#' @param pi_base Baseline event probability for the reference node (the first
#'   treatment found in the network).
#' @param tau Between-study heterogeneity (standard deviation).
#' @param verbose Logical. If FALSE, suppresses JAGS/MCMC console output.
#' @param n.cores Number of cores to use for parallel processing.
#' @param seed Optional non-negative integer used to create reproducible,
#'   independent random-number streams on parallel workers.
#' @param convergence_threshold Maximum acceptable Gelman-Rubin PSRF. Set to
#'   `NULL` to skip convergence checks.
#'
#' @return An object of class "nma_power_sim" containing $summary and $iterations.
#' @export
#'
#' @examples
#' \dontrun{
#' # Define an arbitrary network: A is standard of care.
#' # B, C, and D are active treatments.
#' design <- data.frame(
#'   t1 = c("A", "A", "A", "B"),
#'   t2 = c("B", "C", "D", "C"),
#'   k  = c(6,   6,   2,   3),
#'   OR = c(1.2, 1.4, 1.1, NA) # OR for B vs C is inferred under consistency
#' )
#'
#' res <- nma_power_sim(
#'   S = 100,
#'   target_contrast = c("B", "C"),
#'   network_design = design
#' )
#' print(res)
#' }
#'
nma_power_sim <- function(S = 100,
                             target_contrast = c("B", "C"),
                             network_design,
                             pi_base = 0.5,
                             tau = 0.2,
                             verbose = FALSE,
                             n.cores = NULL,
                             seed = NULL,
                             convergence_threshold = 1.1) {

  # Ensure required packages
  if (!requireNamespace("gemtc", quietly = TRUE)) stop("Package 'gemtc' is required.")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required.")
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")

  S <- .assert_scalar_count(S, "S")
  n.cores <- .validate_n_cores(n.cores)
  seed <- .validate_seed(seed)
  convergence_threshold <- .validate_convergence_threshold(convergence_threshold)
  pi_base <- .assert_probability(pi_base, "pi_base")
  tau <- .assert_nonnegative(tau, "tau")

  if (!is.data.frame(network_design) || nrow(network_design) == 0L) {
    stop("`network_design` must be a non-empty data frame.", call. = FALSE)
  }
  required_columns <- c("t1", "t2", "k", "OR")
  if (!all(required_columns %in% names(network_design))) {
    stop("`network_design` must contain columns `t1`, `t2`, `k`, and `OR`.", call. = FALSE)
  }

  network_design$t1 <- as.character(network_design$t1)
  network_design$t2 <- as.character(network_design$t2)
  network_design$k <- as.numeric(network_design$k)
  network_design$OR <- as.numeric(network_design$OR)

  if (anyNA(network_design[, c("t1", "t2", "k")]) ||
      any(network_design$t1 == "" | network_design$t2 == "")) {
    stop("Treatment names and study counts cannot be missing or empty.", call. = FALSE)
  }
  if (any(network_design$t1 == network_design$t2)) {
    stop("Each network edge must compare two distinct treatments.", call. = FALSE)
  }
  if (any(!is.finite(network_design$k)) || any(network_design$k < 0) ||
      any(network_design$k != floor(network_design$k))) {
    stop("`network_design$k` must contain non-negative integers.", call. = FALSE)
  }
  if (any(!is.na(network_design$OR) &
          (!is.finite(network_design$OR) | network_design$OR <= 0))) {
    stop("Non-missing odds ratios must be finite and greater than zero.", call. = FALSE)
  }

  target_contrast <- as.character(target_contrast)
  if (length(target_contrast) != 2L || anyNA(target_contrast) ||
      target_contrast[1] == target_contrast[2]) {
    stop("`target_contrast` must identify two distinct treatments.", call. = FALSE)
  }

  # Ensure target contrast exists in the network
  all_nodes <- unique(c(network_design$t1, network_design$t2))
  if (!all(target_contrast %in% all_nodes)) {
    stop("Both treatments in `target_contrast` must exist in `network_design`.")
  }

  # -------------------------------------------------------------------------
  # Infer True Parameters for All Nodes (Consistency Check)
  # -------------------------------------------------------------------------
  ref_node <- all_nodes[1]
  true_log_OR <- setNames(rep(NA, length(all_nodes)), all_nodes)
  true_log_OR[ref_node] <- 0

  resolved <- c(ref_node)
  iter_count <- 0
  max_iters <- length(all_nodes) * 2

  while (length(resolved) < length(all_nodes) && iter_count < max_iters) {
    added <- FALSE
    for (i in 1:nrow(network_design)) {
      u <- network_design$t1[i]
      v <- network_design$t2[i]
      or <- network_design$OR[i]

      if (is.na(or)) next

      if (u %in% resolved && !(v %in% resolved)) {
        true_log_OR[v] <- true_log_OR[u] + log(or)
        resolved <- c(resolved, v)
        added <- TRUE
      } else if (v %in% resolved && !(u %in% resolved)) {
        true_log_OR[u] <- true_log_OR[v] - log(or)
        resolved <- c(resolved, u)
        added <- TRUE
      }
    }
    iter_count <- iter_count + 1
    if (!added && length(resolved) < length(all_nodes)) {
      stop("Cannot resolve network consistency. Please provide enough ORs to link all treatments.")
    }
  }

  supplied <- !is.na(network_design$OR)
  implied <- true_log_OR[network_design$t2[supplied]] - true_log_OR[network_design$t1[supplied]]
  if (any(abs(implied - log(network_design$OR[supplied])) > 1e-8)) {
    stop("Supplied odds ratios are inconsistent around at least one network cycle.", call. = FALSE)
  }

  # Calculate True Target OR and Expected Ranking
  point_true <- true_log_OR[target_contrast[2]] - true_log_OR[target_contrast[1]]

  # Filter to only nodes that actually have data generated (k > 0 anywhere)
  active_nodes <- unique(c(
    network_design$t1[network_design$k > 0],
    network_design$t2[network_design$k > 0]
  ))
  # Ensure target nodes are tracked in rank even if no direct data is generated
  active_nodes <- unique(c(active_nodes, target_contrast))
  observed_nodes <- unique(c(
    network_design$t1[network_design$k > 0],
    network_design$t2[network_design$k > 0]
  ))
  if (!all(target_contrast %in% observed_nodes)) {
    stop("Both target treatments must occur on an edge with at least one study.", call. = FALSE)
  }
  observed_edges <- network_design$k > 0
  if (!.network_is_connected(network_design$t1[observed_edges],
                             network_design$t2[observed_edges], observed_nodes)) {
    stop("Edges with `k > 0` must form one connected evidence network.", call. = FALSE)
  }
  true_rank_expected <- names(sort(true_log_OR[active_nodes], decreasing = TRUE))

  # -------------------------------------------------------------------------
  # Parallel Processing Setup
  # -------------------------------------------------------------------------
  cl <- .start_cluster(n.cores, seed)
  on.exit(parallel::stopCluster(cl))

  results <- foreach::foreach(iter = seq_len(S), .combine = .append_foreach_result,
                              .init = list(), .errorhandling = "pass",
                              .packages = c("gemtc", "coda"),
                              .export = ".check_mcmc_convergence") %dopar% {

    # Generate binomial data for a specific edge
    simulate_edge <- function(k, t1_name, t2_name, log_OR_true, base_prob, tau_val, start_id) {
      dat = data.frame(matrix(ncol = 4, nrow = 0))
      for (j in 1:k) {
        n_j = round(runif(1, 100, 500))
        n1 = round(n_j / 2)
        n2 = n_j - n1

        log_OR_j = rnorm(1, log_OR_true, tau_val)

        # Vary baseline prob around true base (U(0.5p, 1.5p)) protecting bounds
        min_pi = max(0.001, 0.5 * base_prob)
        max_pi = min(0.999, 1.5 * base_prob)
        pi_1j = runif(1, min_pi, max_pi)

        odds_1j = pi_1j / (1 - pi_1j)
        odds_2j = odds_1j * exp(log_OR_j)
        pi_2j = odds_2j / (1 + odds_2j)

        dat = rbind(dat,
                    data.frame(study = start_id + j - 1, treatment = t1_name, sampleSize = n1, responders = rbinom(1, n1, pi_1j)),
                    data.frame(study = start_id + j - 1, treatment = t2_name, sampleSize = n2, responders = rbinom(1, n2, pi_2j)))
      }
      return(dat)
    }

    # Generate all data
    all_data <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(all_data) <- c("study", "treatment", "sampleSize", "responders")
    study_idx <- 1

    for (i in 1:nrow(network_design)) {
      if (network_design$k[i] > 0) {
        u <- network_design$t1[i]
        v <- network_design$t2[i]

        # Calculate true baseline pi for treatment 'u' relative to network 'pi_base'
        odds_ref <- pi_base / (1 - pi_base)
        odds_u_true <- odds_ref * exp(true_log_OR[u])
        pi_u_true <- odds_u_true / (1 + odds_u_true)

        true_log_or_uv <- true_log_OR[v] - true_log_OR[u]

        dat_edge <- simulate_edge(network_design$k[i], u, v, true_log_or_uv, pi_u_true, tau, study_idx)
        all_data <- rbind(all_data, dat_edge)
        study_idx <- study_idx + network_design$k[i]
      }
    }

    # Format and fit MCMC model
    all_data[, c('sampleSize', 'responders')] <- lapply(all_data[, c('sampleSize', 'responders')], as.numeric)
    network <- gemtc::mtc.network(data.ab = all_data)
    cons.model <- gemtc::mtc.model(network, type = "consistency", likelihood = "binom", link = "logit", linearModel = "random")

    cons.out <- NULL
    if (!verbose) {
      capture.output(suppressWarnings(suppressMessages({
        cons.out <- gemtc::mtc.run(cons.model, n.adapt = 500, n.iter = 2000, thin = 1)
      })))
    } else {
      cons.out <- gemtc::mtc.run(cons.model, n.adapt = 500, n.iter = 2000, thin = 1)
    }
    .check_mcmc_convergence(cons.out, convergence_threshold)

    # SUCRA Ranking
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    sucra <- gemtc::sucra(round(prob, digits = 3))
    # Extract ranking strictly for active nodes
    rank_order <- rownames(as.matrix(sort(sucra, decreasing = TRUE)))
    expected_target_order <- if (point_true >= 0) target_contrast[2:1] else target_contrast
    observed_target_order <- rank_order[rank_order %in% target_contrast]

    # Target Relative Effect Extraction
    res <- summary(gemtc::relative.effect(cons.out, target_contrast[1], target_contrast[2]))
    s <- res$summaries
    if (is.matrix(s$quantiles)) {
      ci_lower <- s$quantiles[1, "2.5%"]
      ci_upper <- s$quantiles[1, "97.5%"]
      point_est <- s$statistics[1, "Mean"]
    } else {
      ci_lower <- s$quantiles["2.5%"]
      ci_upper <- s$quantiles["97.5%"]
      point_est <- s$statistics["Mean"]
    }

    # Evaluate
    c(
      reject_null = as.numeric(0 > ci_upper || 0 < ci_lower),
      rank_correct = as.numeric(identical(rank_order, true_rank_expected)),
      target_order_correct = as.numeric(identical(observed_target_order, expected_target_order)),
      bias = point_est - point_true,
      abs_bias = abs(point_est - point_true),
      point_est = point_est,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      true_log_or = point_true
    )
  }

  collected <- .collect_simulation_results(results, S)
  if (is.null(collected)) return(NULL)
  df_iterations <- collected$iterations
  summary_stats <- .simulation_summary(df_iterations, S, collected$failed_iterations)

  out <- list(
    summary = summary_stats,
    iterations = df_iterations,
    failure_messages = collected$failure_messages,
    parameters = list(S = S, target_contrast = target_contrast, network_design = network_design, tau = tau, true_log_ORs = true_log_OR)
  )

  class(out) <- "nma_power_sim"
  return(out)
}

#' Print method for nma_power_sim objects
#'
#' @param x An object returned by `nma_power_sim()`.
#' @param ... Additional arguments (currently ignored).
#' @export
print.nma_power_sim <- function(x, ...) {
  cat("\n==============================================\n")
  cat("NMA Power Simulation Results\n")
  cat("==============================================\n")
  cat(sprintf("Iterations Executed : %d / %d\n", x$summary$successful_iterations, x$parameters$S))
  cat(sprintf("Iterations Failed   : %d\n", x$summary$failed_iterations))
  cat(sprintf("Statistical Power   : %.2f%%\n", x$summary$power * 100))
  cat(sprintf("Correct Full Rank   : %.2f%%\n", x$summary$rank_correct_prob * 100))
  cat(sprintf("Correct Target Order: %.2f%%\n", x$summary$target_order_correct_prob * 100))
  cat(sprintf("Average Bias        : %.4f\n", x$summary$avg_bias))
  cat(sprintf("Average Abs. Bias   : %.4f\n", x$summary$avg_abs_bias))
  cat("----------------------------------------------\n")
  cat("Target Contrast Focus:\n")
  cat(sprintf("  %s vs %s (True OR = %.3f)\n",
              x$parameters$target_contrast[2],
              x$parameters$target_contrast[1],
              exp(x$parameters$true_log_ORs[x$parameters$target_contrast[2]] - x$parameters$true_log_ORs[x$parameters$target_contrast[1]])))
  cat("----------------------------------------------\n")
  cat("Network Configuration Input:\n")
  print(x$parameters$network_design[, setdiff(names(x$parameters$network_design), "OR"), drop = FALSE])
  cat("==============================================\n")
}
