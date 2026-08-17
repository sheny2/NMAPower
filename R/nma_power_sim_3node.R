#' Power Simulation for 3-Node Network Meta-Analysis
#'
#' Evaluates statistical power, ranking accuracy, and bias for a Bayesian NMA.
#' Returns both summary statistics and raw iteration data for plotting.
#'
#' @section Note:
#' This function is specifically designed to mimic the exact simulation settings
#' described in the manuscript (e.g., drawing sample sizes from U(100, 500)
#' and using a specific binomial data-generation process for 3 nodes).
#'
#' @param S Number of simulation iterations.
#' @param k_ab Number of A vs B studies (indirect).
#' @param k_ac Number of A vs C studies (indirect).
#' @param k_bc Number of B vs C studies (direct).
#' @param pi_a Baseline event probability for treatment A.
#' @param OR_ab True odds ratio for A vs B.
#' @param OR_ac True odds ratio for A vs C.
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
#' res <- nma_power_sim_3node(S = 100, k_ab = 6, k_ac = 6, k_bc = 3)
#' print(res)
#' hist(res$iterations$point_est, main = "Distribution of Point Estimates", xlab = "Log Odds Ratio (B vs C)")
#' }
nma_power_sim_3node <- function(S = 100,
                          k_ab = 0, k_ac = 0, k_bc = 0,
                          pi_a = 0.5,
                          OR_ab = 1.4, OR_ac = 1.8,
                          tau = 0.01,
                          verbose = FALSE,
                          n.cores = NULL,
                          seed = NULL,
                          convergence_threshold = 1.1) {

  if (!requireNamespace("gemtc", quietly = TRUE)) stop("Package 'gemtc' is required.")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required.")
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")

  S <- .assert_scalar_count(S, "S")
  k_ab <- .assert_scalar_count(k_ab, "k_ab", minimum = 0L)
  k_ac <- .assert_scalar_count(k_ac, "k_ac", minimum = 0L)
  k_bc <- .assert_scalar_count(k_bc, "k_bc", minimum = 0L)
  pi_a <- .assert_probability(pi_a, "pi_a")
  tau <- .assert_nonnegative(tau, "tau")
  if (length(OR_ab) != 1L || is.na(OR_ab) || !is.numeric(OR_ab) || OR_ab <= 0 || !is.finite(OR_ab) ||
      length(OR_ac) != 1L || is.na(OR_ac) || !is.numeric(OR_ac) || OR_ac <= 0 || !is.finite(OR_ac)) {
    stop("`OR_ab` and `OR_ac` must be finite numbers greater than zero.", call. = FALSE)
  }
  n.cores <- .validate_n_cores(n.cores)
  seed <- .validate_seed(seed)
  convergence_threshold <- .validate_convergence_threshold(convergence_threshold)

  edge_counts <- c(k_ab, k_ac, k_bc)
  edge_t1 <- c("A", "A", "B")[edge_counts > 0]
  edge_t2 <- c("B", "C", "C")[edge_counts > 0]
  observed_nodes <- unique(c(edge_t1, edge_t2))
  if (!all(c("B", "C") %in% observed_nodes) ||
      !.network_is_connected(edge_t1, edge_t2, observed_nodes)) {
    stop("The positive study counts must form a connected network containing B and C.", call. = FALSE)
  }

  cl <- .start_cluster(n.cores, seed)
  on.exit(parallel::stopCluster(cl))

  # Run parallel simulations
  results <- foreach::foreach(i = seq_len(S), .combine = .append_foreach_result,
                              .init = list(), .errorhandling = "pass",
                              .packages = c("gemtc", "coda"),
                              .export = ".check_mcmc_convergence") %dopar% {

    # 1. Simulate k_ab studies (indirect)
    dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
    if (k_ab != 0){
      for (j in 1:k_ab) {
        n_j = runif(n = 1, min = 100, max = 500)
        n_aj = n_bj = round(n_j / 2)

        log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab), sd = tau)
        pi_aj = runif(n = 1, min = max(0.001, 1/2 * pi_a), max = min(0.999, 3/2 * pi_a))
        pi_bj = pi_aj * exp(log_OR_ab_j) / (1 - pi_aj + pi_aj * exp(log_OR_ab_j))

        e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
        e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)

        dat_ab = rbind(dat_ab, rbind(cbind(j, "A", n_aj, e_aj), cbind(j, "B", n_bj, e_bj)))
      }
    }

    # 2. Simulate k_ac studies (indirect)
    dat_ac = data.frame(matrix(ncol = 4, nrow = 0))
    if (k_ac != 0){
      for (j in 1:k_ac) {
        n_j = runif(n = 1, min = 100, max = 500)
        n_aj = n_cj = round(n_j / 2)

        log_OR_ac_j = rnorm(n = 1, mean = log(OR_ac), sd = tau)
        pi_aj = runif(n = 1, min = max(0.001, 1/2 * pi_a), max = min(0.999, 3/2 * pi_a))
        pi_cj = pi_aj * exp(log_OR_ac_j) / (1 - pi_aj + pi_aj * exp(log_OR_ac_j))

        e_aj = rbinom(n = 1, size = n_aj, prob = pi_aj)
        e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)

        dat_ac = rbind(dat_ac, rbind(cbind(j + k_ab, "A", n_aj, e_aj), cbind(j + k_ab, "C", n_cj, e_cj)))
      }
    }

    # 3. Simulate k_bc studies (direct)
    OR_bc = exp(log(OR_ac) - log(OR_ab))
    pi_b = (OR_ab * pi_a / (1 - pi_a)) / (1 + OR_ab * pi_a / (1 - pi_a))

    dat_bc = data.frame(matrix(ncol = 4, nrow = 0))
    if (k_bc != 0){
      for (j in 1:k_bc) {
        n_j = runif(n = 1, min = 100, max = 500)
        n_bj = n_cj = round(n_j / 2)

        log_OR_bc_j = rnorm(n = 1, mean = log(OR_bc), sd = tau)
        pi_bj = runif(n = 1, min = max(0.001, 1/2 * pi_b), max = min(0.999, 3/2 * pi_b))
        pi_cj = pi_bj * exp(log_OR_bc_j) / (1 - pi_bj + pi_bj * exp(log_OR_bc_j))

        e_bj = rbinom(n = 1, size = n_bj, prob = pi_bj)
        e_cj = rbinom(n = 1, size = n_cj, prob = pi_cj)

        dat_bc = rbind(dat_bc, rbind(cbind(j + k_ab + k_ac, "B", n_bj, e_bj), cbind(j + k_ab + k_ac, "C", n_cj, e_cj)))
      }
    }

    # Compile dataset and format column types
    colnames(dat_ab) = colnames(dat_ac) = colnames(dat_bc) = c("study", "treatment", "sampleSize", "responders")
    all_data = rbind(dat_ab, dat_ac, dat_bc)
    all_data[, c('sampleSize', 'responders')] <- lapply(all_data[, c('sampleSize', 'responders')], as.numeric)

    # Construct network and MTC model
    network_abc <- gemtc::mtc.network(data.ab = all_data)
    cons.model <- gemtc::mtc.model(network_abc, type = "consistency", likelihood = "binom", link = "logit", linearModel = "random")

    # Run MCMC, suppressing output if verbose = FALSE
    cons.out <- NULL
    if (!verbose) {
      capture.output(suppressWarnings(suppressMessages({
        cons.out <- gemtc::mtc.run(cons.model, n.adapt = 500, n.iter = 2000, thin = 1)
      })))
    } else {
      cons.out <- gemtc::mtc.run(cons.model, n.adapt = 500, n.iter = 2000, thin = 1)
    }
    .check_mcmc_convergence(cons.out, convergence_threshold)

    # Rank order via SUCRA
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits = 3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = TRUE)))
    true_effects <- c(A = 0, B = log(OR_ab), C = log(OR_ac))
    true_rank <- names(sort(true_effects[rownames(prob)], decreasing = TRUE))
    expected_target_order <- if (OR_ac >= OR_ab) c("C", "B") else c("B", "C")
    observed_target_order <- rank_order[rank_order %in% c("B", "C")]

    # Extract the B-to-C effect directly, independent of network composition.
    res <- summary(gemtc::relative.effect(cons.out, "B", "C"))
    s_stats <- res$summaries
    if (is.matrix(s_stats$quantiles)) {
      lower_95 <- s_stats$quantiles[1, "2.5%"]
      upper_95 <- s_stats$quantiles[1, "97.5%"]
      point_est <- s_stats$statistics[1, "Mean"]
    } else {
      lower_95 <- s_stats$quantiles["2.5%"]
      upper_95 <- s_stats$quantiles["97.5%"]
      point_est <- s_stats$statistics["Mean"]
    }

    # Calculate performance metrics
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    point_true = log(OR_bc)

    # Return metrics
    c(
      reject_null = reject_null,
      rank_correct = rank_correct,
      target_order_correct = as.numeric(identical(observed_target_order, expected_target_order)),
      bias = point_est - point_true,
      abs_bias = abs(point_est - point_true),
      point_est = point_est,
      ci_lower = lower_95,
      ci_upper = upper_95,
      true_log_or = point_true
    )
  }

  collected <- .collect_simulation_results(results, S)
  if (is.null(collected)) return(NULL)
  df_iterations <- collected$iterations
  summary_stats <- .simulation_summary(df_iterations, S, collected$failed_iterations)

  # Build final return object
  out <- list(
    summary = summary_stats,
    iterations = df_iterations,
    failure_messages = collected$failure_messages,
    parameters = list(S = S, k_ab = k_ab, k_ac = k_ac, k_bc = k_bc, pi_a = pi_a, OR_ab = OR_ab, OR_ac = OR_ac, tau = tau)
  )

  class(out) <- c("nma_power_sim_3node", "nma_power_sim")
  return(out)
}

#' Print method for nma_power_sim_3node objects
#'
#' @param x An object returned by `nma_power_sim_3node()`.
#' @param ... Additional arguments (currently ignored).
#' @export
print.nma_power_sim_3node <- function(x, ...) {
  cat("\n==============================================\n")
  cat("   NMA Power Simulation Results\n")
  cat("==============================================\n")
  cat(sprintf("Iterations Executed : %d / %d\n", x$summary$successful_iterations, x$parameters$S))
  cat(sprintf("Iterations Failed   : %d\n", x$summary$failed_iterations))
  cat(sprintf("Statistical Power   : %.2f%%\n", x$summary$power * 100))
  cat(sprintf("Correct Full Rank   : %.2f%%\n", x$summary$rank_correct_prob * 100))
  cat(sprintf("Correct Target Order: %.2f%%\n", x$summary$target_order_correct_prob * 100))
  cat(sprintf("Average Bias        : %.4f\n", x$summary$avg_bias))
  cat(sprintf("Average Abs. Bias   : %.4f\n", x$summary$avg_abs_bias))
  cat("----------------------------------------------\n")
  cat("Network Configuration:\n")
  cat(sprintf("  k_ab (indirect) = %d\n", x$parameters$k_ab))
  cat(sprintf("  k_ac (indirect) = %d\n", x$parameters$k_ac))
  cat(sprintf("  k_bc (direct)   = %d\n", x$parameters$k_bc))
  cat("==============================================\n")
}
