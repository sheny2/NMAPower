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
#'
#' @return An object of class "nma_power_sim" containing $summary and $iterations.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- power_sim_nma_3(S = 100, k_ab = 6, k_ac = 6, k_bc = 3)
#' print(res)
#' hist(res$iterations$point_est, main = "Distribution of Point Estimates", xlab = "Log Odds Ratio (B vs C)")
#' }
power_sim_nma_3 <- function(S = 100,
                          k_ab = 0, k_ac = 0, k_bc = 0,
                          pi_a = 0.5,
                          OR_ab = 1.4, OR_ac = 1.8,
                          tau = 0.01,
                          verbose = FALSE,
                          n.cores = parallel::detectCores() - 1) {

  if (!requireNamespace("gemtc", quietly = TRUE)) stop("Package 'gemtc' is required.")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required.")
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")

  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  # Run parallel simulations
  result_matrix <- foreach::foreach(i = 1:S, .combine = rbind, .errorhandling = 'remove', .packages = "gemtc") %dopar% {

    # 1. Simulate k_ab studies (indirect)
    dat_ab = data.frame(matrix(ncol = 4, nrow = 0))
    if (k_ab != 0){
      for (j in 1:k_ab) {
        n_j = runif(n = 1, min = 100, max = 500)
        n_aj = n_bj = round(n_j / 2)

        log_OR_ab_j = rnorm(n = 1, mean = log(OR_ab), sd = tau)
        pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
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
        pi_aj = runif(n = 1, min = 1/2 * pi_a, max = 3/2 * pi_a)
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
        pi_bj = runif(n = 1, min = 1/2 * pi_b, max = 3/2 * pi_b)
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

    # Rank order via SUCRA
    prob <- gemtc::rank.probability(cons.out, preferredDirection = 1)
    prob <- round(prob, digits = 3)
    sucra <- gemtc::sucra(prob)
    rank_order = rownames(as.matrix(sort(sucra, decreasing = TRUE)))

    # Extract relative effects based on network composition
    if (k_ab == 0 & k_ac == 0) {
      res = summary(gemtc::relative.effect(cons.out, "B", c("B", "C")))
      lower_95 = res$summaries$quantiles[2, 1]
      upper_95 = res$summaries$quantiles[2, 5]
      point_est = res$summaries$statistics[2, 1]
      true_rank = c("C", "B")
    } else {
      res = summary(gemtc::relative.effect(cons.out, "B", c("A", "B", "C")))
      lower_95 = res$summaries$quantiles[3, 1]
      upper_95 = res$summaries$quantiles[3, 5]
      point_est = res$summaries$statistics[3, 1]
      true_rank = c("C", "B", "A")
    }

    # Calculate performance metrics
    reject_null = as.numeric(0 > upper_95 || 0 < lower_95)
    rank_correct = as.numeric(identical(rank_order, true_rank))
    point_true = log(OR_bc)

    # Return metrics
    c(
      reject_null = reject_null,
      rank_correct = rank_correct,
      bias = point_est - point_true,
      abs_bias = abs(point_est - point_true),
      point_est = point_est,
      ci_lower = lower_95,
      ci_upper = upper_95,
      true_log_or = point_true
    )
  }

  # Handle case where all iterations failed
  if (is.null(result_matrix) || nrow(result_matrix) == 0) {
    warning("All simulation iterations failed.")
    return(NULL)
  }

  # Format the raw iterations into a dataframe
  df_iterations <- as.data.frame(result_matrix)

  # Calculate summary statistics across all successful iterations
  actual_S <- nrow(df_iterations)
  summary_stats <- list(
    successful_iterations = actual_S,
    power = mean(df_iterations$reject_null, na.rm = TRUE),
    rank_correct_prob = mean(df_iterations$rank_correct, na.rm = TRUE),
    avg_bias = mean(df_iterations$bias, na.rm = TRUE),
    avg_abs_bias = mean(df_iterations$abs_bias, na.rm = TRUE)
  )

  # Build final return object
  out <- list(
    summary = summary_stats,
    iterations = df_iterations,
    parameters = list(S = S, k_ab = k_ab, k_ac = k_ac, k_bc = k_bc, pi_a = pi_a, OR_ab = OR_ab, OR_ac = OR_ac, tau = tau)
  )

  class(out) <- "nma_power_sim"
  return(out)
}

#' Print method for nma_power_sim objects
#' @export
print.nma_power_sim <- function(x, ...) {
  cat("\n==============================================\n")
  cat("   NMA Power Simulation Results\n")
  cat("==============================================\n")
  cat(sprintf("Iterations Executed : %d / %d\n", x$summary$successful_iterations, x$parameters$S))
  cat(sprintf("Statistical Power   : %.2f%%\n", x$summary$power * 100))
  cat(sprintf("Correct Rank Prob.  : %.2f%%\n", x$summary$rank_correct_prob * 100))
  cat(sprintf("Average Bias        : %.4f\n", x$summary$avg_bias))
  cat(sprintf("Average Abs. Bias   : %.4f\n", x$summary$avg_abs_bias))
  cat("----------------------------------------------\n")
  cat("Network Configuration:\n")
  cat(sprintf("  k_ab (indirect) = %d\n", x$parameters$k_ab))
  cat(sprintf("  k_ac (indirect) = %d\n", x$parameters$k_ac))
  cat(sprintf("  k_bc (direct)   = %d\n", x$parameters$k_bc))
  cat("==============================================\n")
}
