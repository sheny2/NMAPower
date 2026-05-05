#' Post-Hoc Power Simulation for Existing Network Meta-Analysis Data
#'
#' Fits a Bayesian NMA to an existing dataset to extract underlying parameters,
#' then performs a simulation-based power calculation preserving the original
#' network topology (sample sizes, study arms, and arm configurations).
#'
#' @param data A data.frame containing the existing NMA data. Expected columns:
#'   \code{study.id}, \code{treatment.id}, \code{sample.size}, \code{response}
#'   (matching the output of generate_NMA).
#' @param target_contrast A vector of length 2 specifying the treatments to compare
#'   (e.g., c(1, 2) or c("TrtA", "TrtB")).
#' @param S Number of simulation iterations.
#' @param verbose Logical. If FALSE, suppresses JAGS output during simulation.
#' @param n.cores Number of cores to use for parallel processing.
#'
#' @return An object of class "nma_posthoc_power" containing summary metrics.
#' @export
nma_power_posthoc <- function(data,
                                    target_contrast,
                                    S = 100,
                                    verbose = FALSE,
                                    n.cores = parallel::detectCores() - 1) {

  if (!requireNamespace("gemtc", quietly = TRUE)) stop("Package 'gemtc' is required.")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required.")
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")

  # Map user column names (e.g., from generate_NMA) to gemtc requirements
  expected_cols <- c("study.id", "treatment.id", "sample.size", "response")
  if (all(expected_cols %in% colnames(data))) {
    dat_gemtc <- data.frame(
      study = as.character(data$study.id),
      treatment = as.character(data$treatment.id),
      sampleSize = as.numeric(data$sample.size),
      responders = as.numeric(data$response)
    )
  } else {
    # Assume it's already in gemtc format
    dat_gemtc <- data
    dat_gemtc$study <- as.character(dat_gemtc$study)
    dat_gemtc$treatment <- as.character(dat_gemtc$treatment)
  }

  target_contrast <- as.character(target_contrast)
  if (!all(target_contrast %in% unique(dat_gemtc$treatment))) {
    stop("Target contrast treatments are not found in the dataset.")
  }

  # 2. Initial Model Fit (The Design Extraction Phase)
  message("Fitting initial model to extract design parameters... This may take a moment.")
  network_init <- gemtc::mtc.network(data.ab = dat_gemtc)
  model_init <- gemtc::mtc.model(network_init, type = "consistency",
                                 likelihood = "binom", link = "logit",
                                 linearModel = "random")

  out_init <- NULL
  capture.output(suppressWarnings(suppressMessages({
    out_init <- gemtc::mtc.run(model_init, n.adapt = 1000, n.iter = 5000, thin = 1)
  })))

  # Extract true relative effects relative to the network's reference treatment
  ref_trt <- network_init$treatments$id[1]
  rel_effects <- summary(gemtc::relative.effect(out_init, t1 = ref_trt))$summaries$statistics

  true_d <- setNames(rep(0, length(network_init$treatments$id)), network_init$treatments$id)
  for (trt in network_init$treatments$id) {
    if (trt != ref_trt) {
      # Handle naming conventions in gemtc relative effect output
      row_name <- paste0("d.", ref_trt, ".", trt)
      if (row_name %in% rownames(rel_effects)) {
        true_d[trt] <- rel_effects[row_name, "Mean"]
      } else {
        # Fallback if names are formatted differently
        true_d[trt] <- rel_effects[grep(trt, rownames(rel_effects))[1], "Mean"]
      }
    }
  }

  # Extract heterogeneity (tau)
  tau_true <- 0.01 # Default fallback
  stats_all <- summary(out_init)$statistics
  if ("sd.d" %in% rownames(stats_all)) tau_true <- stats_all["sd.d", "Mean"]

  # Target parameter for bias calculation
  true_target_logOR <- true_d[target_contrast[2]] - true_d[target_contrast[1]]
  true_rank_expected <- names(sort(true_d, decreasing = TRUE))

  # Calculate empirical study-specific baselines (mu_i) using a smoothed logit
  # to anchor the data generation process in the exact observed prevalence.
  study_baselines <- list()
  for (s in unique(dat_gemtc$study)) {
    s_data <- dat_gemtc[dat_gemtc$study == s, ]
    # Anchor to the first arm in the study (alphabetically/numerically)
    base_arm <- s_data[order(s_data$treatment), ][1, ]
    # Smoothed logit: log((y + 0.5) / (n - y + 0.5))
    y <- base_arm$responders
    n <- base_arm$sampleSize
    study_baselines[[s]] <- log((y + 0.5) / (n - y + 0.5))
  }

  message(sprintf("Parameters extracted. True Target log(OR) = %.3f, tau = %.3f",
                  true_target_logOR, tau_true))
  message("Starting parallel simulation loop...")

  # 3. Parallel Simulation Setup
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  # 4. Simulation Loop
  result_matrix <- foreach::foreach(iter = 1:S, .combine = rbind, .errorhandling = 'remove', .packages = "gemtc") %dopar% {

    # Generate new responses exactly mapping the existing network architecture
    sim_data <- dat_gemtc
    sim_data$responders <- NA

    for (s in unique(sim_data$study)) {
      s_data <- sim_data[sim_data$study == s, ]
      base_trt <- s_data[order(s_data$treatment), "treatment"][1]
      mu_s <- study_baselines[[s]]

      # Draw study-specific random effect vector matching the NMA hierarchical structure
      for (i in 1:nrow(s_data)) {
        trt_k <- s_data$treatment[i]
        n_k <- s_data$sampleSize[i]

        # True expected relative effect between arm k and the study's baseline arm
        delta_expected <- true_d[trt_k] - true_d[base_trt]
        # Incorporate between-study heterogeneity
        delta_jk <- rnorm(1, mean = delta_expected, sd = tau_true)

        # Calculate probability and generate binomial outcome
        logit_p <- mu_s + delta_jk
        p_jk <- exp(logit_p) / (1 + exp(logit_p))

        sim_data$responders[sim_data$study == s & sim_data$treatment == trt_k] <- rbinom(1, n_k, p_jk)
      }
    }

    # Fit the simulated dataset
    sim_network <- gemtc::mtc.network(data.ab = sim_data)
    sim_model <- gemtc::mtc.model(sim_network, type = "consistency",
                                  likelihood = "binom", link = "logit",
                                  linearModel = "random")

    sim_out <- NULL
    if (!verbose) {
      capture.output(suppressWarnings(suppressMessages({
        sim_out <- gemtc::mtc.run(sim_model, n.adapt = 500, n.iter = 2000, thin = 1)
      })))
    } else {
      sim_out <- gemtc::mtc.run(sim_model, n.adapt = 500, n.iter = 2000, thin = 1)
    }

    # Extract Metrics: Rank Probability
    prob <- gemtc::rank.probability(sim_out, preferredDirection = 1)
    sucra <- gemtc::sucra(round(prob, digits = 3))
    rank_order <- rownames(as.matrix(sort(sucra, decreasing = TRUE)))

    # Extract Metrics: Target Effect
    res <- summary(gemtc::relative.effect(sim_out, target_contrast[1], target_contrast[2]))
    s_stats <- res$summaries
    if (is.matrix(s_stats$quantiles)) {
      ci_lower <- s_stats$quantiles[1, "2.5%"]
      ci_upper <- s_stats$quantiles[1, "97.5%"]
      point_est <- s_stats$statistics[1, "Mean"]
    } else {
      ci_lower <- s_stats$quantiles["2.5%"]
      ci_upper <- s_stats$quantiles["97.5%"]
      point_est <- s_stats$statistics["Mean"]
    }

    # Return metrics array
    c(
      reject_null = as.numeric(0 > ci_upper || 0 < ci_lower),
      rank_correct = as.numeric(identical(rank_order, true_rank_expected)),
      bias = point_est - true_target_logOR,
      abs_bias = abs(point_est - true_target_logOR),
      point_est = point_est,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  }

  if (is.null(result_matrix) || nrow(result_matrix) == 0) {
    warning("All simulation iterations failed. Check initial model convergence.")
    return(NULL)
  }

  df_iterations <- as.data.frame(result_matrix)

  out <- list(
    summary = list(
      successful_iterations = nrow(df_iterations),
      power = mean(df_iterations$reject_null, na.rm = TRUE),
      rank_correct_prob = mean(df_iterations$rank_correct, na.rm = TRUE),
      avg_bias = mean(df_iterations$bias, na.rm = TRUE),
      avg_abs_bias = mean(df_iterations$abs_bias, na.rm = TRUE)
    ),
    iterations = df_iterations,
    parameters_used = list(
      target_contrast = target_contrast,
      true_target_logOR = true_target_logOR,
      tau = tau_true,
      true_d = true_d
    )
  )

  class(out) <- "nma_posthoc_power"
  return(out)
}

#' Print method for nma_posthoc_power objects
#' @export
print.nma_posthoc_power <- function(x, ...) {
  cat("\n==============================================\n")
  cat(" Post-Hoc NMA Power Evaluation Results\n")
  cat("==============================================\n")
  cat(sprintf("Iterations Executed : %d\n", x$summary$successful_iterations))
  cat(sprintf("Statistical Power   : %.2f%%\n", x$summary$power * 100))
  cat(sprintf("Correct Rank Prob.  : %.2f%%\n", x$summary$rank_correct_prob * 100))
  cat(sprintf("Average Bias        : %.4f\n", x$summary$avg_bias))
  cat(sprintf("Average Abs. Bias   : %.4f\n", x$summary$avg_abs_bias))
  cat("----------------------------------------------\n")
  cat("Extracted Design Parameters:\n")
  cat(sprintf("  Target Contrast: %s vs %s\n",
              x$parameters_used$target_contrast[2],
              x$parameters_used$target_contrast[1]))
  cat(sprintf("  True Target log(OR): %.3f\n", x$parameters_used$true_target_logOR))
  cat(sprintf("  Estimated Heterogeneity (tau): %.3f\n", x$parameters_used$tau))
  cat("==============================================\n")
}
