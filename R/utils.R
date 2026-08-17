#' NMAPower package
#'
#' Simulation-based power analysis for Bayesian network meta-analysis.
#'
#' @keywords internal
#' @importFrom foreach %dopar%
#' @importFrom stats rbinom rnorm runif setNames
#' @importFrom utils capture.output
"_PACKAGE"

.assert_scalar_count <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x != floor(x) || x < minimum) {
    stop(sprintf("`%s` must be a single integer greater than or equal to %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(x)
}

.assert_probability <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x <= 0 || x >= 1) {
    stop(sprintf("`%s` must be a single number strictly between 0 and 1.", name),
         call. = FALSE)
  }
  as.numeric(x)
}

.assert_nonnegative <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x < 0) {
    stop(sprintf("`%s` must be a single non-negative number.", name), call. = FALSE)
  }
  as.numeric(x)
}

.validate_n_cores <- function(n.cores) {
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected)) detected <- 1L

  if (missing(n.cores) || is.null(n.cores)) {
    return(max(1L, detected - 1L))
  }

  n.cores <- .assert_scalar_count(n.cores, "n.cores")
  min(n.cores, detected)
}

.validate_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("`seed` must be NULL or a single non-negative integer.", call. = FALSE)
  }
  as.integer(seed)
}

.validate_convergence_threshold <- function(threshold) {
  if (is.null(threshold)) return(NULL)
  if (length(threshold) != 1L || is.na(threshold) || !is.numeric(threshold) || threshold <= 1) {
    stop("`convergence_threshold` must be NULL or a number greater than 1.", call. = FALSE)
  }
  as.numeric(threshold)
}

.network_is_connected <- function(t1, t2, nodes) {
  if (length(nodes) < 2L) return(TRUE)
  reached <- nodes[1]
  repeat {
    adjacent <- unique(c(t2[t1 %in% reached], t1[t2 %in% reached]))
    updated <- unique(c(reached, adjacent))
    if (length(updated) == length(reached)) break
    reached <- updated
  }
  all(nodes %in% reached)
}

.start_cluster <- function(n.cores, seed) {
  cl <- parallel::makeCluster(n.cores)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  doParallel::registerDoParallel(cl)
  cl
}

.append_foreach_result <- function(accumulator, value) {
  accumulator[[length(accumulator) + 1L]] <- value
  accumulator
}

.check_mcmc_convergence <- function(fit, threshold) {
  if (is.null(threshold)) return(invisible(TRUE))

  psrf <- tryCatch(
    coda::gelman.diag(fit, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."],
    error = function(e) {
      stop("Unable to calculate MCMC convergence diagnostics: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (any(!is.finite(psrf)) || max(psrf) > threshold) {
    stop(sprintf("MCMC convergence check failed (maximum PSRF %.3f; threshold %.3f).",
                 max(psrf), threshold), call. = FALSE)
  }
  invisible(TRUE)
}

.collect_simulation_results <- function(results, requested_iterations) {
  ok <- vapply(results, function(x) is.numeric(x) && !is.null(names(x)), logical(1))
  failures <- results[!ok]

  if (!any(ok)) {
    messages <- unique(vapply(failures, function(x) {
      if (inherits(x, "condition")) conditionMessage(x) else "Unknown worker failure"
    }, character(1)))
    warning(sprintf("All %d simulation iterations failed. %s",
                    requested_iterations, paste(messages, collapse = " | ")),
            call. = FALSE)
    return(NULL)
  }

  iteration_matrix <- do.call(rbind, results[ok])
  failed <- sum(!ok)
  if (failed > 0L) {
    messages <- unique(vapply(failures, function(x) {
      if (inherits(x, "condition")) conditionMessage(x) else "Unknown worker failure"
    }, character(1)))
    warning(sprintf("%d of %d simulation iterations failed and were excluded. %s",
                    failed, requested_iterations, paste(messages, collapse = " | ")),
            call. = FALSE)
  }

  list(
    iterations = as.data.frame(iteration_matrix),
    failed_iterations = failed,
    failure_messages = unique(vapply(failures, function(x) {
      if (inherits(x, "condition")) conditionMessage(x) else "Unknown worker failure"
    }, character(1)))
  )
}

.simulation_summary <- function(iterations, requested, failed) {
  list(
    requested_iterations = requested,
    successful_iterations = nrow(iterations),
    failed_iterations = failed,
    failure_rate = failed / requested,
    power = mean(iterations$reject_null, na.rm = TRUE),
    rank_correct_prob = mean(iterations$rank_correct, na.rm = TRUE),
    target_order_correct_prob = if ("target_order_correct" %in% names(iterations)) {
      mean(iterations$target_order_correct, na.rm = TRUE)
    } else {
      NA_real_
    },
    avg_bias = mean(iterations$bias, na.rm = TRUE),
    avg_abs_bias = mean(iterations$abs_bias, na.rm = TRUE)
  )
}
