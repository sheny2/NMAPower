test_that("simulation result collection reports failures", {
  results <- list(
    c(reject_null = 1, rank_correct = 1, target_order_correct = 1,
      bias = 0.1, abs_bias = 0.1),
    simpleError("worker failed")
  )

  expect_warning(
    collected <- NMAPower:::.collect_simulation_results(results, 2),
    "1 of 2"
  )
  summary <- NMAPower:::.simulation_summary(
    collected$iterations, 2, collected$failed_iterations
  )
  expect_equal(summary$successful_iterations, 1)
  expect_equal(summary$failed_iterations, 1)
  expect_equal(summary$failure_rate, 0.5)
})

test_that("general and three-node objects have distinct print methods", {
  summary <- list(
    successful_iterations = 1, failed_iterations = 0,
    power = 0.5, rank_correct_prob = 0.4,
    target_order_correct_prob = 0.8, avg_bias = 0, avg_abs_bias = 0
  )
  general <- structure(
    list(
      summary = summary,
      parameters = list(
        S = 1, target_contrast = c("B", "C"),
        network_design = data.frame(t1 = "B", t2 = "C", k = 1, OR = 2),
        true_log_ORs = c(B = 0, C = log(2))
      )
    ),
    class = "nma_power_sim"
  )
  three_node <- structure(
    list(
      summary = summary,
      parameters = list(S = 1, k_ab = 1, k_ac = 1, k_bc = 1)
    ),
    class = c("nma_power_sim_3node", "nma_power_sim")
  )

  expect_output(print(general), "Target Contrast Focus")
  expect_output(print(general), "B\\s+C\\s+1")
  expect_output(print(three_node), "k_ab.*1")
})

test_that("parallel operator is available inside the package namespace", {
  expect_true(exists("%dopar%", envir = asNamespace("NMAPower"), inherits = TRUE))
})
