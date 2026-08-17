test_that("general simulation validates designs before fitting", {
  valid <- data.frame(
    t1 = c("A", "A", "B"),
    t2 = c("B", "C", "C"),
    k = c(1, 1, 0),
    OR = c(1.2, 1.4, 1.4 / 1.2)
  )

  expect_error(
    nma_power_sim(network_design = valid[, -4], n.cores = 1),
    "must contain columns"
  )
  expect_error(
    nma_power_sim(network_design = transform(valid, k = c(1, -1, 0)), n.cores = 1),
    "non-negative integers"
  )
  inconsistent <- valid
  inconsistent$OR[3] <- 9
  expect_error(
    nma_power_sim(network_design = inconsistent, n.cores = 1),
    "inconsistent"
  )
})

test_that("three-node simulation requires an estimable B-C contrast", {
  expect_error(
    nma_power_sim_3node(k_ab = 1, k_ac = 0, k_bc = 0, n.cores = 1),
    "connected network containing B and C"
  )
})

test_that("post-hoc data validation catches invalid arm data", {
  bad <- data.frame(
    study.id = c(1, 1), treatment.id = c(1, 2),
    sample.size = c(10, 10), response = c(4, 11)
  )
  expect_error(
    nma_power_posthoc(bad, c(1, 2), n.cores = 1),
    "valid, non-missing binomial counts"
  )
})
