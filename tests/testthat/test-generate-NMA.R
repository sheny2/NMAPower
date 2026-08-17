test_that("generate_NMA returns valid binomial arm data", {
  set.seed(123)
  dat <- generate_NMA(5, 3, c(0.2, 0.5, 0.8), c(20, 30))

  expect_named(dat, c("study.id", "treatment.id", "sample.size", "response"))
  expect_setequal(unique(dat$study.id), 1:5)
  expect_true(all(table(dat$study.id) >= 2L))
  expect_true(all(dat$sample.size >= 20 & dat$sample.size <= 30))
  expect_true(all(dat$response >= 0 & dat$response <= dat$sample.size))
})

test_that("generate_NMA rejects malformed inputs", {
  expect_error(generate_NMA(2, 1, 0.5), "num_treatments")
  expect_error(generate_NMA(2, 2, c(0.2, NA_real_)), "probabilities")
  expect_error(generate_NMA(2, 2, c(0.2, 0.5), c(30, 20)), "sample_size_range")
})
