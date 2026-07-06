test_that("alpha_hat sums to 1", {
  dataframe <- make_simple_df()
  alpha_hat <- mle_fit(dataframe, D=3)$estimator$alpha
  expect_equal(sum(alpha_hat), 1)
})