dataframe <- make_simple_df()
fit <- mle_fit(dataframe, D=3)

test_that("alpha_hat sums to 1", {
  alpha_hat <- fit$estimator$alpha
  expect_equal(sum(alpha_hat), 1)
})

test_that("rows of P_hat sum to 1", {
  P_hat <- fit$estimator$P
  expect_equal(rowSums(P_hat), rep(1, nrow(P_hat)))
})

test_that("diagonal of P_hat is 0", {
  P_hat <- fit$estimator$P
  expect_equal(diag(P_hat), rep(0, nrow(P_hat)))
})