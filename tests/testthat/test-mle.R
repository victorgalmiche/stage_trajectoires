# Function to make a simple dataframe w/ 2 individuals and 3 different states
make_simple_df <- function() {
  data.frame(
    id    = c(1, 1, 2, 2),
    state = c(1, 2, 1, 3),
    time  = c(5, 3, 2, 4)
  )
}

test_that("alpha_hat sums to 1", {
  dataframe <- make_simple_df()
  alpha_hat <- mle_fit(dataframe, D=3)$estimator$alpha
  expect_equal(sum(alpha_hat), 1)
})