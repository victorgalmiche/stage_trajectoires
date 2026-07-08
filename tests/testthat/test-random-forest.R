n <- 30
D <- 3
law_sojourn <- "exponential"

trajectories <- make_toy_trajectories(n, D)
covariates <- make_toy_covariates(n)
weights <- rep(1, n)

n_trees <- 5
forest <- random_forest(trajectories, covariates, weights,
                        D, law_sojourn, likelihood_ratio_test, 
                        n_trees)

test_that("correct number of trees", {
  expect_length(forest, n_trees)
})

test_that("checking that trees are actual trees", {
  for (i in seq_len(n_trees)){
    expect_true(forest[[i]]$type %in% c('node', 'leaf'))
  }
})

test_that("oob ids are present and greater than 0", {
  for (i in seq_len(n_trees)){
    expect_gt(forest[[i]]$oob_ids, 0)
  }
})
