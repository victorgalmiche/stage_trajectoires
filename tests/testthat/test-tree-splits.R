n <- 30
D <- 3
law_sojourn <- "exponential"

trajectories <- make_toy_trajectories(n, D)
covariates <- make_toy_covariates(n)
weights <- rep(1, n)

min_leaf <- 5

tree <- build_tree(trajectories, covariates, weights, 
                   D, law_sojourn, likelihood_ratio_test, 
                   min_leaf=min_leaf)

tree_with_pop <- attach_node_population(tree, trajectories, covariates)
tree_with_new_estimators <- attach_leaf_estimators(tree_with_pop, trajectories,
                                                   D, weights, law_sojourn)

test_that("no leaf with more than min_leaf", {
  for (i in seq_len(n)){
    leaf <- get_leaf(tree_with_pop, covariates[i, ])
    expect_gte(length(leaf$population), min_leaf)
  }
})

test_that("verification of the estimators in the built tree", {
  for (i in seq_len(n)){
    leaf <- get_leaf(tree, covariates[i, ])
    leaf_2 <- get_leaf(tree_with_new_estimators, covariates[i, ])
    expect_equal(leaf_2$estimator, leaf$estimator)
  }
})



