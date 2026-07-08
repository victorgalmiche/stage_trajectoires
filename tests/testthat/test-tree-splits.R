n <- 30
D <- 3
law_sojourn <- "exponential"

trajectories <- make_toy_trajectories(n)
covariates <- make_toy_covariates(n)
weights <- rep(1, n)

min_leaf <- 5

tree <- build_tree(trajectories, covariates, weights, 
                   D, law_sojourn, likelihood_ratio_test, 
                   min_leaf=min_leaf)

tree_with_pop <- attach_node_population(tree, trajectories, covariates)

test_that("no leaf with more than min_leaf", {
  for (i in seq_len(n)){
    leaf <- get_leaf(tree_with_pop, covariates[i, ])
    expect_lt(length(leaf$population), min_leaf)
  }
})


# tree_with_new_estimators <- attach_leaf