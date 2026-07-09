source_from_root <- function(path) {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/")
  source(file.path(root, path), local = FALSE)
}

# project sources
source_from_root("src/semi_markov/mle_estimation.R")
source_from_root("src/semi_markov/synthesis_data_generation.R")
source_from_root("src/two_samples_test.R")
source_from_root("src/random_forest/tree_construction.R")
source_from_root("src/random_forest/random_forest.R")
source_from_root("src/random_forest/variable_importance.R")

# helper functions to create toy data
make_simple_df <- function() {
  data.frame(
    id    = c(1, 1, 2, 2),
    state = c(1, 2, 1, 3),
    time  = c(5, 3, 2, 4)
  )
}

make_toy_trajectories <- function(n = 30, D = 3, M = 5, 
                                  law_sojourn = "exponential", seed = 1) {
  set.seed(seed)
  theta <- generate_theta(D, law_sojourn)
  generate_dataset_H0(theta, law_sojourn, n %/% 2, n - n %/% 2, M)
}

make_toy_covariates <- function(n) {
  data.frame(
    x_num = rnorm(n),
    x_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
}

make_leaf <- function(D=3, law_sojourn="exponential"){
  list(type = 'leaf', estimator = generate_theta(D, law_sojourn))
}

make_node <- function(n, pval, var, left, right){
  list(type = 'node', n = n, 
       split = list(pval = pval, var = var), 
       left = left, 
       right = right)
}