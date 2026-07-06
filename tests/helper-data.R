source_from_root <- function(path) {
  source(file.path(testthat::test_path("..", ".."), path))
}

source_from_root("src/semi_markov/mle_estimation.R")
source_from_root("src/semi_markov/synthesis_data_generation.R")
source_from_root("src/two_samples_test.R")
source_from_root("src/random_forest/tree_construction.R")

make_toy_trajectories <- function(n = 30, D = 4, M = 5, seed = 1) {
  set.seed(seed)
  theta <- generate_theta(D, "exponential")
  generate_dataset_H0(theta, "exponential", n %/% 2, n - n %/% 2, M)
}

make_toy_covariates <- function(n, signal = FALSE) {
  data.frame(
    x_num = rnorm(n),
    x_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
}