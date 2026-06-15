source('src/semi_markov/mle_estimation.R')

# We begin by using 1-pvalue as a proxy for the impurity
MDI_tree <- function(node, covariate_name, n_root) {
  if (node$type == 'leaf') {
    return(0)
  }
  
  node_contribution <- if (node$split$var == covariate_name){
    (node$n/n_root)*(1-node$split$pval)
  } else {
    0
  }
  
  left_contribution <- MDI_tree(node$left, covariate_name, n_root)
  right_contribution <- MDI_tree(node$right, covariate_name, n_root)
  
  node_contribution + left_contribution + right_contribution
}


# Average MDI across all trees in the forest
MDI <- function(forest, covariate_name) {
  M <- length(forest)
  
  contributions <- vapply(forest, function(tree) {
    n_root <- tree$n
    MDI_tree(tree, covariate_name, n_root)
  }, numeric(1))
  
  mean(contributions)
}

# Rank all covariates by importance
MDI_all <- function(forest, covariates) {
  importance <- vapply(names(covariates), function(cov) {
    MDI(forest, cov)
  }, numeric(1))
  
  sort(importance, decreasing = TRUE)
}