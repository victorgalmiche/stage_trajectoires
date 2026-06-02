library(doParallel)
library(foreach)

source('src/tree_construction.R')

random_forest <- function(dataframe, covariates, pval_algo, 
                          n_trees=100, min_obs, min_leaf, alpha, max_depth,
                          max_features='sqrt', max_samples) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  forest <- foreach(
    i = 1:n_tres,
    .combine = c,
    .export = c("build_tree")
  ) %dopar% { 
    sample_indices <- sample(nrow(dataframe), size = max_samples, replace = TRUE)
    bootstrap_sample <- dataframe[sample_indices, ]
    
    build_tree(bootstrap_sample, covariates, pval_algo, min_obs, min_leaf, 
               alpha, max_depth)
  }
  
  return(forest)
}