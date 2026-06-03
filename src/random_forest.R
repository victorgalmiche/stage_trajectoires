library(doParallel)
library(foreach)

source('src/tree_construction.R')
source('src/two_samples_test.R')

random_forest <- function(dataframe, covariates, pval_algo, 
                          n_trees=100, min_obs, min_leaf, alpha, max_depth,
                          max_features='sqrt', max_samples) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  forest <- foreach(
    i = 1:n_trees,
    .combine = c,
    .export = c("build_tree", "find_best_split", "best_split_categorical",
                "best_split_numeric", "generate_bipartitions")
  ) %dopar% { 
    sample_indices <- sample(nrow(dataframe), size = max_samples, replace = TRUE)
    bootstrap_sample <- dataframe[sample_indices, ]
    
    alg <- function(df1, df2) {
      likelihood_ratio_test(df1, df2, 6, law_sojourn='exponential')
    }
    
    build_tree(bootstrap_sample, covariates[sample_indices, ], alg,
               max_features, min_obs, min_leaf, alpha, max_depth)
  }
  
  return(forest)
}
