library(doParallel)
library(foreach)

source('src/random_forest/tree_construction.R')

random_forest <- function(dataframe, covariates, pval_algo, 
                          n_trees=100, min_obs, min_leaf, alpha, max_depth,
                          max_features='sqrt', max_samples) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)
  
  clusterEvalQ(cl, {
    source('src/semi_markov/synthesis_data_generation.R')
    source('src/semi_markov/mle_estimation.R')
    source('src/random_forest/tree_construction.R')
    source('src/two_samples_test.R')
  })
  
  # Export variables from caller's environment that alg() needs
  clusterExport(cl, c("D", "weights", "law_sojourn"), 
                envir = parent.frame())
  
  ids <- unique(dataframe$id)
  forest <- foreach(
    i = 1:n_trees,
    .combine = list,
    .multicombine = TRUE,
    .maxcombine = n_trees,
    .errorhandling = "remove",
    .export = c(
      "dataframe", "covariates", "pval_algo", "ids", "max_samples",
      "min_obs", "min_leaf", "alpha", "max_depth", "max_features")
  ) %dopar% { 
    boot_ids <- sample(ids, size = max_samples, replace = TRUE)
    # Need to do this to keep multiple-selected ids
    idx <- unlist(lapply(boot_ids, function(id) which(dataframe$id == id)))
    bootstrap_sample <- dataframe[idx, ] 
    
    build_tree(bootstrap_sample, covariates, pval_algo,
               max_features, min_obs, min_leaf, alpha, max_depth)
  }
  
  return(forest)
}

