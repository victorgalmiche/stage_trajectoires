library(doParallel)
library(foreach)

source('src/tree_construction.R')
source('src/two_samples_test.R')

random_forest <- function(dataframe, covariates, pval_algo, 
                          n_trees=100, min_obs, min_leaf, alpha, max_depth,
                          max_features='sqrt', max_samples) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)
  
  clusterEvalQ(cl, {
    source('src/semi_markov/synthesis_data_generation.R')
    source('src/semi_markov/mle_estimation.R')
    source('src/two_samples_test.R')
    source('src/tree_construction.R')
  })
  ids <- unique(dataframe$id)
  forest <- foreach(
    i = 1:n_trees,
    .combine = c,
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
    
    # Wrap into a list for the aggregation
    list(build_tree(bootstrap_sample, covariates, pval_algo,
               max_features, min_obs, min_leaf, alpha, max_depth))
  }
  
  return(forest)
}


# TEST 
library(TraMineR)
data(mvad)
trajectories <- mvad[, 17:86]
covariates <- mvad[, 3:14]
traj_df <- traj_to_df(trajectories)

alg <- function(df1, df2) {
  likelihood_ratio_test(df1, df2, 6, law_sojourn='exponential')
}

rf <- random_forest(traj_df, covariates, alg, 100, 20, 5, 0.05, 5, 'sqrt', 200)


source('src/variable_importance.R')
