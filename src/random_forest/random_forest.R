library(doParallel)
library(foreach)

source('src/random_forest/tree_construction.R')

random_forest <- function(dataframe, covariates, pval_algo, 
                          n_trees=500, max_samples=1, 
                          max_features=1/3, min_leaf=5) {
  
  # registering the CPUs for parallelization
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
  
  # Listing the individuals in the dataframe and their number
  ids <- unique(dataframe$id)
  n <- length(ids)
  
  # Size of the bootstrap samples
  boot_size <- as.integer(floor(max_samples*n))
  
  # Parallelized construction of the trees 
  forest <- foreach(
    i = 1:n_trees,
    .combine = list,
    .multicombine = TRUE,
    .maxcombine = n_trees,
    .errorhandling = "remove",
    .export = c(
      "dataframe", "covariates", "pval_algo",
      "max_features", "min_leaf", "ids", "boot_size")
  ) %dopar% { 
    
    # Bootstrap sample
    boot_ids <- sample(ids, size = boot_size, replace = TRUE)
    # Need to do this to keep multiple-selected ids
    idx <- unlist(lapply(boot_ids, function(id) which(dataframe$id == id)))
    bootstrap_sample <- dataframe[idx, ] 
    
    # Tree construction 
    build_tree(bootstrap_sample, covariates, pval_algo,
               max_features, min_leaf)
  }
  
  return(forest)
}

