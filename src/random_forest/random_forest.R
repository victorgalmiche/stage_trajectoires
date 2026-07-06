library(doParallel)
library(foreach)

source('src/random_forest/tree_construction.R')

random_forest <- function(dataframe, covariates, weights, 
                          D, law_sojourn, pvalue_algo,
                          n_trees=500, max_samples=1, 
                          max_features=1/3, min_leaf=5, alpha=0.1) {
  
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
    .errorhandling = "pass",
    .export = c(
      "dataframe", "covariates", "weights", "D", "law_sojourn", "pvalue_algo",
      "max_features", "min_leaf", "alpha", "ids", "boot_size")
  ) %dopar% { 
    
    # Bootstrap sample
    boot_ids <- sample(ids, size = boot_size, replace = TRUE)
    # Need to do this to keep multiple-selected ids
    idx <- unlist(lapply(boot_ids, function(id) which(dataframe$id == id)))
    bootstrap_sample <- do.call(rbind, lapply(seq_along(boot_ids), function(i) {
      d <- dataframe[dataframe$id == boot_ids[i], ]
      d$id <- i          # new unique id per draw, 1:boot_size
      d
    }))
    
    # Keep only the interesting weights and covariates and reorder
    boot_weights <- if (is.null(weights)) NULL else weights[boot_ids]
    boot_covariates <- covariates[boot_ids, , drop = FALSE]
    
    # Tree construction 
    tree <- build_tree(bootstrap_sample, boot_covariates, boot_weights,
                       D, law_sojourn, pvalue_algo,
                       max_features, min_leaf, alpha)
    tree$oob_ids <- setdiff(ids, boot_ids)
    tree
  }
  
  return(forest)
}

