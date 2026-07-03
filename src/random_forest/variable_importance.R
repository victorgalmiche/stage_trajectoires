source('src/random_forest/tree_construction.R')


### MDI
# We begin by using 1-pvalue as a proxy for the impurity
MDI_tree <- function(node, covariate_name, n_root) {
  if (node$type == 'leaf') {
    return(0)
  }
  
  node_contribution <- if (node$split$var == covariate_name){
    n <- length(node$population)
    (n/n_root)*(1-node$split$pval)
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
    n_root <- length(tree$population)
    MDI_tree(tree, covariate_name, n_root)
  }, numeric(1))
  
  mean(contributions)
}

# Rank all covariates by importance
MDI_all <- function(forest, dataframe, covariates) {
  # Compute the population in each node
  forest <- lapply(forest, attach_node_population, dataframe, covariates)
  
  importance <- vapply(names(covariates), function(cov) {
    MDI(forest, cov)
  }, numeric(1))
  
  sort(importance, decreasing = TRUE)
}



### MDA 
oob_score <- function(tree, oob_ids, dataframe, covariates,
                      D, weights, law_sojourn) {
  scores <- vapply(oob_ids, function(obs_id) {
    tryCatch(
      neg_log_lik(tree, obs_id, dataframe, covariates, D, weights, law_sojourn),
      error = function(e) NA_real_   # guard against empty leaves
    )
  }, numeric(1))
  
  mean(scores, na.rm = TRUE)
}

MDA <- function(forest, covariate_name, dataframe, covariates, 
                D, weights, law_sojourn) {
  all_ids <- unique(dataframe$id)
  
  decreases <- vapply(forest, function(tree) {
    oob_ids <- tree$oob_ids
    if (length(oob_ids) == 0) return(NA_real_)
    
    base_score <- oob_score(tree, oob_ids, dataframe, 
                            covariates, D, weights, law_sojourn)
    
    # Permute the target covariate for OOB observations only
    covariates_permuted <- covariates
    covariates_permuted[oob_ids, covariate_name] <- 
      sample(covariates[oob_ids, covariate_name, drop=TRUE])
    
    perm_score <- oob_score(tree, oob_ids, dataframe, 
                            covariates_permuted, D, weights, law_sojourn)
    
    # Positive = permutation hurt = covariate was useful
    perm_score - base_score
  }, numeric(1))
  
  mean(decreases, na.rm = TRUE)
}


# Rank all covariates
MDA_all <- function(forest, dataframe, covariates, D, weights, law_sojourn) {
  # Then compute importance for each covariate
  importance <- vapply(names(covariates), function(cov) {
    MDA(forest, cov, dataframe, covariates, D, weights, law_sojourn)
  }, numeric(1))
  
  # And rank them
  sort(importance, decreasing = TRUE)
}