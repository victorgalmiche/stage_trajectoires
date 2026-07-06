### REGRESSION TREE CONSTRUCTION
# Split criterion: p-value of the two-sample test

# Data structure:
# - dataframe: data.frame, columns are id, state and time
# - covariates: data.frame, each row are the covariate for an individual
# - weights: vector of weights for each individual
# The id in dataframe refer to a row number in covariates and weights

source('src/semi_markov/mle_estimation.R')

# Find the best split among the different covariates

# Helper function to generate bipartitions 
generate_bipartitions <- function(n_levels) {
  parts <- list()
  for (k in 1:(2^(n_levels-1) - 1)) {
    bits <- as.integer(intToBits(k))[1:n_levels]
    left  <- which(bits==1)
    right <- which(bits==0)
    parts[[length(parts) + 1]] <- list(left = left, right = right)
  }
  parts
}

best_split_categorical <- function(dataframe, covariate, min_leaf, 
                                   pvalue_algo, weights, D, law_sojourn) {
  best <- list(pval=1, left_levels=NULL, right_levels=NULL)
  
  ids <- unique(dataframe$id) # Selecting only the ids of the individuals in the current dataframe
  levs <- levels(droplevels(covariate[ids])) # And the corresponding possible levels
  n_levels <- length(levs) # Counting the number of different levels
  
  if (n_levels > 1){
    parts <- generate_bipartitions(n_levels)
    
    # Iterating through the different partitions
    for (partition in parts){
      left_ids <- which(covariate %in% levs[partition$left])
      
      df_left <- subset(dataframe, id %in% left_ids)
      df_right <- subset(dataframe, !(id %in% left_ids))
      
      # If not enough values, don't take this threshold
      if (length(unique(df_left$id)) < min_leaf || 
          length(unique(df_right$id)) < min_leaf) next
      
      pval <- pvalue_algo(df_left, df_right, D, weights, law_sojourn)
      if (pval < best$pval){
        best <- list(pval=pval, 
                     left_levels=levs[partition$left], 
                     right_levels=levs[partition$right])
      }
    }
  }
  best
  
}

# Find the best split for a numeric variable
best_split_numeric <- function(dataframe, covariate, min_leaf, 
                               pvalue_algo, weights, D, law_sojourn) {
  best <- list(pval=1, threshold=NULL)
  
  ids <- unique(dataframe$id) # Selecting only the ids of the individuals in the current dataframe
  sorted_values <- sort(unique(covariate[ids])) # Sorting the values taken by the covariate
  n_values <- length(sorted_values) # Counting the number of different values
  
  if (n_values > 1){
    # The thresholds are the mean of successive values
    thresholds <- (sorted_values[-n_values] + sorted_values[-1]) / 2
    
    # Iterating through the different thresholds
    for (thresh in thresholds){
      left_ids <- which(covariate < thresh)
      
      df_left <- subset(dataframe, id %in% left_ids)
      df_right <- subset(dataframe, !(id %in% left_ids))
      
      # If not enough values, don't take this threshold
      if (length(unique(df_left$id)) < min_leaf || 
          length(unique(df_right$id)) < min_leaf) next
      
      pval <- pvalue_algo(df_left, df_right, D, weights, law_sojourn)
      if (pval < best$pval){
        best <- list(pval=pval, threshold=thresh)
      }
    }
  }
  best
}

find_best_split <- function(dataframe, covariates, min_leaf, 
                            pvalue_algo, weights, D, law_sojourn){
  best <- list(pval=1)
  for (var in names(covariates)){
    covariate <- covariates[[var]] # Extracting the column of the covariate
    
    # Selecting the good type of covariate
    if (is.numeric(covariate) || is.integer(covariate)){
      best_split <- best_split_numeric(dataframe, covariate, min_leaf, 
                                       pvalue_algo, weights, D, law_sojourn)
      split_type <- 'numeric'
    } else {
      best_split <- best_split_categorical(dataframe, covariate, min_leaf, 
                                           pvalue_algo, weights, D, law_sojourn)
      split_type <- 'categorical'
    }
    
    # Comparing w/ the best curent p-value
    if(best_split$pval < best$pval){
      best <- best_split
      best$var <- var
      best$type <- split_type
    }
  }
  best
} 


# Building the tree recursively 
build_tree <- function(dataframe, covariates, weights, 
                       D, law_sojourn, pvalue_algo,
                       max_features = 1, min_leaf = 5, alpha = 1, 
                       max_depth = Inf, depth = 0) {
  
  # Population of the current node and its size
  population <- unique(dataframe$id)
  pop_size <- length(population)
  
  # Attaining the maximum depth
  if (depth >= max_depth) {
    estimation <- mle_fit(dataframe, D, weights, law_sojourn)
    return(list(type = "leaf", estimator = estimation$estimator))
  }
  
  # Random selection of features among the covariates table
  size_sample <- as.integer(floor(max_features*ncol(covariates)))
  sample_cols <- sample(ncol(covariates), size = size_sample)
  sample_features <- covariates[, sample_cols]
  
  # Finding the best split
  best <- find_best_split(dataframe, sample_features, min_leaf, 
                          pvalue_algo, weights, D, law_sojourn)
  
  # No significant split
  if (is.null(best$var) || best$pval >= alpha) {
    estimation <- mle_fit(dataframe, D, weights, law_sojourn)
    return(list(type = "leaf", estimator = estimation$estimator))
  }
    
  # Split and recursively construct the subtree
  left_ids <- switch(best$type, 
                     categorical= (which(covariates[[best$var]] %in% 
                                           best$left_levels)),
                     numeric = (which(covariates[[best$var]] < best$threshold)))
  right_ids <- switch(best$type, 
                      categorical= (which(covariates[[best$var]] %in% 
                                            best$right_levels)),
                      numeric = (which(covariates[[best$var]] >= best$threshold)))
  
  df_left <- subset(dataframe, id %in% left_ids)
  df_right <- subset(dataframe, id %in% right_ids)
  
  list(
    type = "node",
    split = best,
    left = build_tree(df_left, covariates, weights, 
                      D, law_sojourn, pvalue_algo, 
                      max_features, min_leaf, alpha, 
                      max_depth, depth + 1),
    right = build_tree(df_right, covariates, weights,
                       D, law_sojourn, pvalue_algo, 
                       max_features, min_leaf, alpha, 
                       max_depth, depth + 1)
  )
}


# Function to get the leaf an observation belongs to
# obs is a row of covariates
get_leaf <- function(node, obs) {
  if (node$type == 'leaf') return(node)
  
  val <- obs[[node$split$var]]
  
  goes_left <- if (node$split$type == 'categorical') {
    val %in% node$split$left_levels
  } else {
    val < node$split$threshold       
  }
  
  if (goes_left) get_leaf(node$left, obs) else get_leaf(node$right, obs)
}

attach_node_population <- function(node, dataframe, covariates) {
  node$population <- unique(dataframe$id)
  if (node$type == 'node') {
    left_ids <- switch(node$split$type, 
                       categorical= (which(covariates[[node$split$var]] %in% 
                                             node$split$left_levels)),
                       numeric = (which(covariates[[node$split$var]] < node$split$threshold)))
    right_ids <- switch(node$split$type, 
                        categorical= (which(covariates[[node$split$var]] %in% 
                                              node$split$right_levels)),
                        numeric = (which(covariates[[node$split$var]] >= node$split$threshold)))
    
    df_left <- subset(dataframe, id %in% left_ids)
    df_right <- subset(dataframe, id %in% right_ids)
    node$left <- attach_node_population(node$left, df_left, covariates)
    node$right <- attach_node_population(node$right, df_right, covariates)
  }
  node
}

# Attach the Semi-Markov estimator to each leaf
attach_leaf_estimators <- function(node, dataframe, D, weights, law_sojourn) {
  if (node$type == 'leaf') {
    # Leaf: compute and attach estimator
    estimation <- mle_fit(subset(dataframe, id %in% node$population),
                          D, weights, law_sojourn)
    node$estimator <- estimation$estimator
  } else {
    node$left <- attach_leaf_estimators(node$left,  dataframe, D, weights, law_sojourn)
    node$right <- attach_leaf_estimators(node$right, dataframe, D, weights, law_sojourn)
  }
  node
}

# negative log-likelihood of an observation 
smooth_P <- function(P, epsilon = 1e-6) {
  P_smooth <- P + epsilon
  P_smooth / rowSums(P_smooth)
}

smooth_alpha <- function(alpha, epsilon = 1e-6) {
  alpha_smooth <- alpha + epsilon
  alpha_smooth / sum(alpha_smooth)
}

neg_log_lik <- function(tree, obs_id, dataframe, covariates, 
                        D, weights=NULL, law_sojourn='gamma') {
  obs <- covariates[obs_id, ]
  trajectory_df <- subset(dataframe, id==obs_id)
  
  node <- get_leaf(tree, obs)
  
  # To ensure non-zero entries 
  alpha_smooth <- smooth_alpha(node$estimator$alpha)
  P_smooth     <- smooth_P(node$estimator$P)
  
  ll_alpha <- log_likelihood_alpha(trajectory_df, alpha_smooth, weights)
  ll_P <- log_likelihood_P(trajectory_df, P_smooth, weights)
  ll_omega <- log_likelihood_omega(trajectory_df, node$estimator$omega, 
                                   weights, law_sojourn)
  
  - ll_alpha - ll_P - ll_omega
}


