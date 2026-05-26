### REGRESSION TREE CONSTRUCTION
# Split criterion: p-value of the two-sample test

# Data structure:
# - trajectories: data.frame, each row is a trajectory
# - covariates: data.frame, each row are the covariate for an individual
# For now, we don't need an id column, we give id as the row number for each individual


# Conversion function from trajectories to dataframe format used by two_sample_test.R
traj_to_df <- function(trajectories) {
  res <- list()
  
  for (i in seq_len(nrow(trajectories))) {
    states <- as.integer(trajectories[i, ]) # Extracting the trajectory of i
    changes <- c(TRUE, diff(states) != 0) # Breakpoint detection
    episodes <- cumsum(changes) # Episode numbering
    durations <- table(episodes) # Duration of each episode
    episode_states <- states[changes] # Corresponding state
    
    # Construction of the resulting dataframe
    res[[i]] <- data.frame(
      id = i,
      state = episode_states,
      time = as.integer(durations)
    )
  }
  
  do.call(rbind, res)
}

# Find the best split among the different covariates
# We suppose that each covariates are binary categorical (2 modalities)
best_split_categorical <- function(df, covariate, min_leaf, pvalue_algo) {
  best <- list(pval = 1, var = NULL, left_level = NULL, right_level = NULL)
  
  levs <- levels(covariate) # For now, a max of 2 levels 
  left_ids <- which(covariate==levs[1])
  
  df_left <- subset(df, id %in% left_ids)
  df_right <- subset(df, !(id %in% left_ids))
  
  if (length(unique(df_left$id))>min_leaf &&
      length(unique(df_right$id))>min_leaf) {
    
    pval <- pvalue_algo(df_left, df_right)
    if (pval < best$pval)
      best <- list(pval = pval, 
                   var = var, 
                   left_level = levs[1], 
                   right_level = levs[2])
  }
  best
}

# Find the best split for a numeric variable
best_split_numeric <- function(df, covariate, min_leaf, pvalue_algo) {
  best <- list(pval=1, threshold=NULL)
  
  ids <- unique(df$id) # Selecting only the ids of the individuals in the current dataframe
  sorted_values <- sort(unique(covariate[ids])) # Sorting the values taken by the covariate
  n_values <- length(sorted_values) # Counting the number of different values
  
  if (n_values > 1){
    # The thresholds are the mean of successive values
    thresholds <- (sorted_values[-n_values] + sorted_values[-1]) / 2
    
    # Iterating through the different thresholds
    for (thresh in thresholds){
      left_ids <- which(covariate < thresh)
      
      df_left <- subset(df, id %in% left_ids)
      df_right <- subset(df, !(id %in% left_ids))
      
      # If not enough values, don't take this threshold
      if (length(unique(df_left$id)) < min_leaf || 
          length(unique(df_right$id)) < min_leaf) next
      
      pval <- pvalue_algo(df_left, df_right)
      if (pval < best$pval){
        best <- list(pval=pval, threshold=thresh)
      }
    }
  }
  best
}

find_best_split <- function(df, covariates, min_leaf, pvalue_algo){
  best <- list(pval=1)
  for (var in names(covariates)){
    covariate <- covariates[[var]] # Extracting the column of the covariate
    
    # Selecting the good type of covariate
    if (is.numeric(covariate) || is.integer(covariate)){
      best_split <- best_split_numeric(df, covariate, min_leaf, pvalue_algo)
      split_type <- 'numeric'
    } else {
      best_split <- best_split_categorical(df, covariate, min_leaf, pvalue_algo)
      split_type <- 'categorical'
    }
    
    # Comparing w/ the best curent p-value
    if(best_split$pval <- best$pval){
      best <- best_split
      best$var <- var
      best$type <- split_type
    }
  }
  best
} 


# Building the tree recursively 
build_tree <- function(df, covariates, pvalue_algo, min_obs = 20, min_leaf = 5,
                       alpha = 0.05, max_depth = 5, depth = 0) {
  population <- unique(df$id)
  pop_size <- length(population)
  
  # Stopping criterion
  if (depth >= max_depth || pop_size < min_obs)
    return(list(type = "leaf", population = population, n = pop_size))
  
  best <- find_best_split(df, covariates, min_leaf, pvalue_algo)
  
  # No significant split
  if (is.null(best$var) || best$pval >= alpha)
    return(list(type = "leaf", population = population, n = pop_size))
  
  # Split and recursively construct the subtree
  left_ids <- switch(best$type, 
                     categorical= (which(covariates[[best$var]] == 
                                           best$left_level)),
                     numeric = (which(covariates[[best$var]] < best$threshold)))
  right_ids <- switch(best$type, 
                      categorical= (which(covariates[[best$var]] == 
                                            best$right_level)),
                      numeric = (which(covariates[[best$var]] >= best$threshold)))
  
  df_left <- subset(df, id %in% left_ids)
  df_right <- subset(df, id %in% right_ids)
  
  list(
    type = "node",
    population = population,
    n = pop_size,
    split = best,
    left = build_tree(df_left, covariates, pvalue_algo, min_obs, min_leaf, alpha, max_depth, depth + 1),
    right = build_tree(df_right, covariates, pvalue_algo, min_obs, min_leaf, alpha, max_depth, depth + 1)
  )
}

# Print function 
print_tree <- function(node, indent = 0) {
  pad <- strrep("  ", indent)
  if (node$type == "leaf") {
    cat(pad, "└─ Leaf:  (n =", node$n, ")\n")
  } else {
    cat(pad, "├─ [", node$var, "] p =", formatC(node$pval, digits = 3, format = "e"), "\n")
    cat(pad, "  ├─ "); print_tree(node$left,  indent + 2)
    cat(pad, "  └─ "); print_tree(node$right, indent + 2)
  }
}

# Use example w/ mvad data
source('src/two_samples_test.R')
library(TraMineR)
#data(mvad)
#trajectories <- mvad[, 17:86]
#covariates <- mvad[, 3:14]
traj_df <- traj_to_df(trajectories)


data(biofam)
trajectories <- biofam[, 10:25]
traj_df <- traj_to_df(trajectories) 
traj_df$state <- traj_df$state + 1 # to have state number beginning at 1
covariates <- biofam[, 2:9]

D <- 8

# Function for p_value computation
alg <- function(df1, df2){
  likelihood_ratio_test(df1, df2, D, 'exponential')
  # permutation_test(df1, df2, D, 'exponential')
}

# Tree construction
tree <- build_tree(traj_df, covariates[,1:2], alg, alpha=1.1)

# Affichage
cat("=== ARBRE DE TRAJECTOIRES ===\n\n")
print_tree(tree)

