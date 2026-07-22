library(TraMineR) # To compare with Studer et al. methodology

source('src/semi_markov/mle_estimation.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')
source('src/visualization.R')

# Conversion function from trajectories to dataframe format 
# used by two_samples_test.R
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

# Charging mvad data
data(mvad)
trajectories_mvad <- mvad[, 17:86] # the trajectories 
dataframe <- traj_to_df(trajectories_mvad) # the dataframe associated to the trajectories
covariates <- mvad[, 3:14] # the covariates
weights <- mvad[, 2] # the weights

# Number of states and sojourn time law
D <- 6
law_sojourn <- 'exponential'
min_leaf <- as.integer(floor(nrow(covariates)/20)) # 5% of the total nb of ind

# Tree construction and visualization
tree <- build_tree(dataframe, covariates, weights,
                   D, law_sojourn, likelihood_ratio_test, 
                   min_leaf = min_leaf, alpha = 0.05, max_depth = 5)
plot_tree(tree)

# And a random forest
rf <- random_forest(dataframe, covariates, weights,
                    D, law_sojourn, likelihood_ratio_test,
                    min_leaf = min_leaf, alpha = 0.05)


# Evaluating variable importance
ranking_MDI <- MDI_all(rf, covariates)
ranking_MDA <- MDA_all(rf, dataframe, covariates, D, weights, law_sojourn)

barplot(ranking_MDI, 
        main = "MDIs of the covariates",
        ylab = "MDI", ylim = c(0, 1),
        col = "blue", 
        las = 2)

barplot(ranking_MDA, 
        main = "MDAs of the covariates",
        ylab = "MDA",
        col = "blue", 
        las = 2)




