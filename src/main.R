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

### TEST RANDOM FOREST AND TREES
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')
source('src/two_samples_test.R')

library(doParallel)
library(foreach)
library(TraMineR)


# mvad data
data(mvad)
trajectories <- mvad[, 17:86]
covariates <- mvad[, 3:14]
traj_df <- traj_to_df(trajectories)

# Number of states 
D <- 6

weights <- mvad[, 2]
law_sojourn <- 'exponential'

# Function for p_value computation
alg <- function(df1, df2){
  likelihood_ratio_test(df1, df2, D, weights, law_sojourn)
  #permutation_test(df1, df2, D, weights, law_sojourn)
}


# Tree construction
# tree <- build_tree(traj_df, covariates, alg)

# And a random forest
rf <- random_forest(traj_df, covariates, alg, 100, 20, 5, 0.05, 5, 'sqrt', 200)

system.time({
  ranking_MDA <- MDA_all(rf, traj_df, covariates, D, weights, law_sojourn)
})


system.time({
  ranking_MDI <- MDI_all(rf, covariates)
})



barplot(ranking_MDA, 
        main = "Chi^2 test and Exponential Law",
        ylab = "MDA", 
        col = "blue", 
        las = 2)

barplot(ranking_MDI, 
        main = "Chi^2 test and Exponential Law",
        ylab = "MDI", 
        ylim = c(0,1),
        col = "blue", 
        las = 2)

# # biofam
# data(biofam)
# trajectories <- biofam[, 10:25]
# traj_df <- traj_to_df(trajectories) 
# traj_df$state <- traj_df$state + 1 # to have state number beginning at 1
# covariates <- biofam[, 5:9]
# 
# # actcal
# data(actcal)
# trajectories <- actcal[, 13:24]
# traj_df <- traj_to_df(trajectories)
# traj_df$state <- traj_df$state - 5 # state number beginning at 1
# covariates <- actcal[, 2:12]
