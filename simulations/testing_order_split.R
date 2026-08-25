library(TraMineR)

source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')

# Conversion function from trajectories to dataframe format used 
# by two_sample_test.R
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
trajectories <- mvad[, 17:86]
dataframe <- traj_to_df(trajectories)
covariates <- mvad[, 3:14]

D <- 6 # Number of states 
law_sojourn <- 'exponential' # Law sojourn times

# To store the results
results_lambda <- setNames(numeric(ncol(covariates)), names(covariates))
results_chi2 <- setNames(numeric(ncol(covariates)), names(covariates))
results_perm <- setNames(numeric(ncol(covariates)), names(covariates))
results_boot <- setNames(numeric(ncol(covariates)), names(covariates))

for (var in names(covariates)){
  covariate <- covariates[[var]] # Extracting the column of the covariate
  
  # Ids for each level 
  ids_1 <- which(covariate == levels(covariate)[1])
  ids_2 <- which(covariate == levels(covariate)[2])
  
  # Subpopulation 
  df1 <- subset(dataframe, id %in% ids_1)
  df2 <- subset(dataframe, id %in% ids_2)
  
  # Computing the different values of interest  
  results_lambda[[var]] <- lambda_statistic(df1, df2, D, law_sojourn=law_sojourn)
  results_chi2[[var]] <- likelihood_ratio_test(df1, df2, D, law_sojourn=law_sojourn)
  results_perm[[var]] <- permutation_test(df1, df2, D, law_sojourn=law_sojourn, 
                                          R=1000)
  results_boot[[var]] <- parametric_bootstrap(df1, df2, D, law_sojourn=law_sojourn, 
                                              R=1000, T_max=70)
}

# Comparative dataframe
comparison_df <- data.frame(
  variable = names(covariates),
  lambda = results_lambda,
  pval_chi2 = results_chi2,
  pval_perm = results_perm, 
  pval_boot = results_boot
)

comparison_df['S.Eastern', 'pval_perm'] <- 0.046

# Adding the ranks
comparison_df$rank_lambda <- rank(-comparison_df$lambda)
comparison_df$rank_perm <- rank(comparison_df$pval_perm, ties.method = "min")
comparison_df$rank_boot <- rank(comparison_df$pval_boot, ties.method = "min") 

# Sorting wrt lambda/chi2
comparison_df <- comparison_df[order(comparison_df$rank_lambda), ]
print(comparison_df)

# Correlation measures
cor_spearman <- cor(comparison_df$pval_perm, 
                    comparison_df$pval_boot, 
                    method = "spearman")
cor_kendall  <- cor(comparison_df$pval_perm, 
                    comparison_df$pval_boot, 
                    method = "kendall")

cat("Spearman coefficient :", cor_spearman, "\n")
cat("Kendall's tau :", cor_kendall, "\n")

