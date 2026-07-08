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
source('src/semi_markov/mle_estimation.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')

library(TraMineR)

data(mvad)
trajectories <- mvad[, 17:86]
covariates <- mvad[, 3:14]
dataframe <- traj_to_df(trajectories_mvad)


# Number of states 
D <- 5

weights <- rep(1, nrow(covariates))
law_sojourn <- 'exponential'


results_chi2 <- list()
results_perm <- list()

for (var in names(covariates)){
  covariate <- covariates[[var]] # Extracting the column of the covariate
    
  # Selecting the good type of covariate
  if (is.numeric(covariate) || is.integer(covariate)){
    best_split_chi2 <- best_split_numeric(dataframe, covariate, min_leaf=1, 
                                          likelihood_ratio_test, weights,
                                          D, law_sojourn)
    
    best_split_perm <- best_split_numeric(dataframe, covariate, min_leaf=1,
                                          permutation_test, weights,
                                          D, law_sojourn)
  } else {
    best_split_chi2 <- best_split_categorical(dataframe, covariate, min_leaf=1, 
                                              likelihood_ratio_test, weights, 
                                              D, law_sojourn)
    
    best_split_perm <- best_split_categorical(dataframe, covariate, min_leaf=1,
                                              permutation_test, weights, 
                                              D, law_sojourn)
  }
  results_chi2[[var]] <- best_split_chi2
  results_perm[[var]] <- best_split_perm
}

pval_chi2 <- sapply(results_chi2, function(x) x$pval)
pval_perm <- sapply(results_perm, function(x) x$pval)

# Construction d'un tableau comparatif
comparison_df <- data.frame(
  variable = names(covariates),
  pval_chi2 = pval_chi2,
  pval_perm = pval_perm
)

# Ajout des rangs (1 = p-valeur la plus faible = covariate la plus significative)
comparison_df$rank_chi2 <- rank(comparison_df$pval_chi2)
comparison_df$rank_perm <- rank(comparison_df$pval_perm)

# Tri selon le test chi2 pour visualiser facilement
comparison_df <- comparison_df[order(comparison_df$rank_chi2), ]
print(comparison_df)

# Mesures de concordance entre les deux classements
cor_spearman <- cor(comparison_df$pval_chi2, comparison_df$pval_perm, method = "spearman")
cor_kendall  <- cor(comparison_df$pval_chi2, comparison_df$pval_perm, method = "kendall")

cat("Corrélation de Spearman entre les classements :", cor_spearman, "\n")
cat("Corrélation de Kendall (tau) entre les classements :", cor_kendall, "\n")

