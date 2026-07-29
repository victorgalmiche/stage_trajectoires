
source('src/data/extract_data.R')

source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')

source('src/two_samples_test.R')
source('src/visualization.R')

source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')

source('src/visualization.R')


law_sojourn <- 'exponential'

df_PHD1 <- subset(dataframe, id %in% which(covariates$PHD==1))

system.time({
  min_leaf <- as.integer(floor(length(unique(dataframe$id))/20)) # 5% of the total nb of ind
  tree <- build_tree(dataframe, covariates, weights, 
                     D, law_sojourn, likelihood_ratio_test, 
                     min_leaf = min_leaf, alpha = 0.05, max_depth = 5)
})

system.time({
  forest <- random_forest(dataframe, covariates, weights, 
                          D, law_sojourn, likelihood_ratio_test, 
                          min_leaf = min_leaf, alpha = 0.05)
})

system.time({
  ranking_MDA <- MDA_all(forest, df_PHD1, covariates, 
                         D, weights, law_sojourn)
})

system.time({
  ranking_MDI <- MDI_all(forest, covariates)
})

barplot(ranking_MDA, 
        main = "Permutation test and Exponential Law - PHD=1",
        ylab = "MDA",
        col = "blue", 
        las = 2)

barplot(ranking_MDI, 
        main = "Permutation test and Exponential Law - PHD=1",
        ylab = "MDI", 
        # ylim = c(0,1),
        col = "blue", 
        las = 2)
