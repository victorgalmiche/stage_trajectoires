
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

# Tree construction
system.time({
  min_leaf <- as.integer(floor(length(unique(dataframe$id))/20)) # 5% of the total nb of ind
  tree <- build_tree(dataframe, covariates, weights, 
                     D, law_sojourn, likelihood_ratio_test, 
                     min_leaf = min_leaf, alpha = 0.05, max_depth = 5)
})

# RF construction
system.time({
  forest <- random_forest(dataframe, covariates, weights, 
                          D, law_sojourn, likelihood_ratio_test, 
                          min_leaf = min_leaf, alpha = 0.05)
})

# MDI
system.time({
  ranking_MDI <- MDI_all(forest, covariates)
})

barplot(ranking_MDI, 
        main = "Mean Decrease Impurity for each covariate",
        ylab = "MDI", 
        ylim = c(0,1),
        col = "blue", 
        las = 2)


# MDA
system.time({
  ranking_MDA <- MDA_all(forest, dataframe, covariates, 
                         D, weights, law_sojourn)
})

barplot(ranking_MDA, 
        main = "Mean Decrease Acurracy for each covariate",
        ylab = "MDA",
        col = "blue", 
        las = 2)


