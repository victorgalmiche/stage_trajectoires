
source('src/data/extract_data.R')

source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')

source('src/two_samples_test.R')

source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')

source('src/visualization.R')


law_sojourn <- 'exponential'

# Tree construction
system.time({
  tree <- build_tree(dataframe, covariates, weights=NULL, 
                     D, law_sojourn, permutation_test, 
                     min_leaf = 100, alpha = 0.05, max_depth = 5)
})
plot_tree(tree)


# RF construction
system.time({
  forest <- random_forest(dataframe, covariates, weights=NULL, 
                          D, law_sojourn, permutation_test, 
                          min_leaf = 100, alpha = 0.05)
  saveRDS(forest, file = "forest_temp.rds")
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


