source('src/data/extract_data.R')
source('src/two_samples_test.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')
source('src/visualization.R')

law_sojourn <- 'exponential'

little_df <- subset(dataframe, id<200)
df_PHD1 <- subset(dataframe, id %in% which(covariates$PHD_NEW==1))

# system.time({
#   tree <- build_tree(dataframe, covariates, alg, max_depth=5, alpha=0.1)
# })
# plot_tree(tree)


system.time({
  forest <- random_forest(little_df, covariates, weights, D, law_sojourn)
})

plot_tree(forest[[2]])


system.time({
  ranking_MDA <- MDA_all(forest, little_df, covariates, 
                         D, weights, law_sojourn)
})

system.time({
  ranking_MDI <- MDI_all(forest, covariates)
})

barplot(ranking_MDA, 
        main = "Permutation test and Exponential Law",
        ylab = "MDA",
        col = "blue", 
        las = 2)

barplot(ranking_MDI, 
        main = "Permutation test and Exponential Law",
        ylab = "MDI", 
        # ylim = c(0,1),
        col = "blue", 
        las = 2)
