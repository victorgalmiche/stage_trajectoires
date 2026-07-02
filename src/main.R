source('src/data/extract_data.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')

law_sojourn <- 'exponential'
alg <- function(df1, df2) {
  permutation_test(df1, df2, D, weights, law_sojourn)
}

little_df <- subset(dataframe, id<200)
df_PHD1 <- subset(dataframe, id %in% which(covariates$PHD_NEW==1))

system.time({
  tree <- build_tree(dataframe, covariates, alg, max_depth=5, alpha=0.1)
})

source('src/random_forest/random_forest.R')
system.time({
  forest <- random_forest(df_PHD1, covariates, alg)
})


source('src/visualization.R')
plot_tree(tree)
plot_tree(forest[[1]])


source('src/random_forest/variable_importance.R')
system.time({
  ranking_MDA <- MDA_all(forest, df_PHD1, covariates, 
                         D, weights, law_sojourn)
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
