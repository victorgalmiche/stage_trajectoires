
source('src/data/extract_data.R')

source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')

source('src/two_samples_test.R')
source('src/visualization.R')

source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')


law_sojourn <- 'exponential'

little_df <- subset(dataframe, id<200)
df_PHD1 <- subset(dataframe, id %in% which(covariates$PHD_NEW==1))

system.time({
  forest <- random_forest(df_PHD1, covariates, weights, 
                          D, law_sojourn, likelihood_ratio_test)
})

system.time({
  ranking_MDA <- MDA_all(forest, df_PHD1, covariates, 
                         D, weights, law_sojourn)
})

system.time({
  ranking_MDI <- MDI_all(forest, df_PHD1, covariates)
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
