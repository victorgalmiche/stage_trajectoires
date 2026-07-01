source('src/data/extract_data.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')

law_sojourn <- 'exponential'
alg <- function(df1, df2) {
  permutation_test(df1, df2, D, weights, law_sojourn)
}

df_PHD1 <- subset(dataframe, id %in% which(covariates$PHD_NEW==1))
little_df <- subset(dataframe, id<200)

tree <- build_tree(little_df, covariates, alg, max_depth=2)
