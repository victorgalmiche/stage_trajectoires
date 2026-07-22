source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')
source('src/two_samples_test.R')
source('src/random_forest/tree_construction.R')
source('src/random_forest/random_forest.R')
source('src/random_forest/variable_importance.R')
source('src/visualization.R')

D <- 2
law_sojourn <- 'exponential'

# 3 groups of data
theta1 <- generate_theta(D, law_sojourn)
theta2 <- generate_theta(D, law_sojourn)
theta3 <- generate_theta(D, law_sojourn)

df1 <- generate_dataset(theta1, law_sojourn, n = 50, M = 5)
df2 <- generate_dataset(theta2, law_sojourn, n = 30, M = 5)
df2$id <- df2$id + 50
df3 <- generate_dataset(theta3, law_sojourn, n = 20, M = 5)
df3$id <- df3$id + 80

dataframe <- rbind(df1, df2, df3)

covariates <- data.frame(X1 = rep(c(1, 2), times=c(50, 50)), 
                         X2 = rep(c(1, 2), times=c(80, 20)),
                         X3 = rnorm(100))
weights <- rep(1, 100)

tree <- build_tree(dataframe, covariates, weights, 
                   D, law_sojourn, likelihood_ratio_test, 
                   min_leaf = 5, alpha = 0.05)
plot_tree(tree)

random_forest <- random_forest(dataframe, covariates, weights, 
                               D, law_sojourn, likelihood_ratio_test)

mdi <- MDI_all(random_forest, covariates)
mda <- MDA_all(random_forest, dataframe, covariates, D, weights, law_sojourn)