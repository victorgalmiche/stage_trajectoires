source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')

# Classical LR test (also called chi2)
likelihood_ratio_test <- function(df1, df2, D, weights=NULL, law_sojourn='gamma'){
  # Combining the two groups
  dataframe <- rbind(df1, df2)
  
  # Computing the different estimation
  global_est <- mle_fit(dataframe, D, weights, law_sojourn)
  est1 <- mle_fit(df1, D, weights, law_sojourn)
  est2 <- mle_fit(df2, D, weights, law_sojourn)
  
  # And the lambda statistic
  lambda <- 2*(est1$log_likelihood + 
                 est2$log_likelihood - 
                 global_est$log_likelihood)
  
  # if (!is.finite(lambda) || lambda < 0) return(NA)
  
  # The number of degrees of freedom depends on D
  dof <- switch(law_sojourn, 
                gamma=D^2+D-1,
                weibull=D^2+D-1,
                exponential=D^2-1)
  p_asymp <- 1-pchisq(lambda, df=dof)
  p_asymp
}

# Parametric bootstrap test 
# !!! Not working currently because n1, n2 and M are not defined in generate_dataset_H0
parametric_bootstrap <- function(df1, df2, D, weights=NULL, law_sojourn='gamma', R=100) {
  dataframe <- rbind(df1, df2)
  
  global_est <- mle_fit(dataframe, D, weights, law_sojourn)
  est1 <- mle_fit(df1, D, weights, law_sojourn)
  est2 <- mle_fit(df2, D, weights, law_sojourn)
  
  Tl <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  theta_hat <- global_est$estimator
  
  T_star <- numeric(R)
  for (r in 1:R){
    df_bootstrap <- generate_dataset_H0(theta_hat, law_sojourn, n1, n2, M)
    df1_bootstrap <- subset(df_bootstrap, id<=n1)
    df2_bootstrap <- subset(df_bootstrap, id>n1)
    
    global_est <- mle_fit(df_bootstrap, D, weights, law_sojourn)
    est1 <- mle_fit(df1_bootstrap, D, weights, law_sojourn)
    est2 <- mle_fit(df2_bootstrap, D, weights, law_sojourn)
    
    T_star[r] <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  }
  p_boot <- mean(T_star<=Tl)
  p_boot
}

# Permutation test
permutation_test <- function(df1, df2, D, weights=NULL, law_sojourn='gamma', R=100) {
  dataframe <- rbind(df1, df2)
  n1 <- length(unique(df1$id))
  n2 <- length(unique(df2$id))
  
  global_est <- mle_fit(dataframe, D, weights, law_sojourn)
  est1 <- mle_fit(df1, D, weights, law_sojourn)
  est2 <- mle_fit(df2, D, weights, law_sojourn)
  
  Tl <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  T_star <- numeric(R)
  
  for (r in 1:R) {
    sample1_id   <- sample(unique(dataframe$id), n1)
    df1_permuted <- subset(dataframe, id %in% sample1_id)
    df2_permuted <- subset(dataframe, !(id %in% sample1_id))
    
    est1 <- mle_fit(df1_permuted, D, weights, law_sojourn)
    est2 <- mle_fit(df2_permuted, D, weights, law_sojourn)
    
    T_star[r] <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  }
  
  mean(T_star <= Tl, na.rm=TRUE)
}
