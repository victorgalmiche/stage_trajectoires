library(doParallel)
library(foreach)

run_simulation <- function(cl, D, n, law_sojourn, R,
                           M = 5, nb_datasets=500){
  
  clusterExport(cl, varlist = c("D", "n", "M", "law_sojourn", "R"),
                envir = environment())
  
  results <- foreach(
    i = 1:nb_datasets,
    .combine = rbind
  ) %dopar% {  
    
    theta <- generate_theta(D, law_sojourn)
    df <- generate_dataset_H0(theta, law_sojourn, n, n, M)
    
    df1 <- subset(df, id<=n)
    df2 <- subset(df, id>n)
    
    t_asymp <- system.time(
      p_asymp <- likelihood_ratio_test(df1, df2, D, law_sojourn = law_sojourn)
    )
    
    t_boot <- system.time(
      p_boot <- parametric_bootstrap(df1, df2, D, law_sojourn = law_sojourn, 
                                     M = M, R = R)
    )
    
    t_perm <- system.time(
      p_perm <- permutation_test(df1, df2, D, law_sojourn = law_sojourn, R = R)
    )
    
    c(t_asymp=unname(t_asymp['elapsed']), 
      t_boot=unname(t_boot['elapsed']), 
      t_perm=unname(t_perm['elapsed']))
  }
  
  list(
    D = D, n = n, law_sojourn = law_sojourn, R = R,
    M = M, nb_datasets = nb_datasets,
    t_asymp = results[, "t_asymp"],
    t_boot = results[, "t_boot"],
    t_perm = results[, "t_perm"]
  )
}

plot_times <- function(sim_result) {
  df_time <- data.frame(
    time   = c(sim_result$t_asymp, sim_result$t_boot, sim_result$t_perm),
    algo = rep(c("Chi^2", "Bootstrap", "Permutation"),
                  each = sim_result$nb_datasets)
  )
  
  boxplot(time ~ algo, data = df_time,
          log = 'y',
          main = "Time Executions Comparison",
          xlab = "", ylab = "Time (seconds, in log-scale)",
          col = c("#66c2a5", "#fc8d62", "#8da0cb"))
}


n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  source('src/semi_markov/synthesis_data_generation.R')
  source('src/semi_markov/mle_estimation.R')
  source('src/two_samples_test.R')
})

res <- run_simulation(cl, D = 10, n = 200, 
                      law_sojourn = 'exponential', R = 100)
plot_times(res)

stopCluster(cl)


