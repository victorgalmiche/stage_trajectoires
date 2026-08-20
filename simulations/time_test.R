library(doParallel)
library(foreach)

run_simulation <- function(cl, D, n, law_sojourn, R,
                           T_max = 5, nb_datasets=500){
  
  clusterExport(cl, varlist = c("D", "n", "T_max", "law_sojourn", "R"),
                envir = environment())
  
  results <- foreach(
    i = 1:nb_datasets,
    .combine = rbind
  ) %dopar% {  
    
    theta <- generate_theta(D, law_sojourn)
    df <- generate_dataset_H0(theta, law_sojourn, n, n, T_max)
    
    df1 <- subset(df, id<=n)
    df2 <- subset(df, id>n)
    
    t_asymp <- system.time(
      p_asymp <- likelihood_ratio_test(df1, df2, D, law_sojourn = law_sojourn)
    )
    
    t_boot <- system.time(
      p_boot <- parametric_bootstrap(df1, df2, D, law_sojourn = law_sojourn, 
                                     T_max = T_max, R = R)
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
    T_max = T_max, nb_datasets = nb_datasets,
    t_asymp = results[, "t_asymp"],
    t_boot = results[, "t_boot"],
    t_perm = results[, "t_perm"]
  )
}

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  source('src/semi_markov/synthesis_data_generation.R')
  source('src/semi_markov/mle_estimation.R')
  source('src/two_samples_test.R')
})

res <- run_simulation(cl, D = 10, n = 30, 
                      law_sojourn = 'exponential', R = 500, 
                      nb_datasets = 50)

# Varying n, number of trajectories
n_vals <- c(30, 60, 100, 200)
asymp_D4 <- c(0.03694, 0.05298, 0.05596, 0.102)
asymp_D10 <- c(0.049, 0.0595, 0.0773, 0.14406)
boot_D4 <- c(45.26, 74.92, 120.35, 216.74)
boot_D10 <- c(52.17, 90.37, 154.58, 308.10)
perm_D4 <- c(0.99794, 1.34, 1.82, 2.93)
perm_D10 <- c(1.40, 1.589, 2.2382, 3.53)

plot(n_vals, asymp_D4, col='red', pch=1, type='b',
     xlab='Number of trajectories in each group (n)', 
     ylab='Time in seconds (log-scale)', log='y', ylim=c(1e-2, 5e2),
     main='Execution Time Comparison')
points(n_vals, asymp_D10, col='red', pch=2, type='b')
points(n_vals, boot_D4, col='blue', pch=1, type='b')
points(n_vals, boot_D10, col='blue', pch=2, type='b')
points(n_vals, perm_D4, col='green', pch=1, type='b')
points(n_vals, perm_D10, col='green', pch=2, type='b')


legend("center", 
       legend = c("Chi^2", "Permutation", "Parametric Bootstrap", "D=4", "D=10"),
       col = c("red", "green", "blue", 'grey', 'grey'),
       pch = c(NA, NA, NA, 1, 2), 
       lty = c(1, 1, 1, NA, NA))



# Varying R number of permutation/bootstrap sample
R_vals <- c(100, 500, 1000)
boot_D4 <- c(45.26, 217.64, 438.52)
boot_D10 <- c(52.17, 247.09, 494.27)
perm_D4 <- c(0.99794, 4.96, 10.45)
perm_D10 <- c(1.40, 6.75, 13.56)

plot(R_vals, boot_D4, col='blue', pch=1, type='b',
     xlab='R', 
     ylab='Time in seconds (log-scale)', log='y', ylim=c(1, 5e2),
     main='Execution Time Comparison')
points(R_vals, boot_D10, col='blue', pch=2, type='b')
points(R_vals, perm_D4, col='green', pch=1, type='b')
points(R_vals, perm_D10, col='green', pch=2, type='b')

stopCluster(cl)


