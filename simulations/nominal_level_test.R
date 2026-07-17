library(doParallel)
library(foreach)

run_simulation <- function(cl, D, n1, n2, M, nb_datasets=500, 
                           law_sojourn= 'exponential'){
  
  clusterExport(cl, varlist = c("D", "n1", "n2", "M", "law_sojourn"),
                envir = environment())
  
  results <- foreach(
    i = 1:nb_datasets,
    .combine = rbind
    ) %dopar% {  
    
    theta <- generate_theta(D, law_sojourn)
    df <- generate_dataset_H0(theta, law_sojourn, n1, n2, M)
    
    df1 <- subset(df, id<=n1)
    df2 <- subset(df, id>n1)
    
    p_asymp  <- likelihood_ratio_test(df1, df2, D, law_sojourn = law_sojourn)
    p_boot <- parametric_bootstrap(df1, df2, D, law_sojourn = law_sojourn, M=M)
    p_perm <- permutation_test(df1, df2, D, law_sojourn = law_sojourn)
    
    c(p_asymp=p_asymp, p_boot=p_boot, p_perm=p_perm)
  }
  
  list(
    D = D, n1 = n1, n2 = n2, M = M,
    nb_datasets = nb_datasets,
    p_asymp = results[, "p_asymp"],
    p_boot = results[, "p_boot"],
    p_perm = results[, "p_perm"]
  )
}

plot_pvalues <- function(sim_result, title = NULL) {
  if (is.null(title)) {
    title <- sprintf("n1=%d, n2=%d, D=%d", sim_result$n1, sim_result$n2, sim_result$D)
  }
  
  par(mar = c(4, 4, 2, 1))
  plot(ecdf(sim_result$p_asymp),
       main = title,
       xlab = "p-value",
       ylab = "Cumulative probability",
       col  = "red",
       do.points = FALSE)
  lines(ecdf(sim_result$p_perm), col = "green", do.points = FALSE)
  lines(ecdf(sim_result$p_boot), col = "blue", do.points = FALSE)
  abline(a = 0, b = 1, col = "black")
  legend("topleft",
         legend = c("Chi^2", "Permutation", "Parametric Bootstrap"),
         col    = c("red", "green", "blue"),
         lty    = c(1, 1, 1))
}

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  source('src/semi_markov/synthesis_data_generation.R')
  source('src/semi_markov/mle_estimation.R')
  source('src/two_samples_test.R')
})

res <- run_simulation(cl, D = 10, n1 = 60, n2 = 60, M = 5, nb_datasets = 500)
plot_pvalues(res)

stopCluster(cl)
