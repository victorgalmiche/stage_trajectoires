library(doParallel)
library(foreach)

run_simulation <- function(cl, D, n, M, nb_datasets=500, nb_steps=20,
                           var_parameter, law_sojourn= 'exponential'){
  
  clusterExport(cl, varlist = c("D", "n", "M", "nb_steps", 
                                "law_sojourn", "var_parameter"),
                envir = environment())
  
  results <- foreach(i = 1:nb_datasets) %dopar% {  
    
    # Theta_0 and the associated trajectories
    theta0 <- generate_theta(D, law_sojourn)
    df0 <- generate_dataset(theta0, law_sojourn, n, M)
    
    # Theta_1
    theta1 <- theta0
    if (var_parameter=='alpha'){
      theta1$alpha <- generate_alpha(D)
    } 
    if (var_parameter=='P'){
      theta1$P <- generate_P(D)
    }
    if (var_parameter=='omega'){
      theta1$omega <- generate_omega(D, law_sojourn)
    }
    
    p_asymp <- numeric(nb_steps+1)
    p_boot <- numeric(nb_steps+1)
    p_perm <- numeric(nb_steps+1)
    
    for (k in 0:nb_steps){
      t <- k/nb_steps
      theta <- list(
        alpha = (1-t) * theta0$alpha + t * theta1$alpha,
        P     = (1-t) * theta0$P     + t * theta1$P,
        omega = (1-t) * theta0$omega + t * theta1$omega
      )
      
      df1 <- generate_dataset(theta, law_sojourn, n, M)
      df1$id <- df1$id + n # To avoid collisions w/ df0
      
      p_asymp[k+1] <- likelihood_ratio_test(df0, df1, D, law_sojourn=law_sojourn)
      p_boot[k+1] <- parametric_bootstrap(df0, df1, D, law_sojourn=law_sojourn, M=M)
      p_perm[k+1] <- permutation_test(df0, df1, D, law_sojourn=law_sojourn)
      
    }
    list(p_asymp=p_asymp, p_boot=p_boot, p_perm=p_perm)
  }
  
  list(
    D = D, n = n, M = M,
    nb_datasets = nb_datasets,
    nb_steps = nb_steps,
    var_parameter = var_parameter,
    p_asymp = do.call(rbind, lapply(results, `[[`, "p_asymp")),
    p_boot = do.call(rbind, lapply(results, `[[`, "p_boot")),
    p_perm = do.call(rbind, lapply(results, `[[`, "p_perm"))
  )
}

plot_power <- function(sim_result, title = NULL) {
  if (is.null(title)) {
    title <- sprintf("Varying Parameter: %s", sim_result$var_parameter)
  }
  
  vals <- 0:sim_result$nb_steps/sim_result$nb_steps
  reject_asymp <- colMeans(sim_result$p_asymp < 0.05)
  reject_boot <- colMeans(sim_result$p_boot < 0.05)
  reject_perm  <- colMeans(sim_result$p_perm < 0.05)
  
  plot(vals, reject_asymp, col = "red", pch=1, cex=.8,
       xlab = "t", ylab = "Proportion of p<0.05",
       main = title, 
       ylim = c(0, 1))
  points(vals, reject_boot, col="blue", pch=0, cex=.8)
  points(vals, reject_perm, col = "green", pch=2, cex=.8)
  abline(a = 0.05, b = 0, col = "grey", lty=2, lwd=2)
  
  legend("bottomright",
         legend = c("Chi^2", "Permutation", "Parametric Bootstrap", "Level 0.05"),
         col = c("red", "green", "blue", "grey"),
         pch = c(1, 0, 2, NA), 
         lty = c(NA, NA, NA, 2))
  
}

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  source('src/semi_markov/synthesis_data_generation.R')
  source('src/semi_markov/mle_estimation.R')
  source('src/two_samples_test.R')
})

res <- run_simulation(cl, D = 4, n = 30, M = 5,
                      var_parameter = 'omega')
res_alpha <- run_simulation(cl, D = 4, n = 30, M = 5, 
                            var_parameter = 'alpha')
res_P <- run_simulation(cl, D = 4, n = 30, M = 5, 
                        var_parameter = 'P')
plot_power(res_alpha)
stopCluster(cl)