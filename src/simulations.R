library(doParallel)
library(foreach)

source('src/semi_markov/synthesis_data_generation.R')
source('src/semi_markov/mle_estimation.R')
source('src/two_samples_test.R')

law_sojourn <- 'exponential'
theta1 <- list(alpha = c(1, 0), 
               P = matrix(c(0, 1, 1, 0), nrow=2), 
               omega = cbind(rate=rep(1, times=2)))

nb_sets <- 100
step <- 20

p_asymp <- matrix(0, ncol=step+1, nrow=nb_sets)
p_perm <- matrix(0, ncol=step+1, nrow=nb_sets)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

results <- foreach(
  i = 1:nb_sets,
  .combine = rbind,
  .export = c("generate_dataset", "likelihood_ratio_test", "permutation_test")
) %dopar% {
  df1 <- generate_dataset(theta1, law_sojourn, n=250, M=2)
  
  p_asymp <- numeric(step+1)
  p_perm <- numeric(step+1)
  
  for (k in 0:step){
    theta2 <- list(alpha = c(1, 0), 
                   P = matrix(c(0, 1, 1, 0), nrow=2), 
                   omega = cbind(rate=rep(1+k/step, times=2)))
    df2 <- generate_dataset(theta2, law_sojourn, n=250, M=2)
    df2$id <- df2$id + 250
    
    p_asymp[k+1] <- likelihood_ratio_test(df1, df2, D=2, law_sojourn=law_sojourn)
    p_perm[k+1] <- permutation_test(df1, df2, D=2, law_sojourn=law_sojourn)
    
  }
  c(p_asymp, p_perm)
}
stopCluster(cl)

reject_asymp <- colMeans(results[, 1:step+1] < 0.05)
reject_perm  <- colMeans(results[, step+2:2*(step+1)] < 0.05)

rate_vals <- 1+0:step/step

plot(rate_vals, reject_asymp, type = "l", col = "blue", lwd = 2,
     xlab = "rate", ylab = "Proportion of p<0.05",
     main = "Exponential: rate",
     ylim = c(0, 1))
lines(rate_vals, reject_perm, col = "red", lwd = 2)
legend("topright", legend = c("Chi^2", "Permutation"),
       col = c("blue", "red"), lwd = 2)
