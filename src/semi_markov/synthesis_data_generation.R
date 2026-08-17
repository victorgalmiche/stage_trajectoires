
### Generate parameters defining an SMP
# D is the number of states
# law_sojourn can be Gamma, Weibull or Exponential
generate_alpha <- function(D){
  alpha <- runif(D)
  alpha/sum(alpha) # sum = 1
}

generate_P <- function(D){
  P <- matrix(runif(D^2), nrow=D)
  diag(P) <- 0 # diagonal coefficients are 0
  P/rowSums(P) # each row sum to one 
}

generate_omega <- function(D, law_sojourn) {
  switch(law_sojourn,
         gamma = {
           a <- rexp(D, rate = 0.5) # shape
           lambda <- rexp(D, rate = 3) # rate=1/scale
           cbind(shape = a, rate = lambda)
         },
         weibull = {
           eta <- rexp(D, rate = 0.8) # shape
           beta <- rexp(D, rate = 0.15) # scale
           cbind(shape = eta, scale = beta)
         },
         exponential = {
           lambda <- rexp(D, rate = 1) # rate
           cbind(rate = lambda)
         },
         stop(paste("Unknown sojourn law:", law_sojourn))  # input validation
         )
}


generate_theta <- function(D, law_sojourn='gamma'){
  # Initial probabilities 
  alpha <- generate_alpha(D)
  
  # Transition matrix
  P <- generate_P(D)
  
  # Sojourn times
  omega <- generate_omega(D, law_sojourn)

  return (list(alpha=alpha, P=P, omega=omega))
}


### Generating a dataset of multiple SMP
# theta is the parameter of the SMPs
# n is the number of SMP
# T_max is the total time of observation
generate_dataset <- function(theta, law_sojourn='gamma', n, T_max) {
  D <- length(theta$alpha) # Number of states
  
  # Sojourn times generating function
  fun <- switch(law_sojourn,
                gamma = function(s) {
                  rgamma(1, shape=theta$omega[s, 'shape'], 
                         rate=theta$omega[s, 'rate'])
                },
                weibull = function(s) {
                  rweibull(1, shape=theta$omega[s, 'shape'],
                           scale=theta$omega[s, 'scale'])
                },
                exponential = function(s) {
                  rexp(1, rate=theta$omega[s, 'rate'])
                },
                stop(paste("Unknown sojourn law:", law_sojourn))  # input validation
  )
  
  trajectories <- vector("list", n)
  for (i in 1:n) {
    current_time <- 0
    current_state <- sample.int(D, size=1, replace=TRUE, prob=theta$alpha)
    
    trajectory <- list()
    while (current_time < T_max) {
      sojourn <- fun(current_state)
      observed_sojourn <- min(sojourn, T_max-current_time)
      trajectory[[length(trajectory)+1]] <- data.frame(
        id = i, 
        state = current_state, 
        time = observed_sojourn
      )
      
      probs <- theta$P[current_state, ]
      current_state <- sample.int(D, size=1, replace=TRUE, prob=probs)
      current_time <- current_time + sojourn
    }
    trajectories[[i]] <- do.call(rbind, trajectory)
  }
  
  data.frame(do.call(rbind, trajectories), row.names=NULL)
}


### Generating a dataset under H0 
# n1 and n2 are the number of processes in each group
generate_dataset_H0 <- function(theta, law_sojourn='gamma', n1, n2, T_max) {
  df1 <- generate_dataset(theta, law_sojourn, n1, T_max)
  df2 <- generate_dataset(theta, law_sojourn, n2, T_max)
  df2$id <- df2$id + n1
  
  rbind(df1, df2)
}


