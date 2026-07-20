
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
# M is the number of observed transitions per SMP
generate_dataset <- function(theta, law_sojourn='gamma', n, M) {
  D <- length(theta$alpha) # Number of states
  
  # Matrix describing the succession of states 
  states <- matrix(0L, nrow = n, ncol = M) 
  
  # The first states are sampled using alpha parameter
  states[, 1] <- sample.int(D, size=n, replace=TRUE, prob=theta$alpha)
  
  for (j in 2:M) {
    # Successions of states 
    current_states <- states[, j - 1] 
    
    probs <- theta$P[current_states,]
    
    # Using an inverse-CDF approach
    cum_probs <- t(apply(probs, 1, cumsum))
    u <- runif(n)
    states[, j] <- rowSums(cum_probs < u) + 1L
  }
  
  # Sojourn times
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
  
  # Calculating sojourn times in each state
  sojournTimes <- apply(states[, ], c(1, 2), fun)
  
  # Build a dataframe storing the precedent information
  id <- rep(1:n, each = M)
  state <- as.vector(t(states))
  time <- as.vector(t(sojournTimes))
  
  data.frame(id = id, state=state, time=time)
}


### Generating a dataset under H0 
# n1 and n2 are the number of processes in each group
generate_dataset_H0 <- function(theta, law_sojourn='gamma', n1, n2, M) {
  df1 <- generate_dataset(theta, law_sojourn, n1, M)
  df2 <- generate_dataset(theta, law_sojourn, n2, M)
  df2$id <- df2$id + n1
  
  rbind(df1, df2)
}


