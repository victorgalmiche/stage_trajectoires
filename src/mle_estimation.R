# We dispose of a dataframe with columns id, state, time
# id correspond to an individual
# the rows attached to an id are ordered in the order of the visited states


### LOG-LIKELIHOOD
# Part of the log-likelihood depending on alpha
log_likelihood_alpha <- function(df, alpha) {
  # We compute the initial states of each process
  init_states <- tapply(df$state, df$id, head, 1)
  
  # And then sum the corresponding log(alpha)
  ll <- sum(log(alpha[init_states]))
  if (!is.finite(ll)) warning("Log-likelihood is not finite, check alpha for zero entries")
  ll
}

# Part of the log-likelihood depending on P
log_likelihood_P <- function(df, P) {
  # For each row that has a successor in the same chain
  same_chain <- c(df$id[-1] == df$id[-nrow(df)], FALSE)
  
  # Computing the transitions
  state_i <- df$state[same_chain]
  state_j <- df$state[which(same_chain) + 1]
  
  # Summing the corresponding log(P_ij) 
  sum(log(P[cbind(state_i, state_j)]))
}

# Part of the log-likelihood depending on omega
log_likelihood_omega <- function(df, omega){
  sum(dgamma(df$time,
             shape=omega[df$state, 'shape'],
             rate=omega[df$state.h, 'rate'],
             log=TRUE))
}

### MLE
# MLE for alpha
mle_alpha <- function(df, D) {
  init_states <- tapply(df$state, df$id, head, 1) 
  counts <- tabulate(init_states, nbins = D) # Number for each state
  counts/sum(counts) # Renormalize to have sum=1
}

# MLE for P
mle_P <- function(df, D) {
  # Count transitions
  same_chain <- c(df$id[-1] == df$id[-nrow(df)], FALSE)
  
  state_i <- df$state[same_chain]
  state_j <- df$state[which(same_chain) + 1]
  
  # Encoding in a matrix 
  counts <- matrix(
    tabulate((state_i - 1) * D + state_j, nbins = D * D),
    nrow = D, byrow=TRUE)
  
  # Force diagonal to 0 (semi-Markov: no self-transitions)
  diag(counts) <- 0
  
  # Normalize each row by off-diagonal sum (closed-form MLE under constraint)
  row_sums <- rowSums(counts)
  row_sums[row_sums == 0] <- Inf  # to get 0 for unvisited states
  counts / row_sums
}


# MLE for omega - using Nelder-Mead optimization
mle_omega_nm <- function(df, D){
  # Creating a matrix object to store the result - initialized w/ ones
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c('shape', 'rate')))
  
  # For each state s
  for (s in 1:D){
    # Extracting the corresponding sojourn times
    times_s <- df$time[df$state==s]
    
    # If less than 2 observations, impossible to estimate
    if (length(times_s) < 2) {
      warning(paste('State', s, 'has less than 2 valid observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return(1e10)
      ll <- sum(dgamma(times_s, shape=pars[1], rate=pars[2], log=TRUE))
      if (!is.finite(ll)) return(1e10)
      -ll
    }
    
    # Starting values using method of moments
    mean_x <- mean(times_s)
    var_x <- var(times_s)
    
    shape_start <- mean_x^2 / var_x
    rate_start <- mean_x / var_x
    
    # Test starting values
    if (!is.finite(nll(c(shape_start, rate_start)))) {
      warning(paste('State', s, 'invalid starting values'))
      next
    }
    
    # Optimize w/ Nelder-Mead
    result <- tryCatch(
      optim(c(shape_start, rate_start), nll, 
            method='Nelder-Mead',
            control=list(maxit=5000)),
      error = function(e) {
        warning(paste('State', s, 'optimization error'))
        list(par=c(shape_start, rate_start), convergence=1)
      }
    )
    
    if (result$convergence != 0)
      warning(paste('State', s, 'did not converge'))
    
    omega[i, 'shape'] <- result$par[1]
    omega[i, 'rate'] <- result$par[2]
  }
  omega
}

mle_omega_bfgs <- function(df, D){
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c("a", "lambda")))
  for (i in 1:D){
    times_i <- df$time[df$state.h==i]
    
    if (length(times_i) < 2) {
      warning(paste('State', i, 'has less than 2 valid observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return(1e10)
      ll <- sum(dgamma(times_i, shape=pars[1], rate=pars[2], log=TRUE))
      if (!is.finite(ll)) return(1e10)
      -ll
    }
    
    # Starting values using method of moments
    mean_x <- mean(times_i)
    var_x <- var(times_i)
    
    # Handle invalid variance
    if (var_x <= 0 || !is.finite(var_x)) {
      var_x <- (mean_x / 4)^2
    }
    
    shape_start <- mean_x^2 / var_x
    rate_start <- mean_x / var_x
    
    # Bound starting values to reasonable range
    shape_start <- pmax(0.01, pmin(shape_start, 500))
    rate_start <- pmax(0.001, pmin(rate_start, 500))
    
    # Test starting values
    if (!is.finite(nll(c(shape_start, rate_start)))) {
      warning(paste('State', i, 'invalid starting values'))
      next
    }
    
    # Optimize w/ Nelder-Mead
    result <- tryCatch(
      optim(c(shape_start, rate_start), nll, 
            method='L-BFGS-B',
            lower = c(0.001, 0.001),
            control=list(maxit=5000)),
      error = function(e) {
        warning(paste('State', i, 'optimization error'))
        list(par=c(shape_start, rate_start), convergence=1)
      }
    )
    
    omega[i, 'a'] <- result$par[1]
    omega[i, 'lambda'] <- result$par[2]
  }
  omega
}


### MAXIMUM LIKELIHOOD ESTIMATION ###
mle_fit <- function(df, D){
  alpha_hat <- mle_alpha(df, D)
  ll_alpha <- log_likelihood_alpha(df, alpha_hat)
  
  P_hat <- mle_P(df, D)
  ll_P <- log_likelihood_P(df, P_hat)
  
  omega_hat <- mle_omega_nm(df, D)
  ll_omega <- log_likelihood_omega(df, omega_hat)
  
  theta_hat <- list(alpha=alpha_hat, P=P_hat, omega=omega_hat)
  ll <- ll_alpha + ll_P + ll_omega
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}