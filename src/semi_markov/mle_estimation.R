# We dispose of a dataframe with columns id, state, time
# id correspond to an individual and are 1 to N w/ N the number of individual
# the rows attached to an id are ordered in the order of the visited states

# We also dispose of a weight vector 
# i-th element is the weight assigned to the individual for which id=i

# CENSORING: the LAST row of each id is right-censored, i.e. we only know
# the sojourn in that state lasted AT LEAST dataframe$time for that row,
# not that it ended exactly there.

### HELPER
# Function to flag the last row per id, useful for P and omega with censoring
last_row_per_id <- function(dataframe) {
  c(dataframe$id[-1] != dataframe$id[-nrow(dataframe)], TRUE)
}


### LOG-LIKELIHOOD
# Part of the log-likelihood depending on alpha
log_likelihood_alpha <- function(dataframe, alpha, weights=NULL) {
  # We compute the initial states of each process
  init_states <- tapply(dataframe$state, dataframe$id, head, 1)
  
  # And then sum the corresponding log(alpha) weighted 
  if (is.null(weights)) {
    w <- rep(1, length(init_states))
  } else {
    w <- weights[as.integer(names(init_states))]
  }
  ll <- sum(w*log(alpha[init_states]))
  if (!is.finite(ll)) warning("Log-likelihood is not finite, check alpha for zero entries")
  ll
}

# Part of the log-likelihood depending on P
log_likelihood_P <- function(dataframe, P, weights=NULL) {
  # For each row that has a successor in the same chain 
  # This is the opposite of last row per id
  same_chain <- !last_row_per_id(dataframe)
  
  # Computing the transitions
  state_i <- dataframe$state[same_chain]
  state_j <- dataframe$state[which(same_chain) + 1]
  
  log_p <- log(P[cbind(state_i, state_j)])
  
  # Summing the corresponding log(P_ij) 
  if (is.null(weights)){
    w <- rep(1, length(log_p))
  } else {
    chain_id <- dataframe$id[same_chain]
    w <- weights[chain_id]
  }
  ll <- sum(w*log_p)
  if (!is.finite(ll)) warning("Log-likelihood is not finite, check P for zero entries")
  ll
}

# Part of the log-likelihood depending on omega
log_likelihood_omega <- function(dataframe, omega, 
                                 weights=NULL, law_sojourn='gamma') {
  # censored flagged the censored sojourn time
  censored <- last_row_per_id(dataframe)
  observed <- !censored 
  
  log_f <- numeric(nrow(dataframe))
  
  # For the observed sojourn times 
  s_obs <- dataframe$state[observed]
  t_obs <- dataframe$time[observed]
  log_f[observed] <- switch(law_sojourn, 
                  gamma = dgamma(t_obs,
                            shape=omega[s_obs, 'shape'],
                            rate=omega[s_obs, 'rate'],
                            log=TRUE),
                  weibull = dweibull(t_obs,
                                shape=omega[s_obs, 'shape'],
                                scale=omega[s_obs, 'scale'],
                                log=TRUE),
                  exponential = dexp(t_obs, 
                                rate=omega[s_obs, 'rate'],
                                log=TRUE)
  )
  
  # For the censored sojourn times
  t_cen <- dataframe$time[censored]
  s_cen <- dataframe$state[censored]
  log_f[censored] <- switch(law_sojourn, 
                       gamma = pgamma(t_cen,
                                      shape=omega[s_cen, 'shape'],
                                      rate=omega[s_cen, 'rate'],
                                      lower.tail=FALSE, log.p=TRUE),
                       weibull = pweibull(t_cen,
                                          shape=omega[s_cen, 'shape'],
                                          scale=omega[s_cen, 'scale'],
                                          lower.tail=FALSE, log.p=TRUE),
                       exponential = pexp(t_cen, 
                                          rate=omega[s_cen, 'rate'],
                                          lower.tail=FALSE, log.p=TRUE)
  )
  
  if (is.null(weights)){
    w <- rep(1, length(log_f))
  } else {
    w <- weights[dataframe$id]
  }
  ll <- sum(w*log_f)
  if (!is.finite(ll)) warning("Log-likelihood is not finite, check omega")
  ll
}

### MLE
# MLE for alpha
mle_alpha <- function(dataframe, D, weights=NULL) {
  init_states <- tapply(dataframe$state, dataframe$id, head, 1) 
  
  if (is.null(weights)) {
    w <- rep(1, length(init_states))
  } else {
    w <- weights[as.integer(names(init_states))]
  }
  counts <- tapply(w, init_states, sum)
  
  alpha <- numeric(D)
  alpha[as.integer(names(counts))] <- counts 
  alpha/sum(alpha) # Renormalize to have sum=1
}

# MLE for P
mle_P <- function(dataframe, D, weights=NULL) {
  # Count transitions
  same_chain <- !last_row_per_id(dataframe)
  
  state_i <- dataframe$state[same_chain]
  state_j <- dataframe$state[which(same_chain) + 1]
  
  
  if (is.null(weights)){
    w <- rep(1, length(state_i))
  } else {
    chain_id <- dataframe$id[same_chain]
    w <- weights[chain_id]
  }
  
  # Encoding in a matrix 
  idx <- (state_j - 1) * D + state_i
  counts <- tapply(w, idx, sum)
  
  P <- matrix(0, nrow=D, ncol=D)
  P[as.integer(names(counts))] <- counts
    
  # Normalize each row by off-diagonal sum (closed-form MLE under constraint)
  row_sums <- rowSums(P)
  
  # To get 1/(D-1) for unvisited states
  P[row_sums==0, ] = 1
  diag(P) <- 0 # Force diagonal to be 0
  row_sums[row_sums == 0] <- D-1  
  
  P / row_sums
}


# MLE for omega - using Nelder-Mead optimization
mle_omega_gamma <- function(dataframe, D, weights=NULL){
  censored <- last_row_per_id(dataframe)
  
  # Creating a matrix object to store the result - initialized w/ ones
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c('shape', 'rate')))
  
  # For each state s
  for (s in 1:D){
    # Extracting the corresponding rows
    keep_s <- dataframe$state==s
    df_s <- dataframe[keep_s,]
    cen_s <- censored[keep_s]
    
    # and weights
    if (is.null(weights)){
      w <- rep(1, nrow(df_s))
    } else {
      w <- weights[df_s$id]
    }
    
    # If less than 2 observations, impossible to estimate
    if (nrow(df_s) < 2) {
      warning(paste('State', s, 'has less than 2 valid observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return(1e10)
      obs_s <- !cen_s
      ll_obs <- if(any(obs_s))
        sum(w[obs_s]*dgamma(df_s$time[obs_s], shape=pars[1], rate=pars[2], 
                            log=TRUE))
      else 0
      ll_cen <- if(any(cen_s)) 
        sum(w[cen_s]*pgamma(df_s$time[cen_s], shape=pars[1], rate=pars[2], 
                            lower.tail=FALSE, log.p=TRUE))
      else 0
      ll <- ll_obs + ll_cen
      if (!is.finite(ll)) return(1e10)
      -ll
    }
    
    # Starting values using method of moments
    mean_x <- mean(df_s$time)
    var_x <- var(df_s$time)
    
    if (is.na(var_x) || var_x < .Machine$double.eps) {
      warning(paste('State', s, 'has near-zero or NA variance'))
      var_x <- 1
    }
    
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
    
    omega[s, 'shape'] <- result$par[1]
    omega[s, 'rate'] <- result$par[2]
  }
  omega
}

mle_omega_weibull <- function(dataframe, D, weights=NULL){
  censored <- last_row_per_id(dataframe)

  # Creating a matrix object to store the result - initialized w/ ones
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c('shape', 'scale')))
  
  # For each state s
  for (s in 1:D){
    # Extracting the corresponding rows
    keep_s <- dataframe$state==s
    df_s <- dataframe[keep_s,]
    cen_s <- censored[keep_s]
    
    # and weights
    if (is.null(weights)){
      w <- rep(1, nrow(df_s))
    } else {
      w <- weights[df_s$id]
    }
    
    # If less than 2 observations, impossible to estimate
    if (nrow(df_s) < 2) {
      warning(paste('State', s, 'has less than 2 valid observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return(1e10)
      obs_s <- !cen_s
      ll_obs <- if (any(obs_s))
        sum(w[obs_s]*dweibull(df_s$time[obs_s], shape=pars[1], scale=pars[2], 
                              log=TRUE))
      else 0
      ll_cen <- if (any(cen_s))
        sum(w[cen_s]*pweibull(df_s$time[cen_s], shape=pars[1], scale=pars[2], 
                              lower.tail=FALSE, log.p=TRUE))
      else 0
      ll <- ll_obs + ll_cen
      if (!is.finite(ll)) return(1e10)
      -ll
    }
    
    # Starting values - no easy form w/ method of moments
    shape_start <- 1
    scale_start <- 1
    
    # Test starting values
    if (!is.finite(nll(c(shape_start, scale_start)))) {
      warning(paste('State', s, 'invalid starting values'))
      next
    }
    
    # Optimize w/ Nelder-Mead
    result <- tryCatch(
      optim(c(shape_start, scale_start), nll, 
            method='Nelder-Mead',
            control=list(maxit=5000)),
      error = function(e) {
        warning(paste('State', s, 'optimization error'))
        list(par=c(shape_start, scale_start), convergence=1)
      }
    )
    
    if (result$convergence != 0)
      warning(paste('State', s, 'did not converge'))
    
    omega[s, 'shape'] <- result$par[1]
    omega[s, 'scale'] <- result$par[2]
  }
  omega
}

mle_omega_exponential <- function(dataframe, D, weights=NULL){
  censored <- last_row_per_id(dataframe)
  
  # Creating a matrix object to store the result - initialized w/ ones
  omega <- matrix(1, nrow=D, ncol=1, dimnames=list(1:D, c('rate')))
  
  for (s in 1:D){
    # Extracting the corresponding rows
    keep_s <- dataframe$state==s
    df_s <- dataframe[keep_s,]
    cen_s <- censored[keep_s]
    
    if (nrow(df_s) < 1){
      warning(paste('State', s, 'has no time observations'))
      next
    }
    
    # and weights
    if (is.null(weights)){
      w <- rep(1, nrow(df_s))
    } else {
      w <- weights[df_s$id]
    }
    
    # MLE for exponential distribution
    omega[s, 'rate'] <- sum(w[!cen_s])/sum(w*df_s$time)
  }
  omega
}

### MAXIMUM LIKELIHOOD ESTIMATION ###
mle_fit <- function(dataframe, D, weights=NULL, law_sojourn='gamma'){
  alpha_hat <- mle_alpha(dataframe, D, weights)
  ll_alpha <- log_likelihood_alpha(dataframe, alpha_hat, weights)
  
  P_hat <- mle_P(dataframe, D, weights)
  ll_P <- log_likelihood_P(dataframe, P_hat, weights)
  
  omega_hat <- switch(law_sojourn, 
                      gamma = mle_omega_gamma(dataframe, D, weights),
                      weibull = mle_omega_weibull(dataframe, D, weights), 
                      exponential = mle_omega_exponential(dataframe, D, weights)
  )
  ll_omega <- log_likelihood_omega(dataframe, omega_hat, weights, law_sojourn)
  
  theta_hat <- list(alpha=alpha_hat, P=P_hat, omega=omega_hat)
  ll <- ll_alpha + ll_P + ll_omega
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}