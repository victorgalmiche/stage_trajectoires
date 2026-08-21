source('src/data/extract_data.R')
source('src/semi_markov/mle_estimation.R')
library(dplyr)
library(tidyr)
library(TraMineR)


# Function to transform a dataframe of id, state, time into trajectories
df_to_traj <- function(dataframe){
  # To obtain the trajectories
  trajectories <- dataframe %>%
    group_by(id) %>%
    summarise(t = list(rep(state, times = time)), .groups = "drop") %>%
    tidyr::unnest_wider(t, names_sep = "_")
  trajectories
}

# Creating the dataframe of trajectories
trajectories <- df_to_traj(dataframe)

# Defining the seq object for TraMineR functions
labels <- c("Self-employed", "Permanent Contract", "Subsidized Contract", 
            "Fixed-term Contract", "Temporary Work", "Job Search", 
            "Inactivity", "Training", "Return to School", "Holidays")
seq <- seqdef(trajectories, 2:99, labels=labels)

# Printing 10 randomly selected trajectories w/ legend
par(mfrow = c(1, 2))
seqiplot(seq, idxs = sample(8882, size=10), 
         with.legend = FALSE, border = NA, weighted=FALSE, 
         main = "10 randomly selected trajectories",
         xtstep = 12, yaxis=FALSE)
seqlegend(seq)

# Reinitializing the plot window
dev.off()

# Legend
seqlegend(seq)

# 10 randomly selected trajectories
seqiplot(seq, idxs = sample(8882, size=10), 
         with.legend = FALSE, border = NA, weighted=FALSE, 
         main = "10 randomly selected trajectories",
         xtstep = 12, yaxis=FALSE)

# State distribution over time
seqdplot(seq, with.legend = FALSE, border = NA, 
         main = "State Distribution", 
         xtstep = 12)

# Empirical transition matrix
P_hat <- mle_P(dataframe, D)

# Mean time in each state
seqmtplot(seq, with.legend=FALSE, weighted = FALSE,
          main = "Mean time spent in each state", 
          ylab = "Time in months")

# Modal state sequence
seqmsplot(seq, with.legend=FALSE, border = NA,
          ylab = 'State frequency', 
          xtstep = 12)

# Barplot of the observation in each state
barplot(table(dataframe$state), 
     names.arg = 1:10, 
     col = cpal(seq),
     xlab = 'State',
     ylab = 'Number of observations', 
     main = 'Number of observations for each state')



# Fitting of sojourn times for a specific state - Job Search here
job_search <- dataframe[dataframe$state==6, ]

# Computing MLE for sojourn times
omega_exp <- mle_omega_exponential(job_search, D)
omega_gamma <- mle_omega_gamma(job_search, D)
omega_weibull <- mle_omega_weibull(job_search, D)


hist(job_search$time, freq=FALSE, breaks=100,
     main = "Sojourn time in Job Search (state 6)", 
     xlab = "Time in months", xlim = c(0, 50), ylim = c(0, 0.3))
curve(dexp(x, rate = omega_exp[6, 'rate']), 
      add = TRUE, col = "blue", lwd=2)
curve(dgamma(x, shape = omega_gamma[6, 'shape'], 
             rate = omega_gamma[6, 'rate']), 
      add = TRUE, col = "red", lwd=2)
curve(dweibull(x, shape = omega_weibull[6, 'shape'], 
               scale = omega_weibull[6, 'scale']), 
      add = TRUE, col = "green", lwd=2)
legend("topright",
       legend = c("Gamma", "Weibull", "Exponential"),
       col    = c("red", "green", "blue"),
       lty    = c(1, 1, 1))

