source('src/data/extract_data.R')
source('src/semi_markov/mle_estimation.R')
library(dplyr)
library(tidyr)
library(MASS)
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
labels <- c("Self-employed", "Permanent Contract", "Subsidized Contract", 
            "Fixed-term Contract", "Temporary Work", "Job Search", 
            "Inactivity", "Training", "Return to School", "Holidays")
seq <- seqdef(trajectories, 2:99, weights=weights, labels=labels)

# Printing 10 randomly selected trajectories w/ legend
par(mfrow = c(1, 2))
seqiplot(seq, idxs = sample(8882, size=10), 
         with.legend = FALSE, border = NA, weighted=FALSE, 
         main = "10 randomly selected trajectories",
         xtstep = 12, yaxis=FALSE)
seqlegend(seq)

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
P_hat <- mle_P(dataframe, D, weights)

# Mean time in each state
seqmtplot(seq, with.legend=FALSE, 
          main = "Mean time spent in each state", 
          ylab = "Time in months")

# Modal state sequence
seqmsplot(seq, with.legend=FALSE, border = NA,
          ylab = 'State frequency', 
          xtstep = 12)



# Sum of sojourn times in each state
df_sum <- dataframe %>%
  group_by(state) %>%
  summarise(total_time = sum(time, na.rm = TRUE))

# Barplot of the time in each state
barplot(df_sum$total_time,
        # log = 'y',
        names.arg = state_labels,
        xlab = "State",
        ylab = "Time in months",
        # ylim = c(1e+03, 1e+06),
        main = "Cumulative time spent in each state",
        cex.names = 0.7
)

# Barplot of the observation in each state
barplot(table(dataframe$state), 
     names.arg = state_labels, 
     xlab = 'State', 
     ylab = 'Number of observations', 
     main = 'Number of observations for each state',
     cex.names = 0.7)

# Histogram of sojourn times
hist(dataframe$time,
     xlab = 'Time in months', 
     ylab = 'Number of observations', 
     main = 'Histogram of sojourn times')


# Histogram of sojourn times for a specific state - here CDD
cdd <- dataframe[dataframe$state==4, ]
fit_exp <- fitdistr(cdd$time, "exponential")
fit_gamma <- fitdistr(cdd$time, "gamma")
fit_weibull <- fitdistr(cdd$time, "weibull")

x_vals <- seq(0, max(cdd$time), length.out = 100)
hist(cdd$time, 
     breaks = x_vals,
     freq = FALSE,
     xlab = 'Time in months', 
     ylab = 'Density', 
     main = 'Histogram of sojourn times in CDD')
lines(x_vals, 
      dgamma(x_vals, shape=fit_gamma$estimate[['shape']], rate=fit_gamma$estimate[['rate']]), 
      col = "red")
lines(x_vals, 
      dweibull(x_vals, shape=fit_weibull$estimate[['shape']], scale=fit_weibull$estimate[['scale']]),
      col = "green")
lines(x_vals, 
      dexp(x_vals, rate=fit_exp$estimate),
      col="blue")
legend("top",
       legend = c("Gamma", "Weibull", "Exponential", paste0("n = ", nrow(cdd))),
       col    = c("red", "green", "blue", NA),
       lty    = c(1, 1, 1, NA))

