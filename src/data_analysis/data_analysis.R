source('src/data_analysis/extract_data.R')
library(dplyr)
library(tidyr)
library(MASS)
library(TraMineR)

# Assigning the id (ie the row number in individus)
emplois$id <- match(emplois$IDENT, individus$IDENT)
non_emplois$id <- match(non_emplois$IDENT, individus$IDENT)

# Extracting the states
emplois$state <- as.integer(factor(emplois$CONTRAT_EMB))

mapping <- c(
  '05'=6, '06'=6, '11'=6, '12'=6, # job search
  '07'=7, '08'=7, '13'=7, '14'=7, # inactivity
  '09'=8, '10'=8, '15'=8, '16'=8, # training
  '17'=9, '18'=9, # school
  '21'=10 # holidays
)
non_emplois$state <- as.integer(mapping[non_emplois$CAL])

# And the sojourn time
emplois$time <- emplois$DUREE
non_emplois$time <- non_emplois$DUREE

# Merging the two dataframes to regroup emplois and non_emplois
df_merged <- emplois |> dplyr::select(id, NSEQ, state, time) |>
  bind_rows(non_emplois |> dplyr::select(id, NSEQ, state, time)) |>
  arrange(id, NSEQ) |> 
  dplyr::select(id, state, time)


# Now, creating the trajectory dataframe used for analysis
group <- cumsum(c(TRUE, diff(as.numeric(interaction(df_merged$id, df_merged$state)))!=0))
dataframe <- aggregate(time ~ group + id + state, data = cbind(df_merged, group), sum)
dataframe <- dataframe[order(dataframe$group), c("id", "state", "time")]


# Descriptive statistics
state_labels <- c('Non salarie', 'CDI', 'Contrat aide', 'CDD', 'Interim',
                  'Job search', 'Inactivity', 'Training', 'School', 'Holidays')
df_sum <- dataframe %>%
  group_by(state) %>%
  summarise(total_time = sum(time, na.rm = TRUE))

barplot(df_sum$total_time,
        # log = 'y',
        names.arg = state_labels,
        xlab = "State",
        ylab = "Time in months",
        # ylim = c(1e+03, 1e+06),
        main = "Cumulative time spent in each state",
        cex.names = 0.7
)

barplot(table(dataframe$state), 
     names.arg = state_labels, 
     xlab = 'State', 
     ylab = 'Number of observations', 
     main = 'Number of observations for each state',
     cex.names = 0.7)

hist(dataframe$time,
     xlab = 'Time in months', 
     ylab = 'Number of observations', 
     main = 'Histogram of sojourn times')


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




# To obtain the trajectories
trajectories <- dataframe %>%
  group_by(id) %>%
  summarise(trajectoire = list(rep(state, times = time)), .groups = "drop") %>%
  tidyr::unnest_wider(trajectoire, names_sep = "_t")


# Visualize the trajectories
seq <- seqdef(trajectories, 2:80)
par(mfrow = c(2, 2))
seqiplot(seq, with.legend=FALSE, border=NA)

