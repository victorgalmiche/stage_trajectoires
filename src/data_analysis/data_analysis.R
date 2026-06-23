library(haven)
library(aws.s3)
library(dplyr)
library(tidyr)

BUCKET <- 'victorgalmiche'
FOLDER <- 'stage-trajectoires/lil-1439/lil-1439.dta/Stata'
FILE_INDIV <- paste(FOLDER, 'g107individusvf.dta', sep='/')
FILE_EMPLOI <- paste(FOLDER, 'g107seqentrvf.dta', sep='/')
FILE_NON_EMPLOI <- paste(FOLDER, 'g107nonemplvf.dta', sep='/')


individus <- aws.s3::s3read_using(
  FUN = haven::read_dta,
  object = FILE_INDIV,
  bucket = BUCKET,
  opts = list("region"="")
)

emplois <- aws.s3::s3read_using(
  FUN = haven::read_dta,
  object = FILE_EMPLOI,
  bucket = BUCKET,
  opts = list("region"="")
)

non_emplois <- aws.s3::s3read_using(
  FUN = haven::read_dta,
  object = FILE_NON_EMPLOI,
  bucket = BUCKET,
  opts = list("region"="")
)

# Assigning the id (ie the row number in individus)
emplois$id <- match(emplois$IDENT, individus$IDENT)
non_emplois$id <- match(non_emplois$IDENT, individus$IDENT)

# Extracting the states - We use TYPESEQ and codify into integers
mapping <- c(
  'int'=2, 'asc'=1, 'afa'=1, 'sco'=1, 'slo'=1, 'vac'=6,
  'chc'=3, 'chl'=3, 'inc'=6, 'inl'=6, 'foc'=5, 'fol'=5, 'rep'=4
)
emplois$state <- mapping[emplois$TYPESEQ]
non_emplois$state <- mapping[non_emplois$TYPESEQ]

# And the sojourn time
emplois$time <- emplois$DUREE
non_emplois$time <- non_emplois$DUREE

# Merging the two dataframes to regroup emplois and non_emplois
df_merged <- emplois |> select(id, NSEQ, state, time) |>
  bind_rows(non_emplois |> select(id, NSEQ, state, time)) |>
  arrange(id, NSEQ) |> 
  select(id, state, time)


# Now, creating the trajectory dataframe used for analysis
group <- cumsum(c(TRUE, diff(as.numeric(interaction(df_merged$id, df_merged$state)))!=0))
dataframe <- aggregate(time ~ group + id + state, data = cbind(df_merged, group), sum)
dataframe <- dataframe[order(dataframe$group), c("id", "state", "time")]


# Descriptive statistics
state_labels <- c('Employment', 'Interim', 'Job Search', 
                  'School', 'Training', 'Others')
df_sum <- dataframe %>%
  group_by(state) %>%
  summarise(total_time = sum(time, na.rm = TRUE))

barplot(df_sum$total_time,
        # log = 'y',
        names.arg = state_labels,
        xlab = "State",
        ylab = "Time in months",
        # ylim = c(1e+03, 1e+06),
        main = "Cumulative time spent in each state"
)

barplot(table(dataframe$state), 
     names.arg = state_labels, 
     xlab = 'State', 
     ylab = 'Number of observations', 
     ylim = c(0, 2e+04),
     main = 'Number of observations for each state')

hist(dataframe$time,
     xlab = 'Time in months', 
     ylab = 'Number of observations', 
     main = 'Histogram of sojourn times')


employment <- dataframe[dataframe$state==1, ]
library(MASS)
fit_exp <- fitdistr(employment$time, "exponential")
fit_gamma <- fitdistr(employment$time, "gamma")
fit_weibull <- fitdistr(employment$time, "weibull")

x_vals <- seq(0, max(employment$time), length.out = 100)
hist(employment$time, 
     breaks = x_vals,
     freq = FALSE,
     xlab = 'Time in months', 
     ylab = 'Density', 
     main = 'Histogram of sojourn times in employment')
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
       legend = c("Gamma", "Weibull", "Exponential", paste0("n = ", nrow(employment))),
       col    = c("red", "green", "blue", NA),
       lty    = c(1, 1, 1, NA))




# To obtain the trajectories
trajectories <- dataframe %>%
  group_by(id) %>%
  summarise(trajectoire = list(rep(state, times = time)), .groups = "drop") %>%
  tidyr::unnest_wider(trajectoire, names_sep = "_t")


# Visualize the trajectories
library(TraMineR)
seq <- seqdef(trajectories, 2:80)
par(mfrow = c(2, 2))
seqiplot(seq, with.legend=FALSE, border=NA)

