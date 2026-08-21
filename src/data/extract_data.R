library(haven)
library(aws.s3)
library(dplyr)

# Pointing the files of interests
BUCKET <- 'victorgalmiche'
FOLDER <- 'stage-trajectoires/lil-1439/lil-1439.dta/Stata'
FILE_INDIV <- paste(FOLDER, 'g107individusvf.dta', sep='/')
FILE_EMPLOI <- paste(FOLDER, 'g107seqentrvf.dta', sep='/')
FILE_NON_EMPLOI <- paste(FOLDER, 'g107nonemplvf.dta', sep='/')

# Charging the 3 tables: 
# individus, sequences d'emplois, sequences de non emploi

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


### COVARIATES TABLE ###
chosen_vars <- c('PHD', 'CA13', 'CA20_17', 'CA22', 'CA24', 'ETR1', 'OS1', 
                 'OS3', 'PER1', 'Q1', 'Q16', 'Q2', 'Q31', 'Q53_13', 'Q9',
                 'SITMERE', 'SITPERE', 'ZUS')
individus_clean <- individus[, c('IDENT', chosen_vars)]

individus_clean <- individus_clean[individus_clean$CA20_17 != 4, ]
individus_clean <- individus_clean[individus_clean$CA22 != 5, ]
individus_clean <- individus_clean[individus_clean$CA24 != 3, ]
individus_clean <- individus_clean[individus_clean$Q31 != 3, ]


# PHD
mapping_PHD <- rep(1, times=18)
mapping_PHD[2:5] <- 2
mapping_PHD[6:12] <- 3
mapping_PHD[13:18] <- 4
covariates <- data.frame(
  PHD=mapping_PHD[as.integer(substr(individus_clean$PHD, 1, 2))])

# CA13
covariates$CA13 <- factor(individus_clean$CA13, 
                          levels=c(1,2), 
                          labels=c('Yes', 'No'))

# CA20_17
mapping_CA20 <- c(1, 2, 2)
covariates$CA20_17 <- factor(mapping_CA20[as.integer(individus_clean$CA20_17)], 
                             levels=c(1,2), 
                             labels=c('Yes', 'No'))

# CA22
covariates$CA22 <- as.integer(individus_clean$CA22)

# CA24
covariates$CA24 <- factor(individus_clean$CA24, 
                          levels=c(1, 2), 
                          labels=c('Yes', 'No'))

# ETR1
mapping_ETR1 <- c(1, 1, 2)
covariates$ETR1 <- factor(mapping_ETR1[as.integer(individus_clean$ETR1)], 
                          levels=c(1,2),
                          labels=c('Yes','No'))

# OS1
covariates$OS1 <- factor(individus_clean$OS1, 
                         levels=c(1,2),
                         labels=c('Yes','No'))

# OS3 
mapping_OS3 <- c("4"=0, "1"=1, "2"=1, "3"=1, 
                 "12"=2, "21"=2, "13"=2, "31"=2, "23"=2, "32"=2, 
                 "123"=3, "132"=3, "213"=3, "231"=3, "312"=3, "321"=3)
covariates$OS3 <- unname(mapping_OS3[individus_clean$OS3])

# PER1
covariates$PER1 <- factor(individus_clean$PER1, 
                          levels=c(1,2), 
                          labels=c('Yes','No'))

# Q1
covariates$Q1 <- factor(individus_clean$Q1, 
                        levels=c(1, 2), 
                        labels=c('Man', 'Woman'))

# Q16
covariates$Q16 <- factor(individus_clean$Q16, 
                         levels=c(1,2), 
                         labels=c('Yes','No'))

# Q2
covariates$Q2 <- as.integer(individus_clean$Q2)

# Q31
covariates$Q31 <- factor(individus_clean$Q31, 
                         levels=c(1,2), 
                         labels=c('Yes', 'No'))

# Q53_13
mapping_Q53 <- c(1, 1, 2)
covariates$Q53_13 <- factor(mapping_Q53[as.integer(individus_clean$Q53_13)], 
                          levels=c(1,2),
                          labels=c('Yes','No'))

# Q9
covariates$Q9 <- factor(individus_clean$Q9, 
                        levels=c(1,2), 
                        labels=c('Yes', 'No'))

# SITMERE
covariates$SITMERE <- as.factor(individus_clean$SITMERE)

# SITPERE
covariates$SITPERE <- as.factor(individus_clean$SITPERE)

# ZUS
covariates$ZUS <- factor(individus_clean$ZUS, 
                         levels=c(1,2), 
                         labels=c('Yes', 'No'))



### DATAFRAME SEMI-MARKOV ###

# Assigning the id (ie the row number in individus)
emplois$id <- match(emplois$IDENT, individus_clean$IDENT)
non_emplois$id <- match(non_emplois$IDENT, individus_clean$IDENT)

# Extracting the states
emplois$state <- as.integer(factor(emplois$CONTRAT_EMB))
mapping_CAL2state <- c(
  '05'=6, '06'=6, '11'=6, '12'=6, # job search
  '07'=7, '08'=7, '13'=7, '14'=7, # inactivity
  '09'=8, '10'=8, '15'=8, '16'=8, # training
  '17'=9, '18'=9, # school
  '21'=10 # holidays
)
non_emplois$state <- as.integer(mapping_CAL2state[non_emplois$CAL])

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

# Number of states
D <- 10

