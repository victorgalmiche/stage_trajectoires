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

### DATAFRAME SEMI-MARKOV ###

# Assigning the id (ie the row number in individus)
emplois$id <- match(emplois$IDENT, individus$IDENT)
non_emplois$id <- match(non_emplois$IDENT, individus$IDENT)

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



### COVARIATES TABLE ###
# Creating a new column w/ PHD
mapping_PHD <- rep(1, 18)
mapping_PHD[2:5] <- 2
mapping_PHD[6:12] <- 3
mapping_PHD[13:18] <- 4
individus$PHD_NEW <- mapping_PHD[as.integer(substr(individus$PHD, 1, 2))]


# Columns chosen for the covariates table
# cols_quali <- c('Q1', 'Q16', 'Q31', 'OS1', 'OS3_1', 'OS3_2', 'OS3_3',
#                'ETR1', 'PER1', 'SITPERE', 'SITMERE', 'CA13', 'CA22')
# cols_quanti <- c('PHD_NEW', 'AGE10')

cols_quali <- c('Q1', 'Q16', 'Q31', 'OS1', 'OS3_1', 'OS3_2', 'OS3_3',
                'ETR1', 'PER1', 'CA13', 'CA22')
cols_quanti <- c('PHD_NEW')

# Creating of the covariates table
covariates <- individus[, cols_quanti]
covariates[cols_quali] <- lapply(individus[cols_quali], as.factor)


### WEIGHTS and others ###
weights <- individus$pondef

# Number of states
D <- 10

