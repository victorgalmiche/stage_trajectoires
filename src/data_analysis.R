library(haven)
library(aws.s3)
library(dplyr)

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

# Assigning the id (ie the row number in individus)
emplois$id <- match(emplois$IDENT, individus$IDENT) 
# Extracting the states (the number starting CONTRAT_EMB)
emplois$state <- as.integer(
  substr(iconv(emplois$CONTRAT_EMB, to='UTF-8', sub=''), 1, 2))
# And the sojourn time
emplois$time <- emplois$DUREE

# Now, creating the dataframe used for analysis
group <- cumsum(c(TRUE, diff(as.numeric(interaction(emplois$id, emplois$state)))!=0))
dataframe <- aggregate(time ~ group + id + state, data = cbind(emplois, group), sum)
dataframe <- result[order(dataframe$group), c("id", "state", "time")]


covariates <- data.frame(
  sex=factor(individus$Q1, levels=c(1,2), labels=c('H', 'F'))
)
covariates$etr <- factor(
  individus$FRAETR09, levels=c(1,2), labels=c('FRA', 'ETR')
  )
covariates$bp3 <- factor(
  individus$BP3, levels=c(1,2), labels=c('yes', 'no')
)
covariates$reg <- factor(individus$REGETAB)


source('src/tree_construction.R')
source('src/two_samples_test.R')

D <- 5
alg <- function(df1, df2){
  likelihood_ratio_test(df1, df2, D, law_sojourn='exponential')
}

tree <- build_tree(result, covariates, alg)
