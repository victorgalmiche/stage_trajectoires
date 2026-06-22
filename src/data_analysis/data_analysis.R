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

