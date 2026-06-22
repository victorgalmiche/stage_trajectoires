library(haven)
library(aws.s3)

BUCKET <- 'victorgalmiche'
FOLDER <- 'stage-trajectoires/lil-1439/lil-1439.dta/Stata'
FILE_INDIV <- paste(FOLDER, 'g107individusvf.dta', sep='/')

individus <- aws.s3::s3read_using(
  FUN = haven::read_dta,
  object = FILE_INDIV,
  bucket = BUCKET,
  opts = list("region"="")
)

colSums(is.na(individus))
mean(is.na(individus)) 



library(PCAmixdata)

individus_clean <- individus[, colSums(is.na(individus)) == 0]
cols_quanti <- names(individus_clean)[sapply(individus_clean, is.numeric)]
cols_quali  <- names(individus_clean)[sapply(individus_clean, function(x)
is.character(x) || is.factor(x))]

cat("Quantitatives :", cols_quanti, "\n")
cat("Qualitatives  :", cols_quali,  "\n")

individus_clean[cols_quali] <- lapply(individus_clean[cols_quali], as.factor)


# Colonnes avec variance nulle (inutiles pour PCA) 
zero_var <- cols_quanti[sapply(individus_clean[cols_quanti], 
                               function(x) var(x, na.rm = TRUE) == 0)] 
cat("Variance nulle :", zero_var, "\n") 

# Variables quali avec 1 seule modalité (inutiles) 
unique_modal <- cols_quali[sapply(df_clean[cols_quali], 
                                  function(x) nlevels(x) < 2)] 
cat("Modalité unique :", unique_modal, "\n") 

# Supprimer ces colonnes 
cols_quanti <- setdiff(cols_quanti, zero_var) 
cols_quali <- setdiff(cols_quali, unique_modal)


# Blocs finaux 
X.quanti <- individus_clean[cols_quanti] 
X.quali <- individus_clean[cols_quali] 
cat("Prêt pour PCAmix :\n") 
cat(" -", ncol(X.quanti), "variables quantitatives\n") 
cat(" -", ncol(X.quali), "variables qualitatives\n") 
cat(" -", nrow(X.quanti), "individus\n") 

# Lancer PCAmix 
res <- PCAmix(X.quanti = as.data.frame(X.quanti), 
              X.quali = as.data.frame(X.quali), 
              rename.level = TRUE, graph = TRUE)
