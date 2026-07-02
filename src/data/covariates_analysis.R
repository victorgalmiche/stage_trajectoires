source('src/data/extract_data.R')
library(PCAmixdata)

# Histogramme du nombre de modalités des colonnes quali
hist(sapply(covariates[cols_quali], nlevels))

# Blocs finaux 
X.quanti <- covariates[cols_quanti] 
X.quali <- covariates[cols_quali] 
cat("Prêt pour PCAmix :\n") 
cat(" -", ncol(X.quanti), "variables quantitatives\n") 
cat(" -", ncol(X.quali), "variables qualitatives\n") 
cat(" -", nrow(X.quanti), "individus\n") 

# Lancer PCAmix 
res <- PCAmix(X.quanti = as.data.frame(X.quanti), 
              X.quali = as.data.frame(X.quali), 
              rename.level = TRUE, graph = FALSE)
