source('src/data_analysis/extract_data.R')

library(PCAmixdata)

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
              rename.level = TRUE, graph = FALSE)

