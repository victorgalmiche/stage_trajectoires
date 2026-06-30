source('src/data_analysis/extract_data.R')


# Creating a new column w/ PHD
mapping <- rep(1, 18)
mapping[2:5] <- 2
mapping[6:12] <- 3
mapping[13:18] <- 4
individus$PHD_NEW <- mapping[as.integer(substr(individus$PHD, 1, 2))]


cols_quali <- c('Q1', 'Q2', 'Q16', 'Q31', 'OS1', 'OS3_1', 'OS3_2', 'OS3_3',
                'ETR1', 'PER1', 'SITPERE', 'SITMERE', 'CA13', 'CA22')

cols_quanti <- c('PHD_NEW', 'AGE10')


library(PCAmixdata)
individus_clean <- individus[, cols_quanti]
individus_clean[cols_quali] <- lapply(individus[cols_quali], as.factor)


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

