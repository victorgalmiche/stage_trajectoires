source('src/data/extract_data.R')
library(PCAmixdata)

# Histogramme du nombre de modalités des colonnes quali
hist(sapply(covariates[cols_quali], nlevels))

# Blocs finaux 
X.quanti <- covariates[, cols_quanti, drop=FALSE] 
X.quali <- covariates[, cols_quali, drop=FALSE] 
cat("Prêt pour PCAmix :\n") 
cat(" -", ncol(X.quanti), "variables quantitatives\n") 
cat(" -", ncol(X.quali), "variables qualitatives\n") 
cat(" -", nrow(X.quanti), "individus\n") 

res.pcamix <- PCAmix(X.quanti = X.quanti, X.quali = X.quali, 
                     rename.level = TRUE, graph = TRUE)

# Résumé et variance expliquée
res.pcamix$eig

# Coordonnées des variables (quanti + quali) sur les axes
res.pcamix$quanti$coord
res.pcamix$quali$coord

# Graphique combiné
plot(res.pcamix, choice = "cor")   # variables quantitatives (cercle des corrélations)
plot(res.pcamix, choice = "levels") # modalités des variables qualitatives
plot(res.pcamix, choice = "sqload") # carrés des loadings (quanti + quali ensemble)


library(rcompanion)  # cramerV()
library(corrplot)

# Toutes les variables concernées
all_vars <- c(cols_quali, cols_quanti)
n <- length(all_vars)

mat_assoc <- matrix(NA, n, n, dimnames = list(all_vars, all_vars))

for (i in 1:n) {
  for (j in 1:n) {
    var_i <- all_vars[i]
    var_j <- all_vars[j]
    
    if (i == j) {
      mat_assoc[i, j] <- 1
      
    } else if (var_i %in% cols_quali & var_j %in% cols_quali) {
      # quali - quali : V de Cramér
      mat_assoc[i, j] <- cramerV(covariates[[var_i]], covariates[[var_j]])
      
    } else if (var_i %in% cols_quanti & var_j %in% cols_quali) {
      # quanti - quali : eta² (ANOVA)
      formule <- as.formula(paste(var_i, "~", var_j))
      modele <- aov(formule, data = covariates)
      ss <- summary(modele)[[1]]$`Sum Sq`
      mat_assoc[i, j] <- ss[1] / sum(ss)
      
    } else if (var_i %in% cols_quali & var_j %in% cols_quanti) {
      # symétrique du cas précédent
      formule <- as.formula(paste(var_j, "~", var_i))
      modele <- aov(formule, data = covariates)
      ss <- summary(modele)[[1]]$`Sum Sq`
      mat_assoc[i, j] <- ss[1] / sum(ss)
    } else if (var_i %in% cols_quanti & var_j %in% cols_quanti) {
      # quanti - quanti : Pearson (valeur absolue pour rester cohérent 0-1)
      mat_assoc[i, j] <- abs(cor(covariates[[var_i]], covariates[[var_j]], 
                                 use = "pairwise.complete.obs", method = "pearson"))
    }
  }
}

round(mat_assoc, 2)

corrplot(mat_assoc, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         col = COL1("YlOrRd"), number.cex = 0.7,
         title = "Matrice d'association (Cramér V / eta²)", 
         mar = c(0,0,2,0))
