source('src/data/extract_data.R')
library(dplyr)
library(tidyr)

# Transforming empty spaces into NA
na_strings <- c("", " ", "NA", "N/A")
individus_clean <- individus %>%
  mutate(across(where(is.character), ~ ifelse(. %in% na_strings, NA, .)))

# Counting the proportion of NA for each variable
na_summary <- individus_clean %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_na") %>%
  arrange(desc(pct_na))

# Excluding some variables when too much NA
threshold_na <- 0
vars_to_exclude <- na_summary$variable[na_summary$pct_na > threshold_na]
cat(length(vars_to_exclude), "variables have at least one NA")


# Excluding the variables w/ NA and the weights and id
vars_to_exclude <- c(vars_to_exclude, "IDENT", "pondef")
vars_to_keep <- setdiff(names(individus_clean), vars_to_exclude)

individus_clean <- individus_clean[vars_to_keep]


# Extracting the numeric and categorical variables
detect_type <- function(x) {
  if (is.numeric(x)) {
    return("numeric")
  } else {
    return("categorical")
  }
}

var_types <- sapply(individus_clean, detect_type)
vars_cat <- names(var_types)[var_types == "categorical"]
vars_num <- names(var_types)[var_types == "numeric"]

cat("Numeric variables:", length(vars_num), "\n")
cat("Categorical variables:", length(vars_cat), "\n")

# Converting categorical cols in factor
individus_clean[vars_cat] <- lapply(individus_clean[vars_cat], as.factor)



# Correlation Function for mixt variables

# Cramer's V for categorical vs categorical
cramers_v <- function(x, y) {
  tbl <- table(x, y)
  if (any(dim(tbl) < 2)) return(NA)
  chi2 <- suppressWarnings(chisq.test(tbl, correct = FALSE)$statistic)
  n <- sum(tbl)
  k <- min(dim(tbl)) - 1
  if (k == 0 || n == 0) return(NA)
  sqrt(as.numeric(chi2) / (n * k))
}

# Correlation ratio for categorical vs numeric
correlation_ratio <- function(categories, values) {
  ok <- complete.cases(categories, values)
  categories <- droplevels(as.factor(categories[ok]))
  values <- values[ok]
  if (length(unique(categories)) < 2 || length(values) < 3) return(NA)
  ss_total <- sum((values - mean(values))^2)
  if (ss_total == 0) return(NA)
  moyennes <- tapply(values, categories, mean)
  effectifs <- tapply(values, categories, length)
  ss_between <- sum(effectifs * (moyennes - mean(values))^2)
  sqrt(ss_between / ss_total)
}


assoc <- function(x, y, type_x, type_y) {
  if (type_x == "numeric" && type_y == "numeric") {
    # Pearson correlation for numeric vs numeric - Abs to get between 0 and 1
    return(abs(suppressWarnings(cor(x, y, method = "pearson"))))
  } else if (type_x != "numeric" && type_y != "numeric") {
    return(cramers_v(x, y))
  } else if (type_x == "numeric") {
    return(correlation_ratio(y, x))
  } else {
    return(correlation_ratio(x, y))
  }
}


# Building the correlation matrix
n_vars <- length(vars_to_keep)
mat_corr <- matrix(NA, n_vars, n_vars, 
                   dimnames = list(vars_to_keep, vars_to_keep))

for (i in seq_len(n_vars)) {
  for (j in i:n_vars) {
    if (i == j) {
      mat_corr[i, j] <- 1
    } else {
      v <- assoc(individus_clean[[vars_to_keep[i]]], 
                 individus_clean[[vars_to_keep[j]]],
                 var_types[i], var_types[j])
      mat_corr[i, j] <- v
      mat_corr[j, i] <- v
    }
  }
  if (i %% 50 == 0) cat("Variable", i, "/", n_vars, "done\n")
}

# Hierarchical clustering with the correlation
dist_corr <- as.dist(1 - mat_corr)
hc <- hclust(dist_corr, method = "average")

heights <- sort(hc$height, decreasing = TRUE)
plot(heights[1:100], type = "b",
     xlab = "Fusion (decreasing order)", ylab = "Height", ylim = c(0,1),
     main = "Search of an elbow for optimal number of clusters")

# Cutting at a given threshold
threshold_correlation <- 0.5 
groups <- cutree(hc, h = 1 - threshold_correlation)
cat("Number of clusters:", length(unique(groups)), "\n")

# Arranging the groups in a dataframe
groups_df <- data.frame(variable = names(groups), group = groups) %>%
  arrange(group)

# Correlation plot of the remaining chosen covariates
chosen_vars <- c('PHD', 'CA13', 'CA20_17', 'CA22', 'CA24', 'ETR1', 'OS1', 
                 'OS3', 'PER1', 'Q1', 'Q16', 'Q2', 'Q31', 'Q53_13', 'Q9',
                 'SITMERE', 'SITPERE', 'ZUS')
corrplot(mat_corr[chosen_vars, chosen_vars],
         title = "Correlation plot of the chosen covariates", 
         tl.col = 'black', tl.cex=0.6, mar = c(0, 0, 2, 0))


individus_reduced <- individus_clean[chosen_vars]
