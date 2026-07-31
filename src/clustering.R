# Identify the leaf of an individual
# obs is a row of covariates
get_leaf_path <- function(node, obs) {
  path <- ""
  while (node$type != 'leaf') {
    val <- obs[[node$split$var]]
    goes_left <- if (node$split$type == 'categorical') {
      val %in% node$split$left_levels
    } else {
      val < node$split$threshold
    }
    path <- paste0(path, if (goes_left) "L" else "R")
    node <- if (goes_left) node$left else node$right
  }
  path
}

get_forest_leaves <- function(forest, covariates) {
  n <- nrow(covariates)
  M <- length(forest)
  leaf_ids <- matrix(NA_character_, nrow=n, ncol=M)
  for (t in seq_len(M)){
    for (i in seq_len(n)){
      leaf_ids[i,t] <- get_leaf_path(forest[[t]], covariates[i, ])
    }
    cat('Tree ', t, 'finished')
  }
  leaf_ids
}

compute_consensus_matrix <- function(leaf_ids) {
  n <- nrow(leaf_ids)
  M <- ncol(leaf_ids)
  C <- matrix(0, n, n)

  for (t in seq_len(M)) {
    lp <- leaf_ids[, t]
    same <- outer(lp, lp, "==")
    C <- C + same
  }
  C
}


C <- compute_consensus_matrix(leaf_ids)
C <- C/500
D <- as.dist(1 - C)
hc <- hclust(D, method = "average")
groups <- cutree(hc, k=10)

