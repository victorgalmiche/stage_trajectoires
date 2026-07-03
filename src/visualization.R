library(igraph)

df_to_traj <- function(dataframe){
  # To obtain the trajectories
  trajectories <- dataframe %>%
    group_by(id) %>%
    summarise(trajectoire = list(rep(state, times = time)), .groups = "drop") %>%
    tidyr::unnest_wider(trajectoire, names_sep = "_t")
  trajectories
}



plot_tree <- function(tree, show_pval = TRUE) {
  
  edges <- character(0)
  labels <- character(0)
  colors <- character(0)
  edge_labels <- character(0)
  compteur <- 0
  
  # Libellé du nœud (variable + n, sans le seuil qui passe maintenant sur les arêtes)
  node_label <- function(noeud) {
    if (noeud$type == "leaf") {
      return(sprintf("Feuille"))
    }
    split <- noeud$split
    label <- split$var
    if (show_pval && !is.null(split$pval)) {
      label <- paste0(label, sprintf("\n(p=%.3g)", split$pval))
    }
  }
  
  # Libellés gauche/droite pour les arêtes sortant d'un nœud de split
  branch_labels <- function(split) {
    if (split$type == "numeric") {
      list(
        left  = sprintf("< %.3g", split$threshold),
        right = sprintf("\u2265 %.3g", split$threshold)  # ≥
      )
    } else {
      list(
        left  = paste0("{", paste(split$left_levels, collapse = ", "), "}"),
        right = paste0("{", paste(split$right_levels, collapse = ", "), "}")
      )
    }
  }
  
  parcourir <- function(noeud, id_parent = NULL, sens = NULL) {
    compteur <<- compteur + 1
    id_actuel <- paste0("n", compteur)
    
    labels[id_actuel] <<- node_label(noeud)
    colors[id_actuel] <<- if (noeud$type == "leaf") "lightgreen" else "lightblue"
    
    if (!is.null(id_parent)) {
      edges <<- c(edges, id_parent, id_actuel)
      edge_labels <<- c(edge_labels, sens)
    }
    
    if (noeud$type == "node") {
      bl <- branch_labels(noeud$split)
      parcourir(noeud$left,  id_actuel, bl$left)
      parcourir(noeud$right, id_actuel, bl$right)
    }
    
    id_actuel
  }
  
  parcourir(tree)
  
  g <- make_graph(edges)
  V(g)$label <- labels[V(g)$name]
  V(g)$color <- colors[V(g)$name]
  E(g)$label <- edge_labels
  
  plot(g,
       layout = layout_as_tree(g, root = 1, mode = "out"),
       vertex.label = V(g)$label,
       vertex.shape = "rectangle",
       vertex.size = 70,
       vertex.size2 = 35,
       vertex.color = V(g)$color,
       vertex.label.cex = 0.7,
       vertex.label.family = "sans",
       edge.label = E(g)$label,
       edge.label.cex = 0.75,
       edge.label.color = "darkred",
       edge.arrow.size = 0.4)
}



# 
# # Visualize the trajectories
# seq <- seqdef(trajectories, 2:80)
# par(mfrow = c(2, 2))
# seqiplot(seq, with.legend=FALSE, border=NA)






