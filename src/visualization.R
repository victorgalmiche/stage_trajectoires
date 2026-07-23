library(igraph)

plot_tree <- function(tree, show_pval = TRUE, hgap = 1, vgap = 2,
                      vertex_size = 100, vertex_size2 = 150) {
  
  edges <- character(0)
  labels <- character(0)
  colors <- character(0)
  edge_labels <- character(0)
  compteur <- 0
  depths  <- list()   # profondeur de chaque nœud
  
  node_label <- function(noeud) {
    if (noeud$type == "leaf") {
      return("Leaf")
    }
    split <- noeud$split
    label <- split$var
    if (show_pval && !is.null(split$pval)) {
      label <- paste0(label, sprintf("\n(p=%.3g)", split$pval))
    }
    label
  }
  
  branch_labels <- function(split) {
    if (split$type == "numeric") {
      list(
        left  = sprintf("< %.3g", split$threshold),
        right = sprintf("\u2265 %.3g", split$threshold)
      )
    } else {
      list(
        left  = paste0("{", paste(split$left_levels, collapse = ", "), "}"),
        right = paste0("{", paste(split$right_levels, collapse = ", "), "}")
      )
    }
  }
  
  # 1er parcours : construit les nœuds/arêtes + calcule la profondeur
  parcourir <- function(noeud, id_parent = NULL, sens = NULL, profondeur = 0) {
    compteur <<- compteur + 1
    id_actuel <- paste0("n", compteur)
    
    labels[id_actuel] <<- node_label(noeud)
    colors[id_actuel] <<- if (noeud$type == "leaf") "lightgreen" else "lightblue"
    depths[[id_actuel]] <<- profondeur
    
    if (!is.null(id_parent)) {
      edges <<- c(edges, id_parent, id_actuel)
      edge_labels <<- c(edge_labels, sens)
    }
    
    if (noeud$type == "node") {
      bl <- branch_labels(noeud$split)
      parcourir(noeud$left,  id_actuel, bl$left,  profondeur + 1)
      parcourir(noeud$right, id_actuel, bl$right, profondeur + 1)
    }
    
    id_actuel
  }
  
  parcourir(tree)
  
  g <- make_graph(edges)
  V(g)$label <- labels[V(g)$name]
  V(g)$color <- colors[V(g)$name]
  E(g)$label <- edge_labels
  
  # 2. Calcul manuel des positions x (feuilles espacées régulièrement,
  #    parents centrés sur leurs enfants) et y (= -profondeur)
  xpos <- setNames(rep(NA_real_, vcount(g)), V(g)$name)
  leaf_counter <- 0
  
  assign_x <- function(id) {
    enfants <- neighbors(g, id, mode = "out")
    if (length(enfants) == 0) {
      leaf_counter <<- leaf_counter + 1
      xpos[id] <<- leaf_counter
    } else {
      for (e in enfants) assign_x(V(g)$name[e])
      enfants_x <- xpos[V(g)$name[enfants]]
      xpos[id] <<- mean(enfants_x)
    }
  }
  racine <- V(g)$name[1]
  assign_x(racine)
  
  x <- xpos[V(g)$name] * hgap
  y <- -sapply(V(g)$name, function(n) depths[[n]]) * vgap
  lay <- cbind(x, y)
  
  plot(g,
       layout = lay,
       rescale = FALSE,
       xlim = range(x) + c(-hgap, hgap) * 0.6,
       ylim = range(y) + c(-vgap, vgap) * 0.6,
       asp = 0,
       vertex.label = V(g)$label,
       vertex.shape = "rectangle",
       vertex.size = vertex_size,
       vertex.size2 = vertex_size2,
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






