# **************************************************************************
#
# Functions beyond the cvn package needed for the analysis
#
# **************************************************************************

# Edges and graphs ------------------------------------------------------------------
## convert cvn to igraph
cvn2igraph <- function(cvn){
  
  cvn.igraph <- lapply(cvn$adj_matrices, function(x){
    lapply(x, function(z){
      igraph::graph_from_adjacency_matrix(unlist(z), mode = 'undirected') 
    })
  })
  cvn.igraph
}  

## make list with edges based on an adjacency matrix  
make_edge_list <- function(adj_mat){
  el <- igraph::as_edgelist(
    igraph::graph_from_adjacency_matrix(adj_mat, mode = 'undirected'))
  el <- paste(el[,1], el[,2], sep = "-")
  el
}


# outcome is a list of edges that are shared within and between groups
# just for a single cvn ()
reduce_my_adj <- function(cvn){
  
  if(length(cvn$adj_matrices) > 1) warning("Function is just applied on cvn$adj_matrices[[1]]. Use reduce_cvn()")
  
  combi <- gtools::combinations(cvn$m, 2)
  out <- tmp <- vector(mode = 'list', length = nrow(combi))
  
  for(k in 1:nrow(combi)){
    reduced <- Reduce('+', cvn$adj_matrices[[1]][c(combi[k,1],combi[k,2])])
    reduced[reduced == 2] <- 0
    tmp[[k]] <- reduced
  }
  
  out <- lapply(tmp, make_edge_list)
  namen <- gtools::combinations(cvn$m, 2, c('N00', 'N01', 'N02', 'N10', 'N11', 'N12', 'N20', 'N21', 'N22'))
  namen <- paste(namen[,1], namen[,2], sep = '-')
  names(out) <- namen
  return(out)
}


# adds network descriptives to igraph object
make_my_graph <- function(g, labels = NULL){

  if(!is.null(labels)) V(g)$name = labels

  deg <- igraph::degree(g, v = V(g), mode = "all")
  bet <- igraph::betweenness(g, v = V(g), directed = FALSE, normalized = TRUE)
  clo <- igraph::closeness(g, mode = "total", normalized = TRUE)
  betedge <- igraph::edge_betweenness(g, directed = FALSE)
  
  ### Set graph attributes
  g <- set_vertex_attr(g, "degree", index = V(g), value = deg)
  g <- set_vertex_attr(g, "betweenness", index = V(g), value = bet)
  g <- set_vertex_attr(g, "closeness", index = V(g), value = clo)
  g <- set_edge_attr(g, "betweenness", index = E(g), value = betedge)
  
  ## components
  comps <- decompose.graph(g)
  comps.size <- sapply(comps, gorder)
  gc <- comps[which(comps.size > 1)]
  
  if(length(gc) == 1){
    gc <- gc[[1]]
  } else if(length(gc) > 1){

    gc <- igraph::union(gc[[1]], gc[[2]])
    
    deg <- igraph::degree(gc, v = V(gc), mode = "all")
    bet <- igraph::betweenness(gc, v = V(gc), directed = FALSE, normalized = TRUE)
    clo <- igraph::closeness(gc, mode = "total", normalized = TRUE)
    betedge <- igraph::edge_betweenness(gc, directed = FALSE)
    
    ### Set graph attributes
    gc <- set_vertex_attr(gc, "degree", index = V(gc), value = deg)
    gc <- set_vertex_attr(gc, "betweenness", index = V(gc), value = bet)
    gc <- set_vertex_attr(gc, "closeness", index = V(gc), value = clo)
    gc <- set_edge_attr(gc, "betweenness", index = E(gc), value = betedge)
  }
  return(gc)
}



graph_descriptives <- function(g){
  
  # ---   output matrix   ---
  out <- matrix(NA, nrow = 16, ncol = 2,
                dimnames = list(c("Graph density",
                                  "Graph transitivity",
                                  "Graph reciprocity",
                                  "Longest shortest path (diameter)",
                                  "Avg. shortest path length",
                                  "Articulation points",
                                  "Number of nodes",
                                  "Number of edges",
                                  "Max. Degree and nodes",
                                  "Avg. Degree",
                                  "Max. Betweenness",
                                  "Avg. Betweenness centrality [0-100]",
                                  "Max. Closeness",
                                  "Avg. Closeness centrality  [0-100]",
                                  "Max. Jaccard",
                                  "Avg. Jaccard similarity  [0-100]"), c("value", "node")))
  # graph base

  out[1,1] <- edge_density(g)
  out[2,1] <- transitivity(g)       #~20% connected triples close to form triangles.
  out[3,1] <- reciprocity(g)
  out[4,1] <- length(get_diameter(g, weights = NA))
  out[5,1] <- mean_distance(g, weights = NA)
  out[6,1] <- length(articulation.points(g)) # A single vertex that disconnects the graph
  out[7,1] <- length(V(g))
  out[8,1] <- length(E(g))
  out[9,]  <- c(max(V(g)$degree), which(V(g)$degree == max(V(g)$degree)))
  out[10,1] <- mean(V(g)$degree)
  out[11,] <- c(max(V(g)$betweenness * 100), which(V(g)$betweenness == max(V(g)$betweenness)))
  out[12,1] <- mean(V(g)$betweenness * 100)
  if(length(which(V(g)$closeness == max(V(g)$closeness))) == 1){
    out[13,] <- c(max(V(g)$closeness * 100), which(V(g)$closeness == max(V(g)$closeness))) 
    } else {
    cat('Closeness nodes: ', which(V(g)$closeness == max(V(g)$closeness)))    
  }
  out[14,1] <- mean(V(g)$closeness * 100)
  out[15,] <- c(max(E(g)$betweenness), which(E(g)$betweenness == max(E(g)$betweenness)))
  out[16,1] <- mean(E(g)$betweenness)

  
  # return
  out <- as.data.frame(out)
  out[,1] <- sapply(out[,1], round_custom)
  return(out)
  
}



# Rounds a number to a maximum of 2 decimal places, but omits trailing zeros ------
round_custom <- function(x, digits = 2) {
  # Round the number to the specified number of digits
  rounded <- round(x, digits)
  
  # Format the number to remove trailing zeros
  formatted <- format(rounded, nsmall = 0, scientific = FALSE, trim = TRUE)
  
  return(formatted)
}

