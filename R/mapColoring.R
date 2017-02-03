
## Make adjecency matrix from SpatialPolygons* object
getAM <- function(x) {
  
  requireNamespace("rgeos")
  
  nbmat <- rgeos::gIntersects(x, byid=T)
  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      if (nbmat[i,j]) {
        suppressWarnings(nbmat[i,j] <- ifelse(rgeos::gLength(rgeos::gIntersection(x[i,], x[j,]))>0, T, F))
      }
    }
  }
  diag(nbmat) <- FALSE
  return(nbmat)
} 

## Calculate relative luminance from sRGB matrix 
# See http://www.w3.org/TR/WCAG20/#relativeluminancedef
getRelLum <- function(color, maxColorValue=255) {
  color <- color/maxColorValue
  RGB <- apply(color, 2, function(x) ifelse(x<= 0.03928, x/12.92, ((x+0.055)/1.055) ^ 2.4))
  L <- 0.2126 * RGB[,1] + 0.7152 * RGB[,2] + 0.0722 * RGB[,3]
  return(L)
}

## Calculate color getContrast 
# See http://www.w3.org/TR/WCAG20/#getContrast-ratiodef
getContrast <- function(lum1, lum2) {
  ifelse(lum1 > lum2, (lum1 + 0.05) / (lum2 + 0.05), (lum2 + 0.05) / (lum1 + 0.05))
}

## Count occurrences of color in given nodes
getAmountColor <- function(node_indexes, color_number, coloring) {
  node_colors <- coloring[node_indexes]
  return(sum(node_colors==color_number))
}

## Get neighboring verteces
getNeighbors <- function(adj_mat, node_index) {
  nb <- which(adj_mat[node_index,])
  nb <- nb[!(nb==node_index)]
  return(nb)
}

## Greedy DSATUR graph coloring algorithm
# Reference: D.Brelaz (1979) - New Methods to color the vertices of a graph. Communications of the ACM: 22(4).
# Ported from Python implementation by Andrei Novikov (pyclustering@yandex.ru)
# Under GNU Public license
dsatur <- function(x, coloring=NULL) {
  
  if (is.null(coloring)) {  # Set up vertex coloring from scratch
    color_counter = 1
    adj_mat <- x
    diag(adj_mat) <- FALSE
    degrees = list()
    saturation_degrees = rep(0, nrow(adj_mat))
    
    coloring = rep(0, nrow(adj_mat))
    uncolored_vertices = 1:nrow(adj_mat)
    
    index_maximum_degree = 0
    maximum_degree = 0
    for (index_node in 1:nrow(adj_mat)) {
      # Fill degree of nodes in the input graph
      degrees[[length(degrees)+1]] <- c(sum(adj_mat[index_node,]), index_node)
      
      # And find node with maximal degree at the same time.
      if ((degrees[[index_node]])[1] > maximum_degree) {
        maximum_degree <- (degrees[[index_node]])[1]
        index_maximum_degree <- index_node
      }
    }
    
    # Update saturation
    neighbors = getNeighbors(adj_mat, index_maximum_degree)
    for (index_neighbor in neighbors){
      saturation_degrees[index_neighbor] <- saturation_degrees[index_neighbor] + 1
    }
    
    # Coloring the first node
    coloring[index_maximum_degree] = color_counter
    uncolored_vertices <- uncolored_vertices[-index_maximum_degree]
    
  } else {  # Set up vertex coloring given input coloring
    color_counter = max(coloring)
    adj_mat <- x
    diag(adj_mat) <- FALSE
    degrees = list()
    saturation_degrees = rep(0, nrow(adj_mat))
    
    uncolored_vertices = 1:nrow(adj_mat)
    uncolored_vertices <- uncolored_vertices[coloring==0]
    
    # Fill degree of nodes in the input graph and update saturation
    for (index_node in 1:nrow(adj_mat)) {
      # Set degree
      degrees[[length(degrees)+1]] <- c(sum(adj_mat[index_node,]), index_node)
      # Set saturation
      index_neighbors <- getNeighbors(adj_mat, index_node)
      index_saturation <- 0
      for (number_color in 1:color_counter) {
        if (getAmountColor(index_neighbors, number_color, coloring) > 0) {
          index_saturation <- index_saturation + 1
        }
      }
      saturation_degrees[index_node] <- index_saturation
    }
  }
  
  
  # Color the remaining verteces
  while(length(uncolored_vertices) > 0) {
    # Get maximum saturation degree
    maximum_satur_degree = -1
    for (index in uncolored_vertices) {
      if (saturation_degrees[index] > maximum_satur_degree) {
        maximum_satur_degree = saturation_degrees[index]
      }
    }
    
    # Get list of indexes with maximum saturation degree
    indexes_maximum_satur_degree <- c()
    for (index in uncolored_vertices) {
      if (saturation_degrees[index] == maximum_satur_degree) {
        indexes_maximum_satur_degree <- c(indexes_maximum_satur_degree, index)
      }
    }
    
    coloring_index = indexes_maximum_satur_degree[1]
    if (length(indexes_maximum_satur_degree) > 1) {  # There are more then one node with maximum saturation
      # Find node with maximum degree
      maximum_degree = -1
      for (index in indexes_maximum_satur_degree) {
        degree <- (degrees[[index]])[1]
        node_index <- (degrees[[index]])[2]
        if (degree > maximum_degree) {
          coloring_index = node_index
          maximum_degree = degree
        }
      }
    }
    
    # Coloring
    node_index_neighbors = getNeighbors(adj_mat, coloring_index)
    for (number_color in 1:(color_counter)) {
      if (getAmountColor(node_index_neighbors, number_color, coloring) == 0) {
        coloring[coloring_index] = number_color
        break;
      }
    }
    
    # If it has not been colored then
    if (coloring[coloring_index] == 0) {
      color_counter <- color_counter + 1  # Add new color
      coloring[coloring_index] = color_counter
    }
    
    # Remove node from uncolored set
    uncolored_vertices <- uncolored_vertices[!(uncolored_vertices==coloring_index)]
    
    # Update degree of saturation
    for (index_neighbor in node_index_neighbors) {
      subneighbors = getNeighbors(adj_mat, index_neighbor)
      if (getAmountColor(subneighbors, coloring[coloring_index], coloring) == 1) {
        saturation_degrees[index_neighbor] <- saturation_degrees[index_neighbor] + 1
      }
    }
  }
  # Return final coloring
  return(coloring)
}

## Get number of colors required to color map
getNColors <- function(x) {
  if (is(x, "matrix")) {
    adj_mat <- x
    diag(adj_mat) <- FALSE
  } else if (is(x, "SpatialPolygons")) {
    adj_mat <- getAM(x)
  } else {
    stop("x must be an adjacency matrix or a SpatialPolygons* object.")
  }
  nColors <- length(unique(dsatur(adj_mat)))
  return(nColors)
}

## Get coloring for adjancency matrix or SpatialPolygons* object
getColoring <- function(x) {
  if (is(x, "matrix")) {
    adj_mat <- x
    diag(adj_mat) <- FALSE
  } else if (is(x, "SpatialPolygons")) {
    adj_mat <- getAM(x)
  } else {
    stop("x must be an adjacency matrix or a SpatialPolygons* object.")
  }
  coloring <- dsatur(adj_mat)
  return(coloring)
}

## Cost function for optimal contrast fitting
oc <- function(cand.idx, coloring, adj_mat, cand.colors, beta) {
  
  # Get candidate colors as sRGB matrix
  this.colors <- cand.colors[cand.idx+1]
  cand.rgb <- t(col2rgb(this.colors))
  cand.rgb <- cand.rgb[coloring,]
  
  # Get contrast matrix
  cand.lum <- getRelLum(cand.rgb, 255)
  contrast.mat <- outer(cand.lum, cand.lum, getContrast) - 1
  
  # Calculate score value (average inverse contrast values)
  score <- (contrast.mat[adj_mat] + 1e-3)^(-beta)
  score <- mean(score)
  
  return(score)
}

## Get optimal contrast coloring for adjancency matrix or SpatialPolygons* object
getOptimalContrast <- function(x, colors, beta=1.5, seed=1, ...) {
  
  if (!requireNamespace("CEoptim", quietly = TRUE)) {
    stop("CEoptim needed for this function to work. Please install it.",
         call. = FALSE)
  }
  requireNamespace("CEoptim")
  
  set.seed(seed)
  
  if (is(x, "matrix")) {
    adj_mat <- x
    diag(adj_mat) <- FALSE
  } else if (is(x, "SpatialPolygons")) {
    adj_mat <- getAM(x)
  } else {
    stop("x must be an adjacency matrix or a SpatialPolygons* object.")
  }
  
  coloring <- dsatur(adj_mat)
  ncolors <- length(unique(coloring))
  M <- length(colors)
  if (M < ncolors) {
    stop(paste("Coloring this map requries at least", ncolors, "distinct colors."))
  }
  
  opt <- CEoptim::CEoptim(f=oc,
                          f.arg=list(coloring=coloring, adj_mat=adj_mat, cand.colors=colors, beta=beta),
                          maximize=FALSE,
                          discrete=list(categories=rep(M, ncolors)),
                          ...)
  opt.idx <- opt$optimizer$discrete+1
  opt.colors <- colors[opt.idx]
  opt.coloring <- opt.colors[coloring]
  return(opt.coloring)
}