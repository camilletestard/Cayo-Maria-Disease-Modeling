
# AdjacencyFromEdgelist.R ####

AdjacencyFromEdgelist <- function(elist, check.full=TRUE) {
  # Guardians
  stopifnot(is(elist, "data.frame"),
            3 == ncol(elist),
            is.numeric(elist[,3]))
  
  if(check.full) elist <- EdgelistFill(elist)  # I assume that this is sorted by first col, then second, 
  # so the third col has consecutive rows of the adjacency matrix
  
  nodelist <- sort(unique(elist[,1]))
  n <- length(nodelist)
  weights <- elist[,3]
  
  adjacency <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n) {
    adjacency[i,] <- weights[((i-1)*n+1):(i*n)]
  }
  out <- list(adjacency=adjacency, nodelist=nodelist)
  return(out)
}

EdgelistFill <- function(elist, 
                         fillBlanksWith=0, 
                         nodelist) {
  # Guardians
  if( (is(elist, "matrix") || is(elist, "data.frame")) && 2 == ncol(elist) ) {
    # elist has no weights; perhaps it was the result of get.edgelist(g)
    elist <- cbind(elist, 1)  # put 1s in for each edge
  } else if( (is(elist, "matrix") || is(elist, "data.frame")) && 3 == ncol(elist) ) {
    # elist has a third column, presumably with weights
    stopifnot( is.numeric(elist[,3]) )
  } else {
    stop("elist must have two or three columns")
  }
  
  if( missing(nodelist) ) {
    nodelist <- sort(union(unique(elist[,1]), unique(elist[,2])))
  } else {
    stopifnot(all(elist[,1] %in% nodelist),
              all(elist[,2] %in% nodelist))
    nodelist <- sort(nodelist)
  }
  
  out <- expand.grid(nodelist, nodelist)
  out <- data.frame(out[,2], out[,1], fillBlanksWith)  # ensures sorted by first column, then second column
  if( is.null(names(elist)) ) {
    names(out) <- c("fromnode", "tonode", "weight")
  } else {
    names(out) <- names(elist)
  }
  
  for(i in 1:nrow(elist)) {
    if( is.factor(elist[,1]) ) {
      match1 <- levels(elist[i,1])[as.numeric(elist[i,1])]  # "unfactor" in case factor levels don't agree
    } else {
      match1 <- elist[i,1]
    }
    if( is.factor(elist[,2]) ) {
      match2 <- levels(elist[i,2])[as.numeric(elist[i,2])]  # "unfactor" in case factor levels don't agree
    } else {
      match2 <- elist[i,2]
    }
    out[out[,1]==match1 & out[,2]==match2, 3] <- elist[i,3]
  }
  return( out )
}