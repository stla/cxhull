#' Convex hull
#' @description Computes the convex hull of a set of points.
#' @param points numeric matrix, one point per row
#' @param triangulate logical, whether to triangulate the convex hull
#' @return A list providing a lot of information about the convex hull. See 
#'   the \strong{README} file for details.
#' @export
#' @useDynLib cxhull, .registration = TRUE
#' @examples
#' library(cxhull)
#' points <- rbind(
#'  c(0.5,0.5,0.5),
#'  c(0,0,0),
#'  c(0,0,1),
#'  c(0,1,0),
#'  c(0,1,1),
#'  c(1,0,0),
#'  c(1,0,1),
#'  c(1,1,0),
#'  c(1,1,1)
#' )
#' cxhull(points)
cxhull <- function(points, triangulate = FALSE){
  stopifnot(isBoolean(triangulate))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.")
  }
  dimension <- ncol(points)
  if(dimension < 2L){
    stop("The dimension must be at least 2.")
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.")
  }
  if(any(is.na(points))){
    stop("Missing values are not allowed.")
  }
  if(anyDuplicated(points)){
    stop("There are some duplicated points.")
  }
  storage.mode(points) <- "double"
  errfile <- tempfile(fileext = ".txt")
  hull <- tryCatch({
    .Call("cxhull_", points, as.integer(triangulate), errfile)
  }, error = function(e){
    cat(readLines(errfile), sep = "\n")
    stop(e)
  })
  hull[["volume"]] <- 1/dimension * 
    sum(sapply(hull[["facets"]],
               function(f) crossprod(f[["center"]], f[["normal"]])) *
          sapply(hull[["facets"]], `[[`, "volume"))
  if(dimension == 3L){
    attr(hull, "3d") <- TRUE
  }
  if(triangulate){
    attr(hull, "triangulate") <- TRUE
  }
  hull
}

#' @title Vertices and edges of convex hull
#' @description Computes the vertices and the edges of the convex hull of a set 
#'   of points.
#' @param points numeric matrix, one point per row; it must contain at least 
#'   three columns (the two-dimensional case is not implemented yet)
#' @param adjacencies Boolean, whether to return the vertex adjacencies
#' @param orderEdges Boolean, whether to order the edges in the output
#' @return A list with two fields: \code{vertices} and \code{edges}. The 
#'   \code{vertices} field is a list which provides an id for each vertex and 
#'   its coordinates. If \code{adjacencies=TRUE}, it provides in addition the 
#'   ids of the adjacent vertices for each vertex. The \code{edges} fields is 
#'   an integer matrix with two columns. Each row provides the two ids of the 
#'   vertices of the corresponding edge.
#' @export
#' @useDynLib cxhull, .registration = TRUE
#' @examples library(cxhull)
#' # let's try with the hexacosichoron (see `?hexacosichoron`)
#' #   it is convex so its convex hull is itself
#' VE <- cxhullEdges(hexacosichoron)
#' edges <- VE[["edges"]]
#' random_edge <- edges[sample.int(720L, 1L), ]
#' A <- hexacosichoron[random_edge[1L], ]
#' B <- hexacosichoron[random_edge[2L], ]
#' sqrt(c(crossprod(A - B))) # this is 2/phi
#' # Now let's project the polytope to the H4 Coxeter plane 
#' phi <- (1 + sqrt(5)) / 2
#' u1 <- c(
#'   0, 
#'   2*phi*sin(pi/30), 
#'   0, 
#'   1
#' )
#' u2 <- c(
#'   2*phi*sin(pi/15), 
#'   0, 
#'   2*sin(2*pi/15), 
#'   0
#' )
#' u1 <- u1 / sqrt(c(crossprod(u1)))
#' u2 <- u2 / sqrt(c(crossprod(u2)))
#' # projections to the Coxeter plane
#' proj <- function(v){
#'   c(c(crossprod(v, u1)), c(crossprod(v, u2)))
#' }
#' points <- t(apply(hexacosichoron, 1L, proj))
#' # we will assign a color to each edge  
#' #   according to the norms of its two vertices
#' norms2 <- round(apply(points, 1L, crossprod), 1L)
#' ( tbl <- table(norms2) )
#' #> 0.4 1.6 2.4 3.6 
#' #>  30  30  30  30 
#' values <- as.numeric(names(tbl))
#' grd <- as.matrix(expand.grid(values, values)) 
#' grd <- grd[grd[, 1L] <= grd[, 2L], ]
#' pairs <- apply(grd, 1L, paste0, collapse = "-")
#' colors <- hcl.colors(nrow(grd), palette = "Hawaii", rev = TRUE)
#' if(require("colorspace")) {
#'   colors <- colorspace::darken(colors, amount = 0.3)
#' }
#' names(colors) <- pairs
#' # plot ####
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(
#'   points[!duplicated(points), ], pch = 19, cex = 0.3, asp = 1, 
#'   axes = FALSE, xlab = NA, ylab = NA
#' )
#' for(i in 1L:nrow(edges)){
#'   twopoints <- points[edges[i, ], ]
#'   nrms2 <- round(sort(apply(twopoints, 1L, crossprod)), 1L)
#'   pair <- paste0(nrms2, collapse = "-")
#'   lines(twopoints, lwd = 0.5, col = colors[pair])
#' }
#' par(opar)
cxhullEdges <- function(points, adjacencies = FALSE, orderEdges = FALSE){
  stopifnot(isBoolean(adjacencies), isBoolean(orderEdges))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.")
  }
  dimension <- ncol(points)
  if(dimension < 2L){
    stop("The dimension must be at least 2.")
  }
  if(dimension == 2L){
    stop("This function is not implemented yet for the two-dimensional case.")
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.")
  }
  if(any(is.na(points))){
    stop("Missing values are not allowed.")
  }
  if(anyDuplicated(points)){
    stop("There are some duplicated points.")
  }
  storage.mode(points) <- "double"
  errfile <- tempfile(fileext = ".txt")
  hullEdges <- tryCatch({
    .Call(
      "cxhullEdges_", points, 
      as.integer(orderEdges), 
      errfile, PACKAGE = "cxhull"
    )
  }, error = function(e){
    cat(readLines(errfile), sep = "\n")
    stop(e)
  })
  if(adjacencies){
    edges <- hullEdges[["edges"]]
    vertices <- hullEdges[["vertices"]]
    for(i in 1L:length(vertices)){
      id <- vertices[[i]][["id"]]
      vertices[[i]][["neighbors"]] <- 
        sort(c(edges[edges[, 1L] == id, 2L], edges[edges[, 2L] == id, 1L]))
    }
    hullEdges[["vertices"]] <- vertices
  }
  if(dimension == 3L){
    attr(hullEdges, "3d") <- TRUE
  }
  hullEdges
}



#' @title Convex hull vertices
#' @description The coordinates of the vertices of a 3D convex hull.
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3D points
#'
#' @return A matrix with three columns. Each row represents the coordinates of 
#'   a vertex and the row names are the ids of the vertices.
#' @export
VerticesXYZ <- function(hull){
  if(!isTRUE(attr(hull, "3d"))){
    stop("This function is restriced to the 3D case.")
  }
  t(vapply(hull[["vertices"]], `[[`, numeric(3L), "point"))
}

#' @title Triangles of a triangulated 3D convex hull
#' @description Coordinates of the vertices of the triangles of a 
#'  triangulated 3D convex hull.
#'  
#' @param hull an output of \code{\link{cxhull}} applied to 3D points and 
#'   with the option \code{triangulate=TRUE}
#'
#' @return A matrix with three columns. Each row represents the coordinates of 
#'   a vertex of a triangle.
#' @export
#'
#' @examples library(cxhull)
#' library(rgl)
#' dodecahedron <- t(dodecahedron3d()$vb[-4L, ])
#' hull <- cxhull(dodecahedron, triangulate = TRUE)
#' triangles <- TrianglesXYZ(hull)
#' triangles3d(triangles, color = "firebrick")
TrianglesXYZ <- function(hull){
  Vertices <- VerticesXYZ(hull)
  if(!isTRUE(attr(hull, "triangulate"))){
    stop(
      "You didn't compute the convex hull with the option `triangulate=TRUE`."
    )
  }
  Facets <- hull[["facets"]]
  ntriangles <- length(Facets)
  triangles <- matrix(NA_real_, nrow = 3L*ntriangles, ncol = 3L)
  rownames(triangles) <- paste0(
    rep(c("A", "B", "C"), times = ntriangles),
    rep(1L:ntriangles, each = 3L)
  )
  colnames(triangles) <- c("x", "y", "z")
  #triangles <- cbind(as.data.frame(triangles), family = NA_integer_)
  families <- rep(NA_integer_, ntriangles)
  for(i in 1L:ntriangles){
    facet <- Facets[[i]]
    vertices <- as.character(facet[["vertices"]])
    if(facet[["orientation"]] == -1L){
      vertices <- vertices[c(1L, 3L, 2L)]
    }
    triangles[(3L*i-2L):(3L*i), ] <- Vertices[vertices, ]
    family <- facet[["family"]]
    if(!is.na(family)){
      families[i] <- family
    }
  }
  attr(triangles, "families") <- families
  triangles
}

#' @title Edges of a triangulated 3D convex hull
#' @description Edges of a triangulated 3D convex hull given by the ids of 
#'   the vertices in a matrix, plus a column indicating the border edges.
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3D points and 
#'   with the option \code{triangulate=TRUE}
#'
#' @return A character matrix with three columns. Each row provides the ids of 
#'   the two vertices of an edge, and a yes/no indicator of whether the edge 
#'   is a border edge.
#' @export
#'
#' @examples library(cxhull)
#' library(rgl)
#' dodecahedron <- t(dodecahedron3d()$vb[-4L, ])
#' hull <- cxhull(dodecahedron, triangulate = TRUE)
#' triangles <- TrianglesXYZ(hull)
#' triangles3d(triangles, color = "yellow")
#' edges <- EdgesAB(hull)
#' trueEdges <- edges[edges[, 3L] == "yes", c(1L, 2L)]
#' otherEdges <- edges[edges[, 3L] == "no", c(1L, 2L)]
#' vertices <- VerticesXYZ(hull)
#' for(i in 1:nrow(trueEdges)){
#'   lines3d(vertices[trueEdges[i, ], ], color = "blue", lwd = 3)
#' }
#' for(i in 1:nrow(otherEdges)){
#'   lines3d(vertices[otherEdges[i, ], ], color = "red", lwd = 3)
#' }
EdgesAB <- function(hull){
  if(!isTRUE(attr(hull, "3d"))){
    stop("This function is restriced to the 3D case.")
  }
  if(!isTRUE(attr(hull, "triangulate"))){
    stop(
      "You didn't compute the convex hull with the option `triangulate=TRUE`."
    )
  }
  ridges <- hull[["ridges"]]
  nedges <- length(ridges)
  edges <- matrix("yes", nrow = nedges, ncol = 3L)
  colnames(edges) <- c("A", "B", "border")
  Facets <- hull[["facets"]]
  for(i in 1L:nedges){
    ridge <- ridges[[i]]
    edges[i, c(1L, 2L)] <- as.character(ridge[["vertices"]])
    facets <- ridge[["ridgeOf"]]
    if(!is.na(f1 <- Facets[[facets[1L]]][["family"]])){
      f2 <- Facets[[facets[2L]]][["family"]]
      if(!is.na(f2) && f1 == f2){
        edges[i, 3L] <- "no"
      }
    }
  }
  edges
}

#' @title Edges coordinates
#' @description The coordinates of the extremities of the edges in a matrix, 
#'   plus a column indicating which edges are border edges.
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3D points and 
#'   with the option \code{triangulate=TRUE}
#'
#' @return A numeric matrix with four columns. The first three values of a row 
#'   are the coordinates of a vertex at the extremity of an edge, and the 
#'   fourth column indicates whether the edge is a border edge.
#' @export
EdgesXYZ <- function(hull){
  edges <- EdgesAB(hull)
  nedges <- nrow(edges)
  Vertices <- t(vapply(hull[["vertices"]], `[[`, numeric(3L), "point"))
  Edges <- matrix(NA_real_, nrow = 2L*nedges, ncol = 4L)
  rownames(Edges) <- paste0(
    rep(c("A", "B"), times = nedges),
    rep(1L:nedges, each = 2L)
  )
  colnames(Edges) <- c("x", "y", "z", "border")
  for(i in 1L:nedges){
    edge <- edges[i, ]
    border <- ifelse(edge[3L] == "yes", 1, 0)
    edge <- edge[c(1L, 2L)]
    Edges[(2L*i-1L):(2L*i), ] <- cbind(Vertices[edge, ], border)
  }  
  Edges
}


singleRows <- function(M){
  rownames(M) <- paste0(M[, 1L], "-", M[, 2L])
  tab12 <- table(rownames(M))
  tab12 <- tab12[-which(tab12 == 2L)]
  M[names(tab12), ]
}

polygonize <- function(edges){
  nedges <- nrow(edges)
  vs <- edges[1L, ]
  v <- vs[2L]
  edges <- edges[-1L, ]
  for(. in 1L:(nedges-2L)){
    j <- which(apply(edges, 1L, function(e) v %in% e))
    v <- edges[j, ][which(edges[j, ] != v)]
    vs <- c(vs, v)
    edges <- edges[-j, ]
  }
  cbind(as.character(c(vs)), as.character(c(vs[-1L], vs[1L])))
}


#' @title Summary of 3D convex hull
#' @description Summary of a triangulated 3D convex hull
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3D points and 
#'   with the option \code{triangulate=TRUE}
#'
#' @return A list with the vertices and the facets.
#' @export
#'
#' @examples library(cxhull)
#' # pyramid
#' pts <- rbind(
#'   c(0, 0, 0), 
#'   c(1, 0, 0), 
#'   c(1, 1, 0), 
#'   c(0.5, 0.5, 1), 
#'   c(0.5, 0.5, 0.9),
#'   c(0, 1, 0)
#' )
#' hull <- cxhull(pts, triangulate = TRUE)
#' hullSummary(hull)
hullSummary <- function(hull){
  if(!isTRUE(attr(hull, "3d"))){
    stop("This function is restriced to the 3D case.")
  }
  if(!isTRUE(attr(hull, "triangulate"))){
    stop(
      "You didn't compute the convex hull with the option `triangulate=TRUE`."
    )
  }
  hullFacets <- hull[["facets"]]
  families <- vapply(hullFacets, `[[`, integer(1L), "family")
  for(i in seq_along(hullFacets)){
    attr(hullFacets[[i]][["edges"]], "id") <- i
  }
  otherFacets <- Filter(function(x) !is.na(x[["family"]]), hullFacets)
  Other <- tapply(lapply(lapply(
    otherFacets, `[[`, "edges"
  ),
  identity
  ), families[!is.na(families)], function(x){
    ids <- vapply(x, attr, integer(1L), "id")
    list(
      family   = families[ids[1L]],
      facetids = ids,
      edges    = polygonize(singleRows(do.call(rbind, x)))
    )
  }, simplify = FALSE)
  nOtherFacets <- length(Other)
  
  triFacets <- Filter(function(x) is.na(x[["family"]]), hullFacets)
  nTriangles <- length(triFacets)
  
  Triangles <- tapply(lapply(
    triFacets, `[[`, "edges"
  ), seq_along(triFacets), function(x){
    id <- attr(x[[1L]], "id")
    list(
      facetid = id,
      edges   = polygonize(x[[1L]])
    )
  }, simplify = FALSE)
  
  Vertices <- lapply(unname(hull[["vertices"]]), `[`, c("id", "point"))
  vertices <- t(vapply(Vertices, `[[`, numeric(3L), "point"))
  rownames(vertices) <- as.character(vapply(Vertices, `[[`, integer(1L), "id"))
  
  out <- list(
    vertices    = vertices,
    triangles   = unname(Triangles),
    otherfacets = unname(Other)
  )
  attr(out, "facets") <- sprintf(
    "%d triangular facet%s, %d other facet%s",
    nTriangles, ifelse(nTriangles > 1L, "s", ""),
    nOtherFacets, ifelse(nOtherFacets > 1L, "s", "")
  )
  
  out
}