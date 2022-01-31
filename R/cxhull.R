#' Convex hull
#' @description Computes the convex hull of a set of points.
#' @param points numeric matrix, one point per row
#' @param triangulate logical, whether to triangulate the convex hull
#' @return A list.
#' @export
#' @useDynLib cxhull, .registration = TRUE
#' @examples
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
cxhull <- function(points, triangulate=FALSE){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("`points` must be a numeric matrix")
  }
  dimension <- ncol(points)
  if(dimension < 2L){
    stop("dimension must be at least 2")
  }
  if(nrow(points) <= dimension){
    stop("insufficient number of points")
  }
  if(any(is.na(points))){
    stop("missing values are not allowed")
  }
  errfile <- tempfile(fileext=".txt")
  hull <- tryCatch({
    .Call("cxhull_", points, as.integer(triangulate), errfile)
  }, error = function(e){
    cat(readLines(errfile), sep="\n")
    stop(e)
  })
  hull$volume <- 1/dimension * 
    sum(sapply(hull$facets,
               function(f) crossprod(f[["center"]], f[["normal"]])) *
          sapply(hull$facets, "[[", "volume"))
  if(dimension == 3L){
    attr(hull, "3d") <- TRUE
  }
  if(triangulate){
    attr(hull, "triangulate") <- TRUE
  }
  hull
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
      "You didn't compute the convex hull with the option `triangulate=TRUE`"
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
  for(i in 1L:ntriangles){
    facet <- Facets[[i]]
    vertices <- as.character(facet[["vertices"]])
    if(facet[["orientation"]] == -1L){
      vertices <- vertices[c(1L, 3L, 2L)]
    }
    triangles[(3L*i-2L):(3L*i), ] <- Vertices[vertices, ]
  }
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
      "You didn't compute the convex hull with the option `triangulate=TRUE`"
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
