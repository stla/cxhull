orderFace <- function(face) {
  ids <- face[["vertices"]]
  if(length(ids) == 3L) {
    if(face[["orientation"]] == -1L) {
      ids <- rev(ids)
    }
  } else {
    edges <- face[["edges"]]
    nedges <- nrow(edges)
    ids <- edges[1L, ]
    v <- ids[2L]
    edges <- edges[-1L, ]
    for(. in 1L:(nedges-2L)){
      j <- which(apply(edges, 1L, function(e) v %in% e))
      v <- edges[j, ][which(edges[j, ] != v)]
      ids <- c(ids, v)
      edges <- edges[-j, ]
    }
  }
  ids
}

#' @title Mesh of a 3d convex hull
#' @description Extract the vertices and the faces from a 3d convex hull.
#'
#' @param hull a 3d convex hull, output of \code{\link{cxhull}}
#' @param simplify Boolean, whether to return the faces as a matrix instead 
#'   of a list if possible, i.e. when all faces have the same number of sides
#'
#' @return A list giving the vertices and the faces.
#' @export
#' 
#' @note Unless all faces are triangular, the output does not define a mesh 
#'   with coherently oriented faces. If you want correct face orientations, 
#'   you can use the \strong{PolygonSoup} package. 
#'
#' @examples
#' library(cxhull)
#' hull <- cxhull(daVinciSphere)
#' septuaginta <- hullMesh(hull)
hullMesh <- function(hull, simplify = TRUE) {
  stopifnot(isBoolean(simplify))
  vertices <- VerticesXYZ(hull)
  faces <- lapply(hull[["facets"]], orderFace)
  if(simplify) {
    nsides <- uniqueN(lengths(faces))
    if(nsides == 1L) {
      faces <- do.call(rbind, faces)
    }
  }
  list("vertices" = vertices, "faces" = faces)
}
