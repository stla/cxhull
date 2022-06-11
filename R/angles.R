#' @title Dihedral angles
#' 
#' @description Dihedral angles of a convex hull.
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3D points 
#'
#' @return A dataframe with three columns. The two first columns represent 
#'   the edges, given as a pair of vertex indices. The third column provides 
#'   the dihedral angle in degrees corresponding to the edge, that is the 
#'   angle between the two faces incident to this edge. This is useful to find 
#'   edges between two coplanar faces: if the faces are exactly coplanar then 
#'   the dihedral angle is 180, but because of numerical approximation one can 
#'   consider that there is coplanarity when the dihedral angle is greater 
#'   than 179, for example. This function is used in 
#'   \code{\link{plotConvexHull3d}} to get rid of such edges (if the user 
#'   sets a value to the argument \code{angleThreshold}).
#' @export
#'
#' @examples
#' # a cube ####
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
#' hull <- cxhull(points)
#' dihedralAngles(hull)
dihedralAngles <- function(hull){
  if(!isTRUE(attr(hull, "3d"))){
    stop("Not a 3d convex hull.", call. = TRUE)
  }
  ridges <- hull[["ridges"]]
  edges <- lapply(ridges, `[[`, "vertices")
  edges <- do.call(rbind, edges)
  colnames(edges) <- c("i1", "i2")
  edges <- as.data.frame(edges)
  nedges <- nrow(edges)
  angles <- numeric(nedges)
  ridgeOf <- lapply(ridges, `[[`, "ridgeOf")
  ridgeOf <- do.call(cbind, ridgeOf)
  Facets <- hull[["facets"]]
  for(i in 1L:nedges){
    facets <- ridgeOf[, i]
    f1 <- facets[1L]
    f2 <- facets[2L]
    nrml1 <- Facets[[f1]][["normal"]]
    nrml2 <- Facets[[f2]][["normal"]]
    angles[i] <- acos(min(1, max(-1, c(crossprod(nrml1, nrml2)))))
  }
  edges[["angle"]] <- 180 - angles * 180/pi
  edges
}

