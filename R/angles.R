#' @title
#' 
#' @description 
#'
#' @param hull 
#'
#' @return
#' @export
#'
#' @examples
#' library(cxhull)
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
  edges[["angle"]] <- angles * 180/pi
  edges
}

