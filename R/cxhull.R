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
  if(ncol(points) < 2L){
    stop("dimension must be at least 2")
  }
  if(nrow(points) <= ncol(points)){
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
  hull$volume <- 1/ncol(points) * 
    sum(sapply(hull$facets,
               function(f) crossprod(f[["center"]], f[["normal"]])) *
          sapply(hull$facets, "[[", "volume"))
  hull
}
