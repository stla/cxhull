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
#'   of a list if possible, i.e. if all faces have the same number of edges;
#'   this argument is possibly ignored if \code{rgl=TRUE}, see below
#' @param rgl Boolean, whether to return a \strong{rgl} mesh (class 
#'   \code{mesh3d}) if possible, i.e. if each face has three or four edges; 
#'   if \code{TRUE} and the \strong{rgl} mesh is possible, then the 
#'   \code{simplify} argument has no effect
#'
#' @return A list giving the vertices and the faces, or a \strong{rgl} mesh.
#' @export
#' 
#' @importFrom data.table uniqueN
#' @importFrom rgl mesh3d
#' 
#' @note Unless all faces are triangular, the output does not define a mesh 
#'   with coherently oriented faces. 
#'
#' @examples
#' library(cxhull)
#' hull <- cxhull(daVinciSphere)
#' septuaginta <- hullMesh(hull, rgl = TRUE)
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(septuaginta, color = "darkred")
#' # some quad faces are misoriented:
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(septuaginta, color = "tomato", back = "culled")
hullMesh <- function(hull, simplify = TRUE, rgl = FALSE) {
  stopifnot(isBoolean(simplify))
  stopifnot(isBoolean(rgl))
  vertices <- VerticesXYZ(hull)
  faces <- lapply(hull[["facets"]], orderFace)
  # transform indices to get vertex indices of the mesh
  dict        <- 1L:nrow(vertices)
  names(dict) <- rownames(vertices)
  faces <- lapply(faces, function(x) dict[as.character(x)])
  #
  if(rgl) {
    nedges <- lengths(faces)
    if(all(nedges %in% c(3L, 4L))) {
      sfaces <- split(faces, nedges)
      trgls <- sfaces[["3"]]
      quads <- sfaces[["4"]]
      if(!is.null(trgls)) {
        trgls <- do.call(cbind, trgls)
      }
      if(!is.null(quads)) {
        quads <- do.call(cbind, quads)
      }
      rmesh <- mesh3d(vertices, triangles = trgls, quads = quads)
      return(rmesh)
    } else {
      warning(
        "Cannot do a 'rgl' mesh (there are faces with five edges or more)."
      )
    }
  }
  if(simplify) {
    N <- uniqueN(lengths(faces))
    if(N == 1L) {
      faces <- do.call(rbind, faces)
    } else {
      message(
        "Cannot stack the faces to a matrix - ignoring `simplify=TRUE`."
      )
    }
  }
  list("vertices" = vertices, "faces" = faces)
}
