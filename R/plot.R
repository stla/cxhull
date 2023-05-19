#' @importFrom Rvcg vcgIsotropicRemeshing
#' @noRd
refineMesh <- function(mesh){
  vcgIsotropicRemeshing(
    vcgIsotropicRemeshing(mesh, TargetLen = 0), 
    TargetLen = 0
  )
}

#' @title Plot triangulated 3d convex hull
#' @description Plot a triangulated 3d convex hull with \strong{rgl}.
#'
#' @param hull an output of \code{\link{cxhull}} applied to 3d points and 
#'   with the option \code{triangulate=TRUE}
#' @param angleThreshold a threshold angle in degrees, typically \code{179}, 
#'   to get rid of edges between coplanar faces: edges whose corresponding 
#'   dihedral angle is greater than this threshold are removed; \code{NULL} 
#'   to use another method (see the Leonardo example)
#' @param edgesAsTubes Boolean, whether to draw the edges as tubes
#' @param verticesAsSpheres Boolean, whether to draw the vertices as spheres
#' @param palette a vector of colors to make a color gradient for the faces; 
#'   if \code{NULL}, the colors of the faces are controlled by the 
#'   \code{facesColor} argument
#' @param bias,interpolate if \code{palette} is not \code{NULL}, these arguments are 
#'   passed to \code{\link[grDevices]{colorRamp}}
#' @param g a function defined on [0, 1] and taking its values in [0, 1]; it is 
#'   composed with the function created by \code{\link[grDevices]{colorRamp}}, 
#'   based on \code{palette}
#' @param facesColor the color(s) for the faces; this argument is ignored if 
#'   the argument \code{palette} is not \code{NULL}; otherwise there are three 
#'   possibilities for \code{facesColor}: 
#'   a single color, a vector of colors with length the number of triangles, 
#'   in which case one color is assigned per triangle, or a vector of colors 
#'   with length the number of faces, after merging the triangles, in 
#'   which case one color is assigned per face; use \code{\link{hullSummary}} 
#'   to know the number of faces
#' @param edgesColor the color for the edges 
#' @param tubesRadius the radius of the tubes when \code{edgesAsTubes=TRUE}
#' @param spheresRadius the radius of the spheres when 
#'   \code{verticesAsSpheres=TRUE}
#' @param spheresColor the color of the spheres when 
#'   \code{verticesAsSpheres=TRUE}
#' @param alpha number between 0 and 1 controlling the opacity of the faces
#'
#' @return No value.
#' @export
#'
#' @importFrom rgl triangles3d cylinder3d shade3d lines3d spheres3d as.mesh3d
#' @importFrom Morpho mergeMeshes
#' @importFrom grDevices colorRamp rgb
#'
#' @examples 
#' library(cxhull)
#' library(rgl)
#' cuboctahedron <- t(cuboctahedron3d()$vb[-4L, ])
#' hull <- cxhull(cuboctahedron, triangulate = TRUE)
#' # single color ####
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotConvexHull3d(hull)
#' # gradient ####
#' open3d(windowRect = c(50, 50, 562, 562))
#' if(getRversion() < "4.1.0"){
#'   palette <- "Viridis"
#' }else{
#'   palette <- "Rocket"
#' }
#' plotConvexHull3d(hull, palette = hcl.colors(256, palette), bias = 0.5)
#' 
#' \donttest{
#' library(cxhull)
#' library(rgl)
#' # Leonardo da Vinci's 72-sided sphere ####
#' hull <- cxhull(daVinciSphere, triangulate = TRUE)
#' # there are some undesirable edges:
#' plotConvexHull3d(
#'   hull, tubesRadius = 0.07, spheresRadius = 0.1
#' )
#' # => use `angleThreshold` to get rid of these edges:
#' plotConvexHull3d(
#'   hull, angleThreshold = 179,
#'   tubesRadius = 0.07, spheresRadius = 0.1
#' )
#' }
plotConvexHull3d <- function(
    hull, angleThreshold = NULL,
    edgesAsTubes = TRUE, verticesAsSpheres = TRUE, 
    palette = NULL, bias = 1, interpolate = "linear", g = identity, 
    facesColor = "navy", edgesColor = "gold", 
    tubesRadius = 0.03, spheresRadius = 0.05, spheresColor = edgesColor,
    alpha = 1
){
  if(is.null(angleThreshold)){
    edges <- EdgesAB(hull)
    trueEdges <- edges[edges[, 3L] == "yes", c(1L, 2L)]
  }else{
    edges <- dihedralAngles(hull)
    trueEdges <- 
      as.matrix(subset(edges, angle < angleThreshold)[, c("i1", "i2")])
  }
  if(is.null(palette)){
    ncolors <- length(facesColor) 
    if(ncolors == 1L){
      triangles3d(TrianglesXYZ(hull), color = facesColor, alpha = alpha)
    }else{
      nTriangles <- length(hull[["facets"]])
      trianglesxyz <- TrianglesXYZ(hull)
      triangles <- split(trianglesxyz, gl(nTriangles, 3L))
      if(ncolors == nTriangles){
        for(i in 1L:nTriangles){
          triangles3d(
            matrix(triangles[[i]], nrow = 3L, ncol = 3L), 
            color = facesColor[i], alpha = alpha
          )
        }
      }else{
        families <- as.character(attr(trianglesxyz, "families"))
        families[is.na(families)] <- 
          paste0("NA", seq_along(which(is.na(families))))
        ufamilies <- unique(families)
        if(ncolors == length(ufamilies)){
          names(facesColor) <- ufamilies
          for(i in 1L:nTriangles){
            family <- families[i]
            triangles3d(
              matrix(triangles[[i]], nrow = 3L, ncol = 3L), 
              color = facesColor[family], alpha = alpha
            )
          }
        }else{
          warning("Invalid number of colors.")
        }
      }
    }
  }else{
    nTriangles <- length(hull[["facets"]])
    trianglesxyz <- TrianglesXYZ(hull)
    triangles <- split(trianglesxyz, gl(nTriangles, 3L))
    families <- as.character(attr(trianglesxyz, "families"))
    families[is.na(families)] <- 
      paste0("NA", seq_along(which(is.na(families))))
    ufamilies <- unique(families)
    mergedFaces <- rep(list(list()), length(ufamilies))
    names(mergedFaces) <- ufamilies
    for(i in 1L:nTriangles){
      family <- families[i]
      mesh <- as.mesh3d(matrix(triangles[[i]], nrow = 3L, ncol = 3L))
      mergedFaces[[family]] <- c(mergedFaces[[family]], list(mesh))
    }
    for(family in ufamilies){
      tomerge <- mergedFaces[[family]]
      if(length(tomerge) > 1L){
        mesh <- refineMesh(mergeMeshes(tomerge))
      }else{
        mesh <- refineMesh(tomerge[[1L]])
      }
      vertices <- mesh[["vb"]][-4L, ]
      center <- rowMeans(vertices)
      vertices <- sweep(vertices, 1L, center, `-`)
      dists <- sqrt(apply(vertices, 2L, crossprod))
      dists <- (dists - min(dists)) / diff(range(dists))
      fpalette <- colorRamp(palette, bias = bias, interpolate = interpolate)
      RGB <- fpalette(g(dists))
      colors <- rgb(RGB[, 1L], RGB[, 2L], RGB[, 3L], maxColorValue = 255)
      mesh[["material"]][["color"]] <- colors
      shade3d(mesh, alpha = alpha)
    }
  }
  Vertices <- VerticesXYZ(hull)
  for(i in 1L:nrow(trueEdges)){
    edge <- trueEdges[i, ]
    if(edgesAsTubes){
      tube <- cylinder3d(
        Vertices[edge, ], radius = tubesRadius, sides = 90
      )
      shade3d(tube, color = edgesColor)
    }else{
      lines3d(Vertices[edge, ], color = edgesColor, lwd = 2)
    }
  }
  if(verticesAsSpheres){
    spheres3d(Vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}

