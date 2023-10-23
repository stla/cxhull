utils::globalVariables("angle")

.onLoad <- function(lib, pkg) {
  loadNamespace("rgl") # to use merge.mesh3d
}