cxhull
================
2019-03-09

The purpose of the `cxhull` package is to compute the convex hull of a
set of points in arbitrary dimension. It contains only one function:
`cxhull`.

The output of the `cxhull` function is a list with the following fields.

  - `vertices`: The vertices of the convex hull. Each vertex is given
    with its neighbour vertices, its neighbour ridges and its neighbour
    facets.

  - `edges`: The edges of the convex hull, given as pairs of vertices
    identifiers.

  - `ridges`: The ridges of the convex hull, i.e. the elements of the
    convex hull of dimension `dim-2`. Thus the ridges are just the
    vertices in dimension 2, and they are the edges in dimension 3.

  - `facets`: The facets of the convex hull, i.e. the elements of the
    convex hull of dimension `dim-1`. Thus the facets are the edges in
    dimension 2, and they are the faces of the convex polyhedron in
    dimension 3.

  - `volume`: The volume of the convex hull (area in dimension 2, volume
    in dimension 3, hypervolume in higher dimension).

Let’s look at an example. The points we take are the vertices of a cube
and the center of this cube (in the first row):

``` r
library(cxhull)
points <- rbind(
  c(0.5, 0.5, 0.5),
  c(0, 0, 0),
  c(0, 0, 1),
  c(0, 1, 0),
  c(0, 1, 1),
  c(1, 0, 0),
  c(1, 0, 1),
  c(1, 1, 0),
  c(1, 1, 1)
)
hull <- cxhull(points)
```

Obviously, the convex hull of these points is the cube. We can quickly
see that the convex hull has 8 vertices, 12 edges, 12 ridges, 6 facets,
and its volume is 1:

``` r
str(hull, max=1)
## List of 5
##  $ vertices:List of 8
##  $ edges   : int [1:12, 1:2] 2 2 2 3 3 4 4 5 6 6 ...
##  $ ridges  :List of 12
##  $ facets  :List of 6
##  $ volume  : num 1
```

Each vertex, each ridge, and each facet has an identifier. A vertex
identifier is the index of the row corresponding to this vertex in the
set of points passed to the `cxhull` function. It is given in the field
`id` of the vertex:

``` r
hull$vertices[[1]]
## $id
## [1] 2
## 
## $point
## [1] 0 0 0
## 
## $neighvertices
## [1] 3 4 6
## 
## $neighridges
## [1] 1 2 5
## 
## $neighfacets
## [1] 1 2 4
```

Also, the list of vertices is named with the identifiers:

``` r
names(hull$vertices)
## [1] "2" "3" "4" "5" "6" "7" "8" "9"
```

Edges are given as a matrix, each row representing an edge given as a
pair of vertices identifiers:

``` r
hull$edges
##       [,1] [,2]
##  [1,]    2    3
##  [2,]    2    4
##  [3,]    2    6
##  [4,]    3    5
##  [5,]    3    7
##  [6,]    4    5
##  [7,]    4    8
##  [8,]    5    9
##  [9,]    6    7
## [10,]    6    8
## [11,]    7    9
## [12,]    8    9
```

The ridges are given as a list:

``` r
hull$ridges[[1]]
## $id
## [1] 1
## 
## $ridgeOf
## [1] 1 4
## 
## $vertices
## [1] 2 4
```

The `vertices` field provides the vertices identifiers of the ridge. A
ridge is between two facets; the identifiers of these facets are given
in the field `ridgeOf`.

Facets are given as a list:

``` r
hull$facets[[1]]
## $vertices
## [1] 8 4 6 2
## 
## $edges
##      [,1] [,2]
## [1,]    2    4
## [2,]    2    6
## [3,]    4    8
## [4,]    6    8
## 
## $ridges
## [1] 1 2 3 4
## 
## $neighbors
## [1] 2 3 4 5
## 
## $volume
## [1] 1
## 
## $center
## [1] 0.5 0.5 0.0
## 
## $normal
## [1]  0  0 -1
## 
## $offset
## [1] 0
## 
## $family
## [1] NA
## 
## $orientation
## [1] 1
```

There is no `id` field for the facets: the integer `i` is the identifier
of the `i`-th facet of the list.

The `orientation` field has two possible values, `1` or `-1`, it
indicates the orientation of the facet. See the plotting example below.

Here, the `family` field is `NA` for every facet:

``` r
sapply(hull$facets, `[[`, "family")
## [1] NA NA NA NA NA NA
```

This field has a possibly non-missing value only when one requires the
triangulation of the convex hull:

``` r
thull <- cxhull(points, triangulate = TRUE)
sapply(thull$facets, `[[`, "family")
##  [1]  0  0  2  2  4  4  6  6  8  8 10 10
```

The hull is triangulated into 12 triangles: each face of the cube is
triangulated into two triangles. Therefore one gets six different
families, each one consisting of two triangles: two triangles belong to
the same family mean that they are parts of the same facet of the
non-triangulated hull.

## Plotting a 3-dimensional hull

The triangulation is useful to plot a 3-dimensional hull:

``` r
library(rgl)
# generates 50 points in the unit sphere (uniformly):
set.seed(666)
npoints <- 50
r <- runif(npoints)^(1/3)
theta <- runif(npoints, 0, 2*pi)
cosphi <- runif(npoints, -1, 1)
points <- t(mapply(function(r,theta,cosphi){
  r * c(cos(theta)*sin(acos(cosphi)), sin(theta)*sin(acos(cosphi)), cosphi)
}, r, theta, cosphi))
# computes the triangulated convex hull:
h <- cxhull(points, triangulate = TRUE)
# plot:
open3d(windowRect = c(50,50,450,450))
view3d(10, 80, zoom = 0.75)
for(i in 1:length(h$facets)){
  triangle <- t(sapply(h$facets[[i]]$vertices,
                       function(id) h$vertices[[as.character(id)]]$point))
  triangles3d(triangle, color = "blue")
}
```

![Imgur](https://i.imgur.com/mlbOsco.png)

## Facets orientation

The facets can be differently oriented. This can be seen if we hide the
back side of the triangles:

``` r
open3d(windowRect = c(50,50,450,450))
view3d(10, 80, zoom = 0.75)
for(i in 1:length(h$facets)){
  triangle <- t(sapply(h$facets[[i]]$vertices,
                       function(id) h$vertices[[as.character(id)]]$point))
  triangles3d(triangle, color = "blue", back = "culled")
}
```

![Imgur](https://i.imgur.com/xLSvGtm.png)

The `orientation` field of a facet indicates its orientation (`1` or
`-1`). We can use it to change the facet orientation:

``` r
open3d(windowRect = c(50,50,450,450))
view3d(10, 80, zoom = 0.75)
for(i in 1:length(h$facets)){
  facet <- h$facets[[i]]
  triangle <- t(sapply(facet$vertices,
                       function(id) h$vertices[[as.character(id)]]$point))
  if(facet$orientation == 1){
    triangles3d(triangle, color="blue", back="culled")
  }else{
    triangles3d(triangle[c(1,3,2),], color="red", back="culled")
  }
}
```

![Imgur](https://i.imgur.com/xKr48HX.png)

## Ordering the vertices

Observe the vertices of the first face of the cube:

``` r
hull$facets[[1]]$vertices
## [1] 8 4 6 2
```

They are given as `8-4-6-2`. They are not ordered, in the sense that
`4-6` and `2-8` are not edges of this face:

``` r
( face_edges <- hull$facets[[1]]$edges )
##      [,1] [,2]
## [1,]    2    4
## [2,]    2    6
## [3,]    4    8
## [4,]    6    8
```

One can order the vertices as follows:

``` r
polygon <- function(edges){
  nedges <- nrow(edges)
  vs <- edges[1,]
  v <- vs[2]
  edges <- edges[-1,]
  for(. in 1:(nedges-2)){
    j <- which(apply(edges, 1, function(e) v %in% e))
    v <- edges[j,][which(edges[j,]!=v)]
    vs <- c(vs,v)
    edges <- edges[-j,]
  }
  vs
}
polygon(face_edges)
## [1] 2 4 8 6
```

See a four-dimensional example in [the
vignette](https://cran.r-project.org/web/packages/cxhull/vignettes/truncated_tesseract.html).
