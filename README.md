cxhull
================
2018-06-04

The purpose of the `cxhull` package is to compute the convex hull of a
set of points in arbitrary dimension.

The output of the `cxhull` function is a list with the following fields.

  - `vertices`: The vertices of the convex hull. Each vertex is given
    with its neighbour vertices, its neighbour ridges and its neighbour
    facets.

  - `edges`: The edges of the convex hull, given as pairs of vertices
    identifiants.

  - `ridges`: The ridges of the convex hull, i.e. the elements of the
    convex hull of dimension `dim-2`. Thus the ridges are just the
    vertices in dimension 2, and they are the edges in dimension 3.

  - `facets`: The facets of the convex hull, i.e. the elements of the
    convex hull of dimension `dim-1`. Thus the facets are the edges in
    dimension 2, and they are the faces of the convex polyhedron in
    dimension 3.

<!-- end list -->

``` r
library(cxhull)
vertices <- rbind(
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
hull <- cxhull(vertices)
```

``` r
hull$vertices
## $`2`
## $`2`$id
## [1] 2
##
## $`2`$point
## [1] 0 0 0
##
## $`2`$neighvertices
## [1] 3 4 6
##
## $`2`$neighridges
## [1] 1 2 5
##
## $`2`$neighfacets
## [1] 1 2 4
##
##
## $`3`
## $`3`$id
## [1] 3
##
## $`3`$point
## [1] 0 0 1
##
## $`3`$neighvertices
## [1] 2 5 7
##
## $`3`$neighridges
## [1]  5  7 11
##
## $`3`$neighfacets
## [1] 2 4 6
##
##
## $`4`
## $`4`$id
## [1] 4
##
## $`4`$point
## [1] 0 1 0
##
## $`4`$neighvertices
## [1] 2 5 8
##
## $`4`$neighridges
## [1]  1  4 10
##
## $`4`$neighfacets
## [1] 1 4 5
##
##
## $`5`
## $`5`$id
## [1] 5
##
## $`5`$point
## [1] 0 1 1
##
## $`5`$neighvertices
## [1] 3 4 9
##
## $`5`$neighridges
## [1] 10 11 12
##
## $`5`$neighfacets
## [1] 4 5 6
##
##
## $`6`
## $`6`$id
## [1] 6
##
## $`6`$point
## [1] 1 0 0
##
## $`6`$neighvertices
## [1] 2 7 8
##
## $`6`$neighridges
## [1] 2 3 6
##
## $`6`$neighfacets
## [1] 1 2 3
##
##
## $`7`
## $`7`$id
## [1] 7
##
## $`7`$point
## [1] 1 0 1
##
## $`7`$neighvertices
## [1] 3 6 9
##
## $`7`$neighridges
## [1] 6 7 9
##
## $`7`$neighfacets
## [1] 2 3 6
##
##
## $`8`
## $`8`$id
## [1] 8
##
## $`8`$point
## [1] 1 1 0
##
## $`8`$neighvertices
## [1] 4 6 9
##
## $`8`$neighridges
## [1] 3 4 8
##
## $`8`$neighfacets
## [1] 1 3 5
##
##
## $`9`
## $`9`$id
## [1] 9
##
## $`9`$point
## [1] 1 1 1
##
## $`9`$neighvertices
## [1] 5 7 8
##
## $`9`$neighridges
## [1]  8  9 12
##
## $`9`$neighfacets
## [1] 3 5 6
```

Edges are given as a matrix, each row representing an edge given as a
pair of vertices identifiants:

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
hull$ridges
## [[1]]
## [[1]]$id
## [1] 1
##
## [[1]]$ridgeOf
## [1] 1 4
##
## [[1]]$vertices
## [1] 2 4
##
##
## [[2]]
## [[2]]$id
## [1] 2
##
## [[2]]$ridgeOf
## [1] 1 2
##
## [[2]]$vertices
## [1] 2 6
##
##
## [[3]]
## [[3]]$id
## [1] 3
##
## [[3]]$ridgeOf
## [1] 1 3
##
## [[3]]$vertices
## [1] 6 8
##
##
## [[4]]
## [[4]]$id
## [1] 4
##
## [[4]]$ridgeOf
## [1] 1 5
##
## [[4]]$vertices
## [1] 4 8
##
##
## [[5]]
## [[5]]$id
## [1] 5
##
## [[5]]$ridgeOf
## [1] 2 4
##
## [[5]]$vertices
## [1] 2 3
##
##
## [[6]]
## [[6]]$id
## [1] 6
##
## [[6]]$ridgeOf
## [1] 2 3
##
## [[6]]$vertices
## [1] 6 7
##
##
## [[7]]
## [[7]]$id
## [1] 7
##
## [[7]]$ridgeOf
## [1] 2 6
##
## [[7]]$vertices
## [1] 3 7
##
##
## [[8]]
## [[8]]$id
## [1] 8
##
## [[8]]$ridgeOf
## [1] 3 5
##
## [[8]]$vertices
## [1] 8 9
##
##
## [[9]]
## [[9]]$id
## [1] 9
##
## [[9]]$ridgeOf
## [1] 3 6
##
## [[9]]$vertices
## [1] 7 9
##
##
## [[10]]
## [[10]]$id
## [1] 10
##
## [[10]]$ridgeOf
## [1] 4 5
##
## [[10]]$vertices
## [1] 4 5
##
##
## [[11]]
## [[11]]$id
## [1] 11
##
## [[11]]$ridgeOf
## [1] 4 6
##
## [[11]]$vertices
## [1] 3 5
##
##
## [[12]]
## [[12]]$id
## [1] 12
##
## [[12]]$ridgeOf
## [1] 5 6
##
## [[12]]$vertices
## [1] 5 9
```

The `vertices` field provides the vertices identifiants of the ridge. A
ridge is between two facets. The identifiants of these facets are given
in the field `ridgeOf`.

Facets are given as a list:

``` r
hull$facets
## [[1]]
## [[1]]$vertices
## [1] 2 4 6 8
##
## [[1]]$edges
##      [,1] [,2]
## [1,]    2    4
## [2,]    2    6
## [3,]    4    8
## [4,]    6    8
##
## [[1]]$ridges
## [1] 1 2 3 4
##
## [[1]]$neighbors
## [1] 2 3 4 5
##
## [[1]]$volume
## [1] 1
##
## [[1]]$center
## [1] 0.5 0.5 0.0
##
## [[1]]$normal
## [1]  0  0 -1
##
## [[1]]$offset
## [1] 0
##
## [[1]]$family
## [1] NA
##
##
## [[2]]
## [[2]]$vertices
## [1] 2 3 6 7
##
## [[2]]$edges
##      [,1] [,2]
## [1,]    2    3
## [2,]    2    6
## [3,]    3    7
## [4,]    6    7
##
## [[2]]$ridges
## [1] 2 5 6 7
##
## [[2]]$neighbors
## [1] 1 3 4 6
##
## [[2]]$volume
## [1] 1
##
## [[2]]$center
## [1] 0.5 0.0 0.5
##
## [[2]]$normal
## [1]  0 -1  0
##
## [[2]]$offset
## [1] 0
##
## [[2]]$family
## [1] NA
##
##
## [[3]]
## [[3]]$vertices
## [1] 6 7 8 9
##
## [[3]]$edges
##      [,1] [,2]
## [1,]    6    7
## [2,]    6    8
## [3,]    7    9
## [4,]    8    9
##
## [[3]]$ridges
## [1] 3 6 8 9
##
## [[3]]$neighbors
## [1] 1 2 5 6
##
## [[3]]$volume
## [1] 1
##
## [[3]]$center
## [1] 1.0 0.5 0.5
##
## [[3]]$normal
## [1] 1 0 0
##
## [[3]]$offset
## [1] -1
##
## [[3]]$family
## [1] NA
##
##
## [[4]]
## [[4]]$vertices
## [1] 2 3 4 5
##
## [[4]]$edges
##      [,1] [,2]
## [1,]    2    3
## [2,]    2    4
## [3,]    3    5
## [4,]    4    5
##
## [[4]]$ridges
## [1]  1  5 10 11
##
## [[4]]$neighbors
## [1] 1 2 5 6
##
## [[4]]$volume
## [1] 1
##
## [[4]]$center
## [1] 0.0 0.5 0.5
##
## [[4]]$normal
## [1] -1  0  0
##
## [[4]]$offset
## [1] 0
##
## [[4]]$family
## [1] NA
##
##
## [[5]]
## [[5]]$vertices
## [1] 4 5 8 9
##
## [[5]]$edges
##      [,1] [,2]
## [1,]    4    5
## [2,]    4    8
## [3,]    5    9
## [4,]    8    9
##
## [[5]]$ridges
## [1]  4  8 10 12
##
## [[5]]$neighbors
## [1] 1 3 4 6
##
## [[5]]$volume
## [1] 1
##
## [[5]]$center
## [1] 0.5 1.0 0.5
##
## [[5]]$normal
## [1] 0 1 0
##
## [[5]]$offset
## [1] -1
##
## [[5]]$family
## [1] NA
##
##
## [[6]]
## [[6]]$vertices
## [1] 3 5 7 9
##
## [[6]]$edges
##      [,1] [,2]
## [1,]    3    5
## [2,]    3    7
## [3,]    5    9
## [4,]    7    9
##
## [[6]]$ridges
## [1]  7  9 11 12
##
## [[6]]$neighbors
## [1] 2 3 4 5
##
## [[6]]$volume
## [1] 1
##
## [[6]]$center
## [1] 0.5 0.5 1.0
##
## [[6]]$normal
## [1] 0 0 1
##
## [[6]]$offset
## [1] -1
##
## [[6]]$family
## [1] NA
```

Observe that the `family` field is always `NA`. This field has a
possibly non-missing value only when one requires the triangulation of
the convex hull:

``` r
thull <- cxhull(vertices, triangulate = TRUE)
lapply(thull$facets, `[[`, "family")
## [[1]]
## [1] 0
##
## [[2]]
## [1] 0
##
## [[3]]
## [1] 2
##
## [[4]]
## [1] 2
##
## [[5]]
## [1] 4
##
## [[6]]
## [1] 4
##
## [[7]]
## [1] 6
##
## [[8]]
## [1] 6
##
## [[9]]
## [1] 8
##
## [[10]]
## [1] 8
##
## [[11]]
## [1] 10
##
## [[12]]
## [1] 10
```

The hull is triangulated into 12 triangles: each face of the cube is
triangulated into two triangles. Therefore one gets six different
families, each one consisting of two triangles: two triangles belong to
the same family mean that they are parts of the same facet of the
non-triangulated hull.

The triangulation is useful to plot a 3-dimensional hull:

``` r
library(rgl)
# generates 50 points in the unit sphere (uniformly)
npoints <- 50
r <- runif(npoints)^(1/3)
theta <- runif(npoints, 0, 2*pi)
v <- runif(npoints, -1, 1)
points <- t(mapply(function(r,theta,v){
  r*c(cos(theta)*sin(acos(v)), sin(theta)*sin(acos(v)), v)
}, r, theta, v))
# computes the triangulated convex hull
h <- cxhull(points, triangulate = TRUE)
# plot
for(i in 1:length(h$facets)){
  triangle <- t(sapply(h$facets[[i]]$vertices,
                       function(id) h$vertices[[as.character(id)]]$point))
  triangles3d(triangle, color="blue")
}
```

![Imgur](https://i.imgur.com/9Awcfg7.png)

## Ordonning the vertices

Observe the vertices of the first face of the cube:

``` r
hull$facets[[1]]$vertices
## [1] 2 4 6 8
```

They are given as `2-4-6-8`. They are not ordered, in the sense that
`4-6` and `8-2` are not edges of this face:

``` r
( face_edges <- hull$facets[[1]]$edges )
##      [,1] [,2]
## [1,]    2    4
## [2,]    2    6
## [3,]    4    8
## [4,]    6    8
```

To order the vertices, one can use the `igraph` package as follows:

``` r
library(igraph)
g <- graph_from_edgelist(apply(face_edges,2,as.character))
biconnected_components(g)$components[[1]]
## + 4/4 vertices, named, from d9f762c:
## [1] 6 8 4 2
```

## Volume

To get the volume of the convex hull (area in dimension 2, volume in
dimension 3, hypervolume in higer dimension), do:

``` r
1/d * sum(sapply(hull$facets,
                 function(f) crossprod(f[["center"]], f[["normal"]])) *
            sapply(hull$facets, "[[", "volume"))
```

where `d` is the dimension.
