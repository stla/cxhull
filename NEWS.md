# cxhull 0.4.0.9001

* New function `cxhullEdges`, which computes only the vertices and the edges 
of a convex hull.


# cxhull 0.4.0

* New function `hullSummary`, which returns a summary of a 3D triangulated 
convex hull.

* Now, the `facesColor` argument of `plotConvexHull3d` can be a vector, allowing 
to assign different colors to the faces. See the README file for examples.

* The `plotConvexHull3d` function has a new argument `palette`; it allows to 
decorate the faces with a color gradient.

* Updated the README file.


# cxhull 0.3.0

* New function `plotConvexHull3d`, to easily plot a 3D convex hull, with a 
correct orientation of the facets.

* A couple of other new functions, e.g. `EdgesAB`, which returns a matrix 
representing the edges and indicating which edges are on the border of the 
convex hull. These functions are used by `plotConvexHull3d`, they are not 
of great interest otherwise.


# cxhull 0.2.0

* New component of facets: `orientation`, indicating the facet orientation.


# cxhull 0.1.1

* Maintenance release.


# cxhull 0.1.0

* First release.



