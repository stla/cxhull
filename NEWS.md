# cxhull 0.7.4

- The faces of the mesh returned by the `hullMesh` function are now given by 
the vertex indices of this mesh, while in the previous versions they were 
given by the vertex indices of the original vertices (those for which the 
convex hull is computed).

- The package does no longer depend on the 'Morpho' package.


# cxhull 0.7.3

New argument `alpha` in the `plotConvexHull3d` function, to control the 
transparency of the faces. 


# cxhull 0.7.2

Fixed a small mistake in the C code.


# cxhull 0.7.1

Replaced `sprintf` with `snprintf` in the C code.


# cxhull 0.7.0

* New function `hullMesh` to extract the vertices and the faces of a 3d convex 
hull.

* Updated the vertices of Leonardo da Vinci's 72-sided sphere; the new ones 
are more accurate.


# cxhull 0.6.0

* New function `dihedralAngles` which computes the dihedral angles of a 3d 
convex hull.

* New argument `angleThreshold` in the `plotConvexHull3d` function. Edges whose
corresponding dihedral angle is greater than `angleThreshold` are removed from 
the plot. 


# cxhull 0.5.0

* New function `cxhullEdges`, which computes only the vertices and the edges 
of a convex hull (for speed gain and less memory consumption).


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



