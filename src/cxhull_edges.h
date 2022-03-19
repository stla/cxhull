/* author: Stéphane Laurent */
#ifndef CXHULLEDGESHEADER
#define CXHULLEDGESHEADER

typedef struct Site {
  unsigned  id;
  double*   point;
  unsigned* neighvertices;
  unsigned  nneighvertices;
} SiteT;

typedef struct SimpleSite {
  unsigned  id;
  double*   point;
} SimpleSiteT;

typedef struct SetOfSites {
  SiteT*     sites;
  unsigned   nsites;
  unsigned** edges;
  unsigned   nedges;
} SetOfSitesT;

typedef struct Ridge2 {
  SimpleSiteT*  vertices;
  unsigned      ridgeOf1;
  unsigned      ridgeOf2;
  unsigned      nvertices;
  unsigned      id;
} Ridge2T;

#endif