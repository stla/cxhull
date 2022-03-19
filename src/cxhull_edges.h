/* author: Stéphane Laurent */
#ifndef CXHULLEDGESHEADER
#define CXHULLEDGESHEADER

#ifndef NULL
  #define NULL ((void *)0)
#endif

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

typedef struct SimpleRidge {
  SimpleSiteT*  vertices;
  unsigned      ridgeOf1;
  unsigned      ridgeOf2;
  unsigned      nvertices;
  unsigned      cantorid;
} SimpleRidgeT;

#endif