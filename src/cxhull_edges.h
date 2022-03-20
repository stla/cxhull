/* author: Stéphane Laurent */
#ifndef CXHULLEDGESHEADER
#define CXHULLEDGESHEADER

#ifndef NULL // just to get rid of the RStudio warnings
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

typedef struct SimpleHull {
  SiteT*     sites;
  unsigned   nsites;
  unsigned** edges;
  unsigned   nedges;
} SimpleHullT;

typedef struct SimpleRidge {
  SimpleSiteT*  vertices;
  unsigned      ridgeOf1;
  unsigned      ridgeOf2;
  unsigned      nvertices;
  unsigned      cantorid;
} SimpleRidgeT;

#endif