/* author: Stéphane Laurent */
#ifndef CXHULLEDGESHEADER
#define CXHULLEDGESHEADER

typedef struct Site {
  unsigned  id;
  double*   point;
  unsigned* neighvertices;
  unsigned  nneighvertices;
} SiteT;

typedef struct SetOfSites {
  SiteT*    sites;
  unsigned  nsites;
} SetOfSitesT;

#endif