/* author: Stéphane Laurent */
#include "Rheaders.h"
#include "cxhull_edges.h"
#include "qhull_ra.h"
#include "utils.h"

// to use the qsort function - sort vertices according to their ids ------------
int cmpsites(const void* a, const void* b) {
  return ((*((SiteT*)a)).id - (*((SiteT*)b)).id);
}

// main function ---------------------------------------------------------------
SetOfSitesT cxhullEdges(double* points,
                        unsigned dim,
                        unsigned n,
                        unsigned* exitcode,
                        const char* errfilename) {
  char opts[] = "qhull ";  // option flags for qhull, see qh_opt.htm
  qhT qh_qh;               // Qhull's data structure
  qhT* qh = &qh_qh;
  QHULL_LIB_CHECK
  boolT ismalloc = False;  // True if qhull should free points in qh_freeqhull()
  // or reallocation
  FILE* errfile = fopen(errfilename, "w+");  // NULL;
  FILE* outfile = NULL;
  qh_meminit(qh, errfile);
  qh_zero(qh, errfile);
  exitcode[0] =
      qh_new_qhull(qh, dim, n, points, ismalloc, opts, outfile, errfile);
  fclose(errfile);

  if(!exitcode[0]) {  // 0 if no error from qhull

    printf("no qhull error");
    // all vertices
    unsigned nvertices = qh->num_vertices;
    SiteT* vertices = malloc(nvertices * sizeof(SiteT));
    {
      qh_vertexneighbors(qh);  // make the neighbor facets of the vertices
      vertexT* vertex;
      unsigned i_vertex = 0;
      FORALLvertices {
        // vertex id and coordinates
        unsigned vertex_id = qh_pointid(qh, vertex->point);
        printf("vertexid: %d", vertex_id);
        vertices[i_vertex].id = vertex_id;
        vertices[i_vertex].point = getpoint(points, dim, vertex_id);
        // vertices[i_vertex].nneighsvertices = 0;

        {
          unsigned* neighs = malloc(0);
          unsigned nneighs = 0;
          facetT *neighbor, **neighborp;
          unsigned i_neighbor = 0;
          FOREACHneighbor_(vertex) {
            // vertices[i_vertex].neighfacets[i_neighbor] = neighbor->id;
            printf("neighbor i: %d\n", i_neighbor);
            i_neighbor++;
            if(neighbor->simplicial) {
              printf("neighbor simplicial\n");
              // vertices[i_vertex].nneighsvertices += qh_setsize(qh,
              // neighbor->vertices);
              vertexT *v, **vp;
              FOREACHsetelement_(vertexT, neighbor->vertices, v) {
                unsigned v_id = qh_pointid(qh, v->point);
                if(v_id != vertex_id) {
                  unsigned pushed;
                  appendu(v_id, &neighs, nneighs, &pushed);
                  if(pushed) {
                    printf("pushed");
                    nneighs++;
                  }
                }
              }
            } else {
              qh_makeridges(qh, neighbor);
              // vertexT *v, **vp;
              // FOREACHsetelement_(vertexT, neighbor->vertices, v) {
              ridgeT *r, **rp;
              FOREACHsetelement_(ridgeT, qh_vertexridges(qh, vertex), r) {
                vertexT *v, **vp;
                FOREACHsetelement_(vertexT, r->vertices, v) {
                  unsigned v_id = qh_pointid(qh, v->point);
                  if(v_id != vertex_id) {
                    unsigned pushed;
                    appendu(v_id, &neighs, nneighs, &pushed);
                    if(pushed) {
                      printf("pushed u");
                      nneighs++;
                    }
                  }
                }
              }
              // }
            }
          }
          vertices[i_vertex].neighvertices = neighs;
          vertices[i_vertex].nneighvertices = nneighs;
          i_vertex++;
        }
      }
    }

    qsort(vertices, nvertices, sizeof(SiteT), cmpsites);

    int curlong, totlong;
    qh_freeqhull(qh, !qh_ALL);  // free long memory
    qh_memfreeshort(qh, &curlong,
                    &totlong);  // free short memory and memory allocator
    SetOfSitesT vset = {.sites = vertices, .nsites = nvertices};

    return vset;
  } else {
    int curlong, totlong;
    qh_freeqhull(qh, !qh_ALL);  // free long memory
    qh_memfreeshort(qh, &curlong,
                    &totlong);  // free short memory and memory allocator
    error("Received error code %d from qhull.", exitcode[0]);
  }
}

// -------------------------------------------------------------------------- //
// ----------------------------------- R ------------------------------------ //
// -------------------------------------------------------------------------- //

// SiteT to SEXP ---------------------------------------------------------
SEXP SiteSXP(SiteT vertex, unsigned dim) {
  unsigned nprotect = 0;
  SEXP R_vertex, names, id, point, neighvertices;

  PROTECT(id = allocVector(INTSXP, 1));
  nprotect++;
  INTEGER(id)[0] = 1 + vertex.id;

  PROTECT(point = allocVector(REALSXP, dim));
  nprotect++;
  for(int i = 0; i < dim; i++) {
    REAL(point)[i] = vertex.point[i];
  }

  unsigned nneighvertices = vertex.nneighvertices;
  PROTECT(neighvertices = allocVector(INTSXP, nneighvertices));
  nprotect++;
  for(unsigned i = 0; i < nneighvertices; i++) {
    INTEGER(neighvertices)[i] = 1 + vertex.neighvertices[i];
  }

  PROTECT(R_vertex = allocVector(VECSXP, 3));
  nprotect++;
  SET_VECTOR_ELT(R_vertex, 0, id);
  SET_VECTOR_ELT(R_vertex, 1, point);
  SET_VECTOR_ELT(R_vertex, 2, neighvertices);

  PROTECT(names = allocVector(VECSXP, 3));
  nprotect++;
  SET_VECTOR_ELT(names, 0, mkChar("id"));
  SET_VECTOR_ELT(names, 1, mkChar("point"));
  SET_VECTOR_ELT(names, 2, mkChar("neighvertices"));
  setAttrib(R_vertex, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return R_vertex;
}

// main function ---------------------------------------------------------------
SEXP cxhullEdges_(SEXP p, SEXP errfile) {
  unsigned nprotect = 0;

  unsigned dim = ncols(p);
  unsigned n = nrows(p);

  double* points = (double*)R_alloc(n * dim, sizeof(double));
  for(unsigned i = 0; i < n; i++)
    for(unsigned j = 0; j < dim; j++)
      points[dim * i + j] = REAL(
          p)[i + n * j];  // could have been REAL(p) if p had been transposed

  unsigned exitcode;
  const char* e = R_CHAR(Rf_asChar(errfile));
  SetOfSitesT vset = cxhullEdges(points, dim, n, &exitcode, e);
  if(exitcode) {
    error("Received error code %d from qhull.", exitcode);
  }

  unsigned nvertices = vset.nsites;
  SiteT* vertices = vset.sites;

  SEXP R_vertices, vnames;

  PROTECT(R_vertices = allocVector(VECSXP, nvertices));
  PROTECT(vnames = allocVector(STRSXP, nvertices));
  nprotect += 2;
  for(unsigned i = 0; i < nvertices; i++) {
    SEXP vertex;
    PROTECT(vertex = SiteSXP(vertices[i], dim));
    nprotect++;
    SET_VECTOR_ELT(R_vertices, i, vertex);
    SET_STRING_ELT(vnames, i, Rf_asChar(VECTOR_ELT(vertex, 0)));
  }
  setAttrib(R_vertices, R_NamesSymbol, vnames);

  UNPROTECT(nprotect);
  return R_vertices;
}