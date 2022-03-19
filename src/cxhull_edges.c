/* author: Stéphane Laurent */
#include "Rheaders.h"
#define qh_QHimport
#include "cxhull_edges.h"
#include "qhull_ra.h"
#include "utils.h"

// to use the qsort function - sort vertices according to their ids ------------
int cmpsites(const void* a, const void* b) {
  return ((*((SiteT*)a)).id - (*((SiteT*)b)).id);
}

// all edges from all vertices -------------------------------------------------
unsigned** getEdges1(SiteT* vertices, unsigned nvertices, unsigned outlength) {
  unsigned** out = malloc(outlength * sizeof(unsigned*));
  SiteT v0 = vertices[0];
  unsigned id0 = v0.id;
  unsigned* neighs_v0 = v0.neighvertices;
  unsigned n = v0.nneighvertices;
  for(unsigned i = 0; i < n; i++) {
    out[i] = malloc(2 * sizeof(unsigned));
    out[i][0] = id0;
    out[i][1] = neighs_v0[i];
  }
  for(unsigned k = 1; k < nvertices; k++) {
    SiteT vk = vertices[k];
    unsigned ids0 = vk.id;
    unsigned* neighs_vk = vk.neighvertices;
    unsigned nk = vk.nneighvertices;
    for(unsigned i = 0; i < nk; i++) {
      unsigned ids[2];
      unsigned ids1 = neighs_vk[i];
      if(1) {  //}ids1 > ids0) {
        ids[0] = ids0;
        ids[1] = ids1;
        qsortu(ids, 2);
        unsigned j;
        for(j = 0; j < n; j++) {
          if(ids[0] == out[j][0] && ids[1] == out[j][1]) {
            break;
          }
        }
        if(j == n) {
          out[n] = malloc(2 * sizeof(unsigned));
          out[n][0] = ids[0];
          out[n][1] = ids[1];
          n++;
        }
        if(n == outlength) {
          break;
        }
      }
    }
    if(n == outlength) {
      break;
    }
  }
  return out;
}

unsigned** getEdges2(SiteT* vertices, unsigned nvertices, unsigned outlength) {
  unsigned** out = malloc(outlength * sizeof(unsigned*));
  SiteT v0 = vertices[0];
  unsigned id0 = v0.id;
  unsigned* neighs_v0 = v0.neighvertices;
  unsigned n0 = v0.nneighvertices;
  for(unsigned i = 0; i < n0; i++) {
    out[i] = malloc(2 * sizeof(unsigned));
    out[i][0] = id0;
    out[i][1] = neighs_v0[i];
  }
  for(unsigned k = 1; k < nvertices; k++) {
    SiteT vk = vertices[k];
    unsigned idk = vk.id;
    unsigned* neighs_vk = vk.neighvertices;
    unsigned nk = vk.nneighvertices;
    for(unsigned i = 0; i < nk; i++) {
      unsigned j = n0 + i;
      out[j] = malloc(2 * sizeof(unsigned));
      out[j][0] = idk;
      out[j][1] = neighs_vk[i];
    }
    n0 = n0 + nk;
  }
  return out;
}

unsigned makeSites1(qhT* qh, SiteT* vertices, double* points, unsigned dim) {
  unsigned nalledgesdouble = 0;
  {
    vertexT* vertex;
    unsigned i_vertex = 0;
    FORALLvertices {
      // vertex id and coordinates
      unsigned vertex_id = qh_pointid(qh, vertex->point);
      vertices[i_vertex].id = vertex_id;
      vertices[i_vertex].point = getpoint(points, dim, vertex_id);
      {
        unsigned* neighs = malloc(0);
        unsigned nneighs = 0;
        facetT *neighbor, **neighborp;
        FOREACHneighbor_(vertex) {
          if(neighbor->simplicial) {
            vertexT *v, **vp;
            FOREACHsetelement_(vertexT, neighbor->vertices, v) {
              unsigned v_id = qh_pointid(qh, v->point);
              if(v_id != vertex_id) {
                unsigned pushed;
                appendu(v_id, &neighs, nneighs, &pushed);
                if(pushed) {
                  nneighs++;
                }
              }
            }
          } else {
            qh_makeridges(qh, neighbor);
            ridgeT *r, **rp;
            FOREACHsetelement_(ridgeT, qh_vertexridges(qh, vertex), r) {
              vertexT *v, **vp;
              FOREACHsetelement_(vertexT, r->vertices, v) {
                unsigned v_id = qh_pointid(qh, v->point);
                if(v_id != vertex_id) {
                  unsigned pushed;
                  appendu(v_id, &neighs, nneighs, &pushed);
                  if(pushed) {
                    nneighs++;
                  }
                }
              }
            }
          }
        }
        vertices[i_vertex].neighvertices = neighs;
        vertices[i_vertex].nneighvertices = nneighs;
        i_vertex++;
        nalledgesdouble += nneighs;
      }
    }
  }
  return nalledgesdouble / 2;
}

// the threshold distance to detect neighbor vertices --------------------------
double ridgeThreshold(qhT* qh, ridgeT* ridge, double* point, unsigned dim) {
  unsigned nvertices = qh_setsize(qh, ridge->vertices);
  unsigned point_id = qh_pointid(qh, point);
  double dists[nvertices - 1];
  unsigned count = 0;
  vertexT *vertex, **vertexp;
  FOREACHvertex_(ridge->vertices) {
    unsigned vid = qh_pointid(qh, vertex->point);
    if(vid != point_id) {
      dists[count] = squaredDistance(vertex->point, point, dim);
      count++;
    }
  }
  qsort(dists, count - 1, sizeof(double), cmpfuncdbl);
  return dists[1];
}

// neighbor vertices of a vertex for dim>2 --------------------
unsigned* neighVertices2(qhT* qh,
                         vertexT* vertex,
                         unsigned dim,
                         unsigned* lengthout) {
  unsigned* neighs = malloc(0);
  *lengthout = 0;
  unsigned vertex_id = qh_pointid(qh, vertex->point);
  ridgeT *ridge, **ridgep;
  FOREACHridge_(qh_vertexridges(qh, vertex)) {
    // unsigned nvertices = qh_setsize(qh, ridge->vertices);
    vertexT *vert, **vertp;
    FOREACHsetelement_(vertexT, ridge->vertices, vert) {
      unsigned vert_id = qh_pointid(qh, vert->point);
      // unsigned v = SETindex_(ridge->vertices, vert);
      // if(vert_id == allridges[e].vertices[v].id) {
      // for(unsigned w = 0; w < allridges[e].nvertices; w++) {
      // for dim 3 needless: only two connected vertices
      if(vert_id != vertex_id &&
         (dim == 3 || squaredDistance(vert->point, vertex->point, dim) <=
                          ridgeThreshold(qh, ridge, vertex->point, dim))) {
        unsigned pushed;
        appendu(vert_id, &neighs, *lengthout, &pushed);
        if(pushed) {
          (*lengthout)++;
        }
      }
    }
    // break;
  }
  return neighs;
}

unsigned makeSites2(qhT* qh, SiteT* vertices, double* points, unsigned dim) {
  unsigned nalledges = 0;
  {
    vertexT* vertex;
    unsigned i_vertex = 0;
    FORALLvertices {
      // vertex id and coordinates
      unsigned vertex_id = qh_pointid(qh, vertex->point);
      vertices[i_vertex].id = vertex_id;
      vertices[i_vertex].point = getpoint(points, dim, vertex_id);
      //{
      // unsigned* neighs = malloc(0);
      unsigned nneighs;

      // facetT *neighbor, **neighborp;
      // FOREACHneighbor_(vertex) {
      //   if(neighbor->simplicial) {
      //     vertexT *v, **vp;
      //     FOREACHsetelement_(vertexT, neighbor->vertices, v) {
      //       unsigned v_id = qh_pointid(qh, v->point);
      //       if(v_id > vertex_id) {
      //         unsigned pushed;
      //         appendu(v_id, &neighs, nneighs, &pushed);
      //         if(pushed) {
      //           nneighs++;
      //           nalledges++;
      //         }
      //       }
      //     }
      //   } else {
      // printf("non-simplicial");
      // error("STOP");
      // qh_makeridges(qh, neighbor);
      // ridgeT *ridge, **ridgep;
      // FOREACHridge_(qh_vertexridges(qh, vertex)) {
      //   vertexT *v, **vp;
      //   FOREACHsetelement_(vertexT, ridge->vertices, v) {  //
      //     unsigned v_id = qh_pointid(qh, v->point);
      //     if(v_id > vertex_id) {
      //       unsigned pushed;
      //       appendu(v_id, &neighs, nneighs, &pushed);
      //       if(pushed) {
      //         nneighs++;
      //         nalledges++;
      //       }
      //     }
      //   }
      // }

      //   }
      // }
      vertices[i_vertex].neighvertices =
          neighVertices2(qh, vertex, dim, &nneighs);
      vertices[i_vertex].nneighvertices = nneighs;
      i_vertex++;
      nalledges += nneighs;
      //}
    }
  }
  return nalledges;
}

unsigned cantorPairing(unsigned i, unsigned j){
  return (i+j)*(i+j+1)/2 + j;
}

// main function ---------------------------------------------------------------
SetOfSitesT cxhullEdges(double* points,
                        unsigned dim,
                        unsigned n,
                        unsigned adjacencies,
                        unsigned order_edges,
                        unsigned* exitcode,
                        const char* errfilename) {
  char opts[] = "qhull s FF ";  // option flags for qhull, see qh_opt.htm
  qhT qh_qh;                    // Qhull's data structure
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

    {
      facetT* facet;
      unsigned i_facet = 0;
      FORALLfacets {
        facet->id = i_facet;  // for neighbors and ridgeOf
        i_facet++;
      }
    }

    unsigned* cantorPairs = malloc(0);
    unsigned nCantorPairs = 0;
    {
      facetT* facet;
      unsigned i_facet = 0;
      FORALLfacets {
        {  // face ridges
          qh_makeridges(qh, facet);
          unsigned nridges = qh_setsize(qh, facet->ridges);
          Ridge2T* ridges = malloc(nridges * sizeof(Ridge2T));
          ridgeT *ridge, **ridgep;
          unsigned i_ridge = 0;
          FOREACHridge_(facet->ridges) {
            //ridges[i_ridge].nedges = 0;
            unsigned ridgeSize = qh_setsize(qh, ridge->vertices);  // = dim-1
            ridges[i_ridge].nvertices = ridgeSize;
            unsigned ids[ridgeSize];
            for(unsigned v = 0; v < ridgeSize; v++) {
              ids[v] =
                  qh_pointid(qh, ((vertexT*)ridge->vertices->e[v].p)->point);
            }
            qsortu(ids, ridgeSize);
            ridges[i_ridge].vertices = malloc(ridgeSize * sizeof(SimpleSiteT));
            for(unsigned v = 0; v < ridgeSize; v++) {
              ridges[i_ridge].vertices[v].id = ids[v];
              ridges[i_ridge].vertices[v].point = getpoint(points, dim, ids[v]);
            }
            unsigned ridgeofs[2];
            ridgeofs[0] = ridge->bottom->id;
            ridgeofs[1] = ridge->top->id;
            qsortu(ridgeofs, 2);
            ridges[i_ridge].ridgeOf1 = ridgeofs[0];
            ridges[i_ridge].ridgeOf2 = ridgeofs[1];
            unsigned cantorPair = cantorPairing(ridgeofs[0], ridgeofs[1]);
            unsigned pushed;
            appendu(cantorPair, &cantorPairs, nCantorPairs, &pushed);
            if(pushed) {
              nCantorPairs++;
            }
            
            ////
            i_ridge++;
          }
          // merge triangulated ridges
          if(dim > 3 && !triangulate) {
            unsigned l;
            faces[i_facet].ridges = mergeRidges(ridges, nridges, &l);
            faces[i_facet].nridges = l;
          } else {  // dim 2 or 3, or triangulate option
            faces[i_facet].ridges = ridges;
            faces[i_facet].nridges = nridges;
          }
        }
        ////
        i_facet++;
      }
    }

    qh_vertexneighbors(qh);  // make the neighbor facets of the vertices
    unsigned nvertices = qh->num_vertices;
    SiteT* vertices = malloc(n * sizeof(SiteT));
    unsigned nedges;
    unsigned** edges;
    if(adjacencies) {
      nedges = makeSites1(qh, vertices, points, dim);
      printf("nedges adj: %d\n", nedges);
      qsort(vertices, nvertices, sizeof(SiteT), cmpsites);
      edges = getEdges1(vertices, nvertices, nedges);
    } else {
      nedges = makeSites2(qh, vertices, points, dim);
      printf("nedges xxx: %d\n", nedges);
      qsort(vertices, n, sizeof(SiteT), cmpsites);
      edges = getEdges2(vertices, n, nedges);
    }
    //{
    // vertexT* vertex;
    // unsigned i_vertex = 0;
    // FORALLvertices {
    //   // vertex id and coordinates
    //   unsigned vertex_id = qh_pointid(qh, vertex->point);
    //   vertices[i_vertex].id = vertex_id;
    //   vertices[i_vertex].point = getpoint(points, dim, vertex_id);
    //   // vertices[i_vertex].nneighsvertices = 0;
    //
    //   {
    //     unsigned* neighs = malloc(0);
    //     unsigned nneighs = 0;
    //     facetT *neighbor, **neighborp;
    //     FOREACHneighbor_(vertex) {
    //       if(neighbor->simplicial) {
    //         // vertices[i_vertex].nneighsvertices += qh_setsize(qh,
    //         // neighbor->vertices);
    //         vertexT *v, **vp;
    //         FOREACHsetelement_(vertexT, neighbor->vertices, v) {
    //           unsigned v_id = qh_pointid(qh, v->point);
    //           if(v_id > vertex_id) {
    //             unsigned pushed;
    //             appendu(v_id, &neighs, nneighs, &pushed);
    //             if(pushed) {
    //               printf("simplicial pushed\n");
    //               nneighs++;
    //             }else{
    //               printf("simplicial NOT pushed\n");
    //             }
    //           }
    //         }
    //       } else {
    //         qh_makeridges(qh, neighbor);
    //         ridgeT *r, **rp;
    //         FOREACHsetelement_(ridgeT, qh_vertexridges(qh, vertex), r) {
    //           vertexT *v, **vp;
    //           FOREACHsetelement_(vertexT, r->vertices, v) {
    //             unsigned v_id = qh_pointid(qh, v->point);
    //             if(v_id > vertex_id) {
    //               unsigned pushed;
    //               appendu(v_id, &neighs, nneighs, &pushed);
    //               if(pushed) {
    //                 printf("pushed\n");
    //                 nneighs++;
    //               }else{
    //                 printf("NOT pushed\n");
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //     vertices[i_vertex].neighvertices = neighs;
    //     vertices[i_vertex].nneighvertices = nneighs;
    //     i_vertex++;
    //     nalledgesdouble += nneighs;
    //   }
    // }
    //}

    // edges

    if(order_edges) {
      qsort(edges, nedges, sizeof(unsigned*), cmpedges);
    }
    SetOfSitesT vset = {.sites = vertices,
                        .nsites = nvertices,
                        .edges = edges,
                        .nedges = nedges};

    int curlong, totlong;
    qh_freeqhull(qh, !qh_ALL);  // free long memory
    qh_memfreeshort(qh, &curlong,
                    &totlong);  // free short memory and memory allocator

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
SEXP SiteSXP(SiteT vertex, unsigned dim, unsigned a) {
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

  unsigned nneighvertices = a ? vertex.nneighvertices : 0;
  PROTECT(neighvertices = allocVector(INTSXP, nneighvertices));
  nprotect++;
  for(unsigned i = 0; i < nneighvertices; i++) {
    INTEGER(neighvertices)[i] = 1 + vertex.neighvertices[i];
  }

  size_t l = a ? 3 : 2;
  PROTECT(R_vertex = allocVector(VECSXP, l));
  nprotect++;
  SET_VECTOR_ELT(R_vertex, 0, id);
  SET_VECTOR_ELT(R_vertex, 1, point);
  if(a) {
    SET_VECTOR_ELT(R_vertex, 2, neighvertices);
  }

  PROTECT(names = allocVector(VECSXP, l));
  nprotect++;
  SET_VECTOR_ELT(names, 0, mkChar("id"));
  SET_VECTOR_ELT(names, 1, mkChar("point"));
  if(a) {
    SET_VECTOR_ELT(names, 2, mkChar("neighvertices"));
  }
  setAttrib(R_vertex, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return R_vertex;
}

// main function ---------------------------------------------------------------
SEXP cxhullEdges_(SEXP p, SEXP a, SEXP o, SEXP errfile) {
  unsigned nprotect = 0;

  unsigned dim = ncols(p);
  unsigned n = nrows(p);

  double* points = (double*)R_alloc(n * dim, sizeof(double));
  for(unsigned i = 0; i < n; i++)
    for(unsigned j = 0; j < dim; j++)
      points[dim * i + j] = REAL(
          p)[i + n * j];  // could have been REAL(p) if p had been transposed

  unsigned orderEdges = INTEGER(o)[0];
  unsigned adjacencies = INTEGER(a)[0];

  unsigned exitcode;
  const char* e = R_CHAR(Rf_asChar(errfile));
  SetOfSitesT vset =
      cxhullEdges(points, dim, n, adjacencies, orderEdges, &exitcode, e);
  if(exitcode) {
    error("Received error code %d from qhull.", exitcode);
  }

  unsigned nvertices = vset.nsites;
  SiteT* vertices = vset.sites;
  unsigned nedges = vset.nedges;
  unsigned** edges = vset.edges;

  SEXP out, names, R_vertices, vnames, R_edges;

  PROTECT(R_vertices = allocVector(VECSXP, nvertices));
  PROTECT(vnames = allocVector(STRSXP, nvertices));
  nprotect += 2;
  for(unsigned i = 0; i < nvertices; i++) {
    SEXP vertex;
    PROTECT(vertex = SiteSXP(vertices[i], dim, adjacencies));
    nprotect++;
    SET_VECTOR_ELT(R_vertices, i, vertex);
    SET_STRING_ELT(vnames, i, Rf_asChar(VECTOR_ELT(vertex, 0)));
  }
  setAttrib(R_vertices, R_NamesSymbol, vnames);

  PROTECT(R_edges = allocMatrix(INTSXP, nedges, 2));
  nprotect++;
  for(unsigned i = 0; i < nedges; i++) {
    INTEGER(R_edges)[i] = 1 + edges[i][0];
    INTEGER(R_edges)[i + nedges] = 1 + edges[i][1];
  }

  PROTECT(out = allocVector(VECSXP, 2));
  nprotect++;
  SET_VECTOR_ELT(out, 0, R_vertices);
  SET_VECTOR_ELT(out, 1, R_edges);

  PROTECT(names = allocVector(VECSXP, 2));
  nprotect++;
  SET_VECTOR_ELT(names, 0, mkChar("vertices"));
  SET_VECTOR_ELT(names, 1, mkChar("edges"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}