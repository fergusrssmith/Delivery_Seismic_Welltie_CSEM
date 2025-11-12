#ifndef __VDEFS_H
#define __VDEFS_H

# ifndef MIN
#  define       MIN(a,b)                        (((a) < (b)) ? (a) : (b))
# endif

# ifndef MAX
#  define       MAX(a,b)                        (((a) > (b)) ? (a) : (b))
# endif

# ifndef SQR
#  define       SQR(x)                          ((x) * (x))
# endif


#ifndef NULL
#define NULL 0
#endif

#define DELETED -2

struct EDGE
{
  int i;      /*  vertex index at end 1 */
  int j;      /*  vertex index at end 2 */
  int indx1;  /*  indices of sites which this edge bisects */
  int indx2;
};

struct VERTEX
{
  int indx;
  double x,y;
};



typedef struct tagFreenode
{
  struct tagFreenode * nextfree;
} Freenode ;


typedef struct tagFreelist
{
  Freenode * head;
  int nodesize;
} Freelist ;

typedef struct tagPoint
{
  float x ;
  float y ;
} Point ;

/* structure used both for sites and for vertices */

typedef struct tagSite
{
  Point coord ;
  int sitenbr ;
  int refcnt ;
} Site ;

struct VORONOI_INFO
{
  int nsites;
  int nEdges;
  struct EDGE* edges;       /* array of indices of joined vertices */
  int nVertices;
  struct VERTEX* vertices;   /* array of vertex points */
  int **siteEdges;           /* siteEdges[i][0...k] is list of edges for site i */
  int **edgeLvertex;         /* edgeLVertex[i][0...k] is list of vertex indices for site i, ordered to traverse the polygon */
  int* siteEdgesSize;
  char data_type;           /* f,d, or i */ 
  void *x_data,*y_data;     /* ptrs to original data */
  /* Delaunay triangles if requested */
  int nTriangles;
  int** triangleSites;

  /* internal accounting */
  /* counters */
  int edge_count;
  int vertex_count;
  int* siteEdgesCount;
  int triangle_count;
};

struct NEAREST_REGULAR_GRID_INFO
{
  char info;
  double uhat[2],vhat[2];
  int Nu,Nv;
  double du,dv;
  double u_min,u_max,v_min,v_max;
};


typedef struct tagEdge
{
  float a, b, c ;
  Site * ep[2] ;
  Site * reg[2] ;
  int edgenbr ;
} Edge ;

#define le 0
#define re 1

typedef struct tagHalfedge
{
  struct tagHalfedge * ELleft ;
  struct tagHalfedge * ELright ;
  Edge * ELedge ;
  int ELrefcnt ;
  char ELpm ;
  Site * vertex ;
  float ystar ;
  struct tagHalfedge * PQnext ;
} Halfedge ;

/* edgelist.c */
void ELinitialize(void) ;
Halfedge * HEcreate(Edge *, int) ;
void ELinsert(Halfedge *, Halfedge *) ;
Halfedge * ELgethash(int) ;
Halfedge * ELleftbnd(Point *) ;
void ELdelete(Halfedge *) ;
Halfedge * ELright(Halfedge *) ;
Halfedge * ELleft(Halfedge *) ;
Site * leftreg(Halfedge *) ;
Site * rightreg(Halfedge *) ;
extern int ELhashsize ;
extern Site * bottomsite ;
extern Freelist hfl ;
extern Halfedge * ELleftend, * ELrightend, **ELhash ;

/* geometry.c */
void geominit(void) ;
Edge * bisect(Site *, Site *) ;
Site * intersect(Halfedge *, Halfedge *) ;
int right_of(Halfedge *, Point *) ;
void endpoint(Edge *, int, Site *,struct VORONOI_INFO *info) ;
float dist(Site *, Site *) ;
void makevertex(Site *,struct VORONOI_INFO *);
void deref(Site *) ;
void ref(Site *) ;
extern float deltax, deltay ;
extern int nsites, nedges, sqrt_nsites, nvertices ;
extern Freelist sfl, efl ;

/* heap.c */
void PQinsert(Halfedge *, Site *, float) ;
void PQdelete(Halfedge *) ;
int PQbucket(Halfedge *) ;
int PQempty(void) ;
Point PQ_min(void) ;
Halfedge * PQextractmin(void) ;
void PQinitialize(void) ;
extern int PQmin, PQcount, PQhashsize ;
extern Halfedge * PQhash ;

/* main.c */
extern int sorted, triangulate, plot, _debug, nsites, siteidx ;
extern float xmin, xmax, ymin, ymax ;
extern Site * sites ;
extern Freelist sfl ;

/* getopt.c */
extern int getopt(int, char *const *, const char *);

/* memory.c */
void freeinit(Freelist *, int) ;
char *getfree(Freelist *) ;
void makefree(Freenode *, Freelist *) ;
char *myalloc(unsigned) ;

/* output.c */
void out_bisector(Edge *) ;
void out_ep(Edge *) ;
void out_vertex(Site *);
void out_site(Site *) ;
void out_triple(Site *, Site *, Site *,struct VORONOI_INFO *info);

/* voronoi.c */
void voronoi(Site *(*)(),struct VORONOI_INFO*) ;

/* getdiagram.c */
struct VORONOI_INFO getdiagram(int debug, int doTriangulation, int n, char data_type, const void* x,const void* y);

void orderVoronoiPolygons(struct VORONOI_INFO *info);
struct NEAREST_REGULAR_GRID_INFO nearest_regular_lattice(struct VORONOI_INFO *info,
														  int print,int n, char data_type, const void* x,const void* y,
														  double eps_cosine,double acceptable_off_lattice_fraction);
void nearest_regular_lattice_to_nontriangulated_points(struct VORONOI_INFO *info,
													   struct NEAREST_REGULAR_GRID_INFO *RGinfo,
													   int print,int n, char data_type, 
													   const void* x,const void* y,
													   double eps_cosine,double acceptable_off_lattice_fraction);
void make_unit_vectors(int count,double lenx_j,double leny_j,double uv1[2],double uv2[2],double uv1_RA[2]);
float computeMeanSiteSeparation(struct VORONOI_INFO *info);
void free_voronoi_info(struct VORONOI_INFO *info);
void initialise_info(int nsites,struct VORONOI_INFO *info);
void add_edge(Edge * e,struct VORONOI_INFO *info);
void reAllocEdges(struct VORONOI_INFO *info); 
void add_vertex(Site * v,struct VORONOI_INFO *info);
void reAllocVertices(struct VORONOI_INFO *info);
int edgecomp(const void * vs1, const void * vs2);
int compare_pairs_by_x_and_y(const void *p1, const void *p2);
int compare_doubles(const void *p1, const void *p2);

#endif  


