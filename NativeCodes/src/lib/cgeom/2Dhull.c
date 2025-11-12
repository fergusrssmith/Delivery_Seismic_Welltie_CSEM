/*
  Thread-safe adaptation of Graham 2D convex hull algorithm for library routine.
  J Gunning, CSIRO Petroleum Mar 05. 
  
  Usage:
    given a set of data points data[0...n-1][0,1]
	where 0,1, is x,y on the last index, call

    twoDhull(int n, double** data,char dataType, int print,int *size)

	to get index list of convex hull points.

	example:

	int size;
	int* indx = twoDhull(n,data,'f',0,&size)

	returns an array indx[0....(size-1)] of the indices
	in data[i][0,1] that define the convex hull. Set dataType='f', 'd','i'
	for floats, doubles, ints.

	To see verbose output, set print=1, else set print=0.
	Remember to free *indx when no longer required.

	To make a test binary, set MAKE_TEST_BINARY to 1
	use
	
	% gcc 2Dhull.c -lpthread -lm -o 2Dhull

	To restore/compile library, MAKE_TEST_BINARY 0
	and make.


	Originally from O'Rourke, "Computational Geometry in C" (Second Edition),
	Chapter 3.  
	
	--------------------------------------------------------------------
	This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
	redistributed in its entirety provided that this copyright notice is
	not removed.
	--------------------------------------------------------------------
*/

#include   <stdio.h>
#include   <math.h>
#include   <stdlib.h> 
//#include   <values.h> 
#include   <float.h> 
#include   "2Dhull.h"


#ifndef  MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef  MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif


#define MAKE_TEST_BINARY 0

#define X   0
#define Y   1
typedef enum { FALSE, TRUE }   bool;

#define DIM 2                 /* Dimension of points */
typedef float tPointf[DIM];   /* Type float point */

#define NEW(p, type) if ((p=(type *) malloc (sizeof(type))) == NULL) {\
      fprintf (stdout,"NEW: Out of Memory!\n");\
      exit(EXIT_FAILURE);\
    }

#define FREE(p)  if (p) { free ((char *) p); p = NULL; }

#define SWAP(t,x,y)     { t = x; x = y; y = t; }


/*----------Point(s) Structure-------------*/
typedef struct tPointStructure tsPoint;
typedef tsPoint *tPoint;
struct tPointStructure {
  int     vnum;
  tPointf v;
  bool    delete;
  //  const void *gbs; /* JG addition will be ptr to GRAHAM_GLOBALS */
};

/*----------Stack Structure-------------*/
typedef struct tStackCell tsStack;
typedef tsStack *tStack;
// typedef tsPoint tPointArray[PMAX];

/* Old Global variables in cgeom original */
struct tStackCell 
{
  tPoint   p;
  tStack   next;
};

typedef struct GRAHAM_GLOBALS
{
  int     vnum;
  tPointf v;
  bool    delete;
  int n;                         /* Actual # of points */
  int ndelete;                   /* Number deleted */

  // tPointArray P; // static ??
  tsPoint *P;
  struct tStackCell tsStack;
} GRAHAM_GLOBALS;

/* thread-safe global pointer needed by qsort */

// (1) No thread local storage on OSX. Library 
// will not be thread-safe until OSX supports it
// or we rewrite... sigh.
#ifdef __APPLE__
 #define THREAD
#else
 #define THREAD __thread 
#endif

THREAD struct GRAHAM_GLOBALS *_globals=NULL;

/*----------Internal Function Prototypes-------------*/

tStack	Pop( tStack s );
tStack  Push( tPoint p, tStack top );
tStack  Graham(GRAHAM_GLOBALS*);
void    Squash(GRAHAM_GLOBALS*);
void	Copy(GRAHAM_GLOBALS*, int i, int j );
int     CompareInfo( const void *tp1, const void *tp2,void* globals);
int     Compare( const void *tp1, const void *tp2);
void    FindLowest(GRAHAM_GLOBALS*);
void	Swap(GRAHAM_GLOBALS*, int i, int j );
int     AreaSign( tPointf a, tPointf b, tPointf c );
double  Area2( tPointf a, tPointf b, tPointf c );
bool    Left( tPointf a, tPointf b, tPointf c );
void    FreeStack(tStack top);
int*    getBoundaryIndexArray(GRAHAM_GLOBALS* gbs, tStack t,int* size);


#if MAKE_TEST_BINARY

int main(int argc, char *argv[])
{
  void data_test(void *a);

  int i,j,tmp,NUM_THREADS=1;
  pthread_t thread[NUM_THREADS];
  
  if (argc>1) srand(atoi(argv[1])); // seed rand generator
  
  for (i=0; i<NUM_THREADS; i++ )
	{
	  // spawn threads to make multiple data sets and print convex hulls
	  tmp = pthread_create( &thread[i],NULL, &data_test, (void *)&i);
	}
  for (i=0; i<NUM_THREADS; i++ ) pthread_join (thread[i], NULL);

  return EXIT_SUCCESS;
}

void data_test(void *a) 
{
  int print,getMinBoundingRectangleAlso;
  int i,j,n,tmp;
  double** data;
  double* corners=NULL;
  int offset;
  char dataType='d';
  int size,*list;
  
  // flags
  print = 1;
  getMinBoundingRectangleAlso = 1;

  // data size
  n = 117;

  offset = *((int*)a);
  // printf("Offset=%d\n",offset);

  data = (double **)malloc(n*sizeof(double *));
  data[0] = (double *)malloc((n*2)*sizeof(double));
  for (i=0;i<n;i++)  if (i>0) data[i] =  data[i-1] + 2;

  if (0)
	{
	  // random data
	  for (i=0;i<n;i++)
		{
		  data[i][0] = 2*offset + rand()/(double)RAND_MAX;
		  data[i][1] = data[i][0] + rand()/(double)RAND_MAX;
		}
	}
  else
	{
	  // noisy rotated incomplete regular grid
	  // random data
	  int k,m;
	  float noise=0.2, theta=.4,c,s,x,y;
	  m = (int)sqrt(n);
	  c = cos(theta);
	  s = sin(theta);
	  for (i=0,k=0;k<n;i++)
		{
		  for (j=0;(j<m) && (k<n); j++,k++)
			{
			  x = 1*(i + noise*(rand()/(double)RAND_MAX - 0.5));
			  y =  j + noise*(rand()/(double)RAND_MAX - 0.5);
			  data[k][0] = c*x - s*y;
			  data[k][1] = s*x + c*y;
			}
		}
	}

  list = twoDhull(n,(void**)data,'d',print,&size,getMinBoundingRectangleAlso,&corners);

  free(corners);
  free(list);
  free(data[0]);
  free(data);
  /* terminate the thread */
  pthread_exit(NULL);
}
#endif

/* 
   Called as option to twoDhull(), if bounding rectangle requested.
   See main() example

   Given the convex hull from a call to hullIndices = twoDhull(int n,
   void **data, char dataType, int print,int *size) this returns a
   list of corner corrdinates of a minimum area bounding rectangle,
   corners=[x0,y0,x1,y1,x2,y2,x3,y3], with 4 points in clockwise
   order.  Remember to free double array when finished.

   Algorithm does simple loop over all smallest bounding rectangles
   with one side collinear with an edge of the convex hull, saving the
   smallest. Since convex hull has usual few edges, this is efficient.
*/
double* getMinimumAreaBoundingRectangle( int n, void **data, char dataType,int print, int* indxArr, int indxArrSize)
{
  int i,j,iA,iB;
  double u[2],v[2],uhat[2],vhat[2];
  double norm,dU,dV,Umin,Umax,Vmin,Vmax;
  double min_area,area;
  double* corners;
  int** data_i=NULL;
  float** data_f=NULL;
  double** data_d=NULL;

  switch (dataType) 
	{
	case 'i':
	  data_i = (int**)data; break;
	case 'f':
	  data_f = (float**)data; break;
	case 'd':
	  data_d = (double**)data; break;
	}


  corners = (double*)malloc(8*sizeof(double));
  min_area = DBL_MAX;
  for (i=0;i<indxArrSize;i++)
	{
	  iA = i;
	  iB = (i+1) % indxArrSize; // points to 0 when i = indxArrSize-1
	  switch (dataType) 
		{
		case 'i':
		  u[0] = data_i[indxArr[iA]][0];
		  u[1] = data_i[indxArr[iA]][1];
		  uhat[0] = data_i[indxArr[iB]][0];
		  uhat[1] = data_i[indxArr[iB]][1];
		  break;
		case 'f':
		  u[0] = data_f[indxArr[iA]][0];
		  u[1] = data_f[indxArr[iA]][1];
		  uhat[0] = data_f[indxArr[iB]][0];
		  uhat[1] = data_f[indxArr[iB]][1];
		  break;
		case 'd':
		  u[0] = data_d[indxArr[iA]][0];
		  u[1] = data_d[indxArr[iA]][1];
		  uhat[0] = data_d[indxArr[iB]][0];
		  uhat[1] = data_d[indxArr[iB]][1];
		  break;
		}
	  // fix origin at point u. Get unit vectors uhat, vhat wrt hull segment joining iA, iB
	  // where u is the vector to vertex iA, uhat a unit vector from point iA to point iB along hull.
	  uhat[0] -= u[0]; uhat[1] -= u[1];
	  norm = hypot(uhat[0],uhat[1]);
	  if (norm>0)
		{
		  uhat[0]/=norm; uhat[1]/=norm;
		}
	  vhat[1] = uhat[0];  vhat[0] = -uhat[1];
	  
	  // loop over vertices v of the hull, finding spanning rectangle 
	  // with orientation defined by {uhat,vhat}
	  Umin=Umax=Vmin=Vmax=0.0;
	  for (j=0; j< indxArrSize; j++)
		{
		  switch (dataType) 
			{
			case 'i':
			  v[0] = data_i[indxArr[j]][0];
			  v[1] = data_i[indxArr[j]][1];
			  break;
			case 'f':
			  v[0] = data_f[indxArr[j]][0];
			  v[1] = data_f[indxArr[j]][1];
			  break;
			case 'd':
			  v[0] = data_d[indxArr[j]][0];
			  v[1] = data_d[indxArr[j]][1];
			  break;
			} 
		  // vector from u to v
		  v[0] -= u[0]; v[1] -= u[1];
		  // coords in current {uhat,vhat} system, origin u
		  dU = v[0]*uhat[0] + v[1]*uhat[1];
		  dV = v[0]*vhat[0] + v[1]*vhat[1];
		  Umin = MIN(dU,Umin);
		  Umax = MAX(dU,Umax);
		  Vmin = MIN(dV,Vmin);
		  Vmax = MAX(dV,Vmax);
		}
	  area = (Umax-Umin)*(Vmax-Vmin);
	  if (area < min_area)
		{
		  // save rectangle
		  min_area = area;
		  for (j=0;j<2;j++) // x,y
			{
			  corners[j  ] = u[j] + Umin*uhat[j] + Vmin*vhat[j];
			  corners[2+j] = u[j] + Umax*uhat[j] + Vmin*vhat[j];
			  corners[4+j] = u[j] + Umax*uhat[j] + Vmax*vhat[j];
			  corners[6+j] = u[j] + Umin*uhat[j] + Vmax*vhat[j];
			}
		}
	}

  if (print)
	{
	  printf("#m=2,S=0\n");
	  for (i=0;i<10;i+=2)
		{
		 printf("%lf %lf\n",corners[i%8],corners[(i+1)%8]);
	   }
	 printf("\n");
   }
  return corners;
}

#if 0
/*
  Given minimum bounding rectangle in corners,
  on assumption that supplied data sequence {x,y} is a raster scan in some orientation,
  compute the number of lattice points in the u and v directions, where these
  directions come from the minimum bounding rectangle
  // FIXME NOT FUNCTIONAL YET
*/
void getNearestLatticeSpacing(double* corners,int n, void **data, char dataType, int *Nu, int *Nv)
{
  int i,j;
  int n_inline_incs;
  char inLineIsU,inLineIsV ;
  double len,lenU,lenV;
  double v[2],uhat[2],vhat[2];
  int** data_i=NULL;
  float** data_f=NULL;
  double** data_d=NULL;

  switch (dataType) 
	{
	case 'i':
	  data_i = (int**)data; break;
	case 'f':
	  data_f = (float**)data; break;
	case 'd':
	  data_d = (double**)data; break;
	}

  *Nu = *Nv = -1;

  // get unit vectors
  uhat[0] = corners[2] - corners[0];
  uhat[1] = corners[3] - corners[1];
  vhat[0] = corners[7] - corners[0];
  vhat[1] = corners[8] - corners[1];

  lenU = hypot(uhat[0],uhat[1]);
  if (lenU > 0) {uhat[0]/=lenU; uhat[1]/=lenU;}

  lenV = hypot(vhat[0],vhat[1]);
  if (lenV > 0) {vhat[0]/=lenV; vhat[1]/=lenV;}


  n_inline_incs=0;
  inLineIsU = inLineIsV = 0;
  inline_av_dist = 0.0;

  for (i=0;i<n;i++)
	{
	  switch (dataType) 
		{
		case 'i':
		  v[0] = data_i[i];
		  v[1] = data_i[i];
		  break;
		case 'f':
		  v[0] = data_f[i];
		  v[1] = data_f[i];
		  break;
		case 'd':
		  v[0] = data_d[i];
		  v[1] = data_d[i];
		  break;
		} 
	  if (i>1)
		{
		  dotPu = (v[0] - u[0])*uhat[0] + (v[1] - u[1])*uhat[1];
		  dotPv = (v[0] - u[0])*vhat[0] + (v[1] - u[1])*vhat[1];
		  if (i==1)
			{
			  if (abs(dotPv) < eps * abs(dotPu))
				{
				  inLineIsU = 1;
				}
			  else if (abs(dotPu) < eps * abs(dotPv))
				{
				  inLineIsV = 1;
				}
			  else
				{
				  // Grid not aligned with MBR
				  return -1;
				}
			}
		  if ((inLineIsU && (abs(dotPv) < eps * abs(dotPu)))
			  || (inLineIsV && (abs(dotPu) < eps * abs(dotPv))))
			{
			  len = inLineIsU ? dotPu ; dotPv;
			  // uv vector in uhat direction
			  n_inline_incs++;
			  // update average
			  inline_av_dist = (inline_av_dist*(n_inline_incs-1) + len)/n_inline_incs;
			}
		}
	  u[0] = v[0];
	  u[1] = v[1];
	}
  if (inLineIsU)
	{
	  if (fmod(lenU,abs(inline_av_dist)) > eps)
		{
		  // MBR doesn't divide well into inline increments
		  return -2;
		}
	  else
		{
		  *Nu = (int)(lenU/abs(inline_av_dist) + 0.5);
		}
	}
  else
	{
	  if (fmod(lenV, abs(inline_av_dist)) > eps)
		{
		  // MBR doesn't divide well into inline increments
		  return -2;
		}
	  else
		{
		  *Nv = (int)(lenV/abs(inline_av_dist) + 0.5);
		}
	}
  // FIXME UP TO HERE. DO SOMETHING USEFUL NOW FOR OTHER Nu,Nv....

}
#endif

/**
    twoDhull(int n, double** data,char dataType, int print,int *size,
	        int getMinBoundingRectangleAlso,double **corners)

	Returns index list of convex hull points.

	example:

	int size;
	int* indx = twoDhull(n,data,'f',0,&size,0,NULL)

	returns an array indx[0....(size-1)] of the indices
	in data[i][0,1] that define the convex hull. Set dataType='f', 'd','i'
	for floats, doubles, ints.

	To see verbose output, set print=1, else set print=0.
	Remember to free *indx when no longer required.

	If getMinBoundingRectangleAlso ==1, alloc and returns also
	corners = [x0,y0,x1,y1,x2,y2,x3,y3]
	bounding corners of smallest area bounding rectangle
	with 4 points in clockwise order.
	Remember to free double array when finished. See main() example
*/
int* twoDhull(int n, void **data, char dataType, int print,int *size,
			  int getMinBoundingRectangleAlso,double **corners)
{ 
  int i;
  tStack   top;
  GRAHAM_GLOBALS gbs;

  if (dataType!='i' && dataType!='f' && dataType!='d')
	{
	  fprintf(stderr,"Wrong data type to twoDhull()");
	  *size=0;
	  return NULL;
	}

  _globals = &gbs;
  
   // allocate space for pts
   gbs.P = (tsPoint*)malloc(n*sizeof(tsPoint));
   gbs.n = n;
   gbs.ndelete = 0;

   for (i=0;i<n;i++)
	 {
	   switch (dataType) {
	   case 'i':
		 gbs.P[i].v[0] = ((int**)data)[i][0];
		 gbs.P[i].v[1] = ((int**)data)[i][1];
		 break;
	   case 'f':
		 gbs.P[i].v[0] = ((float**)data)[i][0];
		 gbs.P[i].v[1] = ((float**)data)[i][1];
		 break;
	   case 'd':
		 gbs.P[i].v[0] = ((double**)data)[i][0];
		 gbs.P[i].v[1] = ((double**)data)[i][1];
		 break;
	   }
	   gbs.P[i].vnum = i;
	   gbs.P[i].delete = FALSE;
	   //	   gbs.P[i].gbs = &gbs;      /* needed for sort */
	 }
   
   
   FindLowest(&gbs);

   qsort(
		 &gbs.P[1],         /* pointer to 1st elem */
		 n-1,               /* number of elems */
		 sizeof( tsPoint ), /* size of each elem */
		 Compare           /* -1,0,+1 compare function */			 
		 );

   if (gbs.ndelete > 0) 
	 {
	   Squash(&gbs);
	 }
   
   top = Graham(&gbs);

   int *indxArr;
   indxArr = getBoundaryIndexArray(&gbs,top,size);

   
   if (getMinBoundingRectangleAlso)
	 {
	   *corners = getMinimumAreaBoundingRectangle(n,data,dataType,print,indxArr,*size);
	 }
   else 
	 {
	   *corners=NULL;
	 }
   if (print)
	 {
	   // Output
	   // Points
	   printf("#m=0,S=17\n");
	   switch (dataType) {
	   case 'i':
		 for (i = 0; i < gbs.n; i++) printf("%d %d\n",((int**)data)[i][0],((int**)data)[i][1]);
		 break;
	   case 'f':
		 for (i = 0; i < gbs.n; i++) printf("%f %f\n",((float**)data)[i][0],((float**)data)[i][1]);
		 break;
	   case 'd':
		 for (i = 0; i < gbs.n; i++) printf("%lf %lf\n",((double**)data)[i][0],((double**)data)[i][1]);
		 break;
	   }
	   printf("\n");

	   /* Bounding polygon. */
	   printf("#m=1,S=0\n");
	   switch (dataType) {
	   case 'i':
		 for(i=0;i<*size;i++) printf("%d %d\n",((int**)data)[indxArr[i]][0],((int**)data)[indxArr[i]][1]);
		 // final loop
		 printf("%d %d\n",((int**)data)[indxArr[0]][0],((int**)data)[indxArr[0]][1]);
		 break;
	   case 'f':
		 for(i=0;i<*size;i++) printf("%f %f\n",((float**)data)[indxArr[i]][0],((float**)data)[indxArr[i]][1]);
		 // final loop
		 printf("%f %f\n",((float**)data)[indxArr[0]][0],((float**)data)[indxArr[0]][1]);
		 break;
	   case 'd':
		 for(i=0;i<*size;i++) printf("%lf %lf\n",((double**)data)[indxArr[i]][0],((double**)data)[indxArr[i]][1]);
		 // final loop
		 printf("%lf %lf\n",((double**)data)[indxArr[0]][0],((double**)data)[indxArr[0]][1]);
		 break;
	   }
	 
	 }
   /* free stuff */
	  
   FreeStack(top);
   free(gbs.P);

   _globals =NULL;
   return indxArr;
}

/*---------------------------------------------------------------------
FindLowest finds the rightmost lowest point and swaps with 0-th.
The lowest point has the min y-coord, and amongst those, the
max x-coord: so it is rightmost among the lowest.
---------------------------------------------------------------------*/
void   FindLowest(GRAHAM_GLOBALS *gbs)
{
   int i;
   int m = 0;   /* Index of lowest so far. */

   for ( i = 1; i < gbs->n; i++ )
	 {
	   if ( (gbs->P[i].v[Y] <  gbs->P[m].v[Y]) ||
			((gbs->P[i].v[Y] == gbs->P[m].v[Y]) && (gbs->P[i].v[X] > gbs->P[m].v[X])) ) 
         {
		   m = i;
		 }
	 }
   Swap(gbs,0,m); /* Swap gbs->P[0] and gbs->P[m] */
}

void Swap(GRAHAM_GLOBALS *gbs, int i, int j )
{
   int temp;
   float ftemp;
   /* Uses swap macro. */

   SWAP( temp, gbs->P[i].vnum,   gbs->P[j].vnum );
   SWAP( ftemp, gbs->P[i].v[X],   gbs->P[j].v[X] );
   SWAP( ftemp, gbs->P[i].v[Y],   gbs->P[j].v[Y] );
   SWAP( temp, gbs->P[i].delete, gbs->P[j].delete );

}

/*---------------------------------------------------------------------
Compare: returns -1,0,+1 if p1 < p2, =, or > respectively;
here "<" means smaller angle.  Follows the conventions of qsort.
---------------------------------------------------------------------*/
int   Compare( const void *tpi, const void *tpj)
{
   int a;             /* area */
   int x, y;          /* projections of ri & rj in 1st quadrant */
   tPoint pi, pj;
   pi = (tPoint)tpi;
   pj = (tPoint)tpj;
   GRAHAM_GLOBALS *gbs = (GRAHAM_GLOBALS*)(_globals);

   a = AreaSign(gbs->P[0].v, pi->v, pj->v );
   if (a > 0)
      return -1;
   else if (a < 0)
      return 1;
   else { /* Collinear with P[0] */
      x =  fabs( pi->v[X] -  gbs->P[0].v[X] ) - fabs( pj->v[X] -  gbs->P[0].v[X] );
      y =  fabs( pi->v[Y] -  gbs->P[0].v[Y] ) - fabs( pj->v[Y] -  gbs->P[0].v[Y] );

      gbs->ndelete++;
      if ( (x < 0) || (y < 0) ) {
         pi->delete = TRUE;
         return -1;
      }
      else if ( (x > 0) || (y > 0) ) {
         pj->delete = TRUE;
         return 1;
      }
      else { /* points are coincident */
         if (pi->vnum > pj->vnum)
             pj->delete = TRUE;
         else 
             pi->delete = TRUE;
         return 0;
      }
   }
}



/*---------------------------------------------------------------------
Pops off top elment of stack s, frees up the cell, and returns new top.
---------------------------------------------------------------------*/
tStack   Pop( tStack s )
{
   tStack top;

   top = s->next;
   FREE( s );
   return top;
}

/*---------------------------------------------------------------------
Get a new cell, fill it with p, and push it onto the stack.
Return pointer to new stack top.
---------------------------------------------------------------------*/
tStack   Push( tPoint p, tStack top )
{
   tStack   s;

   /* Get new cell and fill it with point. */
   NEW( s, tsStack );
   s->p = p;
   s->next = top;
   return s;
}

/*--------------------------------------------------------------------
Frees the stack afterwards
---------------------------------------------------------------------*/
void FreeStack(tStack top)
{
   while (top) {
	  top = Pop( top );
   }
}

/*---------------------------------------------------------------------
Performs the Graham scan on an array of angularly sorted points P.
---------------------------------------------------------------------*/
tStack   Graham(GRAHAM_GLOBALS *gbs)
{
   tStack   top;
   int i;
   tPoint p1, p2;  /* Top two points on stack. */

   /* Initialize stack. */
   top = NULL;
   top = Push ( &(gbs->P[0]), top );
   top = Push ( &(gbs->P[1]), top );

   /* Bottom two elements will never be removed. */
   i = 2;

   while ( i < gbs->n ) {
      if( !top->next) {
		fprintf(stdout,"Error\n");
		return NULL;
	  }
      p1 = top->next->p;
      p2 = top->p;
      if ( Left( p1->v , p2->v, (gbs->P)[i].v ) ) {
         top = Push ( &(gbs->P[i]), top );
         i++;
      } else {
		top = Pop( top );
	  }
   }

   return top;

}

/*---------------------------------------------------------------------
Squash removes all elements from P marked delete.
---------------------------------------------------------------------*/
void Squash(GRAHAM_GLOBALS *gbs)
{
   int i, j;
   i = 0; j = 0;

   while ( i < gbs->n ) {
      if ( !(gbs->P[i].delete) ) { /* if not marked for deletion */
         Copy(gbs, i, j ); /* Copy P[i] to P[j]. */
         j++;
      }
      /* else do nothing: delete by skipping. */
      i++;
   }
   gbs->n = j;
}

void Copy( GRAHAM_GLOBALS *gbs, int i, int j )
{
   (gbs->P)[j].v[X] = (gbs->P)[i].v[X];
   (gbs->P)[j].v[Y] = (gbs->P)[i].v[Y];
   (gbs->P)[j].vnum = (gbs->P)[i].vnum;
   (gbs->P)[j].delete = (gbs->P)[i].delete;
}


/*---------------------------------------------------------------------
Returns twice the signed area of the triangle determined by a,b,c.
The area is positive if a,b,c are oriented ccw, negative if cw,
and zero if the points are collinear.
---------------------------------------------------------------------*/
double Area2( tPointf a, tPointf b, tPointf c )
{
   return ((b[X] - a[X]) * (c[Y] - a[Y]) - (c[X] - a[X]) * (b[Y] - a[Y]));
}

/*---------------------------------------------------------------------
Returns true iff c is strictly to the left of the directed
line through a to b.
---------------------------------------------------------------------*/
bool   Left( tPointf a, tPointf b, tPointf c )
{ 
   return  (Area2( a, b, c ) > 0);
}

int* getBoundaryIndexArray(GRAHAM_GLOBALS* gbs, tStack t,int* size)
{
  int i;
  int* indexArray;
  tStack tn;
  
  tn=t;
  *size=0;
  while (tn) {
	(*size)++;
	tn = tn->next;
  }
  tn=t;

  indexArray = malloc((*size)*sizeof(int));
  i=0;
  while (tn) {
	indexArray[i++]= tn->p->vnum;
	tn = tn->next;
  }
  return indexArray;
}

int AreaSign( tPointf a, tPointf b, tPointf c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * ( c[1] - a[1] ) -
            ( c[0] - a[0] ) * ( b[1] - a[1] );

	if      ( area2 > 0) return  1;
    else if ( area2 < 0 ) return -1;
    else                     return  0;
}


#undef X
#undef Y
#undef DIM
#undef NEW
#undef FREE
#undef SWAP

