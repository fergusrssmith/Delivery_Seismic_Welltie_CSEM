#include <stdio.h>
#include <stdlib.h>  /* realloc(), qsort() */
#include <math.h>
#include <float.h>
#include <string.h>

#include "vdefs.h"


Site * readone(void), * nextone(void) ;
void readsites(void) ;

int sorted=0;
int _triangulate, plot, _debug, nsites, siteidx ;
float xmin, xmax, ymin, ymax ;
Site * sites ;
Freelist sfl ;

/* 
   Compute voronoi diagram OR Delaunay triangulation
   of set of n points x[0...n-1], y[0....n-1],
   where x,y can be floats, doubles or ints.  Supply data_type =
   'f','d' or 'i' for each case. All the salient info is then found in
   the returned struct VORONOI_INFO.

   Sometimes you want to traverse each voronoi polygon. To do this call
   orderVoronoiPolygons(info) afterwards, and then the info.edgeLvertex
   has an ordered list of voronoi vertices for you

   Clean up as per the example in main(). (free_voronoi_info()).
*/
#define ELIMINATE_DUPS  0
struct VORONOI_INFO getdiagram(int debug, int triangulate, int n, char data_type, const void* x,const void* y)
{
  struct VORONOI_INFO info;
  int scomp(const void * vs1, const void * vs2);
  int i,j;
  int *xi,*yi;
  float *xf,*yf;
  double *xd,*yd;

  Site *(*next)() ;
  
  nsites = n;
  _triangulate = triangulate;
  _debug = debug;

  freeinit(&sfl, sizeof(Site)) ;

  /* initialise */
  if (data_type=='f')
	{
	  xf = (float *)x;
	  yf = (float *)y;
	  info.data_type = 'f';
	  info.x_data = xf;
	  info.y_data = yf;
	}
  else if (data_type=='d')
	{
	  xd = (double *)x;
	  yd = (double *)y;
	  info.data_type = 'd';
	  info.x_data = xd;
	  info.y_data = yd;
	}
  else if (data_type=='i')
	{
	  xi = (int *)x;
	  yi = (int *)y;
	  info.data_type = 'i';
	  info.x_data = xi;
	  info.y_data = yi;
	}
  else 
	{
	  printf("Invalid data type to getdiagram()");
	  exit(0);
	}

  sites = (Site *) myalloc(n * sizeof(Site));
  for (i=0; i< nsites; i++)
	{
	  if (data_type=='f')
		{
		  sites[i].coord.x = xf[i];
		  sites[i].coord.y = yf[i];
		}
	  else if (data_type=='d')
		{
		  sites[i].coord.x = (float) xd[i];
		  sites[i].coord.y = (float) yd[i];
		}
	  else if (data_type=='i')
		{
		  sites[i].coord.x = (float) xi[i];
		  sites[i].coord.y = (float) yi[i];
		}
	  sites[i].sitenbr = i;
	  sites[i].refcnt = 0;
	}
  qsort((void *)sites, nsites, sizeof(Site), scomp) ;

  xmin = sites[0].coord.x ;
  xmax = sites[0].coord.x ;
  int n_diff;
  n_diff = nsites;
  for (i = 0, j=0 ; i < nsites ; ++i)
	{
	  if (sites[i].coord.x < xmin)
		{
		  xmin = sites[i].coord.x ;
		}
	  if (sites[i].coord.x > xmax)
		{
		  xmax = sites[i].coord.x ;
		}

#if ELIMINATE_DUPS 
	  /* eliminate duplicates if present : not debugged yet!  */
	  if (i<(nsites-1))
		{
		  char dup=0;
		  for (j=i+1;j < nsites ;j++)
			{
			  if ((sites[j].coord.x != sites[i].coord.x)
				  || (sites[j].coord.y != sites[i].coord.y))
				{
				  break;
				}
			  else {dup=1;n_diff--; /* fprintf(stderr,"duplicates...\n"); */}
			}
		  if (dup)
			{
			  sites[i+1].coord.x = sites[j].coord.x;
			  sites[i+1].coord.y = sites[j].coord.y;
			}
		}
#endif
	}
  nsites = n_diff;

  ymin = sites[0].coord.y ;
  ymax = sites[nsites-1].coord.y ;

  next = nextone;
  siteidx = 0 ;
  geominit() ;
  initialise_info(nsites,&info);

  /* compute diagram */
  voronoi(next,&info); 



  /* clean up info */
  info.nEdges = info.edge_count;
  info.nVertices = info.vertex_count;
  if (info.siteEdgesSize != NULL)
	{
	  for (i=0;i<n;i++)
		{
		  info.siteEdgesSize[i] = info.siteEdgesCount[i];
		}
	}

  info.nTriangles = info.triangle_count;


  return info; 
}
#undef ELIMINATE_DUPS 

/*
  Find and order the sequence of vertices info.edgeLvertex defining
  the voronoi polygon, callable after getDiagram() completes
  successfully.
*/

void orderVoronoiPolygons(struct VORONOI_INFO *info)
{
  int i,j,k,vL,vR,tmp;
  struct EDGE *e,*e2;

  if (!_triangulate)
	{ 
	  /* loop over sites */
	  for (i=0;i<info->nsites;i++)
		{
		  /* loop over edges surrounding site */
		  /* sort, in-situ, edges into an ordered traversal of the polygon */
		  info->edgeLvertex[i] = (int *) malloc (info->siteEdgesSize[i] * sizeof(int));
		  for (j=0;j<info->siteEdgesSize[i]-1;j++)
			{
			  e = &(info->edges[info->siteEdges[i][j]]);
			  vL = e->i;
			  vR = e->j;
			  /* store left vertex (the critical one) */
			  info->edgeLvertex[i][j] = vL;
			  // linear search for match, then swap to end...
			  for (k=j+1;k<info->siteEdgesSize[i];k++)
				{
				  e2 = &(info->edges[info->siteEdges[i][k]]);
				  if (e2->j == vR)
					{
					  /* a match, but flipped. Flip it */
					  e2->j = e2->i;  e2->i = vR;
					}
				  if (e2->i == vR)
					{
					  /* a match */
					  tmp = info->siteEdges[i][j+1];
					  info->siteEdges[i][j+1] = info->siteEdges[i][k];
					  info->siteEdges[i][k] = tmp;
 					  break;
					}
				}
			}
		  /* left vertex of last element */
		  info->edgeLvertex[i][j] = vR;
		}
	}
}

/*
  Returns structure containing nearest gridding approximation to data x[0...n-1],y[0..n-1]
  based on examination of Delaunay triangulation contained in VORONOI_INFO *info.
  print=1 spews out stuff for ascii graphing.
  eps_cosine defines a typical epsilon for the cosine of angles that ought to be 90 degrees
  if the points are truly on a lattice.
  If the fraction of points that don't fall on the lattice according to this epsilon
  criterion is still less than acceptable_off_lattice_fraction, the code will
  return a grid based on the ones that look OK.

  Failure signalled by info field of returned NEAREST_REGULAR_GRID_INFO structure:
  info = -1 : wrong data type
         -2 : no triangles - probably collinear points, fall on 1D lattice within epsilon OK
         -3 : no triangles - probably collinear points, exceeds acceptable_off_lattice_fraction
		 -4 : too many points off-lattice wrt acceptable_off_lattice_fraction criterion
		 -10: malloc errors

		 RGinfo.Nu, RGinfo.Nv        now has number of Blocks in lattice
		 RGinfo.u_hat, RGinfo.v_hat  are unit vectors in the principal directions
		 RGinfo.u_min/max,RGinfo.v_min/max, are maximum extents in these directions
		 RGinfo.du/dv are block widths in these directions
*/

struct NEAREST_REGULAR_GRID_INFO nearest_regular_lattice(struct VORONOI_INFO *info, 
														  int print,int n, char data_type, 
														  const void* x,const void* y,
														  double eps_cosine,double acceptable_off_lattice_fraction)
{
  int i,j,j1,j2,pass;
  int cnt;
  int ret_value;
  int signDP;
  int *xi,*yi;
  float *xf,*yf;
  double *xd,*yd;
  double rA[2],rB[2],rC[2],uv1[2],uv2[2];
  double uv1_RA[2],uv2_RA[2]; 
  double duv1[2];
  double len[3],lenx[3],leny[3];
  double hyp,dotP,relCosine,s,s1,s2,s1min,s2min,s1minSq,s2minSq;
  double one_plus_eps,one_minus_eps;
  double u_min,u_max,v_min,v_max;
  char* partOfRightAngledTriangle=NULL;
  char init,aligned=0;
  struct NEAREST_REGULAR_GRID_INFO RGinfo;


  /* If collinear, won't work...*/
  if (info->nTriangles < 1) 
	{
	  nearest_regular_lattice_to_nontriangulated_points(info,&RGinfo,
														  print,n,data_type,x,y, 
														  eps_cosine, acceptable_off_lattice_fraction);
	  return RGinfo;
	} 

  if (data_type=='f')
	{
	  xf = (float *)x;
	  yf = (float *)y;
	}
  else if (data_type=='d')
	{
	  xd = (double *)x;
	  yd = (double *)y;
	}
  else if (data_type=='i')
	{
	  xi = (int *)x;
	  yi = (int *)y;
	}
  else 
	{
	 RGinfo.info=-1;
	 return RGinfo;
	} 

 

  partOfRightAngledTriangle = (char*) calloc(n,sizeof(char));
  if (!partOfRightAngledTriangle) 
	{
	  RGinfo.info=-10;
	  return RGinfo;
	}

  s1min = s2min = DBL_MAX;
  one_plus_eps = 1+eps_cosine;
  one_minus_eps = 1-eps_cosine;
  cnt = 0;
  init=1;
  // Two pass algorithm
  // pass=0: get unit vector directions & increment lengths by examining all right angled triangles
  // pass=1: mark all acceptable nodes which seem part of a lattice because they are part of a smallest right-angled
  //         triangle which is aligned with the principal axes
  for (pass=0; pass<2; pass++)
	{
	  if (pass==1)
		{
		  s1minSq = s1min*s1min;
		  s2minSq = s2min*s2min;
		}
	  for (i=0;i<info->nTriangles;i++) 
		{ 
		  // get triangle ABC
		  int * sv = info->triangleSites[i];
		  if (data_type=='f')
			{
			  rA[0] = xf[sv[0]];
			  rB[0] = xf[sv[1]];
			  rC[0] = xf[sv[2]];
			  rA[1] = yf[sv[0]];
			  rB[1] = yf[sv[1]];
			  rC[1] = yf[sv[2]];
			}
		  else if (data_type=='d')
			{
			  rA[0] = xd[sv[0]];
			  rB[0] = xd[sv[1]];
			  rC[0] = xd[sv[2]];
			  rA[1] = yd[sv[0]];
			  rB[1] = yd[sv[1]];
			  rC[1] = yd[sv[2]];
			}
		  else if (data_type=='i')
			{
			  rA[0] = xi[sv[0]];
			  rB[0] = xi[sv[1]];
			  rC[0] = xi[sv[2]];
			  rA[1] = yi[sv[0]];
			  rB[1] = yi[sv[1]];
			  rC[1] = yi[sv[2]];
			}
		  // triangle has vertices rA,rB,rC
		  // len*[0] is vector BA
		  // len*[1] is vector BC
		  // len*[2] is vector CA
		  lenx[0] = rA[0] - rB[0];
		  leny[0] = rA[1] - rB[1];
		  lenx[1] = rC[0] - rB[0];
		  leny[1] = rC[1] - rB[1];
		  lenx[2] = rA[0] - rC[0];
		  leny[2] = rA[1] - rC[1];

		  for(j=0;j<3;j++) len[j] = SQR(lenx[j])+SQR(leny[j]);
		  // identify permutation with (j+2) on longest side
		  for(j=0;j<3;j++)
			{
			  j2 = (j+2)%3;
			  j1 = (j+1)%3;
			  if ((len[j2] > len[j]) && (len[j2] > len[j1])) break;
			}
		  /*
			Side j2 longest. Triangle now looks like
		             j2
		             | \ len[j2]
		     len[j1] |  \
		             |   \
		             j1---j
		              len[j]
		   Test if corner j1 is 90 degrees
		  */
		  hyp = len[j2];
		  dotP = lenx[j]*lenx[j1] + leny[j]*leny[j1];
		  relCosine = fabs(dotP/hyp);
		  if (relCosine < eps_cosine)
			{					
			  // A right angled triangle
			  // Tests for validly are:
			  // (1) All rightsized right-angled triangles confer 'satisfactory' status on corner points
			  // (2) A rightsized triangles must align with main directions, and have shortest
			  //     sides of matching length with principal directions increments
			  if (pass == 0)
				{
				  // Creat unit grid vectors
				  if (init)
					{
					  // First right-angled triangle. Create unit increment vectors inline, xline
					  make_unit_vectors(1,lenx[j],leny[j],uv1,uv2,uv1_RA);					  
					  init = 0;
					  cnt = 1;
					}
				  // Test current triangle increments against established directions uv1, uv2
				  s1 = sqrt(len[j]);
				  s2 = sqrt(len[j1]);
				  dotP = lenx[j]*uv1[0] + leny[j]*uv1[1];  //  side_j . uv1 
				  signDP = (dotP > 0) ? 1 : -1;
				  dotP = fabs(dotP);  // | side_j . uv1 |
				  aligned =0;
				  if (dotP > s1*(one_minus_eps))
					{
					  // side_j aligned with uv1
					  aligned = 1;
					  duv1[0] = signDP*lenx[j]; duv1[1] = signDP*leny[j];
					}
				  else 
					{
					  dotP = lenx[j1]*uv1[0] + leny[j1]*uv1[1];  //  side_j1 . uv1 
					  signDP = (dotP > 0) ? 1 : -1;
					  dotP = fabs(dotP);                         // | side_j1 . uv1 |
					  if (dotP > s2*(one_minus_eps))
						{
						  // side_j aligned with uv2  <=> side_j1 aligned with uv1
						  aligned = 2;
						  // swap meanings of s1,s2
						  s = s1; s1=s2; s2=s;
						  duv1[0] = signDP*lenx[j1]; duv1[1] = signDP*leny[j1];
						}
					}
				  if (aligned != 0)
					{
					  if ((s1 < s1min) || (s2 < s2min))  
						{
						  if (s1 < s1min)
							{
							  s1min = s1;
							  // re-initialise unit vectors
							  make_unit_vectors(1,duv1[0],duv1[1],uv1,uv2,uv1_RA);	
							  cnt=1;
							}
						  if (s2 < s2min) 
							{
							  s2min = s2;
							  // reform unit vectors
							  make_unit_vectors(1,duv1[0],duv1[1],uv1,uv2,uv1_RA);	
							  cnt=1;
							}
						}
					  else if ((s1 > (s1min*one_minus_eps)) && (s1 < (s1min*one_plus_eps))
							   && (s2 > (s2min*one_minus_eps)) && (s2 < (s2min*one_plus_eps)))
						{
						  // Triangle similar to best smallest right angled triangle found to date
						  // Refine average
						  cnt++;
						  s1min = ((cnt-1)*s1min+s1)/cnt;
						  s2min = ((cnt-1)*s2min+s2)/cnt;
						  // re-average unit vectors
						  make_unit_vectors(cnt,duv1[0],duv1[1],uv1,uv2,uv1_RA);	
						}
					}
				}
			  else
				{
				  // On pass 2, increment vectors are uv1 and uv2, 
				  // of lengths s1min & s2min respectively
				  // Mark this node as having a right angle between shortest two sides
				  partOfRightAngledTriangle[sv[j1]] = 1;
						  
				  dotP = lenx[j]*uv1[0] + leny[j]*uv1[1];  // side_j . uv1
				  if (fabs(dotP) < s1min * eps_cosine)
					{
					  // side_j aligned with uv2
					  s = fabs(len[j])/s2minSq;
					  if ((s < one_plus_eps) && (s > one_minus_eps))
						{
						  // side_j has basic incremental length in uv1 direction
						  // mark node sv[j]
						  partOfRightAngledTriangle[sv[j]] = 1;
						}
					  s = fabs(len[j1])/s1minSq;
					  if ((s < one_plus_eps) && (s > one_minus_eps))
						{
						  // side_(j+1) has basic incremental length in uv2 direction
						  // mark node sv[j]
						  partOfRightAngledTriangle[sv[j2]] = 1;
						}
					}
				  else if (fabs(lenx[j]*uv2[0] + leny[j]*uv2[1]) < s2min * eps_cosine)
					{
					  // side_j aligned with uv1
					  s = fabs(len[j])/s1minSq;
					  if ((s < one_plus_eps) && (s > one_minus_eps))
						{
						  // side_j has basic incremental length in uv1 direction
						  // mark node sv[j]
						  partOfRightAngledTriangle[sv[j]] = 1;
						}
					  s = fabs(len[j1])/s2minSq;
					  if ((s < one_plus_eps) && (s > one_minus_eps))
						{
						  // side_(j+1) has basic incremental length in uv2 direction
						  // mark node sv[j]
						  partOfRightAngledTriangle[sv[j2]] = 1;
						}
					}
				  // else {return -3;}
				} // pass chunks 
			} // detect right-angled triangle
		} // loop on triangles
	} // pass loop

  // Now, all nodes should be marked as valid. If too many aren't 
  // we have to give up
  for (i=0,j=0;i<n;i++)
	{
	  if (!partOfRightAngledTriangle[i]) j++;
	}
  if ((((double)j)/n > acceptable_off_lattice_fraction))
	{
	  free(partOfRightAngledTriangle);
	  RGinfo.info = -4;
	  return RGinfo;	
	}
  

  // if we got to here, everything sits on lattice OK, within the "epsilons" stipulated
  // we have unit-vectors, so determine lattice span
  RGinfo.u_min = RGinfo.v_min = DBL_MAX;
  RGinfo.u_max = RGinfo.v_max = -DBL_MAX;
  for (i=0;i<n;i++)
 	{
	  switch (data_type)
		{
		case 'f':
		  {
			rA[0] = ((float *)x)[i];
			rA[1] = ((float *)y)[i];
			break;
		  }
		case 'd':
		  {
			rA[0] = ((double *)x)[i];
			rA[1] = ((double *)y)[i];
			break;
		  }
		case 'i':
		  {
			rA[0] = ((int *)x)[i];
			rA[1] = ((int *)y)[i]; 
			break;
		  }
		}
	  s1 = rA[0]*uv1[0]+rA[1]*uv1[1];
	  RGinfo.u_min = MIN(RGinfo.u_min,s1);
	  RGinfo.u_max = MAX(RGinfo.u_max,s1);
	  s2 = rA[0]*uv2[0]+rA[1]*uv2[1];
	  RGinfo.v_min = MIN(RGinfo.v_min,s2);
	  RGinfo.v_max = MAX(RGinfo.v_max,s2);
	  if (print) printf("#m=0,S=%d\n%g %g\n", (partOfRightAngledTriangle[i])? 4: 6, rA[0],rA[1]);
	}

  // Set grid...
  // rounding
  RGinfo.Nu = (int)((RGinfo.u_max-RGinfo.u_min)/s1min + 0.5);
  RGinfo.Nv = (int)((RGinfo.v_max-RGinfo.v_min)/s2min + 0.5);
  // unit vectors
  memcpy(&RGinfo.uhat[0],&uv1[0],2*sizeof(double));  
  memcpy(&RGinfo.vhat[0],&uv2[0],2*sizeof(double)); 
  // increments
  RGinfo.du = (RGinfo.Nu > 0) ? ((RGinfo.u_max-RGinfo.u_min)/(RGinfo.Nu)) : 1;
  RGinfo.dv = (RGinfo.Nv > 0) ? ((RGinfo.v_max-RGinfo.v_min)/(RGinfo.Nv)) : 1;

  if (print)
   {
	 fprintf(stderr,"du=%g, dv=%g\n",s1min,s2min);
	 fprintf(stderr,"u_min=%g, u_max=%g\n",RGinfo.u_min,RGinfo.u_max);
	 fprintf(stderr,"v_min=%g, v_max=%g\n",RGinfo.v_min,RGinfo.v_max);
	 fprintf(stderr,"Nu=%d, Nv=%d\n",RGinfo.Nu,RGinfo.Nv);

	 printf("#m=3,S=0\n");
	 for (s1=RGinfo.u_min,i=0; i<=RGinfo.Nu; i++, s1 += RGinfo.du)
	   {
		 printf("%g %g\n",s1*uv1[0]+ RGinfo.v_min*uv2[0], s1*uv1[1] + RGinfo.v_min*uv2[1] );
		 printf("%g %g\n\n",s1*uv1[0] + RGinfo.v_max*uv2[0],s1*uv1[1] + RGinfo.v_max*uv2[1] );
	   }
	 for (s2=RGinfo.v_min,i=0; i<=RGinfo.Nv; i++, s2 += RGinfo.dv)
	   {
		 printf("%g %g\n",RGinfo.u_min*uv1[0] + s2*uv2[0], RGinfo.u_min*uv1[1] + s2*uv2[1] );
		 printf("%g %g\n\n",RGinfo.u_max*uv1[0] + s2*uv2[0], RGinfo.u_max*uv1[1] + s2*uv2[1] );
	   }
   }

  // clean up
  free(partOfRightAngledTriangle);

  RGinfo.info=0; 
   

  return RGinfo;
}

/* 
   Construct best representation of set of points that cannot be triangulated. Probaly a
   straight line.

   RGinfo.Nu, RGinfo.Nv        now has number of Blocks in lattice
   RGinfo.u_hat, RGinfo.v_hat  are unit vectors in the principal directions
   RGinfo.u_min/max,RGinfo.v_min/max, are maximum extents in these directions
   RGinfo.du/dv are block widths in these directions

   For these degenerate case, Nv=0. u_hat is the direction of the line
*/
void nearest_regular_lattice_to_nontriangulated_points(struct VORONOI_INFO *info,
														 struct NEAREST_REGULAR_GRID_INFO *RGinfo,
														 int print,int n, char data_type, 
														 const void* x,const void* y,
														 double eps_cosine,double acceptable_off_lattice_fraction)
{
  int i,count;
  double **r;
  double *dists;
  double s,spacing;

  /* Simplest thing. copy to doubles, sort by x & y, then compute *modal* increment */
  // Signal 1D lattice */
  RGinfo->info = -2;

  /* allocate pointers to rows */
  r = (double**)malloc (n* sizeof(double*));
  if (!r ) 
	{
	  RGinfo->info=-10;
	  return;
	}
  /* allocate rows and set pointers to them */
  r[0] = (double*) malloc((size_t)(2 * n * sizeof(double)));
  if (!r[0]) 
	{
	  free(r); 
	  RGinfo->info=-10;
	  return;
	}
  for (i=1;i<n;i++)  r[i] = r[i-1] + 2;
  for (i=0;i<n;i++)
	{
	  switch (data_type)
		{
		case 'f':
		  r[i][0] = ((float*)x)[i];
		  r[i][1] = ((float*)y)[i];
		  break;
		case 'd':
		  r[i][0] = ((double*)x)[i];
		  r[i][1] = ((double*)y)[i];
		  break;
		case 'i':
		  r[i][0] = ((int*)x)[i];
		  r[i][1] = ((int*)y)[i];
		  break;
		}
	}
  // sort by x,then y
  qsort(&r[0][0],n,2*sizeof(double),compare_pairs_by_x_and_y);
  // list of successive square spacings
  dists = (double*) malloc((size_t)((n-1) * sizeof(double)));
  if (!dists)
	{
	  RGinfo->info=-10;
	  return;
	}
  for (i=0;i<(n-1);i++) dists[i] = SQR(r[i+1][0] - r[i][0]) + SQR(r[i+1][1] - r[i][1]);
  // get median. this would be faster as quickselect, but this'll do for now:
  qsort(dists,n-1,sizeof(double),compare_doubles);
  spacing = sqrt(dists[n/2]);


  // Align uhat with (almost certainly linear) span
  RGinfo->uhat[0] = r[n-1][0] - r[0][0];
  RGinfo->uhat[1] = r[n-1][1] - r[0][1];
  // Total span
  s = hypot( RGinfo->uhat[0], RGinfo->uhat[1]);
  if (s>0) 
	{
	  RGinfo->uhat[0]/=s; 
	  RGinfo->uhat[1]/=s;
	}
  else
	{
	  // no finite span. Probably degenerate pts. Stuff in x,y unit vectors FWIW
	  RGinfo->uhat[0]=1; RGinfo->uhat[1]=0;
	}
  // rotate 90deg for vhat
  RGinfo->vhat[0] = -RGinfo->uhat[1];
  RGinfo->vhat[1] = RGinfo->uhat[0];
  // number of nodes
  RGinfo->Nu = (spacing>0) ?  ( (int)(s/spacing + 0.5)) : 0;  // round to nearest
  RGinfo->Nv = 0; // assume this must be the case
  // principal direction spacings
  RGinfo->du = (RGinfo->Nu > 0) ? s/(RGinfo->Nu) : 0.0;
  RGinfo->dv = 0.0;
  // principal direction bounds
  RGinfo->u_min = r[0][0] * RGinfo->uhat[0] + r[0][1] * RGinfo->uhat[1];
  RGinfo->u_max = r[n-1][0] * RGinfo->uhat[0] + r[n-1][1] * RGinfo->uhat[1];
  RGinfo->v_min = r[0][0] * RGinfo->vhat[0] + r[0][1] * RGinfo->vhat[1];
  RGinfo->v_max = r[n-1][0] * RGinfo->vhat[0] + r[n-1][1] * RGinfo->vhat[1];

  // Off-lattice check
  if (spacing>0)
	{
	  for (i=1,count=0;i<n;i++) 
		{
		  s = sqrt(SQR(r[i][0] - r[0][0]) + SQR(r[i][1] - r[0][1]));
		  if (fmodf(s,spacing) > eps_cosine) count++;
		}
	  if (count > (acceptable_off_lattice_fraction * n))
		{
		  RGinfo->info = -3;
		}
	}

 if (print)
   {
	 fprintf(stderr,"du=%g, dv=%g\n",RGinfo->du,RGinfo->dv);
	 fprintf(stderr,"u_min=%g, u_max=%g\n",RGinfo->u_min,RGinfo->u_max);
	 fprintf(stderr,"v_min=%g, v_max=%g\n",RGinfo->v_min,RGinfo->v_max);
	 fprintf(stderr,"Nu=%d, Nv=%d\n",RGinfo->Nu,RGinfo->Nv);

	 printf("#m=3,S=8\n");
	 double s1;
	 for (s1=RGinfo->u_min,i=0; i<=RGinfo->Nu; i++, s1 += RGinfo->du)
	   {
		 printf("%g %g\n",s1*RGinfo->uhat[0]+ RGinfo->v_min *RGinfo->vhat[0],
				s1*RGinfo->uhat[1]+ RGinfo->v_min *RGinfo->vhat[1]);
	   }
   }

  // clean up
  free(dists); dists = NULL;
  free(r[0]); r[0] = NULL;
  free(r); r = NULL;
}

int compare_doubles(const void *p1, const void *p2)
{
  return (*(double*)p1 > *(double*)p2) ? 1 : ((*(double*)p1 < *(double*)p2) ? -1 : 0);
}

int compare_pairs_by_x_and_y(const void *p1, const void *p2)
{
  double *ra,*rb;
  ra = (double*)p1;
  rb = (double*)p2;
  
  if (ra[0] > rb[0])
	{
	  return 1;
	}
  else if (ra[0] < rb[0])
	{
	  return -1;
	}
  else
	{
	  return (ra[1] > rb[1]) ? 1 : ((ra[1] < rb[1]) ? -1 : 0);	  
	}
}

void make_unit_vectors(int count,double dx,double dy,double uv1[2],double uv2[2],double uv1_RA[2])
{
  double s;

  if (count==1)
	{
	  // initialise
	  uv1[0] = uv1_RA[0] = dx; 
	  uv1[1] = uv1_RA[1] = dy;
	}
  else
	{
	  // accumulate average directions
	  uv1_RA[0] = ((count-1)*uv1_RA[0]+dx)/count;
	  uv1_RA[1] = ((count-1)*uv1_RA[1]+dy)/count;
	  uv1[0] = uv1_RA[0];
	  uv1[1] = uv1_RA[1];
	}
  s = hypot(uv1[0],uv1[1]); 
  if (s > 0)
	{
	  uv1[0]/=s; 
	  uv1[1]/=s;
	}
  uv2[0] = -uv1[1]; 
  uv2[1] = uv1[0];   // rotated 90 deg clockwise
}


/*
 Call only if voronoi has been requested & computed
*/
float computeMeanSiteSeparation(struct VORONOI_INFO *info)
{
  int i,indx1,indx2;
  float dx,dy,dsq;
  struct EDGE *e;

  /* loop over edges */
  dsq=0.0;
  for (i=0;i<info->nEdges;i++)
	{
	  e = &(info->edges[i]);
	  indx1 = e->indx1;  
	  indx2 = e->indx2;
	  switch (info->data_type)
		{
		case 'f':
		  dx = ((float*)(info->x_data))[indx1] -  ((float*)(info->x_data))[indx2];
		  dy = ((float*)(info->y_data))[indx1] -  ((float*)(info->y_data))[indx2];
		  break;

		case 'd':
		  dx = ((double*)(info->x_data))[indx1] -  ((double*)(info->x_data))[indx2];
		  dy = ((double*)(info->y_data))[indx1] -  ((double*)(info->y_data))[indx2];
		  break;

		case 'i':
		  dx = ((int*)(info->x_data))[indx1] -  ((int*)(info->x_data))[indx2];
		  dy = ((int*)(info->y_data))[indx1] -  ((int*)(info->y_data))[indx2];
		  break;
		}
	  // printf("#m=2,S=0\n%g %g\n", ((float*)(info->x_data))[indx1] , ((float*)(info->y_data))[indx1] );
	  // printf("%g %g\n\n", ((float*)(info->x_data))[indx2] , ((float*)(info->y_data))[indx2] );
	  dsq += dx*dx + dy*dy;
	}
  dsq = sqrt(dsq/(info->nEdges));
  // exit(0);
  return dsq;
}

/* 
   Useful test function to see if point is in polygon. INCOMPLETE....FIXME
*/
int isInVoronoiPolygon(double x, double y,struct VORONOI_INFO *info)
{
  /*
  if (i<0 || i>=info->nsites) return -1;
  pts = info->edgeLvertex[i];
  for (j=0; j<info->siteEdgesSize[i]; j++) 
	{
	  if Area2(x,y, 
			   info->vertices[info->edgeLvertex[j]].x, info->vertices[info->edgeLvertex[j]].y,
			   info->vertices[info->edgeLvertex[j+1]].x, info->vertices[info->edgeLvertex[j+1]].y) >=0
			   
	}
  */
  return 0;
}




/*** sort sites on y, then x, coord ***/
int scomp(const void * vs1, const void * vs2)
{
  Point * s1 = (Point *)vs1 ;
  Point * s2 = (Point *)vs2 ;

  if (s1->y < s2->y)
	{
	  return (-1) ;
	}
  if (s1->y > s2->y)
	{
	  return (1) ;
	}
  if (s1->x < s2->x)
	{
	  return (-1) ;
	}
  if (s1->x > s2->x)
	{
	  return (1) ;
	}
  return (0) ;
}

/*** return a single in-storage site ***/

Site *nextone(void)
{
  Site * s ;

  if (siteidx < nsites)
	{
	  s = &sites[siteidx++];
	  return (s) ;
	}
  else
	{
	  return ((Site *)NULL) ;
	}
}

/* Voronoi array additions */
void initialise_info(int nsites,struct VORONOI_INFO *info)
{
  int i;
  int DEFAULT_N_EDGES_GUESS = 8;
  info->nsites = nsites;
  info->nEdges = 0;
  info->nVertices = 0;
  if (!_triangulate)
	{
	  /* edges */
	  /* The wierd numbers are the result of testing large grids on random points in a square */
	  /* The ratios are close to 3 and 2 rspectively, so this is a safety margin */
	  info->nEdges = (int)(3.2*nsites);
	  info->edges = (struct EDGE *) malloc ((info->nEdges) * sizeof(struct EDGE));
	  info->edge_count=0;

	  /* vertices */
	  info->nVertices= (int)(2.2*nsites);
	  info->vertices = (struct VERTEX *) malloc ((info->nVertices) * sizeof(struct VERTEX));
	  info->vertex_count=0;
 
	  /* siteEdges */
	  info->siteEdges = (int **) malloc (nsites * sizeof(int *));
	  info->siteEdgesSize = (int *) malloc (nsites * sizeof(int));
	  info->siteEdgesCount = (int *) calloc (nsites, sizeof(int));  // initialise to 0
	  for (i=0;i<nsites;i++)  
		{
		  info->siteEdgesSize[i] = DEFAULT_N_EDGES_GUESS;
		  info->siteEdges[i] = (int *) malloc(DEFAULT_N_EDGES_GUESS * sizeof(int));
		}
	  /* edgeLvertex */
	  info->edgeLvertex = (int **) malloc (nsites * sizeof(int *));


	  /* Collect no Delaunay */
	  info->nTriangles=0;
	  info->triangleSites = NULL;
	  info->triangle_count=0;
	}
  else 
	{
	  info->edges=NULL;
	  info->vertices=NULL;
	  info->siteEdges =NULL;
	  info->siteEdgesSize =NULL;
	  info->siteEdgesCount = NULL;

	  // JG: I think an upper bound on the # triangles is 2n ’¡Ý 2 ’¡Ý b, where
	  // the convex hull has b vertices, which must be b>1 I think.
	  // So, for safety,
	  info->nTriangles= (int)(2*nsites-2);

	  /* allocate pointers to rows */
	  info->triangleSites = (int**) malloc ((info->nTriangles) * sizeof(int*));
	  /* allocate rows and set pointers to them */
	  info->triangleSites[0] = (int*) malloc((size_t)(info->nTriangles * 3 * sizeof(int)));
	  for (i=1;i<info->nTriangles;i++)  
		{
		  info->triangleSites[i] = info->triangleSites[i-1] + 3;
		}

	  info->triangle_count=0;
	}
}

/*** sort edges on first index, then second ***/
int edgecomp(const void * vs1, const void * vs2)
{
  struct EDGE *s1 = (struct EDGE *)vs1 ;
  struct EDGE *s2 = (struct EDGE *)vs2 ;

  if (s1->i < s2->i)
	{
	  return (-1) ;
	}
  if (s1->i > s2->i)
	{
	  return (1) ;
	}
  if (s1->j < s2->j)
	{
	  return (-1) ;
	}
  if (s1->j > s2->j)
	{
	  return (1) ;
	}
  return (0) ;
}


void free_voronoi_info(struct VORONOI_INFO *info)
{
  int i;
  if (info->edges) free(info->edges);
  if (info->vertices) free(info->vertices);
  /* siteEdges */
  if (info->siteEdgesCount) free (info->siteEdgesCount);
  if (info->siteEdgesSize) free(info->siteEdgesSize);
  if (info->siteEdges)
	{
	  for (i=0;i<info->nsites;i++) 
		{
		  if (info->siteEdges[i]) free(info->siteEdges[i]);
		}
	  free(info->siteEdges);
	}
  if (!_triangulate)
	{
	  if (info->edgeLvertex != NULL) 
		{
		  for (i=0;i<info->nsites;i++) 
			{
			  if (info->edgeLvertex[i]) 
				{
				  free(info->edgeLvertex[i]);
				  info->edgeLvertex[i] = NULL;
				}
			}
		  free(info->edgeLvertex);
		  info->edgeLvertex = NULL;
		}
	}

  if (info->triangleSites)
	{
	  free(info->triangleSites[0]);
	  info->triangleSites[0] = NULL;
	  free(info->triangleSites);
	  info->triangleSites = NULL;
	}

} 


void add_edge(Edge * e,struct VORONOI_INFO *info)
{
  int s1,s2;
  if (!_triangulate)
	{
	  if (info->edge_count == info->nEdges) reAllocEdges(info);
	  struct EDGE * ep = &info->edges[info->edge_count];
	  // printf("e %d", e->edgenbr);
	  ep->i = (e->ep[le] != (Site *)NULL ? e->ep[le]->sitenbr : -1);
	  ep->j = (e->ep[re] != (Site *)NULL ? e->ep[re]->sitenbr : -1);
	  ep->indx1 = s1 = e->reg[le]->sitenbr;
	  ep->indx2 = s2 = e->reg[re]->sitenbr;

	  /* bookkeeping for left site */
	  if (info->siteEdgesCount[s1] >= info->siteEdgesSize[s1])
		{ 
		  info->siteEdgesSize[s1] *= 2;
		  info->siteEdges[s1] = (int *) realloc(info->siteEdges[s1], info->siteEdgesSize[s1] * sizeof(int));
		}

	  info->siteEdges[s1][info->siteEdgesCount[s1]] = info->edge_count;
	  info->siteEdgesCount[s1]++;

	  /* bookkeeping for right site */
	  if (info->siteEdgesCount[s2] >= info->siteEdgesSize[s2])
		{
		  info->siteEdgesSize[s2] *= 2;
		  info->siteEdges[s2] = (int *) realloc(info->siteEdges[s2], info->siteEdgesSize[s2] * sizeof(int));
		}

	  info->siteEdges[s2][info->siteEdgesCount[s2]] = info->edge_count;
	  info->siteEdgesCount[s2]++;


	  info->edge_count++;
	}		  
}

/* Dynamic array control */
void reAllocEdges(struct VORONOI_INFO *info)
{
  info->nEdges *=2;
  info->edges = (struct EDGE *) realloc (info->edges, (info->nEdges) * sizeof(struct EDGE));
}

void add_vertex(Site * v,struct VORONOI_INFO *info)
{
  if (!_triangulate)
	{
	  if (info->vertex_count == info->nVertices) reAllocVertices(info);
	  info->vertices[info->vertex_count].indx = v->sitenbr;
	  info->vertices[info->vertex_count].x = v->coord.x;
	  info->vertices[info->vertex_count].y = v->coord.y;
	  info->vertex_count++;
	}		  
}

/* Dynamic array control */
void reAllocVertices(struct VORONOI_INFO *info)
{
  info->nVertices *=2;
  info->vertices = (struct VERTEX *) realloc (info->vertices, (info->nVertices) * sizeof(struct VERTEX));
}

