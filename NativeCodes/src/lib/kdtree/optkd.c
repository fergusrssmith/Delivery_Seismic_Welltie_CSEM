/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
/*
  Simple kdtree library kludged together from fragments of "Ranger",
  found on the Stonybrook Algorithm repository:
  http://www.cs.sunysb.edu/~algorith/implement/ranger/implement.shtml
  Modified to allow more general 2D and 3D Euclidean
  metrics commonly used in geostatistics.

  James Gunning, CSIRO, Mar 2004.

*/

 
/**************************************************************************/
/* Program to perform orthogonal range searches and nearest neighbor      */
/* querys in a more sophisticated k-d tree.  In this implementation the,  */
/* nodes on any given level of the tree do not have the same              */
/* discriminating dimension as the discrimiator is chosen based on the    */
/* dimension with   the "maxspead."                                       */
/*                                                                        */
/* References:  J.H. Friedman, J.L. Bentley, R.A. Finkel  "An Algorithm   */
/* for Finding Best Matches in Logarithmic Expected Time."                */
/* ACM Transactions on Mathematical Software, Vol 3 No. 3 Sept. 1977      */
/* pp. 209-226.                                                           */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "nn.h"
#include "metric.h"
#include "pqueue.h"
#include "optkd.h"
#include "NR/nrutil.h"
#include "thread.h"

// all these globals below subsumed in optKglobals member

//static int *perm;  /* permutation array */
//int *optfound;
//static int Count=0;
//double *nndist;

//double (*Distance)();

static THREAD int recurs_clean_level=0;

/***************************************************************************/

void Selection(a, l,N, k,discrim, p)

	 /***************************************************************************/
	 /* Makes the perm partition the array Values along the element k.          */
	 /* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
	 /***************************************************************************/

	 double **a;
	 int l,N,k,discrim;
	 optkdNode* p;
	 
{
  double v;
  int t,i,j,r;
  int* perm;

  perm = (*p->Kglobals).perm;

  r=N;

  while(r>l) {
	v=a[perm[r]][discrim]; i=l-1; j=r;
	for (;;) {
	  while (a[perm[++i]][discrim] < v);
	  while (a[perm[--j]][discrim] > v && j>l); 
	  if (i >= j) break;
	  t=perm[i]; perm[i] = perm[j]; perm[j]=t;
	}
	t=perm[i]; perm[i] = perm[r]; perm[r]=t;
	if (i>=k) r=i-1;
	if (i<=k) l=i+1;
  }
}

/****************************************************************************/

int findmaxspread(l,u,dimension,points,p)

/****************************************************************************/

	 int l,u,dimension;
	 double **points;
	 optkdNode* p;

{
  int i,j,maxdim;
  double max       =-DBL_MAX;
  double min       = DBL_MAX;
  double maxspread =-DBL_MAX;

  int *perm;

  perm = (*p->Kglobals).perm;

  for (i=0; i < dimension; i++) {
	max =-DBL_MAX;
	min = DBL_MAX;
	for (j=l; j <= u; j++) {
	  if (max < points[perm[j]][i]) { 
		max = points[perm[j]][i];
	  }
	  if (min > points[perm[j]][i]) { 
		min = points[perm[j]][i];
	  }
	  if (maxspread < fabs(max-min)) {
		maxspread = fabs(max-min);
		maxdim = i;
	  }
	}
  }
  return(maxdim);
}

/*******************************************************************************/

optkdNode *BuildkdTree(points,l,u,dimension,globals)

	 /*******************************************************************************/

	 int l,u;
	 double **points;
	 int dimension;
	 optKglobals *globals;
{
  optkdNode *p;
  int *perm;
  int m;


  if ((p = (optkdNode *)malloc(sizeof(optkdNode)))==NULL) 
	{
	  fprintf(stderr,"memory error\n");
	  return (optkdNode *)NULL;
	}

  p->Kglobals = globals; // JG. Ref back to single copy of globals
  perm =  globals->perm; // JG

  if (u-l+1 <= BUCKETSIZE) {
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
    p->loson = NULL;
    p->hison = NULL;
  } else {
    p->bucket =0;
    p->discrim = findmaxspread(l,u,dimension,points,p);
    m=(l+u)/2;
    Selection(points,l,u,m,p->discrim,p);
    p->cutval = points[perm[m]][p->discrim];
    p->loson = BuildkdTree(points,l,m,dimension,globals);
    p->hison = BuildkdTree(points,m+1,u,dimension,globals);
  }


  return(p);
}

/*******************************************************************************/

optkdNode *BuildOptTree(points,numPoints,dimension)

	 /*******************************************************************************/
	 int dimension,numPoints;
	 double **points;
{

  int j;
  optKglobals *globals=NULL;

  // allocate single instance of globals
  if ((globals = (optKglobals *)malloc(sizeof(optKglobals)))==NULL) 
	{
	  fprintf(stderr,"memory error in optKglobals\n");
	  return (optkdNode *)NULL;
	}
  /* with sensible initialisations */
  globals->nndist = NULL;
  globals->distance_metric = NULL;
  globals->wksp1 = NULL;
  globals->wksp2 = NULL;

  /* initialize perm array */

  (*globals).perm = (int *) malloc(numPoints*sizeof(int));
  if (!(*globals).perm)
	{
	  free(globals);
	  fprintf(stderr,"memory error in optKglobals\n");
	  return (optkdNode *)NULL;
	}

  for (j=0; j < numPoints; j++) (*globals).perm[j]=j;

  return (BuildkdTree(points,0,numPoints-1,dimension,globals));
}


/*******************************************************************************/

void rnnEuclidean(p,querpoint,points,dimension,numpoints,Kglobals)
	 
	 /*******************************************************************************/

	 /* special searching algorithm to take advantage of the fact that square roots
		do not need to be evaluated */

	 optkdNode *p;
	 double *querpoint;
	 double **points;
	 int dimension,numpoints;
	 optKglobals* Kglobals;
{
  int i,j,k;
  int* perm;
  int* optfound;
  double d,thisdist,val,thisx;
  
  perm = Kglobals->perm;
  optfound = Kglobals->optfound;

  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist=0.0;
      for (j=0; j<dimension; j++) {
        d=(querpoint[j]-points[perm[i]][j]);
        thisdist += d*d;
      }        

      if (optfound[0] < numpoints) {
        PQInsert_kd(thisdist,perm[i],Kglobals->nndist,optfound);
      } else {
        PQreplace_kd(thisdist,Kglobals->nndist,optfound,perm[i]);
      }
    }
  } else {
    val = querpoint[p->discrim] - p->cutval;
    if (val < 0) {
      rnnEuclidean(p->loson,querpoint,points,dimension,numpoints,Kglobals);
      if (Kglobals->nndist[1] >= val*val) {
        rnnEuclidean(p->hison,querpoint,points,dimension,numpoints,Kglobals);
      }
    } else {
      rnnEuclidean(p->hison,querpoint,points,dimension,numpoints,Kglobals);
      if (Kglobals->nndist[1] >= val*val) {
        rnnEuclidean(p->loson,querpoint,points,dimension,numpoints,Kglobals);
      }
    }
  }
}

/*
  Method to set up a metric and relevant workspace for a rnnEuclideanGen()
  query, with a metric defined, in 2D, by an azimuth angle and principle
  direction length ratio r = lambda_x/lambda_y.
  Ensure r!=0.
  The azimuth angle IN RADIANS! is measured clockwise.

*/
void setUp2DEuclidenGenMetric(optKglobals* Kglobals,double angle_radians, double r)
{
  int i;
  double st,ct,r1,r2=1.0;

  Kglobals->dim = 2;
  Kglobals->distance_metric = dmatrix(0,Kglobals->dim-1,0,Kglobals->dim-1);
  r1=1.0/fabs(r*r);
  if (r1<1.0)
	{
	  r2=1.0/fabs(r1);
	  r1=1.0;
	} 
  ct = cos(angle_radians);
  st = sin(angle_radians);
  Kglobals->distance_metric[0][0] = r1*ct*ct+r2*st*st;
  Kglobals->distance_metric[0][1] = (r2-r1)*ct*st;
  Kglobals->distance_metric[1][0] = Kglobals->distance_metric[0][1];
  Kglobals->distance_metric[1][1] = r1*st*st+r2*ct*ct;

  /* useful workspace...*/
  Kglobals->wksp1 = dvector(0,Kglobals->dim-1);
  Kglobals->wksp2 = NULL;
}

/*
  Method to set up a metric and relevant workspace for a rnnEuclideanGen()
  query, with a metric defined, in 3D, using angles and lengths 
  as per GSLIB definitions. See section II.3 and section IV.3 in "GSLIB: geostatistical
  software library and users guide, second edition" for documentation on what the
  three angles and ranges mean.

  ang1             Azimuth angle for principal direction
  ang2             Dip angle for principal direction
  ang3             Third rotation angle
  a_hmax           Major horizontal corr length
  a_hmin           Minor horizontal corr length
  a_vert           Vertical corr length

  N.B. Angles are in degrees.
  Usually expects a_hmax > a_hmin I believe.

  Kludged from GSLIB.

  Author: James Gunning CSIRO.
*/
void setUp3DEuclidenGenMetric(optKglobals* Kglobals,double a_hmax,double a_hmin,double a_vert, 
							  double ang1,double ang2,double ang3)
{
  int i,j,k;
  double min,max,sum,aa1,aa2,aa3;
  double afac1,afac2,sina,sinb,sint,cosa,cosb,cost;
  static double DEG2RAD = 0.01745329251994;   /* pi/180 */
  static double EPSLON = 1.e-20;
  double anis1,anis2,alpha,beta,theta;
  double rotmat[3][3];

  /* Allocate space */
  Kglobals->dim = 3;
  Kglobals->distance_metric = dmatrix(0,Kglobals->dim-1,0,Kglobals->dim-1);

  /* For kdtree: scale all lengths so the minimum length is unity */
  min = MIN(MIN(fabs(a_hmax),fabs(a_hmin)) , fabs(a_vert));
  /* equivalent scaled lengths */
  aa1 = fabs(a_hmax)/min;
  aa2 = fabs(a_hmin)/min;
  aa3 = fabs(a_vert)/min;
  max = MAX(MAX(aa1,aa2),aa3);
  /* compute equivalents to GSLIB's anis1 anis2 parameters */
  anis1 = aa2/max;
  anis2 = aa3/max;
  /* now continue on with usual GSLIB stuff.... */
  //   ang1             Azimuth angle for principal direction
  //   ang2             Dip angle for principal direction
  //   ang3             Third rotation angle
  //   anis1            First anisotropy ratio
  //   anis2            Second anisotropy ratio
  //
  // Converts the input angles to three angles which make more
  //  mathematical sense:
  //
  //         alpha   angle between the major axis of anisotropy and the
  //                 E-W axis. Note: Counter clockwise is positive.
  //         beta    angle between major axis and the horizontal plane.
  //                 (The dip of the ellipsoid measured positive down)
  //         theta   Angle of rotation of minor axis about the major axis
  //                 of the ellipsoid.
  //
  if (ang1>=0.0 && ang1<270.0)
	{
	  alpha = (90.0   - ang1) * DEG2RAD;
	}
  else
	{
	  alpha = (450.0  - ang1) * DEG2RAD;
	}
  beta  = -1.0 * ang2 * DEG2RAD;
  theta =        ang3 * DEG2RAD;
//
// Get the required sines and cosines:
//
  sina  = sin(alpha);
  sinb  = sin(beta);
  sint  = sin(theta);
  cosa  = cos(alpha);
  cosb  = cos(beta);
  cost  = cos(theta);
  //
  // Construct the rotation matrix in the required memory:
  //
  afac1 = 1.0 / MAX(anis1,EPSLON);
  afac2 = 1.0 / MAX(anis2,EPSLON);
  rotmat[0][0] =       (cosb * cosa);
  rotmat[0][1] =       (cosb * sina);
  rotmat[0][2] =       (-sinb);
  rotmat[1][0] = afac1*(-cost*sina + sint*sinb*cosa);
  rotmat[1][1] = afac1*(cost*cosa + sint*sinb*sina);
  rotmat[1][2] = afac1*( sint * cosb);
  rotmat[2][0] = afac2*(sint*sina + cost*sinb*cosa);
  rotmat[2][1] = afac2*(-sint*cosa + cost*sinb*sina);
  rotmat[2][2] = afac2*(cost * cosb);

  // Then, the quadratic distance is ||rotmat.dr||^2 = dr^T.rotmat^T.rotmat.dr/aa1^2
  // as per the GSLIB function sqdist. So let's stuff rotmat^T.rotmat into
  // Kglobals->distance_metric, and normalise so minimum diagonal entry =1

  /* hack to make minimum metric entry on diagonal =1 */
  min = DBL_MAX;
  for (i=0;i<Kglobals->dim;i++)
	{
	  for (j=0;j<=i;j++)
		{
		  sum=0;
		  for (k=0;k<Kglobals->dim;k++)
			{
			  sum += rotmat[k][i] * rotmat[k][j];
			}
		  Kglobals->distance_metric[i][j] = sum/(aa1*aa1);
		  Kglobals->distance_metric[j][i] = Kglobals->distance_metric[i][j];		  
		}
	  min=MIN(min,Kglobals->distance_metric[i][i]);
	}
  /* hack to make minimum metric entry on diagonal =1 */
  for (i=0;i<Kglobals->dim;i++)
	{
	  for (j=0;j<Kglobals->dim;j++)
		{
		  Kglobals->distance_metric[i][j] /= min;
		}
	}

  /* useful workspace...*/
  Kglobals->wksp1 = dvector(0,Kglobals->dim-1);
}


















/*
  Generalisation to rotated ellipse case.
  Make sure you set up the metrics and workspace with a call to
  e.g. setUp2DEuclidenGenMetric() (or higher dim equivalents)
  before calling this
 
  Author: James Gunning, CSIRO, 2004.
*/
/*******************************************************************************/

void rnnEuclideanGen(p,querpoint,points,dimension,numpoints,Kglobals)

/*******************************************************************************/
	 
optkdNode *p;
double *querpoint;
double **points;
int dimension,numpoints;
optKglobals* Kglobals;
{
  int i,j,k,index;
  int* perm;
  int* optfound;
  double d,thisdist,val,cutval,thisx;
  double tmp,*dist_vec;  
  perm = Kglobals->perm;
  optfound = Kglobals->optfound;
  dist_vec = Kglobals->wksp1;

  if (p->bucket) 
	{
	  for (i=p->lopt; i <= p->hipt; i++) 
		{
		  // want:
		  // thisdist = Kglobals->Distance(points,perm[i],querpoint,dimension,MinkP);
		  // e.g. ManhattDist(Points,Index,NNQPoint,dimension)
		  // returns the Manhattan distance between Points[Index] and NNQPoint


		  thisdist=0.0;
		  index=perm[i];
		  for (k=0; k<dimension; k++) dist_vec[k] = (querpoint[k]-points[index][k]);
		  for (j=0; j<dimension; j++) 
			{
			  tmp=0.0; 
			  // distance from quadratic form
			  for (k=0; k<dimension; k++)  
				{
				  tmp += Kglobals->distance_metric[j][k]*dist_vec[k];
				}
			  thisdist += dist_vec[j]*tmp;
			}
		  if (optfound[0] < numpoints) 
			{
			  PQInsert_kd(thisdist,perm[i],Kglobals->nndist,optfound);
			} 
		  else 
			{
			  PQreplace_kd(thisdist,Kglobals->nndist,optfound,perm[i]);
			}
		}
	} 
  else 
	{
	  cutval = p->cutval;
	  thisx = querpoint[p->discrim];
	  val = querpoint[p->discrim] - p->cutval;
	  if (val < 0) 
		{
		  // rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP,Kglobals);
		  rnnEuclideanGen(p->loson,querpoint,points,dimension,numpoints,Kglobals);
		  if (Kglobals->nndist[1] >= val*val) 
			{
			  //rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP,Kglobals);
			  rnnEuclideanGen(p->hison,querpoint,points,dimension,numpoints,Kglobals);
			}
		} 
	  else 
		{
		  // rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP,Kglobals);
		  rnnEuclideanGen(p->hison,querpoint,points,dimension,numpoints,Kglobals);
		  if (Kglobals->nndist[1] >= val*val) 
			{
			  // rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP,Kglobals);
			  rnnEuclideanGen(p->loson,querpoint,points,dimension,numpoints,Kglobals);
			}
		}
	}
}


/*******************************************************************************/

void rnnGeneral(p,querpoint,points,dimension,numpoints,MinkP,Kglobals)

/*******************************************************************************/


	 optkdNode *p;
	 double *querpoint;
	 double **points;
	 int dimension,numpoints,MinkP;
	 optKglobals* Kglobals;
{
  int i,j,k;
  int *perm;              // JG
  int* optfound;          // JG
  double d,thisdist,val,thisx;

  perm = Kglobals->perm; // JG
  optfound = Kglobals->optfound; // JG


  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist = Kglobals->Distance(points,perm[i],querpoint,dimension,MinkP);

      if (optfound[0] < numpoints) {
        PQInsert_kd(thisdist,perm[i],Kglobals->nndist,optfound);
      } else {
        PQreplace_kd(thisdist,Kglobals->nndist,optfound,perm[i]);
      }
    }
  } else {
    val = p->cutval;
    thisx=querpoint[p->discrim];
    if (thisx < val) {
      rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP,Kglobals);
      if (thisx + Kglobals->nndist[1] > val) {
        rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP,Kglobals);
      }
    } else {
      rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP,Kglobals);
      if (thisx - Kglobals->nndist[1] < val) {
        rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP,Kglobals);
      }
    }
  }
}


/******************************************************************************/
/*
  Returns list of neighbours in an array arr[1...numNN], which are the indexes
  in points[] of the nearest numNN points to the given point querpoint[]
*/
/*******************************************************************************/

int *kdOptNNQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)
	 optkdNode *root;
	 double *querpoint, **points;
	 int dimension,numNN,Metric,MinkP;

{
  int j;
  int* optfound;
  double *nndist;

  /* set up found array */
  (*root->Kglobals).optfound = (int *) malloc((numNN+1)*sizeof(int)); /* JG */
  optfound = (*root->Kglobals).optfound;  /* JG */

  optfound[0]=1;  /* for now */

  /* nndist is a priority queue of the distances of the nearest neighbors found */
  nndist = (double *)malloc((numNN+1)*sizeof(double));
  (*root->Kglobals).nndist = nndist; // JG
  for (j=0; j < numNN+1; j++) {
    nndist[j] = DBL_MAX;
  }
   
  switch(Metric) 
	{
	case EUCLIDEAN_GEN :
	  rnnEuclideanGen(root,querpoint,points,dimension,numNN,root->Kglobals);
	  break;
	case EUCLIDEAN : 
	  rnnEuclidean(root,querpoint,points,dimension,numNN,root->Kglobals);
	  break;
	case MANHATTAN : 
	  (*root->Kglobals).Distance=ManhattDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	case L_INFINITY: 
	  (*root->Kglobals).Distance=LInfinityDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	case L_P       : 
	  (*root->Kglobals).Distance=LGeneralDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	}
  free(nndist);
  return(optfound);
}

INDEX_DIST_PAIR *kdOptNNDistQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)

/*******************************************************************************/

	 optkdNode *root;
	 double *querpoint, **points;
	 int dimension,numNN,Metric,MinkP;

{
  int j;
  int* optfound;
  double *nndist;
  INDEX_DIST_PAIR *pairs;

  pairs = (INDEX_DIST_PAIR *)malloc(sizeof(INDEX_DIST_PAIR));
  /* set up found array */
  (*root->Kglobals).optfound = (int *) malloc((numNN+1)*sizeof(int)); /* JG */
  optfound = (*root->Kglobals).optfound;                              /* JG */
  pairs->indx = optfound;                                             /* JG */

  optfound[0]=1;  /* for now */

  /* nndist is a priority queue of the distances of the nearest neighbors found */
  nndist = (double *)malloc((numNN+1)*sizeof(double));
  (*root->Kglobals).nndist = nndist; // JG
  pairs->dist = nndist;

  for (j=0; j < numNN+1; j++) {
    nndist[j] = DBL_MAX;
  }
   
  switch(Metric) 
	{
	case EUCLIDEAN_GEN :
	  rnnEuclideanGen(root,querpoint,points,dimension,numNN,root->Kglobals);
	  break;
	case EUCLIDEAN : 
	  rnnEuclidean(root,querpoint,points,dimension,numNN,root->Kglobals);
	  break;
	case MANHATTAN : 
	  (*root->Kglobals).Distance=ManhattDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	case L_INFINITY: 
	  (*root->Kglobals).Distance=LInfinityDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	case L_P       : 
	  (*root->Kglobals).Distance=LGeneralDist;
	  rnnGeneral(root,querpoint,points,dimension,numNN,MinkP,root->Kglobals);
	  break;
	}
  // free(nndist);
  return pairs;
}


/***************************************************************************/

void KillOptTree(optkdNode *P)

/***************************************************************************/

/*  Kills a kd-tree   */
{
  if (P==NULL) 
	{
	  return;
	} 
  /* free permutation array */
  if (P->Kglobals && P->Kglobals->perm) 
	{
	  free(P->Kglobals->perm);
	  P->Kglobals->perm = NULL;
	}  

  /* recursively down, just to be sure */
  if (P->loson != NULL) 
	{
	  recurs_clean_level++;
	  KillOptTree(P->loson);
	  recurs_clean_level--;
	}

  if (P->hison != NULL)  
	{
	  recurs_clean_level++;
	  KillOptTree(P->hison);
	  recurs_clean_level--;
	}

  

  /*
	Plus any generalised Euclidean bits (James Gunning)
  */
  KillMetricSpace(P->Kglobals);
  
  if (P->Kglobals && recurs_clean_level==0)
	{
	  free(P->Kglobals);
	  P->Kglobals = NULL;
	}

  free(P);
  P=NULL;

}

/*
  Kill wkspace for generalised Euclidean bits (James Gunning)
*/
void KillMetricSpace(optKglobals* Kglobals)
{
  if (Kglobals->distance_metric != NULL)
	{
	  free_dmatrix(Kglobals->distance_metric,0,Kglobals->dim-1,0,Kglobals->dim-1);
	  Kglobals->distance_metric = NULL;
	}
  if (Kglobals->wksp1 != NULL)
	{
	  free_dvector(Kglobals->wksp1,0,Kglobals->dim-1);
	  Kglobals->wksp1 = NULL;
	}
  if (Kglobals->wksp2 != NULL)
	{
	  free_dvector(Kglobals->wksp2,0,Kglobals->dim-1);
	  Kglobals->wksp2 = NULL;
	}
}

/***************************************************************************/

void optInRegion(P,Dimension,Points,RectQuery)

	 /***************************************************************************/
	 /* Determines if the treenode P falls inside the rectangular query         */
	 /* RectQuery.  If so, adds the array index of the point to the found       */
	 /* array.                                                                  */
	 /***************************************************************************/

	 optkdNode *P;
	 int Dimension;
	 double **Points, **RectQuery;

{
  int index,dc,InsideRange;
  int* optfound;           // JG UNSURE ABOUT THIS
  int* perm;               // JG UNSURE ABOUT THIS
  perm = (*P->Kglobals).perm;  // JG UNSURE ABOUT THIS
  optfound = (*P->Kglobals).optfound;   // JG UNSURE ABOUT THIS

  for (index=P->lopt;index<=P->hipt;index++) 
	{
	  InsideRange=1;
	  
	  for (dc=0; dc < Dimension; dc ++) 
		{
		  if ((Points[perm[index]][dc] < RectQuery[dc][0] || Points[perm[index]][dc] > 
			   RectQuery[dc][1])) {   /* P is in the region */
			InsideRange=0;
			break;
		  }
		}
	  if (InsideRange) 
		{
		  (*P->Kglobals).Count++;
		  if (((*P->Kglobals).Count % 4000) == 0) 
			{
			  if (realloc(optfound,(4000+(*P->Kglobals).Count)*sizeof(int)) == NULL)
				{
				  printf("realloc() failure in optInRegion()\n");
				}
			}
		  if (optfound == NULL) 
			{
			  printf("we have a memory problem\n");
			}
		  optfound[(*P->Kglobals).Count] = perm[index];
		}
	}
}

/***************************************************************************/

void optAddRegion(P,Dimension,Points,RectQuery)

	 /***************************************************************************/
	 /* Adds the array index of each point in the bucket the point to the found */
	 /* array.  There is no need to check if the points are in it because we    */
	 /* have proven so already.                                                 */
	 /***************************************************************************/

	 optkdNode *P;
	 int Dimension;
	 double **Points, **RectQuery;

{
  int index,dc;

  int* perm; 
  int* optfound;
  optfound = (*P->Kglobals).optfound;
  perm = (*P->Kglobals).perm;

  for (index=P->lopt;index<=P->hipt;index++) {
    (*P->Kglobals).Count++;
    if (((*P->Kglobals).Count % 4000) == 0) {
	  if (realloc(optfound,(4000+(*P->Kglobals).Count)*sizeof(int)) == NULL)
		{
		  fprintf(stderr,"realloc() failure in optAddRegion()\n");
		}
    }
    optfound[(*P->Kglobals).Count] = perm[index];
  }
}

/***************************************************************************/

int optBoundsIntersectRegion(B,RectQuery,Dimension)

	 /***************************************************************************/
	 /* Returns true iff the hyper-rectangle defined by bounds array B          */
	 /* intersects the rectangular query RectQuery.                             */
	 /***************************************************************************/

	 double *B,**RectQuery;
	 int Dimension;

{
  int dc;

  for (dc=0; dc < Dimension; dc++) {
    if (B[2*dc] > RectQuery[dc][1] || B[2*dc+1] < RectQuery[dc][0]) {
      return(0);
    }
  }
  return(1);
}

/***************************************************************************/

int optBoundsContainsRegion(B,RectQuery,Dimension)

	 /***************************************************************************/
	 /* Returns true iff the hyper-rectangle defined by bounds array B          */
	 /* is completely contained inside the rectangular query RectQuery.         */
	 /***************************************************************************/

	 double *B,**RectQuery;
	 int Dimension;

{
  int dc;

  for (dc=0; dc < Dimension; dc++) {
    if (!(B[2*dc]   >= RectQuery[dc][0] &&
		  B[2*dc+1] <= RectQuery[dc][1])) {
      return(0);
    }
  }
  return(1);
}


/***************************************************************************/

void optRangeSearch(P,Points,Dimension,RectQuery,B)

	 /***************************************************************************/

	 optkdNode *P;
	 double **RectQuery, **Points, *B;
	 int Dimension;


{
  int dc, disc;
  double *BHigh,*BLow;
  

  if (P==NULL) {printf("somehow a null pointer got sent here\n");}

  if (P->bucket) {   
    if (optBoundsContainsRegion(B,RectQuery,Dimension)) {
      optAddRegion(P,Dimension,Points,RectQuery);
    } else {
      optInRegion(P,Dimension,Points,RectQuery);
    }
    return;
  }

  /* Claim: P is not a bucket node */ 
  disc=P->discrim;
  BLow =  (double *)(malloc(2*Dimension*sizeof(double)));
  BHigh = (double *)(malloc(2*Dimension*sizeof(double)));

  if (BLow == NULL || BHigh == NULL) {
	printf("we have a memory error\n");
  }

  /* copy the region B into BLow, BHigh */
  for (dc=0; dc < 2*Dimension; dc++) {
    BLow[dc]  = B[dc];
    BHigh[dc] = B[dc];
  }

  /* Improve the Bounds for the subtrees */
  BLow[2*disc+1] = P->cutval;
  BHigh[2*disc] = P->cutval;

  if (optBoundsIntersectRegion(BLow,RectQuery,Dimension)) {
    optRangeSearch(P->loson,Points,Dimension,RectQuery,BLow);
  }
  free(BLow);
  if (optBoundsIntersectRegion(BHigh,RectQuery,Dimension)) {
    optRangeSearch(P->hison,Points,Dimension,RectQuery,BHigh);
  }
  free(BHigh);
}


/***************************************************************************/

int *kdOptRectQuery(root,Points,Dimension,RectQuery)

	 /***************************************************************************/

	 optkdNode *root;
	 double **RectQuery, **Points;
	 int Dimension;

{
  double *B;
  int dc;
  int* optfound;
  optfound = (*root->Kglobals).optfound;


  B =  (double *)(malloc(2*Dimension*sizeof(double)));
  if (B == NULL) {
    printf("We have a memory problem\n");
  }

  for (dc =0; dc < Dimension; dc++) {
    B[2*dc]   = RectQuery[dc][0];
    B[2*dc+1] = RectQuery[dc][1];
  }

  (*root->Kglobals).Count=0;
  optfound = (int *)(malloc(4000*sizeof(int)));
  if (optfound == NULL) {
    printf("We have a memory problem\n");
  }

  optRangeSearch(root,Points,Dimension,RectQuery,B);
  free(B);
  optfound[0] = (*root->Kglobals).Count;
  return(optfound);
}

