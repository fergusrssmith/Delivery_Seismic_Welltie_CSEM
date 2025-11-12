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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ECLIO/eclio.h"
#include "ECLIO/eclutil.h"
#include "ECLIO/blockops.h"
#include "ECLIO/access.h"
#include "NR/nrutil.h"

/* Simple routines for block calculations based on Eclipse records. 
   More advanced calculations are in cornerpoint.c (in prep) 
                    CD White, LSU, 5/04 */

/*************************************************************************/
int topAct(ECL *ecl, int i, int j)
/*************************************************************************/
{ /* topmost active block in col i j */
  int k;
  int *a = (int *) NULL;
  
  a = (ecl->actnum)[i][j];

  for (k=1;k<=ecl->specGrid.nz;k++) if (a[k]) return (k);

  return (I_MISSING);
}
/*************************************************************************/
int botAct(ECL *ecl, int i, int j)
/*************************************************************************/
{ /* bottommost */
  int k;
  int *a = (int *) NULL;
  
  a = (ecl->actnum)[i][j];

  for (k=ecl->specGrid.nz;k>=1;k--) if (a[k]) return (k);

  return (I_MISSING);
}
/*************************************************************************/
double zTopBlock(ECL *ecl, int i, int j, int ic, int jc, 
		 int skipVoids)
/*************************************************************************/
     /* gets a particlar corner of a topmost block ic,jc = 0,0 for the corner
     nearest the origin, 1,1 for the farthest. */
{
  int kt = I_MISSING;
  double ****z = (double ****) NULL;

  if (skipVoids) kt = topAct(ecl, i, j);
  else kt = 1;

  if (kt != I_MISSING) {
    z = ecl->zCorn;
    return (z[i][j][kt][1 + ic + 2*jc]);
  }
  else return (D_MISSING);
}
/*************************************************************************/
double zBotBlock(ECL *ecl, int i, int j, int ic, int jc, 
		 int skipVoids)
/*************************************************************************/
{
  int kb = I_MISSING;
  double ****z = (double ****) NULL;

  if (skipVoids) kb = botAct(ecl, i, j);
  else kb = ecl->specGrid.nz;

  if (kb != I_MISSING) {
    z = ecl->zCorn;
    return (z[i][j][kb][5 + ic + 2*jc]);
  }
  else return (D_MISSING);
}
/*************************************************************************/
int zBarBlock(ECL *ecl, int skipVoids, 
	      int iBlk, int jBlk, int ic, int jc,
	      double *ztbar, int *ntbar, 
	      double *zbbar, int *nbbar)
/*************************************************************************/
/* accumulates averages of top and bottom coordinate */
{
  double z1 = D_MISSING;
  z1 = zTopBlock(ecl,iBlk, jBlk, ic, jc, skipVoids);
  if (!APPROX(z1,D_MISSING,D_TOL)) {
    *ztbar+=z1;
    *ntbar++;
  }
  z1 = zBotBlock(ecl,iBlk, jBlk, ic, jc, skipVoids);
  if (!APPROX(z1,D_MISSING,D_TOL)) {
    *zbbar+=z1;
    *nbbar++;
  }
  return (0);
}
/*************************************************************************/
void interpVert(double z, double *point, double *c)
/*************************************************************************/
{ /* computes x and y for a z along c */
  double xt=D_MISSING, yt=D_MISSING, zt=D_MISSING;
  double xb=D_MISSING, yb=D_MISSING, zb=D_MISSING;
  double x=D_MISSING, y=D_MISSING;
  double fac=D_MISSING;

  if (c) {
    xt = c[1]; yt=c[2]; zt=c[3];  
    xb = c[4]; yb=c[5]; zb=c[6];
  }
  
  if (xt != D_MISSING &&
      yt != D_MISSING &&
      zt != D_MISSING &&
      xb != D_MISSING &&
      yb != D_MISSING &&
      zb != D_MISSING ) {
    
    if (z == D_MISSING) z = AVE(zt,zb);

    if (!APPROX(zb,zt,D_TOL)) {
      fac = (z - zt) / (zb-zt);
      x = xt + fac*(xb-xt);
      y = yt + fac*(yb-yt);
    }
    else {
      x = AVE(xt,xb);
      y = AVE(yt,yb);
    }
  }
   
  point[1] = x;
  point[2] = y;
  point[3] = z;
}
/*************************************************************************/
int applyPinch(ECL *ecl)
/*************************************************************************/
{ /* voids block on the pinch->pVal criteria */
  int ***a = (int ***) NULL;
  double ****z = (double ****) NULL;
  double *z1 = (double *) NULL;
  PROP *zp = (PROP *) NULL;
  PROP *ap = (PROP *) NULL;
  double pval = 0;
  int i, j, k, l, m;
  int nx = 1, ny=1, nz=1;
  int status = 0;
  int **bnd = (int **) NULL;
  const char routine[] = "applyPinch";
  double dz=0.0;

  if (!ecl->pinch.set) return (0);

  status = checkSpecGrid(ecl,&nx, &ny, &nz);

  ap = fetchProp(ecl,"ACTNUM");
  if (!ap) {
    ap = addPropGrid(ecl, "ACTNUM", 'i');
    if (!ap || !ap->p) {
      status = 1;
      (void) eclError(routine,status, 
		      "Could not allocate storage for ACTNUM\n");
      goto bottom;
    }
  }
  a = (int ***)ap->p;
  
  if (!ap->set) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
		for (k=1;k<=nz; k++) a[i][j][k]=1;
      }
    }
    ap->set = 1;
    (void) setAuxPointers(ecl);
  }
  
  pval = 4.*ecl->pinch.pVal; /* Set to 4* value to simplify later test */

  zp = ecl->zCornp;
  if (!zp->set) {
    status = 2;
    eclError(routine,status,"ZCORN must be set before call\n");
    goto bottom;
  }

  bnd = zp->bnd;
  z = ecl->zCorn;
  a = ecl->actnum;
  for (i=bnd[1][1];i<=bnd[1][2];i++) {
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (k=bnd[3][1];k<=bnd[3][2];k++) {
#if JG_FEATURES
		if (a[i][j][k]) {
		  z1 = z[i][j][k];
		  for (l=1,dz=0.;l<=4;l++) dz+=z1[l+4]-z1[l];
		  if (fabs(dz)<=pval) a[i][j][k] = 0;  /* test only makes sense with computed dz value */
		}
#else
		if (a[i][j][k]) {
		  z1 = z[i][j][k];
		  for (l=1,dz=0.;l<=4;l++) dz+=z1[l+4]-z1[l];
		}
		if (fabs(dz)<=pval) a[i][j][k] = 0;
#endif
	  }
	}
  }

 bottom:
  return (status);  
}
/*************************************************************************/
 PROP *prepareTraces(ECL *ecl, int forceCalc, int skipVoids)
/*************************************************************************/
{
  /* computes coordinates of top and bottom of traces;
     see fetchTraces for more detail 
  */
  int i, j;
  double ***c = (double ***) NULL;
  double ***t = (double ***) NULL;
  double ztbar, zbbar;
  int ntbar, nbbar;
  int nx, ny;
  PROP *tP = (PROP *) NULL;
  const char routine[] = "prepareTraces";
  int status = 0;

  if (!ecl->specGrid.set) {
    status = 1;
    eclError(routine,status,
	     "Must set SPECGRID before calling\n");
    goto bottom;
  }
  
  nx=ecl->specGrid.nx; ny=ecl->specGrid.ny;
  
  if (!ecl->coordp || !ecl->coordp->set) {
    status = 2;
    eclError(routine,status, "Must set COORD before calling\n");
    goto bottom;
  }

  if (!ecl->zCornp || !ecl->zCornp->set) {
    status = 3;
    eclError(routine,status,"Must set ZCORN before calling\n");
    goto bottom;
  }

  tP = fetchProp(ecl, "NEWCOORD"); /* exists */
  if (tP && tP->p && tP->set) {
    if (!forceCalc) return (tP); /* no update requested */
  }
  else {
    tP = addProp(ecl, "NEWCOORD", 'd', 3, 0, nx, 0, ny, 1, 6);
    if (!tP || !tP->p) {
      status = 4;
      eclError(routine,status,"Could not get NEWCOORD\n");
      goto bottom;
    }
  }
  
  t = (double ***) tP->p;
  c = ecl->coord;

  for (i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++) {
      ztbar = 0; zbbar = 0;
      ntbar = 0; nbbar = 0;
      if (i>0) {
		if (j>0) 
		  zBarBlock(ecl,skipVoids,i,j,1,1,
					&ztbar,&ntbar,&zbbar,&nbbar);
		if (j<ny)
		  zBarBlock(ecl,skipVoids,i,j+1,1,0,
					&ztbar,&ntbar,&zbbar,&nbbar);
      }
      if (i<nx) {
		if (j>0) 
		  zBarBlock(ecl,skipVoids,i+1,j,0,1,
					&ztbar,&ntbar,&zbbar,&nbbar);
		if (j<ny)
		  zBarBlock(ecl,skipVoids,i+1,j+1,0,0,
					&ztbar,&ntbar,&zbbar,&nbbar);
      }

      if (ntbar > 0) ztbar /= (double) ntbar;
      else ztbar = D_MISSING;

      if (nbbar > 0) zbbar /= (double) nbbar;
      else zbbar = D_MISSING;
      
      interpVert(ztbar, t[i][j], c[i][j]);
      interpVert(zbbar, &(t[i][j][3]), c[i][j]);
    }
  }
 bottom:

  if (!status) tP->set = 1;
  else {
    if (tP) deleteProp(ecl,"NEWCOORD");
    tP = NULL;
  } 

  return (tP);
}

#if JG_FEATURES

/*************************************************************************/
PROP *prepareTracesAndVolumes(ECL *ecl, int forceCalc, int skipVoids, int layerNumberForTraces)
/*************************************************************************/
{
  /* computes coordinates of top and bottom of traces, plus suitable trace-centred volumes;
     see fetchTraces for more detail 
  */
  int i, j;
  double ***c = (double ***) NULL;
  double ***t = (double ***) NULL;
  double ztbar, zbbar;
  int ntbar, nbbar;
  int nx, ny;
  PROP *tP = (PROP *) NULL;
  PROP *aP = (PROP *) NULL;
  double ***area = (double ***) NULL;
  const char routine[] = "prepareTracesAndVolumes";
  int status = 0;

  if (!ecl->specGrid.set) {
    status = 1;
    eclError(routine,status,
			 "Must set SPECGRID before calling\n");
    goto bottom;
  }
  
  nx=ecl->specGrid.nx; 
  ny=ecl->specGrid.ny;
  
  /* protect against realloc's */
  setAuxPointers(ecl);

  if (!ecl->coordp || !ecl->coordp->set) {
    status = 2;
    eclError(routine,status, "Must set COORD before calling\n");
    goto bottom;
  }

  if (!ecl->zCornp || !ecl->zCornp->set) {
    status = 3;
    eclError(routine,status,"Must set ZCORN before calling\n");
    goto bottom;
  }
    

  tP = fetchProp(ecl, "NEWCOORD"); /* exists */
  if (tP && tP->p && tP->set) {
    if (!forceCalc) return (tP); /* no update requested */
  }
  else {
    tP = addProp(ecl, "NEWCOORD", 'd', 3, 0, nx, 0, ny, 1, 6);
    if (!tP || !tP->p) {
      status = 4;
      eclError(routine,status,"Could not get NEWCOORD\n");
      goto bottom;
    }
  }
  t = (double ***) tP->p;

  /* another realloc protect */
  setAuxPointers(ecl);
  c = ecl->coord;

  for (i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++) {
      ztbar = 0; zbbar = 0;
      ntbar = 0; nbbar = 0;
      if (i>0) {
		if (j>0) 
		  zBarNumberedBlock(ecl,layerNumberForTraces, skipVoids,i,j,1,1,
					&ztbar,&ntbar,&zbbar,&nbbar);
		if (j<ny)
		  zBarNumberedBlock(ecl,layerNumberForTraces, skipVoids,i,j+1,1,0,
					&ztbar,&ntbar,&zbbar,&nbbar);
      }
      if (i<nx) {
		if (j>0) 
		  zBarNumberedBlock(ecl,layerNumberForTraces, skipVoids,i+1,j,0,1,
					&ztbar,&ntbar,&zbbar,&nbbar);
		if (j<ny)
		  zBarNumberedBlock(ecl,layerNumberForTraces, skipVoids,i+1,j+1,0,0,
					&ztbar,&ntbar,&zbbar,&nbbar);
      }

      if (ntbar > 0) ztbar /= (double) ntbar;
      else ztbar = D_MISSING;

      if (nbbar > 0) zbbar /= (double) nbbar;
      else zbbar = D_MISSING;
      
      interpVert(ztbar, t[i][j], c[i][j]);
      interpVert(zbbar, &(t[i][j][3]), c[i][j]);

	}
  }
  fprintf(stderr,"...fetching block areas...\n");
  int** active = fetchActCoord(ecl,0,nx,0,ny);
  if (!active)
	{ 
	  printf("ERROR: Returned null active array in TRACEAREAS\n");
	  return (PROP *) NULL;
	}


  aP = fetchProp(ecl, "TRACEAREAS"); /* exists */
  if (aP && aP->p && aP->set) {
    if (!forceCalc) return (aP); /* no update requested */
  }
  else {
    aP = addProp(ecl, "TRACEAREAS", 'd', 3, 0, nx, 0, ny,0,0);
    if (!aP || !aP->p) {
      status = 4;
      eclError(routine,status,"Could not get TRACEAREAS\n");
      goto bottom;
    }
  }
  area = (double ***) aP->p;

  /* re-fetches as protection against realloc */
  tP = fetchProp(ecl, "NEWCOORD");
  t = (double ***) tP->p;


  for (i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++) {
	  // t[i][j][1,2,3] = x,y,z of trace "top"
	  // t[i][j][4,5,6] = x,y,z of trace "bottom"

	  /* trace-centred areas */
	  area[i][j][0] = midPointTraceCentredArea(ecl,t,active,i,j);
	}
  }
  free_imatrix(active,0,nx,0,ny);

 bottom:

  if (!status) 
	{
	  tP->set = 1;
	  aP->set = 1;
	}
  else 
	{
	  if (tP) deleteProp(ecl,"NEWCOORD");
	  tP = NULL;
	  if (aP) deleteProp(ecl,"TRACEAREAS");
	  aP = NULL;
	}

  return (tP);
}  

#define PRINT_TRIANGLES 0

double midPointTraceCentredArea(ECL *ecl,double ***t, int** active, int i, int j)
{
  /* trace-centred areas */
  double A,A1,A2;
  double x1,x2,x3,y1,y2,y3,dx1,dx2,dy1,dy2;
  int ii,jj,triangle_OK;

  A=0;
  // loop over four blocks surrounding this trace
  for (ii=i-1;ii<=i;ii++)
	{
	  for (jj=j-1;jj<=j;jj++)
		{
		  // Add up 1/4 of the area of each of the 4 blocks.
		  // To get area of each block, split block into two triangles, 
		  // get area (A1,A2) of each from vector cross-product
		  A1=A2=0;
		  // triangle 1
		  triangle_OK=1;
		  if (!getMidPoint(ecl,t,active,ii,jj,&x1,&y1)) triangle_OK=0;
		  if (!getMidPoint(ecl,t,active,ii+1,jj,&x2,&y2)) triangle_OK=0;
		  if (!getMidPoint(ecl,t,active,ii+1,jj+1,&x3,&y3)) triangle_OK=0;
		  if (PRINT_TRIANGLES && triangle_OK)
			{
			  printf("Triangle %d %g %g\n",active[ii][jj],x1,y1);
			  printf("Triangle %d %g %g\n",active[ii+1][jj],x2,y2);
			  printf("Triangle %d %g %g\n",active[ii+1][jj+1],x3,y3);
			  printf("Triangle %d %g %g\n",active[ii][jj],x1,y1);
			  printf("Triangle\n");
			}
		  if (triangle_OK)
			{
			  dx1 = x2-x1;
			  dx2 = x3-x1;
			  dy1 = y2-y1;
			  dy2 = y3-y1;
			  A1 = 0.5*fabs(dx1*dy2-dx2*dy1);
			}
		  // triangle 2
		  triangle_OK=1;
		  if (!getMidPoint(ecl,t,active,ii,jj,&x1,&y1)) triangle_OK=0;
		  if (!getMidPoint(ecl,t,active,ii,jj+1,&x2,&y2)) triangle_OK=0;
		  if (!getMidPoint(ecl,t,active,ii+1,jj+1,&x3,&y3)) triangle_OK=0;
		  if (PRINT_TRIANGLES && triangle_OK)
			{
			  printf("Triangle %d %g %g\n",active[ii][jj],x1,y1);
			  printf("Triangle %d %g %g\n",active[ii][jj+1],x2,y2);
			  printf("Triangle %d %g %g\n",active[ii+1][jj+1],x3,y3);
			  printf("Triangle %d %g %g\n",active[ii][jj],x1,y1);
			  printf("Triangle\n");
			}
		  if (triangle_OK)
			{
			  dx1 = x2-x1;
			  dx2 = x3-x1;
			  dy1 = y2-y1;
			  dy2 = y3-y1;
			  A2 = 0.5*fabs(dx1*dy2-dx2*dy1);
			}
		  // printf("AREA %g %g\n",A1,A2);
		  // take a 1/4 of this quadrilateral into the trace-centred area
		  A+=(A1+A2)/4.0;
		}
	}
  return A;
}
#undef PRINT_TRIANGLES 

/*
  Get x,y's of midpoint of trace i,j. Return 1 if OK, 0 otherwise
*/
int getMidPoint(ECL *ecl,double ***t, int **active, int i,int j,double *x,double *y)
{

  if ((i<0) || (j<0) || (i > ecl->specGrid.nx) || (j>ecl->specGrid.ny)) return 0;
  if (!active[i][j]) return 0;

  // t[i][j][1,2,3] = x,y,z of trace "top"
  // t[i][j][4,5,6] = x,y,z of trace "bottom"

  *x = (t[i][j][1] + t[i][j][4])/2.0;
  *y = (t[i][j][2] + t[i][j][5])/2.0;
	  
  return 1;
}

#endif /* JG_FEATURES */
    
/*************************************************************************/
 PROP *prepareActCoord(ECL *ecl, int forceCalc)
/*************************************************************************/
{
  PROP *aP = (PROP *) NULL;
  PROP *aCoordP = (PROP *) NULL;
  int **aCoord = (int **) NULL;
  PROP *aColP = (PROP *) NULL;
  int **aCol = (int **) NULL;
  int i, j;
  int nx, ny, nz, status = 0;
  char routine[] = "prepareActCoord";

  status = checkSpecGrid(ecl,&nx,&ny,&nz);
  if (status) {
    status = 1;
    (void) eclError(routine,status, "Must set SPECGRID before calling\n");
    goto bottom;
  }

  aP = fetchProp(ecl,"ACTNUM");
  if (!aP || !aP->set) {
    status = 2;
    (void) eclError(routine,status, "Must set ACTNUM before calling\n");
    goto bottom;
  }

  aCoordP = fetchProp(ecl, "ACTCOORD"); /* See if it already exists */
  if (aCoordP && aCoordP->p && aCoordP->set) {
    if (!forceCalc) return (aCoordP); /* may be willing to return it w/o calcs */
  }
  else { /* need to allocate */				       
    aCoordP = addProp(ecl, "ACTCOORD", 'i', 2, 0, nx, 0, ny);

    if (!aCoordP || !aCoordP->p) {
      status = 3;
      (void) eclError(routine, status, "Could not add property ACTCOORD\n");
      goto bottom;
    }
    aCoord = (int **) aCoordP->p;
  }

  aColP = prepareActCol(ecl, forceCalc); /* Need the ACTCOL PROP for this calculation. Might invoke a realloc */

#if JG_FEATURES
  aCoordP = fetchProp(ecl, "ACTCOORD");  /* prepareActCol might have forced a realloc, so get correct pointer */
#endif

  if (!aColP || !aColP->p || !aColP->set) {
    status = 4;
    (void) eclError(routine,status, "Could not set ACTCOL\n");
    goto bottom;
  }

  aCol = (int **) aColP->p;
  if (!aCol) {
    status = 5;
    (void) eclError(routine, status, "Could not access array for ACTCOL\n");
    goto bottom;
  }
  /* FInally, the actual calculation */
  for (i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++) {
      aCoord[i][j]=0;
      if (i>0) {
		if (j>0  && aCol[i  ][j  ]) aCoord[i][j] = 1;
		if (j<ny && aCol[i  ][j+1]) aCoord[i][j] = 1;
      }
      if (i<nx) {
		if (j>0  && aCol[i+1][j  ]) aCoord[i][j] = 1;
		if (j<ny && aCol[i+1][j+1]) aCoord[i][j] = 1;
      }
    }
  }
 bottom:
  if (status) {
    if (aCoordP) deleteProp(ecl, "ACTCOORD");
    aCoordP = NULL;
  }
  else aCoordP->set = 1;

  return (aCoordP);
}
/*************************************************************************/
 PROP *prepareActCol(ECL *ecl, int forceCalc)
/*************************************************************************/
{
  int i, j, k;
  PROP *aP = (PROP *) NULL;
  PROP *acP = (PROP *) NULL;
  int ***a = (int ***) NULL, **ac = (int **) NULL;
  int nx, ny, nz;
  char routine[] = "prepareActCol";
  int status = 0;
  
  status = checkSpecGrid(ecl,&nx,&ny,&nz);
  if (status) {
    status = 1;
    eclError(routine, status, 
	     "checkSpecGrid returned error = %d\n",status);
    goto bottom;
  }

  acP = fetchProp(ecl,"ACTCOL"); /* Try to fetch the array */

  if (acP && acP->p && acP->set) {
    if (!forceCalc) return (acP); /*if exists and no update needed, return */
  }
  else { /* Need to allocate ACTCOL */
    acP = addProp(ecl, "ACTCOL", 'i', 2, 1,nx,1,ny);

	if (!acP || !acP->p) {
      status = 2;
      eclError(routine,status, "Could not fetch or create ACTCOL\n");  
      goto bottom;
    }
    ac = acP->p;
  }

  aP = fetchProp(ecl,"ACTNUM");
  if (!aP || !aP->set || !aP->p) {
    status = 3;
    eclError(routine, status, "Must set ACTNUM before calling\n");
    goto bottom;
  }
  a = (int ***) aP->p;

  /* The actual operation */
  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      ac[i][j]=0;
      for (k=1;k<=nz;k++) {
		if (a[i][j][k]) {
		  ac[i][j] = 1;
		  break;
		}
      }
    }
  }
 bottom:
 
  if (!status) acP->set = 1;
  else {
    if (acP) deleteProp(ecl, "ACTCOL");
    acP = NULL;
  }

  return (acP);
}
/*************************************************************************/
int countActive(ECL *ecl)
/*************************************************************************/
{
  int n;
  PROP *aP = (PROP *) NULL;
  int ***a = (int ***) NULL;
  int i, j, k;
  int nx, ny, nz;
  int status = 0;

  aP = fetchProp(ecl, "ACTNUM");
  if (!aP || !aP->set) return (-1);
  a = (int ***)aP->p;
  if (!a) return (-2);

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) return (-3);

  for (n=0,i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      for (k=1;k<=nz;k++) if (a[i][j][k]) ++n;
    }
  }
  return (n);
}
/*************************************************************************/
int countActCol(ECL *ecl, int forceActCol)
/*************************************************************************/
{
  int n=0;
  PROP *aP = (PROP *) NULL;
  int **a = (int **) NULL;
  int i, j;
  int nx, ny, nz;
  int status = 0;

  
  aP = prepareActCol(ecl, forceActCol);
  if (!aP || !aP->set) return (-1);
  a = (int **)aP->p;
  if (!a) return (-2);

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) return (-3);

  for (n=0,i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++)  if (a[i][j]) ++n;
  }
  return (n);
}  
/*************************************************************************/
int countActCoord(ECL *ecl, int forceActCoord)
/*************************************************************************/
{
  int n;
  PROP *aP = (PROP *) NULL;
  int **a = (int **) NULL;
  int i, j;
  int nx, ny, nz;
  int status = 0;

  
  aP = prepareActCoord(ecl, forceActCoord);
  if (!aP || !aP->set) return (-1);
  a = (int **)aP->p;
  if (!a) return (-2);

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) return (-3);

  for (n=0,i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++)  if (a[i][j]) ++n;
  }
  return (n);
}  
/*************************************************************************/
int rotateGridMap(ECL *ecl)
/*************************************************************************/
{
  /* uses a jacobian computed from MAPAXES to transform 
     grid -> map coordinates 

     Detects need for and may apply three types of transforms
        if mapUnits != gridUnit, then rescales
	if off-diagonal elements of Jacobian != 0, then rotates
	if origin != 0, then translates 
  */

  const char routine[] = "rotateGridMap";
  double **a = NULL;
  int didRotate = 0, didTranslate = 0, didScale = 0;
  double **jacobian = NULL;
  double det = 0;
  double *tmp = NULL;
  double x0, y0;
  int status = 0;
  int i,j;

  /* Check units and rescale if needed */
  if (ecl->mapUnits.set && ecl->gridUnit.set) {
    if(!rescaleGridMap(ecl)) didScale = 1;
  }
  else
    (void) eclError(routine,1, "No rescaling applied; units not set\n");

  /* Compute rotation and apply if needed */
  if (ecl->mapAxes.set) {
    if (!rotateCoord(ecl, &jacobian, &det)) didRotate = 1;
  }
  else
    (void) eclError(routine,2,"No rotation applied; mapAxes not set\n");

  /* Is translation needed ? */
  if (ecl->mapAxes.set) {
    if (!translateCoord(ecl)) didTranslate = 1;
  }
  else
    (void) eclError(routine,3,"No translation applied; mapAxes not set\n");


  if (didRotate) { /* rotate mapAxes now */
    tmp = dvector(1,3);
    if (!tmp) {
      if (jacobian) free_dmatrix(jacobian, 1,3,1,3);
      return (1);
    } 
    a = ecl->mapAxes.axes;
    for (i=1;i<3;i++) {
      status = dmvMult(1,2,1,2,jacobian,1,2,a[i],tmp);
      for (j=1;j<=2;j++) {
	a[i][j] = tmp[j];
      }
    }
  }

  if (didTranslate) { /* translate mapAxes */
    if (didRotate && det < 0.) ecl->mapAxes.isReversed = 1;
    if (!a) a = ecl->mapAxes.axes;
    x0 = a[2][1]; y0 = a[2][2];
    for (i=1;i<=3;i++) {
      a[i][1] -= x0;
      a[i][2] -= y0;
    }
    if (ecl->mapAxes.isReversed) a[3][2] += 2. * y0;
  }

 bottom:
  if (jacobian) free_dmatrix(jacobian,1,3,1,3);
  if (tmp) free_dvector(tmp,1,3);  

  return (0);
}
/*************************************************************************/
int rescaleGridMap (ECL *ecl)
/*************************************************************************/
     /* Uses Eclipse MAPUNITS and GRIDUNIT records to 
	rescale grid to map units */
{
  GRIDUNIT *gup = (GRIDUNIT *) NULL;
  MAPUNITS *mup = (MAPUNITS *) NULL; 
  char *gu = (char *) NULL;
  char *mu = (char *) NULL;
  const char routine[] = "rescaleGridMap";
  double g2m = 1.;
  int status = 0;

  mup = &ecl->mapUnits; 
  if (!mup || !mup->set) {
    status = 1;
    (void) eclError(routine, status, "Must set mapUnits first\n");
    goto bottom;
  }

  mu = mup->units;
  if (!mu) {
    status = 2;
    (void) eclError(routine,status, "Could not access mapUnits.units\n");
    goto bottom;
  }

  gup = &ecl->gridUnit; 
  if (!gup || !gup->set) {
    status = 3;
    (void) eclError(routine,status, "Must set gridUnit first\n");
    goto bottom;
  }

  gu = gup->units;
  if (!gu) {
    status = 4;
    (void) eclError(routine,status, "Could not access gridUnit.units\n");
    goto bottom;
  }

  /* Convert units if needed) */
  if (strncmp(mu, gu,strlen(mu))) return (0);

  if (!strcmp(mu, "METRES")) {
    if (!strcmp(gu, "METRES")) return (0);
    if (!strcmp(gu, "FEET")) g2m = 0.3048;
    else { 
      status = 5;
      (void) eclError(routine, status, "Unrecongnized units for grid: %s\n",gu);
      goto bottom;
    }
  }

  else if (!strcmp(mu, "FEET")) {
    if (!strcmp(gu, "FEET")) return (0);
    if (!strcmp(gu, "METRES")) g2m = 1./0.3048;
    else {
      status = 6;
      eclError(routine, status, "Unrecognized units for grid: %s\n", gu);
      goto bottom;
    }
  }
  else {
    status = 7;
    eclError(routine, status, "Unrecognized units for map: %s\n", mu);
    goto bottom;
  }

  (void) transProp(ecl->zCornp, g2m, 0.);
  (void) transProp(ecl->coordp, g2m, 0.);
  
 bottom:
  if (!status) strcpy(gup->units,mu); /* reset grid units if scaled */

  return (status);
}
/*************************************************************************/
int rotateCoord(ECL *ecl, double ***jacP, double *det)
/*************************************************************************/
{
  double d = 0;
  int i,j,k,r;
  double *tmp = (double *) NULL;
  char routine[] = "rotateCoord";
  int nx, ny, nz;
  int status=0;
  double **jac = NULL;
  double ***cp = (double ***) NULL;
  double *c = (double *) NULL;
  MAPAXES *ma = (MAPAXES *) NULL;
  double **a = (double **) NULL;
  double da[4];

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) {
    status = 1;
    (void) eclError(routine, status, "Must set SPECGRID before call\n");
    goto bottom;
  }

  ma = &ecl->mapAxes;
  if (!ma || !ma->set) {
    status = 2;
    (void) eclError(routine, status, "Could not access MAPAXES\n");
    goto bottom;
  }

  a = ma->axes;
  if (!a) {
    status = 3;
    (void) eclError(routine, status, "Could not locate mapAxes.axes\n");
    goto bottom;
  }

  /* First compute factors to normalize the mapaxis lengths 
     This scaling replaces the earlier hard-wiring and rescaling on 
     the inverse of the square root of the Jacobian, and is intended to
     correctly treat strangely scaled mapAxes records */
  da[1]=sqrt((a[3][1]-a[2][1])*(a[3][1]-a[2][1]) +
	     (a[3][2]-a[2][2])*(a[3][2]-a[2][2])  );
  da[2]=sqrt((a[1][1]-a[2][1])*(a[1][1]-a[2][1]) +
	     (a[1][2]-a[2][2])*(a[1][2]-a[2][2])  );

  if (APPROX(da[1],0.0,D_TOL) || APPROX(da[2],0.0,D_TOL)) {
    status = 4;
    (void) eclError(routine,status, "Axes not sufficiently distinct:\n"
		    "XY1: %10.2lf %10.2lf\n"
		    "XY2: %10.2lf %10.2lf\n"
		    "XY3: %10.2lf %10.2lf\n",
		    a[1][1],a[1][2],a[2][1],a[2][2],a[3][1],a[3][2]);
    goto bottom;
  }
  
  da[1] = 1./da[1]; da[2]=1./da[2];

  /* Compute scaled differences needed for Jacobian */
  jac = dmatrix(1,3,1,3);
  for (i=1;i<=2;i++) {
    jac[1][i] = (a[3][i] - a[2][i])*da[1];
    jac[2][i] = (a[1][i] - a[2][i])*da[2];
    jac[3][i] = 0.;
  }
  jac[1][3] = 0.;
  jac[2][3] = 0.;
  jac[3][3] = 1.; /*strictly speaking, should have jac[3][3]=-1
		    to preserve right-handedness (and keep det(jac)
		    positive */
    
  /* Check Jacobian */
  d = det3mat(1, 3, 1, 3, jac);
  *det = d;
  if (APPROX(d, 0.0, D_TOL)) {
    status = 5;
    (void) eclError(routine, status, 
		    "Singular transform, det(J) = %lf\n", d);
    goto bottom;
  }

  /* rotate if jacobian is not diagonal || diagonal is not identity 
     This handles reflections as well as rotations */
  if (APPROX(jac[1][2],0.0,D_TOL) && APPROX(jac[2][1],0.0,D_TOL) &&
      APPROX(jac[1][1],1.0,D_TOL) && APPROX(jac[2][2],1.0,D_TOL)) {
    status = 6;
    (void)eclError(routine,status,  /* no rotation / reflection needed */
		   "mapAxes is loaded but transformation is an identity matrix\n"
		   "\t...transformation skipped\n");
    goto bottom;
  }

  else { /* rotate or reflect */
    tmp = dvector(1,3);
    if (!tmp) {
      status = 7;
      (void)eclError(routine,status, "Could not allocate tmp\n");
      goto bottom;
    }

    cp = ecl->coord;
    if (!cp ) {
      status = 8;
      (void)eclError(routine,status, "Could not locate COORD\n");
      goto bottom;
    }

    for (i=0;i<=nx;i++) { /* rotate and translate COORD lines */
      for (j=0;j<=ny;j++) {
	c = cp[i][j];
	for (k=0;k<=3;k+=3) { /* q-loop handles top and bottom */
	  
	  status = dmvMult(1, 3, 1, 3, jac, 1, 3, &c[k], tmp); /*rotate*/
	  if (status) {
	    status = 9;
	    (void) eclError(routine,status,"Unconfomable matrices in rotation\n");
	    goto bottom;
	  }
	  
	  for(r=1;r<=3;r++) c[k+r]=tmp[r]; /* could use r<=2 with no
					      change in this routine
					      b.c. jac(i,3)=jac(3,i) = 0
					      except jac(3,3) = 1 */
	}
      }
    }
  }

 bottom:
  if (tmp) free_dvector(tmp, 1, 3);
  *jacP = jac;
  return (status);
}
/*************************************************************************/
int translateCoord(ECL *ecl)
/*************************************************************************/
{
  int i, j, k;
  double ***cp = NULL;
  MAPAXES *ma = NULL;
  double **a = NULL;
  double *c = NULL;
  int nx, ny, nz, status = 0;
  char routine[] = "translateCoord";
  double x0, y0;

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) {
    status = 1;
    (void) eclError(routine,1,"Must set SPECGRID before call");
    goto bottom;
  }

  cp = ecl->coord;
  if (!cp) {
    status = 2;
    (void) eclError(routine,status,"Could not access COORD\n");
    goto bottom;
  }
  
  ma = &ecl->mapAxes;
  if (!ma || !ma->set){
    status = 3;
    (void) eclError(routine,status, "Could not access MAPAXES\n");
    goto bottom;
  }

  a = ma->axes;
  if (!a) {
    status = 4;
    (void) eclError(routine,status, "Could not access mapAxes.axes\n");
    goto bottom;
  }
  
  /* Apply translation if origin is not at zero */
  if (!APPROX(a[2][1],0.0,D_TOL) || !APPROX(a[2][2],0.0,D_TOL)) {
    for (i=0;i<=nx;i++) { 
      for (j=0;j<=ny;j++) {
	c = cp[i][j];
	for (k=0;k<=3;k+=3) { 
	  c[k+1]+= a[2][1];  
	  c[k+2]+= a[2][2];
	}
      }
    }
  }
 bottom:
  return (status);
}
/*************************************************************************/
 PROP *prepareSegCol(ECL *ecl, int forceCalc)
/*************************************************************************/
{
  int i, j, k;
  PROP *aP = NULL;
  PROP *acP = NULL;
  int ***a = NULL, **ac = NULL;
  int aij = I_MISSING, aijk = I_MISSING;
  int nx=I_MISSING, ny=I_MISSING, nz=I_MISSING;
  char routine[] = "prepareSegCol";
  int status=0;
  
  status = checkSpecGrid(ecl,&nx,&ny,&nz);
  if (status) {
    status = 1;
    eclError(routine, status, 
	     "checkSpecGrid returned error = %d\n",status);
    goto bottom;
  }

  acP = fetchProp(ecl,"SEGCOL");
  if (acP && acP->p && acP->set) { /* exists */
    if (!forceCalc) goto bottom;  /* no update needed */
  }
  else { /* Allocate */
    acP = addProp(ecl, "SEGCOL", 'i', 2, 1,nx,1,ny);
    if (!acP || !acP->p) {
      status = 2;
      eclError(routine,status, "Could not fetch or create SEGCOL\n");  
      goto bottom;
    }
    ac = (int **) acP->p;
  }

  aP = fetchProp(ecl,"SEGMENTS");
  if (!aP || !aP->set || !aP->p) {
    status = 3;
    eclError(routine, status, "Must set SEGMENTS before calling\n");
    goto bottom;
  }
  a = (int ***) aP->p;

  /* the actual calc */

#if JG_FEATURES

  for (i=1;i<=nx;i++) 
	{
	  for (j=1;j<=ny;j++) 
		{
		  aij=0;
		  for (k=1;k<=nz;k++) 
			{
			  aijk = a[i][j][k];
			  if (aijk > 0) 
				{
				  if (aij == 0) 
					{
					  aij = aijk;
					}
				  else if (aij != aijk && (aij>=0) && (aijk>=0)) 
					{
					  /* CW: 
						 check this - I think we should take the segment 
						 number of the non-negative labels
						 if a mixture is supplied, e.g. {-999,-999,2,2,2,2,}, 
						 assuming -ve means (don't know|too lazy)
					  */
					  ac[i][j]=-1;
					  break;
					}
				}
			}
		  ac[i][j]=aij;
		}
	}

#else  

  for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
      aij=0;
      for (k=1;k<=nz;k++) {
		aijk = a[i][j][k];
		if (aijk) {
		  if (aij == 0) {
			aij = aijk;
		  }
          else if (aij != aijk) {
			ac[i][j]=-1;
			break;
		  }
		}
      }
      ac[i][j]=aij;
    }
  }
#endif

 bottom:
  if (status) {
    if (acP) deleteProp(ecl, "SEGCOL");
    acP = NULL;
  }
  else acP->set = 1;

  return (acP);
}
/*************************************************************************/
 PROP *prepareTrace4(ECL *ecl, int forceCalc, int skipVoids)
/*************************************************************************/
{
  /* computes coordinates of top and bottom of all blocks neighboring traces;
     see fetchTrace4 for more detail 
  */
  int i, j, it, status=0;
  double ****t=NULL, ***c=NULL;
  int nx=I_MISSING, ny=I_MISSING;
  PROP *tP=NULL;
  const char routine[] = "prepareTrace4";

  if (!ecl->specGrid.set) {
    status = 1;
    eclError(routine,status,
	     "Must set SPECGRID before calling\n");
    goto bottom;
  }
  nx=ecl->specGrid.nx; ny=ecl->specGrid.ny;

  if (!ecl->coordp || !ecl->coordp->set) {
    status = 2;
    eclError(routine,status, "Must set COORD before calling\n");
    goto bottom;
  }

  if (!ecl->zCornp || !ecl->zCornp->set) {
    status = 3;
    eclError(routine,status,"Must set ZCORN before calling\n");
    goto bottom;
  }

  c = ecl->coord;

  tP = fetchProp(ecl, "NEWCOOR4"); /* exists */
  if (tP && tP->p && tP->set) {
    if (!forceCalc) return (tP); /* no update requested */
  }
  else {
    tP = addProp(ecl, "NEWCOOR4", 'd', 4, 0, nx, 0, ny, 1, 4, 1, 6);
    if (!tP || !tP->p) {
      status = 4;
      eclError(routine,status,"Could not get NEWCOOR4\n");
      goto bottom;
    }
  }
  
  t = (double ****) tP->p;

  for (i=0;i<=nx;i++) {
    for (j=0;j<=ny;j++) {
      for (it = 1; it<=4;it++)	{
	status = oneBlockEdge(ecl, nx, ny, i, j, it, c, t);
	if (status > 2) {
	  status = 5;
	  eclError(routine, status, "Bounds violation in oneBlockEdge.\n");
	  goto bottom;
	}
	else status = 0;
      }
    }
  }
 bottom:
  if (status) {
    if (tP) deleteProp(ecl,tP->name);
    tP = NULL;
  }
  else tP->set = 1;
 
  return (tP);
}
int oneBlockEdge(ECL *ecl, int nx, int ny, int i, int j, int it, double ***c,
		 double ****t)
{
  double z=D_MISSING;
  int ic=I_MISSING, jc=I_MISSING, iCoord=I_MISSING;
  int iBlk=I_MISSING, jBlk=I_MISSING;

  if (i < 0 || i > nx || j < 0 || j > ny || it < 1 || it > 4) return (2);

  if ((it == 1 && (i ==  0 || j ==  0)) ||
      (it == 2 && (i == nx || j ==  0)) ||
      (it == 3 && (i ==  0 || j == ny)) ||
      (it == 4 && (i == nx || j == ny))   ){
    for (iCoord=1;iCoord<=6;iCoord++) t[i][j][it][iCoord]=0.;
    return (1);
  }

  switch (it) {
  case 1:
    ic = 1; jc = 1; 
    iBlk = i; jBlk = j;
    break;
  case 2:
    ic = 0; jc = 1; 
    iBlk = i+1; jBlk = j;
    break;
  case 3:
    ic = 1; jc = 0; 
    iBlk = i; jBlk = j+1;
    break;
  case 4:
    ic = 0; jc = 0; 
    iBlk = i+1; jBlk = j+1;
    break;
  default:
    break;
  }
  z = zTopBlock(ecl, iBlk, jBlk, ic, jc, 1);
  interpVert(z, t[i][j][it], c[i][j]);
  z = zBotBlock(ecl, iBlk, jBlk, ic, jc, 1);
  interpVert(z, &(t[i][j][it][3]), c[i][j]);
  return (0);
}

#ifdef JG_FEATURES

/*************************************************************************/
int zBarNumberedBlock(ECL *ecl, int kBlock, int skipVoids,
					  int iBlk, int jBlk, int ic, int jc,
					  double *ztbar, int *ntbar, 
					  double *zbbar, int *nbbar)
/*************************************************************************/
/* accumulates averages of top and bottom coordinate of layer of index kBlock */
{
  if (kBlock<0)
	{
	  /* default complete trace if kBlock<0 */
	  return zBarBlock(ecl,skipVoids,iBlk,jBlk,ic,jc,ztbar,ntbar,zbbar,nbbar);
	}
  /* else more precise locations... */

  double z1 = D_MISSING;
  z1 = zNumberedBlock(ecl,iBlk, jBlk, ic, jc, kBlock, skipVoids);
  if (!APPROX(z1,D_MISSING,D_TOL)) {
    (*ztbar)+=z1;
    (*ntbar)++;
  }
  z1 = zNumberedBlock(ecl,iBlk, jBlk, ic, jc, kBlock+1, skipVoids);
  if (!APPROX(z1,D_MISSING,D_TOL)) {
    (*zbbar)+=z1;
    (*nbbar)++;
  }
  return (0);
}


/*************************************************************************/
double zNumberedBlock(ECL *ecl, int i, int j, int ic, int jc, int kBlock,
		 int skipVoids)
/*************************************************************************/
{
  int kb = I_MISSING;
  double ****z = (double ****) NULL;

  if (kBlock < 1)
	{
	  return zBotBlock(ecl,i,j,ic,jc,skipVoids);
	}
  else if (kBlock > ecl->specGrid.nz)
	{
	  return zTopBlock(ecl,i,j,ic,jc,skipVoids);
	}
  else
	{
	  if (skipVoids) 
		{
		  if ((ecl->actnum)[i][j][kBlock]) kb  = kBlock;
		}
	  else
		{
		  kb  = kBlock;
		}
	}
  if (kb != I_MISSING) 
	{
	  z = ecl->zCorn;
	  return (z[i][j][kb][5 + ic + 2*jc]);
	}
  else return (D_MISSING);
}



/* JG utility for units change */

int flipGridAndMapUnits(ECL *ecl)
{
  int i,j,k,c,Nx,Ny,Nz;
  double ****z;
  double METRES_PER_FOOT=0.3048;  // m per foot 
  int doMap=0, doGrid=0;

  // MAPUNITS
  if (ecl->mapUnits.set)
	{
	  if (strncmp(ecl->mapUnits.units,"FEET",4))
		{
		  // not feet. skip
		  doMap=0;
		}
	  else
		{
		  // mapUnits = feet
		  doMap=1;
		}
	}
  else
	{
	  // assume we can blaze away
	  doMap=1;
	}

  // GRIDUNITS
  if (ecl->gridUnit.set)
	{
	  if (strncmp(ecl->gridUnit.units,"FEET",4))
		{
		  // gridUnits not feet. Skip
		  doGrid=0;
		}
	  else
		{
		  // gridUnits = feet
		  doGrid=1;
		}
	}
  else
	{
	  // assume we can blaze away
	  doGrid=1;
	}

  if (doGrid)
	{
	  /* Eclipse dimensions */
	  Nx = getNX(ecl);
	  Ny = getNY(ecl);
	  Nz = getNZ(ecl); 

	  // ZCORN
	  PROP *from=NULL;
	  from = fetchProp(ecl,"ZCORN");
	  if (!from)
		{
		  fprintf(stderr,"flipGridAndMapUnits(): ZCORN not available");
		  return EXIT_FAILURE;
		}

	  z = ((double ****)(from->p));  // should be eclOut->zCorn

	  for(i=1;i<=Nx;i++)
		{
		  for (j=1;j<=Ny;j++)
			{
			  for (k=1;k<=Nz;k++)
				{
				  for (c=1;c<=8;c++)
					{
					  z[i][j][k][c] *= METRES_PER_FOOT;
					}
				}
			}
		}

	  // COORD
	  for(i=0;i<=Nx;i++)
		{
		  for (j=0;j<=Ny;j++)
			{
			  for (k=1;k<=6;k++)
				{
				  ecl->coord[i][j][k] *= METRES_PER_FOOT;
				}
			}
		}

	  // GRIDUNIT
	  sprintf(&(ecl->gridUnit.units[0]),"METRES");
	  ecl->gridUnit.set=1;
	}

  if (doMap)
	{
	  // MAPAXES
	  for (i=1;i<=3;i++) 
		{
		  for (j=1;j<=2;j++) 
			{
			  ecl->mapAxes.axes[i][j] *= METRES_PER_FOOT;
			}
		}
	  // MAPUNITS
	  sprintf(&(ecl->mapUnits.units[0]),"METRES");
	  ecl->mapUnits.set = 1;
	}

  return EXIT_SUCCESS;
}
#endif
