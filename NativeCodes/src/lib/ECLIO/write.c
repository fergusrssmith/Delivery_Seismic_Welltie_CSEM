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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include "ECLIO/eclio.h"
#include "ECLIO/eclutil.h"
#include "ECLIO/access.h"
#include "NR/nrutil.h"
#include "ECLIO/write.h"

/* Reads from a model file fromFile to toFile
		echoes fromFile until a key on keysFromNew is encountered
		inserts record corresponding to keysFromNew[] into toFile
		then resumes echoing until the next key is found
		until an nKeysFromNew are found
		finally echoes remainder of fromFile to toFile
		keysFromNew must be in the same order as keywords are encountered in
			fromFile <--- ****

	Returns various error conditions if keys are not found, etc.
			CD White, LSU 5/04
*/

int wRecord(FILE *f, ECL *ecl, char *name)
     /* Writes a named record or a PROP array headers and other ECLIPSEish formatting */
{
  int status = 0;

  if (!f) status = 3;
  else if (!strcmp(name,"SPECGRID"))      status = wSpecGrid(f, ecl);
  else if (!strcmp(name,"ECHO"))   fprintf(f,"\n%s\n\n",name);
  else if (!strcmp(name,"NOECHO")) fprintf(f,"\n%s\n\n",name);
  else if (!strcmp(name,"GRIDUNIT")) status = wGridUnit(f, ecl);
  else if (!strcmp(name,"PINCH"))    status = wPinch(f, ecl); 
  else if (!strcmp(name,"MAPAXES"))  status = wMapAxes(f, ecl);
  else if (!strcmp(name,"MAPUNITS")) status = wMapUnits(f, ecl);
  else if (!strcmp(name, "NUMRES")) status = wNumRes(f,ecl);
  else if (!strcmp(name, "DIMENS")) status = wDimens(f,ecl);
  else if (!strcmp(name,"COORDSYS")) status = wCoordSys(f, ecl);
  else if (!strncmp(name,"COORD",5)) status = wCoord(f, ecl, name); /* allows COORD4, etc */
  else if (!strcmp(name,"ZCORN")) status = wZCorn(f, ecl);
  else status = wProp(f, ecl, name);
  
  return (status);
}
int wProp(FILE *f,ECL *ecl, char *name)
{
  int status = 0;
  int nD = I_MISSING;
  int **b = NULL;
  char pt = '\0';
  PROP *pP = NULL;

  GENERIC_POINTER x = NULL;

  pP = fetchProp(ecl, name);
  if (!pP) return(1);
  
  x = pP->p;
  if (!x) return (2);
  
  pt = pP->ptype; nD = pP->nInd; b = pP->bnd;

  if (!f) return (3);

  fprintf (f, "%s \n", pP->name);

  switch(pt) {
  case 'i':
    status = wPropArrayI(f, x, nD, b);
    break;
  case 'f':
    status = wPropArrayF(f, x, nD, b);
    break;
  case 'd':
    status = wPropArrayD(f, x, nD, b);
    break;
  default:
    status = 4;
    break;
  }
  
  fprintf (f, RECORD_END_FORMAT);

  return (status);
}
  

/*************************************************************************/
int wSpecGrid(FILE * f, ECL *ecl)
/*************************************************************************/
{
  char buffer[LINE_LENGTH];
  int status=0;
  int nx, ny, nz;
  SPECGRID *rP = NULL;

  rP = &(ecl->specGrid);
  if (!rP) status = 1;
  else if (!f) status = 3;
  else
    if (checkSpecGrid(ecl, &nx, &ny, &nz)) status = 2;
    else {
      fprintf(f,
	      "SPECGRID \n" 
	      "%d %d %d %d %c", 
	      nx,ny,nz,rP->numRes, rP->cylindrical);
      fprintf(f,RECORD_END_FORMAT);
    }
  return (status);
}
/*************************************************************************/
int wDimens(FILE * f, ECL *ecl)
/*************************************************************************/
{
  char buffer[LINE_LENGTH];
  int status=0;
  int nx, ny, nz;
  SPECGRID *rP = NULL;

  rP = &(ecl->specGrid);
  if (!rP) status = 1;
  else if (!f) status = 3;
  else
    if (checkSpecGrid(ecl, &nx, &ny, &nz)) status = 2;
    else {
      fprintf(f,
	      "DIMENS \n" 
	      "%d %d %d", 
	      nx,ny,nz);
      fprintf(f,RECORD_END_FORMAT);
    }
  return (status);
}
/*************************************************************************/
int wNumRes(FILE * f, ECL *ecl)
/*************************************************************************/
{
  char buffer[LINE_LENGTH];
  int status=0;
  int nx, ny, nz;
  SPECGRID *rP = NULL;

  rP = &(ecl->specGrid);
  if (!rP) status = 1;
  else if (!f) status = 3;
  else
    if (checkSpecGrid(ecl, &nx, &ny, &nz)) status = 2;
    else {
      fprintf(f,
	      "NUMRES \n" 
	      "%d", 
	      rP->numRes);
      fprintf(f,RECORD_END_FORMAT);
    }
  return (status);
}
/*************************************************************************/
int wPinch(FILE *f, ECL *ecl)
/*************************************************************************/
{
  PINCH *rP = NULL;
  int status = 0;

  rP = &(ecl->pinch);
  if (!rP) status = 1;
  else if (!f) status = 3;
  else {
    fprintf(f,
	       "PINCH \n"
	       "%lf %12s %12.2lE %12s %12s",
	       rP->pVal, rP->gap, rP->pinchPV, rP->topbot, rP->tMeth);
    fprintf(f,RECORD_END_FORMAT);
  }
  return (status);
}
/*************************************************************************/
int wPropHeader(FILE *f, ECL *ecl, char *name)
/*************************************************************************/
{
  PROP *pp = NULL;
  int status = 0;
  
  pp = fetchProp(ecl,name);

  if (!pp) status = 1;
  else if (!f) status = 3;
  else fprintf (f, "%s \n",name);
  
  return (status);
}
/*************************************************************************/
int wCoordSys(FILE *f, ECL *ecl)
/*************************************************************************/
{
  COORDSYS *cs = NULL;
  int numRes = I_MISSING;
  int i,status = 0;

  if (!ecl->specGrid.set) status = 1;
  if (!f) status = 3;
  else 
	{
	  numRes = ecl->specGrid.numRes;
	  if (!ecl->coordSys.set) status = 2;
	  else 
		{
		  cs = &ecl->coordSys;
		  fprintf(f,"COORDSYS \n");
		  if (numRes <1) status = 4;
		  else 
			{
			  for (i=1;i<=numRes;i++) 
				{
				  fprintf(f,"\t%4d\t%4d\t%10s\t%10s\t/ \n",cs->k1[i],cs->k2[i],cs->comp[i], cs->join[i]);
				}
			}
		}
	}
  return (status);
}
/*************************************************************************/
int wGridUnit(FILE *f, ECL *ecl)
/*************************************************************************/
{
  int status = 0;
  GRIDUNIT *rP = NULL;

  rP = &(ecl->gridUnit);
  if (!rP) status = 1;
  else if (!rP->set) status = 2;
  else if (!f) status = 3;
  else {
    fprintf(f,
	    "GRIDUNIT \n"
	    "\t%10s\t%10s",rP->units, rP->map);
    fprintf(f,RECORD_END_FORMAT);
  }
  return (status);
}
/*************************************************************************/
int wMapUnits(FILE *f, ECL *ecl)
/*************************************************************************/
{
  int status=0;
  MAPUNITS *rP = NULL;

  rP = &(ecl->mapUnits);
  if (!rP) status = 1;
  else if (!rP->set) status = 2;
  else if (!f) status = 3;
  else {
    fprintf(f,
	    "MAPUNITS\n"
	    "\t%10s",
	    rP->units);
    fprintf(f,RECORD_END_FORMAT);
  }
  return (status);
}
/*************************************************************************/
int wMapAxes(FILE *f, ECL *ecl)
/*************************************************************************/
{
  int status = 0;
  MAPAXES *rP = NULL;
  double **axes = NULL;

  rP = &(ecl->mapAxes);
  if (!rP) status = 1;
  else if (!rP->set) status = 2;
  else if (!(axes = rP->axes)) status = 4;
  else if (!f) status = 3;
  else {
    fprintf(f,
	    "MAPAXES\n"
	    "\t%12.3lf\t%12.3lf\n"
	    "\t%12.3lf\t%12.3lf\n"
	    "\t%12.3lf\t%12.3lf",
	    axes[1][1],axes[1][2],
	    axes[2][1],axes[2][2],
	    axes[3][1],axes[3][2]);
    fprintf(f,RECORD_END_FORMAT);
  }
  return (status);
}

int wZCorn(FILE *f, ECL *ecl)
{
  int i, j, k, p, q, r;
  int nx, ny, nz;
  int nOnLine;
  int status = 0;
  double ****z = NULL;
  PROP *zP = NULL;

  int nPerLine = DOUBLE_ON_LINE ;
  char fmt[] = DOUBLE_FORMAT ;

  nOnLine = 0;
  
  if (checkSpecGrid(ecl,&nx,&ny,&nz)) return (1);
 
  fprintf(f, "ZCORN\n");

  zP = fetchProp(ecl,"ZCORN");
  if (!zP) return (2);
  if (!zP->set) return (4);
  if (!f) return (3);

  z = (double ****) zP->p;
  if (!z) return (5);

  for (k=1;k<=nz;k++) {  /*z-direction scans slowest */
    for (p=0;p<=4;p+=4) { /* ...then top and bottom of layer k */
      for (j=1;j<=ny;j++) { /* then along each row */
		for (q=0;q<=2;q+=2) { /* do +- corners in j for each sweep */
		  for (i=1;i<=nx;i++) { /* write corners in i for each block */
			for (r=1;r<=2;r++) { /* get each of two corners on this line */
			  nOnLine++;
			  if (nOnLine > nPerLine) {
				fprintf(f,"\n");
				nOnLine = 1;
			  }
			  fprintf (f,fmt,z[i][j][k][p+q+r]);
			}}}}}}
  
  fprintf (f,RECORD_END_FORMAT);
  return (status);

}
/*************************************************************************/
int wCoord(FILE *f, ECL *ecl, char *name)
/*************************************************************************/
{
  int i, j, k;
  int nx, ny, nz;
  int nOnLine;
  int status = 0;

  PROP *cP = NULL;
  double ***c = NULL;

  int nPerLine = 6 ;
  char fmt[] = " %12.2lf" ;

  if (checkSpecGrid(ecl, &nx, &ny, &nz)) return (1);

  cP = fetchProp(ecl, name);
  if (!cP) return (4);
  else if (!cP->set) return (2);
  else if (!f) return (3);

  c = (double ***) cP->p;
  if (!c) return (5);

  fprintf (f,"%s\n",name);
    
  nOnLine = 0;
  
  for (j=0;j<=ny;j++) {
    for (i=0;i<=nx;i++) {
      for (k=1;k<=6;k++) {
	nOnLine++;
	if (nOnLine > nPerLine) {
	  fprintf(f,"\n");
	  nOnLine = 1;
	}
	fprintf(f,fmt,c[i][j][k]); 
      }
    }
  }
  fprintf (f,RECORD_END_FORMAT);
  return (status);
}

int wPropArrayI(FILE *f, GENERIC_POINTER x, int nD, int **b)
{
  int status = 0;

  if (!f) return (3);

  switch (nD) {
  case 1:
    status = wPropArrayI1(f,x,b);
    break;
  case 2:
    status = wPropArrayI2(f,x,b);
    break;
  case 3:
    status = wPropArrayI3(f,x,b);
    break;
  case 4:
    status = wPropArrayI4(f,x,b);
    break;
  default:
    status = 1;
    break;
  }
  return (status);
}
int wPropArrayI1(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = INT_ON_LINE ;
  char fmt[] = INT_FORMAT ;
  int *y = NULL;
  
  y = (int    *) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (i=bnd[1][1];i<=bnd[1][2];i++) {
      nOnLine++;
      if (nOnLine > nPerLine) {
	fprintf(f,"\n");
	nOnLine = 1;
      }
      fprintf(f,fmt,y[i]);
    }
  }
  return (status);
}
int wPropArrayI2(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = INT_ON_LINE ;
  char fmt[] = INT_FORMAT ;
  int **y = NULL;
  
  y = (int   **) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) {
	nOnLine++;
	if (nOnLine > nPerLine) {
	  fprintf(f,"\n");
	  nOnLine = 1;
	}
	fprintf(f,fmt,y[i][j]);
      }
    }
  }
  return (status);
}    
int wPropArrayI3(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = INT_ON_LINE ;
  char fmt[] = INT_FORMAT ; 
  int ***y = NULL;
  
  y = (int  ***) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
		for (i=bnd[1][1];i<=bnd[1][2];i++) {
		  nOnLine++;
		  if (nOnLine > nPerLine) {
			fprintf(f,"\n");
			nOnLine = 1;
		  }
		  fprintf(f,fmt,y[i][j][k]);
		}
      }
    }
  }
  return (status);
}    
int wPropArrayI4(FILE *f,GENERIC_POINTER x, int **bnd)
{
  /* assume last index should still go fastest !!! */
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = INT_ON_LINE ;
  char fmt[] = INT_FORMAT ;
  int ****y = NULL;
  
  y = (int ****) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) 
	  for (p=bnd[4][1];p<=bnd[4][2];p++) {
	    nOnLine++;
	    if (nOnLine > nPerLine) {
	      fprintf(f,"\n");
	      nOnLine = 1;
	    }
	    fprintf(f,fmt,y[i][j][k][p]);
	  }
      }
    }
  }
  
  return (status);
}    
int wPropArrayF(FILE *f, GENERIC_POINTER x, int nD, int **b)
{
  int status = 0;

  if (!f) return (3);

  switch (nD) {
  case 1:
    status = wPropArrayF1(f,x,b);
    break;
  case 2:
    status = wPropArrayF2(f,x,b);
    break;
  case 3:
    status = wPropArrayF3(f,x,b);
    break;
  case 4:
    status = wPropArrayF4(f,x,b);
    break;
  default:
    status = 1;
    break;
  }
  return (status);
}
int wPropArrayF1(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = FLOAT_ON_LINE ;
  char fmt[] = FLOAT_FORMAT ;
  float *y = NULL;
  
  y = (float    *) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (i=bnd[1][1];i<=bnd[1][2];i++) {
      nOnLine++;
      if (nOnLine > nPerLine) {
	fprintf(f,"\n");
	nOnLine = 1;
      }
      fprintf(f,fmt,y[i]);
    }
  }
  return (status);
}
int wPropArrayF2(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = FLOAT_ON_LINE ;
  char fmt[] = FLOAT_FORMAT ;
  float **y = NULL;
  
  y = (float   **) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) {
	nOnLine++;
	if (nOnLine > nPerLine) {
	  fprintf(f,"\n");
	  nOnLine = 1;
	}
	fprintf(f,fmt,y[i][j]);
      }
    }
  }
  return (status);
}    
int wPropArrayF3(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = FLOAT_ON_LINE ;
  char fmt[] = FLOAT_FORMAT ;
  float ***y = NULL;
  
  y = (float  ***) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
		for (i=bnd[1][1];i<=bnd[1][2];i++) {
		  nOnLine++;
		  if (nOnLine > nPerLine) {
			fprintf(f,"\n");
			nOnLine = 1;
		  }
		  fprintf(f,fmt,y[i][j][k]);
		}
      }
    }
  }

  return (status);
}    
int wPropArrayF4(FILE *f,GENERIC_POINTER x, int **bnd)
{
  /* assume last index should still go fastest !!! */
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = FLOAT_ON_LINE ;
  char fmt[] = FLOAT_FORMAT ;
  float ****y = NULL;
  
  y = (float ****) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) 
	  for (p=bnd[4][1];p<=bnd[4][2];p++) {
	    nOnLine++;
	    if (nOnLine > nPerLine) {
	      fprintf(f,"\n");
	      nOnLine = 1;
	    }
	    fprintf(f,fmt,y[i][j][k][p]);
	  }
      }
    }
  }
  return (status);
}    
int wPropArrayD(FILE *f, GENERIC_POINTER x, int nD, int **b)
{
  int status = 0;

  if (!f) return (3);

  switch (nD) {
  case 1:
    status = wPropArrayD1(f,x,b);
    break;
  case 2:
    status = wPropArrayD2(f,x,b);
    break;
  case 3:
    status = wPropArrayD3(f,x,b);
    break;
  case 4:
    status = wPropArrayD4(f,x,b);
    break;
  default:
    status = 1;
    break;
  }
  return (status);
}
int wPropArrayD1(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = DOUBLE_ON_LINE ;
  char fmt[] = DOUBLE_FORMAT ;
  double *y = NULL;
  
  y = (double    *) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (i=bnd[1][1];i<=bnd[1][2];i++) {
      nOnLine++;
      if (nOnLine > nPerLine) {
	fprintf(f,"\n");
	nOnLine = 1;
      }
      fprintf(f,fmt,y[i]);
    }
  }
  return (status);
}
int wPropArrayD2(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = DOUBLE_ON_LINE ;
  char fmt[] = DOUBLE_FORMAT ;
  double **y = NULL;
  
  y = (double   **) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) {
	nOnLine++;
	if (nOnLine > nPerLine) {
	  fprintf(f,"\n");
	  nOnLine = 1;
	}
	fprintf(f,fmt,y[i][j]);
      }
    }
  }
  return (status);
}    
int wPropArrayD3(FILE *f,GENERIC_POINTER x, int **bnd)
{
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = DOUBLE_ON_LINE ;
  char fmt[] = DOUBLE_FORMAT ;
  double ***y = NULL;
  
  y = (double  ***) x;

  if (!y) status = 1;
  else if (!f) status = 3;  
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) {
	  nOnLine++;
	  if (nOnLine > nPerLine) {
	    fprintf(f,"\n");
	    nOnLine = 1;
	  }
	  fprintf(f,fmt,y[i][j][k]);
	}
      }
    }
  }
  return (status);
}    
int wPropArrayD4(FILE *f,GENERIC_POINTER x, int **bnd)
{
  /* assume last index should still go fastest !!! */
  int status = 0;
  int i, j, k, p;
  int nOnLine;

  int nPerLine = DOUBLE_ON_LINE ;
  char fmt[] = DOUBLE_FORMAT ;
  double ****y = NULL;
  
  y = (double ****) x;

  if (!y) status = 1;
  else if (!f) status = 3;
  else {
    nOnLine = 0;
    for (k=bnd[2][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) 
	  for (p=bnd[4][1];p<=bnd[4][2];p++) {
	    nOnLine++;
	    if (nOnLine > nPerLine) {
	      fprintf(f,"\n");
	      nOnLine = 1;
	    }
	    fprintf(f,fmt,y[i][j][k][p]);
	  }
      }
    }
  }
  return (status);
}    
