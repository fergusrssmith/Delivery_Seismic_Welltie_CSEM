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
#include "ECLIO/read.h"
#include "ECLIO/access.h"
#include "ECLIO/records.h"
#include "ECLIO/set.h"
#include "NR/nrutil.h"

#define PINCH_DEFAULT 0.001
#define GAP_DEFAULT "GAP"
#define PINCHPV_DEFAULT 1.0E20
#define TOPBOT_DEFAULT "TOPBOT"
#define TMETH_DEFAULT "TOP"

#define DEFAULT_L_UNITS "METRES"

#define FREE_ARG char*

/* Reads data from Eclipse-formatted keyword-oriented files and loads
   into ECL */

/*************************************************************************/
int findKey(ECL *ecl, char *keyIn, int toUpper, int rewindF)
/*************************************************************************/
{
  /* Returns 0 on successfully finding keyIn 
             1 if key not found
	     File pointer is positioned on record after keyIn */

  char key[TOKEN_LENGTH+1];
  char buffer[LINE_LENGTH+1];
  size_t i, lKey, lToken;
  char token[TOKEN_LENGTH+1];
  int foundKey =0, nw=0;
  char **words = NULL;
  int iReturn=0;
  const char routine[] = "findKey";
  
  nw = strtokSafe(keyIn, ECL_TOKENS, &words,1);
  if (nw <1) {
    iReturn=1;
    goto bottom;
  }
  
  if (!f) {
    iReturn = 3;
    goto bottom;
  }

  for (i=0;i<TOKEN_LENGTH;i++) key[i]='\0';
  strncpy(key,words[0], TOKEN_LENGTH);
  lKey = strlen(key);
  if (toUpper) for (i=0;i<lKey;i++) key[i]=toupper(key[i]);
  
  if (rewindF) rewind(f);

  foundKey = 0;
  while (fgets(buffer, (int) LINE_LENGTH, f)) { /* work thru file */
    if (strlen(buffer)>1){
      nw = strtokSafe(buffer,ECL_TOKENS, &words, 1);
      if (nw > 0 && words[0]) {	/* grab first word */
	strncpy(token,words[0],TOKEN_LENGTH); 
	lToken = strlen(token);
	if (lToken == 7 && 
	    !strncmp(key, "INCLUDE",lToken) &&
	    nw > 1) {
	  pushStack(words[1],ecl->fStack);
	  fpFromStack(ecl->fStack);
	}
	else if (lToken == lKey) {
	  if (toUpper) for (i=0;i<lToken;i++) token[i] = toupper(token[i]);
	  if (!strncmp(key, token, lKey)) {
	    foundKey = 1;
	    break;
	  }
	}
      }
    }
  }
  if (!foundKey && nInc
 bottom:
  freeWordList(words, nw);
  if (foundKey) return (0);
  return (1);
}
/*************************************************************************/
void freeFindKeyMem(char *key, char* token, char *buffer)
/*************************************************************************/
{
  if (key) free((FREE_ARG) key);
  if (token) free((FREE_ARG) token);
  if (buffer) free((FREE_ARG) buffer);
}
/*************************************************************************/
int findNextKey(ECL *ecl, FILE *f, int toUpper, int nKeys, char **keysIn)
/*************************************************************************/
{
  char **key = NULL;
  int *lKey = NULL;
  char buffer[LINE_LENGTH+1];
  size_t i, j;
  size_t lToken=0;
  char *token = NULL;
  const char routine[] = "findNextKey";
  char *cp = NULL;
  char **words = NULL;
  int nw = 0;
  int status = 0;

  if (!f) {
    status = 1;
    eclError(routine,status,"Input file pointer not set.\n");
    goto bottom;
  }

  token = (char *) malloc((size_t) (TOKEN_LENGTH + 1)*sizeof(char));
  if (!token) {
    status = 2;
    eclError(routine,status,"Could not alloc token\n");
    goto bottom;
  }

  key = (char **) malloc((size_t)(nKeys + 1) * sizeof (char *));
  lKey = ivector(1,nKeys);
  if (!key || !lKey) {
    status = 3;
    eclError(routine,status, "Could not allocate key list");
    goto bottom;
  }

  for (i=1;i<=(size_t)nKeys;i++) {

    nw = strtokSafe(buffer, ECL_WHITESPACE, &words, 1);
    if (nw>0) cp = words[0];
    else cp = (char *) NULL;
    if (cp) strncpy(token, cp, TOKEN_LENGTH);
    else {
      status = 4;
      eclError(routine,status, "Received a NULL key\n");
      goto bottom;
    }

    lKey[i] = (int) strlen(token);
    key[i] = (char *) malloc((size_t) (lKey[i] + 1)*sizeof(char));
    if (!lKey[i]) {
      status = 5;
      eclError(routine,status, "Could not allocate key");
      goto bottom;
    }

    strncpy(key[i],token,(size_t) lKey[i]+1);
    if (toUpper) for (j=0;j<(size_t)lKey[i];j++) key[i][j]=toupper(token[j]);
  }
  
  while (fscanf(f,"%s",buffer) != EOF) {
    nw = strtokSafe(buffer,ECL_WHITESPACE, &words, 1);
    if (nw > 0) cp = words[0];
    else cp = (char *) NULL;
    if (!cp) break;
    strncpy(token,cp,TOKEN_LENGTH);
    lToken = strlen(token);
    for (i=1;i<=(size_t)nKeys;i++) {
      if (lToken != (size_t) lKey[i]) break;
      if (toUpper) for (j=0;j<lToken;j++) token[j] = toupper(token[j]);
      if (!strncmp(token, key[i], lToken)) {
	freeFindNextKeyMem(token, nKeys, key, buffer, words, nw);
	return (0);
      }
    }
  }
 bottom:
  
  freeFindNextKeyMem(token, nKeys, key, buffer, words, nw);
  if (status) status = 0;
  else status = 1;
  return (status);
}
/*************************************************************************/
void freeFindNextKeyMem(char *token, int nKeys, char **key, char *buffer, 
			char **words, int nw)
/*************************************************************************/
{
  int i;

  if (token) free((FREE_ARG) token);
  if (nKeys > 0 && key) {
    for (i=1;i<=nKeys;i++) if (key[i]) free((FREE_ARG)key[i]);
    free((FREE_ARG)key);
  }
  if (buffer) free((FREE_ARG) buffer);
  if (words) freeWordList(words, nw);
}
/*************************************************************************/
int r_specGrid(FILE * f, ECL *ecl)
/*************************************************************************/
{
  int nx, ny, nz, numRes;
  char cylindrical;
  int status = 0;
  char buffer[LINE_LENGTH];
   
  (void) fgets(buffer, LINE_LENGTH, f);
  cleanRecord(buffer);
  status = sscanf(buffer, "%d %d %d %d %c", 
		  &nx, &ny, &nz, &numRes, &cylindrical);

  if (status < 3) {
    nz = 1;
    if (status < 2) {
      ny = 1;
      if (status < 1) nx = 1;
    }
  }

  if (status < 4) numRes = 1;
  if (status < 5) cylindrical = 'F';

  return (s_specGrid(ecl, nx, ny, nz, numRes, cylindrical));
}
/*************************************************************************/
int r_pinch(FILE *f, ECL *ecl)
/*************************************************************************/
{
  double pVal = D_MISSING;
  char gap[TOKEN_LENGTH+1];
  double pinchPV = D_MISSING;
  char topbot[TOKEN_LENGTH+1];
  char tMeth[TOKEN_LENGTH+1];
  int status = 0;
  char buffer[LINE_LENGTH+1];
  char sf[] ="%lf %12s %lf %12s %12s";

  (void) fgets(buffer, LINE_LENGTH, f);
  cleanRecord(buffer);

 
  status = sscanf(buffer,sf, &pVal, gap, &pinchPV, topbot, tMeth);
  if (status<1) pVal = PINCH_DEFAULT;
  if (status<2) strncpy(gap,GAP_DEFAULT, TOKEN_LENGTH);
  if (status<3) pinchPV = PINCHPV_DEFAULT;
  if (status<4) strncpy(topbot, TOPBOT_DEFAULT, TOKEN_LENGTH);
  if (status<5) strncpy(tMeth, TMETH_DEFAULT, TOKEN_LENGTH);
  return (s_pinch(ecl, pVal, gap, pinchPV, topbot, tMeth));
}
/*************************************************************************/
int r_coord(FILE *f, ECL *ecl)
/*************************************************************************/
{
  int nx, ny, nz, status=0;
  const char routine[] = "r_coord";
  PROP *pp = NULL;
  double ***c;
  int i, j, k;

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) {
    status = 1;
    eclError(routine,status, "Must read SPECGRID before COORD\n");
    goto bottom;
  }

  if (!f) {
    status = 2;
    eclError(routine,status, "Pointer to input file not defined.\n");
    goto bottom;
  }

  pp = addProp(ecl, "COORD", 'd', 3, 0, nx, 0, ny, 1, 6);
  if (!pp || !pp->p) {
    status = 3;
    eclError(routine,status, "Could not init COORD\n");
    goto bottom;
  }

  c = (double ***)pp->p;

  for (j=0;j<=ny;j++) {
    for (i=0;i<=nx;i++) {
      for (k=1;k<=6;k++) {
	if (fscanf(f,"%lf",&c[i][j][k]) == EOF) {
	  status = 4;
	  eclError(routine,status,
		   "Incomplete COORD record (i,j,k)=(%d,%d,%d)\n",
		   i, j, k);
	  goto bottom;
	}
      }
    }
  }

 bottom:
  if (!status) {
    ecl->coordp = pp;
    ecl->coord = c;
    pp->set = 1;
  }
  else {
    if (pp) deleteProp(ecl,pp->name);
    pp = NULL;
    ecl->coord = NULL;
  }

  return (status);
}
/*************************************************************************/
int r_zCorn(FILE *f, ECL *ecl)
/*************************************************************************/
{
  int nx, ny, nz;
  const char routine[] = "r_zCorn";
  PROP *pp = NULL;
  double ****z = NULL;
  int i, j, k, p, q;
  int status = 0, iRead;

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status) {
    status = 1;
    (void) eclError(routine,status, "Must read SPECGRID before ZCORN\n");
    goto bottom;
  }

  pp = addProp(ecl, "ZCORN", 'd', 4, 1, nx, 1, ny, 1, nz, 1, 8);
  if (!pp) {
    status = 2;
    (void) eclError(routine,status, "Must read SPECGRID before ZCORN\n");
    goto bottom;
  }

  if (!f) {
    status = 3;
    (void) eclError(routine,status, "Input file pointer not set\n");
    goto bottom;
  }

  z = (double ****)pp->p;
  if (!z) {
    status = 4;
    (void) eclError(routine,status, "Could not access ZCORN array data\n");
    goto bottom;
  }

  for (k=1;k<=nz;k++) {  /*z-direction scans slowest */
    for (p=0;p<=4;p+=4) { /* ...then top and bottom of layer k */
      for (j=1;j<=ny;j++) { /* then along each row */
	for (q=0;q<=2;q+=2) { /* do +- corners in j for each sweep */
	  for (i=1;i<=nx;i++) { /* write +- corners in i for each block */
	    iRead = fscanf(f,"%lf %lf",&z[i][j][k][p+q+1], &z[i][j][k][p+q+2]);
	    if (iRead < 2 || status == EOF) {
	      status = 5;
	      eclError(routine, status, "Incomplete ZCORN record\n");
	      goto bottom;
	    }}}}}}

 bottom:
  if (!status) {
    ecl->zCornp = pp;
    ecl->zCorn = z;
    pp->set = 1;
  }
  else {
    if (pp) deleteProp(ecl, pp->name);
    ecl->zCorn = NULL;
    pp = NULL;
  }
  return (status);
}
/*************************************************************************/
int r_prop(FILE *f, ECL *ecl, char *name, char ptype,
	   int nInd, ...)
/*************************************************************************/
{
  const char routine[] = "r_prop";
  PROP *pp = NULL;
  int **bnd = NULL;
  int status=0;
  va_list ap;  

  va_start(ap, nInd);
  bnd = buildBnd(ecl, nInd, ap, name);
  va_end(ap);

  if (!bnd) {
    status = 1;
    eclError(routine,status, "Could not build bounds: %s\n",name);
    goto bottom;
  }

  pp = addPropN(ecl, name, ptype, nInd, bnd);

  if (!pp) {
    status = 2;
    (void) eclError(routine,status, "Could not access property: %s\n",name);
    goto bottom;
  }

  switch (nInd) {
  case 1:
    status = r_prop1(f, name, pp->p, ptype, nInd, bnd);
    break;
  case 2:
    status = r_prop2(f, name, pp->p, ptype, nInd, bnd);
    break;
  case 3:
    status = r_prop3(f, name, pp->p, ptype, nInd, bnd);
    break;
  case 4:
    status = r_prop4(f, name, pp->p, ptype, nInd, bnd);
    break;
  default:
    status = 3;
    eclError(routine, status, 
	     "PROP %s: Must have 0<nInd<5, but nInd = %d\n",
	     name, nInd);
    goto bottom;
    break;
  }

 bottom:
  if (!status)  pp->set = 1;
  else {
    if (pp) deleteProp(ecl, pp->name);
    pp = NULL;
  }

  return (status);
}
/*************************************************************************/
int r_propGrid(FILE *f, ECL *ecl, char *name, char ptype)
/*************************************************************************/
{
  PROP *pp = NULL;
  int **bnd = NULL;
  int nInd =3;
  const char routine[] = "r_propGrid";
  int status = 0;
  
  bnd = buildBndGrid(ecl);
  if (!bnd) {
    status = 1;
    (void) eclError(routine,1, "Could not allocate bnd array for %s\n",name);
    goto bottom;
  }

  pp = addPropN(ecl, name, ptype, nInd, bnd);
  if (!pp) {
    status = 2;
    (void) eclError(routine,status, "Could not add PROP for %s\n",name);
    goto bottom;
  }

  status = r_prop3(f, name, pp->p, ptype, nInd, bnd);
  
 bottom:
  if (!status)  pp->set = 1;
  else {
    if (pp) deleteProp(ecl, pp->name);
    pp = NULL;
  }

  return (status);
}

/*************************************************************************/
int r_prop1(FILE *f, char *name, GENERIC_POINTER p, char ptype, 
	    int nInd, int **bnd)
/*************************************************************************/
{
  const char routine[] = "r_prop1";
  int *vi = NULL;
  float *vf = NULL;
  double *vd = NULL;
  int i;
  int nReq, nRead=0;

  nReq = calcNinArray(nInd, bnd);
  if (nReq < 1) {
    (void) eclError(routine,1, "No data required for PROP %s\n",name);
    return (1);
  }

  switch (ptype) {
  case 'i':
    vi = (int *) p;
    for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%d",&vi[i]);
    break;
  case 'f':
    vf = (float *) p;
    for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%f",&vf[i]);
    break;
  case 'd':
    vd = (double *) p;
    for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%lf",&vd[i]);
    break;
  default:
    eclError(routine, 2, "Illegal ptype = %c for PROP = %s\n",ptype, name);
    return (2);
    break;
  }
  return (nRead-nReq);
}
/*************************************************************************/
int r_prop2(FILE *f, char *name, GENERIC_POINTER  p, char ptype, 
	    int nInd, int **bnd)
/*************************************************************************/
{
  const char routine[] = "r_prop2";
  int **vi = NULL;
  float **vf = NULL;
  double **vd = NULL;
  int status;

  int i, j;
  int nReq, nRead=0;

  nReq = calcNinArray(nInd, bnd);
  if (nReq < 1) {
    (void) eclError(routine,1, "No data required for PROP %s\n",name);
    return (1);
  }

  switch (ptype) {
  case 'i':
    vi = (int **) p;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%d",&vi[i][j]);
    }
    break;
  case 'f':
    vf = (float **) p;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%f",&vf[i][j]);
    }
    break;
  case 'd':
    vd = (double **) p;
    for (j=bnd[2][1];j<=bnd[2][2];j++) {
      for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%lf",&vd[i][j]);
    }
    break;
  default:
    eclError(routine, 2, "Illegal ptype = %c for PROP = %s\n",ptype, name);
    return (2);
    break;
  }
  return (nRead-nReq);
} 
/*************************************************************************/ 
int r_prop3(FILE *f, char *name, GENERIC_POINTER  p, char ptype, 
	    int nInd, int **bnd)
/*************************************************************************/
{
  const char routine[] = "r_prop3";
  int ***vi = NULL;
  float ***vf = NULL;
  double ***vd = NULL;
  int status = 0;

  int i, j, k;
  int nReq, nRead=0;

  nReq = calcNinArray(nInd, bnd);
  if (nReq < 1) {
    (void) eclError(routine,1, "No data required for PROP %s\n",name);
    return (1);
  }

  switch (ptype) {
  case 'i':
    vi = (int ***) p;
    for (k=bnd[3][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%d",&vi[i][j][k]);
      }
    }
    break;
  case 'f':
    vf = (float ***) p;
    for (k=bnd[3][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) nRead+=fscanf(f, "%f",&vf[i][j][k]);
      }
    }
    break;
  case 'd':
    vd = (double ***) p;
    for (k=bnd[3][1];k<=bnd[3][2];k++) {
      for (j=bnd[2][1];j<=bnd[2][2];j++) {
	for (i=bnd[1][1];i<=bnd[1][2];i++) 
	  nRead+=fscanf(f, "%lf",&vd[i][j][k]);
      }
    }
    break;
  default:
    eclError(routine, 2, "Illegal ptype = %c for PROP = %s\n",ptype, name);
    return (2);
    break;
  }
  return (nRead-nReq);
}
/*************************************************************************/
int r_prop4(FILE *f, char *name, GENERIC_POINTER  p, char ptype, 
	    int nInd, int **bnd)
/*************************************************************************/

{
  const char routine[] = "r_prop4";
  int ****vi = NULL;
  float ****vf = NULL;
  double ****vd = NULL;

  int i, j, k, l;
  int nReq, nRead=0;

  nReq = calcNinArray(nInd, bnd);
  if (nReq < 1) {
    (void) eclError(routine,1, "No data required for PROP %s\n",name);
    return (1);
  }

  switch (ptype) {
  case 'i':
    vi = (int ****) p;
    for (l=bnd[4][1];l<=bnd[4][2];l++) {
      for (k=bnd[3][1];k<=bnd[3][2];k++) {
	for (j=bnd[2][1];j<=bnd[2][2];j++) {
	  for (i=bnd[1][1];i<=bnd[1][2];i++) 
	    nRead+=fscanf(f, "%d",&vi[i][j][k][l]);
	}
      }
    }
    break;
  case 'f':
    vf = (float ****) p;
    for (l=bnd[4][1];l<=bnd[4][2];l++) {
      for (k=bnd[3][1];k<=bnd[3][2];k++) {
	for (j=bnd[2][1];j<=bnd[2][2];j++) {
	  for (i=bnd[1][1];i<=bnd[1][2];i++) 
	    nRead+=fscanf(f, "%f",&vf[i][j][k][l]);
	}
      }
    }
    break;
  case 'd':
    vd = (double ****) p;
    for (l=bnd[4][1];l<=bnd[4][2];l++) {
      for (k=bnd[3][1];k<=bnd[3][2];k++) {
	for (j=bnd[2][1];j<=bnd[2][2];j++) {
	  for (i=bnd[1][1];i<=bnd[1][2];i++) 
	    nRead+=fscanf(f, "%lf",&vd[i][j][k][l]);
	}
      }
    }
    break;
  default:
    eclError(routine, 2, "Illegal ptype = %c for PROP = %s\n",ptype, name);
    return (2);
    break;
  }
  return (nRead-nReq);
}
/*************************************************************************/
int r_coordSys(FILE *f, ECL *ecl)
/*************************************************************************/
{

  const char routine[] = "r_coordSys";
  int *k1 = NULL, *k2 = NULL;
  char **comp = NULL; char **join = NULL;
  int numRes = I_MISSING;
  int i,status = 0, iRead;
  char buffer[LINE_LENGTH];

  if (!ecl->specGrid.set) {
    status = 1;
    eclError(routine, 1, "Must set SPECGRID before reading COORDSYS\n");
    goto bottom;
  }

  numRes = ecl->specGrid.numRes;
  if (numRes < 1) {
    status = 2;
    eclError(routine, status, "Must have numRes > 0\n");
    goto bottom;
  }
  
  k1 = ivector(1,numRes); k2 = ivector(1, numRes);
  comp = cmatrix(1, numRes, 0, TOKEN_LENGTH); 
  join = cmatrix(1, numRes, 0, TOKEN_LENGTH);
  if (!k1 || !k2 ||!comp || !join) {
    status = 3;
    (void) eclError(routine,status,"Could not allocate storage\n");
    goto bottom;
  }
  
  for (i=1;i<=numRes;i++) {
    if (fgets(buffer, LINE_LENGTH, f)) {
      scrubRecord(buffer,ECL_TRINKETS,' ');
      iRead = sscanf(buffer, "%d %d %s %s", 
		      &k1[i],&k2[i],comp[i], join[i]);
      if (iRead < 2){
	status = 4;
	eclError(routine, status, "Cannot default k1, k2, on record %d\n",i);
	goto bottom;
      }
      if (status < 3) strncpy(comp[i],"INCOMP", TOKEN_LENGTH);
      if (status < 4) strncpy(join[i], "SEPARATE", TOKEN_LENGTH);
      status = 0;
    }
    else {
      status = 5;
      eclError(routine,status,
	       "Not enough COORDSYS records (expected %d, found %d)\n",
	       numRes, i-1);
      goto bottom;
    }
  }
 bottom:
  if (status) {
    freeCoordSys(k1,k2,comp,join, numRes,TOKEN_LENGTH);
    return (status);
  }

  return(s_coordSys(ecl, numRes, k1, k2, comp, join));
}
/*************************************************************************/
int r_gridUnit(FILE *f, ECL *ecl)
/*************************************************************************/
{
  const char routine[] = "r_gridUnit";
  char units[TOKEN_LENGTH+1];
  char map[5];
  char buffer[LINE_LENGTH+1];
  int status = 0;
  char sf[]="%12s %12s";

  if (fgets(buffer, LINE_LENGTH, f)) {
    scrubRecord(buffer, ECL_TRINKETS, ' ');
    status = sscanf(buffer, sf, units, map);
    if (status < 1) strncpy(units, DEFAULT_L_UNITS, TOKEN_LENGTH);
    if (status < 2) strncpy(map, "MAP", TOKEN_LENGTH);
    return s_gridUnit(ecl, units, map);
  }
  eclError(routine,1,"Unexpected EOF\n");
  return (1);
}
/*************************************************************************/
int r_mapUnits(FILE *f, ECL *ecl)
/*************************************************************************/
{
  const char routine[] = "r_mapUnits";
  char units[TOKEN_LENGTH+1];
  char buffer[LINE_LENGTH+1];
  int status = 0;
  char sf[]="%12s";

  if (fgets(buffer, LINE_LENGTH, f)) {
    scrubRecord(buffer, ECL_TRINKETS,' ');
    status = sscanf(buffer, sf, units);
    if (status < 1) strncpy(units, DEFAULT_L_UNITS, TOKEN_LENGTH);
    return s_mapUnits(ecl, units);
  }
  eclError(routine,1,"Unexpected EOF\n");
  return (1);
}
/*************************************************************************/
int r_mapAxes(FILE *f, ECL *ecl)
/*************************************************************************/
{
  const char routine[] = "r_mapAxes";
  double **axes = NULL;
  int i, j;
  int status = 0;

  axes = dmatrix(1,3,1,2);
  if (!axes) {
    status = 1;
    (void) eclError(routine,status,"Could not allocate axes\n");
    return (1);
  }

  for (i=1;i<=3;i++) {
    for (j=1;j<=2;j++) {
      if (fscanf(f,"%lf", &axes[i][j])==EOF){
	if (axes) free_dmatrix(axes,1,3,1,2);
	eclError(routine,2, "Unexpected EOF\n");
	return (2);
      }
    }
  }
  return (s_mapAxes(ecl, axes));
}
int pushFStack(MODFILES *fStack, char *fName)
{
  FILE *f;
  
  if (!(f = fopen(fName, "r"))) {
    if (fStack->nInc < N_INC_MAX) {
      fStack->nInc++;
      fStack->fp[nInc] = f;
      strncpy(fn[nInc],fName, NAME_LENGTH);
      return (fStack->nInc);
    }
    else
      return (-1); /* out of space in stack */
  }
  else
    return (-2); /* could not open file */
}
FILE *fpFromStack(MODFILES *fStack)
{
  if (fStack->nInc) return (fStack->fp[fStack->nInc]);
}
int *clearLastFStack(MODFILES *fStack)
{
  FILE *fp;
  
  fp = fStack->fp[fStack->nInc];

  if (fp) fclose(fp);
  fp = NULL;
  fStack->nInc--;

  return (0);
}
  
  
  
  
