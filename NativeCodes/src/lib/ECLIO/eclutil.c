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
#include <stdarg.h>
#include <string.h>
#include "NR/nrutil.h"
#include "ECLIO/eclio.h"
#include "ECLIO/eclutil.h"
#include "ECLIO/access.h"

#define PROP_BUCKET 10
#define INC_BUCKET 10

/* Routines to manipulate  ECL components
   CD White, LSU, 5/04
*/

/*************************************************************************/
int initModel(ECL *e)
/*************************************************************************/
{
  COORDSYS *cs;
  MAPAXES *ma;

  e->specGrid.set = 0;
  e->pinch.set = 0;
  
  cs = &(e->coordSys);
  cs->set = 0;
  cs->k1 = NULL;
  cs->k2 = NULL;
  cs->comp = NULL;
  cs->join = NULL;
  
  e->mapUnits.set = 0;
  e->gridUnit.set = 0;
  
  ma = &(e->mapAxes);
  ma->set  = 0;
  ma->axes = NULL;

  e->zCornp = ( PROP *) NULL;
  e->coordp = ( PROP *) NULL;
  e->actnump = ( PROP *) NULL;
  e->zCorn = (double ****) NULL;
  e->coord = (double ***) NULL;
  e->actnum = (int ***) NULL;

  e->prop = NULL;
  if (!initPropVec(e)) return (1);

  return (0);
}
/*************************************************************************/
 PROP *initPropVec(ECL *ecl)
/*************************************************************************/
{
  PROP *pv = NULL;
  PROP *p = NULL;
  int i;
  /* This piece was causing problems in Fedora; 
     for now will assume that pv is empty going in ...
     will cause a memory leak if this model was not freed properly 
     before calling initModel */
  /*
  if (pv && ecl->nPropAlloc > 1) {
    for (i=1;i<=ecl->nProp;i++) if (pv[i]) free_anyProp(&pv[i]);
    free((FREE_ARG) pv);
  }
  */
  
  pv = ( PROP *) malloc((size_t)(PROP_BUCKET +1)*sizeof( PROP));
  if (!pv) {
    free((FREE_ARG) pv);
    eclError("initPropVec",1,
	     "Could not do initial alloc on PROP vector\n");
    pv = NULL;
  }
  else {
    ecl->prop = pv;
    ecl->nPropAlloc = PROP_BUCKET;
    ecl->nProp = 0;
    for (i=1, p = &pv[1];i<=ecl->nPropAlloc;i++,p++) initProp(p);
  }
  return (pv);
}
/*************************************************************************/
void initProp(PROP *p)
/*************************************************************************/
{
  setProp(p, "NULL", '0', 0, (int **) NULL, (GENERIC_POINTER ) NULL); 
  p->set = 0;
}
/*************************************************************************/
 PROP *addPropN(ECL *ecl, char *name, char ptype, 
		     int nInd, int **bnd)
/*************************************************************************/
{
  PROP *pp = NULL;
  GENERIC_POINTER p = NULL;


  /* set a pointer to the first empty PROP, reallocating if necessary */
  pp = allocProp(ecl, name);

  /* get storage */
  p = anyProp(ptype, nInd, bnd, name);
  /* and intialize property */
  setProp(pp, name, ptype, nInd, bnd, p);

  return(pp);
}
/*************************************************************************/
 PROP *addProp0(ECL *ecl, char *name, char ptype, 
		     int nInd, int **bnd)
/*************************************************************************/
{
  PROP *pp = NULL;
  GENERIC_POINTER p = NULL;

  /* set a pointer to the first empty PROP, reallocating if necessary */
  pp = allocProp(ecl, name);
  if (!pp) return ((GENERIC_POINTER) NULL);

  /* don't get storage */
  p = (GENERIC_POINTER) NULL;
  /* and intialize property */
  setProp(pp, name, ptype, nInd, bnd, p);

  return(pp);
}
/*************************************************************************/
PROP *allocProp(ECL *ecl, char *name)
/*************************************************************************/
{
  PROP *pv=NULL;
  PROP *pp=NULL;
  const char routine[] = "reallocProp";

  /* If this property already exists, fetch its pointer and empty the
     structure */
  pp = fetchProp(ecl, name);
  if (pp) {
	free_anyProp(pp);
    return (pp);
  }
  /* pass 2: if there is an empty prop, use it */
  else if ((pp=fetchProp(ecl,"NULL"))) 
	{
	  return (pp);
	}
  /* Need to add a property */
  else {
    pv = ecl->prop;
    ecl->nProp++;
    /* make the bucket bigger if needed*/
    if (ecl->nProp+1 >= ecl->nPropAlloc) {
      ecl->nPropAlloc += PROP_BUCKET;
	  pv = ( PROP *)
		realloc(pv, (size_t) (ecl->nPropAlloc * sizeof( PROP)));
      if (!pv) {
		eclError(routine,1, "Could not realloc PROP for %s\n",name);  
		return (NULL);
      }
      ecl->prop = pv;
	}
    /* return pointer to the requested property */
    return (&pv[ecl->nProp]);
  }
}
/*************************************************************************/
void setProp(PROP *pp, char *name, char ptype, 
	     int nInd, int **bnd, GENERIC_POINTER p)
/*************************************************************************/
{

  strncpy(pp->name, name, TOKEN_LENGTH);
  pp->ptype = ptype;
  pp->nInd = nInd;
  pp->p = p;
  pp->set = 1;
  pp->bnd = bnd;
}
/*************************************************************************/
int **buildBnd(ECL *ecl, int nInd, va_list ap, char *name)
/*************************************************************************/
{
  int **bnd = NULL;
  int i, j;
  const char routine[] = "buildBnd";
  int status = 0;

  if (nInd < 1 || nInd >4) {
    status = 1;
    eclError(routine, status, 
	     "Number of dimensions not [1..4]\nnInd = %d\nName = %s\n",
	     nInd, name);
    goto bottom;
  }

  bnd = imatrix(1,nInd,1,2);
  if (!bnd) {
    status = 2;
    eclError(routine,status, "Could not allocate pointer to bnd matrix\n");
    goto bottom;
  }

  for (i=1;i<=nInd;i++) {
    for (j=1;j<=2;j++) bnd[i][j] = va_arg(ap, int);
  }
  
  for (i=1;i<=nInd;i++) {
    if (bnd[i][2]-bnd[i][1] < 0) {
      status = 3;
      eclError(routine, status,
	       "Nonsense dimension, direction %d\nName = %s\nLow: %d Hi:%d\n",
	       i, name, bnd[i][1], bnd[i][2]);
      goto bottom;
    }
  }
  
 bottom:
  if (status) {
    if (bnd) free_imatrix(bnd, 1, nInd, 1, 2);
    bnd = NULL;
  }

  return (bnd); 
}
int **copyBnd (int n, int ** b)
{
  int i,j;
  int **new = NULL;
  
  if (!b || n < 1 || n > 4) 
    new = NULL;
  else {
    new = imatrix(1,n,1,2);
    if (!new) new = NULL;
    else
      for (i=1;i<=n;i++) {
	for (j=1;j<=2;j++) new[i][j] = b[i][j];
      }
  }
  return (new);
}
/*************************************************************************/
int **buildBndGrid(ECL *ecl)
/*************************************************************************/
{
  int **bnd = NULL;
  const char routine[] = "buildBndGrid";
  int nx, ny, nz, status=0;
  
  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  if (status){
    status = 1;
    eclError(routine,status, "Must set SPECGRID before calling.");
    goto bottom;
  }
  
  bnd = imatrix(1,3,1,2);
  if (!bnd) {
    status = 2;
    eclError(routine,status, "Could not allocate pointer to bnd matrix\n");
    goto bottom;
  }

  bnd[1][1] = 1;
  bnd[1][2] = nx;
  bnd[2][1] = 1;
  bnd[2][2] = ny;
  bnd[3][1] = 1;
  bnd[3][2] = nz;

 bottom:
  if (status) {
    if (bnd) free_imatrix(bnd,1,3,1,2);
    bnd = NULL;
  }

  return (bnd);
}
/*************************************************************************/
GENERIC_POINTER anyProp(char ptype, int nInd, int **bnd, char *name)
/*************************************************************************/
     /* space for PROP (but not its array) must already be allocated, 
	This routine allocates memory for the property "name" of type "ptype" 
	with rank "nInd". The upper and lower bounds are given 
	pairwise in bnds 
	name is passed only for debugging messages */
{
  GENERIC_POINTER p = NULL;
  const char routine[] = "anyProp";
  
  switch (nInd) {
  case 1:
    p = any1Prop(ptype, nInd, bnd, name);
    break;
  case 2:
    p = any2Prop(ptype, nInd, bnd, name);
    break;
  case 3:
    p = any3Prop(ptype, nInd, bnd, name);
    break;
  case 4:
    p = any4Prop(ptype, nInd, bnd, name);
    break;
  default:
    eclError(routine, 1, 
	     "Illegal value of rank: %d\nProperty: %s\n",
	     nInd, name);
    p = NULL;
    goto bottom;
    break;
  } 

  if (!p) {
    eclError(routine, 2, 
	     "Returned null pointer\nProperty: %s\nnInd: %d\nptype: %c\n",
	     name,nInd, ptype);
    goto bottom;
  }
  
 bottom:
  return (p);
}
/*************************************************************************/
GENERIC_POINTER any1Prop(char ptype, int nInd, int **bnd, char *name)
/*************************************************************************/
{
  const char routine[] = "any1Prop";
  GENERIC_POINTER p = NULL;

  switch (ptype) {
  case 'i':
    p = (GENERIC_POINTER ) ivector(bnd[1][1], bnd[1][2]);
    break;
  case 'f':
    p = (GENERIC_POINTER ) vector(bnd[1][1], bnd[1][2]);
    break;
  case 'd':
    p = (GENERIC_POINTER ) dvector(bnd[1][1], bnd[1][2]);
    break;
  default:
    eclError(routine, 1, 
	     "Called with illegal type = %c\nProperty %s\n",
	     ptype, name);
    return ((GENERIC_POINTER ) NULL);
    break;
  }
  return (p);
}

/*************************************************************************/
GENERIC_POINTER any2Prop(char ptype, int nInd, int **bnd, char *name)
/*************************************************************************/
{
  const char routine[] = "any2Prop";
  GENERIC_POINTER p = NULL;

  switch (ptype) {
  case 'i':
    p = (GENERIC_POINTER) imatrix(bnd[1][1], bnd[1][2],
				   bnd[2][1], bnd[2][2]);
    break;
  case 'f':
    p = (GENERIC_POINTER) matrix(bnd[1][1], bnd[1][2],
				  bnd[2][1], bnd[2][2]);
    break;
  case 'd':
    p = (GENERIC_POINTER) dmatrix(bnd[1][1], bnd[1][2],
				   bnd[2][1], bnd[2][2]);
    break;
  default:
    eclError(routine, 1, 
	     "Called with illegal type = %c\nProperty %s\n",
	     ptype, name);
    return ((GENERIC_POINTER ) NULL);
    break;
  }
  return (p);
}
/*************************************************************************/
GENERIC_POINTER any3Prop(char ptype, int nInd, int **bnd, char *name)
/*************************************************************************/
{
  const char routine[] = "any3Prop";
  GENERIC_POINTER p = NULL;

  switch (ptype) {
  case 'i':
    p = (GENERIC_POINTER ) i3tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2]);
    break;
  case 'f':
    p = (GENERIC_POINTER ) f3tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2]);
    break;
  case 'd':
    p = (GENERIC_POINTER ) d3tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2]);
    break;
  default:
    eclError(routine, 1, 
	     "Called with illegal type = %c\nProperty %s\n",
	     ptype, name);
    return ((GENERIC_POINTER ) NULL);
    break;
  }
  return (p);
}
/*************************************************************************/
GENERIC_POINTER any4Prop(char ptype, int nInd, int **bnd, char *name)
/*************************************************************************/
{
  const char routine[] = "any4Prop";
  GENERIC_POINTER p = NULL;

  switch (ptype) {
  case 'i':
    p = (GENERIC_POINTER ) i4tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2],
				    bnd[4][1], bnd[4][2]);
    break;
  case 'f':
    p = (GENERIC_POINTER ) f4tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2],
				    bnd[4][1], bnd[4][2]);
    break;
  case 'd':
    p = (GENERIC_POINTER ) d4tensor(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2],
				    bnd[3][1], bnd[3][2],
				    bnd[4][1], bnd[4][2]);
    break;
  default:
    eclError(routine, 1, 
	     "Called with illegal type = %c\nProperty %s\n",
	     ptype, name);
    return ((GENERIC_POINTER ) NULL);
    break;
  }
  return (p);
}
/*************************************************************************/
int eclError(const char *caller, int code, char *fmt, ...)
/*************************************************************************/
     /* Error handling routine. Prints an error message
	specified by fmt, inserting an arbitrary number of 
	values.
	if doExit, prints an ERROR andexits with status 1. 
	Else prints a WARNING and returns 0.
     */
{
  va_list ap;
 
  fprintf (stderr,"eclWarning: Issuing WARNING\n");
  fprintf(stderr, "\tCalled by:\t%s\n\tWarning Code:\t%d\n", caller, code);
  
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
 
  fprintf(stderr, "\n");

  return (0);
}
/*************************************************************************/
int reportError(int doExit, const char *caller, int code, char *fmt, ...)
/*************************************************************************/
     /* Error handling routine. Prints an error message
	specified by fmt, inserting an arbitrary number of 
	values.
	if doExit, prints an ERROR andexits with status 1. 
	Else prints a WARNING and returns 0.
     */
{
  va_list ap;
 
  if (doExit) {
    fprintf (stderr, "reportError: Exiting due to FATAL ERROR\n");
    fprintf(stderr, "\tCalled by:\t%s\n\tError Code:\t%d\n", caller, code);
    
  }
  else {
    fprintf (stderr,"reportWarning: Issuing WARNING\n");
    fprintf(stderr, "\tCalled by:\t%s\n\tWarning Code:\t%d\n", caller, code);
  }

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
 
  fprintf(stderr, "\n");

  if (doExit) exit(EXIT_FAILURE);
  else return (0);
}
/*************************************************************************/
void eclExit(ECL *ecl, int status)
/*************************************************************************/
{
  if (ecl)  freeModel(ecl);
  else status = 1;

  if (status) exit(EXIT_FAILURE);
  else exit(EXIT_SUCCESS);
}
/*************************************************************************/
void free_anyProp(PROP *p)
/*************************************************************************/
{
  switch (p->nInd) {
  case 1:
    free_any1Prop(p->nInd, p->bnd, p->p, p->ptype);
    break;
  case 2:
    free_any2Prop(p->nInd, p->bnd, p->p, p->ptype);
    break;
  case 3:
    free_any3Prop(p->nInd, p->bnd, p->p, p->ptype);
    break;
  case 4:
    free_any4Prop(p->nInd, p->bnd, p->p, p->ptype);
    break;
  default:
    break;
  } 
  if (p->bnd) free_imatrix(p->bnd, 1, p->nInd, 1, 2);
  initProp(p);

}
/*************************************************************************/
void free_any1Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype)
/*************************************************************************/
{
  int *iv = NULL;
  float *fv = NULL;
  double *dv = NULL;

  if (ptype){
    switch (ptype) {
    case 'i':
      iv = (int *) p;
      if (iv) free_ivector(iv, bnd[1][1], bnd[1][2]);
      break;
    case 'f':
      fv = (float *) p;
      if (fv) free_vector(fv, bnd[1][1], bnd[1][2]);
      break;
    case 'd':
      dv = (double *) p;
      if (dv) free_dvector(dv, bnd[1][1], bnd[1][2]);
      break;
    default:
      break;
    }
  }
}
/*************************************************************************/
void free_any2Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype)
/*************************************************************************/
{
  int **iv=NULL;
  float **fv=NULL;
  double **dv=NULL;

  switch (ptype) {
  case 'i':
    iv = (int **) p;
    if (iv) free_imatrix(iv,
			 bnd[1][1], bnd[1][2],
			 bnd[2][1], bnd[2][2]);
    break;
  case 'f':
    fv = (float **) p;
    if (fv) free_matrix(fv,
			bnd[1][1], bnd[1][2],
			bnd[2][1], bnd[2][2]);
    break;
  case 'd':
    dv = (double **) p;
    if (dv) free_dmatrix(dv,
			 bnd[1][1], bnd[1][2],
			 bnd[2][1], bnd[2][2]);
    break;
  default:
    break;
  }
}
/*************************************************************************/
void free_any3Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype)
/*************************************************************************/
{
  int ***iv=NULL;
  float ***fv=NULL;
  double ***dv=NULL;

  switch (ptype) {
  case 'i':
    iv = (int ***) p;
    if (iv) free_i3tensor(iv,
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2]);
    break;
  case 'f':
    fv = (float ***) p;
    if (fv) free_f3tensor(fv,
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2]);
    break;
  case 'd':
    dv = (double ***) p;
    if (dv) free_d3tensor(dv, 
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2]);
    break;
  default:
    break;
  }
}
/*************************************************************************/
void free_any4Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype)
/*************************************************************************/
{
  int ****iv=NULL;
  float ****fv=NULL;
  double ****dv=NULL;

  switch (ptype) {
  case 'i':
    iv = (int ****) p;
    if (iv) free_i4tensor(iv, 
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2],
			  bnd[4][1], bnd[4][2]);
    break;
  case 'f':
    fv = (float ****) p;
    if (fv) free_f4tensor(fv, 
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2],
			  bnd[4][1], bnd[4][2]);
    break;
  case 'd':
    dv = (double ****) p;
    if (dv) free_d4tensor(dv, 
			  bnd[1][1], bnd[1][2],
			  bnd[2][1], bnd[2][2],
			  bnd[3][1], bnd[3][2],
			  bnd[4][1], bnd[4][2]);
    break;
  default:
    break;
  }
}
/*************************************************************************/
int transProp (PROP *pp, double factor, double transl)
/*************************************************************************/
/* Multiplies any ECL PROP by factor and adds transl */
{
  char t = '\0';
  const char routine[] = "transProp";
  int i, n=I_MISSING;
  int **bnd = NULL;
  int status = 0;
  
  if (!pp) {
    status = 1;
    (void) eclError(routine, status, "Pointer to PROP %s looks to be null\n",
	     pp->name);
    goto bottom;
  }

  if (!pp->set) { 
    status = 2;
    (void) eclError(routine, status, "Must SET property before call: %s\n", pp->name);
    goto bottom;
  }
  
  n=pp->nInd;
  if (n<1 || n>4) {
    status = 3;
    (void) eclError(routine,3, "Illegal # of dimensions for %s: nInd = %d\n",
	     pp->name, n);
    goto bottom;
  }

  bnd = pp->bnd;
  for (i=1;i<=n;i++) {
    if (bnd[i][1] > bnd[i][2]) {
      status = 4;
      (void) eclError(routine, status, 
		      "Property %s has bad bound in direction %d: %d to %d\n",
		      pp->name, i, bnd[i][1], bnd[i][2]);
      goto bottom;
    }
  }

  if (!(pp->p)) {
    status = 5;
    (void) eclError(routine, status, 
	     "Property %s has a null pointer to data\n", pp->name);
    goto bottom;
  }
  
  t = pp->ptype;
  
  switch (t) {
  case 'i':
    return ((int) transIProp(pp, factor, transl));
    break;
  case 'f':
    return ((int) transFProp(pp, factor, transl));
    break;
  case 'd':
    return ((int) transDProp(pp, factor, transl));
    break;
  default:
    status = 6;
    eclError(routine, status, "Cannot handle unknown type: %c\n for PROP %s\n", 
	     t, pp->name);
    break;
  }
  
 bottom:
  return (status);
}

/*************************************************************************/
int transIProp( PROP *pp, double f, double a)
/*************************************************************************/
{
  int n= I_MISSING, **bnd = NULL;
  GENERIC_POINTER p = NULL;
  
  n=pp->nInd;
  bnd = pp->bnd;
  p = pp->p;

  switch (n) {
  case 1:
    return (transI1Prop(n, bnd, p, f, a));
    break;
  case 2:
    return (transI2Prop(n, bnd, p, f, a));
    break;
  case 3:
    return (transI3Prop(n, bnd, p, f, a));
    break;
  case 4:
    return (transI4Prop(n, bnd, p, f, a));
    break;
  }
  return ((int) 0);
}
int transFProp( PROP *pp, double f, double a)
{
  int n=I_MISSING, **bnd=NULL;
  GENERIC_POINTER p=NULL;

  n=pp->nInd;
  bnd = pp->bnd;
  p = pp->p;  

  switch (n) {
  case 1:
    return (transF1Prop(n, bnd, p, f, a));
    break;
  case 2:
    return (transF2Prop(n, bnd, p, f, a));
    break;
  case 3:
    return (transF3Prop(n, bnd, p, f, a));
    break;
  case 4:
    return (transF4Prop(n, bnd, p, f, a));
    break;
  }
  return ((int) 0);
} 
/*************************************************************************/
int transDProp( PROP *pp, double f, double a)
/*************************************************************************/
{
  int n=I_MISSING, **bnd=NULL;
  GENERIC_POINTER p=NULL;

  n=pp->nInd;
  bnd = pp->bnd;
  p = pp->p;  
  
  switch (n) {
  case 1:
    return (transD1Prop(n, bnd, p, f, a));
    break;
  case 2:
    return (transD2Prop(n, bnd, p, f, a));
    break;
  case 3:
    return (transD3Prop(n, bnd, p, f, a));
    break;
  case 4:
    return (transD4Prop(n, bnd, p, f, a));
    break;
  }
  return ((int) 0);
} 
/*************************************************************************/
int transI1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int *v = NULL;
  int fi=I_MISSING, ai=I_MISSING;
  int i;

  fi = (int) f; ai = (int) a;

  v = (int *) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];

  for (i=nrl;i<=nrh;i++) {
    v[i]=fi*v[i]+ai;
  }
  return ((int) 0);
}
/*************************************************************************/
int transI2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;
  int **v = NULL;
  int fi=I_MISSING, ai=I_MISSING;
  int i,j;

  fi = (int) f; ai = (int) a;

  v = (int **) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      v[i][j]=fi*v[i][j]+ai;
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transI3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  int ***v = NULL;
  int fi=I_MISSING, ai=I_MISSING;
  int i,j,k;

  fi = (int) f; ai = (int) a;

  v = (int ***) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	v[i][j][k]=fi*v[i][j][k]+ai;
      }
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transI4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  int npl=I_MISSING, nph=I_MISSING;
  int ****v = NULL;
  int fi=I_MISSING, ai=I_MISSING;
  int i,j,k, m;

  fi = (int) f; ai = (int) a;

  v = (int ****) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];  
  npl = bnd[4][1]; nph = bnd[4][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	for (m=npl;m<=nph;m++) {
	  v[i][j][k][m]=fi*v[i][j][k][m]+ai;
	}
      }
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transF1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  float *v = NULL;
  int i;
  float ff, af;

  ff = (float) f;
  af = (float) a;

  v = (float *) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];

  for (i=nrl;i<=nrh;i++) {
    v[i]=ff*v[i]+af;
  }
  return ((int) 0);
}
/*************************************************************************/
int transF2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;
  float **v = NULL;
  float ff, af;
  int i,j;

  ff = (float) f; af = (float) a;

  v = (float **) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      v[i][j]=ff*v[i][j]+af;
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transF3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  float ***v = NULL;
  float ff, af;
  int i,j,k;

  ff = (float) f; af = (float) a;

  v = (float ***) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	v[i][j][k]=ff*v[i][j][k]+af;
      }
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transF4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  int npl=I_MISSING, nph=I_MISSING;
  float ****v = NULL;
  float ff, af;
  int i,j,k, m;

  ff = (float) f; af = (float) a;

  v = (float ****) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];  
  npl = bnd[4][1]; nph = bnd[4][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	for (m=npl;m<=nph;m++) {
	  v[i][j][k][m]=ff*v[i][j][k][m]+af;
	}
      }
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transD1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  double *v = NULL;
  int i;
 
 v = (double *) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];

  for (i=nrl;i<=nrh;i++) {
    v[i]=f*v[i]+a;
  }
  return (0);
}
/*************************************************************************/
int transD2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;
  double **v = NULL;
  int i,j;

  v = (double **) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      v[i][j]=f*v[i][j]+a;
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transD3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  double ***v = NULL;
  int i,j,k;

  v = (double ***) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	v[i][j][k]=f*v[i][j][k]+a;
      }
    }
  }
  return ((int) 0);
}  
/*************************************************************************/
int transD4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a)
/*************************************************************************/
{
  int nrl=I_MISSING, nrh=I_MISSING;
  int ncl=I_MISSING, nch=I_MISSING;  
  int ndl=I_MISSING, ndh=I_MISSING;
  int npl=I_MISSING, nph=I_MISSING;
  double ****v = NULL;
  int i,j,k, m;

  v = (double ****) p;
  nrl = bnd[1][1]; nrh = bnd[1][2];
  ncl = bnd[2][1]; nch = bnd[2][2];
  ndl = bnd[3][1]; ndh = bnd[3][2];  
  npl = bnd[4][1]; nph = bnd[4][2];

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	for (m=npl;m<=nph;m++) {
	  v[i][j][k][m]=f*v[i][j][k][m]+a;
	}
      }
    }
  }
  return ((int)0);
}  
/*************************************************************************/
double ddot_ECL(int nrl1, int nrh1, double *v1, 
	    int nrl2, int nrh2, double *v2)
/*************************************************************************/
     /* Dot product of two vectors */
{
  int i;
  double d;
  int dr1, dr2;

  if (!v1 || !v2) 
    reportError(1,"ddot_ECL",1, "Null pointers\n");
  if ((nrh1-nrl1) != (nrh2-nrl2)) 
    reportError(1, "ddot_ECL",2, "Uncomformable\n");
  
  dr1 = nrl1-1; dr2 = nrl2-1;
  for (d=0, i=nrl1;i<=nrh1;i++) d += v1[i+dr1]*v2[i+dr2];

  return (d);
}
/*************************************************************************/
double *dgetCol(int getCol, int nrl, int nrh, int ncl, int nch, double **m)
/*************************************************************************/
{
  int i;
  double *d = NULL;
  if (getCol < ncl || getCol > nch)
    reportError(1, "dgetCol", 1, "Tried to fetch an out of range column\n");

  d = dvector(nrl, nrh);
  for (i=nrl;i<=nrh;i++) d[i] = m[i][getCol];
  return (d);
}
/*************************************************************************/
double *dgetRow(int getRow, int nrl, int nrh, int ncl, int nch, double **m)
/*************************************************************************/
{
  int i;
  double *d=NULL;
  if (getRow < nrl || getRow > nrh) 
    reportError(1, "dgetRow",1, "Tried to fecth an out of range row\n");
  
  d = dvector (ncl, nch);
  for (i=ncl;i<=nch; i++) d[i]=m[getRow][i];
  return (d);
}
  
/*************************************************************************/
double **dmmMult(int nrl1, int nrh1, int ncl1, int nch1, double **m1,
		 int nrl2, int nrh2, int ncl2, int nch2, double **m2)
/*************************************************************************/
{
  /* returns an NR-formatted double resulting from a matrix multiply,
     or NULL if not conformable 
     It is always returned with ranges from m2[1..nr1][1..nc2] and 
     should be shifted by the caller if desired.
  */
  int nr1, nc1, nc2, nr3, nc3;
  int dr1, dc1, dr2, dc2;
  double **d = NULL;
  int i3, j3;
  int p;

  const char routine[] = "dmmMult";


  nr1 = nrh1 - nrl1 + 1;
  nc2 = nch2 - ncl2 + 1;
  
  if (nr1 != nc2) 
    reportError(1,routine, 1, "Uncomformable matrices, nr1 = %d, nc2 = %d/n",
	     nr1, nc1);

  nr3 = nr1;
  nc3 = nc2;
  d = dmatrix(1, nr3, 1, nc3);

  dr1 = nrl1 - 1;
  dc1 = ncl1 - 1;
  dr2 = nrl2 - 1;
  dc2 = ncl2 - 1;

  
  for (i3=1; i3<=nr3;i3++) {
    for (j3=1; j3<=nc3; j3++) {
      for (d[i3][j3]=0.,p=1;p<=nc1;p++) {
	d[i3][j3]+=m1[i3+dr1][p+dc1]*m2[p+dr2][j3+dc2];
      }
    }
  }
  return (d);
}
  
/*************************************************************************/
int dmvMult(int nrl, int nrh, int ncl, int nch, double **m,
	    int nl, int nh, double *v, double *d)
/*************************************************************************/
  /* Computes a (double) matrix vector product

     NOTE CHANGES:
     1. Now returns a vector with same bounds as the input vector,
        rather than the same bounds as the rows of the input matrix --
	seems more natural 
     2. Now accepts a pointer to an *externally* allocated vector d
        for the results rather than internally allocating. Should be 
	allocated NR-like, using dvector.
     3. Now an int; returns zero for success
     4. We will NOT exit on error but will return a code of 1 to allow 
        calling routin to clean up
*/
{
  int nr = nrh-nrl+1, nc = nch-ncl+1, nv = nh-nl+1;
  int dr = nrl-1, dc = ncl-1, dv = nl-1;
  int i, j;

  if (nc != nv) return (1);

  for (i=1; i<=nr;i++) {
    for (d[i+dv]=0.,j=1;j<=nv;j++) {
      d[i+dv]+=m[i+dr][j+dc]*v[j+dv];
    }
  }
  return (0);
}
/*************************************************************************/
double det3mat(int nrl, int nrh, int ncl, int nch, double **m)
/*************************************************************************/
{
  int dr = nrl-1;
  int dc = ncl-1;  
  double d;

  if (nrh -nrl + 1 != 3 || nrh - nrl != nch - ncl)
    (void) reportError(1, "det3mat",1, "Valid only for 3X3 matrices\n");
  
  if (!m)
    (void) reportError(1, "det3mat",2, "Could not access matrix\n");
  
  d = 
    m[1+dr][1+dc]*m[2+dr][2+dc]*m[3+dr][3+dc] +
    m[1+dr][2+dc]*m[2+dr][3+dc]*m[3+dr][1+dc] +
    m[1+dr][3+dc]*m[2+dr][1+dc]*m[3+dr][2+dc] -
    m[1+dr][3+dc]*m[2+dr][2+dc]*m[3+dr][1+dc] -
    m[1+dr][2+dc]*m[2+dr][1+dc]*m[3+dr][3+dc] -
    m[1+dr][1+dc]*m[2+dr][3+dc]*m[3+dr][2+dc] ;

  return (d);
}
/*************************************************************************/
int transD2(int nrl, int nrh, int ncl, int nch, double **v,
	       double f, double a)
/*************************************************************************/
{
  int i,j;
  
  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      v[i][j]=f*v[i][j]+a;
    }
  }
  return (0);
}
/*************************************************************************/
int fetchPropInd( ECL *ecl, char *name)
/*************************************************************************/
{
  int i;
  PROP *p = NULL;
  char *pName = NULL;
  const char routine[] = "fetchPropInd";

  for (i=1;i<=ecl->nProp;i++) {
    p = &(ecl->prop[i]);
    pName = p->name;
    if (!strcmp(pName, name)) return (i);
  }
  
  eclError(routine,1,"Could not find property %s\n",name);
  return (0);
}
/*************************************************************************/
int compressProp( ECL *ecl)
	 /*************************************************************************/
{
  if (ecl->nProp <= ecl->nPropAlloc +1) return (0);
  ecl->prop = ( PROP *) 
    realloc(ecl->prop, (size_t)(ecl->nProp + 1)*sizeof( PROP));
  if (!ecl->prop) return (1);
  ecl->nPropAlloc = ecl->nProp + 1; 
  return (setAuxPointers(ecl));
}
/*************************************************************************/
int setAuxPointers( ECL *ecl)
/*************************************************************************/
{
  PROP *pp = NULL;
  
  pp = fetchProp(ecl, "ZCORN");
  if (!pp) return (1);
  ecl->zCornp = pp;
  ecl->zCorn = (double ****) pp->p;

  pp = fetchProp(ecl, "COORD");
  if (!pp) return (1);
  ecl->coordp = pp;
  ecl->coord = (double ***) pp->p;

  pp = fetchProp(ecl, "ACTNUM");
  if (!pp) return (1);
  ecl->actnump = pp;
  ecl->actnum = (int ***) pp->p;

  return (0);
}
/*************************************************************************/
int checkSpecGrid( ECL *ecl, int *nxp, int *nyp, int *nzp)
/*************************************************************************/
{
  SPECGRID *sg = NULL; 
  
  sg = &ecl->specGrid;
  if (!sg) return (1);
  if (!sg->set) return (2);
  if ((*nxp = sg->nx ) < 1) return (-1);
  if ((*nyp = sg->ny ) < 1) return (-2);
  if ((*nzp = sg->nz ) < 1) return (-3);

  return (0);
}
/*************************************************************************/
int calcNinArray(int nInd, int **bnd)
/*************************************************************************/
{
  int n=0, i;
  
  if (nInd<1) return (0);
  
  for (n=1,i=1;i<=nInd;i++) n=n*(bnd[i][2]-bnd[i][1]+1);

  return (n);
}
/*************************************************************************/
void freeCoordSys(int *k1, int *k2, char **comp, char **join, int n, int l)
/*************************************************************************/
{
  if (k1) free_ivector(k1,1,n);
  if (k2) free_ivector(k2,1,n);
  if (comp) free_cmatrix(comp,1,n,0,l);
  if (join) free_cmatrix(join,1,n,0,l);
}
double trilinearInterp(double *p, double *u)
{
  /* Given value of p at 8 vertices, interpolates to u
     u is a three-vector with each dimension [0,1] */

  double pu, x, y, z, xc, yc, zc;
  x = u[1]; xc = 1.-x;
  y = u[2]; yc = 1.-y;
  z = u[3]; zc = 1.-z;
  
  pu = 
    p[1]*xc*yc*zc +
    p[2]*x *yc*zc +
    p[3]*xc*y *zc +
    p[4]*x *y *zc +
    p[5]*xc*yc*z  +
    p[6]*x *yc*z  +
    p[7]*xc*y *z  +
    p[8]*x *y *z ;
  return (pu);
}
double trilinearNull(double *p, double *u)
{
  return (*p);
}



GENERIC_POINTER manageHeap(GENERIC_POINTER heapP, 
			   size_t need, size_t *haveP, size_t sizeEach, size_t chunk,
			   int condense)
     /* returns a pointer to a heap of sufficient size.
	need = the required allocation
        haveP = pointer to allocated amount
	unless allocation error then -1
        chunk = number of elements to add to the heap if required
	if passed as <= 0, reset to a default of chunkDef
	condense = 
	0 to grow or allocate (realloc to need or *haveP + chunk, need >= *haveP)
	1 to condense (realloc to need, need < *haveP)
     */
{
  const int chunkDef = 50;

  if (chunk < 0) chunk = chunkDef;
  
  if (heapP) {
    if (condense && *haveP != need) {
      *haveP = need;
    }
    else if (need > *haveP) { /* grow the heap */
      if (*haveP == 0) {
		free(heapP);   
		*haveP = need;
		if (need < chunk) *haveP = chunk;
		heapP = (GENERIC_POINTER) malloc(*haveP * sizeEach);
      }
      else {
		*haveP += chunk;
		GENERIC_POINTER newHeap;		  
		newHeap = (GENERIC_POINTER) realloc(heapP, *haveP * sizeEach);
		heapP = newHeap;
      }
    } 
    /* if the pointer is nonNull and need is null, we must want to free */
    /* You better have freed contents outside of this routine!!!! */
    else if (need == 0) { 
	  *haveP = 0;    
      free(heapP);
      heapP = (GENERIC_POINTER) NULL;
    }
    else heapP = heapP;
  }
  
  else { /* heap is null  -- can only allocate sensibly */
    if (condense) { /*can't condense a null heap */
      *haveP = -1;
      heapP = (GENERIC_POINTER) NULL;
    }
    else if (*haveP == 0){ /* assume we want to initialize heap */
      *haveP = need;
      if (*haveP < chunk) *haveP = chunk;
      heapP = (GENERIC_POINTER) malloc(*haveP * sizeEach);
    }
  }
  return (heapP);
}

int copyProp(PROP *to, PROP *from)
     /* copies the array data over ... check for compatability*/

{
  int status = 0;
  GENERIC_POINTER f=NULL;
  GENERIC_POINTER t=NULL;

  if (!from) return (-1);
  if (!to) return (-2);
	      
  f = from->p;
  t = to->p;

  if (!f || !t || !from->bnd || !to->bnd) return (-3);

  if (from->ptype != to->ptype) return (-11);
  status = checkBounds(to->nInd, to->bnd, from->nInd, from->bnd);

  if (!status) return (copyArray(t,to->nInd,to->bnd, 
				 f,from->nInd, from->bnd, from->ptype));
  else return (status);
}


int copyArray (GENERIC_POINTER t, int nT, int **bndT, 
	       GENERIC_POINTER f, int nF, int **bndF,
	       char pt)
{
  long nTotT = 0, nTotF=0;
  GENERIC_POINTER f0 = NULL;
  GENERIC_POINTER t0 = NULL;

  nTotT = getNTotal(nT, bndT);
  nTotF = getNTotal(nF, bndF);

  if (nTotT <1 || nTotF < 1 || nTotT != nTotF) return (-14);

  f0 = locationOne((GENERIC_POINTER) f, nF, bndF, pt);
  if (!f0) return (-15);
  t0 = locationOne((GENERIC_POINTER) t, nT, bndT, pt);
  if (!t0) return (-16);
  
  if (pt == 'd') return (copyArrayD(t0,f0,nTotT));

  else if (pt == 'f') return (copyArrayF(t0,f0,nTotT));

  else if (pt == 'i') return (copyArrayI(t0,f0,nTotT));

  else return (-21);
}

GENERIC_POINTER locationOne(GENERIC_POINTER x, int nb, int **bnd, char pt)
{
  GENERIC_POINTER y = NULL;

  if (pt == 'd') {
    if (nb == 1) 
      y = (GENERIC_POINTER) &(((double *   ) x)[bnd[1][1]]);
    else if (nb == 2)
      y = (GENERIC_POINTER) &(((double **  ) x)[bnd[1][1]][bnd[2][1]]);
    else if (nb == 3)
      y = (GENERIC_POINTER) &(((double *** ) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]]);
    else if (nb == 4)
      y = (GENERIC_POINTER) &(((double ****) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]][bnd[4][1]]);
  }
  if (pt == 'f') {
    if (nb == 1) 
      y = (GENERIC_POINTER) &(((float  *   ) x)[bnd[1][1]]);
    else if (nb == 2)
      y = (GENERIC_POINTER) &(((float  **  ) x)[bnd[1][1]][bnd[2][1]]);
    else if (nb == 3)
      y = (GENERIC_POINTER) &(((float  *** ) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]]);
    else if (nb == 4)
      y = (GENERIC_POINTER) &(((float  ****) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]][bnd[4][1]]);
  }
  if (pt == 'i') {
    if (nb == 1) 
      y = (GENERIC_POINTER) &(((int    *   ) x)[bnd[1][1]]);
    else if (nb == 2)
      y = (GENERIC_POINTER) &(((int    **  ) x)[bnd[1][1]][bnd[2][1]]);
    else if (nb == 3)
      y = (GENERIC_POINTER) &(((int    *** ) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]]);
    else if (nb == 4)
      y = (GENERIC_POINTER) &(((int    ****) x)[bnd[1][1]][bnd[2][1]][bnd[3][1]][bnd[4][1]]);
  }
  return (y);
}
int copyArrayD(GENERIC_POINTER tG, GENERIC_POINTER fG, long n)
{
  /* NOTE, this code assumes allocation in contiguous, e.g., using d4tensor!!! */
  double *t = NULL, *f = NULL;
  long i;

  t = (double *) tG;
  f = (double *) fG;
  if (!t || !f || n<1) return (1);

  for (i=0;i<n;i++) t[i]=f[i];
  return (0);
}
int copyArrayF(GENERIC_POINTER tG, GENERIC_POINTER fG, long n)
{
  /* NOTE, this code assumes allocation in contiguous, e.g., using d4tensor!!! */
  float *t = NULL, *f = NULL;
  long i;

  t = (float *) tG;
  f = (float *) fG;
  if (!t || !f || n<1) return (1);

  for (i=0;i<n;i++) t[i]=f[i];
  return(0);
}
int copyArrayI(GENERIC_POINTER tG, GENERIC_POINTER fG, long n)
{
  /* NOTE, this code assumes allocation in contiguous, e.g., using d4tensor!!! */
  int *t = NULL, *f = NULL;
  long i;

  t = (int *) tG;
  f = (int *) fG;
  if (!t || !f || n<1) return (1);

  for (i=0;i<n;i++) t[i]=f[i];
  return(0);
}
long getNTotal (int nD, int **bnd)
{
  long n = 0;
  int i;

  n=(long) 1;
  if (nD<1 || nD>4) return (-1);
  for (i=1;i<=nD;i++) {
    if (bnd[i][1] > bnd[i][2]) return ((long) -2);
    n *= (long) (bnd[i][2]-bnd[i][1]+1);
  }
  return (n);
}
/*************************************************************************/
GENERIC_POINTER allocAnyArray(char ptype, int nInd, int **bnd)
/*************************************************************************/
{
  const char routine[] = "allocAnyArray";
  GENERIC_POINTER p = NULL;

  switch (nInd) {
  case 1:
    p = allocAnyArray1(ptype, nInd, bnd);
    break;
  case 2:
    p = allocAnyArray2(ptype, nInd, bnd);
    break;
  case 3:
    p = allocAnyArray3(ptype, nInd, bnd);
    break;
  case 4:
    p =allocAnyArray4(ptype, nInd, bnd);
    break;
  default:
    p = (GENERIC_POINTER) NULL;
    break;
  }
  return p;
}
GENERIC_POINTER allocAnyArray1(char ptype, int nInd, int **bnd)
{
  GENERIC_POINTER p = NULL;

  if (nInd != 1 ||
      bnd[1][1]>bnd[1][2]) 
    p = (GENERIC_POINTER) NULL;

  else {
    switch (ptype) {
    case 'i':
      p = (GENERIC_POINTER) ivector(bnd[1][1], bnd[1][2]);
      break;
    case 'f':
      p = (GENERIC_POINTER) vector(bnd[1][1], bnd[1][2]);
      break;
    case 'd':
      p = (GENERIC_POINTER) dvector(bnd[1][1], bnd[1][2]);
      break;
    default:
      p = (GENERIC_POINTER) NULL;
      break;
    }
  } /* end else */
  return (p);
}
GENERIC_POINTER allocAnyArray2(char ptype, int nInd, int **bnd)
{
  GENERIC_POINTER p = NULL;

  if (nInd != 2 ||
      bnd[1][1]>bnd[1][2] || 
      bnd[2][1]>bnd[2][2]) 
    p = (GENERIC_POINTER) NULL;

  else {
    switch (ptype) {
    case 'i':
      p = (GENERIC_POINTER) imatrix(bnd[1][1], bnd[1][2], 
				    bnd[2][1], bnd[2][2]);
      break;
    case 'f':
      p = (GENERIC_POINTER) matrix(bnd[1][1], bnd[1][2],
				   bnd[2][1], bnd[2][2]);
      break;
    case 'd':
      p = (GENERIC_POINTER) dmatrix(bnd[1][1], bnd[1][2],
				    bnd[2][1], bnd[2][2]);
      break;
    default:
      p = (GENERIC_POINTER) NULL;
      break;
    } /* end else */
  }
  return (p);
}
GENERIC_POINTER allocAnyArray3(char ptype, int nInd, int **bnd)
{
  GENERIC_POINTER p = NULL;

  if (nInd != 3 ||
      bnd[1][1]>bnd[1][2] || 
      bnd[2][1]>bnd[2][2] ||
      bnd[3][1]>bnd[3][2]) 
    p = (GENERIC_POINTER) NULL;

  else {
    switch (ptype) {
    case 'i':
      p = (GENERIC_POINTER) i3tensor(bnd[1][1], bnd[1][2], 
				     bnd[2][1], bnd[2][2],
				     bnd[3][1], bnd[3][2]);
      break;
    case 'f':
      p = (GENERIC_POINTER) f3tensor(bnd[1][1], bnd[1][2],
				     bnd[2][1], bnd[2][2],
				     bnd[3][1], bnd[2][2]);
      break;
    case 'd':
      p = (GENERIC_POINTER) d3tensor(bnd[1][1], bnd[1][2],
				     bnd[2][1], bnd[2][2],
				     bnd[3][1], bnd[3][2]);
      break;
    default:
      p = (GENERIC_POINTER) NULL;
      break;
    }
  } /* end else */
  return (p);
}
GENERIC_POINTER allocAnyArray4(char ptype, int nInd, int **bnd)
{
  GENERIC_POINTER p = NULL;
  
  if (nInd != 4 ||
      bnd[1][1]>bnd[1][2] || 
      bnd[2][1]>bnd[2][2] ||
      bnd[3][1]>bnd[3][2] ||
      bnd[4][1]>bnd[4][2]) 
    p = (GENERIC_POINTER) NULL;

  else {
    switch (ptype) {
    case 'i':
      p = (GENERIC_POINTER) i3tensor(bnd[1][1], bnd[1][2], 
				     bnd[2][1], bnd[2][2],
				     bnd[3][1], bnd[3][2]);
      break;
    case 'f':
    p = (GENERIC_POINTER) f3tensor(bnd[1][1], bnd[1][2],
				   bnd[2][1], bnd[2][2],
				   bnd[3][1], bnd[2][2]);
    break;
    case 'd':
      p = (GENERIC_POINTER) d3tensor(bnd[1][1], bnd[1][2],
				     bnd[2][1], bnd[2][2],
				     bnd[3][1], bnd[3][2]);
      break;
    default:
      p = (GENERIC_POINTER) NULL;
      break;
    }
  } /* end else */
  return (p);
}  
int freeAnyArray(GENERIC_POINTER x, char ptype, int nInd, int **bnd)
     /* returns 0 if x is nonNull and free works
	1 if x is null
	negative numbers for various indexing errors
     */
{
  int iReturn = 0;

  if (x) iReturn = 1;
  
  else {
    iReturn = 0;
    
    switch (ptype) {

    case 'i':
      if (nInd == 1) {
		if (bnd[1][1]>bnd[1][2]) iReturn = -11;
		else free_ivector((int *) x, 
						  bnd[1][1],bnd[1][2]);
      }
      else if (nInd == 2) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ) iReturn = -12;
		else free_imatrix((int **) x,
						  bnd[1][1],bnd[1][2],
						  bnd[2][1],bnd[2][2]);
      }
      else if (nInd == 3) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ) iReturn = -13;
		else free_i3tensor((int ***) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2]);
      }
      else if (nInd == 4) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ||
			bnd[4][1]>bnd[4][2] ) iReturn = -14;
		else free_i4tensor((int ****) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2],
						   bnd[4][1],bnd[4][2]);
      }
      else iReturn = -15;
      break;

    case 'f':
      if (nInd == 1) {
		if (bnd[1][1]>bnd[1][2]) iReturn = -21;
		else free_vector((float *) x, 
						 bnd[1][1],bnd[1][2]);
      }
      else if (nInd == 2) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ) iReturn = -22;
		else free_matrix((float **) x,
						 bnd[1][1],bnd[1][2],
						 bnd[2][1],bnd[2][2]);
      }
      else if (nInd == 3) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ) iReturn = -23;
		else free_f3tensor((float ***) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2]);
      }
      else if (nInd == 4) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ||
			bnd[4][1]>bnd[4][2] ) iReturn = -24;
		else free_f4tensor((float ****) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2],
						   bnd[4][1],bnd[4][2]);
      }
      else iReturn = -25;
      break;

    case 'd':
      if (nInd == 1) {
		if (bnd[1][1]>bnd[1][2]) iReturn = -31;
		else free_dvector((double *) x, 
						  bnd[1][1],bnd[1][2]);
      }
      else if (nInd == 2) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ) iReturn = -32;
		else free_dmatrix((double **) x,
						  bnd[1][1],bnd[1][2],
						  bnd[2][1],bnd[2][2]);
      }
      else if (nInd == 3) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ) iReturn = -33;
		else free_d3tensor((double ***) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2]);
      }
      else if (nInd == 4) {
		if (bnd[1][1]>bnd[1][2] ||
			bnd[2][1]>bnd[2][2] ||
			bnd[3][1]>bnd[3][2] ||
			bnd[4][1]>bnd[4][2] ) iReturn = -34;
		else free_d4tensor((double ****) x,
						   bnd[1][1],bnd[1][2],
						   bnd[2][1],bnd[2][2],
						   bnd[3][1],bnd[3][2],
						   bnd[4][1],bnd[4][2]);
      }
      else iReturn = -35;
      break;

    default:
      iReturn = - 41;
      break;
    }
  }
  return (iReturn);
}
int checkBounds(int nNeed, int **bNeed, int nHave, int **bHave)
     /* If there is a mismatch in bnd arrays, trys to straighten it out
	sensibly. The only mismatch that can be corrected is if at least
	nNeed - nHave dimensions in bNeed are 1. This is checked.

	Main use: resetting a two-D bounds for an I,J areal property to an I,J,1 grid
	If successful, then copyArray can be used to move data described by bHave
	to array described by bNeed

	NOTE: no attempt is made to "fix" bNeed or bHave, only "reasonableness"
	and array size are checked!

	Returns the sum of 
	0 on success
	1 if there are too many non-unit dimensions in bNeed
	2 if the array sizes are different
	4 if non-unit dimensions have different lengths
     */
{
  int i,j;
  int nOneInNeed, nDiff;
  long sizeNeed, sizeHave;
  int dimI;
  int ndN, ndH;
  int dimNeed[5], dimHave[5];
  int iReturn = 0;

  nOneInNeed = 0;
  sizeNeed = 1;
  sizeHave = 1;
  ndN = 0; ndH = 0;

  for (i=1;i<=nNeed;i++) {
    dimI = bNeed[i][2] - bNeed[i][1]+1;
    if (dimI>1) dimNeed[++ndN] = dimI;
    else nOneInNeed++;
    sizeNeed *= dimI;
  }
  for (i=1;i<=nHave;i++) {
    dimI = bHave[i][2] - bHave[i][1]+1;
    if (dimI>1) dimHave[++ndH] = dimI;
    else nOneInNeed++;
    sizeHave *= dimI;
  }

  if ((nOneInNeed < (nNeed - nHave)) || 
      ndH != ndN ) iReturn +=  1;
  if (sizeNeed != sizeHave) iReturn += 2;
  
  if (ndH == ndN) {
    for (i=1;i<=ndH;i++) {
      if (dimHave[i]-dimNeed[i]) {
	iReturn += 4;
	break;
      }
    }
  }

  return (iReturn);
}
int checkBoundsStrict(int nNeed, int **bNeed, int nHave, int **bHave)
     /* 
	Must have a match on dimensions and rank
	0 indicates match
	1 mismatch in rank
	2 mismatch in one or more dimensions
     */
{
  int i;

  if (nNeed != nHave) return (1);
  for (i=1;i<=nNeed;i++) {
    if ((bNeed[i][2] - bNeed[i][1]) != (bHave[i][2] - bHave[i][1])) 
      return (2);
  }
  
  return (0);
}
int checkBoundsVeryStrict(int nNeed, int **bNeed, int nHave, int **bHave)
     /* 
	Must have a match on dimensions, rank, and offset
	0 indicates match
	1 mismatch in rank
	2 mismatch in one or more dimensions
	3 mismatch in one or more offsets
     */
{
  int i;

  if (nNeed != nHave) return (1);
  for (i=1;i<=nNeed;i++) {
    if ((bNeed[i][2] - bNeed[i][1]) != (bHave[i][2] - bHave[i][1])) 
      return (2);
  }
  for (i=1;i<=nNeed;i++) {
    if ((bNeed[i][1] != bHave[i][1])) return (3);
  }
  
  return (0);
}
int checkArrayBounds(int nNeed, int **bNeed, int nHave, int **bHave, int levels)
     /* Detects mismatches in array sizes
	level 0 no check
	level 1 check that nonunit array sizes are consistent
	    No requirement that nNeed and nHave are the same, or that offsets are
	    this allows different ranks so long as data area is the same size
	level 2 check that the rank and dimensions are the same
	    Allows different offsets, but rank and data area must be the same
	level 3 check rank, dimensions, and offsets
	   Need and Have arrays must be of identical declaration
	Returns 0 for pass, 1 for failure
     */
{
  int i,j;
  long sizeNeed, sizeHave;
  int dimI;
  int ndN, ndH;
  int dimNeed[5], dimHave[5];
  int iReturn = 0;

  if (levels == 0) iReturn = 0;

  if (levels == 1) {
    sizeNeed = 1;
    sizeHave = 1;
    ndN = 0; ndH = 0;
    /* Computation of nonzero dimensions and array data area size */
    for (i=1;i<=nNeed;i++) {
      dimI = bNeed[i][2] - bNeed[i][1]+1;
      if (dimI>1) dimNeed[++ndN] = dimI;
      sizeNeed *= dimI;
    }
    for (i=1;i<=nHave;i++) {
      dimI = bHave[i][2] - bHave[i][1]+1;
      if (dimI>1) dimHave[++ndH] = dimI;
      sizeHave *= dimI;
    }
    /* check that the number of nonunit dimensions and are the same */
    if (ndH != ndN) iReturn = 1;
    else if (sizeNeed != sizeHave) iReturn = 1;
    else {
      for (i=1;i<=ndH;i++) {
	if (dimHave[i]!=dimNeed[i]) {
	  iReturn = 1;
	  break;
	}
      }
    }
  }
  
  if (levels >= 2) {
    if (nNeed != nHave) iReturn = 1;
    else {
      for (i=1;i<=nNeed;i++) {
	if ((bNeed[i][2] - bNeed[i][1]) != (bHave[i][2] - bHave[i][1])) {
	  iReturn = 1;
	  break;
	}
      }
    }
  }
  
  if (levels == 3) {
    for (i=1;i<=nNeed;i++) {
      if ((bNeed[i][1] != bHave[i][1])) {
	iReturn = 1;
	break;
      }
    }
  }
  
  return (iReturn);
}
