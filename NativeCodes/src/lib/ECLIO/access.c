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
#include <string.h>
#include <stdarg.h>
#include <stdio.h> 
#include <ctype.h>
#include "../NR/nrutil.h"
#include "ECLIO/eclio.h"
#include "ECLIO/access.h"
#include "ECLIO/eclutil.h"
#include "ECLIO/blockops.h"
#include "ECLIO/records.h"
#include "ECLIO/read.h" 
#include "ECLIO/set.h"
#include "ECLIO/copy.h"
#include "ECLIO/write.h"


// No thread local storage on OSX. Library will not be thread-safe until OSX supports it
// or we rewrite... sigh.
#ifdef __APPLE__
 #define THREAD
#else
 #define THREAD __thread
#endif


static THREAD int _skipPinch=0;

/* Relatively high-level routines to load and access
   Eclipse data records CD White, LSU, 5/04 */


/*************************************************************************/
int loadModel(char *fName, ECL *ecl, int transformToMap,
	      int nKeys, ...)
/*************************************************************************/
/* 
     Loads the variable-length keyword/token/record list that
     follows nKeys (the length of the list) onto the  ECL
     pointed to by ecl. Reads from file fName. Optionally applies
     a transform if MAPAXES is present. Will rescale units if MAPUNITS
     and GRIDUNIT are present. The list must be in the SAME order as in
     the file fName.
*/
{
  FILE *f;
  const char routine[] = "loadModel";
  int toUpper = 1;
  int rewindF = 0;
  int status = 0;
  char **keys = (char **) NULL;
  int i;
  size_t tL = TOKEN_LENGTH ;

  va_list ap;

  if (!ecl) {
    (void) eclError(routine,1,"Model pointer is null, input file = %s\n",
		    fName);
    goto bottom;
  }

  keys = cmatrix((long) 1,(long) nKeys, (long) 0, (long) tL);
  if (!keys) {
    (void) eclError(routine,1,"Could not allocate keys vector\n");
    goto bottom;
  }

  va_start(ap,nKeys);
  for (i=1;i<=nKeys;i++) strncpy(keys[i],va_arg(ap,char *),tL);
  va_end(ap);

  (void) initModel(ecl);
			   
  if (!(f = fopen(fName, "r"))){
    (void) eclError(routine,2,"Could not open model file = %s for input\n",
	     fName);
    goto bottom;
  }

  /* scan through keylist */
  for (i=1;i<=nKeys;i++) {
    if (!findKey(ecl,f, keys[i], toUpper, rewindF)) {
      if (!strncmp(keys[i],"SPECGRID",tL)) status += r_specGrid(f,ecl);
      else if (!strncmp(keys[i],"GRIDUNIT",tL)) status += r_gridUnit(f,ecl);
      else if (!strncmp(keys[i],"PINCH",tL)) status += r_pinch(f,ecl);
      else if (!strncmp(keys[i],"MAPAXES",tL)) status += r_mapAxes(f,ecl);
      else if (!strncmp(keys[i],"MAPUNITS",tL)) status += r_mapUnits(f,ecl);
      else if (!strncmp(keys[i],"COORDSYS",tL)) status += r_coordSys(f,ecl);
      else if (!strncmp(keys[i],"COORD",tL)) status += r_coord(f,ecl);
      else if (!strncmp(keys[i],"ZCORN",tL)) status += r_zCorn(f,ecl);
      /* any other keys are assumed to be GRID PROP ARRAYS
		 ACTNUM, SATNUM, PVTNUM, FIPNUM, and SEGMENTS,
		 (JG_FEATURES: and also ->) LAYERS
		 are int, ptype ='i'
         all others are assumed to be double, ptype = 'd' */
      else {
	if (!strncmp(keys[i],"SATNUM",tL) ||
	    !strncmp(keys[i],"PVTNUM",tL) ||
	    !strncmp(keys[i],"FIPNUM",tL) ||
#if JG_FEATURES // hack around Petrel bug, plus addition of LAYERS
        !strncmp(keys[i],"SEGMENTS",8) ||  
		!strncmp(keys[i],"LAYERS",tL)  ||  
#else
	    !strncmp(keys[i],"SEGMENTS",tL) ||  
#endif
        !strncmp(keys[i],"ACTNUM",tL)) 
	  status= r_propGrid(f,ecl,keys[i],'i');
	else
	  status= r_propGrid(f,ecl,keys[i],'d'); 
      }

      if (status) {
	(void) eclError(routine,3, 
			"Could not read data for key = %s\n",keys[i]);
	goto bottom;
      }
    }
    else {
      (void) eclError(routine,4, "Could not find key = %s\n",keys[i]);
      goto bottom;
    }
  }
  strncpy(ecl->name, fName, NAME_LENGTH);
  
  /* some useful pointers in ECL */
  (void) setAuxPointers(ecl);

  /* Additional operations */
  if (transformToMap) (void) rotateGridMap(ecl);
  status += applyPinch(ecl); /*Applies pinchout criteria if pinch->set*/

  /* Clean up */
 bottom:
  freeLoadModelMem(nKeys, tL, keys);
  compressProp(ecl);
  return (status);
}

/*************************************************************************/
int loadModelFromList(char *fName, ECL *ecl, int transformToMap,
		      int nKeys, ...)
/*************************************************************************/
/* 
     Loads the variable-length keyword/token/record list that
     follows nKeys (the length of the list) onto the  ECL
     pointed to by ecl. Reads from file fName. Optionally applies
     a transform if MAPAXES is present. Will rescale units if MAPUNITS
     and GRIDUNIT are present. The list must be in the SAME order as in
     the file fName.
*/
{
  FILE *f;
  const char routine[] = "loadModelFromList";
  int toUpper = 1;
  int status = 0;
  char **keys = (char **) NULL;
  int i;
  size_t tL = TOKEN_LENGTH ;

  va_list ap;

  if (!ecl) {
    status = 1;
    (void) eclError(routine,status,"Model pointer is null, input file = %s\n",
		    fName);
    goto bottom;
  }

  keys = cmatrix((long) 1,(long) nKeys, (long) 0, (long) tL);
  if (!keys){
    status = 2;
    (void) eclError(routine,status,"Could not allocate keys vector\n");
    goto bottom;
  }
    
  va_start(ap,nKeys);
  for (i=1;i<=nKeys;i++) strncpy(keys[i],va_arg(ap,char *),tL);
  va_end(ap);

  status = loadModelFromArray(fName, ecl, transformToMap,
			     nKeys, keys);
  if (status) {
    status = 3;
    (void) eclError(routine,status,
		    "Could not read all %d keys in loadModelFromList\n", nKeys);
    goto bottom;
  }

 bottom:
  return (status);
}

int loadModelFromArraySansPinch(char *fName, ECL *ecl, int transformToMap,
		      int nKeys, char **keys)
{
  int r;
  _skipPinch = 1;
  r = loadModelFromArray(fName,ecl,transformToMap,nKeys,keys);
  _skipPinch =0;

  return r;
}

/*************************************************************************/
int loadModelFromArray(char *fName, ECL *ecl, int transformToMap,
		      int nKeys, char **keys)
/*************************************************************************/
/* 
     Loads the variable-length keyword/token/record list that
     follows nKeys (the length of the list) onto the  ECL
     pointed to by ecl. Reads from file fName. Here the list is an
     array rather than a variable-length list. Optionally applies
     a transform if MAPAXES is present. Will rescale units if MAPUNITS
     and GRIDUNIT are present. The list must be in the SAME order as in
     the file fName.
*/
{
  FILE *f;
  const char routine[] = "loadModelFromArray";
  int toUpper = 1;
  int rewindF = 0;
  int status = 0;
  int i;
  size_t tL = TOKEN_LENGTH ;

  if (!ecl) {
    status = 1;
    (void) eclError(routine, status, "Model pointer is null, input file = %s\n",
		    fName);
    goto bottom;
  }

  (void) initModel(ecl);
			   
  if (!(f = fopen(fName, "r"))){
    status = 2;
    (void) eclError(routine,status,"Could not open model file = %s for input\n",
	     fName);
    goto bottom;
  }

  /* scan through keylist */
  for (i=1;i<=nKeys;i++) {
    if (!findKey(ecl,f, keys[i], toUpper, rewindF)) {
      if (!strncmp(keys[i],"SPECGRID",tL)) status += r_specGrid(f,ecl);
      else if (!strncmp(keys[i],"GRIDUNIT",tL)) status += r_gridUnit(f,ecl);
      else if (!strncmp(keys[i],"PINCH",tL)) status += r_pinch(f,ecl);
      else if (!strncmp(keys[i],"MAPAXES",tL)) status += r_mapAxes(f,ecl);
      else if (!strncmp(keys[i],"MAPUNITS",tL)) status += r_mapUnits(f,ecl);
      else if (!strncmp(keys[i],"COORDSYS",tL)) status += r_coordSys(f,ecl);
      else if (!strncmp(keys[i],"COORD",tL)) status += r_coord(f,ecl);
      else if (!strncmp(keys[i],"ZCORN",tL)) status += r_zCorn(f,ecl);
      /* any other keys are assumed to be GRID PROP ARRAYS
	 ACTNUM, SATNUM, PVTNUM, FIPNUM, and SEGMENTS 
	 (JG_FEATURES: and also ->) LAYERS 
	 are int, ptype ='i'
	 all others are assumed to be double, ptype = 'd' */
      else {
	if (!strncmp(keys[i],"SATNUM",tL) ||
	    !strncmp(keys[i],"PVTNUM",tL) ||
	    !strncmp(keys[i],"FIPNUM",tL) ||
#if JG_FEATURES // hack around Petrel bug
        !strncmp(keys[i],"SEGMENTS",8) ||  
		!strncmp(keys[i],"LAYERS",tL)  ||  
#else
	    !strncmp(keys[i],"SEGMENTS",tL) ||  
#endif
	    !strncmp(keys[i],"ACTNUM",tL)) 
	  status= r_propGrid(f,ecl,keys[i],'i');
	else
	  status= r_propGrid(f,ecl,keys[i],'d');
      }

      if (status) {
	status = 3;
	(void) eclError(routine,status, 
			"Could not read data for key = %s\n",keys[i]);
	goto bottom;
      }
    }
    else {
      status = 4;
      (void) eclError(routine,status, "Could not find key = %s\n",keys[i]);
      goto bottom;
    }
  }
  strncpy(ecl->name, fName, NAME_LENGTH);
  
  /* some useful pointers in ECL */
  (void) setAuxPointers(ecl);

  /* Additional operations */
  if (transformToMap) (void) rotateGridMap(ecl);
  if (!_skipPinch) status += applyPinch(ecl); /* Applies pinchout criteria if pinch->set*/
  if (status) {
    status = 5;
    (void) eclError(routine,status, "Error applying requested pinchout = %s\n");
    goto bottom;
  }

  /* Clean up */
  compressProp(ecl);

 bottom:
  return (status);
}
/*************************************************************************/
    void freeLoadModelMem(int nKeys, size_t l, char **keys)
/*************************************************************************/
{
  if(keys && nKeys>0 && l>=0) free_cmatrix(keys,1,nKeys,0,l);
}

/*************************************************************************/
void freeModel( ECL *ecl)
/*************************************************************************/
  /* Frees all storage allocated with structure pointed to by ecl
     Does not free memory for  ECL ecl itself */
{
  int i, nPropOld;
  if (ecl) {
    if (ecl->coordSys.set) {
      freeCoordSys(ecl->coordSys.k1, ecl->coordSys.k2,
		   ecl->coordSys.comp, ecl->coordSys.join,
		   ecl->specGrid.numRes, TOKEN_LENGTH);
    }
    if (ecl->mapAxes.set && ecl->mapAxes.axes) {
      free_dmatrix(ecl->mapAxes.axes, 1, 3, 1, 2);
    }

    if (ecl->prop) 
	  {
		nPropOld = ecl->nProp;
		for (i=nPropOld;i>=1;i--)
		  {
			if (&ecl->prop[i]) deleteProp(ecl, ecl->prop[i].name);
		  }
		if (ecl->prop) free ((char *) ecl->prop);
	  }
  }
}
/*************************************************************************/
PROP *fetchProp(ECL *ecl, char *name)
/*************************************************************************/
  /* Looks up properties by name. Cannot be used for records
     that are not properties; use fetchRecord or fetchNamedRecord
     for those (SPECGRID, etc.).
     returns a pointer to the requested  PROP if successful
     returns NULL on failure */
{
  int i;
  PROP *p = (PROP *) NULL;

  for (i=1;i<=ecl->nProp;i++) 
	{
	  p = &(ecl->prop[i]);
#if JG_FEATURES // hack around Petrel bug
	  if (!strncmp(name,"SEGMENTS",8))
		{
		  if (!strncmp(p->name, name,8)) return (p);
		}
	  else
		{
		  if (!strcmp(p->name, name)) return (p);
		}
#else
	  if (!strcmp(p->name, name)) return (p);
#endif
	}
  
  return (( PROP *) NULL);
}
/*************************************************************************/
GENERIC_POINTER fetchNamedRecord(ECL *ecl, char *name)
/*************************************************************************/
  /* Looks up a named record (as listed in code) and 
     returns a (GENERIC_POINTER ) pointer to it. Caller must recast
     to the appropriate  type, e.g., (SPECGRID *)
     Returns the (GENERIC_POINTER ) pointer if successful
     Otherwise returns (GENERIC_POINTER ) NULL
     **NOTE** the record may not be set, check
     ( WHATEVER *)fetchNamedRecord(ecl,"WHATEVER")->set
     to determine if record was set.
  */
{
  GENERIC_POINTER  t = (GENERIC_POINTER) NULL;

  if (!strcmp(name,"SPECGRID"))      t = (GENERIC_POINTER) &ecl->specGrid; 
  else if (!strcmp(name,"GRIDUNIT")) t = (GENERIC_POINTER) &ecl->gridUnit; 
  else if (!strcmp(name,"PINCH"))    t = (GENERIC_POINTER) &ecl->pinch; 
  else if (!strcmp(name,"MAPAXES"))  t = (GENERIC_POINTER) &ecl->mapAxes; 
  else if (!strcmp(name,"MAPUNITS")) t = (GENERIC_POINTER) &ecl->mapUnits; 
  else if (!strcmp(name,"COORDSYS")) t = (GENERIC_POINTER) &ecl->coordSys; 
  else if (!strcmp(name,"COORD"))    t = (GENERIC_POINTER) ecl->coordp; 
  else if (!strcmp(name,"ZCORN"))    t = (GENERIC_POINTER) ecl->zCornp; 
  else if (!strcmp(name, "ACTNUM"))  t = (GENERIC_POINTER) ecl->actnump;
  else t = (GENERIC_POINTER) NULL;

  return (t);
}
/*************************************************************************/
GENERIC_POINTER fetchAnyRecord( ECL *ecl, char *name, char *rtype)
/*************************************************************************/
  /* Usage as for fetchNamedRecord, except it will look in 
     the PROP record as well 
     The flag rtype is returned as 'n' for named record or 'p' 
     for a property; '0' if not found
  */
{
  GENERIC_POINTER t = (GENERIC_POINTER) NULL;

  *rtype = 'n';
  t = fetchNamedRecord(ecl, name);
  if (!t) {
    t = (GENERIC_POINTER ) fetchProp(ecl, name);
    if (t) *rtype = 'p';
    else *rtype = '0';
  }
  return ((GENERIC_POINTER) t);
}
/*************************************************************************/
GENERIC_POINTER fetchPropArray(ECL *ecl, char *name)
/*************************************************************************/
/* similar to fetchProp, but returns a pointer to THE DATA
     rather than THE PROP. Returns (GENERIC_POINTER), must be reset to 
     the right level of indirection and type by the caller; e.g.,
     (int ***) for ACTNUM. If uncertain of type, etc., call 
     fetchProp and get nInd and ptype from the returned (PROP *)
*/
{
  GENERIC_POINTER p = (GENERIC_POINTER) NULL;
  PROP *pp = (PROP *) NULL;

  const char routine[] = "fetchPropArray";
  
  if (!(pp = fetchProp(ecl, name))){
    eclError(routine, 1, "Cannot fetch property to get array %s\n", name);
    return ((GENERIC_POINTER) NULL);
  }

  if (!(p = pp->p)) {
    eclError(routine, 2, "Property %s exists but has null pointer to array\n");
    return ((GENERIC_POINTER) NULL);
  }

  return (p);
}

/*************************************************************************/
double ***fetchTraces(ECL *ecl, int skipVoids,
		      int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
	 /* 
		A wrapper for prepareTraces 
	For the model pointed to by ecl,
	Computes the location of the top and bottom of each xy corner
	in the range [i<-nrl..nrh][j<-ncl..nch]
        if (skipVoids) it will not include top/bottom voids (recommend=1)
	returns a (double ***) allocated with 
	NR/nrutil.c:d3tensor(nrl,nrh,ncl,nch,1,6)
	where 
	fetchTraces[i][j][1] = xtop[i][j]
	fetchTraces[i][j][2] = ytop[i][j]
	fetchTraces[i][j][3] = ztop[i][j]
	fetchTraces[i][j][5] = xbot[i][j]
	fetchTraces[i][j][5] = ybot[i][j]
	fetchTraces[i][j][6] = zbot[i][j]

creates also
 
traceP->
*/
{
  return fetchLayerTraces(ecl,-1,skipVoids, nrl, nrh, ncl, nch);
}

/*************************************************************************/
double ***fetchLayerTraces(ECL *ecl, int eclLayerNumber, int skipVoids,
						   int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
/* 
A wrapper for prepareTraces 
	For the model pointed to by ecl,
	Computes the location of the top and bottom of each xy corner
	in the range [i<-nrl..nrh][j<-ncl..nch]
        if (skipVoids) it will not include top/bottom voids (recommend=1)
	returns a (double ***) allocated with 
	NR/nrutil.c:d3tensor(nrl,nrh,ncl,nch,1,6)
	where 
	fetchTraces[i][j][1] = xtop[i][j]
	fetchTraces[i][j][2] = ytop[i][j]
	fetchTraces[i][j][3] = ztop[i][j] 
	fetchTraces[i][j][5] = xbot[i][j]
	fetchTraces[i][j][5] = ybot[i][j]
	fetchTraces[i][j][6] = zbot[i][j]

creates also

traceP->
*/
{
  const char routine[] = "fetchTraces";
  PROP *traceP = (PROP *) NULL;
  double ***trace = (double ***) NULL;
  double ***t = (double ***) NULL;
  int i, j, k;
  int status = 0;

  if (nrl < 0 || nrh > ecl->specGrid.nx || nrl > nrh ||
      ncl < 0 || nch > ecl->specGrid.ny || ncl > nch) {
    (void) eclError(routine,status, 
		    "Nonsense for range of traces\n\tx %d to %d\n\ty %d to %d\n",
		    nrl, nrh, ncl, nch);
    goto bottom;
  }

  /* Get record */
#if JG_FEATURES
  /* force calculation as locations are layer-specific */
  traceP = prepareTracesAndVolumes(ecl, 1, skipVoids, eclLayerNumber);
#else
  traceP = prepareTraces(ecl, 0, skipVoids);
#endif
  if (!traceP) {
    eclError(routine,status,"Could not prepare traces NEWCOORD\n");
    goto bottom;
  }

  trace = (double ***) traceP->p;
  if (!traceP) {
    (void) eclError(routine,status,"Could not access data in NEWCOORD\n");
    goto bottom;
  }

  /* move into storage outside of  ECL */
  t = d3tensor(nrl, nrh, ncl, nch, 1, 6);
  if (!t) eclError(routine,3, "Could not allocate TRACES\n");



  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=1;k<=6;k++) t[i][j][k] = trace[i+nrl][j+ncl][k];
    }
  }

 bottom:
  return (t);
}


#if JG_FEATURES

/*************************************************************************/
double **fetchTraceAreas(ECL *ecl, int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
/* 
	For the model pointed to by ecl,
	Computes the average area of the model region formed by the sum of the 
	"quarter" gridblocks surrounding the trace i,j.
	Location of the top and bottom of each xy corner
	in the range [i<-nrl..nrh][j<-ncl..nch]

	Returns a (double **) allocated with 
	NR/nrutil.c:dmatrix(nrl,nrh,ncl,nch)
	where 
	fetchTraceAreas[i][j] = area around trace i,j

	Usage: call after fetchTraces(...) with exactly the same arguments
*/
{
  const char routine[] = "fetchTraceAreas";
  double **areas = (double **) NULL;
  double ***area = (double ***) NULL;
  PROP *aP = (PROP *) NULL;
  int i, j;
  int status = 0;

  if (nrl < 0 || nrh > ecl->specGrid.nx || nrl > nrh ||
      ncl < 0 || nch > ecl->specGrid.ny || ncl > nch) 
	{
	  (void) eclError(routine,status, 
					  "Nonsense for range of traces\n\tx %d to %d\n\ty %d to %d\n",
					  nrl, nrh, ncl, nch);
	 return (double **)NULL;
	}

  aP = fetchProp(ecl, "TRACEAREAS"); /* exists */
  if (!(aP && aP->p && aP->set)) 
	{
	  status = 1;
      eclError(routine,status,"Could not get TRACEAREAS; call fetchTraces() first\n");
      return (double **)NULL;
    }


  area = (double ***) aP->p;

  /* move into storage outside of  ECL */
  areas = dmatrix(nrl, nrh, ncl, nch);
  if (!areas) eclError(routine,3, "Could not allocate TRACE_AREAS\n");



  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      areas[i][j] = area[i+nrl][j+ncl][0];
    }
  }
  return areas;
}
#endif


/*************************************************************************/
int getNX(ECL *ecl)
/*************************************************************************/
{
  if (ecl->specGrid.set) return (ecl->specGrid.nx);
  else return (0);
}
/*************************************************************************/
int getNY(ECL *ecl)
/*************************************************************************/
{
  if (ecl->specGrid.set) return (ecl->specGrid.ny);
  else return (0);
}
/*************************************************************************/
int getNZ(ECL *ecl)
/*************************************************************************/
{
  if (ecl->specGrid.set) return (ecl->specGrid.nz);
  else return (0);
}
/*************************************************************************/
PROP *addProp(ECL *ecl, char *name, char ptype, 
	      int nInd, ...)
/*************************************************************************/
     /* returns a pointer to a new element in an array of PROP */
{
  int **bnd = (int **) NULL;
  va_list ap; 

  va_start(ap, nInd);  /* Build bnd array */ 
  bnd = buildBnd(ecl, nInd, ap, name);
  va_end(ap);

  return (addPropN(ecl, name, ptype, nInd, bnd));
} 
/*************************************************************************/
PROP *addPropGrid(ECL *ecl, char *name, char ptype)
/*************************************************************************/
     /* returns a pointer to a new element in an array of PROP */

  /*  Adds a grid property, [1..nx][1..ny][1..nz] to the PROP vector 
      of model pointed to by ecl. the new PROP will have name "name" and
      be of type 'ptype' =['i','f','d']
   */

{
  int **bnd = (int **) NULL;
  int nInd;

  bnd = buildBndGrid(ecl);
  nInd = 3;
  return (addPropN(ecl, name, ptype, nInd, bnd));
}
/*************************************************************************/ 
int freeProp(ECL *ecl, char *name)
/*************************************************************************/
     /* Frees the data for a PROP named "name", but does not
	deallocate storage for the PROP itself in ecl->prop
	(see deleteProp)
	returns 0 if freed, 1 if not found */
{
  PROP *p = (PROP *) NULL;
  p = fetchProp(ecl, name);
  if (p) {
    free_anyProp(p);
    return (0);
  }
  return (1);
}
/*************************************************************************/
int deleteProp(ECL *ecl, char *name)
/*************************************************************************/

  /* Frees the data associated with a PROP named "name" in ecl, and
     deletes it from ecl->prop, reallocating ecl->prop.
     returns 0 if deleted, 1 if not found, 2 if realloc failed */
{
  int i;
  int iProp = 0;
  PROP *pv = (PROP *) NULL;
  
  iProp = fetchPropInd(ecl, name);
  if (!iProp) return (1);
  free_anyProp(&ecl->prop[iProp]);
  pv = ecl->prop;
  for (i=iProp;i<=ecl->nProp-1;i++) pv[i]=pv[i+1];
  ecl->nProp--;
  return (compressProp(ecl));
}
/*************************************************************************/
int **fetchActCol(ECL *ecl, int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
{
  PROP *aColP=NULL;
  int **aCol=NULL;
  int **t=NULL;
  char routine[] = "fetchActCol";
  int i, j;
  int di = nrl-1, dj = ncl-1;

  if (nrl < 0 || nrh > ecl->specGrid.nx ||
      ncl < 0 || nch > ecl->specGrid.ny ) {
    (void) eclError(routine,1,
		    "Nonsense bounds [nrl,nrh;ncl,nch] = [%d,%d;%d,%d]\n",
		    nrl,nrh,ncl,nch);
    goto bottom;
  }

  aColP = prepareActCol(ecl, 0);
  if (!aColP || !aColP->p || !aColP->set){
    (void) eclError(routine,2, "Could not set ACTCOL\n");
    goto bottom;
  }

  t = imatrix(nrl, nrh,ncl, nch);
  if (!t) {
    (void) eclError(routine,3, "Could not allocate storage for fetch\n");
    goto bottom;
  }

  aCol = (int **) aColP->p;

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) t[i][j] = aCol[i+di][j+dj];
  }

 bottom:
  return (t);
}
/*************************************************************************/
int ***fetchActnum(ECL *ecl, 
		   int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/*************************************************************************/
{
  PROP *aP = NULL; 
  int ***a = NULL;
  int ***t = NULL;
  char routine[] = "fetchActnum";
  int i, j, k;
  int di=nrl-1, dj=ncl-1, dk=ndl-1;

  if (nrl < 0 || nrh > ecl->specGrid.nx ||
      ncl < 0 || nch > ecl->specGrid.ny ||
      ndl < 0 || ndh > ecl->specGrid.nz ) {
    eclError(routine,1,
	     "Nonsense bounds [nrl,nrh;ncl,nch;ndl,ndh] = "
	     "[%d,%d;%d,%d;%d,%d]\n",
	     nrl,nrh,ncl,nch, ndl, ndh);
    goto bottom;
  }

  aP = fetchProp(ecl, "ACTNUM");
  if (!aP || !aP->p || !aP->set){
    (void) eclError(routine,2, "Could not get ACTNUM\n");
    goto bottom;
  }

  t = i3tensor(nrl, nrh, ncl, nch, ndl, ndh);
  if (!t) {
    (void) eclError(routine,3, "Could not allocate fetch array\n");
    goto bottom;
  }

  a = (int ***) aP->p;

  for (i=nrl;i<=nrh; i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) t[i][j][k] = a[i+di][j+dj][k+dk];
    }
  }
 bottom:
  return (t);
}

/*************************************************************************/
int **fetchActCoord(ECL *ecl, int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
{
  /* Returns an imatrix indicating whether each COORD has any active blocks
     along. 
     Will create new PROPs ACTCOORD and ACTCOL if needed.
     Returns pointer to non-NULL for success, else NULL
  */

  PROP *actCoord;
  char routine[] ="fetchActCoord";
  int **aC;
  int **t;
  int i, j;
  int di = nrl-0, dj=ncl-0;

  if (nrl < 0 || nrh > ecl->specGrid.nx ||
      ncl < 0 || nch > ecl->specGrid.ny ) {
    eclError(routine, 1,
	     "Nonsense bounds [nrl,nrh;ncl,nch] = [%d,%d;%d,%d]\n",
	     nrl,nrh,ncl,nch);
    goto bottom;
  }

  actCoord = prepareActCoord(ecl, 0);
  
  if (!actCoord || !actCoord->p || !actCoord->set){
    eclError(routine, 2,
	     "Could not prepare the array.\n");
    goto bottom; 
  }

  aC = (int **) actCoord->p;
  if (!aC) {
    eclError(routine, 3,
	     "NULL data array.\n");
    goto bottom;
  }


  t = imatrix(nrl,nrh,ncl,nch);
  if (!t) {
    eclError(routine, 4,
	     "NULL results array.\n");
    goto bottom;
  }

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) t[i][j]=aC[i+di][j+dj]; 
  }

 bottom:
  return (t);
}    
/*************************************************************************/
int ***fetchTopBot(ECL *ecl, 
		   int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
{
  PROP *aP = NULL;
  int ***a = NULL;
  int ***t = NULL;
  char routine[] = "fetchTopBot";
  int i, j;
  int di=nrl-1, dj=ncl-1;

  if (nrl < 1 || nrh > ecl->specGrid.nx ||
      ncl < 1 || nch > ecl->specGrid.ny){
    eclError(routine,1,
	     "Nonsense bounds [nrl,nrh;ncl,nch] = "
	     "[%d,%d;%d,%d]\n",
	     nrl,nrh,ncl,nch);
    goto bottom;
  }

  aP = fetchProp(ecl, "ACTNUM");
  if (!aP || !aP->p || !aP->set) {
    (void) eclError(routine,2, "Could not get ACTNUM\n");
    goto bottom;
  }


  a = (int ***) aP->p;
  if (!a) {
    (void) eclError(routine,3, "Null data area for ACTNUM\n");
    goto bottom;
  }

  t = i3tensor(nrl, nrh, ncl, nch, 1, 2);
  if (!t) {
    (void) eclError(routine,4, "Could not allocate fetch array\n");
    goto bottom;
  }

  for (i=nrl;i<=nrh; i++) {
    for (j=ncl;j<=nch;j++) {
      t[i+di][j+dj][1] = topAct(ecl,i,j);
      t[i+di][j+dj][2] = botAct(ecl,i,j);
    }
  }

 bottom:
  return (t);
}
/*************************************************************************/
int ***fetchSegments(ECL *ecl, 
		     int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/*************************************************************************/
{
  PROP *aP = NULL;
  int ***a = NULL;
  int ***t = NULL;
  char routine[] = "fetchSegments";
  int i, j, k;
  int di=nrl-1, dj=ncl-1, dk = ndl-1;
  int nx, ny, nz;
  int status;

  status = checkSpecGrid(ecl, &nx, &ny, &nz);
  
  if (status ||
      nrl < 1 || nrh > nx || nrl > nrh ||
      ncl < 1 || nch > ny || ncl > nch ||
      ndl < 1 || ndh > nz || ndl > ndh ) {
    eclError(routine,1,
	     "Nonsense bounds [nrl,nrh;ncl,nch;ndl,ndh] = "
	     "[%d,%d;%d,%d;%d,%d]\n"
	     "[nx,ny,nz]=[%d,%d,%d]",
	     nrl,nrh,ncl,nch, ndl, ndh,nx,ny,nz);
    goto bottom;
  }

  aP = fetchProp(ecl, "SEGMENTS");
  if (!aP || !aP->p || !aP->set) {
    (void) eclError(routine,2, "Could not get SEGMENTS\n");
    goto bottom;
  }
  a = aP->p;

  t = i3tensor(nrl, nrh, ncl, nch, ndl, ndh);
  if (!t) {
    (void) eclError(routine,3, "Could not allocate fetch array\n");
    goto bottom;
  }

  for (i=nrl;i<=nrh; i++) {
    for (j=ncl;j<=nch;j++) {
      for (k=ndl;k<=ndh;k++) {
	t[i+di][j+dj][k+dk]=a[i][j][k];
      }
    }
  }

 bottom:
  return (t);
}
/*************************************************************************/
int **fetchSegCol(ECL *ecl, 
		   int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
{
  PROP *aP = NULL;
  int **a = NULL;
  int **t= NULL;
  char routine[] = "fetchSegCol";
  int i, j;
  int di=nrl-1, dj=ncl-1;

  if (nrl < 1 || nrh > ecl->specGrid.nx ||
      ncl < 1 || nch > ecl->specGrid.ny) {
    eclError(routine,1,
 	     "Nonsense bounds [nrl,nrh;ncl,nch] = "
	     "[%d,%d;%d,%d]\n",
	     nrl,nrh,ncl,nch);
    goto bottom; 
  }

  aP = prepareSegCol(ecl,0);

  if (!aP || !aP->p || !aP->set) {
    (void) eclError(routine,2, "Could not get SEGMENTS\n");
    goto bottom;
  }

  t = imatrix(nrl, nrh, ncl, nch);
  if (!t) {
    (void) eclError(routine,3, "Could not allocate fetch array\n");
    goto bottom;
  }

  a = (int **) aP->p;

  for (i=nrl;i<=nrh; i++) {
    for (j=ncl;j<=nch;j++) {
      t[i+di][j+dj] = a[i][j];
    }
  }

 bottom:
  return (t);
}
/*************************************************************************/
double ****fetchTrace4(ECL *ecl, int skipVoids,
		       int nrl, int nrh, int ncl, int nch)
/*************************************************************************/
/* 
A wrapper for prepareTrace4 
	For the model pointed to by ecl,
	Computes the location of the top and bottom of each xy corner
	in the range [i<-nrl..nrh][j<-ncl..nch]
        if (skipVoids) it will not include top/bottom voids (recommend=1)
	returns a (double ****) allocated with 
	NR/nrutil.c:d4tensor(nrl,nrh,ncl,nch,1,4,1,6)
	where 
	fetchTrace4[i][j][1][1] = xtop[i-1][j-1]
	fetchTrace4[i][j][1][2] = ytop[i-1][j-1]
	fetchTrace4[i][j][1][3] = ztop[i-1][j-1]
	fetchTrace4[i][j][1][5] = xbot[i-1][j-1]
	fetchTrace4[i][j][1][5] = ybot[i-1][j-1]
	fetchTrace4[i][j][1][6] = zbot[i-1][j-1]
	...and the third index is [2] for [i][j-1]
	                          [3] for [i-1][j]
				              [4] for [i][j]
        values off the grid are returned as zeroes */
{
  const char routine[] = "fetchTrace4";
  PROP *traceP = NULL;
  double ****trace = NULL;
  double ****t = NULL;
  int i, j, k, l;

  if (nrl < 0 || nrh > ecl->specGrid.nx || nrl > nrh ||
      ncl < 0 || nch > ecl->specGrid.ny || ncl > nch) {
    eclError(routine, 1, 
	     "Nonsense for range of traces\n\tx %d to %d\n\ty %d to %d\n",
	     nrl, nrh, ncl, nch);
    goto bottom;
  }

  /* Get record */
  traceP = prepareTrace4(ecl, 0, skipVoids);
  if (!traceP) {
    eclError(routine,2,"Could not prepare traces NEWCOOR4\n");
    goto bottom;
  }

  trace = (double ****) traceP->p;
  if (!traceP) {
    eclError(routine,3,"Could not access data in NEWCOOR4\n");
    goto bottom;
  }

  /* move into storage outside of  ECL */
  t = d4tensor(nrl, nrh, ncl, nch, 1,4,1, 6);
  if (!t) {
    eclError(routine,3, "Could not allocate TRACE4\n");
    goto bottom;
  }

  for (i=nrl;i<=nrh;i++) {
    for (j=ncl;j<=nch;j++) {
      for (l=1;l<=4;l++) {
		for (k=1;k<=6;k++) t[i][j][l][k] = trace[i+nrl][j+ncl][l][k];
      }
    }
  }

 bottom:
  return (t);
}
char contFlags(double *t1, double *t2)
{
  int i;
  char c='T';
  for (i=1;i<=6;i++) {
    if (t1[i]!=t2[i]) {
      c='F';
      break;
    }
  }
  return (c);
}

int fetchKeyList(char *fName, char ***keyWordsP, char *commentString)
     /* returns pointer to a char **[0:getKeyList()-1] list of 
		(probable) encountered keyWords  ***in keyWordsP***
		returns number of keyWords as return value UNLESS
		-1 for allocation error
		-2 for file open error
     */
{
  int nKeyWords, nw, i, wl;
  char **keyWords = NULL, **words = NULL;
  const char routine[] = "fetchKeyList";
  int toUpper = 1;
  size_t nKeyAlloc;
  int nKeyChunk = 100;
  FILE *f = NULL;
  int isFirst = 1;
  int wasSlash = 0;
  int lastWasMulti = 0;
  char buffer[LINE_LENGTH+1];
  char b[LINE_LENGTH+1];
  size_t lKey;
  int iReturn;
  char padWord[TOKEN_LENGTH+3];
  
  nKeyWords = 0; *keyWordsP = NULL;
  nw = 0;

  if (!(f = fopen(fName, "r"))) {
    iReturn = -1;
    eclError(routine,iReturn, "Could not open file for read:\n\t%s\n",fName);    
    goto bottom;
  }

  nKeyAlloc = 0;
  /* Set up keyword array */
  keyWords = (char **) manageHeap((GENERIC_POINTER) NULL, (size_t) nKeyChunk, &nKeyAlloc, 
								  sizeof(char *), (size_t) nKeyChunk, 0);

  if (!keyWords) {
    iReturn = -2;
    keyWordsP = (char ***) NULL;
    eclError(routine,iReturn, "Could not initialize keyWords array.\n");    
    goto bottom;
  }
  
  while (fgets(buffer, (int) LINE_LENGTH, f)) { /* work thru file */
    
    (void) removeComment(buffer,b,commentString);
    
    if (strlen(b)<1) goto botLoop;
    
    nw = strtokSafe(b,ECL_TOKENS, &words, -1); /* get all words */

    if (nw > 0) {	/* found some words */

      /* found a terminating slash so we know this ain't a keyword, and */
      /* expect the next non-null buffer may be a keyword */
      if (words[nw-1][0]=='/') { 
		wasSlash = 1;            
		goto botLoop;            
      }

      if (!(wasSlash || isFirst)) goto botLoop; /* no reason to suspect a keyword */

      isFirst = 0; wasSlash = 0;

      /* handle single-word directives ... add to this def in "eclio.h" if needed */
      padWord[0] = ' ';  wl = strlen(words[0]);
      for (i=1;i<=wl;i++) padWord[i]=words[0][i-1];
      padWord[wl+1] = ' '; padWord[wl+2]='\0';
      if (strstr(ONE_LINE_RECORDS, padWord )) {
		wasSlash = 1;
      }

      nKeyWords++; /* this follows a slash or is first record, is non-null after
					  stripping comments --> probably a keyword */

      keyWords = (char **) manageHeap((GENERIC_POINTER) keyWords, 
									  (size_t) nKeyWords+1, &nKeyAlloc, 
									  sizeof(char *), (size_t) nKeyChunk, 0);
      if (!keyWords) {
		iReturn = -3;
		keyWordsP = (char ***) NULL;
		eclError(routine,iReturn, "Could not reallocate keyWords array.\n");    
		goto bottom;
      }

      keyWords[nKeyWords] = (char *) malloc((size_t)(strlen(words[0])+1) * sizeof(char));
      if (!keyWords[nKeyWords]) {
		iReturn = -4;
		eclError(routine,iReturn, "Could not initialize keyWord string %d\n:%s\n",
				 nKeyWords, words[0]);    
		keyWordsP = (char ***) NULL;
		goto bottom;
      }
      
      /* move keyword to list */
      strncpy (keyWords[nKeyWords],words[0],strlen(words[0])+1);
    }
   
  botLoop:
    ;
	if (words) 
	  {
		freeWordList(words, nw); // JG_MEMORY
		words=NULL;
	  }
  }

  iReturn = nKeyWords;

 bottom:
  if (f) fclose(f);
  if (words) freeWordList(words, nw); // JG_MEMORY
  keyWords = (char **) manageHeap((GENERIC_POINTER) keyWords, 
								  (size_t) nKeyWords+1, &nKeyAlloc, 
								  sizeof(char *), (size_t) nKeyChunk, 1);
  if (keyWords!=NULL && !keyWords[nKeyWords]) 
	{
	  iReturn = -4;
	  eclError(routine,iReturn, "Could not initialize keyWord string %d\n:%s\n",
			   nKeyWords, words[0]);    
	  keyWordsP = (char ***) NULL;
	  goto bottom;
	}
  else *keyWordsP = keyWords;

  return (iReturn);
}


int editKeyList (int nKeys, char ***keysP)
     /* Eliminates keywords that are not to be processed 
		returns -1 for a reallocation failure */
{
  int offset = 0;
  int i, j, nk;
  char **keys = NULL;
  char routine[] = "editKeyList"; 

  
  nk = nKeys; keys = *keysP;

  if (!keys) {
    nk = 0;
    eclError(routine,nk, "Null keys array,\n");
    goto bottom;
  }

  for (i=1;i<=nk;i++) {
    if (!keys[i]){
      eclError(routine,1,"Null key (# %d) ignored.\n");
    }
    else if (strstr(keys[i],"NOECHO") ||
	strstr(keys[i],"ECHO") ) {
      offset++;
      if (keys[i]) free (keys[i]);
      for (j=i+offset;j<=nk;j++) {
	keys[j-offset]=keys[j];
      }
      nk--;
    }
  }

  if (nk != nKeys) {
    keys = realloc(keys,(size_t) (nk + 1) * sizeof(char **));
    if (!keys) {
      nk = -1;
      eclError(routine, 2, "Could not realloc key array.\n");
    }
    *keysP = keys;
  }

 bottom:
  return (nk);
}

void checkKeyWord(char *key, char *weirdOK, char *msg, int lmsg)
     /* issues a status string for a keyword...additional checks easily appended */
{
  int i, l, lw;

  l = strlen(key); lw = strlen(weirdOK);

  for (i=0;i<l;i++) {
    if (!isalnum(key[i])) {
      if (lw) if (!strchr(weirdOK,key[i])) {
	strncpy(msg,"Contains nonalphanumeric/unexcepted characters",lmsg);
	return;
      }
    }
  }
  strncpy(msg,"",lmsg);
  return;
}

int trace4ToCornerpoints(ECL *eclSrc, ECL *eclDest, int update, int skipVoids)
     /* Moves a trace4 array (NEWCOOR4) from eclSrc to eclDest
	   Will create NEWCOOR4 if not already done
	   Reassigns size of array as needed
	   Creates/overwrites ZCORN and COORD in eclDest
	   eclDest should have x,y dimension equal to those of eclSrc
	      if SPECGRID is not set on call, it is set on return
	   eclSrc and eclDest can be identical -- but watch out on later calls!
        Returns 0 for success
     XXX*/
{
  PROP *zcDP=NULL, *cDP=NULL;
  PROP *t4SP=NULL;
  double ****t4S=NULL;
  double ****zcD=NULL, ***cD=NULL;
  double *c=NULL, *z=NULL, *t=NULL;
  int i, j, k, l;
  int nxS, nyS, nzS;
  int nxD, nyD, nzD;
  int ic, jc;
  char routine[] = "trace4ToCornerPoints";

  int status = 0;

  /* Check the SPECGRIDs */
  status = checkSpecGrid(eclSrc, &nxS, &nyS, &nzS);
  if (status) {
    status = 1;
    eclError(routine, status, "Problems with source specgrid.\n");
    goto bottom;
  }

  nxD = nxS; nyD = nyS; nzD = 1;
  copySpecGrid(eclDest, eclSrc);
  eclDest->specGrid.nz = 1; /* identical records except for nz */
  

  /* Get pointers to the arrays we need */
  /* trace4 or NEWCOOR4 in the Src ecl */
  t4SP = prepareTrace4(eclSrc, update, skipVoids);
  if (!t4SP) {
    status = 2;
    eclError(routine, status, "Couldn't get source TRACE4 property.\n");
    goto bottom;
  }

  t4S = (double ****) t4SP->p;
  if (!t4S) {
    status = 3;
    eclError(routine, status, "Couldn't get source TRACE4 data.\n");
    goto bottom;
  }

  /* Two results arrays, ZCORN and COORD, in the Dest ecl */
  zcDP = fetchProp(eclDest, "ZCORN");
  if (zcDP) deleteProp(eclDest, "ZCORN");
  zcDP = addProp(eclDest,"ZCORN", 'd',
		 4,1,nxD,1,nyD,1,1,1,8); 

  if (!zcDP){ 
    status = 4;
    eclError(routine, status, "Couldn't get dest ZCORN prop.\n");
    goto bottom;
  }

  zcD = (double ****) zcDP->p;
  if (!zcD) {
    status = 5;
    eclError(routine, status, "Couldn't get dest ZCORN data.\n");
    goto bottom;
  }

  cDP = fetchProp(eclDest, "COORD");
  if (cDP) deleteProp(eclDest, "COORD");
  cDP = addProp(eclDest,"COORD",'d',3,0,nxD,0,nyD,1,6);

  if (!cDP){
    status = 6;
    eclError(routine, status, "Couldn't get dest COORD prop.\n");
    goto bottom;
  }
  cD = (double ***) cDP->p;
  if (!cD) {
    status = 7;
    eclError(routine, status, "Couldn't get dest COORD data.\n");
    goto bottom;
  }

  /* With traces and dest arrays in hand, it's just a matter of permutation */
  
  /* Do COORD first */
  for (j=0;j<=nyD;j++) {
    if (j == 0) jc = 2;
    else jc = 0;
    for (i=0;i<=nxD;i++) {
      if (i == 0) ic = 1;
      else ic = 0;
      c = cD[i][j]; t = t4S[i][j][ic+jc+1];
      if (!c) {
	status = 8;
	eclError(routine, status, "Couldn't get dest COORD data for IJ = %d %d.\n",i,j);
	goto bottom;
      }
      if (!t) {
	status = 9;
	eclError(routine, status, "Couldn't get dest TRACE data for IJ = %d %d.\n",i,j);
	goto bottom;
      }
      /* all traces must be colinear; differing t/b don't matter */
      for (l=1;l<=6;l++) c[l]= t[l];
    }
  }
  
  /* ZCORN is a bit tricker; must get the corner from the right trace
     For gridblock IJ* (because K==1) we must get the 
        t4S[I-1][J-1][4][3,6] for the 1st and 5th corners
	t4S[I  ][J-1][3][3,6] for the 2nd and 6th corners
	t4S[I-1][J  ][2][3,6] for the 3rd and 7th corners, and
	t4S[I  ][J  ][1][3,6] for the 4th and 8th corners.
  */
  for (j=1;j<=nyD;j++){
    for (i=1;i<=nxD;i++) {

      z = zcD[i][j][1];
      if (!z) {status = 10; goto bottom;}

      t = t4S[i-1][j-1][4];
      if (!t) {status = 11; goto bottom;}
      z[1]=t[3]; z[5]=t[6];
      
      t = t4S[i][j-1][3];
      if (!t) {status = 12; goto bottom;}
      z[2]=t[3]; z[6]=t[6];
      
      t = t4S[i-1][j][2];
      if (!t) {status = 13; goto bottom;}
      z[3]=t[3]; z[7]=t[6];

      t = t4S[i][j][1];
      if (!t) {status = 14; goto bottom;}
      z[4]=t[3]; z[8]=t[6];
    }
  }

 bottom:
  if (status >= 10) {
    eclError(routine, status, "Couldn't get trace or zCorn data in loop.\n");
  }

  return (0);
}
int copyPropECL(ECL *eclDest, char *nameDest, ECL *eclSrc, char *nameSrc,
		int asGridProp)
{
  /* Copies to  PROP storage from fromX, resetting PROP attributes as needed.
     if (toProp) on entry, ecl can be NULL */
  PROP *from, *to;

  int nD, nDest;
  int **bndDest = NULL;
  int iReturn = 0;
  char routine[] = "copyPropECL";
  
  from = fetchProp(eclSrc, nameSrc);
  if (!from) {
    iReturn = 1;
    eclError(routine, iReturn, "Couldn't get prop %s from %s.\n",nameSrc,eclSrc->name);
    goto bottom;
  }
  
  nD = from->nInd;

  if (asGridProp) {
    bndDest = buildBndGrid(eclDest);
    nDest = 3;
    if (!bndDest) {
      iReturn = 2;
      eclError(routine, iReturn, 
	       "Could not build bounds array\nSource:\t%s\t%s\n",
	       nameSrc, nameDest);
      goto bottom;
    }
    iReturn = checkArrayBounds(nDest, bndDest, nD, from->bnd,1);
    if (iReturn) {
      iReturn = 3;
      eclError(routine, iReturn, "Array bounds don't check\n\t%s\t%s\n",
	       nameSrc, nameDest);
      goto bottom;
    }
  }
  else {
    bndDest = copyBnd(nD, from->bnd);
    if (!bndDest) {
      iReturn = 4;
      eclError(routine, iReturn, "Could not copy bounds.\n");
      goto bottom;
    }
    nDest = nD;
  }

  if (1<= nD && nD <=4) {

    to = fetchProp(eclDest, nameDest);

    if (to) { /* get rid of old storage and reset description and pointer to store */
      freeAnyArray(to->p, to->ptype, to->nInd, to->bnd);
      setProp(to, nameDest, from->ptype, nDest, bndDest, (GENERIC_POINTER) NULL);
      to->p = anyProp(to->ptype, nDest, bndDest, nameDest);
    }

    else {    /*Add array-ful PROP .. note the storage at bndDest is used*/
      to = addPropN(eclDest, nameDest, from->ptype, nDest, bndDest); 
    } 

    if (!to || !to->p) {
      iReturn = 5;
      eclError(routine, iReturn, "Could not allocate dest data area\n");
      goto bottom;
    }
  }
  else goto bottom;

  if (!copyProp(to,from)) to->set = 1;
  else {
    iReturn = 6;
    eclError (routine, iReturn, "copyProp failed\n");
    goto bottom;
  }

 bottom:
  if (iReturn && bndDest) free_imatrix(bndDest,1,nDest,1,2);

#if 0 && JG_FEATURES
  if (!strcmp(nameDest,"ZCORN"))
	{
	  eclDest->zCornp = to->p;
	  eclDest->zCorn = (double ****) eclDest->zCornp->p;
	}
#endif 


  return (0);
}
int assignArrayToProp(ECL *ecl, PROP *toProp, char *name,  
		      GENERIC_POINTER fromX,
		      int nInd, int **bnd, char ptype,
		      int asGridProp)
{
  /* Just points a PROP storage at fromX, resetting PROP attributes as needed.
     if (toProp) on entry, ecl can be NULL */
  int i,j;
  int status = 0;
  int nDest;
  int **bndDest = NULL;
  char routine[] = "assignArrayToProp";

  if (asGridProp) {

    bndDest = buildBndGrid(ecl);
    nDest = 3;
    if (!bndDest) {
      status = 1;
      eclError(routine, status, "Could not build bounds\n");
      goto bottom;
    }

    /* 3 = strict check ... identical rank, dimensions and offset */
    status = checkArrayBounds(nDest, bndDest, nInd, bnd, 3); 
    if (status) {
      status = 2;
      eclError(routine, status, "Failed bounds check.\n");
      goto bottom;
    }
  }

  else {
    bndDest = copyBnd(nInd,bnd);
    if (!bndDest) {
      status = 3;
      eclError(routine, status, "Could not copy bounds.\n");
      goto bottom;
    }
    nDest = nInd;
  }

  if (!fromX) {
    status = 4;
    eclError(routine, status, "Source array is NULL..will prodeed\n");
  }

  if (nInd < 1 || nInd > 4) {
    status = 5;
    eclError(routine, status, "Bad ranks, not [1,4]\n");
    goto bottom;
  }

  if (!(ptype == 'i' || ptype == 'f' || ptype == 'd')) {
    status = 6;
    eclError(routine, status, "Bad type, not [i,f,d].\n");
    goto bottom;
  }
  
  if (toProp) { /* get rid of old contents and reinitialize if toProp */
    freeAnyArray(toProp->p, toProp->ptype, toProp->nInd, toProp->bnd);
    setProp(toProp, name, ptype, nDest, bndDest, (GENERIC_POINTER) NULL);
  }
  else {
    toProp = addProp0(ecl, name, ptype, nDest, bndDest); /*Add array-less PROP */
  }

  if (!toProp) {
    status = 7;
    eclError(routine, status, "Could not add PROP.\n");
    goto bottom;
  }

  toProp->p = fromX;
  toProp->set = 1;
  
 bottom:

  if (status && status != 4 && bndDest) free_imatrix(bndDest,1,nDest, 1, 2);

  return (status);
}	   


int copyArrayToProp(ECL *ecl, PROP *toProp, char *name,
					GENERIC_POINTER fromX,
					int nInd, int **bnd, char ptype,
					int asGridProp)
{
  int i,j;
  int status = 0;
  int nDest;
  int **bndDest = NULL;
  char routine[] = "copyArrayToProp";

  if (asGridProp) {
    bndDest = buildBndGrid(ecl);
    nDest = 3;
    if (!bndDest) {
      status = 1;
      eclError(routine, status, "Could not build bounds.\n");
      goto bottom;
    }

    /* 1 = loose check, just need consistent data area */
    status = checkArrayBounds(nDest, bndDest, nInd, bnd,1);
    if (status) {
      status = 2;
      eclError(routine, status, "Failed bounds check.\n");
      goto bottom;
    }
  }

  else {
    bndDest = copyBnd(nInd, bnd);
    if (!bndDest) {
      status = 3;
      eclError(routine, status, "Could not copy bounds.\n");
      goto bottom;
    }
    nDest = nInd;
  }

  if (!fromX) {
    status = 4;
    eclError(routine, status, "Null array to copy!\n");
    goto bottom;
  }

  if (nInd < 1 || nInd > 4) {
    status = 5;
    eclError(routine, status, "Bad rank, not [1,4]\n");
    goto bottom;
  }

  if (!(ptype == 'i' || ptype == 'f' || ptype == 'd')) {
    status = 6;
    eclError(routine, status, "Bad type, not [i,f,p]\n");
    goto bottom;
  }

  if (toProp) { /* get rid of old contents and reinitialize if toProp */
    freeAnyArray(toProp->p, toProp->ptype, toProp->nInd, toProp->bnd);
    setProp(toProp, name, ptype, nDest, bndDest, (GENERIC_POINTER) NULL);
    toProp->p = anyProp(ptype, nDest, bndDest, name);
  }
  else {
    toProp = addPropN(ecl, name, ptype, nDest, bndDest);  /*Add array-ful PROP */
  }

  if (!toProp || !toProp->p) {
    status = 7;
    eclError(routine, status, "Could not get source prop.\n");
    goto bottom;
  }

  if (!copyArray(toProp->p,toProp->nInd, toProp->bnd,
		 fromX, nInd, bnd, 
		 ptype)) toProp->set = 1;
  
  else {
    toProp->set = 0;
    status = 8;
    eclError(routine, status, "copyArray failed.\n");
  }

 bottom:
  if (status && bndDest) free_imatrix(bndDest,1,nDest,1,2);
  return (status);
}     		   
int writeRecordList (FILE *f, ECL *ecl, int nKeys, char **keyList)
{
  int i;
  int status = 0;
  char routine[] = "writeRecordList";
  
  if (!ecl) {
    status = 1;
    eclError(routine, status, "NULL ecl\n");
  }

  if (nKeys < 1) {
    status = 2;
    eclError(routine, status, "Nonpositive nKeys\n");
  }

  if (!f) {
    status = 3;
    eclError(routine, status, "NULL output file pointer\n");
  }
  
  if (!keyList) {
    status = 4;
    eclError(routine, status, "NULL keyList\n");
  }
  
  if (status) goto bottom;

  for (i=1;i<=nKeys;i++) {
    if (!keyList[i]) {
      status = 5;
      eclError(routine, status, "NULL key, keyword # %d ... skipping\n", i);
    }
    else {
      status = wRecord(f, ecl, keyList[i]);
      if (status) {
		status = 6;
		eclError(routine, status, "Error writing keyWord: %s ... continuing\n",keyList[i]); // JG_FEATURES
      }
    }
  }
 bottom:
  return (status);
}

int squashZCorn(ECL *eclDest, ECL *eclSrc)

{
  int i, j, k, p, q, r;
  double ****zS=NULL, ****zD= NULL;
  int ***aS = NULL;
  PROP *zSP = NULL;
  PROP *zDP = NULL;
  PROP *aSP = NULL;

  double *zT = NULL, *zB = NULL, *zN = NULL;

  int nxS, nyS, nzS;
  int nxD, nyD, nzD;
  int kTop, kBot;
  char routine[] = "squashZCorn";

  int status = 0;

  /* Check the SPECGRIDs */
  status = checkSpecGrid(eclSrc, &nxS, &nyS, &nzS);
  if (status) {
    status = 1;
    eclError(routine, status, "Failed specgrid check on src ecl.\n");
    goto bottom;
  }

  status = checkSpecGrid(eclSrc, &nxD, &nyD, &nzD);
  if (status) {
    status = 2;
    eclError(routine, status, "Failed specgrid check on dest ecl.\n");
    goto bottom;
  }


  /* Get pointers to the arrays we need */
  /* trace4 or NEWCOOR4 in the Src ecl */
  zSP = fetchProp(eclSrc, "ZCORN");
  if (!zSP || !zSP->p) {
    status = 3;
    eclError(routine, status, "Could not get ZCORN from src ecl\n");
    goto bottom;
  }
  zS = (double ****) zSP->p;

  aSP = fetchProp(eclSrc,"ACTNUM");
  if (!aSP || !aSP->p) {
    status = 4;
    eclError(routine, status, "Could not get ACTNUM from src ecl\n");
    goto bottom;
  }
  aS = (int ***) aSP->p;

  /*  results arrays, ZCORN in the Dest ecl */
  zDP = fetchProp(eclDest, "ZCORN");
  if (zDP) deleteProp(eclDest, "ZCORN");
  zDP = addProp(eclDest,"ZCORN", 'd',4,1,nxD,1,nyD,1,1,1,8); 

  if (!zDP || !zDP->p) {
    status = 5;
    eclError(routine, status, "Could not fetch or set ZCORN in dest ecl\n");
    goto bottom;
  }
  zD = (double ****) zDP->p;

  /* With traces and dest arrays in hand, it's just a matter of permutation */
  for (i=1;i<=nxD;i++) {
    for (j=1;j<=nyD;j++) {
      kTop = I_MISSING; kBot = I_MISSING;
      zN = zD[i][j][1]; 
      if (!zN) return (-51);
      for (k=1;k<=nzS;k++) {
	if (aS[i][j][k]) {
	  kTop = k;
	  break;
	}
      }
      if (kTop != I_MISSING) {
	for (k=nzS;k>=1;k--) {
	  if (aS[i][j][k]) {
	    kBot = k;
	    break;
	  }
	}
      }
      if (kTop == I_MISSING) 
	for (p=1;p<=8;p++) zN[p] = D_MISSING;
      else {
	zT = zS[i][j][kTop]; 
	if (!zT) return (-52);
	zB = zS[i][j][kBot];
	if (!zB) return (-53);
	for (p=1;p<=4;p++) {
	  zN[p] = zT[p];
	  zN[p+4] = zB[p+4];
	}
      }
    }
  }
 bottom:
  return (status);
}
 
