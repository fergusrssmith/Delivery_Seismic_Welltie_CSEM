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
#include <ctype.h> 
#include "eclio.h"
#include "ECLIO/set.h"
#include "ECLIO/access.h"
#include "ECLIO/eclutil.h"
#include "NR/nrutil.h"

#define PINCH_DEFAULT 0.001
#define GAP_DEFAULT "GAP"
#define PINCHPV_DEFAULT 1.0E20
#define TOPBOT_DEFAULT "TOPBOT"
#define TMETH_DEFAULT "TOP"

#define DEFAULT_L_UNITS "METRES"


/*************************************************************************/
int s_specGrid(ECL *ecl, int nx, int ny, int nz, int numRes, 
	       char cylindrical)
/*************************************************************************/
{
  SPECGRID *sg = NULL;
  const char routine[] = "s_specGrid";

  if (nx < 1 || ny < 1 || nz < 1) {
    eclError(routine,1, 
	     "Must have all dimensions > 1, but (nx,ny,nz) = (%d, %d, %d)\n",
	     nx, ny, nz);
    return (1);
  }
  if (numRes < 1 || numRes > nz) {
    eclError(routine,2, 
	     "Must have 1<=numRes<=nz, but set to %d, nz= %d\nReset to 1\n", 
	     numRes, nz);
    numRes = 1;
  }
  if (cylindrical != 'T' && cylindrical != 'F') {
    eclError(routine,2, 
	     "Must have cylindrical set to T or F but set to %c\nSet to F", 
	     cylindrical);
    cylindrical = 'F';
  }

  sg = &ecl->specGrid;

  sg->nx = nx;
  sg->ny = ny;
  sg->nz = nz;
  sg->numRes = numRes;
  sg->cylindrical = cylindrical;
  
  sg->set = 1;

  return (0);
}
/*************************************************************************/
 int s_pinch (ECL *ecl, double pVal, char *gap, 
	      double pinchPV, char *topbot, char *tMeth)
/*************************************************************************/
{
  PINCH *p = NULL;
  const char routine[] = "s_pinch";
  
  p = &ecl->pinch;

  if (pVal < 0.0) pVal = PINCH_DEFAULT;
  p->pVal = pVal;

  if (strcmp(gap, "GAP") && strcmp(gap,"NOGAP")){
    eclError(routine,1, "Cannot parse GAP = %s\n",gap);
    return (1);
  }

  strncpy(p->gap, gap, TOKEN_LENGTH);

  if (pinchPV <0.0) pinchPV = PINCHPV_DEFAULT;
  p->pinchPV = pinchPV;

  if (strcmp(topbot, "TOPBOT") && strcmp(topbot, "ALL")) {
    eclError(routine,2, "Cannot parse TOPBOT = %s\n",topbot);
    return (2);
  }

  strncpy(p->topbot, topbot, TOKEN_LENGTH);

  if (strcmp(tMeth, "TOP") && strcmp(tMeth, "ALL")) {
    eclError(routine,3, "Cannot parse TMETH = %s\n",tMeth);
    return (3);
  }

  strncpy(p->tMeth, tMeth, TOKEN_LENGTH);

  p->set = 1;

  return (0);
}
/*************************************************************************/
int s_coordSys(ECL *ecl, int numRes, int *k1, int *k2, 
	       char **comp, char ** join)
/*************************************************************************/
{
  COORDSYS *cs = NULL;
  
  cs = &ecl->coordSys;

  cs->numRes = numRes;

  cs->k1 = k1;
  cs->k2 = k2;
  cs->join = join;
  cs->comp = comp;
  cs->set = 1;

  return (0);
}
/*************************************************************************/
int s_gridUnit(ECL *ecl, char *units, char *map)
/*************************************************************************/
{
  GRIDUNIT *gu = NULL;
  
  gu = &ecl->gridUnit;
  strncpy(gu->units, units, TOKEN_LENGTH);
  strncpy(gu->map, map, TOKEN_LENGTH);
  gu->set = 1;
  return (0);
}
/*************************************************************************/
int s_mapUnits(ECL *ecl, char *units)
/*************************************************************************/
{
  MAPUNITS *mu = NULL;
  
  mu = &ecl->mapUnits;
  strncpy(mu->units, units, TOKEN_LENGTH);
  mu->set = 1;
  return (0);
}
/*************************************************************************/
int s_mapAxes(ECL *ecl, double **axes)
/*************************************************************************/
{
  MAPAXES *ma = NULL;
  
  ma = &ecl->mapAxes;
  ma->axes = axes;
  ma->set = 1;
  ma->isReversed = 0;
  return (0);  
}

