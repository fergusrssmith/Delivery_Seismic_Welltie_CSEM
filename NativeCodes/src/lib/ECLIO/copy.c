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
#include "ECLIO/copy.h" 
#include "NR/nrutil.h"

/*************************************************************************/
int copySpecGrid(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  SPECGRID *d;
  SPECGRID *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->specGrid);
  s = &(sEcl->specGrid);
		
  if (!d || !s || !s->set) return (1);

  d->nx = s->nx;
  d->ny = s->ny;
  d->nz = s->nz;
  
  d->numRes = s->numRes;
  d->cylindrical = s->cylindrical;
 
  d->set = 1;

  return (0);
}
/*************************************************************************/
 int copyPinch(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  PINCH *d;
  PINCH *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->pinch);
  s = &(sEcl->pinch);

  if (!d || !s || !s->set) return (1);

  d->pVal = s->pVal;
  strcpy(d->gap, s->gap);
  d->pinchPV = s->pinchPV;
  strncpy(d->topbot, s->topbot, TOKEN_LENGTH+1);
  strncpy(d->tMeth, s->tMeth, TOKEN_LENGTH+1);

  d->set = 1;

  return (0);
}
/*************************************************************************/
int copyCoordSys(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  int set;
  int i;
  int numRes;

  COORDSYS *d;
  COORDSYS *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->coordSys);
  s = &(sEcl->coordSys);

  if (!d || !s || !s->set || s->numRes < 1) return (1);
  
  numRes = s->numRes;
  d->numRes = numRes;

  d->k1 = ivector(1,numRes);
  d->k2 = ivector(1,numRes);

  d->comp = cmatrix(1,numRes,0,TOKEN_LENGTH);
  d->join = cmatrix(1,numRes,0,TOKEN_LENGTH);

  if (!d->k1 || !d->k2 || !d->comp || !d->join) return (1);

  for (i=1;i<=numRes;i++) {
    d->k1[i] = s->k1[i];
    d->k2[i] = s->k2[i];
    if (d->k2[i] > dEcl->specGrid.nz) d->k2[i] = dEcl->specGrid.nz;
    strncpy(d->comp[i],s->comp[i], TOKEN_LENGTH);
    strncpy(d->join[i],s->comp[i], TOKEN_LENGTH);
  }

  d->set = 1;

  return (0);
}
/*************************************************************************/
int copyGridUnit(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  GRIDUNIT *d;
  GRIDUNIT *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->gridUnit);
  s = &(sEcl->gridUnit);

  if (!d || !s || !s->set) return (1);

  strncpy(d->units, s->units, TOKEN_LENGTH);
  strncpy(d->map, s->map, TOKEN_LENGTH);

  d->set = 1;
  
  return (0);
}
/*************************************************************************/
int copyMapUnits(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  MAPUNITS *d;
  MAPUNITS *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->mapUnits);
  s = &(sEcl->mapUnits);

  if (!d || !s || !s->set) return (1);

  strncpy(d->units, s->units, TOKEN_LENGTH);
  d->set = 1;

  return (0);
}
/*************************************************************************/
int copyMapAxes(ECL *dEcl, ECL *sEcl)
/*************************************************************************/
{
  int i, j;
  MAPAXES *d;
  MAPAXES *s;

  if (!dEcl || !sEcl) return (1);
  d = &(dEcl->mapAxes);
  s = &(sEcl->mapAxes);

  if (!d || !s || !s->set || !s->axes) return (1);

  d->axes = dmatrix(1,3,1,2);
  if (!d->axes) return (1);

  for (i=1;i<=3;i++){
    for (j=1;j<=2;j++) d->axes[i][j] = s->axes[i][j];
  }
  d->isReversed = s->isReversed;
  d->set = 1;

  return (0);  
}
