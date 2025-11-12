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
#ifndef _ecl_blockops_h_
#define _ecl_blockops_h_
int topAct(ECL *ecl, int i, int j);
int botAct(ECL *ecl, int i, int j);
double zTopBlock(ECL *ecl, int i, int j, int ic, int jc, 
		 int skipVoids);
double zBotBlock(ECL *ecl, int i, int j, int ic, int jc, 
		 int skipVoids);
int zBarBlock(ECL *ecl, int skipVoids, 
	      int iBlk, int jBlk, int ic, int jc,
	      double *ztbar, int *ntbar, 
	      double *zbbar, int *nbbar);
void interpVert(double z, double *point, double *c);
int rescaleGridMap (ECL *ecl);
int rotateGridMap(ECL *ecl);
int rotateCoord(ECL *ecl, double ***jacP, double *det);
int translateCoord(ECL *ecl);
int applyPinch(ECL *ecl);    
PROP *prepareTraces(ECL *ecl, int forceCalc, int skipVoids);
PROP *prepareActCoord(ECL *ecl, int forceCalc);
PROP *prepareActCol(ECL *ecl, int forceCalc);
PROP *prepareSegCol(ECL *ecl, int forceCalc);
int countActive(ECL *ecl);
int countActCol(ECL *ecl, int forceCalc);
int countActCoord(ECL *ecl, int forceCalc);
void freeRotateMem(double **jac, double *tmp);
PROP *prepareTrace4(ECL *ecl, int forceCalc, int skipVoids);
int oneBlockEdge(ECL *ecl, int nx, int ny, int i, int j, int it, double ***c,
		 double ****t);

#if JG_FEATURES
int zBarNumberedBlock(ECL *ecl, int kBlock, int skipVoids,
					  int iBlk, int jBlk, int ic, int jc,
					  double *ztbar, int *ntbar, 
					  double *zbbar, int *nbbar);
double zNumberedBlock(ECL *ecl, int i, int j, int ic, int jc, int kBlock,int skipVoids);
PROP *prepareTracesAndVolumes(ECL *ecl, int forceCalc, int skipVoids, int layerNumberForTraces);
double midPointTraceCentredArea(ECL *ecl,double ***t,int **active,int i,int j);
int getMidPoint(ECL *ecl,double ***t,int** active, int i,int j,double *x,double *y);
int flipGridAndMapUnits(ECL *ecl);
#endif

#endif
