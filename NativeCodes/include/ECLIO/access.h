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
#ifndef _ecl_access_h_
#define _ecl_access_h_
int loadModel(char *fName, ECL *ecl, int transformToMap, int nKeys, ...);
int loadModelFromList(char *fName, ECL *ecl, int transformToMap,
		      int nKeys, ...);
int loadModelFromArray(char *fName, ECL *ecl, int transformToMap,
		       int nKeys, char **keys);
int loadModelFromArraySansPinch(char *fName, ECL *ecl, int transformToMap,
								int nKeys, char **keys);
void freeModel(ECL *ecl);
void freeLoadModelMem(int nKeys, size_t l, char **keys);
PROP *fetchProp(ECL *ecl, char *name);
GENERIC_POINTER fetchNamedRecord(ECL *ecl, char *name);
GENERIC_POINTER fetchAnyRecord(ECL *ecl, char *name, char *rtype);
GENERIC_POINTER fetchPropArray(ECL *ecl, char *name);
double ***fetchTraces(ECL *ecl, int skipVoids,
		      int nrl, int nrh, int ncl, int nch);
int **fetchActCoord(ECL *ecl, int nrl, int nrh, int ncl, int nch);
int getNX(ECL *ecl);
int getNY(ECL *ecl);
int getNZ(ECL *ecl);
PROP *addProp(ECL *ecl, char *name, char ptype, 
		     int nInd, ...);
int deleteProp(ECL *ecl, char *name);
int freeProp(ECL *ecl, char *name);
PROP *addPropGrid(ECL *ecl, char *name, char ptype);
int ***fetchActnum(ECL *ecl, 
		   int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
int **fetchActCol(ECL *ecl, int nrl, int nrh, int ncl, int nch);
int ***fetchTopBot(ECL *ecl, 
		   int nrl, int nrh, int ncl, int nch);
int ***fetchSegments(ECL *ecl, 
		     int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
int **fetchSegCol(ECL *ecl, int nrl, int nrh, int ncl, int nch);
double ****fetchTrace4(ECL *ecl, int skipVoids,
		       int nrl, int nrh, int ncl, int nch);
char contFlags(double *t1, double *t2);
int fetchKeyList(char *fName, char ***keyWordsP, char *commentString);
int editKeyList(int nKeys, char ***keysP);
void checkKeyWord(char *key, char *weirdOK, char *msg, int lmsg);
int trace4ToCornerpoints(ECL *eclSrc, ECL *eclDest, int update, int skipVoids);
int copyPropECL(ECL *eclDest, char *nameDest, ECL *eclSrc, char *nameSrc, 
		int asGridProp);
int assignArrayToProp(ECL *ecl, PROP *toProp, char *name, 
		      GENERIC_POINTER fromX,
		      int nInd, int **bnd, char ptype,
		      int asGridProp);
int copyArrayToProp(ECL *ecl, PROP *toProp, char *name,
		    GENERIC_POINTER fromX,
		    int nInd, int **bnd, char ptype,
		    int asGridProp);
int writeRecordList (FILE *f, ECL *ecl, int nKeys, char **keyList);
int squashZCorn(ECL *eclDest, ECL *eclSrc);

#if JG_FEATURES
double **fetchTraceAreas(ECL *ecl, int nrl, int nrh, int ncl, int nch);
double ***fetchLayerTraces(ECL *ecl, int eclLayerNumber, int skipVoids,
						   int nrl, int nrh, int ncl, int nch);
#endif 

#endif

