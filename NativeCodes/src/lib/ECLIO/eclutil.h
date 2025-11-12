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
#ifndef _eclutil_h_
#define _eclutil_h_
#include <stdarg.h>
int initModel(ECL *e);
PROP *initPropVec();
void initProp(PROP *p);
PROP *addPropN(ECL *ecl, char *name, char ptype, 
		      int nInd, int **bnd);
PROP *addProp0(ECL *ecl, char *name, char ptype, 
	       int nInd, int **bnd);
PROP *allocProp(ECL *ecl, char *name);
void setProp(PROP *pp, char *name, char ptype, 
	     int nInd, int **bnd, GENERIC_POINTER p);
int **buildBnd(ECL *ecl, int nInd, va_list ap, char *name);
int **copyBnd (int n, int ** b);
int **buildBndGrid(ECL *ecl);
GENERIC_POINTER anyProp(char ptype, int nInd, int **bnd, char *name);
GENERIC_POINTER any1Prop(char ptype, int nInd, int **bnd, char *name);
GENERIC_POINTER any2Prop(char ptype, int nInd, int **bnd, char *name);
GENERIC_POINTER any3Prop(char ptype, int nInd, int **bnd, char *name);
GENERIC_POINTER any4Prop(char ptype, int nInd, int **bnd, char *name);
int eclError(const char *caller, int code, char *fmt, ...);
int reportError(int doExit, const char *caller, int code, char *fmt, ...);
void eclExit(ECL *ecl, int status);
void free_anyProp(PROP *p);
void free_any1Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype);
void free_any2Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype);
void free_any3Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype);
void free_any4Prop(int nInd, int **bnd, GENERIC_POINTER p, char ptype);
int transProp (PROP *pp, double factor, double transl);
int transIProp(PROP *pp, double f, double a);
int transFProp(PROP *pp, double f, double a);
int transDProp(PROP *pp, double f, double a);
int transI1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transI2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transI3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transI4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transF1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transF2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transF3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transF4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transD1Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transD2Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transD3Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
int transD4Prop(int n, int **bnd, GENERIC_POINTER p, double f, double a);
double ddot_ECL(int nrl1, int nrh1, double *v1, 
	    int nrl2, int nrh2, double *v2);
double *dgetCol(int getCol, int nrl, int nrh, int ncl, int nch, double **m);
double *dgetRow(int getRow, int nrl, int nrh, int ncl, int nch, double **m);
double **dmmMult(int nrl1, int nrh1, int ncl1, int nch1, double **m1,
		 int nrl2, int nrh2, int ncl2, int nch2, double **m2);
int dmvMult(int nrl, int nrh, int ncl, int nch, double **m,
	    int nl, int nh, double *v, double *d);
double det3mat(int nrl, int nrh, int ncl, int nch, double **m);
int transD2(int nrl, int nrh, int ncl, int nch, double **v,
	    double f, double a);
int fetchPropInd(ECL *ecl, char *name);
int compressProp(ECL *ecl);
int setAuxPointers(ECL *ecl);
int checkSpecGrid(ECL *ecl, int *nxp, int *nyp, int *nzp);
int calcNinArray(int nInd, int **bnd);
void freeCoordSys(int *k1, int *k2, char **comp, char **join, int n, int l);
double trilinearInterp(double *p, double *u);
double trilinearNull(double *p, double *u);
GENERIC_POINTER manageHeap(GENERIC_POINTER heap, 
			   size_t need, size_t *haveP, size_t sizeEach, size_t chunk,
			   int condense);
int copyProp(PROP *to, PROP *from);
int copyArray (GENERIC_POINTER t, int nT, int **bndT, 
	       GENERIC_POINTER f, int nF, int **bndF,
	       char pt);
GENERIC_POINTER locationOne(GENERIC_POINTER x, int nb, int **b, char pt);
int copyArrayD(GENERIC_POINTER tG, GENERIC_POINTER fG, long n);
int copyArrayF(GENERIC_POINTER tG, GENERIC_POINTER fG, long n);
int copyArrayI(GENERIC_POINTER tG, GENERIC_POINTER fG, long n);
long getNTotal (int nD, int **bnd);
GENERIC_POINTER allocAnyArray(char ptype, int nInd, int **bnd);
GENERIC_POINTER allocAnyArray1(char ptype, int nInd, int **bnd);
GENERIC_POINTER allocAnyArray2(char ptype, int nInd, int **bnd);
GENERIC_POINTER allocAnyArray3(char ptype, int nInd, int **bnd);
GENERIC_POINTER allocAnyArray4(char ptype, int nInd, int **bnd);
int freeAnyArray(GENERIC_POINTER x, char ptype, int nInd, int **bnd);
int checkBounds(int nNeed, int **bNeed, int nHave, int **bHave);
int checkBoundsStrict(int nNeed, int **bNeed, int nHave, int **bHave);
int checkBoundsVeryStrict(int nNeed, int **bNeed, int nHave, int **bHave);
int checkArrayBounds(int nNeed, int **bNeed, int nHave, int **bHave, int levels);
#endif
