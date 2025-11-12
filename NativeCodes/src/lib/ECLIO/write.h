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
#ifndef _ecl_write_h_
#define _ecl_write_h_

// JG_FEATURES: one more digit precision on floats
#define FLOAT_FORMAT   " %10.3f"
#define DOUBLE_FORMAT  " %15.6lf"
#define INT_FORMAT     " %7d"   
#define FLOAT_ON_LINE 8
#define DOUBLE_ON_LINE 5
#define INT_ON_LINE 10
/* Handle record termination options here ... 
   set RECORD_END_FORMAT as desired and recompile */
#define RECORD_END_FORMAT_NO_NEW " /\n"
#define RECORD_END_FORMAT_NEW " \n /\n"
// #define RECORD_END_FORMAT RECORD_END_FORMAT_NO_NEW
// JG_FEATURES
#define RECORD_END_FORMAT RECORD_END_FORMAT_NEW
 
int wRecord(FILE *f, ECL *ecl, char *name);
int wProp(FILE *f, ECL *ecl, char *name);
int wSpecGrid(FILE * f, ECL *ecl);
int wCoordSys(FILE * f, ECL *ecl);
int wNumRes(FILE * f, ECL *ecl);
int wDimens(FILE * f, ECL *ecl);
int wPinch(FILE *f, ECL *ecl);
int wPropHeader(FILE *f, ECL *ecl, char *name);
int wCoordSys(FILE *f, ECL *ecl);
int wGridUnit(FILE *f, ECL *ecl);
int wMapUnits(FILE *f, ECL *ecl);
int wMapAxes(FILE *f, ECL *ecl);
int wZCorn(FILE *f, ECL *ecl);
int wCoord(FILE *f, ECL *ecl, char *name);
int wPropArrayI(FILE *f, GENERIC_POINTER x, int nD, int **b);
int wPropArrayI1(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayI2(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayI3(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayI4(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayF(FILE *f, GENERIC_POINTER x, int nD, int **b);
int wPropArrayF1(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayF2(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayF3(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayF4(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayD(FILE *f, GENERIC_POINTER x, int nD, int **b);
int wPropArrayD1(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayD2(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayD3(FILE *f,GENERIC_POINTER x, int **b);
int wPropArrayD4(FILE *f,GENERIC_POINTER x, int **b);
#endif
