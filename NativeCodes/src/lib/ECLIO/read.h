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
#ifndef _ecl_read_h_
#define _ecl_read_h_
int findKey(ECL *ecl, FILE *f, char *keyIn, int toUpper, int rewindF);
void freeFindKeyMem(char *key, char* token, char *buffer);
int findNextKey(ECL *ecl, FILE *f, int toUpper, int nKeys, char **keysIn);
void freeFindNextKeyMem(char *token, int nKeys, char **key, char *buffer, 
			char **words, int nw);
int r_specGrid(FILE * f, ECL *ecl);
int r_pinch(FILE *f, ECL *ecl);
int r_coord(FILE *f, ECL *ecl);
int r_zCorn(FILE *f, ECL *ecl);
int r_prop(FILE *f,  ECL *ecl, char *name, char ptype,
	   int nInd, ...);
int r_propGrid(FILE *f, struct ECL *ecl, char *name, char ptype);
int r_prop1(FILE *f, char *name, GENERIC_POINTER p, char ptype, int nInd, int **bnd);
int r_prop2(FILE *f, char *name, GENERIC_POINTER  p, char ptype, int nInd, int **bnd);
int r_prop3(FILE *f, char *name, GENERIC_POINTER  p, char ptype, int nInd, int **bnd);
int r_prop4(FILE *f, char *name, GENERIC_POINTER  p, char ptype, int nInd, int **bnd);
int r_coordSys(FILE *f, struct ECL *ecl);
int r_gridUnit(FILE *f, struct ECL *ecl);
int r_mapUnits(FILE *f, struct ECL *ecl);
int r_mapAxes(FILE *f, struct ECL *ecl);
#endif
