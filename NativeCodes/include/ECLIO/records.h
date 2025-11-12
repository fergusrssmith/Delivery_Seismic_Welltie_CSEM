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
#ifndef _ecl_records_h_
#define _ecl_records_h_
int getNSlashReqd(char *key,  ECL *ecl);
int cleanRecord(char *buffer);
int seekToken1(char *buffer, char *key, int toUpper);
int scrubRecord(char *buffer, char *elimChars, char replChar);
int makeTokenPristine(char *src, char *dest, size_t maxDestSize, 
		      char *separators, int forceUp);
int seekToken1(char *buffer, char *key, int toUpper);
int findEndOfRecord( ECL *ecl, char *key, FILE *f);
int findSlash(FILE *f, char *commentString);
int removeTrailingComment(char *b1, char *b2, size_t lineLen, char *commentString);
int removeComment (char *b1, char *b2, char *commentString);
int strtokSafe(char *buffer, char *tokens, char ***wordp, int nToGet);
int strtokSafe2(char *buffer, char *tokens, char ***wordp, int nToGet);
void freeWordList(char ** words, int nWords);
#endif

