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
int buildECLOutList(int nKeyWords, char **keyWords, 
		    int nHeaderWords, char **headerWords,
		    int *nOutWordsP, char **outWords)
{
  int *isHeader = NULL, *isKey = NULL;
  int nOutWords = *nOutWordsP;
  int i,j, li;
  char *tWord;

  if (!(isHeader = (int) malloc ((size_t) (nKeyWords + 1) * sizeof(int))) ||
      !(isKey = (int) malloc ((size_t) (nHeaderWords + 1) * sizeof(int))))
    return (1);

  for (i=1;i<=nKeyWords;i++) {
    tWord = keyWords[i];
    li = strlen(tWord);
    isHeader[i] = 0;
    for (j=1;j<=nHeaderWords; j++) {
      if (!strncmp(tWord,headerWord[j],li)) {
	isHeader[i] = 1;
	j = nHeaderWords +1;
      }
    }
  }

  for (i=1;i<=nHeaderWords;i++) {
    tWord = headerWords[i];
    li = strlen(tWord);
    isKey[i] = 0;
    for (j=1;j<=nKeyWords; j++) {
      if (!strncmp(tWord,headerWord[j],li)) {
	isKey[i] = 1;
	j = nKeyWords +1;
      }
    }
  }
  /* Alloc outWords */
  i = 0;
  i++; strncpy(headerWords[i], "RUNSPEC",5);
  i++; strncpy(headerWords[i], "TITLE",5);
  i++; strncpy(headerWords[i], "DIMENS",6);
  i++; strncpy(headerWords[i], "OIL",3);
  i++; strncpy(headerWords[i], "GAS",3);
  i++; strncpy(headerWords[i], "WATER",5);
   i++; strncpy(headerWords[i], "DISGAS",6);
 i++; strncpy(headerWords[i], "METRIC",6);
