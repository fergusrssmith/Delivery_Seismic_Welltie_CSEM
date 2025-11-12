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
#include "ECLIO/eclio.h"
#include "ECLIO/records.h"
#include "ECLIO/eclutil.h"
/***************************************************************************/
int getNSlashReqd(char *key,  ECL *ecl)
/***************************************************************************/
{
  if (strncmp(key,"COORDSYS", TOKEN_LENGTH)!= 0) {
    if ( ecl->coordSys.set == 1) 
      return (ecl->specGrid.numRes);
    else return (1);
  }
  return (1);
}
/*************************************************************************/
int cleanRecord(char *buffer)
/*************************************************************************/
{
#if JG_FEATURES
  int j;
  int len=strlen(buffer)-1;

  for (j=len;j>=0;j--) 
	{
	  if (buffer[j] == '/') 
		{
		  buffer[j]=' ';
		  break;
		}
	}
  return ((int)j);
#else

  size_t i = 0;

  for (i=strlen(buffer)-1;i>=0;i--) 
	{
	  if (buffer[i] == '/') 
		{
		  buffer[i]=' ';
		  break;
		}
	}
  return ((int)i);

#endif
}

/*************************************************************************/
int scrubRecord(char *buffer, char *elimChars, char replChar)
/*************************************************************************/
{
  size_t i = 0;
  size_t lb = 0;
  
  lb = strlen(buffer);

  for (i=0;i<lb;i++) 
    if (strchr(elimChars,buffer[i])) buffer[i]=replChar;
  return (0);
}
/*************************************************************************/
int seekToken1(char *buffer, char *key, int toUpper)
/*************************************************************************/
{
  /*
	Assumes key is already a clean, compact, string in the correct case
		(has already had makeTokenPristin run on it) 
	Returns 0 on successfully finding keyIn 
             1 if key not found
			 <0 for null buffer or key */
  
  char token[TOKEN_LENGTH+1];
  int status = 0;
  
  if (!key) return (-1);
  
  status=       makeTokenPristine(buffer,token, (size_t) TOKEN_LENGTH,
				  ECL_TOKENS, toUpper);
  
  if (status) return (-2);
  if (!strcmp(key, token)) return (0);
  return (1);
}
/*************************************************************************/
int makeTokenPristine(char *src, char *dest, size_t maxDestSize, 
		      char *separators, int forceUp)
/*************************************************************************/

/* Tries to extract a "pristine" keyword from the buffer, returns in dest
	leaving buffer unaffected 
	Returns 0 for success; 1 or 2 for nul strings, 3 for null separators, and
	4 for maxDestSize <= 1
	*/
{
  int i, nw = 0, iReturn=0;
  char **words = (char **) NULL;
  
  if (!src) {
    iReturn = 1;
    goto bottom;
  }
  if (!separators) {
    iReturn=3;
    goto bottom;
  }
  if (maxDestSize <= 1) {
    iReturn=4;
    goto bottom;
  }
	
  nw = strtokSafe(src,separators, &words, 1);
  if (nw<1 || words[0] == (char *) NULL) {
    iReturn = 2;
    goto bottom;
  }
  
  strncpy(dest, words[0], maxDestSize);
  cleanRecord(dest);
  
  if (forceUp) for (i=0;i<strlen(dest);i++) dest[i]=toupper(dest[i]);
  
 bottom:
  freeWordList(words,nw);
  return (0);
}

/***************************************************************************/
/* Finds the end of an ECL record, by detecting the requisite number of "/"'s */
int findEndOfRecord( ECL *ecl, char *key, FILE *f)
/***************************************************************************/
{
  int nSlashReqd, nSlashRead=0, isSlash=0;

  nSlashReqd = getNSlashReqd(key, ecl);
  while (nSlashRead < nSlashReqd) {
    isSlash = findSlash(f, "--");
    if (isSlash == 1) nSlashRead++;
    if (isSlash < 0) {
      eclError("findEndOfRecord", 1,
	       "\tCould not find terminating / \n");
      return (1);
   }
  }	
  return (0);
}
/***************************************************************************/
int findSlash(FILE *f, char *commentString)
/***************************************************************************/
{
  char buffer1[LINE_LENGTH+1], buffer2[LINE_LENGTH+1];
  char key[LINE_LENGTH+1];
  int status = 0;
  int forceUp = 1;
  int lkey;
  int nw = 0;
  char **words = (char **) NULL;
  int iReturn = 0;

  iReturn = 0;

  while (fgets(buffer1, (int) LINE_LENGTH, f)) {        /* work thru file */

    status = removeTrailingComment(buffer1, buffer2, (size_t) LINE_LENGTH, commentString);
    if (status < 0) {
      iReturn = status; /* screwup */
      goto bottom;
    }
    
    nw = strtokSafe(buffer2,ECL_TOKENS,&words,-1);
    if (nw < 1) {
      iReturn = 0; /* no words on line */
      goto bottom;
    }

    /* Look at last item */
    status = makeTokenPristine(words[nw-1],key, (size_t) LINE_LENGTH, 
				  ECL_TOKENS, forceUp);
    if (status) {
      iReturn =  -10*status; /* error making key pristine */    
      goto bottom;
    }
    if ((lkey = strlen(key)) > 0 && (key[lkey-1] == '/')) iReturn = 1;
  }
 bottom:
  freeWordList(words,nw);
  return (iReturn);
}
int removeComment (char *b1, char *b2, char *commentString)
{
  /* removes comments NONDESTRUCTIVELY, returning b1 less padding whitespace
     at ends and less any comment strings */

  int l1=0, l2=0, i, iBegin=0;
  int offSet = 0;
  char *beginComment = (char *) NULL;
  
  l1 = strlen(b1);
  if (!l1) return(0);

  /* remove padded spaces at beginning and end*/
  for (i=0,l2=0;i<l1;i++) {
    if (!strchr(ECL_WHITESPACE,b1[i])) {
      l2 = l1 - i;
      iBegin = i;
      break;
    }
  }
  // for (i=0;i<l2;i++) b2[i]=b1[i+iBegin]; // JG
  memcpy(b2,&b1[iBegin],l2);

  for (i=l2-1;i>=0;i--) {
    if (strchr(ECL_WHITESPACE,b2[i])) {
      b2[i]='\0';
      l2--;
    }
    else break;
  }
      
  if ((beginComment = strstr(b2,commentString))) 
	{
	  iBegin = (int)(beginComment - b2)/sizeof(char);
	  //for (i=iBegin;i<l2;i++) b2[i] = '\0';  JG
	  memset(&b2[iBegin],'\0',(l2-iBegin));
	  //l2 = strlen(b2); JG
	}
  return (0); 
}
int removeTrailingComment(char *b1, char *b2, size_t lineLen, char *commentString)
{
  int lenIn=0, lenOut=0, i;
  char *beginComment = (char *) NULL, *b1p = (char *) NULL;
  
  lenIn = strlen(b1);  
  if (!lenIn) return(0); 

  if ((beginComment = strstr(b1,commentString))) {
    for (b1p = b1, i=0, lenOut=0; b1p < beginComment && i < lenIn; b1p++,i++) {
      b2[i]=b1[i];
      lenOut++;
    }
    b2[lenOut]='\0';
  }
  else strncpy(b2,b1,lenIn);

  return (0);
}

void freeWordList(char ** words, int nWords)
{
  int i;

  for (i=0;i<nWords;i++) 
    {
	  if (words[i])
		{
		  free ((char *) words[i]);
		  // JG_MEMORY
		  words[i] = NULL;
		}
	}
  if (words) free((char *) words); // JG_MEMORY
}

int strtokSafe2(char *buffer, char *tokens, char ***wordp, int nToGet)
     /* Similar functionality to strtok, but no statics and
	less dangerous. Returns the number of token-delimited words;
	words[0:return() -1] is an array of character strings, each
	null-terminated.
	Returns -1 for failure to allocate
	0 if no words found
     Top routine*/
{
  int WORDCHUNK = 50;
  char **words = (char **) NULL;
  char *b = (char *) NULL, *bb = (char *) NULL, *be = (char *) NULL;
  char t0 = ' ';
  int nw = 0, nwAlloc = 0, nt = 0;
  int lbuffer=0, lb=0, lw=0;
  int ib,it,iw;

  if (nToGet > 0) WORDCHUNK = nToGet;
  else nToGet = I_BIG; 
  
  /* set up words */
  nwAlloc = WORDCHUNK;
  nw = 0;
  words = (char **) malloc((size_t) nwAlloc * sizeof(char *));
  if (!words) return (-1);

  /* set up token list */
  nt = strlen(tokens);
  if (nt<1) return (-2);
  t0 = tokens[0];

  /* set up temp buffer b */
  lbuffer = strlen(buffer);
  b = (char *) malloc((size_t) (lbuffer + 1) * sizeof(char));
  if (!b) return (-1);
  for (ib=0;ib<lbuffer;ib++) b[ib]='\0';
  strncpy(b,buffer,lbuffer); b[lbuffer]='\0';
  lb = strlen(b);

  /* before actually breaking it up, we will change all tokens to 
     tokens[0] */
  for (ib=0;ib<lb;ib++) {
    for (it=1;it<nt;it++) {
      if (b[ib] == tokens[it]) {
	b[ib]=t0;
	break;
      }
    }
  }
  for (ib=lbuffer-1;ib>=0 && b[ib]==t0;ib--) {b[ib]='\0';}
  lb = strlen(b);

  /* Initialize pointer to beginning of b, and step through word-by-word */
  nw = 0; 
  be = b;
  iw = 0;
  while (iw < lb && nw < nToGet){

    for (bb=be;*bb==t0 && iw<lb;bb++){iw++;} /*find beginning*/
    if (*bb == t0) {bb++;iw++;}
    for (be=bb;*be != t0 && iw<lb;be++){iw++;} /*find end */
    if (*be== t0) {be--;iw--;}
    nw++;

    if (nw > nwAlloc) { /*make sure there's space in the heap */
      nwAlloc += WORDCHUNK;
      words = (char **) realloc(words,(size_t) nwAlloc * sizeof(char *));
      if (!words) return (-1);
    }

    lw = (be-bb+1)/sizeof(char); /*allocate and copy this word*/
    words[nw-1] = (char *) malloc((size_t)(lw+1)*sizeof(char));
    if (!words[nw-1]) return (-1);
    
    for (ib=0;ib<=lw;ib++) words[nw-1][ib]='\0';
    strncpy(words[nw-1],bb,lw);
    be++; iw++;
  }

  if (nwAlloc > nw) 
    words = (char **) realloc(words, (size_t) (nw) * sizeof(char *));
  free ((char *) b);
  *wordp = words;
  return (nw);
}



/* JG's verson for speed... */
#define MAX_WORDS 50
int strtokSafeJG(char *buffer, char *tokens, char ***wordp, int nToGet)
//int strtokSafe(char *buffer, char *tokens, char ***wordp, int nToGet)
{
  int asprintf(char **ret, const char *format, ...);

  size_t ptr_static[MAX_WORDS]; 
  size_t *ptr,*ptr_dynamic=NULL;
  int i,max,ptr_size,nw=0;
  char *brkt, *buffer_cp;

  char *base_buffer;
  char* word;
  size_t base;

  ptr_size = MAX_WORDS;
  ptr = &ptr_static[0];

  // loop to cope with really long lines, if encountered
  for(;;)
	{
	  // copy to work array
	  asprintf(&buffer_cp,"%s",buffer); 	  
	  base_buffer = buffer_cp;
	  base = (size_t)base_buffer;

	  if (nToGet<0) 
		{
		  max = ptr_size;
		}
	  else
		{
		  max =  (nToGet<ptr_size) ? nToGet: ptr_size;
		}



	  nw=0;
	  for (word = strtok_r(buffer_cp, tokens, &brkt);
		   word && (nw<max);
		   word = strtok_r(NULL, tokens, &brkt))
		{
		  ptr[nw] = (size_t)(word - base);
		  nw++;
		}
	  if (nw==0) 
		{
		  free(buffer_cp);
		  return -1; 
		}

	  if (nw==ptr_size && word!=NULL)
		{
		  /* buffer filled, with extra tokens unaccounted for */
		  /* try again, with bigger buffer */
		  free(buffer_cp);
		  if (ptr_size != MAX_WORDS) 
			{
			  /* this pass was already a retry: clean up before next try */
			  free(ptr_dynamic);
			} 
		  ptr_size *= 2;

		  ptr_dynamic = (size_t*)malloc(ptr_size*sizeof(size_t));
		  if (!ptr_dynamic) return -1;

		  ptr = ptr_dynamic;

		}
	  else
		{
		  break;
		}
	}

  *wordp = (char**)malloc(nw * sizeof(char*)); 
  if (!(*wordp)) return -1;

  for(i=0;i<nw;i++)
	{
	  max = asprintf(&((*wordp)[i]),"%s",&(base_buffer[ptr[i]])); 
	  if (max==-1) return -1;
	}
  free(buffer_cp);

  if (ptr_size != MAX_WORDS) free(ptr_dynamic);

  return nw;
}
#undef MAX_WORDS


int strtokSafe(char *buffer, char *tokens, char ***wordp, int nToGet)
     /* Similar functionality to strtok, but no statics and
	less dangerous. Returns the number of token-delimited words;
	words[0:return() -1] is an array of character strings, each
	null-terminated.
	Returns -1 for failure to allocate
	0 if no words found
	** New version using manageHeap **
     Top routine*/
{
  int WORDCHUNK = 50;
  char **words = (char **) NULL;
  char *b = (char *) NULL;
  char *bb = (char *) NULL, *be = (char *) NULL;
  char t0 = ' ';
  // OLD (breaks 64bit binaries) 
  //  int nw=0, nwAlloc=0, nt=0;
  int nw=0, nt=0;
  size_t nwAlloc=0;
  int lbuffer=0, lb=0, lw=0;
  int ib,it,iw;

  *wordp = words;
  
  if (nToGet > 0) WORDCHUNK = nToGet;
  else nToGet = I_BIG; 
  
  /* set up words */
  nwAlloc = 0;
  words = (char **) manageHeap((GENERIC_POINTER) NULL,WORDCHUNK, &nwAlloc, 
 			       sizeof(char **), WORDCHUNK,0);
  if (!words) return (-1);

  /* set up token list */
  nt = strlen(tokens);
  if (nt<1) return (-2);
  t0 = tokens[0];

  /* set up temp buffer b ; we don't change buffer in this routine */
  lbuffer = strlen(buffer);
  b = (char *) malloc((size_t) (lbuffer + 1) * sizeof(char));
  if (!b) return (-1);

  for (ib=0;ib<lbuffer;ib++) b[ib]='\0';
  strncpy(b,buffer,lbuffer); b[lbuffer]='\0';
  lb = strlen(b);

  /* before actually breaking it up, we will change all tokens to 
     tokens[0] */
  for (ib=0;ib<lb;ib++) {
    for (it=1;it<nt;it++) {
      if (b[ib] == tokens[it]) {
	b[ib]=t0;
	break;
      }
    }
  }
  /* and remove trailing tokens */
  for (ib=lbuffer-1;ib>=0 && b[ib]==t0;ib--) {b[ib]='\0';}
  lb = strlen(b);

  /* Initialize pointer to beginning of b, and step through word-by-word */
  nw = 0; 
  be = b;
  iw = 0;
  while (iw < lb && nw < nToGet){

    for (bb=be;*bb==t0 && iw<lb;bb++){iw++;} /*find beginning*/
    if (*bb == t0) {bb++;iw++;}
    for (be=bb;*be != t0 && iw<lb;be++){iw++;} /*find end */
    if (*be== t0) {be--;iw--;}
    nw++;

    if (nw) {
      words = (char **) manageHeap((GENERIC_POINTER) words, nw, &nwAlloc, 
				   sizeof(char **), WORDCHUNK,0);
      if (!words) return(-1);

      lw = (be-bb+1)/sizeof(char); /*allocate and copy this word*/
      words[nw-1] = (char *) malloc((size_t)(lw+1)*sizeof(char));
      if (!words[nw-1]) return (-1);
      
      for (ib=0;ib<=lw;ib++) words[nw-1][ib]='\0';
      strncpy(words[nw-1],bb,lw);
      be++; iw++;
    }
  }

  if (nw) words = (char **) manageHeap((GENERIC_POINTER) words, nw, &nwAlloc, 
				       sizeof(char **), WORDCHUNK,1);
  else {
    freeWordList(words,nw);
    words = (char **) NULL;
  }
  if (b) free ((char *) b);
  *wordp = words;
  return (nw);
}
