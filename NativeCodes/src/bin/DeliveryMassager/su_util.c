/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "su_util.h" 

void swap_to_BigEndian_if_needed(segy *tr)
{
  int i, ns, machine_is_little_endian;

  i=1;
  if (*(char *)&i == 1) machine_is_little_endian = 1;
  else machine_is_little_endian = 0;

  if (machine_is_little_endian)
	{
	  ns = tr->ns;
	  /* swap to foreign (BE) header values */
	  for (i = 0; i < SU_NKEYS; ++i) swaphval(tr,i);
	  /* swap to foreign (BE) data values */
	  for (i = 0; i < ns ; ++i) swap_float_4(&(tr->data[i]));
	}
}    



/** 
	Get index of header word (linear search)
*/
int getindex(const char *key)   /* get index for this key */
{
  register int i;
  
  for (i = 0; i < SU_NKEYS; i++)
	{
	  if (!strcmp(hdr[i].key, key)) return i;       /* key found */
	}
  /* not found */
  return -1;
}


/*
  Put header value at known index value.
  value is pointer to either short or int value
*/
void puthval(segy *tr, int index, void *valp)
{
  char *tp = (char*) tr;
  
  switch(*(hdr[index].type)) 
	{
	case 'h': *((short*)  (tp + hdr[index].offs)) = *((short*)valp); break;
	case 'u': *((unsigned short*) (tp + hdr[index].offs)) = *((unsigned short*)valp); break;
	case 'i': *((int*)   (tp + hdr[index].offs)) = *((int*)valp); break;
	case 'p': *((unsigned int*)   (tp + hdr[index].offs)) = *((unsigned int*)valp); break;
	case 'l': *((long*)   (tp + hdr[index].offs)) = *((long*)valp); break;
	case 'v': *((unsigned long*)  (tp + hdr[index].offs)) = *((unsigned long*)valp); break;
	case 'f': *((float*)  (tp + hdr[index].offs)) = *((float*)valp); break;
	case 'd': *((double*) (tp + hdr[index].offs)) = *((double*)valp); break;
	case 's': (void) strcpy(tp + hdr[index].offs, ((char*)valp));  break; /* JG: never occurs I think */
	default: fprintf(stderr,"%s: %d: mysterious data type", __FILE__,__LINE__); break;
	}  
  return;
}

void swaphval(segy *tr, int index)
{
  register char *tp= (char *) tr;

  switch(*(hdr[index].type)) 
	{
	case 'h': swap_short_2((short*)(tp + hdr[index].offs)); 
	  break;
	case 'u': swap_u_short_2((unsigned short*)(tp + hdr[index].offs));
	  break;
	case 'i': swap_int_4((int*)(tp + hdr[index].offs));
	  break;
	case 'p': swap_u_int_4((unsigned int*)(tp + hdr[index].offs));
	  break;
	case 'l': swap_long_4((long*)(tp + hdr[index].offs));
	  break;
	case 'v': swap_u_long_4((unsigned long*)(tp + hdr[index].offs));
	  break;
	case 'f': swap_float_4((float*)(tp + hdr[index].offs));
	  break;
	case 'd': swap_double_8((double*)(tp + hdr[index].offs));
	  break;
	default: fprintf(stderr,"%s: %d: unsupported data type", __FILE__, __LINE__);
	  break;
	}
  return;
}

void swap_short_2(short *tni2)
/**************************************************************************
swap_short_2            swap a short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));  
}

void swap_u_short_2(unsigned short *tni2)
/**************************************************************************
swap_u_short_2          swap an unsigned short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));  
}

void swap_int_4(int *tni4)
/**************************************************************************
swap_int_4              swap a 4 byte integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_u_int_4(unsigned int *tni4)
/**************************************************************************
swap_u_int_4            swap an unsigned integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_long_4(long *tni4)
/**************************************************************************
swap_long_4             swap a long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_u_long_4(unsigned long *tni4)
/**************************************************************************
swap_u_long_4           swap an unsigned long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_float_4(float *tnf4)
/**************************************************************************
swap_float_4            swap a float
***************************************************************************/
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_double_8(double *tndd8)
/**************************************************************************
swap_double_8           swap a double
***************************************************************************/
{
  char *tnd8=(char *)tndd8;
  char tnc;

  tnc= *tnd8;
  *tnd8= *(tnd8+7);
  *(tnd8+7)=tnc;

  tnc= *(tnd8+1);
  *(tnd8+1)= *(tnd8+6);
  *(tnd8+6)=tnc;

  tnc= *(tnd8+2);
  *(tnd8+2)= *(tnd8+5);
  *(tnd8+5)=tnc;

  tnc= *(tnd8+3);
  *(tnd8+3)= *(tnd8+4);
  *(tnd8+4)=tnc;
}


