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
/*
** deliveryMassager.h
**
** Automatically generated header file.
** Source file:   deliveryMassager.c
** Last modified: Monday, March 29, 2004 at 15:26:43
*/

#ifndef __deliveryMassager_h__
# define __deliveryMassager_h__

int     main                            (int , char* []);
int     checkOutputGridisInterior       (struct MASSAGER*);
int     setUpGridMapping                (struct MASSAGER*);
int     readOutputGrids                 (struct MASSAGER*);
int     readMassageAnalysisData         (struct MASSAGER*);
int     setUpWellObservations           (struct MASSAGER *massager);
int     setUpOutputs                    (struct MASSAGER *massager);
struct PT2D getBoundingIndices(struct INPUTGRID *ig,double x,double y);
int     setUpUTMmapping                 (struct MASSAGER *massager); 
struct PT2D getUTMSfromInlineXline(struct UTM_CONVERSION_INFO info,double inLine,double xLine);
int writeEclipseCache(struct MASSAGER *massager,int i);
int readEclipseCache(struct MASSAGER *massager,int i);
int isADesiredECLIPSEKeyword(char* s);
int isABlockCentredECLIPSEKeyword(char* s);
int cleanUpGrids                        (struct MASSAGER *massager);
int cleanUpSundryMemory                 (struct MASSAGER *massager);
int addSystematicDepthUncertainty       (struct MASSAGER *massager);
int freeMassageAnalysisCovariances      (struct MASSAGER *massager);
int freeMassageAnalysisData             (struct MASSAGER *massager);
int isEclipseCommonSource(struct MASSAGER *massager);
int determineEclipseOutputMode(struct MASSAGER *massager);
void checkXsection(char *suffix,int iXsection,int nch, int ndh, double ****zCorn);
void checkYsection(char *suffix,int jXsection,int nrh, int ndh,double ****zCorn);

#endif
