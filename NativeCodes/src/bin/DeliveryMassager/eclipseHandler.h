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
#ifndef __eclipseHandler_h__
# define __eclipseHandler_h__

int getDepthRangeInputGrid(struct MASSAGER *massager, double *zMin, double *zMax);
int makeEclipseGrid(struct MASSAGER *massager);
int getEclipseGriddingFromBoundingRectangle(int n, double *x,double *y, struct NEAREST_REGULAR_GRID_INFO *RGinfo);
int writeEclipseRealisation(struct MASSAGER *massager,int real_no,float ***r_properties);
struct SEG_ENCODING **fetchSegLists(ECL *ecl,int *nNodes, int verbosity);
void search_in_seg_array(GPtrArray *segArray,int indx, int val, int *segCount);
void free_segInfo(struct MASSAGER *massager,struct SEG_ENCODING **m,int Nx,int Ny);
void allocateStochEclipseMem(struct MASSAGER *massager);
void freeStochEclipseMem(struct MASSAGER *massager);
int isInEclipseBlock(ECL *ecl,int ecl_i,int ecl_j,double z,double* testPt);
int getSegOfEclipseBlock(ECL *ecl,int **segcol,int ecl_i,int ecl_j,double z,double* testPt);
int flipGridAndMapUnits(ECL *ecl);
int doBilinearInterp(double* pt_cnrs,double* r,double *u,double *v);
void debug_plotTestPt(double* pt_cnrs,double* r);


#endif
