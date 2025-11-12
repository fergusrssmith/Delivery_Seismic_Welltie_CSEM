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
** solver.h
*/

#ifndef __solver_h__
# define __solver_h__


int      setUpMapping(struct MASSAGER*);
int      solveScaledSmoothSparseProblems(struct MASSAGER *massager);
int      dumpPostMassagedSUFiles(struct MASSAGER *massager,float ***krigMeans,float ***krigStdDev);
int      intPairCompare(const void *a, const void *b);
void     writeAsciiHeaderInfo(int layerNumber,struct MASSAGER *massager,int isRealisationsFile);
int      computeMovingAverage(struct MASSAGER *massager);
void     formMultigridTraverseOrder(struct INTPAIR* mgIndices, int size);
void     solveOKsystem(double** Ckrig,double* RHS,int n);
double **getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(struct MASSAGER *massager, 
																 int layer_index, int output_grid_index, 
																  int debug, double eps_last);
int      cleanUpMappingSetup(struct MASSAGER *massager);
void     showStatistics(struct MASSAGER *massager, float **stats,int isComposite,int layerNumber,int Nsamples, int Nprops,char **names);

int pairCmp(const void *a, const void *b);
void print_dmatrix(double **m,int n);
void print_matrix(void **m,int n,char type,char* name);

#endif /* end of header file */
