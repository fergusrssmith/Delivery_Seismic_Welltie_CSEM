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
#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <search.h>

#include <glib.h>

#include "NR/nrutil.h"
#include "NR/nr.h"
#include "bayes/bayes.h"

#include "kdtree/nn.h"
#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"

#include "std/std.h"
#include "std/error.h"


#include "../../lib/ATLAS/include/atlas_misc.h"
#include "../../lib/ATLAS/include/clapack.h"


#include "massager_types.h"
#include "variogram.h"
#include "solver.h"
#include "util.h"
#include "deliveryMassager.h"


#define PRIOR_STANDARD_DEVIATION_BOOST_FACTOR 2            // Prior model is taken to have a flat mean = data mean, 
                                                           // std dev = PRIOR_STANDARD_DEVIATION_BOOST_FACTOR * data std dev...
#define DETRENDED_PRIOR_STANDARD_DEVIATION_BOOST_FACTOR 3  // Prior model is taken to have mean = smoothed data map
#define DEFAULT_MAX_NEIGHBOURS_IN_SPARSE_PRIOR 80          // Number of neighbours used in local Bayesian kriging systems, *per layer, per property*
#define DEFAULT_MAX_NEIGHBOURS_IN_SIMULATIONS  30          // Number of neighbours used in sequential simulation
#define DEFAULT_N_NEIGHBOURS_IN_INTERPOLATION_KRIGING 8    // Number of neighbours used in interpolation
#define KRIG_WEIGHT_THRESHOLD 0.001                        // interpolation will drop neighbours with weaker contributions than this

/*
  Set up mappings between structured input grids and unstructured output grid using kriging
*/

int setUpMapping(struct MASSAGER *massager)
{
  int i,j,jj,k,kk,dimension;
  int info;   
  double** Points;
  double** Ckrig;
  double* RHSkrig;
  int *ipiv_krig;
  double x,y;
  optkdNode *OptkdRoot;
  int* neigh;
  int dim=2;
  int   nNeighborsDesired = DEFAULT_N_NEIGHBOURS_IN_INTERPOLATION_KRIGING;
  int   outputgrid_indx_offset;
  int  input_grid_indx;
  int  output_grid_indx;
  int count;
  double **A=NULL;
  double pt[2];
  int totalNSoftDataObservations;
  int totalNWellDataObservations=0;
  struct PT2D pt_utm;

  struct SOLVER_INFO *sI;


  sI = &massager->solverInfo;  // Useful shortcut

  // size of interpolation matrix
  // a row for each for (from "delivery") data point
  totalNSoftDataObservations = massager->nLayers * massager->nPropertyNames * massager->inputgrid.nPoints;
  // plus one for each well observation
  for (i=0; i<  massager->nLayers ; i++)
	{
	  totalNWellDataObservations += massager->well_observations[i].n;
	}
  sI->nRowsA = totalNSoftDataObservations + totalNWellDataObservations;
  sI->nColumnsA = massager->totalNOutputGridPoints * massager->nPropertyNames;


  // Allocate
  if (massager->useDenseMethod)
	{
	  // dense  (hope there's room :-)
	  sI->A = dmatrix(0, sI->nRowsA - 1, 0, sI->nColumnsA - 1);
	}
  else
	{
	  // Store array of kriging weights for each grid, and coupled information
	  sI->IGkrigInfoSoft = (struct PAIR_ARRAY**)malloc((size_t) ((massager->nLayers) * sizeof(struct PAIR_ARRAY*)));
	  sI->WellkrigInfo = (struct PAIR_ARRAY**)malloc((size_t) ((massager->nLayers) * sizeof(struct PAIR_ARRAY*)));
	  sI->OGkrigLookupInfoSoft = (struct INDXTRIPLE**)malloc((size_t) ((massager->nLayers) * sizeof(struct INDXTRIPLE*)));
	  sI->OGkrigInfoWell = (struct INDXTRIPLE**)malloc((size_t) ((massager->nLayers) * sizeof(struct INDXTRIPLE*)));
	  sI->OGkrigInfoSoft = (struct PAIR_ARRAY**)malloc((size_t) ((massager->nLayers) * sizeof(struct PAIR_ARRAY*)));

	  for (i=0; i< massager->nLayers ; i++)
		{
		  sI->IGkrigInfoSoft[i] = (struct PAIR_ARRAY*)malloc((size_t) 
								  ((massager->inputgrid.nPoints) * sizeof(struct PAIR_ARRAY)));
		  sI->WellkrigInfo[i] = (struct PAIR_ARRAY*)malloc((size_t) 
								((massager->well_observations[i].n) * sizeof(struct PAIR_ARRAY)));
		  // calloc these, so the size elements are initialised to zero
		  sI->OGkrigLookupInfoSoft[i] = (struct INDXTRIPLE*)calloc((size_t) (massager->outputgrid[i].nPoints),
																				  sizeof(struct INDXTRIPLE));
		  // FIXME later on to be more memory efficient, this could be a hash table...
		  sI->OGkrigInfoWell[i] = (struct INDXTRIPLE*)calloc((size_t) (massager->outputgrid[i].nPoints),
																				  sizeof(struct INDXTRIPLE));
		  sI->OGkrigInfoSoft[i] = (struct PAIR_ARRAY*)malloc((size_t) 
								   ((massager->outputgrid[i].nPoints) * sizeof(struct PAIR_ARRAY)));
		}
									   
	}
  if (massager->useDenseMethod)
	{
	  // shortcut to increase legibility later on
	  A = sI->A;
	  // zero
	  memset(&A[0][0],0,sI->nRowsA * sI->nColumnsA * sizeof(double));
	}

  input_grid_indx=0;
     
  // kriging stuff etc needs to be performed in UTM coordinates if at all possible
  // massager->outputgrid[i] will be in seismic coords if massager->duplicate applies
  if (massager->duplicate && massager->utmInfo.useUTMconversion)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}


  // loop over layers
  for (i=0; i<massager->nLayers; i++)
	{
	  // (1) form kdtree of outputgrid, in UTM's if at all possible
	  massager->outputgrid[i].Points = dmatrix(0,massager->outputgrid[i].nPoints-1,0,1);
	  Points = massager->outputgrid[i].Points; // useful shortcut
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  if (massager->duplicate && massager->utmInfo.useUTMconversion)
			{
			  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
										  massager->outputgrid[i].xValues[j],
										  massager->outputgrid[i].yValues[j]);
			  Points[j][0] = pt_utm.x;
			  Points[j][1] = pt_utm.y;
			}
		  else
			{
			  Points[j][0] = massager->outputgrid[i].xValues[j];
			  Points[j][1] = massager->outputgrid[i].yValues[j];
			}
		}
	  // build kdtree
	  massager->outputgrid[i].OptkdRoot = BuildOptTree(Points,massager->outputgrid[i].nPoints,2);
	  OptkdRoot = massager->outputgrid[i].OptkdRoot;
	  // set up appropriate metric for the kriging neighbour search using a kdtree, 
	  // suitable for ith output grid
	  setUp2DEuclidenGenMetric(OptkdRoot->Kglobals,
							   M_PI*massager->variograms[i].azimuth/180.0,
							   massager->variograms[i].rangeX/massager->variograms[i].rangeY);
	  // Be careful about number of neighbours
	  if (nNeighborsDesired > massager->outputgrid[i].nPoints ) 
		{
		  nNeighborsDesired = massager->outputgrid[i].nPoints-1;   
		}
	  // workspace for the OK kriging 
	  Ckrig = dmatrix(0,nNeighborsDesired,0,nNeighborsDesired);
	  RHSkrig = dvector(0,nNeighborsDesired);
	  ipiv_krig = ivector(0,nNeighborsDesired);
	  // set up variogram for kriging on this ith outputgrid
	  setUpcorr2D(&massager->variograms[i]);

 
	  // (2) run over well observations first
	  for (j=0; j<massager->well_observations[i].n; j++)
		{
		  // x,y of jth observation in this layer. In UTMs unless there is no conversion info
		  pt[0] = massager->well_observations[i].xValues[j];
		  pt[1] = massager->well_observations[i].yValues[j];
		  // get list of outputgrid neighbours of this point from kdtree
		  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,OptkdRoot,1);
		  // form OK kriging system
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  // Matrix. Need form only lower triangle
			  for (kk=1; kk <= k ; kk++) 
				{
				  // last index is "dimension" => 0 = x , 1 = y
				  Ckrig[k-1][kk-1] = corr2D(Points[neigh[k]][0],  Points[neigh[k]][1],
											Points[neigh[kk]][0], Points[neigh[kk]][1],
											&(massager->variograms[i]));
				  Ckrig[kk-1][k-1] = Ckrig[k-1][kk-1]; // FIXME needed?
				}
			  // RHS
			  RHSkrig[k-1] = corr2D(pt[0],pt[1],Points[neigh[k]][0],Points[neigh[k]][1],&(massager->variograms[i]));				             
			}
		  // edges
		  for (kk=0; kk < nNeighborsDesired ; kk++) Ckrig[nNeighborsDesired][kk] = Ckrig[kk][nNeighborsDesired] = 1.0;
		  Ckrig[nNeighborsDesired][nNeighborsDesired]=0.0;
		  RHSkrig[nNeighborsDesired] = 1.0;

		  // Now solve kriging system for the weights mapping the output grid to the well-observation point
		  // print_dmatrix0(Ckrig,nNeighborsDesired+1,"junk.txt");
		  info = clapack_dgetrf(CblasColMajor,nNeighborsDesired+1,nNeighborsDesired+1,
								&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0]);
		  if (info>0)
			{
			  warning ("lapack decomposition failed in kriging step\n");
			  return EXIT_FAILURE;
			}   

		  /* Get weights via dge solution solution */
		  // Cholesky: clapack_dpotrs(CblasColMajor,CblasLower,nNeighborsDesired,1,&Ckrig[0][0],nNeighborsDesired,&RHSkrig[0],nNeighborsDesired);
		  clapack_dgetrs(CblasColMajor,CblasNoTrans,nNeighborsDesired+1,1,&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired+1);
		  
		  /*
			Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1]. 
			
		  */
		  
		  // Wells: store kriging weights etc.
		  // eliminate negligible contributions
		  if (!massager->useDenseMethod)
			{
			  int krigCount=0;
			  for (k=0 ; k <nNeighborsDesired; k++)
				{
				  if (fabs(RHSkrig[k]) > KRIG_WEIGHT_THRESHOLD) krigCount++;
				}
			  // allocate room and store
			  sI->WellkrigInfo[i][j].size = krigCount;
			  sI->WellkrigInfo[i][j].values = 
				(struct PAIR*)malloc((size_t) ((krigCount) * sizeof(struct PAIR)));
			  for (k=0,kk=0 ; k <nNeighborsDesired; k++)
				{
				  if (fabs(RHSkrig[k]) > KRIG_WEIGHT_THRESHOLD) 
					{
					  sI->WellkrigInfo[i][j].values[kk].val =  RHSkrig[k];  // kriging weight
					  sI->WellkrigInfo[i][j].values[kk].indx = neigh[k+1];  // output grid index
					  /*
						printf("Well Obs layer %d, number %d: OGindx=%d, wt=%g\n",i,j,
							 sI->WellkrigInfo[i][j].values[kk].indx,
							 sI->WellkrigInfo[i][j].values[kk].val);
					  */
					  // increment contributions count for output grid node neigh[k+1]
					  sI->OGkrigInfoWell[i][neigh[k+1]].size++;
					  kk++;
					}
				}
			}

		 if (massager->useDenseMethod)
		   {
			 /* Kriging weights are the non-zero entries in the interpolation operator A 
				mapping the output grid(s) to the well observation point.
				Store these weights in the correct position in A:
				Ordering of (y_I - A m_O) vector is as per input grid,
				(1) properties fastest, then 
				(2) nodes
				(3) layers. 
				Number of nodes on input grid
				is constant per layer, as the input grid is the same.
				Output grid ordering has the same sequence, but there may be different numbers of points per
				layer
			 */
			 /* Precompute a useful output grid index offset */
			 for (k = outputgrid_indx_offset = 0; k<i; k++) outputgrid_indx_offset += massager->outputgrid[k].nPoints;
			 outputgrid_indx_offset *= massager->nPropertyNames;

			 // we are looking at an observation of (in layer i) massager->well_observations[i].property_name[j]
			 // which is formed by kriging neighbours of the same property name in this this.
			 // We need to find the index of this property name
			 // loop over property names
			 for (k=0; k< massager->nPropertyNames; k++)
			   {
				 if (!strcmp(massager->well_observations[i].property_name[j], massager->propertyNames[k]))
				   {
					 // ok, this observation is of property k.
					 // Make an entry in the
					 // observations linear algebra
					 // loop over neighbours contributing to kriging
					 for (kk=1 ; kk <=nNeighborsDesired; kk++)
					   {
						 // grid index of neighbour = neigh[kk], kk=1...nNeighborsDesired				  
						 output_grid_indx = k + massager->nPropertyNames * neigh[kk] + outputgrid_indx_offset;
						 // set entry equal to kriging weight
						 A[input_grid_indx][output_grid_indx] = RHSkrig[kk-1];
						 // printf("A_entry: %d %d %g\n",input_grid_indx,output_grid_indx,A[input_grid_indx][output_grid_indx]);
					   }
					 input_grid_indx++;
					 break; // go on to next observation
				   }
			   }
		   }

		 free(neigh);  // neighbourhood of jth observation
		}


	  // (3) Run over all "soft" observations
	  // Run through seismic input grid and set up kriging weights, which
	  // populate the "interpolation" matrix A mapping the input grid to the
	  // ith output grid

	  
	  // now loop over seismic points
	  for (j=0; j < massager->inputgrid.nPoints; j++)
		{
		  // point on input grid where we want to form the kriging estimate. In utm's if at all possible
		  // Recall that Points is in UTMs if at all possible too
		  if (massager->duplicate && massager->utmInfo.useUTMconversion)
			{
			  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
										  massager->inputgrid.xValues[j],
										  massager->inputgrid.yValues[j]);
			  pt[0] = pt_utm.x;
			  pt[1] = pt_utm.y;
			}
		  else
			{
			  pt[0] = massager->inputgrid.xValues[j];
			  pt[1] = massager->inputgrid.yValues[j];
			}

		  // get list of outputgrid neighbours of this point from kdtree
		  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,OptkdRoot,1);
		  // form kriging system
		  // the list of indexes in Points[] of the nearest neighbours is in neigh[1...nNeighborsDesired]
		 
		  
		  // printf("%g %g: ",pt[0],pt[1]);
		  /*
			for (k=1; k <= nNeighborsDesired; k++) 
			{
			  printf("NODES %d %g %g\n",j,massager->outputgrid[i].xValues[neigh[k]],
					 massager->outputgrid[i].yValues[neigh[k]]);
			}
		  */

		  // form OK kriging system
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  // Matrix. Need form only lower triangle
			  for (kk=1; kk <= k ; kk++) 
				{
				  // last index is "dimension" => 0 = x , 1 = y
				  Ckrig[k-1][kk-1] = corr2D(Points[neigh[k]][0],  Points[neigh[k]][1],
											Points[neigh[kk]][0], Points[neigh[kk]][1],
											&(massager->variograms[i]));
				  Ckrig[kk-1][k-1] = Ckrig[k-1][kk-1]; // FIXME needed?
				}
			  // RHS
			  RHSkrig[k-1] = corr2D(pt[0],pt[1],Points[neigh[k]][0],Points[neigh[k]][1],&(massager->variograms[i]));				             
			}
		  // edges
		  for (kk=0; kk < nNeighborsDesired ; kk++) Ckrig[nNeighborsDesired][kk] = Ckrig[kk][nNeighborsDesired] = 1.0;
		  Ckrig[nNeighborsDesired][nNeighborsDesired]=0.0;
		  RHSkrig[nNeighborsDesired] = 1.0;

		  // Now solve kriging system for the weights. OK is not pos. def, so use general solver.  
		  // info = clapack_dpotrf(CblasColMajor,CblasLower,nNeighborsDesired,&Ckrig[0][0],nNeighborsDesired);
		  info = clapack_dgetrf(CblasColMajor,nNeighborsDesired+1,nNeighborsDesired+1,
								&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0]);
		  if (info>0) 
			{
			  warning ("lapack  decomposition failed in kriging step\n");
			  return EXIT_FAILURE;
			}

		  /* If Ckrig=L.L^T, L^T is now in (c-sense) upper triangle of C, including the diagonal */
		  /* Get weights via Cholesky solution */
		  /* void clapack_spotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
			 const int N, const int NRHS, const float *A, const int lda,
			 float *B, const int ldb)
		  */
		  // clapack_dpotrs(CblasColMajor,CblasLower,nNeighborsDesired,1,&Ckrig[0][0], 
		  //                nNeighborsDesired,&RHSkrig[0],nNeighborsDesired);
		  clapack_dgetrs(CblasColMajor,CblasNoTrans,nNeighborsDesired+1,1,&Ckrig[0][0],
						 nNeighborsDesired+1,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired+1);
		  /*
			Solution (kriging weights) are now in RHSkrig[0...nNeighborsDesired-1]. These are the
			non-zero entries in the interpolation operator A mapping the output grid(s)
			to the input grid (for this layer).
			Store these weights in the correct position in A:
			Ordering of (y_I - A m_O) vector is as per input grid,
			(1) properties fastest, then 
			(2) nodes
			(3) layers. 
			Number of nodes on input grid
			is constant per layer, as the input grid is the same.
			Output grid ordering has the same sequence, but there may be different numbers of points per
			layer
		  */
		  /*
			printf("%g:  ",pt[0]);
			for (k=1; k <= nNeighborsDesired; k++) 
			{
			printf("(%g %4.2f)  ",Points[neigh[k]][0], RHSkrig[k-1]);
			}
			printf("\n\n");
		  */
		  if (!massager->useDenseMethod)
			{
			  // Soft data: store kriging weights etc.
			  // eliminate negligible contributions
			  int krigCount=0;
			  for (k=0 ; k <nNeighborsDesired; k++)
				{
				  if (fabs(RHSkrig[k]) > KRIG_WEIGHT_THRESHOLD) krigCount++;
				}
			  // allocate room for storage of kriging weights and store
			  sI->IGkrigInfoSoft[i][j].size = krigCount;
			  sI->IGkrigInfoSoft[i][j].values = 
				(struct PAIR*)malloc((size_t) (krigCount * sizeof(struct PAIR)));
			  for (k=0,kk=0 ; k <nNeighborsDesired; k++)
				{
				  if (fabs(RHSkrig[k]) > KRIG_WEIGHT_THRESHOLD) 
					{
					  sI->IGkrigInfoSoft[i][j].values[kk].val =  RHSkrig[k];  // kriging weight
					  sI->IGkrigInfoSoft[i][j].values[kk].indx = neigh[k+1];  // output grid index
					  sI->OGkrigLookupInfoSoft[i][neigh[k+1]].size++;
					  kk++;
					}
				}
			}
		 if (massager->useDenseMethod)
		   {
			 /* Precompute a useful output grid index offset */
			 for (k = outputgrid_indx_offset = 0; k<i; k++) outputgrid_indx_offset += massager->outputgrid[k].nPoints;
			 outputgrid_indx_offset *= massager->nPropertyNames;

			 /* loop over properties at this input grid location (x,y), of index j */
			 for (k=0; k< massager->nPropertyNames; k++)
			   {		
				 // loop over neighbours contributing to kriging
				 for (kk=1 ; kk <=nNeighborsDesired; kk++)
				   {
					 // grid index of neighbour = neigh[k]				  
					 output_grid_indx = k + massager->nPropertyNames * neigh[kk] + outputgrid_indx_offset;
					 // set entry equal to kriging weight
					 /*
					   if (A[input_grid_indx][output_grid_indx]>0.2) 
					   {
					   printf("A_entry: %d %d %g\n",input_grid_indx,output_grid_indx,A[input_grid_indx][output_grid_indx]);
					   }
					 */
					 A[input_grid_indx][output_grid_indx] = RHSkrig[kk-1];
					 if (massager->ignore_soft) A[input_grid_indx][output_grid_indx] = 0;
					 // printf("A_entry: %d %d %g\n",input_grid_indx,output_grid_indx,A[input_grid_indx][output_grid_indx]);
				   }
				 input_grid_indx++;
			   }
		   }

		  free(neigh);
		  // just for cleanliness
		  OptkdRoot->Kglobals->nndist = NULL;
		}




	  if (!massager->useDenseMethod)
		{

		  // Set up entries in OGkrigInfo structures
		  // (1) Wells
		  // alloc. FIXME could make this an efficient hashmap later...
		  for (j=0; j<massager->outputgrid[i].nPoints; j++)
			{
			  if (sI->OGkrigInfoWell[i][j].size > 0)
				{
				  sI->OGkrigInfoWell[i][j].values = ivector(0,sI->OGkrigInfoWell[i][j].size-1);
				}
			}
		  // fill
		  for (j=0; j<massager->well_observations[i].n; j++)
			{
			  for (k=0; k<sI->WellkrigInfo[i][j].size ;k++)
				{
				  output_grid_indx = sI->WellkrigInfo[i][j].values[k].indx;
				  count = sI->OGkrigInfoWell[i][output_grid_indx].indx;
				  sI->OGkrigInfoWell[i][output_grid_indx].values[count] = j;
				  sI->OGkrigInfoWell[i][output_grid_indx].indx++;
				  /*
					printf("OGkrigInfoWell layer %d, well obs %d OGindx %d, x,y=%g, %g\n",
						 i,j,output_grid_indx,massager->outputgrid[i].xValues[output_grid_indx],
						 massager->outputgrid[i].yValues[output_grid_indx]);
				  */
				}
			}
		  // (2) Soft data
		  // alloc
		  for (j=0; j<massager->outputgrid[i].nPoints; j++)
			{
			  if (sI->OGkrigLookupInfoSoft[i][j].size > 0)
				{
				  sI->OGkrigLookupInfoSoft[i][j].values = ivector(0,sI->OGkrigLookupInfoSoft[i][j].size-1);
				}
			}
		  // fill
		  for (j=0; j<massager->inputgrid.nPoints; j++)
			{
			  for (k=0; k<sI->IGkrigInfoSoft[i][j].size ;k++)
				{
				  output_grid_indx = sI->IGkrigInfoSoft[i][j].values[k].indx;
				  count = sI->OGkrigLookupInfoSoft[i][output_grid_indx].indx;
				  sI->OGkrigLookupInfoSoft[i][output_grid_indx].values[count] = j;
				  sI->OGkrigLookupInfoSoft[i][output_grid_indx].indx++;
				}
			}
		}



	  // Clean up kdtree
	  KillMetricSpace(OptkdRoot->Kglobals);
	  // Leave tree in place for subsequent searches though

	  if (massager->useDenseMethod)
		{
		  free_dmatrix(massager->outputgrid[i].Points,0,massager->outputgrid[i].nPoints-1,0,1);
		}
	  
	  free_dmatrix(Ckrig,0,nNeighborsDesired,0,nNeighborsDesired);
	  free_dvector(RHSkrig,0,nNeighborsDesired);
	  free_ivector(ipiv_krig,0,nNeighborsDesired);
	}

  // ---------------------------------------------------------------------------------
  //  Reverse interpolation structure from input grid to output grid
  // ---------------------------------------------------------------------------------

  // Form kdtree of inputgrid. Map to utm's if desired
  massager->inputgrid.Points = dmatrix(0,massager->inputgrid.nPoints-1,0,1);
  Points = massager->inputgrid.Points; // useful shortcut
  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  /* add utm stuff if need be */
	  if (massager->utmInfo.useUTMconversion)
		{
		  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
										  massager->inputgrid.xValues[j],
										  massager->inputgrid.yValues[j]);
		  Points[j][0] = pt_utm.x;
		  Points[j][1] = pt_utm.y;
		}
	  else
		{			  
		  Points[j][0] = massager->inputgrid.xValues[j];
		  Points[j][1] = massager->inputgrid.yValues[j];
		}
	}
  // build kdtree
  massager->inputgrid.OptkdRoot = BuildOptTree(Points,massager->inputgrid.nPoints,2);
  OptkdRoot = massager->inputgrid.OptkdRoot;

  // Be careful about number of neighbours
  if (nNeighborsDesired > massager->inputgrid.nPoints) nNeighborsDesired = massager->inputgrid.nPoints-1;  
  // workspace for the OK kriging 
  Ckrig = dmatrix(0,nNeighborsDesired,0,nNeighborsDesired);
  RHSkrig = dvector(0,nNeighborsDesired);
  ipiv_krig = ivector(0,nNeighborsDesired);

  // set up kriging structures for interpolation
  // loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  // set up appropriate metric for the kriging neighbour search using a kdtree
	  setUp2DEuclidenGenMetric(OptkdRoot->Kglobals,
							   M_PI*massager->variograms[i].azimuth/180.0,
							   massager->variograms[i].rangeX/massager->variograms[i].rangeY);

	  // set up variogram for kriging on this ith outputgrid
	  setUpcorr2D(&massager->variograms[i]);

 	  // Run over all outputgrid points
	  // Set up kriging weights and indices
	  for (j=0; j < massager->outputgrid[i].nPoints; j++)
		{
		  // point on output grid where we want to form the kriging estimate. In UTMs if at all possible
		  // Recall that Points is in UTMs if at all possible too
		  if (massager->utmInfo.useUTMconversion)
			{
			  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
											  massager->outputgrid[i].xValues[j],
											  massager->outputgrid[i].yValues[j]);
			  pt[0] = pt_utm.x;
			  pt[1] = pt_utm.y;
			}
		  else
			{
			  pt[0] = massager->outputgrid[i].xValues[j];
			  pt[1] = massager->outputgrid[i].yValues[j];
			}

		  // get list of outputgrid neighbours of this point from kdtree
		  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,OptkdRoot,1);

		  // form OK kriging system
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  // Matrix. Need form only lower triangle
			  for (kk=1; kk <= k ; kk++) 
				{
				  // last index is "dimension" => 0 = x , 1 = y
				  Ckrig[k-1][kk-1] = corr2D(Points[neigh[k]][0],  Points[neigh[k]][1],
											Points[neigh[kk]][0], Points[neigh[kk]][1],
											&(massager->variograms[i]));
				  Ckrig[kk-1][k-1] = Ckrig[k-1][kk-1]; // FIXME needed?
				}
			  // RHS
			  RHSkrig[k-1] = corr2D(pt[0],pt[1],Points[neigh[k]][0],Points[neigh[k]][1],&(massager->variograms[i]));				             
			}
		  // edges
		  for (kk=0; kk < nNeighborsDesired ; kk++) Ckrig[nNeighborsDesired][kk] = Ckrig[kk][nNeighborsDesired] = 1.0;
		  Ckrig[nNeighborsDesired][nNeighborsDesired]=0.0;
		  RHSkrig[nNeighborsDesired] = 1.0;

		  // Now solve kriging system for the weights. OK is not pos. def, so use general solver.  
		  // info = clapack_dpotrf(CblasColMajor,CblasLower,nNeighborsDesired,&Ckrig[0][0],nNeighborsDesired);
		  info = clapack_dgetrf(CblasColMajor,nNeighborsDesired+1,nNeighborsDesired+1,
								&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0]);
		  if (info>0) 
			{
			  warning ("lapack decomposition failed in kriging step\n");
			  return EXIT_FAILURE;
			}

		  clapack_dgetrs(CblasColMajor,CblasNoTrans,nNeighborsDesired+1,1,&Ckrig[0][0],
						 nNeighborsDesired+1,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired+1);
		  /*
			Solution (kriging weights) are now in RHSkrig[0...nNeighborsDesired-1].
		  */
		  // Store kriging weights etc.
		  // allocate room for storage of kriging weights and store
		  sI->OGkrigInfoSoft[i][j].size = nNeighborsDesired;
		  sI->OGkrigInfoSoft[i][j].values =	(struct PAIR*)malloc((size_t) (nNeighborsDesired * sizeof(struct PAIR)));
		  for (k=0; k <nNeighborsDesired; k++)
			{
			  sI->OGkrigInfoSoft[i][j].values[k].val =  RHSkrig[k];  // kriging weight
			  sI->OGkrigInfoSoft[i][j].values[k].indx = neigh[k+1];  // input grid index
			}
		  free(neigh);
		  // just for cleanliness
		  OptkdRoot->Kglobals->nndist = NULL;
		}
	  KillMetricSpace(OptkdRoot->Kglobals);
	}

  free_dmatrix(Ckrig,0,nNeighborsDesired,0,nNeighborsDesired);
  free_dvector(RHSkrig,0,nNeighborsDesired);
  free_ivector(ipiv_krig,0,nNeighborsDesired);
  
  return EXIT_SUCCESS;
}



int setUpDensePriorModel(struct MASSAGER *massager)
{
  /*
	Set up a suitably benign prior model for the properties structure
  */
  int i,j,k,output_grid_indx;
  int ii,jj,kk,output_grid_indx_kk;
  int properties_indx_k;
  int breakOut;
  double val,sum,sumSq;
  double *values;
  double **averages;
  double **variances;
  double *mu_p;   // short cut
  double **C_p;   // short cuts

  // sparse stuff
  double pt[2];
  optkdNode *OptkdRoot;
  double var,max,min;
  int i_n;
  int outputgrid_indx_offset;
  int nNeighborsDesired,nNeighInShell;
  int* neigh;





  massager->solverInfo.mu_p  = dvector(0, massager->solverInfo.nColumnsA - 1);
  // shortcut to increase legibility later on
  mu_p = massager->solverInfo.mu_p;

  massager->solverInfo.C_p  = dmatrix(0, massager->solverInfo.nColumnsA - 1,
										  0, massager->solverInfo.nColumnsA - 1);
  // shortcut to increase legibility later on
  C_p = massager->solverInfo.C_p;

  // set the mean of each property to be the average seen in each layer

  // form averages.
  averages = dmatrix(0,massager->nLayers-1,0,massager->nPropertyNames-1);
  variances = dmatrix(0,massager->nLayers-1,0,massager->nPropertyNames-1);
  // loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  for (k=0; k< massager->nPropertyNames; k++) 
		{
		  // need this to extract the right properties index from the unrolled array massager->properties
		  properties_indx_k = k + i*massager->nPropertyNames;
		  sum=sumSq=0.0;
		  values = massager->properties[properties_indx_k].values;
		  for (j=0; j < massager->inputgrid.nPoints; j++)
			{
			  val = values[j];
			  sum += val;
			  sumSq += val*val;
			}
		  averages[i][k] = sum/massager->inputgrid.nPoints;
		  variances[i][k] = (sumSq - massager->inputgrid.nPoints * averages[i][k]*averages[i][k])/(massager->inputgrid.nPoints-1);
		  // scale variances by suitable "safety factor"
		  variances[i][k] *= PRIOR_STANDARD_DEVIATION_BOOST_FACTOR * PRIOR_STANDARD_DEVIATION_BOOST_FACTOR;
		  if (variances[i][k] <=0.0) variances[i][k] = VARIANCE_EPS;
		  printf("average of property %s, layer %d= %g, var = %g\n", 
				 massager->propertyNames[k],massager->layerNumbers[i],averages[i][k], variances[i][k]);
		}
	}

  // Now populate model with these estimates
  // the prior covariance for each property comes from the variogram, scaled to the size of these variances
  // loop over standard ordering of output grid
  output_grid_indx=0;
  for (i=0; i<  massager->nLayers ; i++)
	{
	  /* Precompute a useful output grid index offset */
	  for (k = outputgrid_indx_offset = 0; k<i; k++) outputgrid_indx_offset += massager->outputgrid[k].nPoints;
	  outputgrid_indx_offset *= massager->nPropertyNames;

	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  /* loop over properties at this output grid location j*/
		  for (k=0; k< massager->nPropertyNames; k++)
			{
			  // set mean
			  mu_p[output_grid_indx] = averages[i][k];


			  // Dense filling of prior

			  // loop along comparable nesting to set covariance entries for this
			  // row (output_grid_indx) of the prior covariance Cp[][]
			  // this will run along columns 0,1.....output_grid_indx
			  output_grid_indx_kk=0;
			  breakOut=0;
			  for (ii=0; ii<  massager->nLayers ; ii++)
				{
				  for (jj=0; jj < massager->outputgrid[ii].nPoints; jj++) 
					{
					  /* loop over properties at this output grid location jj*/
					  for (kk=0; kk< massager->nPropertyNames; kk++)
						{
						  // printf("%d %d\n",output_grid_indx,output_grid_indx_kk);
						  // Set covariance entries
						  // Assume: properties in distinct layers are independent, 
						  // also: properties of different index (aka different properties) are independent
						  // Their mutual dependence is what the Bayesian update produces
						  if ((k!=kk) || (i!=ii))
							{
							  C_p[output_grid_indx][output_grid_indx_kk] = 0.0;
							}
						  else
							{

							  C_p[output_grid_indx][output_grid_indx_kk] = variances[ii][kk]
								* corr2D(massager->outputgrid[ii].xValues[j], massager->outputgrid[ii].yValues[j],
										 massager->outputgrid[ii].xValues[jj], massager->outputgrid[ii].yValues[jj],
										 &(massager->variograms[ii]));
							}
						  // symmetry
						  C_p[output_grid_indx_kk][output_grid_indx] = C_p[output_grid_indx][output_grid_indx_kk];
						  // only do up to & including the diagonal. Symmetry has done everything else
						  if (output_grid_indx_kk == output_grid_indx) 
							{
							  breakOut=1;
							  break;
							};
						  output_grid_indx_kk++;
						}
					  if (breakOut) break;
					}
				  if (breakOut) break;
				}				
			  output_grid_indx++;
			}
		}
	}
  // clean up
  free_dmatrix(averages,0,massager->nLayers-1,0,massager->nPropertyNames-1);
  free_dmatrix(variances,0,massager->nLayers-1,0,massager->nPropertyNames-1);


  // FIXME. Uncomment for ascii view of matrices
  // print_dmatrix0(C_p, massager->solverInfo.nColumnsA, massager->solverInfo.nColumnsA, "Cp.txt");

  return EXIT_SUCCESS;
}


int computePosteriorModelFromFullDenseSystems(struct MASSAGER *massager)
{
  double *mu;
  int info;
  int properties_indx_k;
  struct PT2D pt_utm;

  // ascii stuff
  if (massager->writeAsciiOut && massager->utmInfo.useUTMconversion)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}

  // luxurious output for extreme verbosity
  if (massager->verbosity>3)
	{
	  if (massager->useDenseMethod)
		{
		  printf("Cp: rows=%d, cols=%d\n",massager->solverInfo.nColumnsA,massager->solverInfo.nColumnsA);
		  print_dmatrix0(massager->solverInfo.C_p, massager->solverInfo.nColumnsA, massager->solverInfo.nColumnsA, "Cp.txt");
		  
		  printf("mu_p: dim=%d\n",massager->solverInfo.nColumnsA);
		  print_dvector0(massager->solverInfo.mu_p, massager->solverInfo.nColumnsA,"mu.txt");
		  
		  printf("A: rows=%d, cols=%d\n",massager->solverInfo.nRowsA,massager->solverInfo.nColumnsA);
		  print_dmatrix0(massager->solverInfo.A, massager->solverInfo.nRowsA, massager->solverInfo.nColumnsA, "A.txt");
		  
		  printf("C_D: rows=%d, cols=%d\n",massager->solverInfo.nRowsA,massager->solverInfo.nRowsA);
		  print_dmatrix0(massager->solverInfo.C_D,massager->solverInfo.nRowsA,massager->solverInfo.nRowsA,"CD.txt");
		  
		  printf("y: dim=%d\n",massager->solverInfo.nRowsA);
		  print_dvector0(massager->solverInfo.y, massager->solverInfo.nRowsA,"y.txt");
		}
	}

  if (massager->verbosity>3)
	{
	  // print the prior
	  mu = massager->solverInfo.mu_p;
	  int i,j,k;
	  int output_grid_indx=0;
	  for (i=0; i<  massager->nLayers ; i++)
		{
		  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
			{
			  /* loop over properties at this output grid location j*/
			  for (k=0; k< massager->nPropertyNames; k++)
				{
				  // print new mean
				  printf("prior:layer:%d %d %g\n",i,k,
						 mu[output_grid_indx]);
				  output_grid_indx++;
				}
			}
		  printf("\n");
		}

	  // print the soft observations
	  for (i=0; i<  massager->nLayers ; i++)
		{
		  for (j=0; j < massager->inputgrid.nPoints; j++)
			{
			  for (k=0; k< massager->nPropertyNames; k++)
				{
				  // need this to extract the right properties index from the unrolled array massager->properties
				  properties_indx_k = k + i*massager->nPropertyNames;
				  // observation
				  printf("obs:layer:%d %d %g %g\n",i,k,
						 massager->properties[properties_indx_k].values[j],
						 sqrt(massager->property_covariance[j][properties_indx_k][properties_indx_k]));
				}
			}
		  printf("\n");   
		}
	}

  // Solve it all in one fell swoop. Hurrah.
  /*
	Syntax of bayes library calls....
	updateBayes(double **Cp, double *mu, int Np, double **X, double **Cy, double *y, int Ny, 
	int CyIsAlreadyInverted, int CyIsIdentityMatrix)
  */			
  // call the bayes library to do what we want
  massager->solverInfo.update = updateBayes(massager->solverInfo.C_p,
												massager->solverInfo.mu_p,
												massager->solverInfo.nColumnsA, 
												massager->solverInfo.A,
												massager->solverInfo.C_D,
												massager->solverInfo.y,
												massager->solverInfo.nRowsA,
												0,0);
  if (massager->solverInfo.update.info > 0)
	{
	  printf("FAILURE: Full bayes update failed with code %d\n",massager->solverInfo.update.info);
	  return EXIT_FAILURE;
	}
  
  // print the updated mean
	  
  mu = massager->solverInfo.update.mu_tilde;
  int i,j,k,output_grid_indx_tmp;
  int output_grid_indx=0;
  for (i=0; i<  massager->nLayers ; i++)
	{
	  if (massager->writeAsciiOut) writeAsciiHeaderInfo(i,massager,0);

	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  /* loop over properties at this output grid location j*/
		  /*
			for (k=0; k< massager->nPropertyNames; k++)
			{
			printf("post:layer:%d %d %g %g\n",i,k,
			mu[output_grid_indx],
			sqrt(massager->solverInfo.update.C_tilde[output_grid_indx][output_grid_indx]));
			output_grid_indx++;
			}
		  */
		  if (massager->writeAsciiOut)
			{
			  // these are the solutions we want
			  // This is the set of properties at the central grid node...
			  if (massager->writeAsciiOut)
				{
				  // write out x,y
				  fprintf(massager->ascii_out_fp[i],"%g %g ",
						  massager->outputgrid[i].xValues[j],
						  massager->outputgrid[i].yValues[j]);
				  /* add utm stuff if need be */
				  if (massager->utmInfo.useUTMconversion)
					{
					  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
													  massager->outputgrid[i].xValues[j],
													  massager->outputgrid[i].yValues[j]);
					  fprintf(massager->ascii_out_fp[i],"%10.2f %10.2f ",pt_utm.x,pt_utm.y);
					} 
				  // write out means of all properties
				  output_grid_indx_tmp = output_grid_indx;
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  /*
						printf("layer %d node %d mu_tilde=%g stddev=%g\n",i,j,
						localUpdate.mu_tilde[k+ *hashIndx],
						sqrt(localUpdate.C_tilde[k+ *hashIndx][k+ *hashIndx]));
					  */
					  if (mu[output_grid_indx]<0)
						{
						  warning("Resetting negative property!");
						  mu[output_grid_indx] = 0;
						}
					  if (massager->properties[k].normalised && (mu[output_grid_indx]>1))
						{
						  mu[output_grid_indx] = 1;
						}
					  fprintf(massager->ascii_out_fp[i],"%g ",mu[output_grid_indx]);
					  output_grid_indx++;
					}
				  // and std deviations
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  fprintf(massager->ascii_out_fp[i],"%g ",
							  sqrt(massager->solverInfo.update.C_tilde[output_grid_indx_tmp][output_grid_indx_tmp]));
					  output_grid_indx_tmp++;
					}
				  fprintf(massager->ascii_out_fp[i],"\n");
				}
			}
		}
	  // close output files as soon as possible
	  if (massager->writeAsciiOut)
		{
		  fclose(massager->ascii_out_fp[i]);
		  printf("Wrote: %s\n",massager->ascii_out_filenames[i]);
		  free(massager->ascii_out_filenames[i]);
		}
	  if (massager->Nrealisations>0)
		{
		  fclose(massager->ascii_out_realisations_fp[i]);
		  printf("Wrote: %s\n",massager->ascii_out_realisations_filenames[i]);
		  free(massager->ascii_out_realisations_filenames[i]);
		}
	}

  return EXIT_SUCCESS;
}
  
/*
 Approximate solution of problem using scaling from smoothed covariances
*/
int solveScaledSmoothSparseProblems(struct MASSAGER *massager)
{
  int i,j,k,kp,m,mm,sizeP,output_grid_indx;
  int properties_indx_k,properties_indx_kp;
  double val,sum,sumSq,minUnscaledCorr,unscaledCorr,maxWellCov;
  double *values;

  // sparse stuff
  int n,n_y,n_i,n_j;
  int ii,jj,kk;
  int i_p,j_p;
  int i_y,j_y;
  int *hashIndx;
  int breakFlag;
  double wt,weight_sum,wt2;
  double pt[2];
  optkdNode **OptkdRoot;  // array of kdtree roots
  int nNeighInShell;
  int* neigh;
  int *indxs;
  GHashTable** collapsedOG_hashtable;  // a hash table array from output grid index to collapsed local index
  struct SOLVER_INFO* sI = &massager->solverInfo;  // Useful shortcut


  // SK kriging stuff
  double** Ckrig;
  double* RHSkrig,*Crhs;
  float ***krigVariances,***krigMeans;
  int *ipiv_krig;
  int size,IGindex;
  double var,mean,variance,krigVar,krigCorrection;
  double interpolatedSmoothPropCovar;
  struct PT2D pt_utm;
  struct PAIR* krigPairs;

  // dense entries
  double **Cp,**A,**C_D;
  double *y,*mu_p;

  // realisations stuff
  int info,i_sample,real_no;
  long seed = -1;             // FIXME promote
  double condMean,condVar,deviate,deviate_scale;
  float ****simulationValues;    // will be huge, so alloc as a float (not numerically sensitive)
  float ***posteriorMeans;       // mean for grid i,node j, property k 
  double **knownValues;          // working array of known neighbour values in local simulation machinery
  double **Lscale;
  float mu_neigh;
  int **path;
  int nNeighborsDesired;
  dynamicTree *dTree;

  if (massager->utmInfo.useUTMconversion)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}

  // ASCII IO: loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  if (massager->writeAsciiOut) writeAsciiHeaderInfo(i,massager,0);
	  if (massager->Nrealisations > 0) writeAsciiHeaderInfo(i,massager,1);
	}
  // Presmooth means and std deviations
  computeMovingAverage(massager,1);

  krigVariances = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
  krigMeans = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
  
  // precompute and store multigrid paths
  path = (int**)malloc((size_t) ((massager->nLayers * sizeof(int*))));
  for (i=0; i< massager->nLayers ; i++)
	{
	  size = massager->outputgrid[i].nPoints;
	  /* loop on all OG grid nodes in this layer */
	  struct INTPAIR* mgIndices = (struct INTPAIR*)malloc((size_t) ((size * sizeof(struct INTPAIR))));
	  formMultigridTraverseOrder(mgIndices,size);
	  // store path for ith layer
	  path[i] = (int*)malloc((size_t) ((size * sizeof(int))));
	  for (k=0; k < size; k++) path[i][k] = mgIndices[k].val;
	  free(mgIndices);
	}

  if (massager->Nrealisations > 0)
	{
	  simulationValues = (float ****)malloc((size_t) ((massager->nLayers) * sizeof(float ***)));
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  simulationValues[i] = f3tensor(0,massager->outputgrid[i].nPoints - 1,   // points
										 0,massager->nPropertyNames-1,            // props
										 0,massager->Nrealisations-1);            // realisation #
		}	  
	}


  if ((massager->verbosity>2)) printf("Starting well-kriging scaling loop...\n");
  for (i=0; i< massager->nLayers ; i++)
	{
	  /* loop on all OG grid nodes in this layer */
	  // Do pseudo multigrid by traversing binary tree breadth first???
	  size = massager->outputgrid[i].nPoints;

	  // (1) Kriging Loop
	  // Now traverse outputgrid in this order, solving SK kriging system
	  krigMeans[i] = matrix(0,size,0,massager->nPropertyNames-1);
	  krigVariances[i] = matrix(0,size,0,massager->nPropertyNames-1);

	  for (j=0; j < size; j++) 
		{
		  output_grid_indx = path[i][j];   // multigrid path
		  //output_grid_indx = j; // standard raster scan

		  // We solve a SK system for each property at this grid location
		  for (m=0; m< massager->nPropertyNames; m++)
			{
			  nNeighborsDesired = massager->well_observations[i].n_obs[m];

			  if (nNeighborsDesired>0)
				{
				  Ckrig = dmatrix(0,nNeighborsDesired-1,0,nNeighborsDesired-1);
				  RHSkrig = dvector(0,nNeighborsDesired-1);
				  Crhs = dvector(0,nNeighborsDesired-1);
				  ipiv_krig = ivector(0,nNeighborsDesired-1);

				  // form SK kriging system for all wells about this node
				  for (k=0; k < nNeighborsDesired; k++) 
					{
					  // Matrix. Need form only lower triangle
					  // well observations are always in absolute coordinates
					  for (kk=0; kk <= k ; kk++) 
						{
						  Ckrig[k][kk] = corr2D(massager->well_observations[i].obs[m][k].x, massager->well_observations[i].obs[m][k].y,
												massager->well_observations[i].obs[m][kk].x, massager->well_observations[i].obs[m][kk].y,
												&(massager->variograms[i]));
						  // normalise to guarantee unit variance
						  Ckrig[k][kk] /= Ckrig[0][0];
						  Ckrig[kk][k] = Ckrig[k][kk]; // FIXME needed?
						}
					  // RHS. output grid values may be in seismic coordinates....... Convert to UTMs if so
					  if (massager->utmInfo.useUTMconversion)
						{
						  pt_utm = getUTMSfromInlineXline(massager->utmInfo,
														  massager->outputgrid[i].xValues[output_grid_indx],
														  massager->outputgrid[i].yValues[output_grid_indx]);

						  RHSkrig[k] = corr2D(pt_utm.x, pt_utm.y,
											  massager->well_observations[i].obs[m][k].x, massager->well_observations[i].obs[m][k].y,
											  &(massager->variograms[i]));				             
						}
					  else
						{
						  RHSkrig[k] = corr2D(massager->outputgrid[i].xValues[output_grid_indx],massager->outputgrid[i].yValues[output_grid_indx],
											  massager->well_observations[i].obs[m][k].x, massager->well_observations[i].obs[m][k].y,  
											  &(massager->variograms[i]));				             
						}
					  RHSkrig[k] /= Ckrig[0][0]; // normalise to guarantee unit variance
					}
				  memcpy(Crhs,RHSkrig,nNeighborsDesired*sizeof(double));

				  // Now solve kriging system for the weights mapping the output grid to the well-observation point
				  // print_dmatrix0(Ckrig,nNeighborsDesired+1,"junk.txt");
				  info = clapack_dgetrf(CblasColMajor,nNeighborsDesired,nNeighborsDesired,
										&Ckrig[0][0],nNeighborsDesired,&ipiv_krig[0]);
				  if (info>0)
					{
					  warning ("lapack decomposition failed in kriging step\n");
					  return EXIT_FAILURE;
					}   
				  
				  /* Get weights via dge solution solution */
				  clapack_dgetrs(CblasColMajor,CblasNoTrans,nNeighborsDesired,1,&Ckrig[0][0],
								 nNeighborsDesired,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired);
					  
				  /*
					Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1].
					Compute kriged mean and variance
				  */

				  mean = massager->solverInfo.smoothMeans[i][output_grid_indx][m];
				  krigCorrection=0;

				  krigVar = 1.0; // =Ckrig[0][0] ... will be 1
				  for (k=0; k <nNeighborsDesired; k++)
					{
					  // residual at data point: construct from interpolation
					  // interpolate smoothedMeans to Well observation point
					  krigCorrection += RHSkrig[k] * (massager->well_observations[i].obs[m][k].val - 
													  massager->well_observations[i].obs[m][k].smoothVal);
					  krigVar -= RHSkrig[k]*Crhs[k]; // standard SK kriging variance: GSLIB II.14
					}
				  mean += krigCorrection;

				  free_dmatrix(Ckrig,0,nNeighborsDesired-1,0,nNeighborsDesired-1);
				  free_dvector(RHSkrig,0,nNeighborsDesired-1);
				  free_dvector(Crhs,0,nNeighborsDesired-1);
				  free_ivector(ipiv_krig,0,nNeighborsDesired-1);
				}
			  else
				{
				  mean = massager->solverInfo.smoothMeans[i][output_grid_indx][m];
				  krigVar = massager->variograms[i].sill + massager->variograms[i].nugget;
				}
			  krigMeans[i][output_grid_indx][m] = mean;
			  krigVariances[i][output_grid_indx][m] = krigVar;
			  /* krigVariances[i][output_grid_indx][m] is now normalised kriging variance of property m, conditioned only to wells */
			  /* We will have krigVariances[i][output_grid_indx][m] -> 1 far from the wells */
#if 0
			  // All this subsumed below........
			  /* Now scale variance by smooth covariance interpolated from input grid */
			  properties_indx_k = m + i*massager->nPropertyNames;
			  // interpolate smoothpropertyCovariance from input grid to this output grid site, using precomputed kriging weights
			  smoothPropCovar[m] = 0;
			  krigPairs = sI->OGkrigInfoSoft[i][output_grid_indx].values;
			  for (k=0; k< sI->OGkrigInfoSoft[i][output_grid_indx].size; k++)
				{
				  IGindex = krigPairs[k].indx;
				  wt = krigPairs[k].val;
				  smoothPropCovar[m] += massager->solverInfo.smoothPropertyCovariance[IGindex][properties_indx_k][properties_indx_k];
				}
			  // scale this by the well-kriged variance to recover true variance
			  krigVariances[i][output_grid_indx][m] *= smoothPropCovar[m];
#endif
			}
		}

	  // Kriging I/O
	  if (massager->writeAsciiOut) 
		{
		  // write out surfaces in natural order
		  for (j=0; j < size; j++) 
			{
			  // write out x,y
			  fprintf(massager->ascii_out_fp[i],"%g %g ",
					  massager->outputgrid[i].xValues[j],
					  massager->outputgrid[i].yValues[j]);
			  /* add utm stuff if need be */
			  if (massager->utmInfo.useUTMconversion)
				{
				  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
												  massager->outputgrid[i].xValues[j],
												  massager->outputgrid[i].yValues[j]);
				  fprintf(massager->ascii_out_fp[i],"%10.2f %10.2f ",pt_utm.x,pt_utm.y);
				}
			  // run over properties in this layer
			  for (m=0; m< massager->nPropertyNames; m++) 
				{
				  // write out kriged mean
				  fprintf(massager->ascii_out_fp[i],"%g ",krigMeans[i][j][m]);
				  /* Compute variance: scale variance by smooth covariance interpolated from input grid */
				  properties_indx_k = m + i*massager->nPropertyNames;
				  // interpolate smoothpropertyCovariance from input grid to this output grid site, using precomputed kriging weights
				  interpolatedSmoothPropCovar = 0;
				  krigPairs = sI->OGkrigInfoSoft[i][j].values;
				  for (k=0; k< sI->OGkrigInfoSoft[i][j].size; k++)
					{
					  IGindex = krigPairs[k].indx;
					  wt = krigPairs[k].val;
					  interpolatedSmoothPropCovar += massager->solverInfo.smoothPropertyCovariance[IGindex][properties_indx_k][properties_indx_k];
					}
				  // scale this by the well-kriged variance to recover true variance
				  variance = krigVariances[i][j][m] * interpolatedSmoothPropCovar;
				  // write out kriged variance
				  fprintf(massager->ascii_out_fp[i],"%g ",sqrt(MAX(variance,0)));
				  // fprintf(massager->ascii_out_fp[i],"%g ",sqrt(MAX(krigVariances[i][j][m],0)));
				}
			  fprintf(massager->ascii_out_fp[i],"\n");
			}
		  // close output files as soon as possible
		  fclose(massager->ascii_out_fp[i]);
		  printf("Wrote: %s\n",massager->ascii_out_filenames[i]);
		  free(massager->ascii_out_filenames[i]);
		}
	}

 
  // (2) Realisations Loop
  if (massager->Nrealisations>0)
	{
	  if ((massager->verbosity>2)) printf("Starting realisations loop...\n");
	  for (i=0; i< massager->nLayers ; i++)
		{
		  size = massager->outputgrid[i].nPoints;
		  
		  // Dynamic KD tree for output grid if realisations required
		  if (massager->Nrealisations > 0)
			{
			  dTree = buildDynamicTree(massager->outputgrid[i].Points, size, 2, path[i]);
			  // set up appropriate metric for the kriging neighbour search using the dynamic kdtree, 
			  // suitable for ith output grid
			  setUp2DEuclidenGenMetric(dTree->OptkdRoot->Kglobals,
									   M_PI*massager->variograms[i].azimuth/180.0,
									   massager->variograms[i].rangeX/massager->variograms[i].rangeY);
			}

		  for (j=0; j < size; j++) 
			{
			  if (massager->verbosity>2)
				{
				  if (j%1000==0) printf("Layer %d: Solving node %d...\n",massager->layerNumbers[i],j);
				}
			  output_grid_indx = path[i][j];   // multigrid path

			  // Construct SK system from previously simulated nodes
			  // We solve a SK system for each property at this grid location
			  // The structure of this system is the same for all properties, since
			  // the accumulating neighbours are the same for all properties, the path is the same etc,
			  // and the variogram is taken to be the same
			  nNeighborsDesired = MIN(j,DEFAULT_MAX_NEIGHBOURS_IN_SIMULATIONS);
		  
			  // find previsited nearest neighbours
			  neigh = nextDynamicKdOptNNQuery(nNeighborsDesired, EUCLIDEAN_GEN, dTree);
			  // neigh[1....nNeighborsDesired] are now the indexes 
			  // of the already simulated points on massager->outputgrid[i]
			  if (nNeighborsDesired>0)
				{
				  Ckrig = dmatrix(0,nNeighborsDesired-1,0,nNeighborsDesired-1);
				  RHSkrig = dvector(0,nNeighborsDesired-1);
				  Crhs = dvector(0,nNeighborsDesired-1);
				  ipiv_krig = ivector(0,nNeighborsDesired-1);

				  // form SK kriging system over these
				  for (k=0; k < nNeighborsDesired; k++) 
					{
					  // Matrix. Need form only lower triangle
					  for (kk=0; kk <= k ; kk++) 
						{
						  // massager->outputgrid[i].Points is in UTM's if at all possible
						  Ckrig[k][kk] = corr2D(massager->outputgrid[i].Points[neigh[k+1]][0], massager->outputgrid[i].Points[neigh[k+1]][1],
												massager->outputgrid[i].Points[neigh[kk+1]][0], massager->outputgrid[i].Points[neigh[kk+1]][1],
												&(massager->variograms[i]));
						  Ckrig[k][kk] /= Ckrig[0][0]; // normalise to guarantee unit variance
						  Ckrig[kk][k] = Ckrig[k][kk]; // FIXME needed?
						}
					  // RHS.
					  RHSkrig[k] = corr2D(massager->outputgrid[i].Points[output_grid_indx][0],massager->outputgrid[i].Points[output_grid_indx][1],
										  massager->outputgrid[i].Points[neigh[k+1]][0], massager->outputgrid[i].Points[neigh[k+1]][1],  
										  &(massager->variograms[i]));
					  RHSkrig[k] /= Ckrig[0][0]; // normalise to guarantee unit variance
					}
				  memcpy(Crhs,RHSkrig,nNeighborsDesired*sizeof(double));
				  condVar = 1.0;  // store this before we mes up Ckrig

				  // Now solve kriging system for the weights
				  info = clapack_dgetrf(CblasColMajor,nNeighborsDesired,nNeighborsDesired,
										&Ckrig[0][0],nNeighborsDesired,&ipiv_krig[0]);
				  if (info>0)
					{
					  warning ("lapack decomposition failed in kriging step\n");
					  return EXIT_FAILURE; 
					}   
				  
				  /* Get weights via dge solution solution */
				  clapack_dgetrs(CblasColMajor,CblasNoTrans,nNeighborsDesired,1,&Ckrig[0][0],
								 nNeighborsDesired,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired);
					  
				  /*
					Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1].
					Generate unconditional realisations with normalised variogram
				  */
				  for (k=0; k <nNeighborsDesired; k++) condVar -= RHSkrig[k]*Crhs[k]; // standard SK kriging variance: GSLIB II.14
				  double deviate_scale = sqrt(MAX(condVar,0));
				  for (m=0; m< massager->nPropertyNames; m++)
					{
					  //mean = krigMeans[i][output_grid_indx][m];
					  //for (k=0; k <nNeighborsDesired; k++) mean -= RHSkrig[k] * krigMeans[i][neigh[k+1]][m];
					  for (real_no=0; real_no < massager->Nrealisations; real_no++)
						{
						  condMean = 0;
						  for (k=0; k <nNeighborsDesired; k++)
							{
							  // residual at data point: construct from interpolation
							  // interpolate smoothedMeans to Well observation point
							  condMean += RHSkrig[k] * simulationValues[i][neigh[k+1]][m][real_no];
							}
						  //condMean += mean;
						  deviate = deviate_scale * gasdev(&seed);
						  // Generate realisations and store
						  simulationValues[i][output_grid_indx][m][real_no] = condMean + deviate;
						}
					}
				  free_dmatrix(Ckrig,0,nNeighborsDesired-1,0,nNeighborsDesired-1);
				  free_dvector(RHSkrig,0,nNeighborsDesired-1);
				  free_dvector(Crhs,0,nNeighborsDesired-1);
				  free_ivector(ipiv_krig,0,nNeighborsDesired-1);
				}
			  else
				{
				  // first point: no neighbours
				  deviate_scale = 1;
				  for (m=0; m< massager->nPropertyNames; m++)
					{
					  mean = 0;
					  for (real_no=0; real_no < massager->Nrealisations; real_no++)
						{
						  deviate = deviate_scale * gasdev(&seed);
						  // store
						  simulationValues[i][output_grid_indx][m][real_no] = mean + deviate;
						}
					}
				}
			}
		  // Cleanup dynamic tree
		  killDynamicTree(dTree); 
		}

	  // Realisations I/O
	  for (i=0; i< massager->nLayers ; i++)
		{	 
		  for (j=0;j<massager->outputgrid[i].nPoints;j++)
			{
			  // write out x,y
			  fprintf(massager->ascii_out_realisations_fp[i],"%g %g ",
					  massager->outputgrid[i].xValues[j],
					  massager->outputgrid[i].yValues[j]);
			  /* add utm stuff if need be */
			  if (massager->utmInfo.useUTMconversion)
				{
				  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
												  massager->outputgrid[i].xValues[j],
												  massager->outputgrid[i].yValues[j]);
				  fprintf(massager->ascii_out_realisations_fp[i],"%10.2f %10.2f ",pt_utm.x,pt_utm.y);
				}

			  // Now scale to couple in property correlations
			  Lscale = getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(massager,i,j);
			  krigVariances[i][j][m];
			  for (real_no=0; real_no < massager->Nrealisations; real_no++)
				{
				  for (m=0; m< massager->nPropertyNames; m++)
					{
					  mean = 0 * krigMeans[i][j][m]; // FIXME No 0 !!!
					  properties_indx_k = m + i*massager->nPropertyNames;
					  sum=0;
					  // Multiply across row of Cholesky prefactor. Looks messy but that's all it is
					  for (ii=0; ii< massager->nLayers ; ii++)
						{
						  for (mm=0; mm< massager->nPropertyNames; mm++)
							{
							  k = mm + ii*massager->nPropertyNames;   // properties_indx_kk
							  if (k <= properties_indx_k)
								{
								  sum += Lscale[properties_indx_k][k] * simulationValues[ii][j][mm][real_no];
								}
							}
						}
					  sum *= sqrt(MAX(krigVariances[i][j][m],0));   // HACK to get well data right
					  // the realisation is now sum+ krigMeans[i][j][m]
					  val = sum + krigMeans[i][j][m];
					  // val = simulationValues[i][j][m][real_no];
					  // Properties have to be positive. Truncate at 0
					  fprintf(massager->ascii_out_realisations_fp[i],"%12.3f ",MAX(0,val));
					}
				}
			  fprintf(massager->ascii_out_realisations_fp[i],"\n");
			  free_dmatrix(Lscale,0,massager->nProperties-1,0,massager->nProperties-1);
			}
		  fclose(massager->ascii_out_realisations_fp[i]);
		  printf("Wrote: %s\n",massager->ascii_out_realisations_filenames[i]);
		  free(massager->ascii_out_realisations_filenames[i]);
		}
	}

  // Cleanup
  // 1 smoothed means and variances
  for (i=0; i< massager->nLayers ; i++) 
	{
	  free_ivector(path[i],0,size-1);
	  free_matrix(krigMeans[i],0,size,0,massager->nPropertyNames-1);
	  free_matrix(krigVariances[i],0,size,0,massager->nPropertyNames-1);
	  free_matrix(sI->smoothMeans[i], 0, massager->outputgrid[i].nPoints-1, 0, massager->nPropertyNames-1);
	}
  free(sI->smoothMeans);
  free_f3tensor(massager->solverInfo.smoothPropertyCovariance, 0,(massager->inputgrid.nPoints)-1,
				0,(massager->nProperties)-1,0,(massager->nProperties)-1);

  if (massager->Nrealisations > 0)
	{
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  free_f3tensor(simulationValues[i],0,massager->outputgrid[i].nPoints - 1,0,massager->nPropertyNames-1,0,massager->Nrealisations-1);
		}
	  free(simulationValues);
	}


  return EXIT_SUCCESS;
}

double **getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(struct MASSAGER *massager, 
																 int layer_index, int output_grid_index)
{
  int i,j,k,m,mm,info,IGindex,N;
  double wt,eps,max;
  struct PAIR* krigPairs;
  struct SOLVER_INFO* sI;

  i = layer_index;
  j = output_grid_index;
  N = massager->nProperties;
  sI = &massager->solverInfo;  // Useful shortcut

  // Full covariance scaling
  double **interpolatedPropCovar = dmatrix(0,N-1,0,N-1);

  // interpolate the input grid full covariance into the current output-grid point
  // linear combination implicit in kriging guarantees positive definiteness
  for (m=0; m<N; m++)
	{
	  for (mm = 0; mm<=m; mm++)
		{
		  // interpolate smooth propertyCovariance from 
		  // input grid to this output grid site, using precomputed kriging weights
		  interpolatedPropCovar[m][mm] = 0;
		  krigPairs = sI->OGkrigInfoSoft[i][j].values;
		  for (k=0; k< sI->OGkrigInfoSoft[i][j].size; k++)
			{
			  IGindex = krigPairs[k].indx;
			  wt = krigPairs[k].val;
			  interpolatedPropCovar[m][mm]  += wt * sI->smoothPropertyCovariance[IGindex][m][mm];
			}
		  // symmetry
		  interpolatedPropCovar[mm][m] = interpolatedPropCovar[m][mm]; 
		}
	}
  // Now get eigenfactor U such that U.U^T = interpolatedPropCovar
  // Ansatz: Cholesky is a decent surrogate
  // Make this robust against zeros:
  max=0;
  for (m=0; m<N; m++) max = MAX(max,interpolatedPropCovar[m][m]);
  eps = 1e-6 * max;
  for (m=0; m<N; m++) interpolatedPropCovar[m][m] += eps;

  info = clapack_dpotrf(CblasRowMajor,CblasLower,N,&interpolatedPropCovar[0][0],N);
  
  if (info>0) 
	{
	  warning ("lapack Cholesky decomposition failed in interpolatedPropCovar()\n");
	  return NULL;
	}
  // free_dmatrix(interpolatedPropCovar,0,N-1,0,N-1);
  return interpolatedPropCovar;
}

void formMultigridTraverseOrder(struct INTPAIR* mgIndices, int size)
{
  int j,k,v,sizeP,kstart,level;

  v = size-1;
  sizeP=0;
  while ((v >>= 1) > 0 ) sizeP++;
  // sizeP is now the smallest power of 2 s.t. 2^sizeP >= size
  kstart=(2 << sizeP);   // 2^sizeP
  for (j=0; j < size; j++) 
	{
	  k=kstart;
	  v=j+1;
	  level=sizeP+1;
	  // keep dividing by 2 to find what level this index is in a binary tree
	  while (v%k != 0) {k >>= 1;level--;}
	  mgIndices[j].indx = level;
	  mgIndices[j].val = j;
	  // printf("powers %d %d %d\n",kstart,j,level);
	}
  // sort index list into level-based order
  qsort(&mgIndices[0],(size_t)size,sizeof(struct INTPAIR),intPairCompare);
  /*
	for (j=0; j < size; j++) 
	{
	  printf("powers %d %d\n",mgIndices[j].val,mgIndices[j].indx);
	}
  */
}


/*
  Here, we form the local dense approximate problems using local searching at each node, without
  using the sparse matrices
*/
int solveSequentialSparseProblems(struct MASSAGER *massager)
{
  // averages
  double **averages;
  double **variances;
  int i,j,k,kp;
  int properties_indx_k,properties_indx_kp;
  double val,sum,sumSq,minUnscaledCorr,unscaledCorr,maxWellCov;
  double *values;

  // sparse stuff
  int n,n_y,n_i,n_j;
  int ii,jj,kk;
  int i_p,j_p;
  int i_y,j_y;
  int *hashIndx;
  int breakFlag;
  double wt,weight_sum,wt2;
  double pt[2];
  optkdNode **OptkdRoot;  // array of kdtree roots
  int* nNeighborsDesired;
  int nNeighInShell;
  int** neigh;
  int *indxs;
  GHashTable** collapsedOG_hashtable;  // a hash table array from output grid index to collapsed local index
  struct SOLVER_INFO* sI = &massager->solverInfo;  // Useful shortcut
  struct BAYES_UPDATE localUpdate;
  struct PT2D pt_utm;

  // dense entries
  double **Cp,**A,**C_D;
  double *y,*mu_p;

  // realisations stuff
  int info,i_sample;
  float ****simulationValues; // will be huge, so alloc as a float (not numerically sensitive)
  float ***posteriorMeans;    // mean for grid i,node j, property k 
  double **knownValues;        // working array of known neighbour values in local simulation machinery
  float mu_neigh;
  long seed=-1;

  // ascii stuff
  if (massager->writeAsciiOut && massager->utmInfo.useUTMconversion)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}

  // set the mean of each property to be the average seen in each layer
  // form averages.
  averages = dmatrix(0,massager->nLayers-1,0,massager->nPropertyNames-1);
  variances = dmatrix(0,massager->nLayers-1,0,massager->nPropertyNames-1);
  // loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  if (massager->writeAsciiOut) writeAsciiHeaderInfo(i,massager,0);
	  if (massager->Nrealisations > 0) writeAsciiHeaderInfo(i,massager,1);
	}
  // Presmooth to estimate mean if required: if picking out local trend
  if (massager->useMAtrend)
	{
	  computeMovingAverage(massager,0);
	}
  else
	{
	  // simple very rough guesses (constants)
	  for (i=0; i<  massager->nLayers ; i++)
		{
		  for (k=0; k< massager->nPropertyNames; k++) 
			{
			  // need this to extract the right properties index from the unrolled array massager->properties
			  properties_indx_k = k + i*massager->nPropertyNames;
			  sum=sumSq=0.0;
			  values = massager->properties[properties_indx_k].values;
			  for (j=0; j < massager->inputgrid.nPoints; j++)
				{
				  val = values[j];
				  sum += val;
				  sumSq += val*val;
				}
			  averages[i][k] = sum/massager->inputgrid.nPoints;
			  variances[i][k] = (sumSq - massager->inputgrid.nPoints * averages[i][k]*averages[i][k])/(massager->inputgrid.nPoints-1);
			  // scale variances by suitable "safety factor"
			  variances[i][k] *= PRIOR_STANDARD_DEVIATION_BOOST_FACTOR * PRIOR_STANDARD_DEVIATION_BOOST_FACTOR;
			  // safety for zero properties
			  if (variances[i][k] <=0.0) variances[i][k] = VARIANCE_EPS;
			  printf("average of property %s, layer %d= %g, var = %g\n", 
					 massager->propertyNames[k],massager->layerNumbers[i],averages[i][k], variances[i][k]);
			}
		}
	}

  // Now loop over grids and outputgrid nodes
  int output_grid_indx=0;
  // array of ptrs to neighbour lists, for each layer
  neigh = (int**)malloc((size_t) ((massager->nLayers) * sizeof(int*)));
	  
  nNeighborsDesired = ivector(0, (massager->nLayers) -1); 

  // useful indexing array for local neighbourhoods
  // allocate from biggest possible size
  int nIndxs = massager->nPropertyNames * DEFAULT_MAX_NEIGHBOURS_IN_SPARSE_PRIOR * massager->nLayers;
  indxs = ivector(0,nIndxs);
  for (i=0; i<=nIndxs; i++) indxs[i]=i;
  // array of local output grid hash tables for each layer FIXME free me
  collapsedOG_hashtable = (GHashTable**)malloc((size_t) ((massager->nLayers) * sizeof(GHashTable*)));
  GTree** OGindxTreeWells = (GTree**)malloc((size_t) ((massager->nLayers) * sizeof(GTree*))); // FIXME free me
  GTree** OGindxTreeSoft = (GTree**)malloc((size_t) ((massager->nLayers) * sizeof(GTree*))); // FIXME free me


  // realisations
  if (massager->Nrealisations > 0)
	{
	  sI->simFlag = (int ***)malloc((size_t) ((massager->nLayers) * sizeof(int **)));
	  simulationValues = (float ****)malloc((size_t) ((massager->nLayers) * sizeof(float ***)));
	  posteriorMeans = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  sI->simFlag[i] = imatrix(0,massager->outputgrid[i].nPoints - 1,0,massager->nPropertyNames-1);
		  memset(&sI->simFlag[i][0][0],0,massager->outputgrid[i].nPoints * massager->nPropertyNames * sizeof(int));
		  simulationValues[i] = f3tensor(0,massager->outputgrid[i].nPoints - 1,   // points
										 0,massager->nPropertyNames-1,            // props
										 0,massager->Nrealisations-1);            // realisation #
		  posteriorMeans[i] = matrix(0,massager->outputgrid[i].nPoints - 1,       // points
									 0,massager->nPropertyNames-1);               // props
		}	  
	}

  if ((massager->verbosity>2)) printf("Starting sequential loop...\n");

  for (i=0; i< massager->nLayers ; i++)
	{
	  /* loop on all OG grid nodes in this layer */
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  if (massager->verbosity>2) 
			{
			  if (j%((massager->verbosity>4)?1:100) == 0) printf("Solving layer %d, node %d...\n",i,j);
			}
		  // x,y of OG point in this layer
		  pt[0] = massager->outputgrid[i].xValues[j];
		  pt[1] = massager->outputgrid[i].yValues[j];
		  n=0;
		  // form list of neighbours to this x,y point in all layers 
		  for (ii=0; ii<  massager->nLayers ; ii++)
			{
			  // get list of outputgrid neighbours of this point from kdtree
			  // Initial crude form: get nearest neighbours
			  nNeighborsDesired[ii] = MIN(massager->outputgrid[ii].nPoints,
										  DEFAULT_MAX_NEIGHBOURS_IN_SPARSE_PRIOR); 
			  // get neighbour list in ii-th grid from iith-kdtree
			  neigh[ii] = kdOptNNQuery(massager->outputgrid[ii].Points,
									   2,&pt[0],nNeighborsDesired[ii],
									   EUCLIDEAN_GEN,massager->outputgrid[ii].OptkdRoot,1);
			  // sort neighbour list into increasing order
			  qsort(&neigh[ii][1],(size_t)nNeighborsDesired[ii],sizeof(int),intCompare);
			  n += nNeighborsDesired[ii];  // accumulate size of Gaussian prior system
			}
		  /*
		  for (k=1;k<=n;k++)
			{
			  printf("NODE %d %g %g\n",j,massager->outputgrid[i].xValues[neigh[i][k]],
					 massager->outputgrid[i].yValues[neigh[i][k]]);
			}
		  */



		  n*= massager->nPropertyNames;  // size of local system embraces properties, grids, nodes

		  // Allocate prior matrix and mean
		  // Local system is union of all nodes in all layer-neighborhoods, for all properties
		  Cp = dmatrix(0,n-1,0,n-1);
		  mu_p = dvector(0,n-1);
		  // Now run over output nodes, setting prior Cp and mu as we go. i_p is main row index of local Linear system
		  minUnscaledCorr=DBL_MAX;
		  for (ii=0,i_p=0; ii<  massager->nLayers ; ii++)           // layers
			{
			  double** Points_ii = massager->outputgrid[ii].Points;
			  // and form a hash table mapping output grid indices for this layer to the collapsed local model index
			  collapsedOG_hashtable[ii] = g_hash_table_new(g_int_hash,g_int_equal);
			  for (n_i=1; n_i <= nNeighborsDesired[ii]; n_i++) // neighbours
				{
				  // Make hash table map output grid indexes to collapsed indexes, as if for property k=0
				  g_hash_table_insert(collapsedOG_hashtable[ii],(gpointer)&neigh[ii][n_i],(gpointer)&indxs[i_p]);
				  /*
				  printf("NODE %d %g %g\n",j,massager->outputgrid[ii].xValues[neigh[ii][n_i]],
					 massager->outputgrid[ii].yValues[neigh[ii][n_i]]);
				  */
				  /* printf("Hash insert layer %d, OGindx = %d, reduced = %d\n", ii,neigh[ii][n_i], indxs[i_p]); */
				  for (k=0; k< massager->nPropertyNames; k++,i_p++) // properties
					{
					  if (massager->useMAtrend)
						{
						  mu_p[i_p] = massager->solverInfo.smoothMeans[ii][neigh[ii][n_i]][k];
						}
					  else
						{
						  mu_p[i_p] = averages[ii][k]; 
						}
					  // now repeat nested loop for covariance off-diagonals 
					  breakFlag=0;
					  for (jj=0,j_p=0; jj<massager->nLayers ; jj++)           // layers
						{
						  double** Points_jj = massager->outputgrid[jj].Points;
						  for (n_j=1; n_j <= nNeighborsDesired[jj]; n_j++) // neighbours
							{
							  for (kk=0; kk< massager->nPropertyNames; kk++,j_p++) // properties
								{
								  // prior correlated only if same layer, same property
								  if (ii==jj && k==kk)
									{
									  unscaledCorr = corr2D(Points_ii[neigh[ii][n_i]][0], Points_ii[neigh[ii][n_i]][1],
															Points_jj[neigh[jj][n_j]][0], Points_jj[neigh[jj][n_j]][1],
															&(massager->variograms[ii]));
									  Cp[i_p][j_p] = variances[ii][k] * unscaledCorr;
									  minUnscaledCorr = MIN(minUnscaledCorr,fabs(unscaledCorr));
									}
								  else
									{
									  Cp[i_p][j_p] = 0.0;
									}
								  // Symmetry
								  Cp[j_p][i_p] = Cp[i_p][j_p];
								  // efficiency
								  if (j_p == i_p) {breakFlag=1;break;} 
								}
							  if (breakFlag) break; 
							}
						  if (breakFlag) break; 
						}
					} // ...nested loop for covariance off-diagonals
				}
			}
		  // print_dmatrix0(Cp,n,n,"Cp.txt");

		  // Now form observations y and C_D for this node 
		  // Accumulate list of pertinent observations by running over all observations
		  // linked to nodes in the local OG neighbourhoods, using a tree

		  n_y=0;
		  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
			{
			  OGindxTreeWells[ii] = g_tree_new(intCompare);
			  OGindxTreeSoft[ii] = g_tree_new(intCompare);

			  // loop over local suite of output grid points involved in
			  // condensed system, and find which observations are linked to these
			  // via kriging
			  for (n_i=1; n_i <= nNeighborsDesired[ii]; n_i++) // neighbours in layer ii
				{
				  // absolute outputgrid index of this n_i-th neighbour
				  int OGindx = neigh[ii][n_i];
				  // loop over kriging contributors to observations at this point and assemble their
				  // indices in an ordered tree
				  // (1) first, wells.  This loop will probably be mainly null (no entries)
				  for (n_j=0; n_j < sI->OGkrigInfoWell[ii][OGindx].size; n_j++)
					{
					  // Well observation index of an observation contributed to (via kriging weights) by this
					  // n_i-th neighbour. We get a contribution for each distinct property
					  int WellObsIndex = sI->OGkrigInfoWell[ii][OGindx].values[n_j];
					  // printf("WELLOBSCONTRIB %d %d %g %g\n",WellObsIndex,OGindx,massager->outputgrid[ii].xValues[neigh[ii][n_i]],massager->outputgrid[ii].yValues[neigh[ii][n_i]]);
					  
					  // add to tree
					  g_tree_insert(OGindxTreeWells[ii], (gpointer)&sI->OGkrigInfoWell[ii][OGindx].values[n_j],
									(gpointer)&sI->OGkrigInfoWell[ii][OGindx].values[n_j]);		 
					} 
				  // ( 2) Then soft data
				  for (n_j=0; n_j < sI->OGkrigLookupInfoSoft[ii][OGindx].size; n_j++)
					{ 
					  // Absolute soft observation index of an observation contributed to (via kriging weights) by this
					  // n_i-th neighbour. We get only one entry per input grid point.
					  int IGindxNeigh = sI->OGkrigLookupInfoSoft[ii][OGindx].values[n_j];		
					  //  add to tree
					  g_tree_insert(OGindxTreeSoft[ii], (gpointer)&sI->OGkrigLookupInfoSoft[ii][OGindx].values[n_j],
									(gpointer)&sI->OGkrigLookupInfoSoft[ii][OGindx].values[n_j]);
					}
				}
			  // y size is sum over wells + soft info over all layers
			  n_y += g_tree_nnodes(OGindxTreeWells[ii]) + massager->nPropertyNames*g_tree_nnodes(OGindxTreeSoft[ii]);
			}

		  // Now convert the tree into more directly useful index arrays
		  struct INDXTRIPLE *treeOrdererWells,*treeOrdererSoft;
		  treeOrdererWells = (struct INDXTRIPLE*)malloc((size_t) ((massager->nLayers) * sizeof(struct INDXTRIPLE)));
		  treeOrdererSoft = (struct INDXTRIPLE*)malloc((size_t) ((massager->nLayers) * sizeof(struct INDXTRIPLE)));
		  // Create ordering with all well data first
		  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
			{
			  treeOrdererWells[ii].indx=0;
			  treeOrdererWells[ii].size = g_tree_nnodes(OGindxTreeWells[ii]);
			  treeOrdererWells[ii].values=ivector(0,treeOrdererWells[ii].size - 1);
			  // printf("Well tree %d:\n",ii);
			  g_tree_foreach(OGindxTreeWells[ii], fillIndexTripleArray, (gpointer)&treeOrdererWells[ii]);
			}
		  // then soft data
		  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
			{
			  treeOrdererSoft[ii].indx=0;
			  treeOrdererSoft[ii].size = g_tree_nnodes(OGindxTreeSoft[ii]);
			  treeOrdererSoft[ii].values=ivector(0, treeOrdererSoft[ii].size - 1);
			  //printf("Soft tree %d:\n",ii);
			  g_tree_foreach(OGindxTreeSoft[ii], fillIndexTripleArray, (gpointer)&treeOrdererSoft[ii]);
			}
		  // Now we can assemble y and C_D
		  if (n_y>0)
			{
			  C_D = dmatrix(0,n_y-1,0,n_y-1);
			  y = dvector(0,n_y-1);
			  // zero C_D
			  memset(&C_D[0][0],0,n_y*n_y * sizeof(double));
			  i_y=0;                 // "observation" index
			  // Well data. This loop is often null
			  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
				{
				  for (n_j=0; n_j < treeOrdererWells[ii].size; n_j++)
					{
					  int WellObsIndex = treeOrdererWells[ii].values[n_j];
					  y[i_y] = massager->well_observations[ii].values[WellObsIndex];
					  C_D[i_y][i_y] = HARD_OBSERVATIONS_EPS*HARD_OBSERVATIONS_EPS*((y[i_y]!=0) ? y[i_y]*y[i_y] : 1.0);
					  // printf("WELLS %d y[%d]=%g\n",j,i_y,y[i_y]);
					  i_y++;
					}
				}
			  // Soft data
			  int inputGridIndex_ii, inputGridIndex_jj;
			  int offset = i_y;
			  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
				{
				  for (n_i=0; n_i < treeOrdererSoft[ii].size; n_i++)
					{
					  inputGridIndex_ii = treeOrdererSoft[ii].values[n_i];
					  for (k=0; k<massager->nPropertyNames; k++)
						{
						  // need this to extract the right properties index from the unrolled array massager->properties
						  properties_indx_k = k + ii*massager->nPropertyNames;
						  y[i_y] = massager->properties[properties_indx_k].values[inputGridIndex_ii];
						  /* test: only valid for -d option
							 if (inputGridIndex_ii==j) 
							{
							  printf("YASS %g %g %g %g\n",massager->outputgrid[ii].xValues[j],massager->outputgrid[ii].yValues[j],
									 y[i_y], massager->solverInfo.smoothMeans[ii][j][k]); 
							}
						  */
						  // mirrored loop for C_D entries
						  j_y = offset;
						  for (jj=0; jj<  massager->nLayers ; jj++)  // layers
							{
							  for (n_j=0; n_j < treeOrdererSoft[jj].size; n_j++)
								{
								  inputGridIndex_jj = treeOrdererSoft[jj].values[n_j];
								  for (kp=0; kp<massager->nPropertyNames; kp++)
									{
									  // entries only for properties from same input grid "trace"
									  if (inputGridIndex_jj == inputGridIndex_ii)
										{
										  properties_indx_kp = kp + jj*massager->nPropertyNames;
										  C_D[i_y][j_y] = massager->property_covariance[inputGridIndex_jj][properties_indx_k][properties_indx_kp];
										  // Symmetry and efficiency
										  C_D[j_y][i_y] = C_D[i_y][j_y];
										  if (i_y==j_y && C_D[j_y][i_y]==0) C_D[j_y][i_y] = VARIANCE_EPS;  // FIXME.. SENSIBLE????
										  // here, can detect if C_D is non-diagonal....
										}
									  // only need to fill lower triangle...(col<=row)
									  // but this is messy to achieve without screwing up j_y ordering
									  j_y++;
									}
								}
							}
						  i_y++;
						}
					}
				}


			  // Now the design matrix A
			  A = dmatrix(0,n_y-1,0,n-1);
			  // zero A
			  memset(&A[0][0],0,n_y * n * sizeof(double));

			  // fill rows according to
			  // (1) Wells
			  int j_model;
			  i_y=0;
			  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
				{ 
				  // printf("Awell: layer %d\n",ii);
				  // loop over wells that contribute to local A
				  for (n_j=0; n_j < treeOrdererWells[ii].size; n_j++)
					{
					  int WellObsIndex = treeOrdererWells[ii].values[n_j];
					  // printf("A: observation %d\n",WellObsIndex);
					  // loop over list of output grid nodes linked to this well
					  // observation via kriging
					  weight_sum=0.0;
					  for (k=0; k< sI->WellkrigInfo[ii][WellObsIndex].size; k++)
						{
						  // weight of kriging contribution
						  wt = sI->WellkrigInfo[ii][WellObsIndex].values[k].val;
					  
						  // Use hash table to retrieve collapsed indexes corresponding to this
						  // absolute output grid index, and as if for property k=0
						  hashIndx = g_hash_table_lookup(collapsedOG_hashtable[ii],
														 (gpointer) &sI->WellkrigInfo[ii][WellObsIndex].values[k].indx);
						  if (hashIndx!=NULL)
							{
							  // this has to be bumped along by the property index of the property we are observing
							  j_model = *((int*)hashIndx) + massager->well_observations[ii].propertyIndex[WellObsIndex];
							  A[i_y][j_model] = wt;
							  weight_sum += wt;
							}
						}
					  // renormalize so weight sum is 1 (only some kriging neighbours of the well may have 
					  // found their way into the set of output grid nodes in the local prior-model neighbourhood)
					  if (weight_sum==0) weight_sum=1.0; // probly never happens but be safe
					  maxWellCov=0.0;
					  for (k=0; k< sI->WellkrigInfo[ii][WellObsIndex].size; k++)
						{
						  // weight of kriging contribution
						  wt = sI->WellkrigInfo[ii][WellObsIndex].values[k].val;
						  hashIndx = g_hash_table_lookup(collapsedOG_hashtable[ii],
														 (gpointer) &sI->WellkrigInfo[ii][WellObsIndex].values[k].indx);
						  int OGindx = sI->WellkrigInfo[ii][WellObsIndex].values[k].indx; // OG grid index
						  if (hashIndx!=NULL)
							{
							  j_model = *((int*)hashIndx) + massager->well_observations[ii].propertyIndex[WellObsIndex];
							  A[i_y][j_model] /= weight_sum;
							  /*
								printf("AWELLENTRY node %d  OG indx %d, collapsed indexes %d,%d, wt=%g\n",j,OGindx,i_y,j_model,A[i_y][j_model]);
							  */
							  maxWellCov = MAX(maxWellCov, corr2D(pt[0], pt[1],massager->outputgrid[ii].Points[OGindx][0],
																  massager->outputgrid[ii].Points[OGindx][1],&(massager->variograms[ii])));
							}
						}
					  /* And, in case we don't have all the kriging contributors to the well, soften 
						 the C_D constraint proportionally to the correlation distance to the well: 
						 as if sigma <- sigma/weight_sum.

						 If the kriging contributors are just on the edge of the system, we should have maxWellCov/minUnscaledCorr == 1
						 So, soften the equation's strength by the factor (maxWellCov/minUnscaledCorr-1)
					  */

					  
					  double wellCorr = corr2D(pt[0], pt[1],massager->well_observations[ii].xValues[WellObsIndex],
										massager->well_observations[ii].yValues[WellObsIndex],&(massager->variograms[ii]));
					  wt2 = (wellCorr/corr2D(0,0,0,0,&(massager->variograms[ii])));
					  // wt = MIN(1,(maxWellCov/minUnscaledCorr - 1.0));					  
					  wt = fabs(MAX((wellCorr - minUnscaledCorr),0)/((massager->variograms[ii].sill+massager->variograms[ii].nugget) - minUnscaledCorr));
					  
					  printf("WELL-EDGE %g %g %g %g %g\n",pt[0],pt[1],wellCorr,minUnscaledCorr,wt);



					  C_D[i_y][i_y] /= ((wt>0)? wt*wt: 1);

					  i_y++;
					}
				}
			  // (2) Soft observations...
			  for (ii=0; ii<  massager->nLayers ; ii++)  // layers
				{ 
				  // printf("Asoft: layer %d\n",ii);
				  // loop over input grid subset, then properties at each input grid subset point
				  for (n_j=0; n_j < treeOrdererSoft[ii].size; n_j++)
					{
					  int SoftObsIndex = treeOrdererSoft[ii].values[n_j];
					  //printf("A: observation %d\n",SoftObsIndex);
					  // properties loop at input-grid subset point
					  for (k=0; k<massager->nPropertyNames; k++)
						{
						  // Loop over list of output grid nodes linked to this soft observation
						  // via kriging.
						  // This is the run across the row of A
						  weight_sum=0.0;
						  for (kp=0; kp< sI->IGkrigInfoSoft[ii][SoftObsIndex].size; kp++)
							{
							  // weight of kriging contribution
							  wt = sI->IGkrigInfoSoft[ii][SoftObsIndex].values[kp].val;
							  // Output grid index
							  int OGindx = sI->IGkrigInfoSoft[ii][SoftObsIndex].values[kp].indx; // OG grid index
							  // Use hash table to retrieve collapsed indexes corresponding to this
							  // absolute output grid index, and as if for property kp=0
							  hashIndx = g_hash_table_lookup(collapsedOG_hashtable[ii],
															 (gpointer) &sI->IGkrigInfoSoft[ii][SoftObsIndex].values[kp].indx);
							  if (hashIndx!=NULL)
								{
								  // this has to be bumped along by the property index of the property we are observing
								  j_model = *((int*)hashIndx) + k;
								  // printf("Asoft: OG indx %d, collapsed index %d, wt=%g\n",OGindx,j_model,wt);
								  // A[i_y][j_model] = wt;  raw weights may not necc sum to 1 since not all the OK neighborhood may be involved
								  weight_sum += wt;   // acquire sum of weight contributions along row
								}
							}
						  // same loop, but with normalised weights wt <- wt/weight_sum. Row sum is then 1.
						  for (kp=0; kp< sI->IGkrigInfoSoft[ii][SoftObsIndex].size; kp++)
							{
							  wt = sI->IGkrigInfoSoft[ii][SoftObsIndex].values[kp].val;
							  hashIndx = g_hash_table_lookup(collapsedOG_hashtable[ii],
															 (gpointer) &sI->IGkrigInfoSoft[ii][SoftObsIndex].values[kp].indx);
							  if (hashIndx!=NULL)
								{
								  j_model = *((int*)hashIndx) + k;
								  A[i_y][j_model] = wt/weight_sum;
								  if (massager->ignore_soft) A[i_y][j_model] = 0;
								}
							}
						  i_y++;
						}
					}
				}
			}
		  

		  // Now solve the local Bayesian update
#if 0
		  {
			printf("Cp: rows=%d, cols=%d\n",n,n);
			print_dmatrix0(Cp, n, n, "Cp.loc.txt");
		  
			printf("mu_p: dim=%d\n",n);
			print_dvector0(mu_p, n,"mu.loc.txt");
		  
			if (n_y>0)
			  {
				printf("A: rows=%d, cols=%d\n",n_y,n);
				print_dmatrix0(A, n_y, n, "A.loc.txt");
		  
				printf("C_D: rows=%d, cols=%d\n",n_y,n_y);
				print_dmatrix0(C_D,n_y,n_y,"CD.loc.txt");
		   
				printf("y: dim=%d\n",n_y);
				print_dvector0(y, n_y,"y.loc.txt");
			  }
		  }
#endif   
		  // Now HERE we can solve the local Bayesian dense problem..
		  if (n_y>0) 
			{
			  localUpdate = updateBayes(Cp, mu_p, n, A, C_D, y, n_y, 0, 0);
			  if (localUpdate.info > 0)
				{
				  printf("FAILURE: local Bayes update failed with code %d\n",localUpdate.info);
				  return EXIT_FAILURE;
				}
			}
		  else
			{
			  localUpdate.mu_tilde = mu_p;
			  localUpdate.C_tilde = Cp;
			}
		  // Solution of interest is that at central node on main layer i
		  // hashIndx points to its collapsed index.
		  hashIndx = (int*)g_hash_table_lookup(collapsedOG_hashtable[i],(gpointer) &j);
		  if (hashIndx==NULL)
			{
			  printf("Found a hashing bug in solveSequentialSparseProblems(). FIXME\n");
			  return EXIT_FAILURE;
			}
		  else
			{
			  // possible negativity issues (should do QP : FIXME!!!)
			  for (k=0; k<massager->nPropertyNames; k++)
				{
				  // positivity
				  if (localUpdate.mu_tilde[k+ *hashIndx] < 0)
					{
					  warning("Resetting negative property!");
					  localUpdate.mu_tilde[k+ *hashIndx] = 0.0;
					}
				  if (massager->properties[k].normalised)
					{
					  if (localUpdate.mu_tilde[k+ *hashIndx] >1.0)
						{
						  warning("Thresholding normalised property!");
						  localUpdate.mu_tilde[k+ *hashIndx] = 1.0;
						}
					}
				}
			  // outputs
			  if (massager->writeAsciiOut)
				{
   				  // these are the solutions we want
				  // This is the set of properties at the central grid node...
				  // write out x,y
				  fprintf(massager->ascii_out_fp[i],"%g %g ",
						  massager->outputgrid[i].xValues[j],
						  massager->outputgrid[i].yValues[j]);
				  /* add utm stuff if need be */
				  if (massager->utmInfo.useUTMconversion)
					{
					  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
												  massager->outputgrid[i].xValues[j],
												  massager->outputgrid[i].yValues[j]);
					  fprintf(massager->ascii_out_fp[i],"%10.2f %10.2f ",pt_utm.x,pt_utm.y);
					} 
				  // write out means of all properties
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  /*
						printf("layer %d node %d mu_tilde=%g krigVar=%g\n",i,j,
						localUpdate.mu_tilde[k+ *hashIndx],
						sqrt(localUpdate.C_tilde[k+ *hashIndx][k+ *hashIndx]));
					  */
					  fprintf(massager->ascii_out_fp[i],"%12.3f ",localUpdate.mu_tilde[k+ *hashIndx]);
					}
				  // and std deviations
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  fprintf(massager->ascii_out_fp[i],"%12.3f ",sqrt(localUpdate.C_tilde[k+ *hashIndx][k+ *hashIndx]));
					}
				  fprintf(massager->ascii_out_fp[i],"\n"); // FIXME REDEEM!
				}
			  // simulations
			  if (massager->Nrealisations > 0)
				{
				  // Form submatrix of known points and current unknown points & properties at this output grid node
				  int i_L,j_L,k_L,nKnownNbrs,Lsize;
				  int *requiredRows = ivector(0,n-1);
				  double **CtildeSub;
				  double sum;

				  // save this very useful array
				  for (k=0; k<massager->nPropertyNames; k++) posteriorMeans[i][j][k] = localUpdate.mu_tilde[k+ *hashIndx];

				  // count known neighbours & find their index positions
				  nKnownNbrs=0;
				  for (ii=0,i_p=0; ii<massager->nLayers ; ii++)           // layers
					{
					  for (n_i=1; n_i <= nNeighborsDesired[ii]; n_i++) // neighbours
						{
						  for (k=0; k< massager->nPropertyNames; k++,i_p++) // properties
							{
							  if (sI->simFlag[ii][neigh[ii][n_i]][k] == 1)
								{
								  // yes we have simulated this point already
								  requiredRows[nKnownNbrs++]=i_p;
								}
							}
						}
					}
				  knownValues = dmatrix(0, nKnownNbrs-1, 0, massager->Nrealisations-1);
				  nKnownNbrs=0;
				  for (ii=0,i_p=0; ii<massager->nLayers ; ii++)           // layers
					{
					  for (n_i=1; n_i <= nNeighborsDesired[ii]; n_i++) // neighbours
						{
						  for (k=0; k< massager->nPropertyNames; k++,i_p++) // properties
							{
							  if (sI->simFlag[ii][neigh[ii][n_i]][k] == 1)
								{
								  // yes we have simulated this point already
								  // kriged mean at this neighbour point
								  mu_neigh = posteriorMeans[ii][neigh[ii][n_i]][k];
								  for (i_sample=0 ; i_sample < massager->Nrealisations; i_sample++)
									{
									  //printf("Sim node %d, neigh %d: val=%g, mean=%g\n",j,nKnownNbrs,simulationValues[ii][neigh[ii][n_i]][k][i_sample],mu_neigh);
									  knownValues[nKnownNbrs][i_sample] = simulationValues[ii][neigh[ii][n_i]][k][i_sample] - mu_neigh;
									}
								  nKnownNbrs++;
								}
							}
						}
					}

				  // printf("Known neighbor size = %d\n",nKnownNbrs);
				  Lsize = nKnownNbrs + massager->nPropertyNames;
				  CtildeSub = dmatrix(0, Lsize-1, 0, Lsize-1);
				  // populate submatrix
				  // known neighbours
				  for (i_L=0; i_L<nKnownNbrs;i_L++) 
					{
					  for (j_L=0; j_L<=i_L;j_L++) 
						{
						  CtildeSub[i_L][j_L] = localUpdate.C_tilde[requiredRows[i_L]][requiredRows[j_L]];
						  CtildeSub[j_L][i_L] = CtildeSub[i_L][j_L];
						}
					}
				  // submatrix of properties at this point correlated with neighbours
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  for (j_L=0; j_L<nKnownNbrs;j_L++) 
						{
						  CtildeSub[nKnownNbrs+k][j_L] = localUpdate.C_tilde[k+ *hashIndx][requiredRows[j_L]];
						  CtildeSub[j_L][nKnownNbrs+k] = CtildeSub[nKnownNbrs+k][j_L]; // symmetry
						}
					}
				  // final corner block
				  for (k=0; k<massager->nPropertyNames; k++)
					{
					  for (k_L=0; k_L<massager->nPropertyNames; k_L++)
						{
						  CtildeSub[nKnownNbrs+k][nKnownNbrs+k_L] = localUpdate.C_tilde[k+ *hashIndx][k_L+*hashIndx];
						  CtildeSub[nKnownNbrs+k_L][nKnownNbrs+k] = CtildeSub[nKnownNbrs+k][nKnownNbrs+k_L];
						}
					}
				  // Cholesky decompose CtildeSub
				  //print_dmatrix0_mma(CtildeSub, Lsize, Lsize, "Ctildesub.txt");
				  info = clapack_dpotrf(CblasColMajor,CblasUpper,Lsize,&CtildeSub[0][0],Lsize);
				  if (info>0) 
					{
					  warning("lapack Cholesky decomposition failed in random simulation\n");
					  return EXIT_FAILURE;
					}
				  // alloc row for RHS's in the algebra
				  double **qTrans = dmatrix(0, massager->nPropertyNames,0,nKnownNbrs-1);

				  // And solve for the requisite simulation pieces
				  /* int clapack_dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
					 const int N, const int NRHS, const double *A, const int lda,
					 double *B, const int ldb);
				  */
				  //printf("Ctildesubafter: rows=%d, cols=%d\n",Lsize,Lsize);
				  //print_dmatrix0_mma(CtildeSub, Lsize, Lsize, "Ctildesubafter.txt");

				  // backsolve ala NR
				  // FIXME GET clapack_dpotrs to do this !!!!!!!!!!!!!!
				  // double *x = dvector(0,nKnownNbrs-1);
				  for (k_L=0; k_L< massager->nPropertyNames; k_L++)
					{
					  for (i_L=nKnownNbrs-1 ; i_L>=0 ; i_L--)
						{
						  sum = CtildeSub[k_L+nKnownNbrs][i_L];
						  for (kk=i_L+1; kk<nKnownNbrs; kk++) sum -= CtildeSub[kk][i_L] * qTrans[k_L][kk];
						  qTrans[k_L][i_L] = sum/CtildeSub[i_L][i_L];
						}
					  // overwrite solution into qTrans
					  // for (i_L=0; i_L< nKnownNbrs; i_L++) qTrans[k_L][i_L] = x[i_L];
					}
				  // generate samples !!
				  double *L_dm;
				  L_dm = dvector(0,massager->nPropertyNames-1);
				  double *gaussRands = dvector(0,massager->nPropertyNames-1);
				  // FIXME do this as a dgemm!!!!!!!!!!!!
				  for (i_sample=0 ; i_sample < massager->Nrealisations; i_sample++)
					{
					  for (k_L=0; k_L< massager->nPropertyNames; k_L++)
						{
						  L_dm[k_L] = posteriorMeans[i][j][k_L];
						  for (i_L=0; i_L<nKnownNbrs; i_L++) L_dm[k_L] += qTrans[k_L][i_L] * knownValues[i_L][i_sample];
						  gaussRands[k_L] = gasdev(&seed);
						}
					  for (k_L=0; k_L< massager->nPropertyNames; k_L++)
						{
						  for (i_L=0; i_L<=k_L; i_L++)
							{
							  L_dm[k_L] += CtildeSub[nKnownNbrs+k_L][nKnownNbrs+i_L] * gaussRands[i_L];
							}
						  // here's the realisation
						  simulationValues[i][j][k_L][i_sample] = L_dm[k_L];
						  // It has to be thresholded positive
						  fprintf(massager->ascii_out_fp[i],"%10.3g ", MAX(0,simulationValues[i][j][k_L][i_sample]));
						}
					}
				  fprintf(massager->ascii_out_fp[i],"\n");







				  free_dvector(L_dm,0,massager->nPropertyNames-1);
				  free_dvector(gaussRands,0,massager->nPropertyNames-1);



				  //free_dvector(x,0,nKnownNbrs-1);
				  // mark grid point as simulated
				  for (k=0; k<massager->nPropertyNames; k++) sI->simFlag[i][j][k] = 1;

				  free_dmatrix(qTrans, 0, massager->nPropertyNames,0,nKnownNbrs-1);
				  free_dmatrix(CtildeSub,0, Lsize-1, 0, Lsize-1);
				  free_ivector(requiredRows,0,n-1);
				  free_dmatrix(knownValues, 0, nKnownNbrs-1, 0, massager->Nrealisations-1);
				}
			}

		  // free stuff
		  // free up local dense systems
		  free_dmatrix(Cp,0,n-1,0,n-1);
		  free_dvector(mu_p,0,n-1);
		  if (n_y>0)
			{
			  free_dmatrix(A,0,n_y-1,0,n-1);
			  free_dmatrix(C_D,0,n_y-1,0,n_y-1);
			  free_dvector(y,0,n_y-1); 
			}
		  for (ii=0; ii<massager->nLayers ; ii++)  // layers
			{
			  g_hash_table_destroy(collapsedOG_hashtable[ii]);
			  g_tree_destroy(OGindxTreeSoft[ii]);
			  g_tree_destroy(OGindxTreeWells[ii]);
			  free(neigh[ii]);
			  free_ivector(treeOrdererWells[ii].values, 0, treeOrdererWells[ii].size - 1);
			  free_ivector(treeOrdererSoft[ii].values, 0, treeOrdererSoft[ii].size - 1);
			}

		  free(treeOrdererWells);
		  free(treeOrdererSoft);

		}
	  // close output files as soon as possible
	  if (massager->writeAsciiOut)
		{
		  fclose(massager->ascii_out_fp[i]);
		  printf("Wrote: %s\n",massager->ascii_out_filenames[i]);
		  free(massager->ascii_out_filenames[i]);
		}
	}

  // big cleanup
  free_dmatrix(averages,0,massager->nLayers-1,0,massager->nPropertyNames-1);
  free_dmatrix(variances,0,massager->nLayers-1,0,massager->nPropertyNames-1);
  free(neigh);
  free_ivector(nNeighborsDesired,0, (massager->nLayers) -1);
  free_ivector(indxs,0,nIndxs);

  free(collapsedOG_hashtable);
  free(OGindxTreeWells);
  free(OGindxTreeSoft);

  if (massager->Nrealisations > 0)
	{
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  free_imatrix(sI->simFlag[i],0,massager->outputgrid[i].nPoints - 1, 0 ,massager->nPropertyNames-1);
		  free_f3tensor(simulationValues[i],0,massager->outputgrid[i].nPoints - 1,0,massager->nPropertyNames-1,0,massager->Nrealisations-1);
		  free_matrix(posteriorMeans[i],0,massager->outputgrid[i].nPoints - 1,0,massager->nPropertyNames-1);
		}
	  free(sI->simFlag);
	  free(simulationValues);
	  free(posteriorMeans);
	}

  if (massager->useMAtrend)
	{
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  free_matrix(sI->smoothMeans[i], 0, massager->outputgrid[i].nPoints-1, 0, massager->nPropertyNames-1);
		}
	  free(sI->smoothMeans);
	}

  return EXIT_SUCCESS;
}












gboolean showOrder(gpointer key, gpointer value, gpointer data)
{
  printf("key= %d\n",*((int*)key));

  // FALSE indicates to perform a full tree loop
  return FALSE;
}

gboolean fillIndexTripleArray(gpointer key, gpointer value, gpointer indxTriple)
{
  struct INDXTRIPLE *tO = (struct INDXTRIPLE*)indxTriple;
  int *indx =  &(tO->indx);
  int *values = (tO->values);

  // makes the values store the keys in order
  values[*indx] = *((int*)key);

  //printf("indx=%d value=%d\n",*indx, values[*indx]);
  (*indx)++;

  // FALSE indicates full tree loop
  return FALSE;
}



gboolean fillIndexArray(gpointer key, gpointer value, gpointer treeOrderer)
{
  struct TREE_ORDERER *tO = (struct TREE_ORDERER*)treeOrderer;
  GHashTable* hashTable= tO->hash_table;
  int *index_array = tO->index_array;
  int *indx =  &(tO->indx);

  // makes the values simply increment 0...n-1
  g_hash_table_insert(hashTable,key,(gpointer)&index_array[*indx]);
  (tO->reverse_indx)[*indx] = *((int*)key);

  // printf("key=%d indx=%d\n",*((int*)key),index_array[*indx]);
  (*indx)++;

  // FALSE indicates full tree loop
  return FALSE;
}

gint intCompare(gconstpointer a, gconstpointer b)
{
  return ( (*((int*)a) < *((int*)b)) ? -1 : ( (*((int*)a) > *((int*)b)) ? 1 : 0 ));
}

// reverse level based sort
int intPairCompare(const void *a, const void *b)
{
  struct INTPAIR *A,*B;
  A = (struct INTPAIR*)a;
  B = (struct INTPAIR*)b;
  return ( (*A).indx < (*B).indx ? 1 : ( (*A).indx > (*B).indx ? -1 : 0 ));
}


void writeAsciiHeaderInfo(int layerNumber,struct MASSAGER *massager,int isRealisationsFile)
{
  int i,k;

  if (isRealisationsFile==0)
	{
	  fprintf(massager->ascii_out_fp[layerNumber],
			  "Bayesian kriged summary file for layer %d, generated by massager run from XML file %s\n",
			  massager->layerNumbers[layerNumber],massager->massage_xml_file);
	  if (massager->utmInfo.useUTMconversion)
		{
		  fprintf(massager->ascii_out_fp[layerNumber],"%d\n%s\n%s\nx\ny\n",4+2*massager->nPropertyNames,
				  massager->inputgrid.xWord,massager->inputgrid.yWord);
		}
	  else
		{
		  fprintf(massager->ascii_out_fp[layerNumber],"%d\nx\ny\n",2+2*massager->nPropertyNames);
		}
	  for (k=0; k< massager->nPropertyNames; k++) 
		{
		  fprintf(massager->ascii_out_fp[layerNumber],"%s\n",massager->propertyNames[k]);
		  fprintf(massager->ascii_out_fp[layerNumber],"std_dev_%s\n",massager->propertyNames[k]);
		}
	}
  else
	{
	  fprintf(massager->ascii_out_realisations_fp[layerNumber],
			  "Bayesian kriged posterior realisations file for layer %d, generated by massager run from XML file %s\n",
			  massager->layerNumbers[layerNumber],massager->massage_xml_file);
	  if (massager->utmInfo.useUTMconversion)
		{
		  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\n%s\n%s\nx\ny\n",
				  4+massager->Nrealisations*massager->nPropertyNames,
				  massager->inputgrid.xWord,massager->inputgrid.yWord);
		}
	  else
		{
		  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\nx\ny\n",2+massager->Nrealisations*massager->nPropertyNames);
		}
	  for (i=0; i<massager->Nrealisations; i++)
		{
		  for (k=0; k< massager->nPropertyNames; k++) 
			{
			  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%s_sample_%d\n",massager->propertyNames[k],i);
			}
		}
	}
}

int computeMovingAverage(struct MASSAGER *massager, int smoothStdDeviations)
{
  int i,ii,j,k,i_y,kp,kpp,properties_indx_kp,properties_indx_kpp,OGindx,size;
  // Smooth
  int IGindx;
  double wt,val,sum,sumSq,pSum;
  double **Points;
  double pt[2];
  struct PT2D pt_utm;
  int nNeighborsDesired;
  int *neigh;
  struct SOLVER_INFO *sI = &massager->solverInfo;  // Useful shortcut

  printf("Computing moving trend inference\n");

  if (massager->utmInfo.useUTMconversion)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}

  // Use prebuilt kdtree of inputgrid from setUpMapping()

  // alloc room for means (and maybe std deviations)
  sI->smoothMeans = (float***)malloc((size_t) ((massager->nLayers) * sizeof(float**)));
  if (smoothStdDeviations) 
	{
	  sI->smoothPropertyCovariance = f3tensor(0,(massager->inputgrid.nPoints)-1,0,(massager->nProperties)-1,0,(massager->nProperties)-1);
	  memset(&(sI->smoothPropertyCovariance[0][0][0]), 0, 
			 massager->inputgrid.nPoints * massager->nProperties * massager->nProperties);
	}
  for (i=0; i< massager->nLayers ; i++) 
	{
	  sI->smoothMeans[i] = matrix(0,massager->outputgrid[i].nPoints-1, 0, massager->nPropertyNames-1);
	  // zero it
	  memset(&(sI->smoothMeans[i][0][0]), 0,massager->outputgrid[i].nPoints *  massager->nPropertyNames * sizeof(float));
	}
  nNeighborsDesired = MIN(massager->inputgrid.nPoints,
						  DEFAULT_MAX_NEIGHBOURS_IN_SPARSE_PRIOR);
  Points = massager->inputgrid.Points;
  for (i=0; i< massager->nLayers ; i++)
	{
	  // set up appropriate metric for the kriging neighbour search using a kdtree
	  setUp2DEuclidenGenMetric(massager->inputgrid.OptkdRoot->Kglobals,
							   M_PI*massager->variograms[i].azimuth/180.0,
							   massager->variograms[i].rangeX/massager->variograms[i].rangeY);
	  /* loop on all OG grid nodes in this layer */
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  // x,y of OG point in this layer
		  if (massager->utmInfo.useUTMconversion)
			{
			  pt_utm = getUTMSfromInlineXline(massager->utmInfo, massager->outputgrid[i].xValues[j], massager->outputgrid[i].yValues[j]);
			  pt[0] = pt_utm.x;
			  pt[1] = pt_utm.y;
			}
		  else
			{
			  pt[0] = massager->outputgrid[i].xValues[j];
			  pt[1] = massager->outputgrid[i].yValues[j];
			}
		  // find set of neighbours on input grid
		  // get neighbour list from kdtree
		  neigh = kdOptNNQuery(Points,
								  2,&pt[0],nNeighborsDesired,
								  EUCLIDEAN_GEN,massager->inputgrid.OptkdRoot,1);



		  // Form weighted average on this neighborhood, weights from variogram
		  // for (kp=0; kp<massager->nPropertyNames; kp++) massager->solverInfo.smoothMeans[i][j][kp]=0;
		  for (i_y=1,sum=0; i_y <= nNeighborsDesired; i_y++)
			{
			  IGindx = neigh[i_y];
			  // form weight from correlation function
			  wt = corr2D(pt[0],pt[1],
						  Points[IGindx][0],Points[IGindx][1],
						  &(massager->variograms[i]));
			  sum += wt;
			  for (kp=0; kp<massager->nPropertyNames; kp++)
				{
				  properties_indx_kp = kp + i*massager->nPropertyNames;
				  massager->solverInfo.smoothMeans[i][j][kp] += wt * massager->properties[properties_indx_kp].values[IGindx];
				}
			}

		  // renormalise
		  for (kp=0; kp<massager->nPropertyNames; kp++) massager->solverInfo.smoothMeans[i][j][kp]/=sum;
		  if (massager->verbosity>4)
			{
			  printf("MA %g %g ",massager->outputgrid[i].xValues[j],massager->outputgrid[i].yValues[j]);
			  for (kp=0; kp<massager->nPropertyNames; kp++)
				{				  
				  properties_indx_kp = kp + i*massager->nPropertyNames;	
				  printf("%g %g ",massager->properties[properties_indx_kp].values[j],
						 massager->solverInfo.smoothMeans[i][j][kp]);
				}
			  printf("\n");
			}
		}
	  KillMetricSpace(massager->inputgrid.OptkdRoot->Kglobals);
	}



  // Useful pre-computational step: map smooth surfaces to well observation points
  sI = &massager->solverInfo;  // Useful shortcut

  int *obsCounters = ivector(0,massager->nPropertyNames-1); 
  for (i=0; i<  massager->nLayers ; i++)
	{
	  for (kp=0; kp< massager->nPropertyNames; kp++) obsCounters[kp]=0;
	  for (j=0; j<massager->well_observations[i].n; j++)
		{
		  kp = massager->well_observations[i].propertyIndex[j];
		  // interpolate smooth map of kp-th property onto well-observation point for this jth observation
		  sum=0;
		  for (k=0; k <sI->WellkrigInfo[i][j].size; k++)
			{
			  wt = sI->WellkrigInfo[i][j].values[k].val;       // kriging weight
			  OGindx = sI->WellkrigInfo[i][j].values[k].indx;  // output grid index
			  sum += wt * massager->solverInfo.smoothMeans[i][OGindx][kp];
			}
		  // sum is now the interpolated smooth mean of the property corresponding to this jth observation
		  // at the well location
		  massager->well_observations[i].obs[kp][obsCounters[kp]].smoothVal = sum;
		  obsCounters[kp]++;
		}
	}
  free_ivector(obsCounters,0,massager->nPropertyNames-1);

  // Smooth property variances on input grid also. 
  // FIXME: TEMPORARY measure:
  // Do this based on first layer's variogram
  // set up appropriate metric for the kriging neighbour search using a kdtree
  setUp2DEuclidenGenMetric(massager->inputgrid.OptkdRoot->Kglobals,
						   M_PI*massager->variograms[0].azimuth/180.0,
						   massager->variograms[0].rangeX/massager->variograms[0].rangeY);

  size = massager->nPropertyNames * massager->nLayers;

  /* loop on all IG grid nodes in this layer */
  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  // x,y of IG point
	  if (massager->utmInfo.useUTMconversion)
		{
		  pt_utm = getUTMSfromInlineXline(massager->utmInfo, massager->inputgrid.xValues[j], massager->inputgrid.yValues[j]);
		  pt[0] = pt_utm.x;
		  pt[1] = pt_utm.y;
		}
	  else
		{
		  pt[0] = massager->inputgrid.xValues[j];
		  pt[1] = massager->inputgrid.yValues[j];
		}
	  // find set of neighbours on input grid
	  // get neighbour list from kdtree
	  neigh = kdOptNNQuery(Points,
						   2,&pt[0],nNeighborsDesired,
						   EUCLIDEAN_GEN,massager->inputgrid.OptkdRoot,1);


	  // Form weighted average on this neighborhood, weights from variogram

	  for (i_y=1,sum=0; i_y <= nNeighborsDesired; i_y++)
		{
		  IGindx = neigh[i_y];
		  // form weight from correlation function
		  wt = corr2D(pt[0],pt[1],
					  Points[IGindx][0],Points[IGindx][1],
					  &(massager->variograms[0]));
		  sum += wt;		  
		  for (i=0; i<size; i++)
			{
			  for (ii=i; ii< size; ii++)
				{
				  massager->solverInfo.smoothPropertyCovariance[j][i][ii] += wt * 
					massager->property_covariance[IGindx][i][ii];
				}
			}
		}
	  if (massager->verbosity>4)
		{
		  printf("MA covariance %g %g ",massager->inputgrid.xValues[j],massager->inputgrid.yValues[j]);
		}
	  // renormalise at this trace
	  for (i=0; i<size; i++)
		{
		  for (ii=i; ii< size; ii++)
			{
			  massager->solverInfo.smoothPropertyCovariance[j][i][ii] /=sum;
			  massager->solverInfo.smoothPropertyCovariance[j][ii][i] = massager->solverInfo.smoothPropertyCovariance[j][i][ii];
			  if (massager->verbosity>4)
				{
				  printf("%g ",massager->solverInfo.smoothPropertyCovariance[j][i][ii]);
				}
			}
		}
	  if (massager->verbosity>4) printf("\n");

	  // Could precompute and store Cholesky factor here for realisations.....
	}

  KillMetricSpace(massager->inputgrid.OptkdRoot->Kglobals);


  free_dmatrix(massager->inputgrid.Points,0,massager->inputgrid.nPoints-1,0,1);
  KillOptTree(massager->inputgrid.OptkdRoot);
  printf("Finished smoothing of means and covariances\n");

  return EXIT_SUCCESS;
}
