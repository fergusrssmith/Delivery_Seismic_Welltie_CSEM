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

#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <search.h>
// #include <malloc.h>

#include <glib.h>

#include "ECLIO/eclio.h"  /* structure definition and a little miscellany */
#include "ECLIO/access.h" /* access high-level routines */
#include "ECLIO/eclutil.h" /* for model initialization only...my move to access */
#include "ECLIO/blockops.h" /* also added some lower-level routines, flagged below */
#include "ECLIO/copy.h"    /* ways to set records for moving data between models */

#include "NR/nrutil.h"


#include "kdtree/nn.h"
#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"

#include "cgeom/vdefs.h"

#include "errorHandler.h"





#include "massager_types.h" 
#include "variogram.h"
#include "modified_left_Cholesky.h" 
#include "util.h"
#include "solver.h"

  
#include "upperTriangleMatrices.h"
#include "clapack_cblas_wrappers.h"
#include "massager.h"
#include "decorator.h"
#include "eclipseHandler.h"
#include "su_util.h"

#define DEFAULT_MAX_NEIGHBOURS_IN_MOVING_AVERAGE 50       // Number of neighbours used in moving averages (MA), *per layer, per property*
#define MIN_SMOOTHING_WEIGHT 1e-10                         // underflow-prevention eps for smoothing
#define MA_WEIGHT_SUM_EPS 1e-2                             // sum of MA weights must exceed this for weighted sums to operate
#define DEFAULT_MAX_NEIGHBOURS_IN_SIMULATIONS    25        // Number of neighbours used in sequential simulation
#define DEFAULT_N_NEIGHBOURS_IN_INTERPOLATION_KRIGING 8    // Number of neighbours used in interpolation
#define KRIG_WEIGHT_THRESHOLD 0.001                        // interpolation will drop neighbours with weaker contributions than this
#define POSITIVE_DEFINITE_MATRIX_EPS 1e-5                  // Should ensure matrices are pos. def.
#define USE_MODIFIED_CHOLESKY 1
#define USE_MODIFIED_CHOLESKY_IN_INTERPOLATION 1


/*
  Set up mappings between structured input grids and unstructured output grid using kriging
*/ 

int setUpMapping(struct MASSAGER *massager)
{
  int i,j,jj,k,kk,dimension;
  int info;   
  double** Points;
  double** OGPoints;
  float oldNugget;
  double** Ckrig;
  double* RHSkrig;
  int *ipiv_krig;
  double x,y;
  optkdNode *OptkdRoot;
  int* neigh;
  int dim=2;
  int   nNeighborsDesired;
  int   outputgrid_indx_offset;
  int  output_grid_indx;
  int count;
  int OKfailcount=0;

  double pt[2];
  struct PT2D pt_utm;

  struct SOLVER_INFO *sI;


  sI = &massager->solverInfo;  // Useful shortcut

  // Allocate
  // Store array of kriging weights for each grid, and coupled information
  sI->WellkrigInfo = (struct FLPAIR_ARRAY**)malloc((size_t) ((massager->nLayers) * sizeof(struct FLPAIR_ARRAY*)));

  sI->OGkrigInfoSoft = (struct FLPAIR_ARRAY**)malloc((size_t) ((massager->nLayers) * sizeof(struct FLPAIR_ARRAY*)));


  for (i=0; i< massager->nLayers ; i++)
	{
	  sI->WellkrigInfo[i] = (struct FLPAIR_ARRAY*)malloc((size_t) 
													   ((massager->well_observations[i].n) * sizeof(struct FLPAIR_ARRAY)));
	  sI->OGkrigInfoSoft[i] = (struct FLPAIR_ARRAY*)malloc((size_t) 
														 ((massager->outputgrid[i].nPoints) * sizeof(struct FLPAIR_ARRAY)));
	}
     
  // kriging stuff etc needs to be performed in UTM coordinates if at all possible
  // N.B. massager->outputgrid[i] will be in seismic coords if massager->duplicate applies
  if (massager->utmInfo.performInputGridUTMmapping)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	} 


  // loop over layers
  for (i=0; i<massager->nLayers; i++)
	{
	  printf("Setting up output grid mapping (kriging), layer %d\n",massager->layerNumbers[i]);
	  // (1) form kdtree of outputgrid, in UTM's if at all possible
	  massager->outputgrid[i].Points = dmatrix(0,massager->outputgrid[i].nPoints-1,0,1);
	  Points = massager->outputgrid[i].Points; // useful shortcut
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  // outputgrid needs to be transformed if duplicated and 
		  // transformation been done on the input grid
		  if (massager->utmInfo.performOutputGridUTMmapping)
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
	  if (!massager->outputgrid[i].OptkdRoot) return EXIT_FAILURE;

	  OptkdRoot = massager->outputgrid[i].OptkdRoot;
	  // set up appropriate metric for the kriging neighbour search using a kdtree, 
	  // suitable for ith output grid
	  setUp2DEuclidenGenMetric(OptkdRoot->Kglobals,
							   M_PI*massager->variograms[i].azimuth/180.0,
							   massager->variograms[i].rangeX/massager->variograms[i].rangeY);
	  // Be careful about number of neighbours
	  nNeighborsDesired = DEFAULT_N_NEIGHBOURS_IN_INTERPOLATION_KRIGING;
	  if (nNeighborsDesired > massager->outputgrid[i].nPoints ) 
		{
		  nNeighborsDesired = massager->outputgrid[i].nPoints-1;   
		}
	  oldNugget = massager->variograms[i].nugget;
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
		  pt[0] = g_array_index( massager->well_observations[i].xValues, gdouble, j);
		  pt[1] = g_array_index( massager->well_observations[i].yValues, gdouble, j);
		  // get list of outputgrid neighbours of this point from kdtree
		  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,OptkdRoot,1);
		  OKfailcount=0;
OK:
		  // form OK kriging system
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  // Matrix. Need form only lower triangle
			  for (kk=1; kk <= k ; kk++) 
				{
				  // last index is "dimension" => 0 = x , 1 = y
				  // use safe exp variogram for interpolation. Metric same as specified variogram
				  Ckrig[k-1][kk-1] = corr2Dexp(Points[neigh[k]][0],  Points[neigh[k]][1],
											   Points[neigh[kk]][0], Points[neigh[kk]][1],
											   &(massager->variograms[i]));
				  Ckrig[kk-1][k-1] = Ckrig[k-1][kk-1]; // FIXME needed?
				}
			  // RHS
			  RHSkrig[k-1] = corr2Dexp(pt[0],pt[1],Points[neigh[k]][0],Points[neigh[k]][1],&(massager->variograms[i]));				             
			}
		  // edges
		  for (kk=0; kk < nNeighborsDesired ; kk++) Ckrig[nNeighborsDesired][kk] = Ckrig[kk][nNeighborsDesired] = 1.0;
		  Ckrig[nNeighborsDesired][nNeighborsDesired] = 0.0;
		  RHSkrig[nNeighborsDesired] = 1.0;

		  // Now solve kriging system for the weights mapping the output grid to the well-observation point
		  // print_dmatrix0(Ckrig,nNeighborsDesired+1,nNeighborsDesired+1,"junk.txt");
		  // print_dvector0(RHSkrig,nNeighborsDesired+1,"junk2.txt");
		  if (OKfailcount==0)
			{
			  info = clapack_wrapper_dgetrf(nNeighborsDesired+1,nNeighborsDesired+1,
									&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0]);
			  if (info>0)
				{
				  OKfailcount++; // use robust Cholesky next time....
				  warning("LU decomposition failed in kriging step: using safe Cholesky\n");
				  goto OK;
				  // restore
				  // if (OKfailcount!=0) massager->variograms[i].nugget = oldNugget;
				}
			  /* Get weights via LU solution solution: Must solve OK systems with LU methods (matrix not PD)  */
			  clapack_wrapper_dgetrs('N',nNeighborsDesired+1,1,&Ckrig[0][0],
							 nNeighborsDesired+1,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired+1);

			}
		  else
			{
			  solveOKsystem(Ckrig,RHSkrig,nNeighborsDesired);
			}  

		  /*
			Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1]. 
			
		  */
		  
		  // Wells: store kriging weights etc.
		  // eliminate negligible contributions
		  int krigCount=0;
		  for (k=0 ; k <nNeighborsDesired; k++)
			{
			  if (fabs(RHSkrig[k]) > KRIG_WEIGHT_THRESHOLD) krigCount++;
			}
		  // allocate room and store
		  sI->WellkrigInfo[i][j].size = krigCount;
		  sI->WellkrigInfo[i][j].values = 
			(struct FLPAIR*)malloc((size_t) ((krigCount) * sizeof(struct FLPAIR)));
		  if (!sI->WellkrigInfo[i][j].values) 
			{
			  fprintf(stderr,"alloc error in setUpMapping()");
			  return EXIT_FAILURE;
			}
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
				  kk++;
				}
			}
			

		  free(neigh);  // neighbourhood of jth observation
		  neigh=NULL;
		}


	  // Clean up kdtree metric
	  KillMetricSpace(OptkdRoot->Kglobals);
	  // Leave tree in place for subsequent searches though

	  free_dmatrix(Ckrig,0,nNeighborsDesired,0,nNeighborsDesired); Ckrig=NULL;
	  free_dvector(RHSkrig,0,nNeighborsDesired); RHSkrig=NULL;
	  free_ivector(ipiv_krig,0,nNeighborsDesired); ipiv_krig=NULL;
	}


  // ---------------------------------------------------------------------------------
  //  Reverse interpolation structure from input grid to output grid
  // ---------------------------------------------------------------------------------
  printf("Setting up input grid mapping (kriging)\n");
  // Form kdtree of inputgrid. Map to utm's if desired
  massager->inputgrid.Points = dmatrix(0,massager->inputgrid.nPoints-1,0,1);
  Points = massager->inputgrid.Points; // useful shortcut
  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  /* add utm stuff if need be */
	  if (massager->utmInfo.performInputGridUTMmapping)
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
  nNeighborsDesired = DEFAULT_N_NEIGHBOURS_IN_INTERPOLATION_KRIGING;
  if (nNeighborsDesired > massager->inputgrid.nPoints) 
	{
	  nNeighborsDesired = massager->inputgrid.nPoints;
	}
  // workspace for the OK kriging 
  Ckrig = dmatrix(0,nNeighborsDesired,0,nNeighborsDesired);
  RHSkrig = dvector(0,nNeighborsDesired);
  ipiv_krig = ivector(0,nNeighborsDesired);

  // Count usage of input grid points. Major time/mem savings can be made in some cases.
  // The idea here is that probably all input grid points are used to form the smoothed means/covars,
  // but only a  subset of them for subsequent interpolation onto outputgrid. Therefore
  // we form smoothed averages  only on the needed subset.

  massager->inputgrid.isUsedByInterpolation = cvector(0,massager->inputgrid.nPoints-1);
  memset(&(massager->inputgrid.isUsedByInterpolation[0]),0,massager->inputgrid.nPoints);
  unsigned char *isUsedByInterpolation = massager->inputgrid.isUsedByInterpolation;
  int nVisitCnt=0;

  // set up kriging structures for interpolation
  // loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  OGPoints = massager->outputgrid[i].Points; // useful shortcut
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
		  pt[0] = OGPoints[j][0];
		  pt[1] = OGPoints[j][1];
		  // get list of outputgrid neighbours of this point from kdtree
		  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,OptkdRoot,1);
		  if (nNeighborsDesired>0)
			{
			  // form OK kriging system
			  for (k=1; k <= nNeighborsDesired; k++) 
				{
				  // Matrix. Need form only lower triangle
				  for (kk=1; kk <= k ; kk++)
					{
					  // last index is "dimension" => 0 = x , 1 = y
					  // use safe exp variogram for interpolation. Metric same as specified variogram
					  Ckrig[k-1][kk-1] = corr2Dexp(Points[neigh[k]][0],  Points[neigh[k]][1],
												   Points[neigh[kk]][0], Points[neigh[kk]][1],
												   &(massager->variograms[i]));
					  Ckrig[kk-1][k-1] = Ckrig[k-1][kk-1]; // FIXME needed?
					}
				  // RHS
				  RHSkrig[k-1] = corr2Dexp(pt[0],pt[1],Points[neigh[k]][0],Points[neigh[k]][1],&(massager->variograms[i]));				             
				}
			  // edges
			  for (kk=0; kk < nNeighborsDesired ; kk++) Ckrig[nNeighborsDesired][kk] = Ckrig[kk][nNeighborsDesired] = 1.0;
			  Ckrig[nNeighborsDesired][nNeighborsDesired]=0.0;
			  RHSkrig[nNeighborsDesired] = 1.0;

#if USE_MODIFIED_CHOLESKY_IN_INTERPOLATION
			  // Modified-cholesky style OK solver....
			  solveOKsystem(Ckrig,RHSkrig,nNeighborsDesired);
#else
			  // Now solve ordinary kriging system for the weights. OK is not pos. def, so use general solver. 
			  info = clapack_wrapper_dgetrf(nNeighborsDesired+1,nNeighborsDesired+1,
									&Ckrig[0][0],nNeighborsDesired+1,&ipiv_krig[0]);
			  if (info>0) 
				{
				  warning ("lapack decomposition failed in interpolation kriging step\n");
				  goto OKinterp;
				  // return EXIT_FAILURE;
				}
			  clapack_wrapper_dgetrs('N',nNeighborsDesired+1,1,&Ckrig[0][0],
							 nNeighborsDesired+1,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired+1);
#endif
			}

		  /*
			Solution (kriging weights) are now in RHSkrig[0...nNeighborsDesired-1].
		  */
		  // Store kriging weights etc.
		  // allocate room for storage of kriging weights and store
		  sI->OGkrigInfoSoft[i][j].size = nNeighborsDesired;
		  sI->OGkrigInfoSoft[i][j].values =	(struct FLPAIR*)malloc((size_t) (nNeighborsDesired * sizeof(struct FLPAIR)));

		  for (k=0; k <nNeighborsDesired; k++)
			{
			  sI->OGkrigInfoSoft[i][j].values[k].val =  RHSkrig[k];  // kriging weight
			  sI->OGkrigInfoSoft[i][j].values[k].indx = neigh[k+1];  // input grid index
			  if (!isUsedByInterpolation[neigh[k+1]])
				{
				  isUsedByInterpolation[neigh[k+1]] = 1;
				  nVisitCnt++;
				}
			}

		  free(neigh);
		  // just for cleanliness
		  OptkdRoot->Kglobals->nndist = NULL;
		}
	  KillMetricSpace(OptkdRoot->Kglobals);
	}

  printf("Number of used points on inputgrid = %d\n",nVisitCnt);

  
  free_dmatrix(Ckrig,0,nNeighborsDesired,0,nNeighborsDesired); Ckrig=NULL;
  free_dvector(RHSkrig,0,nNeighborsDesired); RHSkrig=NULL;
  free_ivector(ipiv_krig,0,nNeighborsDesired); ipiv_krig=NULL;

  // Seismic trace area calculations for pure seismic stats
  if (massager->acquireStats && massager->duplicate) 
	{
	  double dx,dy,dx1,dy1,dotP,L1,L2,area,minArea =DBL_MAX ;
	  double RIGHT_ANGLE_EPS=0.01;
	  int nAreas=0;
	  double *areas=NULL;

	  if (massager->inputgrid.nPoints<5)
		{
		  printf("Too few traces to get seismic trace area, assuming 1.0\n");
		  massager->inputgrid.traceArea = 1.0;
		}
	  else
		{
		  nNeighborsDesired = 5;

		  Points = massager->inputgrid.Points; // useful shortcut
		  OptkdRoot = massager->inputgrid.OptkdRoot;

		  // Run over all inputgrid points. Find the min area associated with a seismic grid
		  areas = dvector(0,massager->inputgrid.nPoints-1);
		  printf("Determining seismic trace area...");
		  for (j=0; j < massager->inputgrid.nPoints; j++)
			{
			  // This pt will be in UTMs if at all possible
			  pt[0] = Points[j][0];
			  pt[1] = Points[j][1];
			  // printf("TRIGRIG %f %f\n",pt[0],pt[1]);
			  // get list of nearest neighbours of this point from kdtree, which are neigh[1...5]
			  // one of these is itself
			  neigh = kdOptNNQuery(Points,dim,&pt[0],nNeighborsDesired,
								   EUCLIDEAN,OptkdRoot,1);
			  // swap so pt is always indexed neigh[1]
			  for (k=2; k<=5; k++) 
				{
				  if (neigh[k]==j) 
					{
					  neigh[k]=neigh[1]; neigh[1]=j; break;
					}
				}
			  // vector to one of the neighbours, the first indexed one
			  dx = Points[neigh[2]][0] - pt[0];
			  dy = Points[neigh[2]][1] - pt[1];
			  L1 = dx*dx+dy*dy;
			  // find another neighbour which is at right angles (i.e. dot product seems small)
			  for (k=3; k<=5; k++)
				{
				  dx1 = Points[neigh[k]][0] - pt[0];
				  dy1 = Points[neigh[k]][1] - pt[1];
				  dotP = fabs(dx1*dx+dy1*dy);
				  // printf("DOTP: %f\n",dotP/L1);
				  if (dotP < RIGHT_ANGLE_EPS * L1)
					{
					  // area = sqrt(dx1*dx1+dy1*dy1) * sqrt(L1);
					  area = sqrt((dx1*dx1+dy1*dy1) * L1);
#if 0
					  printf("TRI %f %f\n",pt[0],pt[1]);
					  printf("TRI %f %f\n",Points[neigh[2]][0],Points[neigh[2]][1]);
					  printf("TRI %f %f\n",Points[neigh[k]][0],Points[neigh[k]][1]);
					  printf("TRI %f %f\n",pt[0],pt[1]);
					  printf("TRI \n");
#endif
					  // printf("AREA %g\n",area);
					  minArea = MIN(minArea,area);
					  areas[nAreas] = area;
					  nAreas++;
					  break;
					}
				}

			  free(neigh);
			}

		  if (minArea == DBL_MAX)
			{
			  printf("Seismic x,y trace positions do not form an orthogonal grid. "
					 "Check header words, UTM triple in XML file, inline/xline header words etc\n");
			  exit(0);
			}
		  // second scan over areas to diagnose jitter - test for up to 50% area-jitter
		  area = 1.5*minArea;
		  for (i=0,j=0; i<nAreas; i++)
			{
			  if (areas[i] < area)
				{ // collect rectangles within 50% of minimum area at start of array
				  areas[j++] = areas[i];
				}
			}
		  nAreas = j;
		  // get median of these nAreas rectangles
		  area = quick_select(areas,nAreas);
		  printf("minimum seismic trace area = %.2f\n",minArea);
		  printf("Using median seismic trace area for bins within 50%% of minimum: %.2f\n",area);
		  massager->inputgrid.traceArea = area;
		  free_dvector(areas,0,massager->inputgrid.nPoints-1);
		}
	}

  return EXIT_SUCCESS;
}

/* Clean up */
int cleanUpMappingSetup(struct MASSAGER *massager)
{
  int i,j;
  struct SOLVER_INFO *sI;
  sI = &massager->solverInfo;  // Useful shortcut 
  
  for (i=0; i< massager->nLayers ; i++)
	{
	  // Clean up kdtrees of output grids
	  if (massager->outputgrid[i].Points)
		{
		  free_dmatrix(massager->outputgrid[i].Points,0,massager->outputgrid[i].nPoints-1,0,1);
		  massager->outputgrid[i].Points = NULL;
		}
	  KillOptTree(massager->outputgrid[i].OptkdRoot);

	  for (j=0; j<massager->well_observations[i].n; j++)
		{
		  free(sI->WellkrigInfo[i][j].values);
		  sI->WellkrigInfo[i][j].values = NULL;
		}
	  free(sI->WellkrigInfo[i]);
	  sI->WellkrigInfo[i] = NULL;

	  freeObservations(massager,i);
	  for (j=0; j < massager->outputgrid[i].nPoints; j++)
		{
		  free(sI->OGkrigInfoSoft[i][j].values);
		  sI->OGkrigInfoSoft[i][j].values = NULL;
		}
	  free(sI->OGkrigInfoSoft[i]);
	  sI->OGkrigInfoSoft[i] = NULL;
	}
  
  if (massager->inputgrid.Points)
	{
	  free_dmatrix(massager->inputgrid.Points,0,massager->inputgrid.nPoints-1,0,1);
	  massager->inputgrid.Points=NULL;
	}

  free(sI->WellkrigInfo);
  sI->WellkrigInfo= NULL;
  free(sI->OGkrigInfoSoft);
  sI->OGkrigInfoSoft = NULL;

  return EXIT_SUCCESS;  
}



/*
  Approximate solution of problem using scaling from smoothed covariances
*/
int solveScaledSmoothSparseProblems(struct MASSAGER *massager)
{
  int i,j,k,kp,m,mm,sizeP,output_grid_indx;
  int properties_indx_k,properties_indx_kp;
  double val,thickness_val,net_sand_val,sum,sumSq,minUnscaledCorr,unscaledCorr,maxWellCov;
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
  struct SOLVER_INFO* sI = &massager->solverInfo;  // Useful shortcut


  // SK kriging stuff
  double** Ckrig;
  double* RHSkrig,*Crhs;
  float ***krigVariances,***krigMeans,***krigStdDev;
  int *ipiv_krig;
  int size=0,IGindex;
  double var,mean,variance,krigVar,krigCorrection;
  double interpolatedSmoothPropCovar;
  struct PT2D pt_utm;
  struct FLPAIR* krigPairs;
  double **OGPoints;

  // ECLIO stuff
  char makeStochasticEclipse = 0;
  char *name;
  struct ECL_OUT_INFO *eclOutInfo=NULL;

  // realisations stuff
  int info,coLocatedNbrIndx=0,i_sample,real_no;
  long seed = massager->seed;
  double condMean,condVar,deviate,deviate_scale;
  double z=0,area,*y_sum,*sim_values,*w_tmp,tmp[2];
  float ****simulationValues=NULL; // will be huge, so alloc as a float (not numerically sensitive)
  float ***posteriorMeans;         // mean for grid i,node j, property k 
  double **knownValues;            // working of array known neighbour values in local simulation machinery
  double **Lscale;
  double zStoch=0,qty;
  int dimStats=0;
  float **stats=NULL;              // stats for net-sand and thickness
  float **layer_sum_stats=NULL;    // stats for composite net-sand and thickness
  float mu_neigh;
  int **path;
  int nNeighborsDesired;
  dynamicTree *dTree;

  // decoration stuff
  int i_dec,decorate_this_layer;
  double dec_NG=0.5, dec_thick=0,dec_d=0,error;

  if (massager->verbosity>1)
	{
	  printf("Solving model sequentially\n");
	}

  if (massager->utmInfo.utmTripleAvailable)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	} 

  // Eclipse realisations flag
  makeStochasticEclipse = ((massager->Nrealisations > 0) && 
						   massager->eclInfo.writeECLOut && 
						   massager->outputgrid[0].eclOutInfo.outputFilename != NULL) ? 1 : 0;

  // ASCII IO: loop over layers
  for (i=0; i<  massager->nLayers ; i++)
	{
	  if (massager->writeAsciiOut) 
		{
		  writeAsciiHeaderInfo(i,massager,0);
		  if (massager->Nrealisations > 0) writeAsciiHeaderInfo(i,massager,1);
		}
	}

  // Presmooth means and std deviations
  computeMovingAverage(massager);

  krigVariances = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
  krigMeans = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
  krigStdDev = (float ***)malloc((size_t) ((massager->nLayers) * sizeof(float **)));
  // safe initialisation
  for (i=0; i< massager->nLayers ; i++) krigMeans[i] = krigVariances[i] = krigStdDev[i] = NULL;

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
	  // large storage needed for all realisations
	  simulationValues = (float ****)malloc((size_t) ((massager->nLayers) * sizeof(float ***)));
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  simulationValues[i] = f3tensor(0,massager->outputgrid[i].nPoints - 1,   // points
										 0,massager->nPropertyNames-1,            // props
										 0,massager->Nrealisations-1);            // realisation #
		}	  
	  if (massager->acquireStats)
		{
		  dimStats = massager->nPropertyNames;
		  if (dimStats>0)
			{
			  stats = matrix(0,massager->Nrealisations-1,0,massager->nPropertyNames-1);
			  layer_sum_stats = matrix(0,massager->Nrealisations-1,0,massager->nPropertyNames-1);
			}
		}
	}

  // (1) Well kriging loop
  if ((massager->verbosity>2)) printf("Starting well-kriging scaling loop...\n");
  for (i=0; i< massager->nLayers ; i++)
	{
	  if (massager->verbosity>2)
		{
		  printf("Kriging layer %d...\n",massager->layerNumbers[i]);
		}

	  /* loop on all OG grid nodes in this layer */
	  // Do pseudo multigrid by traversing binary tree indexing of points.
	  size = massager->outputgrid[i].nPoints;

	  // (1) Kriging Loop 
	  // Now traverse outputgrid in this order, solving SK kriging system
	  krigMeans[i] = matrix(0,size-1,0,massager->nPropertyNames-1);
	  krigVariances[i] = matrix(0,size-1,0,massager->nPropertyNames-1);
	  krigStdDev[i] = matrix(0,size-1,0,massager->nPropertyNames-1);

	  OGPoints = massager->outputgrid[i].Points; // useful shortcut

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
#if 0
					  if (massager->utmInfo.utmTripleAvailable)
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
#else
					  RHSkrig[k] = corr2D(OGPoints[output_grid_indx][0], OGPoints[output_grid_indx][1],
										  massager->well_observations[i].obs[m][k].x, massager->well_observations[i].obs[m][k].y,
										  &(massager->variograms[i]));	
#endif

					  RHSkrig[k] /= Ckrig[0][0]; // normalise to guarantee unit variance
					}
				  memcpy(Crhs,RHSkrig,nNeighborsDesired*sizeof(double));

#if !USE_MODIFIED_CHOLESKY
				  // Now solve kriging system for the weights mapping the output grid to the well-observation point
				  // print_dmatrix0(Ckrig,nNeighborsDesired+1,"junk.txt");
				  info = clapack_wrapper_dgetrf(nNeighborsDesired,nNeighborsDesired,
										&Ckrig[0][0],nNeighborsDesired,&ipiv_krig[0]);
				  if (info>0)
					{
					  warning ("lapack decomposition failed in well kriging step\n");
					  return EXIT_FAILURE;
					}   
				  
				  /* Get weights via dge solution solution */
				  clapack_wrapper_dgetrs('N',nNeighborsDesired,1,&Ckrig[0][0],
								 nNeighborsDesired,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired);
					  
#else
				  /* Modified Cholesky solve of SK system */
				  int *perm = ivector(0,nNeighborsDesired-1);
				  modified_Cholesky_left_decomposition(nNeighborsDesired,Ckrig,perm,0);
				  modified_Cholesky_solve(nNeighborsDesired,Ckrig,perm,RHSkrig);
				  free_ivector(perm,0,nNeighborsDesired-1); perm=NULL;
#endif


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
				  

				  free_dmatrix(Ckrig,0,nNeighborsDesired-1,0,nNeighborsDesired-1); Ckrig=NULL;
				  free_dvector(RHSkrig,0,nNeighborsDesired-1); RHSkrig=NULL;
				  free_dvector(Crhs,0,nNeighborsDesired-1); Crhs=NULL;
				  free_ivector(ipiv_krig,0,nNeighborsDesired-1); ipiv_krig=NULL;
				}
			  else
				{
				  mean = massager->solverInfo.smoothMeans[i][output_grid_indx][m];
				  krigVar = massager->variograms[i].sill + massager->variograms[i].nugget;
				}

			  // truncate net_sand kriging if both thickness + net_sand available
			  if ((m == massager->net_sand_index) && (massager->thickness_index>=0))
				{
				  thickness_val = krigMeans[i][output_grid_indx][massager->thickness_index];
				  mean = MIN(mean,thickness_val);
				}

			  // override kriging if NG computable from thickness & net_sand:
			  // NG_kriged \def net_sand_kriged/thickness_kriged;
			  if ((m == massager->NG_index) && (massager->thickness_index>=0) && (massager->net_sand_index>=0))
				{
				  thickness_val = krigMeans[i][output_grid_indx][massager->thickness_index];
				  net_sand_val = krigMeans[i][output_grid_indx][massager->net_sand_index];
				  if (thickness_val!=0) mean = net_sand_val/thickness_val;
				}
			  
			  // override kriging if vp_eff computable from times & thickness :
			  // FIXME needs to couple to next layer (i+1)...
			  if (0)
				{
				  if ((m == massager->vp_eff_index) && (massager->thickness_index >= 0) && (massager->t_index >= 0))
					{
					  double dTWT;
					  thickness_val = krigMeans[i][output_grid_indx][massager->thickness_index];
					  dTWT = krigMeans[i+1][output_grid_indx][massager->t_index]  // but i+1 not there yet FIXME (must be done in second pass)
						- krigMeans[i][output_grid_indx][massager->t_index]; // TWT in ms
					  if (dTWT > 0 ) mean = 2.0e3*thickness_val/dTWT;        // vp_eff in normal units
					}
				}


			  // Allow no negative kriged values:
			  mean = MAX(0.0, mean);

			  // also normalise to [0,1] if required.
			  if (massager->properties[m].normalised) mean = MIN(1.0, mean);

			  if (isnan(mean)) printf("NAN BUG!\n");

			  krigMeans[i][output_grid_indx][m] = mean;
			  krigVariances[i][output_grid_indx][m] = krigVar;
			  /* krigVariances[i][output_grid_indx][m] is now normalised kriging variance of property m, conditioned only to wells */
			  /* We will have krigVariances[i][output_grid_indx][m] -> 1 far from the wells */
			}
		}
	  if (massager->verbosity>2)
		{
		  printf("Finished kriging layer %d\n",massager->layerNumbers[i]);
		}

	  // Kriging I/O & variance calculation
	  if (massager->writeAsciiOut || massager->eclInfo.writeECLOut || massager->writeSUOut) 
		{
		  if (massager->verbosity>2 && massager->writeAsciiOut) 
			{
			  printf("Writing to: %s\n",massager->ascii_out_filenames[i]);
			}
		  // write out surfaces in natural order
		  for (j=0; j < size; j++) 
			{
			  if (massager->writeAsciiOut)
				{
				  // write out x,y
				  fprintf(massager->ascii_out_fp[i],"%12.4f %12.4f ",
						  massager->outputgrid[i].xValues[j],
						  massager->outputgrid[i].yValues[j]);
				  /* add utm stuff if need be */ 
				  if (massager->utmInfo.performOutputGridUTMmapping)
					{
					  fprintf(massager->ascii_out_fp[i],"%12.4f %12.4f ",OGPoints[j][0], OGPoints[j][1]);
					}
				  if (!massager->duplicate)
					{
					  // Eclipse
					  fprintf(massager->ascii_out_fp[i],"%6d %6d ",
							  massager->outputgrid[i].ecl_i[j], massager->outputgrid[i].ecl_j[j]);
					  if (massager->outputgrid[i].segnum)
						{
						  fprintf(massager->ascii_out_fp[i],"%4d ",massager->outputgrid[i].segnum[j]);
						}
					}
				}
			  // run over properties in this layer
			  for (m=0; m< massager->nPropertyNames; m++) 
				{
				  // write out +ve truncated kriged mean
				  if (massager->writeAsciiOut) fprintf(massager->ascii_out_fp[i],"%g ",MAX(krigMeans[i][j][m],0));

				  /* Compute variance: scale variance by smooth covariance interpolated from input grid */
				  properties_indx_k = m + i*massager->nPropertyNames;
				  // interpolate smoothpropertyCovariance from input grid to this output grid site, using precomputed kriging weights
				  interpolatedSmoothPropCovar = 0;
				  krigPairs = sI->OGkrigInfoSoft[i][j].values;
				  for (k=0; k< sI->OGkrigInfoSoft[i][j].size; k++)
					{
					  IGindex = krigPairs[k].indx;
					  wt = krigPairs[k].val;
					  interpolatedSmoothPropCovar += wt * getDiag(sI->smoothPropertyCovariance_tri[IGindex],
																  properties_indx_k);
					}
				  // scale this by the well-kriged variance to recover true variance
				  variance = krigVariances[i][j][m] * interpolatedSmoothPropCovar;
				  // store it
				  krigStdDev[i][j][m] = sqrt(MAX(variance,0));



				  // write out kriged std dev
				  if (massager->writeAsciiOut) fprintf(massager->ascii_out_fp[i],"%g ",krigStdDev[i][j][m]);
				}
			  if (massager->writeAsciiOut) fprintf(massager->ascii_out_fp[i],"\n");
			}
		  // close output files as soon as possible
		  if (massager->writeAsciiOut) 
			{
			  fclose(massager->ascii_out_fp[i]);
			  if (massager->verbosity>2) printf("Wrote: %s\n",massager->ascii_out_filenames[i]);
			  free(massager->ascii_out_filenames[i]);
			}
		}
	  
	  // ECLIO I/O
	  if (!massager->duplicate && massager->outputgrid[i].eclOutInfo.outputFilename != NULL)
		{
		  // kriging data is 'accumulated'
		  writeECLIO_InSituOutput(-1,1,i,massager,krigMeans,krigStdDev);
		}
	}

  // (1b) SU version of trends as above, rolled out by location slowest (i.e. tracewise)
  //      Only possible under duplicated grids at present
  if (massager->writeSUOut)  
	{
	  if (dumpPostMassagedSUFiles(massager,krigMeans,krigStdDev) == EXIT_FAILURE)
		{
		 return EXIT_FAILURE; 
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
		  dTree = buildDynamicTree(massager->outputgrid[i].Points, size, 2, path[i]);
		  // set up appropriate metric for the kriging neighbour search using the dynamic kdtree, 
		  // suitable for ith output grid
		  setUp2DEuclidenGenMetric(dTree->OptkdRoot->Kglobals,
								   M_PI*massager->variograms[i].azimuth/180.0,
								   massager->variograms[i].rangeX/massager->variograms[i].rangeY);

		  decorate_this_layer = (massager->nDecorations && (i==massager->decoratorInfo.decoration_layer_indx)) ? 1 : 0; 

		  if (decorate_this_layer) setUp3Dsimulation(massager,path[i]);

		  double Cscale = massager->variograms[i].nugget + massager->variograms[i].sill;

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
			  if (neigh==NULL && (j!=0))
				{
				  printf("BUG!\n"); return EXIT_FAILURE; 
				}
			  // neigh[1....nNeighborsDesired] are now the indexes 
			  // of the already simulated points on massager->outputgrid[i]
			  if (nNeighborsDesired>0)
				{
				  Ckrig = dmatrix(0,nNeighborsDesired-1,0,nNeighborsDesired-1);
				  RHSkrig = dvector(0,nNeighborsDesired-1);
				  Crhs = dvector(0,nNeighborsDesired-1);
				  ipiv_krig = ivector(0,nNeighborsDesired-1);
				  coLocatedNbrIndx=-1;

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
						  Ckrig[k][kk] /= Cscale; // normalise to guarantee unit variance
						  Ckrig[kk][k] = Ckrig[k][kk]; // FIXME needed?
						} 
					  // RHS.
					  RHSkrig[k] = corr2D(massager->outputgrid[i].Points[output_grid_indx][0],massager->outputgrid[i].Points[output_grid_indx][1],
										  massager->outputgrid[i].Points[neigh[k+1]][0], massager->outputgrid[i].Points[neigh[k+1]][1],  
										  &(massager->variograms[i]));
					  RHSkrig[k] /= Cscale; // normalise to guarantee unit variance
					  if (RHSkrig[k]>=1.0) coLocatedNbrIndx = k;
					}

				  // SOLVE SK system
				  if (coLocatedNbrIndx>=0)
					{
					  // one of the neighbours is co-located with the current point. The associated
					  // kriging weight must be one
					  memset(RHSkrig,0,nNeighborsDesired*sizeof(double));
					  // for (k=0; k < nNeighborsDesired; k++) RHSkrig[k]=0;
					  RHSkrig[coLocatedNbrIndx]=1.0;
					  condVar=0.0;
					}
				  else
					{
					  cleanKrigingSystemOfDuplicateNeighbours(nNeighborsDesired,Ckrig,RHSkrig);
					  memcpy(Crhs,RHSkrig,nNeighborsDesired*sizeof(double));
					  condVar = 1.0;  // store this before we mess up Ckrig
#if   !USE_MODIFIED_CHOLESKY
					  // Now solve kriging system for the weights
					  info = clapack_wrapper_dgetrf(nNeighborsDesired,nNeighborsDesired,
											&Ckrig[0][0],nNeighborsDesired,&ipiv_krig[0]);
					  if (info>0)
						{
						  warning ("lapack decomposition failed in kriging step\n");
						  return EXIT_FAILURE; 
						}   
				  
					  /* Get weights via dge solution solution */
					  clapack_wrapper_dgetrs('N',nNeighborsDesired,1,&Ckrig[0][0],
									 nNeighborsDesired,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired);
					  /*
						Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1].
						Generate unconditional realisations with normalised variogram
					  */
#else
					  /* Modified Cholesky solve of SK system */
					  int *perm = ivector(0,nNeighborsDesired-1);
					  modified_Cholesky_left_decomposition(nNeighborsDesired,Ckrig,perm,0);
					  modified_Cholesky_solve(nNeighborsDesired,Ckrig,perm,RHSkrig);
					  free_ivector(perm,0,nNeighborsDesired-1); perm=NULL;
					  /*
						Solution (kriging weights) are now in RHSkrig[0.....nNeighborsDesired-1].
						Generate unconditional realisations with normalised variogram
					  */
#endif

					  for (k=0; k <nNeighborsDesired; k++) condVar -= RHSkrig[k]*Crhs[k]; // standard SK kriging variance: GSLIB II.14
					  if (!finite(condVar))
						{
						  printf("MCBUG1\n");
						}
					}

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
							  if (!finite(simulationValues[i][neigh[k+1]][m][real_no]))
								{
								  printf("MCBUG2\n");
								}
							}
						  //condMean += mean;
						  deviate = deviate_scale * Gaussian_deviate(&seed);
						  // Generate realisations and store
						  simulationValues[i][output_grid_indx][m][real_no] = condMean + deviate;
						  if (!finite(simulationValues[i][output_grid_indx][m][real_no]))
							{
							  printf("BUG\n");
							}
						}
					}
				  free_dmatrix(Ckrig,0,nNeighborsDesired-1,0,nNeighborsDesired-1); Ckrig=NULL;
				  free_dvector(RHSkrig,0,nNeighborsDesired-1); RHSkrig=NULL;
				  free_dvector(Crhs,0,nNeighborsDesired-1); Crhs=NULL;
				  free_ivector(ipiv_krig,0,nNeighborsDesired-1); ipiv_krig=NULL;
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
						  deviate = deviate_scale * Gaussian_deviate(&seed);
						  // store
						  simulationValues[i][output_grid_indx][m][real_no] = mean + deviate;
						}
					}
				}
			  // Decoration....
			  if (decorate_this_layer) 
				{
				  decSetUp3DcolumnSim(massager,i,j);
				}
			  free(neigh);
			  neigh=NULL;
			}
		  // Cleanup dynamic tree
		  killDynamicTree(dTree);		  
		}

	  // Realisations I/O

	  if (massager->writeAsciiOut || makeStochasticEclipse || dimStats)
		{
		  y_sum = dvector(0,massager->nProperties-1);
		  sim_values = dvector(0,massager->nProperties-1);
		  w_tmp = dvector(0,massager->nProperties-1);

		  // Eclipse realisations


		  if (makeStochasticEclipse)
			{ 
			  eclOutInfo = &(massager->outputgrid[0].eclOutInfo);
			  allocateStochEclipseMem(massager);
			}

		  if (dimStats) memset(&layer_sum_stats[0][0],0,sizeof(float)*dimStats*massager->Nrealisations);

		  for (i=0; i< massager->nLayers ; i++)
			{
			  decorate_this_layer = (massager->nDecorations && (i==massager->decoratorInfo.decoration_layer_indx)) ? 1 : 0;
			  if (dimStats) 
				{
				  // reset individual layer stats cumulators
				  memset(&stats[0][0],0,sizeof(float)*dimStats*massager->Nrealisations);
				  // if simulating depth, use stochastic depth in stats. It's always property 0 by construction
				  if (!strcmp(massager->propertyNames[0],"d")) 
					{
					  massager->outputgrid[i].eclOutInfo.statsConstr.useStochasticDepth = 1;
					}
				}
			  OGPoints = massager->outputgrid[i].Points; // useful shortcut
			  for (j=0;j<massager->outputgrid[i].nPoints;j++)
				{
				  // write out x,y
				  if (massager->writeAsciiOut) 
					{
					  fprintf(massager->ascii_out_realisations_fp[i],"%12.4f %12.4f ",
							  massager->outputgrid[i].xValues[j],
							  massager->outputgrid[i].yValues[j]);
					  /* add utm stuff if need be */
					  if (massager->utmInfo.performOutputGridUTMmapping)
						{
						  fprintf(massager->ascii_out_realisations_fp[i],"%12.4f %12.4f ",OGPoints[j][0], OGPoints[j][1]);
						}
					  if (!massager->duplicate)
						{
						  /* and eclio indexes for eclipse */
						  // Use:  ecl_i, ecl_j, seg
						  fprintf(massager->ascii_out_realisations_fp[i],"%6d %6d ",
								  massager->outputgrid[i].ecl_i[j], massager->outputgrid[i].ecl_j[j]); 
						  if (massager->outputgrid[i].segnum)
							{
							  // Add segments: seg
							  fprintf(massager->ascii_out_realisations_fp[i],"%4d ",
									  massager->outputgrid[i].segnum[j]);
							}
						}
					}

				  // Now scale to couple in property correlations within and across layers
				  Lscale = getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(massager,i,j,0,0);

				  for (real_no=0; real_no < massager->Nrealisations; real_no++)
					{
					  // Form realisation by y = P.L.w
					  // y_sum = L.w bit
					  for (m=0; m< massager->nPropertyNames; m++)
						{
						  properties_indx_k = m + i*massager->nPropertyNames;
						  // Multiply across row of Cholesky prefactor. Looks messy but that's all it is
						  // w = simulationValues = correlated vector of normal-deviates
						  sum=0;
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
						  // Debugging: to skip Cross layer/property Cholesky, do...
						  // y_sum[properties_indx_k] = simulationValues[i][j][m][real_no];
						  // else (correct)
						  y_sum[properties_indx_k] = sum;
						}

					  // Add back in smoothed means, and modulate by well kriging scaling factors
					  for (m=0; m< massager->nPropertyNames; m++)
						{
						  mean = krigMeans[i][j][m];
						  properties_indx_k = m + i*massager->nPropertyNames;

						  y_sum[properties_indx_k] *= sqrt(MAX(krigVariances[i][j][m],0));   // To get well data right

						  // the realisation is now y_sum[properties_indx_k]+ krigMeans[i][j][m]
						  val = y_sum[properties_indx_k]+ mean;
						  val = MAX(0,val);
						  // normalisation
						  if (massager->properties[m].normalised) val = MIN(1.0,val);
						  // possible net-sand <= thickness truncation
						  if ((m == massager->net_sand_index) && (massager->thickness_index>=0))
							{
							  // we are simulating net_sand: thickness_index simulation is in:
							  thickness_val = sim_values[massager->thickness_index + i*massager->nPropertyNames];
							  // truncate:
							  val = MIN(thickness_val,val);
							}
						  // hard-coded override on NG if thickness+net-sand available:
						  if ((m == massager->NG_index) && (massager->thickness_index>=0) && (massager->net_sand_index>=0))
							{
							  thickness_val = sim_values[massager->thickness_index + i*massager->nPropertyNames];
							  net_sand_val = sim_values[massager->net_sand_index + i*massager->nPropertyNames];
							  val = (thickness_val>0) ? net_sand_val/thickness_val : 0.0;
							}
						  sim_values[properties_indx_k] = val;

						  // write out/store the actual realisation value 
						  if (isnan(val))
							{
							  printf("BUG2\n");
							}
						  if (makeStochasticEclipse) eclOutInfo->r_properties[real_no][i][j][m] = val;

						  if (massager->writeAsciiOut) 
							{
							  fprintf(massager->ascii_out_realisations_fp[i],"%12.3f ",val);
							}
						  // acquire stats. Sums 
						  if (dimStats)
							{
							  if (massager->acquireStats==1)
								 {
								   // depth is available
								   if (massager->outputgrid[i].eclOutInfo.statsConstr.useStochasticDepth)
									 {
									   if (m==0) z = zStoch = val;
									   else z = zStoch;
									 }
								   else
									 {
									   z = massager->outputgrid[i].zValues[j];
									 }
								 }
							  // Either no depth info available (acquireStats==2), or depth is in-range
							  if ((massager->acquireStats==2) || (( z>= massager->outputgrid[i].eclOutInfo.statsConstr.Zmin)
								  && (z <= massager->outputgrid[i].eclOutInfo.statsConstr.Zmax)))
								{
								  if (massager->duplicate) 
									{
									  area = massager->inputgrid.traceArea;    
									}
								  else
									{
									  area = massager->outputgrid[i].areas[j];
									}
								  if (val > massager->statsCutoff)
									{
									  qty = val * area;
									  stats[real_no][m] += qty;              // individual layer stats
									  layer_sum_stats[real_no][m] += qty;    // composite layer stats
									}
								}
							}
						} // property loop
					  if (decorate_this_layer)
						{
						  dec_NG = sim_values[massager->decoratorInfo.NG_index];
						  dec_thick = sim_values[massager->decoratorInfo.thickness_index];
						  dec_d = sim_values[massager->decoratorInfo.d_index];

						  // save realisation target for later. Note this is in jth column by natural 
						  // (not multigrid ordering)
						  // ensure no NG values >1 get fed in
						  massager->decoratorInfo.NGtarget[real_no][j] = MIN(1,dec_NG);
						  // test wedge....
						  // massager->decoratorInfo.NGtarget[real_no][j] = (float)j/massager->outputgrid[i].nPoints;
						  massager->decoratorInfo.thickness[real_no][j] = dec_thick;
						  massager->decoratorInfo.depth[real_no][j] = dec_d;
						}
					} // realisations loop
				  if (massager->writeAsciiOut) fprintf(massager->ascii_out_realisations_fp[i],"\n");
				  free_dmatrix(Lscale,0,massager->nProperties-1,0,massager->nProperties-1);
				  Lscale=NULL;
				}  // loop over nodes j

			  if (massager->writeAsciiOut) 
				{
				  fclose(massager->ascii_out_realisations_fp[i]);
				  printf("Wrote: %s\n",massager->ascii_out_realisations_filenames[i]);
				  free(massager->ascii_out_realisations_filenames[i]);
				  massager->ascii_out_realisations_filenames[i]=NULL;
				}

			  // Statistics output
			  if (dimStats) 
				{
				  // individual layer stats
				  showStatistics(massager,stats,0,massager->layerNumbers[i],massager->Nrealisations,dimStats,massager->propertyNames);
				  // composite layer stats on last layer
				  if (i==(massager->nLayers-1))
					{
					  showStatistics(massager,layer_sum_stats,1,0,massager->Nrealisations,dimStats,massager->propertyNames);
					}
				}


			  // Decoration for each realisation....
			  if (decorate_this_layer) 
				{
				  for (real_no =0; real_no < massager->Nrealisations; real_no++)
					{
					  for (i_dec=0; i_dec<massager->nDecorations; i_dec++)
						{
						  for (j=0;j<massager->outputgrid[i].nPoints;j++)
							{
							  // get target from multigrid index 
							  output_grid_indx = path[i][j];  
							  error = dec3DSimulate(massager,i,j,massager->decoratorInfo.NGtarget[real_no][output_grid_indx]);
							  // massager->decoratorInfo.NGerror[real_no][output_grid_indx] = error;
							}
						  decOutputSimulation(massager,real_no,i_dec);
						}
					}
				}
			} // loop over layers 
		  free_dvector(y_sum,0,massager->nProperties-1); y_sum=NULL;
		  free_dvector(sim_values,0,massager->nProperties-1); sim_values=NULL;
		  free_dvector(w_tmp,0,massager->nProperties-1); w_tmp=NULL;
		}

	  // After all layers simulated, we can write eclipse realisations
	  if (makeStochasticEclipse)
		{
		  for (real_no=0; real_no < massager->Nrealisations; real_no++)
			{
			  writeEclipseRealisation(massager,real_no,eclOutInfo->r_properties[real_no]);
			}
		}
	}

  // Cleanup
  // 1 ECLIO
  if (makeStochasticEclipse)
	{
	  freeStochEclipseMem(massager);
	}


  // 2 smoothed means and variances
  for (i=0; i< massager->nLayers ; i++) 
	{
	  free(path[i]); path[i]=NULL;
	  size = massager->outputgrid[i].nPoints;
	  free_matrix(krigMeans[i],0,size-1,0,massager->nPropertyNames-1);  krigMeans[i]=NULL;
	  free_matrix(krigVariances[i],0,size-1,0,massager->nPropertyNames-1); krigVariances[i]=NULL;
	  free_matrix(krigStdDev[i],0,size-1,0,massager->nPropertyNames-1); krigStdDev[i]=NULL;
	  free_matrix(sI->smoothMeans[i], 0, massager->outputgrid[i].nPoints-1, 0, massager->nPropertyNames-1);
	  sI->smoothMeans[i]=NULL;
	}
  free(krigMeans);
  krigMeans=NULL;
  free(krigVariances);
  krigVariances=NULL;
  free(krigStdDev);
  krigStdDev=NULL;
  free(sI->smoothMeans);
  sI->smoothMeans = NULL;
  free(path);
  path=NULL;


#if !USE_IN_PLACE_SMOOTHING
  // Free new space with smoothed covariances
  for (i=0; i < (massager->inputgrid.nPoints); i++)
	{
	  free (massager->solverInfo.smoothPropertyCovariance_tri[i]);
	  massager->solverInfo.smoothPropertyCovariance_tri[i] = NULL;
	}
#endif
  free(massager->solverInfo.smoothPropertyCovariance_tri);
  massager->solverInfo.smoothPropertyCovariance_tri = NULL;



  if (massager->Nrealisations > 0)
	{
	  for (i=0; i< massager->nLayers ; i++) 
		{
		  free_f3tensor(simulationValues[i],0,massager->outputgrid[i].nPoints - 1,0,
						massager->nPropertyNames-1,0,massager->Nrealisations-1);
		  simulationValues[i]=NULL;
		}
	  free(simulationValues);
	  simulationValues = NULL;
	  if (dimStats>0)
		{
		  free_matrix(stats,0,dimStats-1,0,massager->Nrealisations-1);
		  stats=NULL;
		  free_matrix(layer_sum_stats,0,dimStats-1,0,massager->Nrealisations-1);
		  layer_sum_stats=NULL;
		}
	}
  if (massager->nDecorations) free_decorator_info(massager);

  return EXIT_SUCCESS;
}


int dumpPostMassagedSUFiles(struct MASSAGER *massager,float ***krigMeans,float ***krigStdDev)
{
  int i,j,k,m,mp,size,ns,n_properties,nbytes,nbytes_w;
  int Xindx, Yindx, InlineIndx=0, XlineIndx=0;
  int i_x, i_y, in_line, x_line;

  double thickness, vp_eff; 
  double** OGPoints;
  segy trace;
  struct SU_OUT_INFO *suinfo;
  struct OUTPUTGRID *mog_i;
  char *su_filename;
  FILE* fSU;

  printf("Writing SU trend volumes\n");
  asprintf(&su_filename,"post_massaged_properties.su");

  if ((fSU = fopen(su_filename,"wb")) == NULL)
	{
	  warning ("Could not write SU post-massaging file %s. Abandoning SU output.",su_filename);
	  goto FAIL;
	}

  if (dumpBHPviewerPropertyKeyFile(massager,"post_massaged_property_keys.txt",1) == EXIT_FAILURE)
	{
	  goto FAIL;
	} 
  printf("Wrote SU property-key file post_massaged_property_keys.txt\n");

  suinfo = massager->su_out_info;
  /* Zero trace data */
  memset( (void *) trace.data, (int) '\0',SU_NFLTS*sizeof(float));

  n_properties = suinfo->nBodyPropertiesPossibleForSUOutput;  // sensible body properties
  n_properties += massager->nPropertyNames;                   // + all std deviations
  ns = (suinfo->nLayersPossibleForSUOutput * n_properties)    // blocks for each layer
	  + (suinfo->nLayersPossibleForSUOutput+1);               // final "time" block


  // Since grids are duplicated, have guaranteed that
  // massager->inputgrid.nPoints = massager->outputgrid[i].nPoints;
  size = massager->inputgrid.nPoints;
	  

  // Header words to try and set
  if (massager->utmInfo.performOutputGridUTMmapping)
	{
	  // in_line/x_line supplied, to get x,y via UTM triple
	  InlineIndx = getindex(massager->inputgrid.xWord);
	  XlineIndx = getindex(massager->inputgrid.yWord);
	  // set gx,gy for want of better choice
	  Xindx = getindex("gx");
	  Yindx = getindex("gy");
	}
  else
	{
	  // have x, y directly. Nothing can really be inferred about in_line/x_line
	  Xindx = getindex(massager->inputgrid.xWord);
	  Yindx = getindex(massager->inputgrid.yWord);
	}

  // write out in trace (location) order
  for (j=0; j < size; j++)
	{
	  mog_i = &(massager->outputgrid[0]);

	  // Set header words
	  trace.ns = ns;
	  // many SU programs like this set, even if it's meaningless
	  trace.dt = 1000;
	  // FIXME gx should be massager->inputgrid.xWord etc
	  /* add utm stuff if need be */ 
	  if (massager->utmInfo.performOutputGridUTMmapping)
		{
		  i_x = (int)floor(mog_i->Points[j][0] + 0.5); // round
		  i_y = (int)floor(mog_i->Points[j][1] + 0.5); // round
		  in_line = (int) (mog_i->xValues[j]);
		  x_line = (int) (mog_i->yValues[j]);

		  puthval(&trace, Xindx, &i_x);
		  puthval(&trace, Yindx, &i_y);
		  puthval(&trace, InlineIndx, &in_line);
		  puthval(&trace, XlineIndx, &x_line);
		}
	  else
		{
		  i_x = (int) floor(mog_i->xValues[j] + 0.5); // round
		  i_y = (int) floor(mog_i->yValues[j] + 0.5);  // round
		  puthval(&trace, Xindx, &i_x);
		  puthval(&trace, Yindx, &i_y);
		}

	  // MEANS blocks
	  k=0;
	  // run over sensible body properties in this layer
	  for (m=0; m< massager->nPropertyNames; m++) 
		{
		  if (indexIsADepictableSUBodyProperty(massager,m))
			{
			  // A body SU property
			  // loop on layers
			  for (i=0; i< suinfo->nLayersPossibleForSUOutput ; i++)
				{
				  // write out +ve truncated kriged mean
				  trace.data[k] = MAX(krigMeans[i][j][m],0);
				  k++;
				} // layers
			} // Body property
		} // loop on properties


	  // STDDEV blocks
	  // run over body properties in this layer
	  for (m=0; m< massager->nPropertyNames; m++) 
		{
		  for (i=0; i< suinfo->nLayersPossibleForSUOutput ; i++,k++)
			{
			  // write out +ve truncated kriged stddev, all properties
			  trace.data[k] = MAX(krigStdDev[i][j][m],0);
			} // layers
		} // properties

	  // Now fill in last block of times etc
	  m = (suinfo->inTime) ? massager->t_index : massager->d_index;
	  for (i=0; i< suinfo->nLayersPossibleForSUOutput; i++,k++)
		{
		  // layer top time or depth
		  trace.data[k] = MAX(krigMeans[i][j][m],0);
		}
	  if (suinfo->getLastSurfacefromTbase)
		{
		  // Straight from t_base, last layer. Must be in time, by definition
		  mp = massager->t_base_index;
		  i = suinfo->nLayersPossibleForSUOutput - 1;
		  trace.data[k] = MAX(krigMeans[i][j][mp],0);
		}
	  else
		{
		  // Computed final bottom time or depth
		  if (!suinfo->getLastSurfacefromThickness)
			{
			  // Omitted last layer, so to get layer top time or depth - simply read from last layer
			  i = suinfo->nLayersPossibleForSUOutput;
			  trace.data[k] = MAX(krigMeans[i][j][m],0);
			}
		  else
			{
			  // compute from previous layer
			  i = suinfo->nLayersPossibleForSUOutput - 1;
			  if (suinfo->inTime)
				{
				  // in time. Expect to have thickness & vp_eff available
				  vp_eff = MAX(krigMeans[i][j][massager->vp_eff_index],0);
				  thickness = MAX(krigMeans[i][j][massager->thickness_index],0);
				  // Factor of 2 since times are two-way times.
				  // Delivery's times are written out as ms, so massaging will preserve this
				  // Thus factor of 1000 so produce ms.
				  // velocities have no Greek prefixes I think (i.e. presume lengths in thickness/vp_eff cancel)
				  trace.data[k] = trace.data[k-1] + ((vp_eff > 0) ? (2.0e3*thickness/vp_eff ) : 0.0);
				}
			  else
				{
				  // In depth. Have thickness as well
				  thickness = MAX(krigMeans[i][j][massager->thickness_index],0);
				  // Factor of 2 since times are two-way times
				  trace.data[k] = trace.data[k-1] + thickness;				  
				}
			}
		}

	  // Dump traces. 
	  nbytes = trace.ns*sizeof(float);
	  // Flip bytes if needed
	  swap_to_BigEndian_if_needed(&trace);

	  // Write out..Poor man's version of SU stream handlers etc.
	  nbytes_w = fwrite(&trace,1, HDRBYTES, fSU);           // header
	  nbytes_w += fwrite(&(trace.data[0]),1,nbytes, fSU);	// data
	  if (nbytes_w != (HDRBYTES+nbytes))
		{
		  warning ("Could not write complete SU data. Abandoning SU output.");
		  return EXIT_FAILURE;
		}		  	  
	} // loop over grid/locations/traces

  fclose(fSU);
  printf("Wrote SU volume %s\n",su_filename);

  return EXIT_SUCCESS;

 FAIL:
  free(su_filename);
  return EXIT_FAILURE;

}

double **getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(struct MASSAGER *massager,
																  int layer_index, int output_grid_index, 
																  int useRobustCholesky,double eps_last)
{
  int i,j,k,m,mm,info=0,IGindex,N;
  double wt,eps=0.0,max,sum;
  struct FLPAIR* krigPairs;
  struct SOLVER_INFO* sI;

  i = layer_index;
  j = output_grid_index;
  N = massager->nProperties; 
  sI = &massager->solverInfo;  // Useful shortcut

  // Full covariance scaling
  double **interpolatedPropCovar = dmatrix(0,N-1,0,N-1);

  // interpolate the input grid full covariance into the current output-grid point
  // linear combination implicit in kriging guarantees positive definiteness, if weights
  // are +ve. Some kriging structures induce -ve weights, so, to be sure, we use wt^2.
  // This is the correct coefficient for the sum of n "random processes" anyway
  sum=0;
  krigPairs = sI->OGkrigInfoSoft[i][j].values;
  for (k=0; k< sI->OGkrigInfoSoft[i][j].size; k++) sum += krigPairs[k].val * krigPairs[k].val;

  for (m=0; m<N; m++)
	{
	  for (mm = 0; mm<=m; mm++)
		{
		  // interpolate smooth propertyCovariance from 
		  // input grid to this output grid site, using precomputed kriging weights
		  interpolatedPropCovar[m][mm] = 0;
		  // krigPairs = sI->OGkrigInfoSoft[i][j].values;
		  for (k=0; k< sI->OGkrigInfoSoft[i][j].size; k++)
			{
			  IGindex = krigPairs[k].indx;
			  wt = krigPairs[k].val;
			  // simplest hack for -ve kriging weights....
			  wt = wt*wt/sum;
			  interpolatedPropCovar[m][mm]  += wt * get(sI->smoothPropertyCovariance_tri[IGindex],mm,m);
			  // symmetry
			  interpolatedPropCovar[mm][m] = interpolatedPropCovar[m][mm]; 
			}
		}
	}
  // Now get eigenfactor U such that U.U^T = interpolatedPropCovar
  // Ansatz: Cholesky is a decent surrogate

  // Make this robust against zeros:
  // iterated call
#if 0
  if (eps_last>0) 
	{
	  // Don't know what to do for sure here. Keep winding up eps perhaps???
	  eps=2*eps_last;
 	  //print_matrix((void**)interpolatedPropCovar,N,'d',"dodgy_matrix.txt");
	  //exit(0);
	}
  else
	{
	  max=0;
	  for (m=0; m<N; m++) max = MAX(max,interpolatedPropCovar[m][m]);
	  eps = POSITIVE_DEFINITE_MATRIX_EPS * max;
	  if (eps<=0)
		{
		  warning ("Zero covariance matrix. Setting eps to +ve small.\n");	
		  eps=1e-20;
		}
	}
#endif

  // preclean matrices to ensure no floating point fluff
  for (i=0;i<N;i++) 
	{
	  if (interpolatedPropCovar[i][i]<=0)
		{
		  for (j=0;j<N;j++) 
			{
			  interpolatedPropCovar[i][j] = interpolatedPropCovar[j][i] = 0.0;
			}
		}
	}

#if 0
  // debugging  
  if (layer_index==1)
	{
	  printf("NODE %f %f ",massager->inputgrid.xValues[output_grid_index],massager->inputgrid.yValues[output_grid_index]);
	  for (i=0;i<N;i++) {
		for (j=i;j<N;j++) printf ("%g ", interpolatedPropCovar[i][j]);
	  }
	  printf("\n");
	}
#endif

#if 0
  // debugging: approximate Cholesky by sqrts of diagonal....
  for (i=0;i<N;i++) 
	{
	  interpolatedPropCovar[i][i] = sqrt(interpolatedPropCovar[i][i]);
	  for (j=0;j<i;j++) interpolatedPropCovar[i][j]=0;
	  for (j=i+1;j<N;j++) interpolatedPropCovar[i][j]=0;
	}
  return interpolatedPropCovar;
#endif

  if (useRobustCholesky)
	{
	  getCholeskySafe(interpolatedPropCovar,N);
	} 
  else
	{
	  // fast but fragile.....
	  // info = clapack_wrapper_dpotrf(CblasLower,N,&interpolatedPropCovar[0][0],N);

	  // Use bullet-proof method every time....
	  /* Modified Cholesky solve of SK system */
	  getCholeskySafe(interpolatedPropCovar,N);
	} 

  if (info>0) 
	{
	  // Fail detected. Use safe Cholesky
	  if (!useRobustCholesky)
		{
		  warning ("ATLAS lapack Cholesky decomposition failed in interpolatedPropCovar(). Using safe Cholesky\n");
		}
	  free_dmatrix(interpolatedPropCovar,0,N-1,0,N-1);
	  interpolatedPropCovar=NULL;

	  // try again, with doubled eps, using useRobustCholesky=1
	  return getMatrixFactorOfCompleteSmoothedLayerPropertyCovariance(massager,layer_index,output_grid_index,1,eps);
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
  
#if 0
  for (j=0; j < size; j++) 
	{
	  printf("MG traverse: %d %d\n",mgIndices[j].val,mgIndices[j].indx);
	}
#endif
  
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
			  "Smoothed+kriged summary file for layer %d, generated by massager run from XML file %s\n",
			  massager->layerNumbers[layerNumber],massager->massage_xml_file);
	  if (!massager->duplicate)
		{
		  // Eclipse.
		  if (massager->outputgrid[layerNumber].segnum)
			{
			  // Use:  x, y, ecl_i, ecl_j, seg
			  fprintf(massager->ascii_out_fp[layerNumber],"%d\nx\ny\necl_i\necl_j\nseg\n",5+2*massager->nPropertyNames);
			}
		  else
			{
			  // Use:  x, y, ecl_i, ecl_j
			  fprintf(massager->ascii_out_fp[layerNumber],"%d\nx\ny\necl_i\necl_j\n",4+2*massager->nPropertyNames);
			}			  
		}
	  else
		{
		  // Duplicated
		  if (massager->utmInfo.utmTripleAvailable)
			{
			  // Duplicated, and UTM's available. Use:  inline, xline, x, y,
			  fprintf(massager->ascii_out_fp[layerNumber],"%d\n%s\n%s\nx\ny\n",4+2*massager->nPropertyNames,
					  massager->inputgrid.xWord,massager->inputgrid.yWord);
			}
		  else
			{
			  // Duplicated, but UTM's not available. Use:  xword, yword,
			  fprintf(massager->ascii_out_fp[layerNumber],"%d\n%s\n%s\n",2+2*massager->nPropertyNames,
					  massager->inputgrid.xWord,massager->inputgrid.yWord);
			}
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
			  "Posterior realisations file for layer %d, generated by massager run from XML file %s\n",
			  massager->layerNumbers[layerNumber],massager->massage_xml_file);
	  if (!massager->duplicate)
		{
		  // Eclipse. 
		  if (massager->outputgrid[layerNumber].segnum)
			{
			  // Use:  x, y, ecl_i, ecl_j, seg
			  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\nx\ny\necl_i\necl_j\nseg\n",
					  5+massager->Nrealisations*massager->nPropertyNames);
			}
		  else 
			{
			  // Use:  x, y, ecl_i, ecl_j
			  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\nx\ny\necl_i\necl_j\n",
					  4+massager->Nrealisations*massager->nPropertyNames);
			}
		}
	  else
		{
		  // Duplicated
		  if (massager->utmInfo.utmTripleAvailable)
			{
			  // Use inline, xline, x, y
			  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\n%s\n%s\nx\ny\n",
					  4+massager->Nrealisations*massager->nPropertyNames,
					  massager->inputgrid.xWord,massager->inputgrid.yWord);
			}
		  else
			{
			  // Use xword, yword
			  fprintf(massager->ascii_out_realisations_fp[layerNumber],"%d\n%s\n%s\n",
					  2+massager->Nrealisations*massager->nPropertyNames,
					  massager->inputgrid.xWord,massager->inputgrid.yWord);
			}
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

/*
  Compute smoothed mean and covariances for each property on the 
  output grid (from input-grid (delivery) points, using
  the trend smoothing variograms. No well data constraints applied here.
  Smoothing of means respects segments if required.
*/
int computeMovingAverage(struct MASSAGER *massager)
{
  int i,ii,j,jj,k,i_y,kp,kpp,properties_indx_kp,properties_indx_kpp,OGindx,size;
  int ecl_i,ecl_j;
  int ecl_ref_layer,massager_depth_layer_indx=-1;
  double zSeg=0;
  // Smooth
  int IGindx,segnum=0,segSafe=1,pass,useSegWeights,useExtrapWeights;
  int block_segnum=0;
  int inputGridHasSomeDataInSameSegment = 1;
  int willUseSegments=0;
  int MAsegWarnings=0;
  int isNAN;
  float *fullKrigSum=NULL;
  double wt,segWt,weight,extrapWtval,sum,segSum,sumSq,pSum;
  double **Points;
  double pt[2];
  struct PT2D pt_utm;
  int nNeighborsDesired;
  int *neigh;
  float cov_val,s_val;
  struct OUTPUTGRID *mog_i;

  /* Special variogram smoothing variables */
  char handle_special_variogram;
  double wt_special=0, segWt_special=0 ,sum_special=0 ,segSum_special=0;


  struct SOLVER_INFO *sI = &massager->solverInfo;  // Useful shortcut

  printf("Computing moving trends for mean and covariance\n");

  if (massager->utmInfo.utmTripleAvailable)
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}

  // Use prebuilt kdtree of inputgrid from setUpMapping()

  // alloc room for means (and maybe std deviations)
  // struct mallinfo mall_info = mallinfo();
  show_malloc_stats();
  
  sI->smoothMeans = (float***)malloc((size_t) ((massager->nLayers) * sizeof(float**)));



  sI->smoothPropertyCovariance_tri = (float**)malloc((size_t) ((massager->inputgrid.nPoints) * sizeof(float*)));
  if (!sI->smoothPropertyCovariance_tri) fprintf(stderr,"allocation failure of smoothPropertyCovariance_tri");
  sI->smoothPropertyCovariance_tri_size = massager->inputgrid.nPoints;

#if !USE_IN_PLACE_SMOOTHING
  // careful use of space here!

  for (i=0; i < (massager->inputgrid.nPoints); i++)
	{
	  if (massager->inputgrid.isUsedByInterpolation[i])
		{			  
		  sI->smoothPropertyCovariance_tri[i] = getUpperTriMatrix(massager->nProperties,massager->diagonalCovariancesOnly);
		}
	  else
		{
		  sI->smoothPropertyCovariance_tri[i] = NULL;
		}
	} 
#endif // USE_IN_PLACE_SMOOTHING


  for (i=0; i< massager->nLayers ; i++)
	{ 
	  sI->smoothMeans[i] = matrix(0,massager->outputgrid[i].nPoints-1, 0, massager->nPropertyNames-1);
	  // zero it
	  memset(&(sI->smoothMeans[i][0][0]), 0,massager->outputgrid[i].nPoints *  massager->nPropertyNames * sizeof(float));
	}

  show_malloc_stats();  

  if (massager->eclInfo.useECLIOsegments) /* will only enter if not duplicating */
	{
	  if (massager->verbosity>2) printf("Segmenting Delivery input grid from Eclipse grid...\n");
	  massager->inputgrid.segnum = ivector(0,massager->inputgrid.nPoints-1);
	  massager->inputgrid.availableSegments = g_hash_table_new(g_int_hash,g_int_equal);

	  willUseSegments = 0;
	  for (kp=0; kp<massager->nPropertyNames; kp++) 
		{
		  if (!massager->properties[kp].smooth_across_faults) willUseSegments = 1;
		}

	  // set segmentation by segments of reference depth layer if possible
	  if (massager->d_index >= 0)  // are massaging depth
		{
		  ecl_ref_layer = massager->outputgrid[0].ECLIPSE_reference_depth_layer_number;
		  for (ii=0; ii<massager->nLayers; ii++)
			{
			  if (massager->outputgrid[ii].ECLIPSE_layer_number == ecl_ref_layer)
				{
				  massager_depth_layer_indx = ii;
				}
			}
		}

	  // set up inputgrid segments by in-polygon? tests of output grid.
	  nNeighborsDesired = 1;

	  // perform snapping on reference layer i if at all possible, or if not available, layer 0	  
	  i = (massager_depth_layer_indx>=0) ? massager_depth_layer_indx : 0;
	  
	  mog_i = &(massager->outputgrid[i]);
	  Points = massager->outputgrid[i].Points;

	  /* loop on all IG grid nodes in this layer */
	  /* by this stage massager->inputgrid.Points will be in UTM's if at all possible */
	  for (j=0; j < massager->inputgrid.nPoints; j++) 
		{
		  // get neighbour list from kdtree
		  neigh = kdOptNNQuery(Points,2,&(massager->inputgrid.Points[j][0]),nNeighborsDesired,
							   EUCLIDEAN,mog_i->OptkdRoot,1);

#if DUPLICATE_AT_FAULTS
		  // Make assigment based on (1) eclipse block segment of eclipse block containing (x,y,z), where
		  // x,y are from delivery grid, z is from ref layer (or layer 0 if ref layer not available)
		  // (2) if this fails, use fallback as follows:
		  // Get (arbitrarily) lowest seg number in (possible) set of duplicates containing this OG point neigh[1]
		  // mog_i->segnum[jj] = segInfo[ecl_i[jj]][ecl_i[jj]].values[i_dup];
		  ecl_i = mog_i->ecl_i[neigh[1]];
		  ecl_j = mog_i->ecl_j[neigh[1]];
		  // search back to start of group of duplicated points, if any
		  for (jj=neigh[1]-1;jj>=0;jj--) if (mog_i->ecl_i[jj] != ecl_i || mog_i->ecl_j[jj] != ecl_j) break;
		  if (jj!=0) jj++; 
		  // assign to minimum of group
		  segnum = mog_i->segnum[jj];
		  for (;jj< massager->inputgrid.nPoints; jj++)
			{
			  if (mog_i->ecl_i[jj] != ecl_i || mog_i->ecl_j[jj] != ecl_j) break;
			  segnum = IMIN(segnum,mog_i->segnum[jj]);
			}
		  if (jj != massager->inputgrid.nPoints-1) jj--;

		  // find segment of eclipse block, as obtained at depth of ref layer if possible
		  if (massager_depth_layer_indx>=0)
			{
			  // depth of middle of ref. layer at this nearest point
			  zSeg = massager->outputgrid[massager_depth_layer_indx].zValues[neigh[1]];
			}
		  else
			{
			  // use current layer depth
			  zSeg = mog_i->zValues[jj];
			}
		  block_segnum = getSegOfEclipseBlock(&(mog_i->ecl), mog_i->segcol, ecl_i, ecl_j, zSeg, &(massager->inputgrid.Points[j][0]));
		  
		  if (block_segnum != I_MISSING) 
			{
			  massager->inputgrid.segnum[j] = block_segnum;
			}
		  else
			{
			  massager->inputgrid.segnum[j] = segnum;
			}
#else
		  // simple nearest-nbr snapping
		  massager->inputgrid.segnum[j] = mog_i->segnum[neigh[1]];
#endif

		  // record available segment numbers in inputgrid via hashmap
		  if (!g_hash_table_lookup(massager->inputgrid.availableSegments, &(massager->inputgrid.segnum[j])))
			{
			  g_hash_table_insert(massager->inputgrid.availableSegments,
								  &(massager->inputgrid.segnum[j]),&(massager->inputgrid.segnum[j]));
			}
		  if (massager->verbosity>4) 
			{
			  printf("IGSEG (%d, ecl_i,j= %d %d) %15f %15f %d\n", j, ecl_i,ecl_j,
					 massager->inputgrid.Points[j][0], massager->inputgrid.Points[j][1],massager->inputgrid.segnum[j]);
			}
		  free(neigh);
		  neigh=NULL;
		}
	  fullKrigSum = vector(0,massager->nPropertyNames-1);
	}
  else  
	{
	  massager->inputgrid.segnum = NULL;
	}


  nNeighborsDesired = MIN(massager->inputgrid.nPoints,
						  DEFAULT_MAX_NEIGHBOURS_IN_MOVING_AVERAGE);
  Points = massager->inputgrid.Points;


  // (A) First, smooth input grid properties to outputgrid means, respecting segments

  

  // Test if exceptional treatment of smoothing of some trend properties is required
  handle_special_variogram = 0;  
  for (kp=0; kp<massager->nPropertyNames; kp++)
	{
	  if (massager->properties[kp].use_trend_smoothing_special_variogram)  
		{
		  handle_special_variogram = 1; 
		  setUpcorr2D(&massager->trend_smoothing_special_variogram);
		  break;
		}
	}

  setUpcorr2D(&massager->trend_smoothing_variogram);

  for (i=0; i< massager->nLayers ; i++)
	{
	  // set up appropriate metric for the neighbour search using a kdtree
	  // the neighbour search metric is defined by the trend-smoothing variogram
	  setUp2DEuclidenGenMetric(massager->inputgrid.OptkdRoot->Kglobals,
							   M_PI*massager->trend_smoothing_variogram.azimuth/180.0,
							   massager->trend_smoothing_variogram.rangeX/massager->trend_smoothing_variogram.rangeY);

	  int properties_indx_offset = i*massager->nPropertyNames;

	  /* loop on all OG grid nodes in this layer */
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  // x,y of OG point in this layer
		  pt[0] = massager->outputgrid[i].Points[j][0];
		  pt[1] = massager->outputgrid[i].Points[j][1];

		  /* segment of this OG point */
		  if (massager->eclInfo.useECLIOsegments) 
			{
			  segnum = massager->outputgrid[i].segnum[j];
			  // is this a point with segment label which can never be matched by any input grid segments?
			  if (g_hash_table_lookup(massager->inputgrid.availableSegments, (gpointer)&segnum) != NULL)
				inputGridHasSomeDataInSameSegment = 1; 
			  else inputGridHasSomeDataInSameSegment = 0;
			} 

		  // find set of neighbours on input grid 
		  // get neighbour list from kdtree
		  neigh = kdOptNNQuery(Points,2,&pt[0],nNeighborsDesired,
							   EUCLIDEAN_GEN,massager->inputgrid.OptkdRoot,1);

	
		  // Form weighted average on this neighborhood, weights from variogram
		  if (massager->eclInfo.useECLIOsegments)
			{
			  memset(&(fullKrigSum[0]), 0, massager->nPropertyNames * sizeof(float));
			}
		  if (handle_special_variogram) 
			{
			  sum_special=0.0;
			  segSum_special = 0.0;
			}
		  for (i_y=1,sum=0,segSum=0; i_y <= nNeighborsDesired; i_y++)
			{
			  IGindx = neigh[i_y];
			  // form weight from correlation function
			  wt = corr2Dexp(pt[0],pt[1],
							 Points[IGindx][0],Points[IGindx][1],
							 &(massager->trend_smoothing_variogram));
			  // preclude numerical underflows on smoothing...
			  wt = MAX(wt,MIN_SMOOTHING_WEIGHT);
			  segWt = wt;
			  for (kp=0,isNAN=0; kp<massager->nPropertyNames; kp++)
				{
				  if (massager->property_NAN[IGindx]) 
					{
					  isNAN = 1;
					  segWt = wt = 0.0;
					}
				}
			  sum += wt;

			  /* special smoothing variogram */
			  if (handle_special_variogram)
				{
				  // Make approximation here that main smoothing variogram 
				  // neighbourhood list is adequate for special variogram
				  wt_special = corr2Dexp(pt[0],pt[1],
								 Points[IGindx][0],Points[IGindx][1],
								 &(massager->trend_smoothing_special_variogram));
				  // preclude numerical underflows on smoothing...
				  wt_special = MAX(wt_special,MIN_SMOOTHING_WEIGHT);
				  segWt_special = wt_special;
				  if (isNAN) segWt_special = wt_special = 0.0;
				  sum_special += wt_special;
				}


			  /* zero-weight neighbour if from wrong segment */
			  if (massager->eclInfo.useECLIOsegments)  /* will only enter if not duplicating */
				{
				  if ((massager->inputgrid.segnum[IGindx] != segnum) && inputGridHasSomeDataInSameSegment) 
					{
					  segWt=0.0;
					  if (handle_special_variogram) segWt_special=0.0;
					}


				  segSum += segWt;				  
				  if (handle_special_variogram) segSum_special += segWt_special;

				  // compute full krig average as backup
				  if (!isNAN)
					{
					  for (kp=0; kp<massager->nPropertyNames; kp++)
						{
						  properties_indx_kp = kp + properties_indx_offset;
						  if (!(massager->properties[kp].use_trend_smoothing_special_variogram))
							{
							  fullKrigSum[kp] += wt * massager->properties[properties_indx_kp].values[IGindx];
							}
						  else
							{
							  fullKrigSum[kp] += wt_special * massager->properties[properties_indx_kp].values[IGindx];
							}
						}
					}
				}
			  // segment-adjusted kriging
			  if (!isNAN)
				{
				  // simplest case
				  for (kp=0; kp<massager->nPropertyNames; kp++)
					{
					  properties_indx_kp = kp + properties_indx_offset;
					  if (massager->properties[kp].smooth_across_faults)
						{
						  weight = (massager->properties[kp].use_trend_smoothing_special_variogram) ? wt_special : wt ; 
						  massager->solverInfo.smoothMeans[i][j][kp] += weight * massager->properties[properties_indx_kp].values[IGindx];
						}
					  else 
						{
						  weight = (massager->properties[kp].use_trend_smoothing_special_variogram) ? segWt_special : segWt ; 
						  massager->solverInfo.smoothMeans[i][j][kp] += weight * massager->properties[properties_indx_kp].values[IGindx];
						}
					}
				  
				}
			} // loop over neighbours

		  if (massager->eclInfo.useECLIOsegments)  /* will only enter if not duplicating */
			{
			  if ((segSum==0) && willUseSegments)
				{ 
				  // Big trouble. No neighbours. Don't bother handling special variogram 
				  // smoothing here, as this is so pathological...
				  if (massager->verbosity>1 && MAsegWarnings <= MAXWARNINGS)
					{
					  printf("WARNING: no moving average neighbours in same segment!!! "
							 "Moving average is suspect. (x,y,seg) = %g %g %d\n",
							 pt[0],pt[1],segnum);
					  if (MAsegWarnings == MAXWARNINGS)
						{
						  printf("Suppressing further moving average warnings...\n"); 
						}
					  MAsegWarnings++;
					}
				  // segment-adjusted kriging had no decent neighbours. Restore full kriging as stopgap if this is ok
				  if (sum>MA_WEIGHT_SUM_EPS)
					{
					  for (kp=0; kp<massager->nPropertyNames; kp++)
						{
						  massager->solverInfo.smoothMeans[i][j][kp] = fullKrigSum[kp];
						}
					}
				  else
					{
					  // probably here we have extrapolation underflow issues in the variograms. Even using all the
					  // neighbours didn't work. Default to
					  // average over all neighbours with equal weight
					  sum = nNeighborsDesired;
					  for (kp=0; kp<massager->nPropertyNames; kp++)
						{
						  fullKrigSum[kp]=0.0;
						  properties_indx_kp = kp + properties_indx_offset;
						  for (i_y=1; i_y <= nNeighborsDesired; i_y++)
							{
							  IGindx = neigh[i_y];
							  fullKrigSum[kp] += massager->properties[properties_indx_kp].values[IGindx];
							}
						  massager->solverInfo.smoothMeans[i][j][kp] = fullKrigSum[kp];
						}
					} 
				  segSum = sum;
				  if (handle_special_variogram) segSum_special = sum;
				  segSafe=0;
				}
			  else segSafe=1;
			}
		  else 
			{
			  segSum = sum;
			  if (handle_special_variogram) segSum_special = sum_special;
			}
		  // renormalise
		  for (kp=0; kp<massager->nPropertyNames; kp++) 
			{
			  if (massager->properties[kp].smooth_across_faults)
				{
				  weight = (massager->properties[kp].use_trend_smoothing_special_variogram) ? sum_special : sum; 
				  massager->solverInfo.smoothMeans[i][j][kp] /= weight;
				}
			  else 
				{
				  weight = (massager->properties[kp].use_trend_smoothing_special_variogram) ? segSum_special : segSum; 
				  massager->solverInfo.smoothMeans[i][j][kp] /= weight;
				}
			}
		  if (massager->verbosity>4) // 4
			{
			  printf("MA %g %g ",massager->outputgrid[i].xValues[j],massager->outputgrid[i].yValues[j]);
			  for (kp=0; kp<massager->nPropertyNames; kp++)
				{		
				  properties_indx_kp = kp + properties_indx_offset;	
				  printf("%g ", massager->solverInfo.smoothMeans[i][j][kp]);
				}
			  printf("%d %d\n",segSafe,segnum);
			}
		  free(neigh);  // neighbourhood of jth observation
		  neigh=NULL;
		}
	  KillMetricSpace(massager->inputgrid.OptkdRoot->Kglobals);
	}

  if (massager->eclInfo.useECLIOsegments) 
	{
	  free_vector(fullKrigSum,0,massager->nPropertyNames-1);
	  fullKrigSum =NULL;
	}

  // (A - footnote) At this point, another useful pre-computational step is to 
  // compute and store the smooth surfaces at the well observation points. This has nothing to do with
  // smoothing the means or the variances. Well data is not used in the smoothing
  sI = &massager->solverInfo;  // Useful shortcut

  int *obsCounters = ivector(0,massager->nPropertyNames-1); 
  for (i=0; i<  massager->nLayers ; i++)
	{
	  for (kp=0; kp< massager->nPropertyNames; kp++) obsCounters[kp]=0;
	  for (j=0; j<massager->well_observations[i].n; j++)
		{
		  kp = g_array_index(massager->well_observations[i].propertyIndex, gint,j);
		  if (kp==-1) continue; /* a non-modelled property */
		  // interpolate smooth map of kp-th property onto well-observation point for this jth observation
		  sum=0;
		  /* a modelled property observation */		  
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
  obsCounters = NULL;





  // (B) Smooth property covariances on input grid also. 

  // set up appropriate metric for the kriging neighbour search using a kdtree
  setUp2DEuclidenGenMetric(massager->inputgrid.OptkdRoot->Kglobals,
						   M_PI*massager->trend_smoothing_variogram.azimuth/180.0,
						   massager->trend_smoothing_variogram.rangeX/massager->trend_smoothing_variogram.rangeY);

  size = massager->nPropertyNames * massager->nLayers;

  /* loop on all IG grid nodes in this layer */
  /* IMPORTANT: Don't segment here (hence the 0) as (1) it's not necessary 
	 for covariances and (2) if you do, weighting across properties has to 
	 be done carefully to preserve positive definiteness
  */
  int useSegmentationForVariances = 0;


	  



#if USE_IN_PLACE_SMOOTHING

  /************ In-place code here assumes useSegmentationForVariances = 0  ******************/

  // don't alloc here, just point at originals. Smoothed versions will get copied in place
  // careful use of space here!
  for (i=0; i < (massager->inputgrid.nPoints); i++)
	{		  
	  sI->smoothPropertyCovariance_tri[i] = massager->property_covariance_tri[i];
	}

  // Lookups

  struct FLPAIR **smLookups;
  smLookups = (struct FLPAIR **)malloc((size_t) (massager->inputgrid.nPoints * sizeof(struct FLPAIR*)));
  if (!smLookups) 
	{
	  fprintf(stderr,("smLookups allocation failure in computeMovingAverage()"));
	  return EXIT_FAILURE;
	}
  smLookups[0] = (struct FLPAIR *)malloc((size_t) (massager->inputgrid.nPoints * nNeighborsDesired * sizeof(struct FLPAIR)));
  if (!smLookups[0]) 
	{
	  fprintf(stderr,("smLookups[0] allocation failure in computeMovingAverage()"));
	  return EXIT_FAILURE;
	}
  for(j=1;j<massager->inputgrid.nPoints ;j++) 
	{
	  smLookups[j] = smLookups[j-1]+nNeighborsDesired;
	}

  show_malloc_stats();

  if (massager->verbosity>2)
	{
	  printf("Getting covariance smoother lookups\n");
	}
  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  // skip smoothing if point won't be used (and is not alloced anyway)
	  if (!massager->inputgrid.isUsedByInterpolation[j]) continue;

	  // x,y of IG point 
	  pt[0] = Points[j][0];
	  pt[1] = Points[j][1];

	  // find set of neighbours on input grid
	  // get neighbour list from kdtree
	  neigh = kdOptNNQuery(Points,
						   2,&pt[0],nNeighborsDesired,
						   EUCLIDEAN_GEN,massager->inputgrid.OptkdRoot,1);

	  // Form weighted average on this neighborhood, weights from variogram (squared)
	  // use second pass only if no neighbours of the same segment are found
	  for (i_y=1,sum=0,segSum=0; i_y <= nNeighborsDesired; i_y++)
		{
		  IGindx = neigh[i_y];
		  // form weight from correlation function
		  wt = corr2Dexp(pt[0],pt[1],
						 Points[IGindx][0],Points[IGindx][1],
						 &(massager->trend_smoothing_variogram));
		  wt=wt*wt;   // summing covariances must be squared (see notes for sum of indep random vectors)
		  segWt = wt;
		  sum += wt;

		  // simply store what we need....
		  smLookups[j][i_y-1].indx = IGindx;
		  smLookups[j][i_y-1].val = segWt;
		}
	  segSum=sum;
	  free(neigh);  // neighbourhood of jth input grid point
	}  // Lookup setups.

  // Order-flipped in-place smoothing, looping over properties slowest
  // Allocate temp storage for [i][ii] matrix element, for all grid nodes
  float*  tmp_covar = vector(0,massager->inputgrid.nPoints-1);

  for (i=0; i<size; i++)
	{
	  if (massager->verbosity>2)
		{
		  printf("Smoothing covariance row %d...\n",i);
		}
	  for (ii=i; ii< (massager->diagonalCovariancesOnly ? (i+1) : size); ii++)
		{
		  // set up fast matrix access
		  set_quick_access(i,ii);
		  // zero accumulator
		  memset(tmp_covar,0,(size_t)(massager->inputgrid.nPoints*sizeof(float)));
		  // average
		  for (j=0; j < massager->inputgrid.nPoints; j++) 
			{
			  // skip smoothing if point won't be used (and is not alloced anyway)
			  
			  if (!massager->inputgrid.isUsedByInterpolation[j]) continue;

			  // Form weighted average on this neighborhood, weights from variogram (squared)
			  for (i_y=1,sum=0; i_y <= nNeighborsDesired; i_y++)
				{
				  IGindx = smLookups[j][i_y-1].indx;
				  wt = smLookups[j][i_y-1].val;				  
				  // next line accesses [i][ii] th element via set_quick_access(i,ii) call above
				  cov_val = get_quick(massager->property_covariance_tri[IGindx]); 
				  // cov_val = get(massager->property_covariance_tri[smLookups[j][i_y-1].indx],i,ii); // [i][ii] th element
				  // if (isnan(cov_val)) {wt=cov_val=0.0;}  // FIXME Would MUCH rather clean up NANs in input reading!!!!
				  if (massager->property_NAN[IGindx]) {wt=cov_val=0.0;} 
				  sum += wt;
				  tmp_covar[j] += wt * cov_val;
				}
			  // renormalise at this trace
			  tmp_covar[j] /= ((sum!=0) ? sum : 1);
			}
		  // now copy [i][ii] th property variance back, since it's no longer needed
		  for (j=0; j < massager->inputgrid.nPoints; j++) 
			{
			  if (!massager->inputgrid.isUsedByInterpolation[j]) continue;
			  if (i!=ii) {
				set_quick(massager->property_covariance_tri[j], tmp_covar[j]); // [i][ii] th element
			  } else {
				set_quick_roundup(massager->property_covariance_tri[j], tmp_covar[j]); // [i][i] th element
			  }
			}
		  if (massager->ignore_soft)
			{
			  for (j=0; j < massager->inputgrid.nPoints; j++) 
				{
				  set(massager->property_covariance_tri[j], (i==ii) ? 1 : 0, i, ii);
				}
			}
		}
	}
  free_vector(tmp_covar,0,massager->inputgrid.nPoints-1);
  tmp_covar= NULL;
  free(smLookups[0]); smLookups[0] = NULL;
  free(smLookups); smLookups = NULL;



#else  // #if USE_IN_PLACE_SMOOTHING

  // Working version where new covariances are alloced separately

  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  // skip smoothing if point won't be used (and is not alloced anyway)
	  if (!massager->inputgrid.isUsedByInterpolation[j]) continue;

	  if (massager->verbosity>2)
		{
		  if (j%1000==0) 
			{
			  printf("Smoothing covariance at node %d...\n",j);
			  // show_malloc_stats();
			}
		}
	  // x,y of IG point
	  pt[0] = Points[j][0];
	  pt[1] = Points[j][1];

	  // find set of neighbours on input grid
	  // get neighbour list from kdtree
	  neigh = kdOptNNQuery(Points,
						   2,&pt[0],nNeighborsDesired,
						   EUCLIDEAN_GEN,massager->inputgrid.OptkdRoot,1);

	  /* segment of this IG point */
	  /* IMPORTANT: Don't segment here (hence the 0) as (1) it's not necessary 
		 for covariances and (2) if you do, weighting across properties has to 
		 be done carefully to preserve positive definiteness
	  */
	  if (useSegmentationForVariances) 
		{
		  segnum = massager->inputgrid.segnum[j];
		  useSegWeights = 1;
		}
	  useExtrapWeights = 0;
	  // Form weighted average on this neighborhood, weights from variogram (squared)
	  // use second pass only if no neighbours of the same segment are found
	  for (pass=0; pass<=1; pass++)
		{
		  for (i_y=1,sum=0,segSum=0; i_y <= nNeighborsDesired; i_y++)
			{
			  IGindx = neigh[i_y];
			  // form weight from correlation function
			  wt = corr2Dexp(pt[0],pt[1],
							 Points[IGindx][0],Points[IGindx][1],
							 &(massager->trend_smoothing_variogram));
			  wt=wt*wt;   // summing covariances must be squared (see notes for sum of indep random vectors)
			  segWt = wt;
			  sum += wt;
			  /* zero-weight neighbour if from wrong segment */
			  if (useSegmentationForVariances)
				{
				  if ((massager->inputgrid.segnum[IGindx] != segnum) && useSegWeights) segWt=0.0;
				  segSum += segWt;
				}
			  if (useExtrapWeights) segWt = 1.0; // bulletproofing if required

			  for (i=0; i<size; i++)
				{
				  for (ii=i; ii< (massager->diagonalCovariancesOnly ? (i+1) : size); ii++)
					{
					  cov_val = get(massager->property_covariance_tri[IGindx],i,ii);
					  accumulate(massager->solverInfo.smoothPropertyCovariance_tri[j],segWt * cov_val,i,ii);
					}
				}
			}
		  // underflow protection for wild extrapolation
		  if (pass==0 && sum==0) useExtrapWeights = 1;
		  if (useSegmentationForVariances)
			{
			  if ((pass=0) && (segSum==0.0)) useSegWeights=0;  // will force recalculation
			  else break;
			}
		  else 
			{
			  segSum=sum;
			  break;  // non segmented mode
			}
		}
	  if (massager->verbosity>4)
		{
		  printf("MovAvCovariance %g %g ",pt[0],pt[1]);
		}

	  if (useExtrapWeights) segSum = nNeighborsDesired;

	  // renormalise at this trace
	  for (i=0; i<size; i++)
		{
		  for (ii=i; ii< size; ii++)
			{
			  divide(massager->solverInfo.smoothPropertyCovariance_tri[j],segSum,i,ii);
			  if (massager->verbosity>4)
				{
				  printf("%g ",get(massager->solverInfo.smoothPropertyCovariance_tri[j],i,ii));
				}
			}
		}
	  if (massager->ignore_soft)
		{
		  for (i=0; i<size; i++)
			{
			  for (ii=i; ii< size; ii++)
				{
				  set( massager->solverInfo.smoothPropertyCovariance_tri[j], (i==ii) ? 1 : 0, i, ii);
				}
			}
		}
#if 0
	  //test fragment for pos definiteness. Must comment out! Only works for conventional storage!!!
	  double max=0;
	  float** cov = massager->solverInfo.smoothPropertyCovariance[j];
	  // print_matrix(cov,massager->nProperties,'f',NULL);
	  for (i=0;i<(massager->nProperties);i++) max=MAX(max,cov[i][i]);
	  for (i=0;i<(massager->nProperties);i++) cov[i][i] += 1e-6*max;
	  int info = clapack_wrapper_spotrf(CblasLower,size,&(massager->solverInfo.smoothPropertyCovariance[j][0][0]),size);
  
	  if (info>0) 
		{
		  // fail detected. 
		  printf("j=%d: smoothPropertyCovariance not pos definite!\n",j);
		  exit(0);
		}
#endif

	  if (massager->verbosity>4) printf("\n");

	  // Could precompute and store Cholesky factor here for realisations.....
	  free(neigh);  // neighbourhood of jth input grid point
    }

#endif // #if USE_IN_PLACE_SMOOTHING


  KillMetricSpace(massager->inputgrid.OptkdRoot->Kglobals);

  KillOptTree(massager->inputgrid.OptkdRoot);

  if (massager->eclInfo.useECLIOsegments)
	{
	  free_ivector(massager->inputgrid.segnum, 0, massager->inputgrid.nPoints-1);
	  massager->inputgrid.segnum = NULL;
	  g_hash_table_destroy(massager->inputgrid.availableSegments);
	}


  if (massager->duplicate)
	{
#if !USE_IN_PLACE_SMOOTHING
	  // can clean up covariances, but must leave x,y's in place
	  freeMassageAnalysisCovariances(massager);
#endif
	}
  else
	{ 

#if USE_IN_PLACE_SMOOTHING
	  // dispose of unneeded original covariances at certain traces to save memory
	  // fold all memory up front so cleaning leaves contig. array
	  // pointless to do this if duplicating
	  printf("Defragging memory...");
	  j = massager->inputgrid.nPoints-1;
	  for (i=0; i< massager->inputgrid.nPoints; i++)
		{
		  // find an empty slot
		  if (!massager->inputgrid.isUsedByInterpolation[i])
			{
			  // search from end for a required matrix, and copy it into the empty i slot.
			  // adjust pointers, free original, and go on
			  for (;j>=0;j--)
				{
				  if (massager->inputgrid.isUsedByInterpolation[j]) break;
				}
			  if (j>=0 && j>i)
				{
				  // copyUpperTriMatrix( float* Adest,float* Asrc,int n)
				  copyUpperTriMatrix(massager->property_covariance_tri[i], massager->property_covariance_tri[j],massager->nProperties);
				  freeUpperTriMatrix(massager->property_covariance_tri[j],massager->nProperties);
				  massager->property_covariance_tri[j] = NULL;
				  massager->property_covariance_tri[j] = massager->property_covariance_tri[i];	  
				  sI->smoothPropertyCovariance_tri[j] = massager->property_covariance_tri[i];		
				  j--; 
				}
			} 
		}
	  massager->property_covariance_tri_size = j+1;
	  sI->smoothPropertyCovariance_tri_size = massager->property_covariance_tri_size;
	  printf("done\n");

#else

	  freeMassageAnalysisData(massager);

#endif
	}

  // free up isUsedByInterpolation flag as this is the last we'll need of it...
  free_cvector(massager->inputgrid.isUsedByInterpolation ,0,massager->inputgrid.nPoints-1);
  massager->inputgrid.isUsedByInterpolation = NULL;

  printf("Finished smoothing of means and covariances\n");

#if 0  
  // Debugging. Print smoothed covars to file

  FILE *fsm;
  fsm = fopen("smooth_covars.txt","w");
  if (!fsm) {printf("Cannot open output file smooth_covars.txt for writing\n");exit(0);}
  for (j=0; j < massager->inputgrid.nPoints; j++) 
	{
	  fprintf(fsm,"%f %f ",massager->inputgrid.xValues[j],massager->inputgrid.yValues[j]);

	  for (i=0; i<size; i++)
		{
		  for (ii=i; ii< size; ii++)
			{
			  fprintf(fsm,"%10.4g ",get(massager->solverInfo.smoothPropertyCovariance_tri[j], i, ii));
			}
		}
	  fprintf(fsm,"\n"); 
	}
  printf("Wrote smooth_covars.txt\n");
  fclose(fsm);	 

#endif

  show_malloc_stats(); 

  return EXIT_SUCCESS;
}

/*
  Compute means and std deviations of the properties. Print CDF, sort
  files by property 0. Supply the argument 
  isComposite==1 for composite statistics
*/
void showStatistics(struct MASSAGER *massager, float **stats, int isComposite, int layerNumber,int Nsamples, int Nprops,char **names)
{
  int i,j,k;
  float mean,var;
  struct PAIR *pairs;
  FILE* fp;
  char *s;

  if (isComposite)
	{
	  asprintf(&s,"Realisation_Summary_Stats_Composite.txt");
	}
  else
	{
	  asprintf(&s,"Realisation_Summary_Stats_layer_%d.txt",layerNumber);
	}

  fp = fopen(s,"w");
  if (fp == NULL)
	{
	  printf("Cannot open output file %s for writing\n",s);
	  free(s);
	  return;
	}

  pairs = (struct PAIR*)malloc((size_t) (Nsamples * sizeof(struct PAIR)));

  for (i = 0; i<Nprops; i++)
	{
	  mean = 0;
	  for (j = 0;j<Nsamples;j++)
		{
		  if (!strcmp(names[i], massager->statsSortPropName))
			{
			  pairs[j].indx = j;
			  pairs[j].val = stats[j][i];
			}
		  mean+=stats[j][i];
		}
	  mean/=Nsamples;
	  var = 0;
	  for (j = 0;j<Nsamples;j++)
		{
		  var+=SQR(stats[j][i]-mean);
		}
	  var = (Nsamples==1) ? 0 : var/(Nsamples-1);
	  if (strcmp(names[i],"d") && strcmp(names[i],"NG")) 
		{
		  fprintf(fp,"Stats for area-weighted, depth-constrained sum of property %s: mean=%.4g, std dev=%.4g, ratio %.4g\n",
				  names[i],mean,sqrt(var),sqrt(var)/mean);
		}
	}

  // sort and rename realisations files
  qsort(pairs,(size_t)Nsamples,sizeof(struct PAIR),pairCmp);

  // write summary
  fprintf(fp,"Cumulative CDF realisation indices, ranked by %s\n",massager->statsSortPropName);
  fprintf(fp,"Percentile        sum Realisation_index");
  // other properties
  for (k=0; k<Nprops; k++)
	{
	  if (strcmp(names[k],"net-sand") && strcmp(names[k],"net-hydrocarbon") 
		  && strcmp(names[k],"d") && strcmp(names[k],"NG"))
		{
		  fprintf(fp,"  %10s",names[k]);
		}
	}
  fprintf(fp,"\n");
  for (j=0;j<Nsamples;j++)
	{
	  fprintf(fp,"%10.2f %10.5g %17d",(100.0*j)/(float)Nsamples, pairs[j].val, pairs[j].indx);
	  for (k=0; k<Nprops; k++)
		{
		  if (strcmp(names[k],"net-sand") && strcmp(names[k],"net-hydrocarbon")
 			  && strcmp(names[k],"d") && strcmp(names[k],"NG")) 
			{
			  fprintf(fp," %10.5lg",stats[pairs[j].indx][k]);
			}
		}
	  fprintf(fp,"\n");
	} 
  fclose(fp);

  printf("Wrote summary stats to %s\n",s); 
  free(s);
  free(pairs);
}

int pairCmp(const void *a, const void *b)
{
  struct PAIR *A = (struct PAIR *)a;
  struct PAIR *B = (struct PAIR *)b;

  return (A->val > B->val)? 1 : ((A->val < B->val)? -1 : 0 );
}

/*
  Solve OK system of form {{Ckrig,1},{1^T,0}}.{{w},{mu}} = {{RHS},{1}}
  using safe Cholesky tricks. Dimensions of Ckrig are nxn.
  Returns solution in-place in RHS. Solution returns only first n components (w)
  of answer.
  Ckrig and RHS can be dimensioned either (n+1)x(n+1) or n x n.
  They must be zero offset.
*/
void solveOKsystem(double** Ckrig,double* RHS,int n)
{
  int i,*perm;
  double sum=0.0, sumR=0.0, mu, *CkrigOne=NULL;

  perm = ivector(0,n-1);
  CkrigOne = dvector(0,n-1);

  for (i=0;i<n;i++) CkrigOne[i]=1.0;

  modified_Cholesky_left_decomposition(n,Ckrig,perm,0);
  modified_Cholesky_solve(n,Ckrig,perm,RHS);
  modified_Cholesky_solve(n,Ckrig,perm,CkrigOne);

  for (i=0;i<n;i++) 
	{
	  sum += CkrigOne[i];
	  sumR += RHS[i];
	}
  mu = (sumR - 1.0)/sum;
  for (i=0;i<n;i++) 
	{
	  RHS[i] -= mu * CkrigOne[i];
	}
  free_ivector(perm,0,n-1); perm=NULL;
  free_dvector(CkrigOne,0,n-1); CkrigOne=NULL;
}

/* A useful debugging func */
void print_matrix(void **m,int n,char type,char* name)
{
	int i,j;
	FILE *out;
	char *filename;
	
	if (name==NULL) asprintf(&filename,"junk.txt");
	else filename=name;
	printf("Dumping matrix to %s\n",filename);
	if ((out = fopen(filename, "w")) == NULL) 
	  {
		printf("Can't open %s\n",filename);
		exit(-1);
	  }
	for (i=0;i<n;i++)
	  {
		for (j=0;j<n;j++)
		  {
			fprintf(out,"%15.15f\n",(double)((type=='f') ? ((float**)m)[i][j] : ((double**)m)[i][j]));
		  }
	  }
	fclose(out);
	if (name==NULL) free(filename);
}
