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
#include <sys/time.h>

#include <glib.h>

#include "NR/nrutil.h"
// #include "NR/nr.h"


 
#include "kdtree/nn.h"
#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h" 

#include "errorHandler.h"


#include "massager_types.h"
#include "variogram.h"
#include "solver.h"
#include "decorator.h" 
#include "util.h"  
#include "modified_left_Cholesky.h"  
#include "clapack_cblas_wrappers.h"
#include "massager.h"

#define DEFAULT_MAX_NEIGHBOURS_IN_TGS_SIMULATIONS    20  // Number of neighbours used in sequential simulation
#define DEFAULT_N_COLUMN_TRIES 100                       // Number of repeat column tries used in sequential simulation
#define USE_STRAT_COORDS 0                               // Do all decoration in strat coordinates (lengths rescaled to unit box)
#define TEST_TREE_TIMING 1                               // To test 3Dtree timing
#define TREAT_Z_COORD_AS_CORNER_PT 1                     // To treat k values at corner point grid locations. Last layer is then not used in NG test


void initialiseDecorator(struct DECORATOR_INFO *decoratorInfo)
{
  decoratorInfo->Nz=0;                 // grid subdivision
  decoratorInfo->decoration_layer=0;   // which massager layer to decorate

  decoratorInfo->z_stretch=0;
  decoratorInfo->Nattempts=0;
  decoratorInfo->t_tree_timing=0; // 3D kdtree timing
  
  memset(&decoratorInfo->variogram.rotmat[0][0],0,9*sizeof(double));

  // set this to flag bad parsing...
  decoratorInfo->decoration_layer_indx=-1;
  decoratorInfo->NG_index=0;
  decoratorInfo->thickness_index=0;
  decoratorInfo->d_index=0; 
  decoratorInfo->NGtarget=NULL;

  // float** NGerror;
  decoratorInfo->thickness=NULL;
  decoratorInfo->depth=NULL;
  decoratorInfo->seed=-1;             // Seed used for realisations
  
  decoratorInfo->col_weights=NULL;

  decoratorInfo->threeDpath=NULL;
  decoratorInfo->simulation=NULL;
  decoratorInfo->facies_simulation=NULL;
  
  // Stuff for use in kd trees
  decoratorInfo->dTree=NULL;
  decoratorInfo->Points=NULL;       // 3D grid for use in KD tree
  decoratorInfo->flag=NULL;         // 3D flag for known hard data
  decoratorInfo->OptkdRoot=NULL;    // in KD tree
}



int setUp3Dsimulation(struct MASSAGER *massager,int *path)
{
  int i,j,k,m,size; 
  int *threeDpath;
  double z_scale;
  dynamicTree *dTree;
  struct OUTPUTGRID *outputgrid;
  struct DECORATOR_INFO *decoratorInfo; 
  struct VARIOGRAM_3D *variogram;

 
  // shortcuts
  decoratorInfo = &massager->decoratorInfo;
  outputgrid = &massager->outputgrid[decoratorInfo->decoration_layer_indx];
  variogram = &decoratorInfo->variogram;

  massager->decoratorInfo.Nattempts = DEFAULT_N_COLUMN_TRIES;

  // set up variogram
  setUpcorr3D(variogram);


  // build points array for this grid: last index 0,1,2 = x,y,z
  size = outputgrid->nPoints * decoratorInfo->Nz;
  decoratorInfo->simulation = dvector(0,size - 1);
  decoratorInfo->facies_simulation = (short*)malloc((size_t) (size * sizeof(short)));
  decoratorInfo->Points = dmatrix(0,size- 1, 0,2);
  decoratorInfo->flag = (short*)malloc((size_t) (size * sizeof(short)));


  outputgrid->xMin = outputgrid->xMax = outputgrid->Points[0][0];
  outputgrid->yMin = outputgrid->yMax = outputgrid->Points[0][1];

  for (i=0; i<outputgrid->nPoints; i++)
	{
	  outputgrid->xMin = MIN(outputgrid->xMin, outputgrid->Points[i][0]);
	  outputgrid->xMax = MAX(outputgrid->xMax, outputgrid->Points[i][0]);
	  outputgrid->yMin = MIN(outputgrid->yMin, outputgrid->Points[i][1]);
	  outputgrid->yMax = MAX(outputgrid->yMax, outputgrid->Points[i][1]);
	}

  // For efficiency in searches, we stretch z so the variogram ellipsoid 
  // becomes approximately a sphere
  // (1) Z-stretch scaled from grid to make approximate cube
  massager->decoratorInfo.z_stretch = sqrt(SQ(outputgrid->xMax - outputgrid->xMin)
										   +SQ(outputgrid->yMax - outputgrid->yMin));

#if USE_STRAT_COORDS

  double DX,DY;
  DX = outputgrid->xMax -outputgrid->xMin;
  DY = outputgrid->yMax -outputgrid->yMin;
  DX = (DX>0)?DX:1;
  DY = (DY>0)?DY:1;
  z_scale = ((decoratorInfo->Nz>1)? decoratorInfo->Nz-1 : 1);
  variogram->a_hmax /= massager->decoratorInfo.z_stretch;
  variogram->a_hmin /= massager->decoratorInfo.z_stretch;

#else

  z_scale = massager->decoratorInfo.z_stretch /
	((decoratorInfo->Nz>1)? decoratorInfo->Nz-1 : 1);

  // grid points are now approx equally spaced in this cube.
  // massager->decoratorInfo.variogram.a_vert *= massager->decoratorInfo.z_stretch;

  // (2) Then magnify by correlation ratio coeff 0.01????
  // massager->decoratorInfo.z_stretch *= .01* variogram->a_hmax/variogram->a_vert; 

#endif

  for (i=0,j=0; i<outputgrid->nPoints; i++)
	{
	  for (k=0;k<decoratorInfo->Nz; k++,j++)
		{
#if USE_STRAT_COORDS
		  // x, y, from strats
		  decoratorInfo->Points[j][0] = (outputgrid->Points[i][0] - outputgrid->xMin)/DX;
		  decoratorInfo->Points[j][1] = (outputgrid->Points[i][1] - outputgrid->yMin)/DY;
		  // stratigraphic depth w in [0,1]
		  decoratorInfo->Points[j][2] =((double)k)*z_scale;
#else
		  // x, y, inherited
		  decoratorInfo->Points[j][0] = outputgrid->Points[i][0];
		  decoratorInfo->Points[j][1] = outputgrid->Points[i][1];
		  // stratigraphic depth w in [0,1]
		  decoratorInfo->Points[j][2] =((double)k)*z_scale;
#endif
		  // scan hard data and set to 1 if known... FIXME (fill in later)
		  decoratorInfo->flag[j]=0;
		}
	}

  // Precompute and store multigrid paths
  // 1D path for column
  struct INTPAIR* mgIndices = (struct INTPAIR*)malloc((size_t) ((decoratorInfo->Nz * sizeof(struct INTPAIR))));
  formMultigridTraverseOrder(mgIndices,decoratorInfo->Nz);

  // Compute interlaced 3Dpath 
  // this is a sequence of chunks on the same ordering as the 2dpath for the columns, but
  // each chunk corresponds to a multigrid ordering of the z gridblocks
  decoratorInfo->threeDpath = (int*)malloc((size_t) ((size * sizeof(int))));

  threeDpath = decoratorInfo->threeDpath;
  for (k=0; k < size; k++) 
	{
	  // index = 2D offset (multigrid) + column index (multigrid)....
	  threeDpath[k] = decoratorInfo->Nz*path[k/decoratorInfo->Nz] + mgIndices[k%decoratorInfo->Nz].val;
	}
  free(mgIndices);

  // Dynamic KD tree for decorator tree
  decoratorInfo->dTree = buildDynamicTree(decoratorInfo->Points, size, 3, threeDpath);

  setUp3DEuclidenGenMetric(decoratorInfo->dTree->OptkdRoot->Kglobals, 
						   variogram->a_hmax, variogram->a_hmin, variogram->a_vert, 
						   variogram->ang1, variogram->ang2, variogram->ang3);

  decoratorInfo->col_weights = (struct COLUMN_TGS_KRIG_WEIGHTS **)
	malloc((size_t)((outputgrid->nPoints * sizeof(struct COLUMN_TGS_KRIG_WEIGHTS*))));

  for (i=0; i<outputgrid->nPoints; i++)
	{
	  decoratorInfo->col_weights[i] = (struct COLUMN_TGS_KRIG_WEIGHTS *)
		malloc((size_t)((decoratorInfo->Nz * 
						 sizeof(struct COLUMN_TGS_KRIG_WEIGHTS))));
	}


  /* Lookup NG targets for realisations */
  decoratorInfo->NGtarget = matrix(0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  // decoratorInfo->NGerror = matrix(0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  decoratorInfo->thickness = matrix(0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  decoratorInfo->depth = matrix(0,massager->Nrealisations-1,0,outputgrid->nPoints-1);

  decoratorInfo->t_tree_timing=0;
  
  return EXIT_SUCCESS; 
}


int decSetUp3DcolumnSim(struct MASSAGER *massager,int layer_number, int column_number)
{
  int i,j,k,kk,k_s,kk_s,nKnown,firstTry;
  int output_grid_indx,nNeighborsDesired;
  int deleteRow=0,coLocatedNbrIndx=0;
  int* neigh,info;
  double** Ckrig;
  double *RHSkrig,*Crhs;
  int *ipiv_krig;
  double condVar;
  struct DECORATOR_INFO *decoratorInfo; 
  struct VARIOGRAM_3D *variogram;
  struct COLUMN_TGS_KRIG_WEIGHTS *col_weights;
  double** Points;
  double dt_timing=0;
  double Cscale;

  i = layer_number;
  j = column_number;

  decoratorInfo = &massager->decoratorInfo;
  Points = decoratorInfo->Points;

  k_s = j*decoratorInfo->Nz;

  nKnown = k_s; // +hard data FIXME

  // decoratorInfo->variogram.nugget=0;
  Cscale = decoratorInfo->variogram.nugget + decoratorInfo->variogram.sill;

  for (kk_s=0; kk_s<decoratorInfo->Nz; k_s++,kk_s++,nKnown++)
	{
	  if (decoratorInfo->flag[j] == 0)
		{
		  // an unknown point. Must simulate/set up simulation

		  output_grid_indx = decoratorInfo->threeDpath[k_s];   // multigrid path

		  // Construct SK system from previously simulated nodes
		  // We solve a SK system for each property at this grid location
		  // The structure of this system is the same for all properties, since
		  // the accumulating neighbours are the same for all properties, the path is the same etc,
		  // and the variogram is taken to be the same
		  nNeighborsDesired = MIN(nKnown,DEFAULT_MAX_NEIGHBOURS_IN_TGS_SIMULATIONS);
		  col_weights = &decoratorInfo->col_weights[column_number][kk_s];
		  col_weights->N = nNeighborsDesired;

		  // find previsited nearest neighbours 
		  struct timeval tv;
		  struct timezone tz;
		  double told;
		  gettimeofday(&tv,&tz);
#if TEST_TREE_TIMING
		  told=(tv.tv_sec*1e6 + tv.tv_usec);
#endif
		  neigh = nextDynamicKdOptNNQuery(nNeighborsDesired, EUCLIDEAN_GEN, decoratorInfo->dTree);
#if TEST_TREE_TIMING
	  // time accesses
		  gettimeofday(&tv,&tz);
		  dt_timing += (double)((tv.tv_sec*1e6+tv.tv_usec) - told);
#endif
		  if (neigh==NULL && (nKnown!=0))
			{
			  printf("TGS BUG!\n"); 
			}
		  if (nNeighborsDesired > 0)
			{
			  col_weights->wts = (double*)malloc((size_t) (nNeighborsDesired * sizeof(double)));
			  // col_weights->nbr_indices = (int*)malloc((size_t) (nNeighborsDesired * sizeof(int)));
			  col_weights->nbr_pointers = (double**)malloc((size_t) (nNeighborsDesired * sizeof(double*)));
			  // neigh[1....nNeighborsDesired] are now the indexes 
			  // of the already simulated points on Points
			  Ckrig = dmatrix(0,nNeighborsDesired-1,0,nNeighborsDesired-1);
			  RHSkrig = decoratorInfo->col_weights[column_number][kk_s].wts;
			  Crhs = dvector(0,nNeighborsDesired-1);
			  ipiv_krig = ivector(0,nNeighborsDesired-1);
			  coLocatedNbrIndx=-1;
			  firstTry=1;

			SK:  // start of assembly of SK matrix and RHS

			  // form SK kriging system over these
			  for (k=0; k < nNeighborsDesired; k++) 
				{
				  // decoratorInfo->col_weights[column_number][kk_s].nbr_indices[k] = neigh[k+1];
				  col_weights->nbr_pointers[k] = &(decoratorInfo->simulation[neigh[k+1]]);
				  // Matrix. Need form only lower triangle
				  for (kk=0; kk <= k ; kk++) 
					{
					  // Points is in UTM's if at all possible
					  Ckrig[k][kk] = corr3D(Points[neigh[k+1]][0], 
											Points[neigh[k+1]][1], Points[neigh[k+1]][2],
											Points[neigh[kk+1]][0], Points[neigh[kk+1]][1], 
											Points[neigh[kk+1]][2],
											&(decoratorInfo->variogram));
					  Ckrig[k][kk] /= Cscale;      // normalise to guarantee unit variance
					  Ckrig[kk][k] = Ckrig[k][kk]; // FIXME needed?
					}
				  // RHS.
				  RHSkrig[k] = corr3D(Points[output_grid_indx][0],
									  Points[output_grid_indx][1],
									  Points[output_grid_indx][2],
									  Points[neigh[k+1]][0], 
									  Points[neigh[k+1]][1], 
									  Points[neigh[k+1]][2], 
									  &(decoratorInfo->variogram));
				  RHSkrig[k] /= Cscale;            // normalise to guarantee unit variance
				  if (RHSkrig[k]>=1.0) coLocatedNbrIndx = k;
				}

			  cleanKrigingSystemOfDuplicateNeighbours(nNeighborsDesired, Ckrig, RHSkrig);

			  memcpy(Crhs,RHSkrig,nNeighborsDesired*sizeof(double));
			  condVar = 1.0;  // store this before we mess up Ckrig

			  // Now solve kriging system for the weights
			  if (coLocatedNbrIndx>=0)
				{
				  // one of the neighbours is co-located with the current point. The kriging weight must be one
				  for (k=0; k < nNeighborsDesired; k++) RHSkrig[k]=0;
				  RHSkrig[coLocatedNbrIndx]=1.0;
				  condVar=0.0;
				}
			  else
				{
				  if (firstTry)
					{
					  info = clapack_wrapper_dgetrf(nNeighborsDesired,nNeighborsDesired,
											&Ckrig[0][0],nNeighborsDesired,&ipiv_krig[0]);
					  if (info>0)
						{
						  warning ("lapack decomposition failed in TGS kriging step: using safe Cholesky\n");
						  firstTry=0;
						  goto SK;
						}
					  /* Get weights via dge solution solution */
					  clapack_wrapper_dgetrs('N',nNeighborsDesired,1,&Ckrig[0][0],
									 nNeighborsDesired,&ipiv_krig[0],&RHSkrig[0],nNeighborsDesired);
					}
				  else
					{
					  /* Modified Cholesky solve of SK system */
					  int *perm = ivector(0,nNeighborsDesired-1);
					  modified_Cholesky_left_decomposition(nNeighborsDesired,Ckrig,perm,0);
					  modified_Cholesky_solve(nNeighborsDesired,Ckrig,perm,RHSkrig);
					  free_ivector(perm,0,nNeighborsDesired-1);perm=NULL;
					}
				  /*
					Solution (kriging weights) are now in RHSkrig[0...nNeighborsDesired-1], which 
					points to col_weights->wts.
					I.e. they are stored.
				  */
				  for (k=0; k <nNeighborsDesired; k++) condVar -= RHSkrig[k]*Crhs[k]; // standard SK kriging variance: GSLIB II.14
				  //			  decoratorInfo->col_weights[column_number][kk_s].deviate_scale = sqrt(MAX(condVar,0));
				}
			  col_weights->deviate_scale = sqrt(MAX(condVar,0));

			  free_dmatrix(Ckrig,0,nNeighborsDesired-1,0,nNeighborsDesired-1);Ckrig=NULL;
			  // free_dvector(RHSkrig,0,nNeighborsDesired-1); these are stored
			  free_dvector(Crhs,0,nNeighborsDesired-1);Crhs=NULL;
			  free_ivector(ipiv_krig,0,nNeighborsDesired-1);ipiv_krig=NULL;

			  free(neigh);  // neighbourhood of jth gridblock
			  neigh=NULL;
			}
		  else
			{
			  // first point: no neighbours. Nothing to do.
			  // decoratorInfo->col_weights[column_number][kk_s].deviate_scale=1.0;
			  col_weights->deviate_scale=1.0;
			}
		}
	  else
		{
		  col_weights->N=0;
		}
	}
#if TEST_TREE_TIMING
  decoratorInfo->t_tree_timing += dt_timing;
  if (column_number%100==0)  printf("Time per 3D neighbour search = %g uS\n",dt_timing/(double)decoratorInfo->Nz);  
#endif

  return EXIT_SUCCESS; 
}




/**
   on calling, column_number is the increasing sequence of column labels (1,2....n),
*/
float dec3DSimulate(struct MASSAGER *massager,int layer_number, int column_number,double NGtarget)
{
  int i,j,k,k_s,kk_s,nKnown,n_try;
  int output_grid_indx,nNeighborsDesired;
  int nBest;
  int Nz_minus_1;

  struct DECORATOR_INFO *decoratorInfo; 
  struct COLUMN_TGS_KRIG_WEIGHTS *c_wts;
  struct COLUMN_TGS_KRIG_WEIGHTS *col_weights;

  double* col_sims;
  double condMean,deviate;
  double NGthresh,NG,best_NG_error;

  // safeguard inverseG:
  //NGtarget=0.8;
  //if(column_number==0)printf("FIX NGTARGET\n");
  NGthresh = inverseG(MAX(0.0001,MIN(0.9999,1-NGtarget)));
  best_NG_error = DBL_MAX; 

  i = layer_number;
  j = column_number;

  decoratorInfo = &massager->decoratorInfo;


  col_sims = dvector(0,decoratorInfo->Nz-1);

  // loop over attempts
  Nz_minus_1 = decoratorInfo->Nz-1;

  col_weights = decoratorInfo->col_weights[column_number];
  for (n_try=0; n_try<decoratorInfo->Nattempts; n_try++)
	{	  
	  // loop over nodes in column
	  k_s = j*decoratorInfo->Nz; // 3Dpath offset
	  NG=0;
	  for (kk_s=0; kk_s<decoratorInfo->Nz; k_s++,kk_s++)
		{
		  output_grid_indx = decoratorInfo->threeDpath[k_s];   // 3D index in multigrid path

		  if (decoratorInfo->flag[output_grid_indx] == 0)
			{
			  // an unknown point. Must simulate			  
			  c_wts = &col_weights[kk_s];
			  nNeighborsDesired = c_wts->N;
			  condMean = 0;
			  for (k=0; k < c_wts->N; k++)
				{
				  // residual at data point: construct from interpolation
				  // interpolate smoothedMeans to Well observation point
				  // condMean += c_wts->wts[k] * decoratorInfo->simulation[c_wts->nbr_indices[k]];
				  condMean += c_wts->wts[k] * (*(c_wts->nbr_pointers[k]));
				}
			  deviate = c_wts->deviate_scale * Gaussian_deviate(&massager->seed);
			  // Generate realisation and store
			  // for neighbours on current try....Overwrite only current column values
			  decoratorInfo->simulation[output_grid_indx] = condMean + deviate;
			  // col_sims[n_try][kk_s] = decoratorInfo->simulation[output_grid_indx];
			}
		  else
			{
			  // Hard values.
			  // col_sims[n_try][kk_s] = decoratorInfo->simulation[output_grid_indx];
			}
#if TREAT_Z_COORD_AS_CORNER_PT
		  if (kk_s < Nz_minus_1) {
			NG += (decoratorInfo->simulation[output_grid_indx] > NGthresh)? 1 : 0;
		  }
#else
		  NG += (decoratorInfo->simulation[output_grid_indx] > NGthresh)? 1 : 0;
#endif
		}

#if TREAT_Z_COORD_AS_CORNER_PT
	  NG/=Nz_minus_1;
#else
	  NG/=decoratorInfo->Nz;
#endif

	  // save best attempt
	  if (fabs(NG-NGtarget) < best_NG_error)
		{
		  nBest = n_try;
		  best_NG_error = fabs(NG-NGtarget);
		  // save best one
		  k_s = j*decoratorInfo->Nz; // path offset
		  for (kk_s=0; kk_s<decoratorInfo->Nz; k_s++,kk_s++)
			{
			  output_grid_indx = decoratorInfo->threeDpath[k_s];   // multigrid path
			  col_sims[kk_s] = decoratorInfo->simulation[output_grid_indx];
			}
		}
	}

  // printf("BEST NG %f\n",best_NG_error);

  // evaluate attempts to find best
  k_s = j*decoratorInfo->Nz; // path offset
  for (kk_s=0; kk_s<decoratorInfo->Nz; k_s++,kk_s++)
	{
	  output_grid_indx = decoratorInfo->threeDpath[k_s];   // multigrid path
	  if (decoratorInfo->flag[output_grid_indx] == 0)
		{
		  decoratorInfo->simulation[output_grid_indx] = col_sims[kk_s];
		  decoratorInfo->facies_simulation[output_grid_indx] = (decoratorInfo->simulation[output_grid_indx] > NGthresh)? 1 : 0;
		}
	}

  free_dvector(col_sims,0,decoratorInfo->Nz-1);
  col_sims=NULL;

  return best_NG_error;
}

/* 
  debugging test call
*/
int decOutputSimulation(struct MASSAGER *massager,int real_no,int dec_no)
{
  void print_postscript_xsection(struct MASSAGER *massager,int real_no,int dec_no,int show_grid_boxes);

  int i,j,j_col,k,size;
  int nCols;
  float zmin,frac;
  FILE *fp;
  char* filename;
  double** OGPoints;

  struct DECORATOR_INFO *decoratorInfo; 
  struct OUTPUTGRID *outputgrid;

  // shortcuts
  decoratorInfo = &massager->decoratorInfo;
  i = decoratorInfo->decoration_layer_indx;
  outputgrid = &massager->outputgrid[i];
  OGPoints = massager->outputgrid[i].Points; // useful shortcut

  size = outputgrid->nPoints * decoratorInfo->Nz;


#if USE_STRAT_COORDS
  double DX,DY;
  DX = outputgrid->xMax -outputgrid->xMin;
  DY = outputgrid->yMax -outputgrid->yMin;
#endif

#if TEST_TREE_TIMING
  if (decoratorInfo->t_tree_timing>0)
	{
	  // print only once
	  printf("mean time per 3D neighbour search = %g uS\n",decoratorInfo->t_tree_timing/(double)size);  
	}
  decoratorInfo->t_tree_timing=0;
#endif

  asprintf(&filename,"post_massaged_decorations_layer_%d.%d.%d.txt",
		   decoratorInfo->decoration_layer,real_no,dec_no);

  fp = fopen(filename,"w");
  if (fp == NULL)
	{ 
	  printf("Cannot open output file %s for writing\n",filename);
	  return EXIT_FAILURE;
	}
  nCols=8;
  if (!massager->duplicate) nCols+=2; // ecl indexes

  // geoEAS header
  fprintf(fp,
		  "Truncated Gaussian decoration layer %d simulation %d\n"
		  "%d\n"
		  "x\n"
		  "y\n",
		  decoratorInfo->decoration_layer,real_no,nCols);
  if (!massager->duplicate) fprintf(fp,"ecl_i\necl_j\n");
  fprintf(fp,"k\n"
		  "z_strat\n"
		  "z\n"
		  "NGtarget\n"
		  "Gauss\n"
		  "truncated_Gauss\n");
  
  for (j=0;j<size;j++)
	{
	  j_col = j/decoratorInfo->Nz;
	  k = j%decoratorInfo->Nz;
	  zmin = decoratorInfo->depth[real_no][j_col];
	  frac = decoratorInfo->Points[j][2]/decoratorInfo->z_stretch;
#if 0
	  // realisations output looks like this.....
	  // can add this in later if desired...

	  // write out x,y
	  fprintf(fp,"%g %g ",
			  outputgrid->xValues[j_col],
			  outputgrid->yValues[j_col]);
	  /* add utm stuff if need be */
	  if (massager->utmInfo.performOutputGridUTMmapping)
		{
		  fprintf(fp,"%10.2f %10.2f ",OGPoints[j_col][0], OGPoints[j_col][1]);
		}
#endif

#if USE_STRAT_COORDS
	  fprintf(fp,"%f %f ",
			  decoratorInfo->Points[j][0]*DX + outputgrid->xMin,
			  decoratorInfo->Points[j][1]*DY + outputgrid->yMin);
#else
	  fprintf(fp,"%f %f ",
			  decoratorInfo->Points[j][0],
			  decoratorInfo->Points[j][1]);
#endif
	  if (!massager->duplicate)
		{
		  /* and eclio indexes for eclipse :Use:  ecl_i, ecl_j, seg */
		  fprintf(fp,"%6d %6d ",
				  outputgrid->ecl_i[j_col], 
				  outputgrid->ecl_j[j_col]);
		}
	  fprintf(fp,"%4d %f %f %f %f %d\n",
			  k,
			  frac,
			  zmin + frac*decoratorInfo->thickness[real_no][j_col],
			  decoratorInfo->NGtarget[real_no][j_col],
			  decoratorInfo->simulation[j],
			  decoratorInfo->facies_simulation[j]);
	} 
  fclose(fp);
  printf("Wrote: %s\n",filename);
  free(filename);

  if (decoratorInfo->dumpPSdecoration)
	{
	  print_postscript_xsection(massager,real_no,dec_no,1);
	}
  return EXIT_SUCCESS;
}

void free_decorator_info(struct MASSAGER *massager)
{
  int i,kk_s,layer_number,size;
  struct DECORATOR_INFO *decoratorInfo; 
  struct OUTPUTGRID *outputgrid;

  // shortcuts
  decoratorInfo = &massager->decoratorInfo; 
  // get the right layer number 
  layer_number = decoratorInfo->decoration_layer_indx;

  outputgrid = &massager->outputgrid[layer_number];

  /* free 3D kdtree */
  killDynamicTree(decoratorInfo->dTree);		

  /* free column stuff */ 
  for (i=0; i<outputgrid->nPoints; i++)
	{
	  for (kk_s=0; kk_s<decoratorInfo->Nz; kk_s++)
		{
		  if (decoratorInfo->col_weights[i][kk_s].N > 0)
			{
			  if (decoratorInfo->col_weights[i][kk_s].wts) 
				{
				  free(decoratorInfo->col_weights[i][kk_s].wts);
				  decoratorInfo->col_weights[i][kk_s].wts = NULL;
				}
			  if (decoratorInfo->col_weights[i][kk_s].nbr_pointers)
				{
				  free(decoratorInfo->col_weights[i][kk_s].nbr_pointers);
				  decoratorInfo->col_weights[i][kk_s].nbr_pointers = NULL;
				}
			}
		}
	  free(decoratorInfo->col_weights[i]);
	  decoratorInfo->col_weights[i] = NULL;
	}
  free(decoratorInfo->col_weights);
  decoratorInfo->col_weights = NULL;

  size = outputgrid->nPoints * decoratorInfo->Nz;

  free_dvector(decoratorInfo->simulation,0,size - 1);
  decoratorInfo->simulation = NULL;
  free(decoratorInfo->facies_simulation);
  decoratorInfo->facies_simulation = NULL;

  free_dmatrix(decoratorInfo->Points,0,size- 1, 0,2);
  decoratorInfo->Points = NULL;
  free(decoratorInfo->flag);
  decoratorInfo->flag = NULL;

  free(decoratorInfo->threeDpath);
  decoratorInfo->threeDpath = NULL;

  free_matrix(decoratorInfo->NGtarget,0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  decoratorInfo->NGtarget = NULL;
  //  free_matrix(decoratorInfo->NGerror,0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  free_matrix(decoratorInfo->thickness,0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  decoratorInfo->thickness = NULL;
  free_matrix(decoratorInfo->depth,0,massager->Nrealisations-1,0,outputgrid->nPoints-1);
  decoratorInfo->depth = NULL;
}



#define GRAY_0 0.2
#define GRAY_1 0.8
#define HEIGHT_INCHES 4
#define WIDTH_INCHES 6
#define DEPTH_SENSE -1

void print_postscript_xsection(struct MASSAGER *massager,int real_no,int dec_no,int show_grid_boxes)
{
  int i,j,j_col,j_col_last,k,size,width,height,pass,iSection;
  double zmin,zminL,zminR,frac,X_g,Y_g,Color_g;
  double xL,xR,yL,yR,dYL,dYR;
  double x0,x1,x2,z0,z1,z2;
  double xBl,xBr,yBlb,yBlt,yBrb,yBrt;
  double xmin, ymin, xmax, ymax;
  double bxmin, bymin, bxmax, bymax;
  float gBfact,xscale,yscale;
  FILE *fp;
  fpos_t pos;
  char* filename;
  GArray *col_sequence,*xvals_XS,*depth_XS,*thick_XS;
  int Nz_minus_1;

  int showsAsStratigraphic=0;  // debugging; switch on to remove thickness and depth control on image
  float zero=0.0f,one=1.0f;

  struct DECORATOR_INFO *decoratorInfo; 
  struct OUTPUTGRID *outputgrid;
  struct INPUTGRID *inputgrid;

  // shortcuts
  decoratorInfo = &massager->decoratorInfo;
  i = decoratorInfo->decoration_layer_indx;
  outputgrid = &massager->outputgrid[i];
  inputgrid = &massager->inputgrid;
  size = outputgrid->nPoints * decoratorInfo->Nz;


  asprintf(&filename,"post_massaged_decorations_layer_%d.%d.%d.Xsection.eps",
		   decoratorInfo->decoration_layer,real_no,dec_no);

  fp = fopen(filename,"w");
  if (fp == NULL)
	{ 
	  printf("Cannot open output file %s for writing\n",filename);
	  return;
	}

  iSection = decoratorInfo->dumpPSdecorationVal;

  void* testArray;
  if (massager->duplicate)
	{
	  if (decoratorInfo->dumpPSdecorationDir == 'x')
		testArray = inputgrid->xValues;
	  else
		testArray = inputgrid->yValues;
	}
  else
	{
	  if (decoratorInfo->dumpPSdecorationDir == 'x')
		testArray = outputgrid->ecl_i;
	  else
		testArray = outputgrid->ecl_j;
	}
  

 
  // range pass
  ymin = xmin = DBL_MAX;
  ymax = xmax= -DBL_MAX;
  bymin = bxmin = DBL_MAX;
  bymax = bxmax= -DBL_MAX;

  float* th = decoratorInfo->thickness[real_no];
  float* dep = decoratorInfo->depth[real_no];
  col_sequence = g_array_new(FALSE, FALSE, sizeof(int));
  xvals_XS = g_array_new(FALSE, FALSE, sizeof(double));
  depth_XS = g_array_new(FALSE, FALSE, sizeof(float));
  thick_XS = g_array_new(FALSE, FALSE, sizeof(float));

  int newCol;
  int nCol=-1;
  Nz_minus_1 = decoratorInfo->Nz - 1;

  for (j = 0,j_col_last=-1;j<size;j++)
	{
	  j_col = j/decoratorInfo->Nz;
	  k = j%decoratorInfo->Nz;
	  if (massager->duplicate)
		{
		  if ((int)((double*)testArray)[j_col] != iSection) continue;
		}
	  else
		{
		  if (((int*)testArray)[j_col] != iSection) continue;
		}
	  if (j_col!=j_col_last)
		{
		  nCol++;
		  j_col_last=j_col;
		  newCol=1;
		} else newCol=0;

	  frac = decoratorInfo->Points[j][2]/decoratorInfo->z_stretch;

	  // get x,y
	  if (decoratorInfo->dumpPSdecorationDir == 'x')
		{
		 X_g = outputgrid->yValues[j_col];
		}
	  else
		{
		  X_g = outputgrid->xValues[j_col];
		}
	  if (showsAsStratigraphic)
		{
		  zmin=0;
		  Y_g = frac;
		}
	  else
		{
		  zmin = dep[j_col];
		  Y_g = (zmin + frac*th[j_col]);
		}
#if TREAT_Z_COORD_AS_CORNER_PT
	  // if (k==Nz_minus_1) continue;
#endif

	  // scale to [0,1]  
	  xmin = MIN(xmin,X_g);
	  xmax = MAX(xmax,X_g);
	  ymin = MIN(ymin,Y_g);
	  ymax = MAX(ymax,Y_g);
	  if (newCol)
		{
		  // save thickness, x's, depths
		  g_array_append_val(col_sequence,j_col);
		  g_array_append_val(xvals_XS,X_g);
 
		  if (showsAsStratigraphic)
			{
			  g_array_append_val(depth_XS,zero);
			  g_array_append_val(thick_XS,one);
			}
		  else
			{
			  g_array_append_val(depth_XS, dep[j_col]);
			  g_array_append_val(thick_XS, th[j_col]);
			}

		  j_col_last=j_col;
		  nCol++;
		}
	}
  int xSectionNcols = xvals_XS->len;

  // drawing pass
  width=WIDTH_INCHES*72;
  height=HEIGHT_INCHES*72;

   /* PostScript header */
  fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp,"%%%%Creator: decorator.c (James Gunning)\n");
  fgetpos(fp,&pos);
  fprintf(fp,"%%%%BoundingBox: %4d %4d %4d %4d\n",0,0,width,height);
  fprintf(fp,"                           \n"); // protection against fsetpos hack below
  fprintf(fp,"%%%%EndComments\n");
  // fprintf(fp,"%d %d translate\n", 72, 72 );
  fprintf(fp,".05 setlinewidth\n");

  /* The +72 shifts the figure one inch from the lower left corner */

  gBfact = (decoratorInfo->Nz>1)? 1.0*(decoratorInfo->Nz-1) : 1;
  xscale=(xmax>xmin) ? width/(xmax-xmin): 1;
  yscale = DEPTH_SENSE* ((ymax>ymin) ? height/(ymax-ymin): 1);
  for (pass=0;pass<=(show_grid_boxes?1:0);pass++)
	{
	  if (pass==1) 
		{
		  fprintf(fp,"\n\n%%%%Grid boxes from here till Showpage\n\n");
		  fprintf(fp,"0.0 setgray\n");
		}
	  for (j = 0,nCol=-1,j_col_last = -1;j<size;j++)
		{
		  j_col = j/decoratorInfo->Nz;
		  k = j%decoratorInfo->Nz;

		  // pick out Xsection by skipping points off-xection
		  if (massager->duplicate)
			{
			  if ((int)((double*)testArray)[j_col] != iSection) continue;
			}
		  else
			{
			  if (((int*)testArray)[j_col] != iSection) continue;
			}
		  if (j_col!=j_col_last)
			{
			  nCol++;
			  j_col_last=j_col;
			}

		  frac = decoratorInfo->Points[j][2]/decoratorInfo->z_stretch;

		  Color_g =  decoratorInfo->facies_simulation[j]; 
		  // Color_g =  0.5+decoratorInfo->simulation[j]/6.0; 

		  if (nCol>0)
			{
			  x1 = g_array_index(xvals_XS,double,nCol);
			  x0 = g_array_index(xvals_XS,double,nCol-1);
			  xL = (x1+x0)/2.0;
			  z1 = g_array_index(depth_XS,float,nCol);
			  z0 = g_array_index(depth_XS,float,nCol-1);
			  zminL = (x1!=x0) ? 0.5*(z1 + z0) : z1;
			  yL = zminL + frac*(g_array_index(thick_XS,float,nCol) +g_array_index(thick_XS,float,nCol-1))/2.0;
			  dYL = (g_array_index(thick_XS,float,nCol)+g_array_index(thick_XS,float,nCol-1))/(2.0*gBfact);
			}
		  else
			{
			  zminL = g_array_index(depth_XS,float,nCol);
			  xL = g_array_index(xvals_XS,double,nCol);
			  yL = zminL + frac*g_array_index(thick_XS,float,nCol);
			  dYL = g_array_index(thick_XS,float,nCol)/gBfact;
			}
		  if (nCol<(xSectionNcols-1))
			{
			  x1 = g_array_index(xvals_XS,double,nCol);
			  x2 = g_array_index(xvals_XS,double,nCol+1);
			  xR = (x1+x2)/2.0;
			  z1 = g_array_index(depth_XS,float,nCol);
			  z2 = g_array_index(depth_XS,float,nCol+1);
			  zminR = (x2!=x1) ? 0.5*(z1+z2) : z1;

			  yR = zminR + frac*(g_array_index(thick_XS,float,nCol) +g_array_index(thick_XS,float,nCol+1))/2.0;
			  dYR = (g_array_index(thick_XS,float,nCol) +g_array_index(thick_XS,float,nCol+1))/(2.0*gBfact);
			}
		  else
			{
			  zminR = g_array_index(depth_XS,float,nCol);
			  xR = g_array_index(xvals_XS,double,nCol);
			  yR = zminR + frac*g_array_index(thick_XS,float,nCol);
			  dYR = g_array_index(thick_XS,float,nCol)/gBfact;
			}

		  // rescale
		  xBl = (xL-xmin)*xscale;
		  xBr = (xR-xmin)*xscale;
		  yBlb = (yL-dYL-ymin)*yscale;
		  yBlt = (yL+dYL-ymin)*yscale;
		  yBrt = (yR+dYR-ymin)*yscale;
		  yBrb = (yR-dYR-ymin)*yscale;

          // don't draw last box if working in corner pt geom
#if TREAT_Z_COORD_AS_CORNER_PT
		  // if (k==Nz_minus_1) continue;
#endif
		  // draw
		  fprintf(fp,"newpath\n");
		  fprintf(fp,"%f %f moveto\n",xBl,yBlb);
		  fprintf(fp,"%f %f lineto\n",xBl,yBlt);
		  fprintf(fp,"%f %f lineto\n",xBr,yBrt);
		  fprintf(fp,"%f %f lineto\n",xBr,yBrb);
		  fprintf(fp,"closepath\n");

		  if (pass==0) {
			fprintf(fp,"%f setgray\nfill\n",(Color_g==1)?GRAY_1:GRAY_0);
		  } else {
			fprintf(fp,"stroke\n");
		  }
		  if (pass==0)
			{
			  // keep track of bounding box...
			  bxmin = MIN(bxmin,xBl);
			  bxmax = MAX(bxmax,xBr);
			  bymin = MIN(MIN(MIN(MIN(bymin,yBrb),yBlb),yBrt),yBlt);
			  bymax = MAX(MAX(MAX(MAX(bymax,yBrb),yBlb),yBrt),yBlt);
			}
		}
	}
  fprintf(fp,"showpage\n%%%%EOF\n");
  
  // rewrite bounding box
  fsetpos(fp,&pos);
  fprintf(fp,"%%%%BoundingBox: %4d %4d %4d %4d\n",(int)floor(bxmin-1),(int)floor(bymin-1),
		  (int)ceil(bxmax+1),(int)ceil(bymax+1));
  fclose(fp);
  printf("Wrote: %s\n",filename);
  free(filename);
  g_array_free (col_sequence,TRUE);
  g_array_free (xvals_XS,TRUE);
  g_array_free (depth_XS,TRUE);
  g_array_free (thick_XS,TRUE);
}

#undef GRAY_0
#undef GRAY_1
#undef WIDTH_INCHES
#undef HEIGHT_INCHES
#undef DEPTH_SENSE 
#undef TREAT_Z_COORD_AS_CORNER_PT
