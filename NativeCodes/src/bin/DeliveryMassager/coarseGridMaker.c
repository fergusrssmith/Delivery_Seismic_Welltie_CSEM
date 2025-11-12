#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>

#include "ECLIO/eclio.h"  /* structure definition and a little miscellany */
#include "ECLIO/access.h" /* access high-level routines */
#include "ECLIO/eclutil.h" /* for model initialization only...my move to access */
#include "ECLIO/blockops.h" /* also added some lower-level routines, flagged below */
#include "ECLIO/copy.h"    /* ways to set records for moving data between models */


#include "errorHandler.h"

#include "NR/nrutil.h"


#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"

#include "cgeom/vdefs.h"

#include "platform_specifics.h"
#include "massager_types.h"
#include "clapack_cblas_wrappers.h"
#include "massager.h"

#include "eclipseHandler.h" 
#include "coarseGridMaker.h"

#define SHELL_SEARCH 1

typedef struct GLBS_INFO
{ 
  int ks_delta; 
  int* k_search;
  GArray* ws;
  float n_eff_depths;
  // eclipse stuff
  char* filename;
  ECL ecl;
  double ****zCorn;
  int nKeys_ecl;
  char** keyWords_ecl; 

  // other useful info
#if SHELL_SEARCH
  int n_minima;
  int shell_search_plateau_mode;
  GHashTable* plateau_hashtable;
  GPtrArray* plateau_index_list;
#endif

  double** z_intersections;
  double s_best;
  int* k_best;
} GLBS_INFO;

typedef struct DEPTH_INFO 
{
  // observations
  int i;
  double x,y,z,weight;
  float guess_weight;
  // ecl info
  int i_ecl,j_ecl,k_ecl,vrtx_ecl;
  int kmin_ecl,kmax_ecl;
  double *coord_ecl;
  // other useful stuff
  double dmin;
} DEPTH_INFO;

//THREAD GLBS_INFO info;
GLBS_INFO info;



int makeOptimalCoarseLayeredECL(struct MASSAGER *massager)
{
  int i,i_g,j,k,kmin,kmax,iB,jB,kB,Nx,Ny,Nz;
  int *kFine,*kFineInit;
  float *kFineNorm,*kFineAvg;
  
  printf("Making optimal subgrid...\n");
  // get refined depth array from XML
  if (massager->d_index<0)
	{
	  fprintf(stderr,"No depth observations available\n");
	}
  getDepthArray(massager);
  if ((info.ws)->len==0)
	{
	  fprintf(stderr,"No depth data available in xml to subdivide Eclipse grid\n");
	  return EXIT_FAILURE;
	}

  // Read eclipse file
  if (readEclipseGridToSubgrid(massager)==EXIT_FAILURE)
	{
	  fprintf(stderr,"Failed to read Eclipse grid\n");
	  return EXIT_FAILURE;
	}

  /* Eclipse dimensions */
  Nx = getNX(&(info.ecl));
  Ny = getNY(&(info.ecl));
  Nz = getNZ(&(info.ecl));  

  // loop over ECl to get nearest points to well depth info, 
  // examining top 4 corners of ECL blocks
  int i_c,j_c,vrtx,found;
  double x,y,z,xMin,xMax,yMin,yMax,ds;
  double pt[4];
  DEPTH_INFO *d;

  for (i=0;i<(info.ws)->len;i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  d->dmin = DBL_MAX;
	}
  printf("Start column search\n");

  // exhaustive search for initial guesses
  // lookup table
  info.z_intersections = dmatrix(1,Nz,0,(info.ws)->len - 1);
  memset(&info.z_intersections[1][0],0,sizeof(double)*Nz*(info.ws)->len);

  double ***c,*zc;
  c = info.ecl.coord;
  found=0;
  for (iB=1; iB<=Nx; iB++)
	{
	  for (jB=1; jB<=Ny; jB++) 
		{		 		  
		  for (kB=1; kB<=Nz ; kB++)
			{
			  // use 4 bottom corners for bilinear....
			  double pt_cnrs[12];
			  double ptr[4];
			  int c_pt;
			  zc = (info.zCorn)[iB][jB][kB];
			  xMin=yMin=DBL_MAX;
			  xMax=yMax=-DBL_MAX;
			  for (vrtx=1,c_pt=0; vrtx<=4; vrtx++,c_pt+=3)
				{ 
				  i_c = (iB-1) + ((vrtx&1) ? 0 : 1);
				  j_c = (jB-1) + ((vrtx==3 || vrtx==4 || vrtx==7 || vrtx==8) ? 1: 0);
				  interpVert(zc[vrtx],ptr,c[i_c][j_c]);
				  pt_cnrs[c_pt] = ptr[1];
				  pt_cnrs[c_pt+1] = ptr[2];
				  pt_cnrs[c_pt+2] = ptr[3];
				  // get rough bounds
				  x = ptr[1];
				  y = ptr[2];
				  xMin=MIN(xMin,x);
				  xMax=MAX(xMax,x);
				  yMin=MIN(yMin,y);
				  yMax=MAX(yMax,y);
				}
			  // now {x1,y1,z1, x2,y2,z2, ...} are in pt_cnrs[0...11]
			  // check all candidate data
			  for (i=0;i<(info.ws)->len;i++)
				{
				  double testPt[3];
				  double u,v;
				  d = g_array_index((info.ws),DEPTH_INFO*,i);
				  if ((d->x < xMin) || (d->x > xMax) || (d->y < yMin) || (d->y>yMax))
					{
					  // can't possibly be here
					  continue;
					}
				  testPt[0] = d->x;
				  testPt[1] = d->y;
				  // see if it's OK
				  // debug_plotTestPt(pt_cnrs,testPt);
				  doBilinearInterp(pt_cnrs,testPt,&u,&v);
				  if (u<0 || u>1 || v<0 || v>1) continue; // data pt doesn't fit in containing rectangle
				  // If here, point is "transversely" locatable
				  // and surface height z = testPt[2]
				  z = testPt[2];
				  // count it first time
				  if (d->dmin==DBL_MAX) found++;				  


				  // save useful loookup
				  info.z_intersections[kB][i] = z;

				  ds = fabs(z - d->z); 
				  if (ds < d->dmin) 
					{
					  d->i_ecl = iB;
					  d->j_ecl = jB;
					  d->k_ecl = kB;
					  d->coord_ecl = c[i_c][j_c];
					  d->vrtx_ecl = vrtx;
					  d->dmin = ds;					 
					}
				}
			}
		}
	}

  if (found==0)
	{
	  // no data in eclipse grid. Give up
	  fprintf(stderr,"No depth data in xml locatable in Eclipse grid. Giving up. Check units?\n");
	  return EXIT_FAILURE;
	}
  else if (found != info.ws->len)
	{
	  fprintf(stderr,"Warning! Only %d of %d depth data in xml locatable in Eclipse grid. Check units?\n",found,info.ws->len);
	  // return EXIT_FAILURE; 
	}

  // second pass reweighting to soften influence for matches in multi-pinchout regions
  for (i=0;i<(info.ws)->len;i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if (d->dmin==DBL_MAX) continue; // data not inside grid, skip

	  iB = d->i_ecl;
	  jB = d->j_ecl;
	  kB = d->k_ecl;
	  ds = d->dmin;
	  vrtx = d->vrtx_ecl;
	  j=1;
	  d->kmin_ecl = d->kmax_ecl = kB;
	  // search up
	  for (k=kB+1;k<=Nz;k++)
		{
		  pt[3] = info.z_intersections[k][i]; 
		  ds = fabs(pt[3] - d->z); 
		  if (ds <= d->dmin) d->kmax_ecl++;
		  else break;
		}
	  // search down
	   for (k=kB-1;k>=1;k--)
		{
		  pt[3] = info.z_intersections[k][i]; 
		  ds = fabs(pt[3] - d->z); 
		  if (ds <= d->dmin) d->kmin_ecl--;
		  else break;
		}
	   // now, all [d->kmin_ecl...d->kmax_ecl] points are equally likely
	   // pick midpoint
	   d->k_ecl = (int)floor((d->kmax_ecl+d->kmin_ecl)/2.0 + 0.5);
	   // and downweight
	   d->guess_weight = 1.0/(d->kmax_ecl - d->kmin_ecl+1.0);
	}
  printf("End column search\n");

  // Form initial estimates of fine grid indices, to seed global search
  kFine = ivector(0,massager->nLayers);
  kFineInit = ivector(0,massager->nLayers);
  kFineNorm = vector(0,massager->nLayers);
  kFineAvg = vector(0,massager->nLayers);
  for (i=0; i<=massager->nLayers; i++)
	{
	  kFine[i]=0;
	  kFineNorm[i]=0.0;
	  kFineAvg[i]=0.0;
	}
  for (i=0; i<(info.ws)->len; i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if (d->dmin==DBL_MAX)
		{
		  // data not inside grid, skip
		  continue;
		}
	  iB = d->i_ecl;
	  jB = d->j_ecl;
	  kB = d->k_ecl;

	  z = info.z_intersections[kB][i]; 
	  printf("Nearest points: well-obs=(%g %g %g,  coarse-layer %d): grid (z=%g, fine-layer k:[%d-%d]), |dist| = %g\n",
			 d->x,d->y,d->z,d->i,z,d->kmin_ecl,d->kmax_ecl,d->dmin);
	  kFineNorm[d->i] += d->guess_weight;
	  kFineAvg[d->i] += d->guess_weight * kB;
	}
  for (i=0; i<=massager->nLayers; i++)
	{
	  if (kFineNorm[i]>0) 
		{ 
		  kFine[i] = floor(kFineAvg[i]/kFineNorm[i]+0.5);
		}
	  kFineInit[i] = kFine[i];
	  printf("Estimated index, massager layer %d, is fine-index %d\n",i, (int)kFine[i]);
	}

  // Recursive search 
  /* max and min matching fine layer indx vectors for coarse laer j */
  int *kspan_min,*kspan_max;
  kspan_min = ivector(0,massager->nLayers);
  kspan_max = ivector(0,massager->nLayers);

  // initialise, then lookup over subgroups
  for (i=0,j=0;i<(info.ws)->len;i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if (d->dmin==DBL_MAX) continue; // data not inside grid, skip
	  i_g = d->i;
	  kspan_max[j] = d->kmax_ecl;
	  kspan_min[j] = d->kmin_ecl;
	  // group starts at i
	  for (k=i+1; k<(info.ws)->len; k++)
		{
		  d = g_array_index((info.ws),DEPTH_INFO*,k);		  
		  if (d->dmin==DBL_MAX) continue; // data not inside grid, skip
		  if (d->i == i_g)
			{ 
			  kspan_max[j] = IMAX(d->kmax_ecl,kspan_max[j]);
			  kspan_min[j] = IMIN(d->kmin_ecl,kspan_min[j]);
			}
		  else
			{
			  j++; // start next group
			  i=k-1; // starting one back...
			  break;
			}
		}
	}
  for (i=0; i<=massager->nLayers; i++)
	{
	  printf("Coarse layer indx=%d:, likely optimum subgrid layer range %d,%d\n",i, kspan_min[i],kspan_max[i]);
	}

  // IDEA: could use ranges should above kspan_min[i],kspan_max[i]), to set hypercube search span
  // in below (not implemented yet)

  // we can have in principle massager->nLayers+1 useful depth controls: all tops + a bottom
  (info.k_search) = ivector(0,massager->nLayers);
  (info.k_best) = ivector(0,massager->nLayers);
  (info.s_best) = DBL_MAX;
  for (i=0; i<=massager->nLayers; i++)  info.k_best[i]=-1;
  int ks_indx=-1;

#if SHELL_SEARCH
  // greedy neighbour search
  int hunt,*list;
  double s_min=DBL_MAX;

  // shell searches are confined to adjacent lattice sites plus center site:
  info.ks_delta = 1; // searches in nearest shell...
  info.plateau_index_list = g_ptr_array_new();
  info.plateau_hashtable = g_hash_table_new_full(g_str_hash,g_str_equal,string_hash_free,NULL);
  info.shell_search_plateau_mode = 0;

  printf("Shell search:\n");
  hunt=1;
  while (hunt)
	{ 
	  if (massager->verbosity>4 && s_min<DBL_MAX) 
		{	  
		  printf("Hypercube-steepest-descent search solution: ");
		  for (i=0; i<=massager->nLayers; i++)  printf("%4d ",info.k_best[i]);
		  printf("\n");
		}
	  // list length can vary in plateau mode, so make it recheckable
	  for (k=0; k< ((info.shell_search_plateau_mode==1)? info.plateau_index_list->len : 1); k++)
		{
		  // do onion shell search about kth location in plateau list
		  if (info.shell_search_plateau_mode==1)
			{
			 list = (int*) g_ptr_array_index(info.plateau_index_list,k);
			 for (i=0; i<=massager->nLayers; i++)  
			   {
				 kFine[i] = list[i]; 
			   }
			}
		  // shell neighbour search
		  info.n_minima=0;
		  hypercube_loop(&ks_indx,Nz,kFine,massager);
		  if (info.shell_search_plateau_mode==0)
			{
			  // last attempt was descent mode, 			  
			  // but was it unique?
			  if (info.n_minima>1)
				{
				  // Not unique. Search the plateau..
				  info.shell_search_plateau_mode = 1;
				  if (info.s_best < s_min) s_min = info.s_best;  // remember best match so far
				  break;
				}
			  else
				{
				  // unique minima. Exit if minima is center of shell
				  // assume final solution (i.e. at shell-center..)
				  hunt=0;
				  for (i=0; i<=massager->nLayers; i++)  
					{
					  // till shown otherwise...
					  if (kFine[i] != info.k_best[i]) {hunt=1;break;}
					}
				}
			  if (info.s_best < s_min)
				{
				  s_min = info.s_best;  // remember best match so far
				  for (i=0; i<=massager->nLayers; i++)  kFine[i] = info.k_best[i]; // steepest descent step
				}
			  if (!hunt) break;
			}
		  else
			{
			  // still in plateau mode. No unique solution if on last list member . Bail out
			  if (k == (info.plateau_index_list->len - 1))
				{
				  hunt=0;
				  break;
				}
			}
		}
	}
  // cleanup
  flush_plateau_list();
  g_ptr_array_free(info.plateau_index_list,TRUE);
  g_hash_table_destroy(info.plateau_hashtable);

#else

  // Exhaustive search +/- 80 layers of initial guess
  // (can be very expensive)
  (info.ks_delta) = 80; 

  // exhaustive search
  hypercube_loop(&ks_indx,Nz,kFine,massager);

#endif

  int isOK=1;
  for (i=0; i<massager->nLayers; i++)
	{
	  printf("Optimum subdivision index, massager layer %d, is fine-index %d\n",i, info.k_best[i]);
	  if (info.k_best[i]==-1) isOK=0;
	}
  printf("Optimum subdivision index, base-massager layer %d, is fine-index %d\n",
		 massager->nLayers-1, info.k_best[massager->nLayers]);

  if (!isOK) 
	{
	  printf("Cannot identify depths of all massager layers. Coarse grid creation not possible.\n"
			 "Add more well observation depth info to XML to help. Giving up\n");
	  return EXIT_FAILURE;
	}
  printf("Optimisation minima = %g\n", info.s_best); 

#if !SHELL_SEARCH
  // check if search terminated on initial subbox wall...
  isOK=1;
  for (i=0; i<=massager->nLayers; i++)
	{
	  if (((int)fabs(info.k_best[i] - kFineInit[i]) ==  info.ks_delta) && 
		  (info.k_best[i]!=1) && (info.k_best[i]!=Nz)) 
		{
		  isOK=0;
		}
	}
  if (!isOK) printf("WARNING: optimal subgridding solution limited by initially guessed index bounds!\n");
#endif



  printf("Optimum fit vertical mismatches:\n");
  printf("%12s %12s %12s %24s %10s %12s\n","x_obs","y_obs","z_obs","coarse-layer-index_obs","grid_z","z_obs-z_surf");
  for (i=0; i<(info.ws)->len; i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if (d->dmin==DBL_MAX)
		{
		  // data not inside grid, skip
		  continue;
		}
	  iB = d->i_ecl;
	  jB = d->j_ecl;
	  kB = info.k_best[d->i];
	  z = info.z_intersections[kB][i];
 
	  printf("%12.6g %12.6g %12.6g %24d %10.6g %12.2g\n",
			 d->x,d->y,d->z,d->i, z, d->z - z);
	}



  // Now make smaller eclipse using surfaces at kBest :-)
  if (exportCoarseGrid(massager)==EXIT_FAILURE)
	{
	  fprintf(stderr,"Failed to write Eclipse grid\n");
	  return EXIT_FAILURE;
	}






  // clean up
  /* input eclipse grid */
  freeModel(&info.ecl);

 
  
  /* clean ups */ 
  free_dmatrix(info.z_intersections,1,Nz,0,(info.ws)->len-1);
  free_ivector(kFine,0,massager->nLayers); kFine=NULL;
  free_ivector(kFineInit,0,massager->nLayers); kFineInit=NULL;
  free_vector(kFineNorm,0,massager->nLayers); kFineNorm=NULL;
  free_vector(kFineAvg,0,massager->nLayers); kFineAvg = NULL;
  
  free_ivector(info.k_search, 0, massager->nLayers);
  free_ivector(info.k_best, 0, massager->nLayers);

  free_ivector(kspan_min,0,massager->nLayers);
  free_ivector(kspan_max,0,massager->nLayers);

  for (i=0;i<(info.ws)->len;i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if(d) {free(d); d = NULL;}
	}

  g_array_free((info.ws), TRUE);
  (info.ws)=NULL;

  return EXIT_SUCCESS;
}


double hypercube_loop(int *ks_indx,int Nz, int *kFine,struct MASSAGER *massager)
{
  int k,kmin=1,kmax;
  double s=0;

  (*ks_indx)++;

  if ((*ks_indx) == (massager->nLayers+1))
	{
	  s = vertical_mismatch();

#if SHELL_SEARCH
	 if (info.shell_search_plateau_mode==1)
	   {
		 // plateau mode;
		 if (s < (info.s_best))
		   {
			 // possible optimum unvisited. try this
			 for (k=0;k<=massager->nLayers; k++) (info.k_best)[k] = (info.k_search)[k]; 
			 info.s_best = s;			 
			 info.n_minima=1;
			 flush_plateau_list();
			 append_to_plateau_list(massager->nLayers);
			 // flip to descent mode
			 info.shell_search_plateau_mode = 0;
		   }
		 else if (s == (info.s_best))
		   {
			 // add new one
			 append_to_plateau_list(massager->nLayers);
			 info.n_minima++; 
		   }
	   }
	 else
	   {
		 // descent mode, or subsequent nodes after breaking from plateau mode..
		 if (s < (info.s_best))
		   {
			 // possible optimum unvisited. try this
			 for (k=0;k<=massager->nLayers; k++) (info.k_best)[k] = (info.k_search)[k]; 
			 info.n_minima=1;
			 info.s_best = s;
			 flush_plateau_list();
			 // add new one
			 append_to_plateau_list(massager->nLayers);
		   }
		 else if (s == info.s_best)
		   {
			 info.n_minima++;
			 append_to_plateau_list(massager->nLayers);
		   }
	   }
#else
	  if (s < (info.s_best))
		{	  
		  // assign	only for improvements			  
		  for (k=0;k<=massager->nLayers; k++) (info.k_best)[k] = (info.k_search)[k]; 
		  info.s_best = s;
		  if (massager->verbosity>4) 
			{
			  printf("Hypercube-search best solution: ");
			  for (k=0;k<=massager->nLayers; k++) printf("%5d ",(info.k_search)[k]);
			  printf(" mean-mismatch=%g\n",sqrt(s)/info.n_eff_depths); 
			}
		}
#endif
	}
  else
	{
	  // search must start at k>=1, k > k_previous_index but also near optimum.
	  if (*ks_indx > 0) 
		{
		  // to insist on separate surfaces...
		  // kmin = IMAX(1,(info.k_search)[*ks_indx-1]+1);
		  // to allow degeneracy...
		  kmin = IMAX(1,(info.k_search)[*ks_indx-1]);
		}
	  kmin = IMAX(kmin, kFine[*ks_indx]-(info.ks_delta));
	  kmax = IMIN(Nz, kFine[*ks_indx]+(info.ks_delta));
	  for ((info.k_search)[*ks_indx] = kmin;
		   (info.k_search)[*ks_indx] <= kmax;
		   (info.k_search)[*ks_indx]++)
		{
		  s = hypercube_loop(ks_indx,Nz,kFine,massager);
		}
	}
  (*ks_indx)--;
  return s; 
}

double vertical_mismatch()
{
  int i,i_c,j_c,kB,vrtx;
  double ds,dsq;
  double z;
  DEPTH_INFO* d; 
 
  dsq=0.0;
  for (i=0; i<(info.ws)->len; i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  if (d->dmin==DBL_MAX)
		{
		  // data not inside grid, skip
		  continue;
		}
	  kB = (info.k_search)[d->i];
	  z = info.z_intersections[kB][i]; 

	  // weight mismatch by number of duplicates
	  ds = d->weight * fabs(z - d->z);

	  dsq += ds;
	}
  return dsq;
}



void getDepthArray(struct MASSAGER *massager)
{
  int i,j,k,ii,jj,n_dup,depthObsCount; 
  int canAdd;
  double x,y,z,z_orig,z_acc;
  DEPTH_INFO *d,*ds;
  struct WELL_OBSERVATIONS *wo,*woSearch;

  
  (info.ws) = g_array_new(FALSE, FALSE, sizeof(DEPTH_INFO*));


  // (1) assemble list ofabsolute depths from well info


  for (i=0; i<massager->nLayers; i++)
	{
	  wo = &massager->well_observations[i];
	  depthObsCount=0; 

  
	  // loop over observations
	  for (j=0; j<wo->n; j++)
		{
		  // is this a depth observation?
		  if ((massager->d_index>=0 && !strcmp(g_array_index(wo->property_name, gchar*,j), 
											   massager->propertyNames[massager->d_index])))
			{
			  // Yes. here are the x,y's
			  d = (DEPTH_INFO*)malloc(sizeof(DEPTH_INFO));
			  d->i = i;
			  d->x = g_array_index (wo->xValues, gdouble, j);
			  d->y = g_array_index (wo->yValues, gdouble, j);
			  d->z = g_array_index(wo->values, gdouble,j);
			  g_array_append_val((info.ws),d);
			  z_orig = z_acc = d->z;

			  if (massager->verbosity>4) printf("Depths obs: layer %d x=%g y=%g z=%g\n",d->i,d->x,d->y,d->z);
			  depthObsCount++;

			  // search for any thickness in other layers with same x,y that
			  // are addable from this ref.
			  if (massager->thickness_index < 0) continue; // no thickness available

			  // (1) go down
			  for (ii=i; ii<massager->nLayers; ii++)
				{
				  woSearch = &massager->well_observations[ii];
				  canAdd=0;
				  // search other observation for this layer
				  for (jj=0; jj<woSearch->n; jj++)
					{
					  // is this a thickness observation?
					  if (!strcmp(g_array_index(woSearch->property_name, gchar*,jj), 
								  massager->propertyNames[massager->thickness_index]))
						{
						  // Yes. here are the x,y's
						  x = g_array_index (woSearch->xValues, gdouble, jj);
						  y = g_array_index (woSearch->yValues, gdouble, jj);
						  if (APPROX(d->x,x,WELL_XY_MATCH_EPS) && APPROX(d->y,y,WELL_XY_MATCH_EPS))
							{
							  // its in the same location. Add to database
							  canAdd=1;
							  ds = (DEPTH_INFO*)malloc(sizeof(DEPTH_INFO));
							  ds->i = ii+1;
							  ds->x = d->x;
							  ds->y = d->y;
							  // add thickness
							  z_acc += g_array_index(woSearch->values, gdouble,jj);
							  ds->z = z_acc;
							  g_array_append_val((info.ws),ds);
							  if (massager->verbosity>4) printf("Depths obs, downward search: %d %g %g %g\n",ds->i,ds->x,ds->y,ds->z);
							  depthObsCount++;
							}
						}						  
					 }
				  if (!canAdd) break;
				}
			  // (2) go up
			  z_acc = z_orig; // restore starting depth
			  for (ii=i-1; ii>=0; ii--)
				{
				  woSearch = &massager->well_observations[ii];
				  canAdd=0;
				  // search other observation for this layer
				  for (jj=0; jj<woSearch->n; jj++)
					{
					  // is this a thickness observation?
					  if (!strcmp(g_array_index(woSearch->property_name, gchar*,jj), 
								  massager->propertyNames[massager->thickness_index]))
						{
						  // Yes. here are the x,y's
						  x = g_array_index (woSearch->xValues, gdouble, jj);
						  y = g_array_index (woSearch->yValues, gdouble, jj);
						  if (APPROX(d->x,x,WELL_XY_MATCH_EPS) && APPROX(d->y,y,WELL_XY_MATCH_EPS))
							{
							  // its in the same location. Add to database
							  canAdd=1;
							  ds = (DEPTH_INFO*)malloc(sizeof(DEPTH_INFO));
							  ds->i = ii;
							  ds->x = d->x;
							  ds->y = d->y;
							  // subtract thickness
							  z_acc -= g_array_index(woSearch->values, gdouble,jj);
							  ds->z = z_acc;
							  g_array_append_val((info.ws),ds);
							  if (massager->verbosity>4) printf("Depths obs, upward search : %d %g %g %g\n",ds->i,ds->x,ds->y,ds->z);
							  depthObsCount++;
							}
						}						  
					 }
				  if (!canAdd) break;
				}
			}
		}
	}

  // 2 deal with possible duplicates....
  // first sort, by x,y,z  
  g_array_sort(info.ws, (GCompareFunc)depth_info_compareLXYZ);

  // then de-weight duplicates
  // initialise locatability while we're here
  for (i=0; i<(info.ws)->len; i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  x = d->x;
	  y = d->y;
	  ii = d->i;
	  d->weight = 1.0;
	  for (j=i+1,n_dup=1; j<(info.ws)->len; j++)
		{
		  ds = g_array_index((info.ws),DEPTH_INFO*,j);
		  if ((x==ds->x) && (y==ds->y) && (ii == ds->i))
			{
			  // a duplicate
			  n_dup++;
			}
		  else break;
		}
	  if (n_dup>1)
		{
		  // de-weight the data in range [i,j-1]
		  for (k=i;k<=(j-1); k++)
			{
			  ds = g_array_index((info.ws),DEPTH_INFO*,k);
			  ds->weight=1.0/n_dup;
			}
		  // and proceed from j (expect ++)
		  i=j-1;
		}
	}
  info.n_eff_depths=0.0;
  for (i=0; i<(info.ws)->len; i++)
	{
	  d = g_array_index((info.ws),DEPTH_INFO*,i);
	  info.n_eff_depths += d->weight;
	  printf("Reweighted depths obs, %d %12.7g %12.7g %12.7g wt=%g\n",d->i,d->x,d->y,d->z,d->weight);
	}
  printf("Effective number of depth observations = %d\n",(int)info.n_eff_depths);
}


/*
  Sort compare to sort coarse layer index, then x,y, depth measurement z
*/
gint depth_info_compareLXYZ(gconstpointer a, gconstpointer b)
{
  DEPTH_INFO *A,*B;

  A = *((DEPTH_INFO **)a);
  B = *((DEPTH_INFO **)b);

  if (A->i > B->i) return 1;
  else if (A->i < B->i) return -1;
  else
	{
	  if (A->x > B->x) return 1;
	  else if (A->x < B->x) return -1;
	  else
		{
		  if (A->y > B->y) return 1;
		  else if (A->y < B->y) return -1;
		  else 
			{
			  if (A->z > B->z) return 1;
			  else if (A->z < B->z) return -1;
			  else return 0; // don't care
			}
		}
	}

}

int readEclipseGridToSubgrid(struct MASSAGER *massager)
{
  int i,j,k;


  int haveSegmentsKeyword = 0;
  int haveLayersKeyword = 0;
  int status,transformToMap;
  struct OUTPUTGRID *mog_i;
  char msg[LINE_LENGTH];
  PROP *from=NULL;
  char* filename;

  /* read Eclipse grids */
  info.keyWords_ecl = NULL; 
  info.ecl.coord = (double ***) NULL;  // coord lines

  transformToMap = 1; /* number of records to read and
						 force transform grid->map coords */
  mog_i = &(massager->outputgrid[0]);

  filename = (massager->outputgrid[0]).filename;

  if (!filename)
	{
	  fprintf(stderr,"No Eclipse grid specified, or duplicate (-d) mode used\n");
	  return EXIT_FAILURE;
	}


  if (massager->verbosity>1) printf("...Loading Eclipse model from file %s\n", filename);
  info.nKeys_ecl = fetchKeyList(filename, &info.keyWords_ecl, "--");

  if (info.nKeys_ecl <= 0) return EXIT_FAILURE;


  /* List possible keywords */
  if (massager->verbosity>2)
	{
	  fprintf(stdout, 
			  "\nFound %d *probable* keys in file %s\n"
			  "Key Number   Key Value                        Message\n"
			  "------------ -------------------------------- ------- \n",
			  info.nKeys_ecl,filename);
	}
  for (j=1;j<=info.nKeys_ecl;j++) 
	{
	  checkKeyWord(info.keyWords_ecl[j],"(){}[]-_",msg,LINE_LENGTH);
	  if (massager->verbosity>2) fprintf(stdout,"%32s %s\n",info.keyWords_ecl[j], msg);
	}
  // reduce list to useful ones... 
  for (j=1,k=1;j<=info.nKeys_ecl;j++) 
	{ 
	  if (isADesiredECLIPSEKeyword(info.keyWords_ecl[j])) 
		{
		  info.keyWords_ecl[k++]=info.keyWords_ecl[j];
		}
	}
  info.nKeys_ecl = k-1;
					  
  /* After a bit of cleaning -- if you get a message here, program will probably crash 
	 when loading ... will provide a cushier exit later*/
  info.nKeys_ecl = editKeyList(info.nKeys_ecl, &info.keyWords_ecl);
  if (massager->verbosity>2)
	{
	  printf("\nAfter editKeyList: %d *probable* keys in file %s\n"
			 "Key Number   Key Value                        Message\n"
			 "------------ -------------------------------- ------- \n",
			 info.nKeys_ecl,filename);
	}
  /* Check we have segments & layers too */ 
  haveSegmentsKeyword=0;
  haveLayersKeyword=0;
  for (j=1;j<=info.nKeys_ecl;j++) 
	{
	  checkKeyWord(info.keyWords_ecl[j],"(){}[]-_",msg,LINE_LENGTH);
	  if (massager->verbosity>2) printf("%12d %32s %s\n",j,info.keyWords_ecl[j], msg);
	  /* check if keywords contains segments */
	  if (massager->eclInfo.useECLIOsegments && 
		  !strncmp(info.keyWords_ecl[j],(*mog_i).ECLIPSE_segments_keyword,
				   strlen((*mog_i).ECLIPSE_segments_keyword))) 
		{
		  haveSegmentsKeyword=1;
		}
	  /* check if keywords contains layers */
	  if (!strcmp(info.keyWords_ecl[j],DESIRED_ECLIPSE_KEYWORDS[9])) haveLayersKeyword=1;
	}

  /* General method */
  /* load model */
  status = loadModelFromArray(filename, &info.ecl, transformToMap, info.nKeys_ecl, info.keyWords_ecl);  
  if (status) 
	{
	  printf("ERROR: Could not load model, exiting\n");
	  return EXIT_FAILURE; 
	}

  if (massager->makeThinBlockActive) 
	{
	  printf("Activating thin blocks in Eclipse grid\n");
	  if (activateThinActnums(&info.ecl,haveSegmentsKeyword) == EXIT_FAILURE)
		{
		  return EXIT_FAILURE;
		}
	}

  // set up coord lines and zCorns
  from = fetchProp(&info.ecl,"ZCORN");
  info.zCorn = ((double ****)(from->p));  // should be eclOut->zCorn


  // metricize if required
  if (massager->flipGridunits) flipGridAndMapUnits(&info.ecl);

  return EXIT_SUCCESS;
}

int exportCoarseGrid(struct MASSAGER *massager)
{
  int i,j,k,kk,kkk,c;
  int k_top,k_bot;
  int k_new_s,k_new_e;
  int status,exitStatus;
  ECL eclOut;
  double ****zD=NULL;
  PROP *zDP = NULL;
  int ***a = (int ***) NULL;
  int ***a_orig = (int ***) NULL;
  PROP *ap = (PROP *) NULL;
  PROP *ap_orig = (PROP *) NULL;


  exitStatus = EXIT_SUCCESS;
  
  /*
	Output of new grid
	Expect these:	 
	SPECGRID 
	PINCH 
	COORDSYS 
	GRIDUNIT 
	MAPAXES 
	COORD 
	ZCORN 
	ACTNUM 
	SEGMENTS 
	LAYERS ?
  */

  status =  initModel(&eclOut);
  if (status) {
    fprintf(stderr, "ERROR: Could not initialize eclOut\n");
    exitStatus = EXIT_FAILURE;
	goto bottom;
  }

  /* set (copy.c) some records , blithely discarding status flags */
  status = 0;
  status += copySpecGrid(&eclOut, &info.ecl); 

  // set new number of layers to the best we could extract from the well observations info
  k_new_s = k_new_e = -1;
  for (k=0;k<=massager->nLayers;k++)
	{
	  if (k_new_s==-1 && info.k_best[k]>0) 
		{
		  k_new_s = k;
		}
	  if (k_new_s!=-1)
		{
		  if (k_new_e!=-1 && info.k_best[k]<0) 
			{
			  // exiting good sequence
			  k_new_e = k-1;
			  break;
			}
		  if (info.k_best[k]>0)	k_new_e = k;
		}
	}

  // k_best has always one more entry than the number of layers we can define, so...
  eclOut.specGrid.nz = k_new_e - k_new_s;

  if (eclOut.specGrid.nz<1 || k_new_e==-1 || k_new_s==-1) 
	{
	  printf("Fewer than 1 layer's depth info available to write to new Eclipse grid.\n");
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
	}
  printf("New Eclipse grid spans massager layers indexed %d %d.\n",k_new_s,k_new_e);

  status += copyPinch(&eclOut, &info.ecl);
  status += copyCoordSys(&eclOut, &info.ecl);
  status += copyGridUnit(&eclOut, &info.ecl);
  status += copyMapAxes(&eclOut, &info.ecl); 
  if (status) 
	{
	  fprintf(stderr, "ERROR(s): Could not copy  %1d records to new Eclipse grid.\n",status);
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
	}

  // Copy remaining properties
  if (copyPropECL(&eclOut, "COORD", &info.ecl, "COORD", 0)) 
	{
	  fprintf(stderr, "Could not copy COORD to new Eclipse grid\n");
	  exitStatus = EXIT_FAILURE; 
	  goto bottom;
	}  

  // Add new ZCORN
  zDP = addProp(&eclOut,"ZCORN", 'd',4, 1, eclOut.specGrid.nx, 
				1,eclOut.specGrid.ny, 1,eclOut.specGrid.nz, 1,8); 
  if (!zDP || !zDP->p) 
	{
	  fprintf(stderr,"Could not fetch or set ZCORN in new Eclipse grid\n");
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
	}

  zD = (double ****) zDP->p;
  /* loop and set */
  for (i=1;i<=eclOut.specGrid.nx; i++)
	{
	  for(j=1; j<=eclOut.specGrid.ny; j++)
		{
		  for (k=1, kk=k_new_s; k<=eclOut.specGrid.nz; k++,kk++)
			{
			  k_top = info.k_best[kk];
			  // layer tops (block bottoms)
			  for (c=1;c<=4;c++)
				{
				  zD[i][j][k][c] =  (info.zCorn)[i][j][k_top][c];
				}
			  // layer bottoms (block tops)
			  k_bot = info.k_best[kk+1];
			  for (c=5;c<=8;c++)
				{
				  zD[i][j][k][c] =  (info.zCorn)[i][j][k_bot][c-4];
				}
			}
		}
	}

  // Add ACTNUM
  ap = addPropGrid(&eclOut, "ACTNUM", 'i');
  if (!ap || !ap->p) 
	{
	  fprintf(stderr,"Could not set ACTNUM in new Eclipse grid\n");
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
    }
  a = (int ***)ap->p;

  ap_orig = fetchProp(&info.ecl, "ACTNUM");
  if (!ap_orig || !ap_orig->p || !ap_orig->set)
	{
	  fprintf(stderr,"Could not get ACTNUM in source Eclipse grid\n");
	  exitStatus = EXIT_FAILURE;
	  goto bottom;;
	}
  a_orig= (int ***) ap_orig->p;

  for (i=1;i<=eclOut.specGrid.nx;i++) 
	{
	  for (j=1;j<=eclOut.specGrid.ny;j++) 
		{
		  for (k=1, kk=k_new_s; k<=eclOut.specGrid.nz; k++,kk++)
			{
			  // layer tops
			  k_top = info.k_best[kk];   //e.g. 90
			  // layer bottoms (block tops)
			  k_bot = info.k_best[kk+1]-1;  //e.g. (PO case) 90-1=89
			  // if info.k_best[kk+1] == info.k_best[kk], protect this..
			  // i.e. pinched out cases probably should be made active
			  k_bot = MAX(k_bot,k_top);
			  int isAct=0;
			  // scan subsumed blocks to see if * at least one * is active
			  for (kkk=k_top; kkk<=k_bot; kkk++)
				{
				  if (a_orig[i][j][kkk] != 0)
					{
					  isAct=1;
					  break;
					}
				}
			  a[i][j][k] = isAct;
			}
		}
	}
  ap->set = 1;
  (void) setAuxPointers(&eclOut);

  // Add SEGMENTS if possible
  ap_orig = fetchProp(&info.ecl, "SEGMENTS");
  // if not present in orig. skip this lot
  if (ap_orig)
	{
	  if (!ap_orig->p || !ap_orig->set)
		{
		  fprintf(stderr,"Could not get SEGMENTS in source Eclipse grid\n");
		  exitStatus = EXIT_FAILURE;
		  goto bottom;;
		}
	  a_orig= (int ***) ap_orig->p;


	  ap = addPropGrid(&eclOut, "SEGMENTS", 'i');
	  if (!ap || !ap->p) 
		{
		  fprintf(stderr,"Could not set SEGMENTS in new Eclipse grid\n");
		  exitStatus = EXIT_FAILURE;
		  goto bottom;
		}
	  a = (int ***)ap->p;

	  GHashTable *seg_hash = g_hash_table_new_full(g_int_hash,g_int_equal,NULL,free);
	  int majority_seg,max;
	  gint *value_ptr,value;
	  for (i=1;i<=eclOut.specGrid.nx;i++) 
		{
		  for (j=1;j<=eclOut.specGrid.ny;j++) 
			{
			  for (k=1, kk=k_new_s; k<=eclOut.specGrid.nz; k++,kk++)
				{
				  // layer tops
				  k_top = info.k_best[kk];
				  // layer bottoms (block tops)
				  k_bot = info.k_best[kk+1]-1;
				  // scan subsumed blocks to see if segment is constant
				  max=0;
				  majority_seg=a[i][j][k_top];
				  for (kkk=k_top; kkk<=k_bot; kkk++)
					{
					  if (!(value_ptr = (gint*)g_hash_table_lookup(seg_hash, (gint*)&a_orig[i][j][kkk])))
						{
						  value_ptr = (gint*)malloc(sizeof(gint));
						  *value_ptr = 0;
						  g_hash_table_insert(seg_hash,(gint*)&a_orig[i][j][kkk],value_ptr);
						}
					  (*value_ptr)++;
					  if ((*value_ptr) > max) 
						{
						  max = *value_ptr; 
						  majority_seg = a_orig[i][j][kkk];
						}
					}
				  a[i][j][k] = majority_seg;
				  // printf("SEG upscale: i=%d j=%d maj_seg=%d, fraction of sites=%g\n",
				  // i,j,majority_seg, (k_bot>k_top ? max/(k_bot-k_top+1.0) : 1));
				  // clear hash table
				  g_hash_table_foreach_remove (seg_hash, GHRCleanSegmentsFunc,NULL);
				}
			}
		}
	  ap->set = 1;
	  g_hash_table_destroy(seg_hash);
	}
 

  // Add LAYERS if possible
  ap_orig = fetchProp(&info.ecl, "LAYERS");
  // if not present in orig. skip this lot
  if (ap_orig)
	{
	  if (!ap_orig->p || !ap_orig->set)
		{
		  fprintf(stderr,"Could not get LAYERS in source Eclipse grid\n");
		  exitStatus = EXIT_FAILURE;
		  goto bottom;;
		}
	  a_orig= (int ***) ap_orig->p;


	  ap = addPropGrid(&eclOut, "LAYERS", 'i');
	  if (!ap || !ap->p) 
		{
		  fprintf(stderr,"Could not set LAYERS in new Eclipse grid\n");
		  exitStatus = EXIT_FAILURE;
		  goto bottom;
		}
	  a = (int ***)ap->p;

  
	  // just set as 1...Nz for now....
	  for (i=1;i<=eclOut.specGrid.nx;i++) 
		{
		  for (j=1;j<=eclOut.specGrid.ny;j++) 
			{
			  for (k=1; k<=eclOut.specGrid.nz; k++)
				{
				  a[i][j][k] = k;
				}
			}
		}
	  ap->set = 1;
	}



  /* now we write to the eclOut file */
  FILE *outEcl = NULL;
  char *filename_out=NULL;

  // use input grid name with new suffix
  asprintf(&filename_out,"%s.reduced.ecl",(massager->outputgrid[0]).filename);
  if (!(outEcl = fopen(filename_out,"w")))
	{
	  fprintf(stderr, "Couldn't open output file %s",filename_out);
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
	}


  int nKeyWords = info.nKeys_ecl; /* truncate the list for what we have available in the smashed model */
  status =  writeRecordList (outEcl, &eclOut, nKeyWords, info.keyWords_ecl);
  if (status) 
	{
	  fprintf(stderr, "Error encountered in writeRecordList: status = %d.\n", status);
	  exitStatus = EXIT_FAILURE;
	  goto bottom;
	}


  fclose(outEcl);

  printf("Wrote new Eclipse grid %s\n",filename_out);

  /* now get rid of eclOut data */
  freeModel(&eclOut);

  // should free up: info.keyWords_ecl

  // Cleanups
  free(filename_out);

 bottom:
  return exitStatus;
}

gboolean GHRCleanSegmentsFunc(gpointer key,gpointer value, gpointer user_data)
{
  return TRUE;
}

#if SHELL_SEARCH
void flush_plateau_list()
{
  int k;
  if (info.plateau_index_list->len > 0)
	{
	  // clean up plateau array. go from end to preserve order
	  for (k=info.plateau_index_list->len-1; k>=0; k--)
		{
		  int* list;
		  list = (int*)g_ptr_array_index(info.plateau_index_list,k);
		  free(list);
		  g_ptr_array_remove_index(info.plateau_index_list,k);
		}
	}
  // clear hash table
#if 0
  /* Simpler version available in newer glibs */
  g_hash_table_remove_all(info.plateau_hashtable);
#else
  /* old glib2 (pre 2.12 g_hash_table_remove_all) */
  g_hash_table_foreach_remove(info.plateau_hashtable, GHRCleanSegmentsFunc, NULL);
#endif
}

void append_to_plateau_list(int Nlayers)
{
  int k,*list;
  GString *gs;
  gpointer orig_key;
  gpointer value;
  gboolean found;

  // check if location is already on list
  // form key of type indx1.indx2.indx3.... 
  gs = g_string_new("");
  for (k=0; k<Nlayers; k++) 
	{
	  g_string_append_printf(gs,"%d.",info.k_search[k]);
	}
  g_string_append_printf(gs,"%d",info.k_search[Nlayers]);

  if (!g_hash_table_lookup_extended(info.plateau_hashtable,gs,&orig_key,&value))
	{
	  // Not in list. Insert location into hashtable
	  g_hash_table_insert(info.plateau_hashtable,gs,(gpointer)NULL);
	  // and store lattice site
	  list = (int*)malloc((Nlayers+1)*sizeof(int));
	  for (k=0;k<=Nlayers; k++) list[k] = info.k_search[k];
	  g_ptr_array_add(info.plateau_index_list,(gpointer)list);
	}
  else
	{
	  // already in list, free lookup key 
	  g_string_free(gs,TRUE);
	}
}

void string_hash_free(gpointer gs)
{
  g_string_free((GString*)gs,TRUE);
}

#endif
