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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>


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
#include "xmlParser.h"
#include "solver.h"
#include "util.h" 
#include "upperTriangleMatrices.h"
#include "eclipseHandler.h" 
#include "coarseGridMaker.h"
#include "su_util.h"


/* 
   Initialiser/constructor for struct MASSAGER massager
*/
void init_massager(struct MASSAGER *massager)
{
	  massager->verbosity=3;
	  massager->makeThinBlockActive=0;
	  massager->patchSegmentLabels=0;
	  massager->duplicate=0;
	  massager->writeAsciiOut = 0;
	  massager->writeSUOut = 0;
	  massager->su_out_info = NULL;
	  massager->make_ecl_grid = 0;
	  massager->Nrealisations=0;
	  massager->acquireStats=0;
	  massager->nDecorations=0;
	  massager->decoratorInfo.dumpPSdecoration=0;
	  massager->decoratorInfo.dumpPSdecorationVal=0;
	  massager->decoratorInfo.dumpPSdecorationDir='x';
	  massager->statsCutoff = -FLT_MAX;
	  massager->ignore_soft = 0;
	  massager->ignore_well_data = 0;
	  massager->eclInfo.writeECLOut  = 0;
	  massager->seed = -1;
	  massager->systematicDepthUncertainty=0;
	  massager->inputgrid.filt.apply=0;
	  massager->inputgrid.filt.propName = NULL;
	  massager->inputgrid.filt.min = 0;
	  massager->inputgrid.filt.max = 0;
	  massager->make_ecl_grid_options.make_ecl_options=NULL;
	  massager->make_ecl_grid_options.use_default=1;
	  massager->make_ecl_grid_options.Nu=-1;
	  massager->make_ecl_grid_options.Nv=-1;
	  massager->make_ecl_grid_options.MBR=0;

	  /* Memory safety */
	  massager->properties = NULL;
}

int parse_make_ecl_grid_options(struct MASSAGER *massager)
{
  int i;
  char r,*opts,*token;
  
  opts = massager->make_ecl_grid_options.make_ecl_options;

  for (i=0;;i++)
	{
	  token = strtok((i==0)?opts:NULL,",");
	  if (token == NULL) break;
	  r=0;
	  if (!strcmp(token,"d")) 
		{
		  massager->make_ecl_grid_options.use_default=1;
		  r=1;
		}
	  if (!strcmp(token,"MBR")) 
		{
		  massager->make_ecl_grid_options.use_default=0;
		  massager->make_ecl_grid_options.MBR=1;
		  printf("make-eclipse: Using minimum bounding rectangle\n");
		  r=1;
		}
	  if (!strncmp(token,"Nu=",3)) 
		{
		  massager->make_ecl_grid_options.Nu = atoi(&token[3]);
		  if (massager->make_ecl_grid_options.Nu <1)
			{
			  warning("Supplied make-eclipse option %s not a number>=1\n",token);
			  return EXIT_FAILURE;
			}
		  printf("make-eclipse: Using Nu=%d\n",massager->make_ecl_grid_options.Nu);
		  r=1;
		}
	  if (!strncmp(token,"Nv=",3)) 
		{
		  massager->make_ecl_grid_options.Nv = atoi(&token[3]);
		  if (massager->make_ecl_grid_options.Nv < 1)
			{
			  warning("Supplied make-eclipse option %s not a number>=1\n",token);
			  return EXIT_FAILURE;
			}
		  printf("make-eclipse: Using Nv=%d\n",massager->make_ecl_grid_options.Nv);
		  r=1;
		}
	  if (!r)
		{
		  warning("unrecognised make-eclipse option %s\n",token);
		  return EXIT_FAILURE;
		}
	}
  return EXIT_SUCCESS;
}
// read output grid using Chris' library
int checkOutputGridisInterior(struct MASSAGER *massager)
{
  int i,j;
  struct PT2D pt;
  // first run over all input points and see if they are 'interior'
  for (i=0; i < massager->nLayers;i++)
	{
	  for (j=0; j < massager->outputgrid[i].nPoints; j++) 
		{
		  pt = getBoundingIndices(&massager->inputgrid,
								  massager->outputgrid[i].xValues[j],
								  massager->outputgrid[i].yValues[j]);
		  if (pt.x<0)
			{
			  warning("Output grid contains x coord point < input grid\n");
			  // return EXIT_FAILURE;
			}
		  if (pt.x > massager->inputgrid.Nx-1)
			{
			  warning("Output grid contains x coord point > input grid\n");
			  // return EXIT_FAILURE;
			}
		  if (pt.y<0)
			{
			  warning("Output grid contains y coord point < input grid\n");
			  // return EXIT_FAILURE;
			}
		  if (pt.y > massager->inputgrid.Ny-1)
			{
			  warning("Output grid contains y coord point > input grid\n");
			  // return EXIT_FAILURE;
			}
		  // if we got here, model this point as obtainable from bilinear interpolation
		}
	}
  return EXIT_SUCCESS;
}

// Get input coordinate bounding indices of any arbitrary point
// Here, this is the corner of the box (i,j),(i+1,j),(i+1,j+1),(i,j+1)
// which encloses x,y: we get i,j by floor(indices.x,indices.y)
struct PT2D getBoundingIndices(struct INPUTGRID *ig,double x,double y)
{
  struct PT2D indices;
  double dotP;

  indices.x = (x - ig->rFirst.x)*ig->dRx.x 
	        + (y - ig->rFirst.y)*ig->dRx.y;

  indices.y = (x - ig->rFirst.x)*ig->dRy.x 
	        + (y - ig->rFirst.y)*ig->dRy.y;
  return indices;
}

// Set up mapping of input grid to output grid using kriging weights derived from
// variogram
int setUpGridMapping(struct MASSAGER *massager)
{
  int i,j;
  struct PT2D pt;

  if (massager->verbosity>1)
	printf("Setting up grid mapping (kriging)\n");

  return setUpMapping(massager);
}

int setUpOutputs(struct MASSAGER *massager)
{
  int i;
  char filename[100];

  /* set up possible output files */
  if (massager->make_ecl_grid)
	{
	  struct stat sb;
	  if (stat(massager->outputgrid[0].filename, &sb) == 0)
		{
		  printf("Making eclipse grid: file %s already exists: move it to safety before running!\n",
				 massager->outputgrid[0].filename);
		   printf("deliveryMassager will not overwrite ECLIPSE 'receptacle' files\n");
		   return EXIT_FAILURE;
		}
	  if (massager->outputgrid[0].ECLIPSE_segments_keyword == NULL)
		{
		  printf("Making eclipse grid: ECLIPSE_segments_keyword must be set: setting to SEGMENTS\n");
		  asprintf(&(massager->outputgrid[0].ECLIPSE_segments_keyword),"SEGMENTS");
		  massager->eclInfo.useECLIOsegments = 1;
		}
	}

  if (massager->writeAsciiOut)
	{
	  for (i=0; i < massager->nLayers;i++)
		{
		  sprintf(filename,"post_massaged_properties_layer_%d.txt",massager->layerNumbers[i]);
		  massager->ascii_out_filenames[i] = strdup(filename);
		  massager->ascii_out_fp[i] = fopen(massager->ascii_out_filenames[i],"w");
		  if (massager->ascii_out_fp[i] == NULL)
			{
			  printf("Cannot open output file %s for writing\n",massager->ascii_out_filenames[i]);
			  return EXIT_FAILURE;
			}
		  if (massager->Nrealisations>0)
			{
			  sprintf(filename,"post_massaged_realisations_layer_%d.txt",massager->layerNumbers[i]);
			  massager->ascii_out_realisations_filenames[i] = strdup(filename);
			  massager->ascii_out_realisations_fp[i] = fopen(massager->ascii_out_realisations_filenames[i],"w");
			  if (massager->ascii_out_realisations_fp[i] == NULL)
				{
				  printf("Cannot open output file %s for writing\n",massager->ascii_out_realisations_filenames[i]);
				  return EXIT_FAILURE;
				}			  
			}
		}
	}
  return EXIT_SUCCESS;
}

// Read output grid(s) using ECLIO library

int readOutputGrids(struct MASSAGER *massager)
{
  int i,j,k,i_dup,cnt,cnt_b;
  int usedCache=0;
  char c;
  ECL *ecl;
  int status,transformToMap,nKeys=0;
  int skipVoids,nrl=0,Nx=0,ncl=0,Ny=0,Nz=0,active_size,dup_active_size;
  double ***traces, **blockAreas,x,y,z,qtrA,A;
  int **active=NULL;
  int **segcol=NULL,***segments=NULL;
  char **keyWords = NULL;
  int haveSegmentsKeyword = 0;
  int haveLayersKeyword = 0;
  char msg[LINE_LENGTH];
  char *s=NULL;
  struct OUTPUTGRID *mog_i;

#if DUPLICATE_AT_FAULTS
  struct SEG_ENCODING **segInfo=NULL;
#endif

  // read and parse here

  massager->totalNOutputGridPoints=0;
  
 if (massager->make_ecl_grid)
   {
	 // Make eclipse grid, so it's ready to read in below
	 if (makeEclipseGrid(massager) == EXIT_FAILURE) return EXIT_FAILURE;
	 // Overwrite any exotic XML settings to do with eclipse file
	 for (i=0; i < massager->nLayers; i++)
	   {
		 massager->outputgrid[i].ECLIPSE_layer_number = (i+1);
		 massager->outputgrid[i].ECLIPSE_reference_depth_layer_number = 1;
		}
	 // and other general cleanups, as per XML parser
	 if (isEclipseCommonSource(massager) == EXIT_FAILURE) return EXIT_FAILURE;
   }


 if (massager->duplicate)  // never occurs with massager->make_ecl_grid
	{
	  for (i=0; i< massager->nLayers; i++)
		{
		  mog_i = &(massager->outputgrid[i]);

		  // Example of a simple case of 2 rectangular grids, mapping input to output
		  (*mog_i).isRectangular=1;
		  (*mog_i).nPoints = massager->inputgrid.nPoints;
		  (*mog_i).xValues = massager->inputgrid.xValues;
		  (*mog_i).yValues = massager->inputgrid.yValues;
		  (*mog_i).areas = NULL;
		  massager->totalNOutputGridPoints = (*mog_i).nPoints; // sum over all output grids
		}
	}
  else
	{
	  // read eclipse files etc

	  for (i=0; i< massager->nLayers; i++)
		{
		  mog_i = &(massager->outputgrid[i]);

		  if (massager->eclInfo.commonEclipseGrid)
			{
			  // struct mallinfo mall_info = mallinfo();

			  // common eclipse grid. Only case coded so far....
			  usedCache  = readEclipseCache(massager,i);
			  if (!usedCache)
				{
				  if (i==0)
					{
					  // If common, only need to read grid once
					  
					  /* read Eclipse grids */
					  transformToMap = 1; /* number of records to read and
											 force transform grid->map coords */
					  if (massager->verbosity>1) printf("...Layer %d, loading Eclipse model from file %s\n", 
														massager->layerNumbers[i],(*mog_i).filename);
					  nKeys = fetchKeyList((*mog_i).filename, &keyWords, "--");

					  if (nKeys<=0) return EXIT_FAILURE;

					  /* List possible keywords */
					  if (massager->verbosity>2)
						{
						  fprintf(stdout, 
								  "\nFound %d *probable* keys in file %s\n"
								  "Key Number   Key Value                        Message\n"
								  "------------ -------------------------------- ------- \n",
								  nKeys,(*mog_i).filename);
						}
					  for (j=1;j<=nKeys;j++) 
						{
						  checkKeyWord(keyWords[j],"(){}[]-_",msg,LINE_LENGTH);
						  if (massager->verbosity>2) fprintf(stdout,"%12d %32s %s\n",j,keyWords[j], msg);
						}
					  // reduce list to useful ones...
					  for (j=1,k=1;j<=nKeys;j++) 
						{ 
						  if (isADesiredECLIPSEKeyword(keyWords[j])) 
							{
							  //asprintf(&s,"%s",keyWords[j]);
							  keyWords[k++]=keyWords[j];
							}
						  else
							{
							  free(keyWords[j]);keyWords[j]=NULL;
							}
						}
					  // required keywords are now indexed j=1...(k-1)
					  nKeys=k-1;
					  
					  /* After a bit of cleaning -- if you get a message here, program will probably crash 
						 when loading ... will provide a cushier exit later*/
					  nKeys = editKeyList(nKeys, &keyWords);
					  if (massager->verbosity>2)
						{
						  fprintf (stdout, 
								   "\nAfter editKeyList: %d *probable* keys in file %s\n"
								   "Key Number   Key Value                        Message\n"
								   "------------ -------------------------------- ------- \n",
								   nKeys,(*mog_i).filename);
						}
					  /* Check we have segments too */ 
					  haveSegmentsKeyword=0;
					  haveLayersKeyword=0;
					  for (j=1;j<=nKeys;j++) 
						{
						  checkKeyWord(keyWords[j],"(){}[]-_",msg,LINE_LENGTH);
						  if (massager->verbosity>2) fprintf(stdout,"%12d %32s %s\n",j,keyWords[j], msg);
						  if (massager->eclInfo.useECLIOsegments && 
							  !strncmp(keyWords[j],(*mog_i).ECLIPSE_segments_keyword,strlen((*mog_i).ECLIPSE_segments_keyword))) 
							{
							  haveSegmentsKeyword=1;
							}
						  /* check if keywords contains layers */
						  if (!strcmp(keyWords[j],DESIRED_ECLIPSE_KEYWORDS[9])) haveLayersKeyword=1;
						} 
					  if (massager->eclInfo.useECLIOsegments && !haveSegmentsKeyword)
						{
						  printf("ERROR: Eclipse file %s does not contain segments keyword %s, exiting\n",
								 (*mog_i).filename,
								 (*mog_i).ECLIPSE_segments_keyword);
						  return EXIT_FAILURE;
						}

					  /* General method */
					  /* load model */
					  status = loadModelFromArray((*mog_i).filename, &(*mog_i).ecl, transformToMap, nKeys, keyWords);  
					  if (status) 
						{
						  printf("ERROR: Could not load model, exiting\n");
						  return EXIT_FAILURE; 
						}
					  ecl = &(*mog_i).ecl;


					  // metricize if required
					  if (massager->flipGridunits) flipGridAndMapUnits(ecl);

					  // fix actives immediately if desired
					  // Activate thin deactivated blocks
					  if (massager->patchSegmentLabels) 
						{
						  printf("Patching segment labels in Eclipse grid\n");
						  if (patchEclSegmentLabels(massager,ecl,haveSegmentsKeyword) == EXIT_FAILURE)
							{
							  return EXIT_FAILURE;
							}
						}
					  if (massager->makeThinBlockActive) 
						{
						  printf("Activating thin blocks in Eclipse grid\n");
						  if (activateThinActnums(ecl,haveSegmentsKeyword) == EXIT_FAILURE)
							{
							  return EXIT_FAILURE;
							}
						}
					}
				  else
					{
					  // subsequent loads: use first grid. This is for (i>0 && massager->eclInfo.commonEclipseGrid)
					  ecl = &(massager->outputgrid[0].ecl);
					  mog_i->ecl = massager->outputgrid[0].ecl;
					}

				  skipVoids = 1; /* do include void blocks in gross count */
				  nrl = 0; 
				  Nx = getNX(ecl); /* subranges n*l > 0 and n*h < nx are allowed */
				  ncl = 0; 
				  Ny = getNY(ecl);  
				  Nz = getNZ(ecl);
				  massager->outputgrid[i].ecl_Nx = Nx;
				  massager->outputgrid[i].ecl_Ny = Ny; 

				  // Simple sanity check
				  if ((Nz < massager->nLayers) && massager->eclInfo.writeECLOut)
					{
					  printf("Too few layers (%d) in Eclipse grid %s to match to %d layers in XML file\n",
							 Nz,(*mog_i).filename,massager->nLayers);
					  return EXIT_FAILURE;
					}
				  if ((*mog_i).ECLIPSE_layer_number > Nz)
					{
					 printf("Too few layers (%d) in Eclipse grid %s to match required Eclipse_layer_number %d\n",
							 Nz,(*mog_i).filename,(*mog_i).ECLIPSE_layer_number);
					  return EXIT_FAILURE; 
					}


				  traces = fetchLayerTraces(ecl, (*mog_i).ECLIPSE_layer_number, skipVoids, 0, Nx, 0, Ny);
				  if (!traces) 
					{ 
					  printf("ERROR: Returned null trace array\n");
					  return EXIT_FAILURE;
					}

				  if (i>0 && massager->eclInfo.commonEclipseGrid)
					{
					  // active column list from first grid can be reused
					  massager->outputgrid[i].ecl_active = massager->outputgrid[0].ecl_active;
					  active = massager->outputgrid[i].ecl_active;
					}
				  else
					{					  
					  printf("...fetching active trace indicators for model %s\n",(*ecl).name); 
					  active = fetchActCoord(ecl, 0, Nx, 0, Ny);
					  if (!active)
						{
						  printf("ERROR: Returned null active array\n");
						  return EXIT_FAILURE;
						}
					  printf("...fetched active array for model %s\n",(*ecl).name); 
					  massager->outputgrid[i].ecl_active = active;
					}
					
				  blockAreas = fetchTraceAreas(ecl, 0, Nx, 0, Ny); // areas coupled to fetchLayerTraces() call
				  // FIXME MAKE sure areas alter with layer number!!!!!!!!!!!!!!!
				  if (!blockAreas) 
					{ 
					  printf("ERROR: Returned null trace-centered block area array\n");
					  return EXIT_FAILURE;
					}

				  /* get grid size */
				  active_size = countActCoord(ecl,0);
				  printf("Active grid size = %d\n",active_size);

#if DUPLICATE_AT_FAULTS
				  if (massager->eclInfo.useECLIOsegments)
					{
					  if (i==0)
						{
						  // segInfo common if eclInfo.commonEclipseGrid
						  segInfo = fetchSegLists(ecl,&dup_active_size,massager->verbosity);
						  segcol = fetchSegCol(ecl, 1, Nx, 1, Ny);
						}
					  printf("With fault duplication, active grid size = %d\n",dup_active_size);
					  active_size = dup_active_size;
					  (*mog_i).segInfo = segInfo;
					  (*mog_i).segcol = segcol;
					}
#endif

				  /* alloc and load */ 
				  (*mog_i).isAMirror=0;
				  (*mog_i).isRectangular=0; 
				  (*mog_i).nPoints = active_size;
				  (*mog_i).xValues = dvector(0,active_size-1);
				  (*mog_i).yValues =  dvector(0,active_size-1);
				  (*mog_i).zValues = NULL;
				  (*mog_i).zValues =  dvector(0,active_size-1);
				  (*mog_i).areas = NULL;
				  (*mog_i).areas =  vector(0,active_size-1);
				  (*mog_i).ecl_i = NULL;
				  (*mog_i).ecl_i =  ivector(0,active_size-1);
				  (*mog_i).ecl_j = NULL;
				  (*mog_i).ecl_j =  ivector(0,active_size-1);

				  massager->totalNOutputGridPoints += (*mog_i).nPoints; // sum over all output grids
				  /* segments */
				  (*mog_i).segnum = NULL;
				  (*mog_i).segnum_corner_encoding = NULL;
				  if (massager->eclInfo.useECLIOsegments)
					{
					  (*mog_i).segnum = ivector(0,active_size-1);					  
#if DUPLICATE_AT_FAULTS
					  (*mog_i).segnum_corner_encoding = cvector(0,active_size-1);
#else
					  /* fetch segments */
					  segments = fetchSegments(ecl, 1, Nx, 1, Ny, 1, Nz); 
					  segcol = fetchSegCol(ecl, 1, Nx, 1, Ny);
#endif
					}
				  else (*mog_i).segnum = NULL;

				  cnt=0;
				  for (j=0;j<=Nx;j++)
					{
					  for (k=0;k<=Ny;k++)
						{
						  if (active[j][k]) 
							{
							  if (massager->eclInfo.useECLIOsegments)
								{
#if DUPLICATE_AT_FAULTS
								  /* midpoints of traces: 6 values are x,y,z_top,x,y,z_bottom */
								  x = 0.5*(traces[j][k][1]+traces[j][k][4]);
								  y = 0.5*(traces[j][k][2]+traces[j][k][5]);
								  z = 0.5*(traces[j][k][3]+traces[j][k][6]);
								  qtrA = 0.25*blockAreas[j][k];  // 1/4 of area centred at this trace

								  // Divvy up the area centred at this trace to the duplicated nodes, based
								  // on segment representation. Total area of duplicated points will sum to 
								  // 4*qtrA == blockAreas[j][k]
								  // In most cases, there are no duplicates (segInfo[j][k].size==1)
								  // so all the area will go to the single point. The corner encoding will of
								  // course be 2^1+2^2+2^3+2^4 == 30, so we get 4 contributions to the area

								  for (i_dup=0; i_dup < segInfo[j][k].size; i_dup++)
									{
									  (*mog_i).xValues[cnt] = x;
									  (*mog_i).yValues[cnt] = y;
									  (*mog_i).zValues[cnt] = z;
									  (*mog_i).ecl_i[cnt] = j;
									  (*mog_i).ecl_j[cnt] = k;
									  (*mog_i).segnum[cnt] = segInfo[j][k].values[i_dup];
									  (*mog_i).segnum_corner_encoding[cnt] = segInfo[j][k].corner_encoding[i_dup];
									  // allocate area based on number of segments
									  for (c=1,A=0.0; c<=4; c++)
										{
										  if ((char)((*mog_i).segnum_corner_encoding[cnt] & (char)(1<<c))) A+=qtrA;
										}
									  (*mog_i).areas[cnt] = A;
									  cnt++;
									}
#else
								  /* midpoints of traces: 6 values are x,y,z_top,x,y,z_bottom */
								  (*mog_i).xValues[cnt] = 0.5*(traces[j][k][1]+traces[j][k][4]);
								  (*mog_i).yValues[cnt] = 0.5*(traces[j][k][2]+traces[j][k][5]);
								  (*mog_i).zValues[cnt] = 0.5*(traces[j][k][3]+traces[j][k][6]);	
								  (*mog_i).ecl_i[cnt] = j;
								  (*mog_i).ecl_j[cnt] = k;
								  (*mog_i).areas[cnt] = blockAreas[j][k];
								  (*mog_i).segnum[cnt] = (j>0 && k>0) ? segcol[j][k] : -1;
								  if (massager->verbosity>4)
									{
									  printf("OGSEG %15.10g %15.10g %3d\n",(*mog_i).xValues[cnt],
											 (*mog_i).yValues[cnt],(*mog_i).segnum[cnt]);
									}
								  cnt++;
#endif
								}
							  else 
								{
								  (*mog_i).ecl_i[cnt] = j;
								  (*mog_i).ecl_j[cnt] = k;
								  (*mog_i).areas[cnt] = blockAreas[j][k];
								  if (massager->verbosity>3)
									{
									  printf("GRID %15.10g %15.10g\n",(*mog_i).xValues[cnt],(*mog_i).yValues[cnt]);
									}
								  cnt++;
								}
							}
						}
					}
				  
				  free_dmatrix(blockAreas, 0, Nx, 0, Ny);
				  blockAreas = NULL;
				  free_d3tensor(traces, 0, Nx, 0, Ny,1,6);
				  traces = NULL;
				  /* Clean ups */
				  if (segments!=NULL) 
					{
					  free_i3tensor(segments,1,Nx, 1, Ny, 1, Nz);
					  segments=NULL;
					}

				  // write cache to save time on next read !!!
				  writeEclipseCache(massager,i);
				}
			  if (massager->verbosity>1) printf("...finished loading model from file %s\n", (*mog_i).filename); 
			}
		  else
			{
			  // HERE! Code case of multiple eclipse input grids
			  // For now, signal failure: (XML parser forces same)
			  return EXIT_FAILURE;
			}
		}


	  // Set up output eclipse grids for mean model

	  if (determineEclipseOutputMode(massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  for (i=0; i< massager->nLayers; i++)
		{ 
		  mog_i = &(massager->outputgrid[i]);

		  /* Set up ECLIPSE output grid */
		  if (massager->eclInfo.writeECLOut)
			{
			  if ((i==0) || !massager->eclInfo.writeInSitu) 
				{
				  if ((*mog_i).eclOutInfo.outputFilename != NULL)
					{
					  ecl = &massager->outputgrid[i].ecl;

					  int status;
					  /* Now some write stuff ... this uses mainly high-level routines
						 from write.c, but also uses some low-level routines in set.c 
						 A "squashModel()" or "copyModel()" routine could be added to 
						 access.c
					  */
					  status =  initModel(&(*mog_i).eclOutInfo.eclOut);
					  if (status) 
						{
						  fprintf(stderr, "ERROR: Could not initialize eclOut\n");
						  return EXIT_FAILURE; /* could clean up here.... */
						}

					  /* set (copy.c) some records , blithely discarding status flags */
					  status = 0;
					  status += copySpecGrid(&(*mog_i).eclOutInfo.eclOut, ecl); 
					  if (!massager->eclInfo.writeInSitu)
						{
						  // NOT in-situ: a new, nz=1 collapsed grid
						  (*mog_i).eclOutInfo.eclOut.specGrid.nz = 1;
						}



					  status += copyPinch(&(*mog_i).eclOutInfo.eclOut, ecl);
					  status += copyCoordSys(&(*mog_i).eclOutInfo.eclOut, ecl);
					  status += copyGridUnit(&(*mog_i).eclOutInfo.eclOut, ecl);
					  status += copyMapAxes(&(*mog_i).eclOutInfo.eclOut, ecl);
					  if (status) 
						{
						  fprintf(stderr, "ERROR(s): Could not copy  %1d records.\n",status);
						  return EXIT_FAILURE; /* could clean up here.... */
						}
				  
					  if (massager->eclInfo.writeInSitu)
						{
						  /* In situ grids */

						  /* Copy ACTNUM. This routine grabs data from the other model, forcing to new grid (last arg is flag) */
						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "ACTNUM",  ecl, "ACTNUM",1)) 
							{
							  fprintf(stderr, "Could not copy ACTNUM to mean model\n");
							  return EXIT_FAILURE;
							}
						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "COORD", ecl, "COORD", 0)) 
							{
							  fprintf(stderr, "Could not copy COORD to mean model\n"); 
							  return EXIT_FAILURE; 
							}
						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "SEGMENTS", ecl, "SEGMENTS", 0)) 
							{
							  fprintf(stderr, "Could not copy SEGMENTS to mean model\n");
							  return EXIT_FAILURE; 
							}

						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "ZCORN", ecl, "ZCORN", 0)) 
							{
							  fprintf(stderr, "Could not copy ZCORN to mean model.\n");
							  return EXIT_FAILURE;
							} 
						  if (haveLayersKeyword)
							{
							  (*mog_i).eclOutInfo.haveLayers=1;
							  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "LAYERS", ecl, "LAYERS", 0))
								{
								  fprintf(stderr, "Could not copy LAYERS to mean model.\n");
								  return EXIT_FAILURE; 
								}
							}
					  
						  /* set up bounds on array */
						  (*mog_i).eclOutInfo.bounds = imatrix(1,3,1,2);
						  int** bounds = (*mog_i).eclOutInfo.bounds;
						  bounds[1][1] = 0; bounds[1][2] = Nx-1;
						  bounds[2][1] = 0; bounds[2][2] = Ny-1;
						  bounds[3][1] = 0; bounds[3][2] = Nz-1;

						  /* Now we set up keywords and file I/O to the eclOut file */
						  /* truncate the list for what we have available in the reduced model */
						  /* I think 8 is the minimum we need */
						  /* Write out  the kriged properties plus their std deviations */
						  (*mog_i).eclOutInfo.nDuplicateKeys = 9;

						  // plus layers if present
						  if ((*mog_i).eclOutInfo.haveLayers)  (*mog_i).eclOutInfo.nDuplicateKeys++;

						  (*mog_i).eclOutInfo.nGridBlockKeys = 2*massager->nPropertyNames;
						  // thickness and std_dev thickness now supplied
						  /*
							if (massager->thickness_index >= 0)
							{
							  // omit thickness: it's a geometry property
							  (*mog_i).eclOutInfo.nGridBlockKeys -= 2;
							}
						  */
						  if (massager->d_index >= 0)
							{
							  // Omit depth: it's a geometry property. But we *will* supply STDDEV_d
							  (*mog_i).eclOutInfo.nGridBlockKeys -= 1;
							}

						  (*mog_i).eclOutInfo.nKeys = (*mog_i).eclOutInfo.nDuplicateKeys + (*mog_i).eclOutInfo.nGridBlockKeys;

						  (*mog_i).eclOutInfo.keyWords = 
							(char**)malloc((size_t)(((*mog_i).eclOutInfo.nKeys+1) * sizeof(char*)));

						  /* duplicated keys */
						  if (massager->verbosity>2)
							{
							  fprintf(stdout,"Keys written to ECLIPSE grids:\n");
							  fprintf (stdout, 
									   "Key Number   Key Value                       \n"
									   "------------ --------------------------------\n");
							}

						  char* s;
						  for (j=1;j<=(*mog_i).eclOutInfo.nDuplicateKeys; j++)
							{
							  asprintf(&s,"%s",keyWords[j]);
							  (*mog_i).eclOutInfo.keyWords[j] = s;
							  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
				  
							}

						  /* concatenate properties. We'll write the kriged means plus std deviations */

						  for (k=0;k<massager->nPropertyNames; k++)
							{
							  /* skip depth as we don't need it. But we'd like to see the std deviation */
							  if (strcmp(massager->propertyNames[k],"d"))
								{
								  /* skip depth (Already easily viewed in Petrel). But write its uncertainty (below) */
								  /* property */
								  (*mog_i).eclOutInfo.keyWords[j] = massager->propertyNames[k];
								  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
								  j++;
								}

							  /* std deviation */
							  asprintf(&s,"STDDEV_%s",massager->propertyNames[k]);
							  (*mog_i).eclOutInfo.keyWords[j] = s;
							  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
							  j++;
							}
						  // allocate block-centred values.
						  // first index is property, then regular contiguous cubes i,j,k */
						  if ((*mog_i).eclOutInfo.nGridBlockKeys > 0)
							{
							  // Nx = getNX(ecl) etc
							  massager->eclInfo.property = f4tensor(0,(*mog_i).eclOutInfo.nGridBlockKeys-1, 0,Nx-1, 0,Ny-1, 0,Nz-1);
							  f4tensor_fill(D_MISSING,massager->eclInfo.property,(*mog_i).eclOutInfo.nGridBlockKeys,Nx,Ny,Nz);
							}
						  else massager->eclInfo.property = NULL;

						  /* keep track of actives */
						  (*mog_i).eclOutInfo.active = active;

						  /* open file for writing */
						  (*mog_i).eclOutInfo.ecl_out_fp = fopen((*mog_i).eclOutInfo.outputFilename,"w");
						  if ((*mog_i).eclOutInfo.ecl_out_fp == NULL)
							{
							  printf("Cannot open output file %s for writing\n",(*mog_i).eclOutInfo.outputFilename);
							  return EXIT_FAILURE;
							}
						}
					  else // else squashed
						{					  
						  /* SQUASHED GRIDS: Not complete */

						  /* redo the zcorn array for the squashed grid (must set specGrid, as above) */
						  if (squashZCorn(&(*mog_i).eclOutInfo.eclOut, ecl)) 
							{
							  fprintf(stderr, "Could not squash the original to a 2D eclipse grid.\n");
							  return EXIT_FAILURE; /* could clean up here.... */
							}

						  /* This routine grabs data from the other model, forcing to new grid (last arg is flag) */
						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "ACTNUM",  ecl, "ACTCOL",1)) 
							{
							  fprintf(stderr, "Could not copy ACTNUM to ACTCOL,\n");
							  return EXIT_FAILURE; 
							}
						  if (copyPropECL(&(*mog_i).eclOutInfo.eclOut, "COORD", ecl, "COORD", 0)) 
							{
							  fprintf(stderr, "Could not copy COORD to COORD (orig->squashed).\n");
							  return EXIT_FAILURE; 
							}

						  /* Now copy the smashed segments arrays over to ECLOUT */
						  /* NOTE: copies are done via vector operations assuming contiguous storage */
						  /* 
							 Here is the syntax for an alternative form, where the data reside in an 
							 array not in a ecl->prop * ... still a vector contiguous copy
							 NOTE that the array is passed by its name (not the first data location), 
							 with indexing data in nInd (here, 2) and bnd
						  */
					  
						  /* set up bounds on array */
						  (*mog_i).eclOutInfo.bounds = imatrix(1,3,1,2);
						  int** bounds = (*mog_i).eclOutInfo.bounds;
						  bounds[1][1] = 0; bounds[1][2] = Nx-1;
						  bounds[2][1] = 0; bounds[2][2] = Ny-1;
						  bounds[3][1] = 0; bounds[3][2] = 0;

					  
						  /* Now we set up keywords and file I/O to the eclOut file */
						  /* truncate the list for what we have available in the smashed model */
						  /* I think 8 is the minimum we need */
						  /* Write out  the kriged properties plus their std deviations */
						  (*mog_i).eclOutInfo.nDuplicateKeys = 8;
						  (*mog_i).eclOutInfo.nKeys = (*mog_i).eclOutInfo.nDuplicateKeys + 2*massager->nPropertyNames;
						  (*mog_i).eclOutInfo.keyWords = 
							(char**)malloc((size_t)(((*mog_i).eclOutInfo.nKeys+1) * sizeof(char*)));
						  /* duplicated keys */
						  if (massager->verbosity>2)
							{
							  fprintf(stdout,"Keys written to ECLIPSE grids:\n");
							  fprintf (stdout, 
									   "Key Number   Key Value                       \n"
									   "------------ --------------------------------\n");
							} 

						  char* s;
						  for (j=1;j<=(*mog_i).eclOutInfo.nDuplicateKeys; j++)
							{
							  asprintf(&s,"%s",keyWords[j]);
							  (*mog_i).eclOutInfo.keyWords[j] = s;
							  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
				  
							}
						  /* concatenate properties. We'll write the kriged means plus std deviations */

						  for (k=0;k<massager->nPropertyNames; k++,j++)
							{
							  /* property */
							  (*mog_i).eclOutInfo.keyWords[j] = massager->propertyNames[k];
							  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
							  /* std deviation */
							  asprintf(&s,"STDDEV_%s",(*mog_i).eclOutInfo.keyWords[j]);
							  (*mog_i).eclOutInfo.keyWords[++j] = s;
							  if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,(*mog_i).eclOutInfo.keyWords[j]);
							}
						  /* allocate. first index is property, then regular grids i,j */
						  // (*mog_i).eclOutInfo.property = f3tensor(0,2*massager->nPropertyNames-1,0,massager->nLayers-1,0,Nx,0,Ny);

						  /* keep track of actives */
						  (*mog_i).eclOutInfo.active = active;

						  /* open file for writing */
						  (*mog_i).eclOutInfo.ecl_out_fp = fopen((*mog_i).eclOutInfo.outputFilename,"w");
						  if ((*mog_i).eclOutInfo.ecl_out_fp == NULL)
							{
							  printf("Cannot open output file %s for writing\n",(*mog_i).eclOutInfo.outputFilename);
							  return EXIT_FAILURE;
							}
						}
					}
				} // alloc new eclio output grids
			}  // write ECLIO		  
		}  // layers
	  if (massager->eclInfo.commonEclipseGrid)
		{	
		  // usual clean up
		  if (keyWords)
			{
			  for (j=1;j<=nKeys;j++) 
				{
				  free(keyWords[j]);
				}
			  keyWords=NULL;	  
			}
		}
	}
  if (massager->verbosity>1)
	{
	  for (i=0; i< massager->nLayers; i++)
		{
		  printf("Total number of output grid points on layer %d = %d\n",i,massager->outputgrid[i].nPoints);
		}
	}


  return EXIT_SUCCESS;
} 

int isEclipseCommonSource(struct MASSAGER *massager)
{
  int i;
  char* word=NULL;
  struct OUTPUTGRID *mog_i;

  // is the eclipse source grid the same?
  if (massager->duplicate)
	{
	  // shouldn't matter
	  massager->eclInfo.commonEclipseGrid = 0;
	}
  else
	{
	  massager->eclInfo.commonEclipseGrid = 1;
	  if (massager->outputgrid[0].filename == NULL)
		{
		  printf("<ECLIPSE_grid_filename> unspecified in XML file %s for layer %s\n",
				 massager->massage_xml_file, massager->outputgrid[0].layer_name);
		  return EXIT_FAILURE;
		}
	  for (i=1; i< massager->nLayers; i++)
		{
		  if (massager->outputgrid[i].filename == NULL)
			{
			  // not supplied in xml => assume duplicate
			  massager->outputgrid[i].filename = massager->outputgrid[i-1].filename;
			}
		  else
			{
			  if (strcmp(massager->outputgrid[i].filename, massager->outputgrid[i-1].filename))
				{
				  massager->eclInfo.commonEclipseGrid = 0;
				  break;
				}
			}
		}

	  // checks
	  if (massager->eclInfo.commonEclipseGrid)
		{
		  // (1) reference layer check
		  // ECLIPSE_reference_depth_layer_number. Required for segmentations, even if no 
		  // ECLIPSE output => Check it

		  // simulating with thickness only (no depths) call
		  int layer=0,layerRef,last_layer;

		  layerRef = massager->outputgrid[0].ECLIPSE_reference_depth_layer_number;
		  if (layerRef == INT_MIN)
			{
			  printf("Unspecified ECLIPSE_reference_depth_layer_number value in first layer of XML\n");
			  return EXIT_FAILURE;
			}
		  if (layerRef <= 0)
			{
			  printf("Illegal ECLIPSE_reference_depth_layer_number value in first layer of XML\n");
			  return EXIT_FAILURE;
			}
		  for (i=0; i< massager->nLayers; i++)
			{
			  layer = massager->outputgrid[i].ECLIPSE_reference_depth_layer_number;
			  if (layer < 0)
				{
				  // not specified: force match..
				  massager->outputgrid[i].ECLIPSE_reference_depth_layer_number = layerRef;
				  layer = layerRef;
				}
			  if (layer != layerRef)
				{
				  printf("ECLIPSE_reference_depth_layer_number %d in XML doesn't match %d in first layer "
						 "for common ECLIPSE grid\n",layer,layerRef);
				  return EXIT_FAILURE;
				}
			}

		  
		  // (2) reference seg keyword check
		  if (massager->eclInfo.useECLIOsegments)
			{
			  word = massager->outputgrid[0].ECLIPSE_segments_keyword;
			  if (!word)
				{
				  printf("Error: ECLIPSE_segments_keyword not specified in first layer of XML\n");
				  return EXIT_FAILURE;
				}
			  for (i=1; i< massager->nLayers; i++)
				{
				  mog_i = &(massager->outputgrid[i]);
				  if (mog_i->ECLIPSE_segments_keyword == NULL )
					{
					  mog_i->ECLIPSE_segments_keyword = word;
					}
				  else
					{
					  if (strcmp(mog_i->ECLIPSE_segments_keyword,word))
						{
						  printf("ECLIPSE_segments_keyword %s in layer %d in XML doesn't match first layer "
								 "for common ECLIPSE grid\n",mog_i->ECLIPSE_segments_keyword ,layer+1);
						  return EXIT_FAILURE;
						}
					}
				}
			}


		  // Eclipse *Output* Checks		 
		  if (massager->eclInfo.writeECLOut)
			{
			  // (a) Eclipse layer numbers
			  if (massager->outputgrid[0].ECLIPSE_layer_number < 0)
				{
				  // not specified: error.
				  printf("ECLIPSE_layer_number not specified in XML file %s for layer %s\n",
							 massager->massage_xml_file,massager->outputgrid[0].layer_name);
				  return EXIT_FAILURE;
				}
			  last_layer = massager->outputgrid[0].ECLIPSE_layer_number;

			  for (i=1; i< massager->nLayers; i++)
				{
				  layer = massager->outputgrid[i].ECLIPSE_layer_number;
				  if (layer < 0)
					{
					  // Not specified
					  // can just increment
					  massager->outputgrid[i].ECLIPSE_layer_number = 
						massager->outputgrid[i-1].ECLIPSE_layer_number + 1;
					  layer = massager->outputgrid[i].ECLIPSE_layer_number;
					  printf("Warning: ECLIPSE_layer_number not specified in XML file %s\n"
							 "         for layer %s. Assuming ECLIPSE_layer_number=%d\n",
							 massager->massage_xml_file,massager->outputgrid[i].layer_name,layer);
					}
				  else
					{
					  if (layer != (last_layer+1))
						{
						  printf("ECLIPSE_layer_numbers are not contiguous (ascending by 1's) in XML file, "
								 "for layer named %s\n",massager->outputgrid[i].layer_name);
						  return EXIT_FAILURE;
						}
					}
				  last_layer = layer;
				}

			  // (b) Eclipse outputgrid filename
			  word = massager->outputgrid[0].eclOutInfo.outputFilename;
			  if (!word)
				{
				  printf("Error: ECLIPSE_output_grid_filename not specified in first layer of XML\n");
				  return EXIT_FAILURE;
				}
			  for (i=1; i< massager->nLayers; i++)
				{
				  mog_i = &(massager->outputgrid[i]);
				  if (mog_i->eclOutInfo.outputFilename == NULL )
					{
					  mog_i->eclOutInfo.outputFilename = word;
					}
				  else
					{
					  if (strcmp(mog_i->eclOutInfo.outputFilename,word))
						{
						  printf("ECLIPSE_output_grid_filename %s in layer %d in XML doesn't match first layer "
								 "for common ECLIPSE grid\n",mog_i->eclOutInfo.outputFilename,layer+1);
						  return EXIT_FAILURE;
						}
					}
				}
			}
		  
		}
	}
  return EXIT_SUCCESS;
}

int determineEclipseOutputMode(struct MASSAGER *massager)
{
  int i,Nz,layer;
  struct OUTPUTGRID *mog_i;

  if (massager->eclInfo.writeECLOut)
	{
	  // is the source grid the same?
	  // Write in-situ if all layers accessible and the ecl source is common
	  if (massager->eclInfo.commonEclipseGrid)
		{
		  // assume we can, till disproved
		  massager->eclInfo.writeInSitu = 1;
		  Nz = getNZ(&(massager->outputgrid[0].ecl));
		  for (i=0; i< massager->nLayers; i++)
			{
			  // check of named source layer
			  layer = massager->outputgrid[i].ECLIPSE_layer_number;
			  if ((layer<=0) || ( layer > Nz))
				{
				  massager->eclInfo.writeInSitu = 0;
				  break;
				}			  
			  if (massager->d_index==-1)
				{
				  // Simulating without depths.

				  // check of named reference layer
				  layer = massager->outputgrid[i].ECLIPSE_reference_depth_layer_number;
				  if ((layer<=0) || ( layer > Nz))
					{
					  printf("ECLIPSE_reference_depth_layer_number value in XML outside scope of ECLIPSE grid");
					  return EXIT_FAILURE;
					}
				}
			}
		}
	}
  printf("Writing eclipse file in-situ = %d\n",massager->eclInfo.writeInSitu);
  return EXIT_SUCCESS;
}

/* 
   Method that sets blocks to have homogeneous segments up a column with common i,j

*/
int patchEclSegmentLabels(struct MASSAGER *massager,ECL *ecl,int haveSegments)
{
  int i,j,k,Nx,Ny,Nz,isThin;
  int seg,cnt=0,cnt_max,seg_max=0;
  int ***segs= NULL;
  PROP *aPS = NULL; 

  if (!haveSegments) return EXIT_SUCCESS;

  Nx = getNX(ecl);
  Ny = getNY(ecl);
  Nz = getNZ(ecl);  


  /* Get segments */
  aPS = fetchProp(ecl,"SEGMENTS");
  if (!aPS || !aPS->set || !aPS->p) 
	{
	  fprintf(stderr,"activateThinActnums(): Could not get SEGMENTS\n");
	  return EXIT_FAILURE;
	}
  segs = (int ***) aPS->p;

  for (i=1; i<=Nx; i++)
	{
	  for (j=1; j<=Ny; j++)
		{
		  // determine majority segment label
		  seg=cnt_max=-1;
		  for (k=1; k<=Nz; k++)
			{
			  if (segs[i][j][k] >= 0)
				{				  
				  if (seg==-1) 
					{
					  seg = segs[i][j][k];
					  cnt=1;
					  seg_max = seg;
					}
				  else if (segs[i][j][k] == seg) 
					{ 
					  cnt++;
					}
				  else 
					{
					  if (cnt>cnt_max)
						{
						  if (massager->verbosity>1)
							{
							  printf("Warning! non-constant segment label in column (%d,%d) of eclipse file\n",i,j);
							}
						  cnt_max = cnt;
						  seg_max = seg;
						}
					  seg = segs[i][j][k];
					  cnt=1;
					}
				}
			}
		  // apply this uniformly down column if we found a useful seg label
		  if (seg >= 0)
			{
			  for (k=1; k<=Nz; k++) segs[i][j][k] = seg_max;
			}
		}
	}
  return EXIT_SUCCESS;
}



/* 
   Method that sets blocks to active if they seem to be inactive
   because of small thickness.  These are switched on if they belong
   to a column of blocks that has at least one active block.
   Approximate "applyPinch()" - like test only. See applyPinch() in
   ECLIO
*/
int activateThinActnums(ECL *ecl,int haveSegments)
{
  int i,j,k,Nx,Ny,Nz,isThin,eligible;
  int ***active= NULL;
  int ***segs= NULL;
  PROP *aP = NULL; 
  PROP *aPS = NULL; 
  PROP *from=NULL;
  double ****z=NULL;

  double eps = 0.25 * ecl->pinch.pVal;

  Nx = getNX(ecl);
  Ny = getNY(ecl);
  Nz = getNZ(ecl);  

  aP = fetchProp(ecl, "ACTNUM");
  if (!aP || !aP->p || !aP->set)
	{
	  fprintf(stderr,"activateThinActnums(): Could not get ACTNUM\n");
	  return EXIT_FAILURE;
	}
  active = (int ***) aP->p;

  from = fetchProp(ecl,"ZCORN");
  if (!from || !from->p)
	{
	  fprintf(stderr,"activateThinActnums(): Could not get ZCORN\n");
	  return EXIT_FAILURE;
	}
  z = ((double ****)(from->p)); 
  double *zC;

  // don't activate a block unless it has a decent segment label, if available
  if (haveSegments)
	{
	  /* Get segments */
	  aPS = fetchProp(ecl,"SEGMENTS");
	  if (!aPS || !aPS->set || !aPS->p) 
		{
		  fprintf(stderr,"activateThinActnums(): Could not get SEGMENTS\n");
		  return EXIT_FAILURE;
		}
	  segs = (int ***) aPS->p;
	}

  for (i=1; i<=Nx; i++)
	{
	  for (j=1; j<=Ny; j++)
		{
		  // determine if column of blocks is eligible for switching
		  eligible=0;
		  for (k=1; k<=Nz; k++)
			{
			  // determine if column of blocks is eligible for switching
			  if (active[i][j][k])
				{
				  // at least one is active.
				  eligible=1;
				  break;
				}
			}
		  if (eligible)
			{
			  for (k=1; k<=Nz; k++)
				{
				  if (!active[i][j][k])
					{
					  // is it thin ???? Test all corners. 
					  // Assume thin till disproved
					  isThin = 1;
					  zC = z[i][j][k];

					  if ( fabs(zC[1] - zC[5]) > eps) isThin=0;
					  if ( fabs(zC[2] - zC[6]) > eps) isThin=0;
					  if ( fabs(zC[3] - zC[7]) > eps) isThin=0;
					  if ( fabs(zC[4] - zC[8]) > eps) isThin=0;
				  
					  if (isThin)
						{
						  // Switch it on: we suspect it's inactive b/c of thinness
						  // But check if it has a decent segment label first
						  if (haveSegments)
							{
							  // restrict reactivation to labelled blocks
							  if (segs[i][j][k]>=0)
								{
								  active[i][j][k] = 1;
								}
							}
						  else
							{
							  // assume we can blaze away....
							  active[i][j][k] = 1;
							}

						}
					}
				}
			}
		}
	}

  return EXIT_SUCCESS;
}

/* 
   Method that sets blocks to active if they seem to be inactive because of small thickness.
   Approximate "applyPinch()" - like test only. See applyPinch() in ECLIO
*/
void debugPrintActnums(ECL *ecl,char* filename)
{
  int i,j,k,Nx,Ny,Nz;
  int ***active= NULL;
  PROP *aP = NULL; 
  FILE *f;


  Nx = getNX(ecl);
  Ny = getNY(ecl);
  Nz = getNZ(ecl);  

  aP = fetchProp(ecl, "ACTNUM");
  active = (int ***) aP->p;

  printf("Dumping actives to %s\n",filename);
  f = fopen(filename,"w");

  for (i=1; i<=Nx; i++)
	{
	  for (j=1; j<=Ny; j++)
		{
		  fprintf(f,"%d %d %d\n",i,j,active[i][j][2]);
		}
	}
  fclose(f);
}

int cleanUpSundryMemory(struct MASSAGER *massager)
{
  int i;
  struct WELL_OBSERVATIONS *wo;

  for (i=0; i<massager->nLayers; i++)
	{
	  wo = &massager->well_observations[i];
	  free(wo->obs);
	  wo->obs = NULL;
	  free_ivector(wo->n_obs,0,massager->nPropertyNames-1);
	  wo->n_obs = NULL;
	}

  if (massager->propertyNames)
	{
	  for (i=0; i < massager->nPropertyNames; i++) 
		{
		  if (massager->propertyNames[i])
			{
			  free(massager->propertyNames[i]); 
			  massager->propertyNames[i]=NULL; 
			} 
		}
	  free(massager->propertyNames);
	  massager->propertyNames = NULL;
	}

  free(massager->properties); 
  massager->properties = NULL;

  free(massager->layerNumbers); 
  massager->layerNumbers = NULL;

  free(massager->variograms); 
  massager->variograms = NULL;

  free(massager->propertyNames); 
  massager->propertyNames = NULL;

  free(massager->well_observations);
  massager->well_observations = NULL;

  free(massager->massage_analysis_file);
  massager->massage_analysis_file = NULL;

  free(massager->inputgrid.xWord);
  massager->inputgrid.xWord=NULL;
  free(massager->inputgrid.yWord);
  massager->inputgrid.yWord=NULL;

  return EXIT_SUCCESS;
}

int cleanUpGrids(struct MASSAGER *massager)
{
  int i;

  freeMassageAnalysisData(massager);

  for (i=0; i< massager->nLayers; i++)
	{
	  /* output grids */
	  if (!(massager->duplicate))
		{
		  if (!massager->outputgrid[i].isAMirror)
			{
			  if (massager->outputgrid[i].xValues)
				{
				  free_dvector(massager->outputgrid[i].xValues,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].xValues = NULL;
				}
			  if (massager->outputgrid[i].yValues)
				{
				  free_dvector(massager->outputgrid[i].yValues,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].yValues = NULL;
				}
			  if (massager->outputgrid[i].zValues)
				{
				  free_dvector(massager->outputgrid[i].zValues,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].zValues = NULL;
				}
			  if (massager->outputgrid[i].areas)
				{
				  free_vector(massager->outputgrid[i].areas,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].areas = NULL;
				}
			  if (massager->outputgrid[i].ecl_i)
				{
				  free_ivector(massager->outputgrid[i].ecl_i,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].ecl_i = NULL;
				}
			  if (massager->outputgrid[i].ecl_j)
				{
				  free_ivector(massager->outputgrid[i].ecl_j,0,massager->outputgrid[i].nPoints-1);
				  massager->outputgrid[i].ecl_j = NULL;
				}
			  if ((i==0) && massager->eclInfo.commonEclipseGrid)
				{				  
				  if (massager->outputgrid[i].ecl_active != NULL)
					{
					  free_imatrix(massager->outputgrid[i].ecl_active,0,massager->outputgrid[i].ecl_Nx,
								   0,massager->outputgrid[i].ecl_Ny);
					  massager->outputgrid[i].ecl_active = NULL;
					}
				}
			}
		}
	  if (massager->eclInfo.useECLIOsegments)
		{
		  if (massager->outputgrid[i].segnum)
			{
			  free_ivector(massager->outputgrid[i].segnum,0,massager->outputgrid[i].nPoints-1);
			  massager->outputgrid[i].segnum = NULL;
			}
		  if (massager->outputgrid[i].segnum_corner_encoding)
			{
			  free_cvector(massager->outputgrid[i].segnum_corner_encoding,0,massager->outputgrid[i].nPoints-1);
			  massager->outputgrid[i].segnum_corner_encoding = NULL;
			}
		}
	  // Simple eclipse cleanup. Will be much messier with multiple input grids
	  if (!massager->duplicate) 
		{
		  if (i==0 && massager->eclInfo.commonEclipseGrid)
			{
#if DUPLICATE_AT_FAULTS
			  free_segInfo(massager,massager->outputgrid[i].segInfo,
						   massager->outputgrid[i].ecl_Nx,massager->outputgrid[i].ecl_Ny);
			  massager->outputgrid[i].segInfo = NULL;
			  if (massager->outputgrid[i].segcol)
				{
				  free_imatrix(massager->outputgrid[i].segcol, 1, massager->outputgrid[i].ecl_Nx ,1, massager->outputgrid[i].ecl_Ny);
				  massager->outputgrid[i].segcol = NULL;
				}
#endif
			  freeModel(&massager->outputgrid[0].ecl);
			}
		}
	}

  // downwards free for copies
  for (i=massager->nLayers-1; i>=0 ;i--)
	{
	  if (massager->outputgrid[i].layer_name)
		{
		  free(massager->outputgrid[i].layer_name);
		  massager->outputgrid[i].layer_name = NULL;
		}
	  if (massager->outputgrid[i].filename)
		{
		  if (i==0 || (massager->outputgrid[i].filename != massager->outputgrid[0].filename))
			{
			  free(massager->outputgrid[i].filename);
			  massager->outputgrid[i].filename = NULL;
			}
		}
	  if (massager->outputgrid[i].ECLIPSE_segments_keyword)
		{
		  if (i==0 || (massager->outputgrid[i].ECLIPSE_segments_keyword != massager->outputgrid[0].ECLIPSE_segments_keyword))
			{
			  free(massager->outputgrid[i].ECLIPSE_segments_keyword);
			  massager->outputgrid[i].ECLIPSE_segments_keyword = NULL;
			}
		}
	}


  free(massager->outputgrid);
  massager->outputgrid=NULL;



  return EXIT_SUCCESS;
}


int isADesiredECLIPSEKeyword(char* s)
{
  int i;

#ifdef PETREL_HACKS
  // hack around Petrel bug
  if (!strncmp(s,"SEGMENTS",8))
	{
	  for (i=0;i<N_DESIRED_ECLIPSE_KEYWORDS;i++)
		{
		  if (!strncmp(s,DESIRED_ECLIPSE_KEYWORDS[i],8)) return 1;
		}
	}
  else
	{
	  for (i=0;i<N_DESIRED_ECLIPSE_KEYWORDS;i++)
		{
		  if (!strcmp(s,DESIRED_ECLIPSE_KEYWORDS[i])) return 1;
		}
	}
#else
  for (i=0;i<N_DESIRED_ECLIPSE_KEYWORDS;i++)
	{
	  if (!strcmp(s,DESIRED_ECLIPSE_KEYWORDS[i])) return 1;
	}
#endif

  /* if here, no match */
  return 0;
}

int isABlockCentredECLIPSEKeyword(char* s)
{
  int i;

  // hack around Petrel bug
  if (!strncmp(s,"thickness",9))
	{
	  return 0;
	}
  else
	{
	  return 1;
	}
}


int writeEclipseCache(struct MASSAGER *massager,int i)
{
  char *fncache;
  FILE *fp;
  char	line_buffer[LINE_LENGTH];
  struct stat file_stat;
  time_t ecl_time=0;
  long int cache_time=0;
  int have_segments=0;
  int active_size;

  if (massager->eclInfo.commonEclipseGrid)
	{
	  // cache for multiple layers not yet working...
	  return 0;
	}
	
  asprintf(&fncache,".%s.%d.cache",massager->outputgrid[i].filename, massager->outputgrid[i].ECLIPSE_layer_number);
  printf("Writing file %s\n",fncache);
  fp = fopen(fncache,"w");
  if (fp==NULL) 
	{
	  printf("Can't write cache!\n");
	  free(fncache);
	  return EXIT_FAILURE;
	}
  else
	{
	  // get filestat date for full eclipse file...
	  stat(massager->outputgrid[i].filename, &file_stat);
	  ecl_time = file_stat.st_mtime;
	  
	  // Compare to data from cache
	  active_size = massager->outputgrid[i].nPoints;
	  fprintf(fp,"%s_date: %ld\n",massager->outputgrid[i].filename,ecl_time);
	  fprintf(fp,"active_size: %d\n",active_size);
	  fprintf(fp,"have_segments: %d\n", (massager->outputgrid[i].segnum == NULL) ? 0:1);

	  fwrite(&massager->outputgrid[i].xValues[0], sizeof(double), active_size, fp);
	  fwrite(&massager->outputgrid[i].yValues[0], sizeof(double), active_size, fp);
	  fwrite(&massager->outputgrid[i].zValues[0], sizeof(double), active_size, fp);

	  fwrite(&massager->outputgrid[i].ecl_i[0], sizeof(int), active_size, fp);
	  fwrite(&massager->outputgrid[i].ecl_j[0], sizeof(int), active_size, fp);

	  fwrite(&massager->outputgrid[i].areas[0], sizeof(float), active_size, fp);

	  if (massager->outputgrid[i].segnum != NULL)
		{
		  fwrite(&massager->outputgrid[i].segnum[0], sizeof(int), active_size, fp);
		}
	  fclose(fp);		
	}
  printf("Finished cache write\n");
  free(fncache);
  return 1;
}


int readEclipseCache(struct MASSAGER *massager,int i)
{
  char *fncache;
  FILE *fp;
  char	line_buffer[LINE_LENGTH];
  struct stat file_stat;
  time_t ecl_time=0;
  long int cache_time=0;
  int have_segments=0;
  int active_size=0;

  // caching doesn't yet work for full ECLIPSE output....
  if (massager->eclInfo.writeECLOut) return 0;

  asprintf(&fncache,".%s.%d.cache", massager->outputgrid[i].filename, massager->outputgrid[i].ECLIPSE_layer_number);

  printf("Reading cache file %s\n",fncache);
  fp = fopen(fncache,"r");
  free(fncache);
  if (fp==NULL) 
	{
	  return 0;
	}
  else
	{
	  // get filestat date for full eclipse file...
	  stat(massager->outputgrid[i].filename, &file_stat);
	  ecl_time = file_stat.st_mtime;
	  
	  // Compare to data from cache
	  if ((fscanf(fp,"%s %ld\n",line_buffer,&cache_time)!=2) || (cache_time != ecl_time)) return 0;
	  else
		{
		  // cache date is OK. Read sizes
		  if (fscanf(fp,"%s %d\n",line_buffer,&active_size)!=2) return 0;
		  if (fscanf(fp,"%s %d\n",line_buffer,&have_segments)!=2) return 0;
		  printf("Cache ecl read: Active grid size = %d\n",active_size);

		  /* alloc and load */ 
		  massager->outputgrid[i].isAMirror=0;
		  massager->outputgrid[i].isRectangular=0; 
		  massager->outputgrid[i].nPoints = active_size;
		  massager->outputgrid[i].xValues = dvector(0,active_size-1);
		  massager->outputgrid[i].yValues =  dvector(0,active_size-1);
		  massager->outputgrid[i].zValues = NULL;
		  massager->outputgrid[i].zValues =  dvector(0,active_size-1);
		  massager->outputgrid[i].areas = NULL;
		  massager->outputgrid[i].areas =  vector(0,active_size-1);
		  massager->outputgrid[i].ecl_i = NULL;
		  massager->outputgrid[i].ecl_i =  ivector(0,active_size-1);
		  massager->outputgrid[i].ecl_j = NULL;
		  massager->outputgrid[i].ecl_j =  ivector(0,active_size-1);
		  massager->totalNOutputGridPoints += massager->outputgrid[i].nPoints; // sum over all output grids

		  /* segments */
		  if (have_segments)
			{
			  massager->outputgrid[i].segnum = ivector(0,active_size-1);
			}
		  else massager->outputgrid[i].segnum = NULL;

		  /* midpoints of traces: 6 values are x,y,z_top,x,y,z_bottom */
		  fread(&massager->outputgrid[i].xValues[0], sizeof(double), active_size, fp);
		  fread(&massager->outputgrid[i].yValues[0], sizeof(double), active_size, fp);
		  fread(&massager->outputgrid[i].zValues[0], sizeof(double), active_size, fp);

		  fread(&massager->outputgrid[i].ecl_i[0], sizeof(int), active_size, fp);
		  fread(&massager->outputgrid[i].ecl_j[0], sizeof(int), active_size, fp);

		  fread(&massager->outputgrid[i].areas[0], sizeof(float), active_size, fp);

		  if (have_segments)
			{
			  fread(&massager->outputgrid[i].segnum[0], sizeof(int), active_size, fp);
			}
		  fclose(fp);
		}
	}
  return 1;
}
/*
  Writes ECLIO in two modes: 
  accumulateMode = 1 for kriged means and
                   stddevs. Will be called for each layer 
				   - then writes on final layer
  accumulateMode = 0 for realisations. Called once with 
                   everything in place in means.  (real_no 
				   only used for debugging)
*/
int writeECLIO_InSituOutput(int real_no, int accumulateMode, int i_layer,struct MASSAGER *massager,
							float*** means,float*** stddevs)
{
  int i,j,iB,jB,k,i_dup,ii,m,mm,m_prop;
  int ecl_layer;
  int i_layer_st,i_layer_end;
  int errCount=0;
  int Nx,Ny,Nz,cnt,status;
  double sum;
  char* name;
  char c_cmp;
  char mightBeRefSurface;
  struct ECL_OUT_INFO *eclOutInfo;
  int ecl_ref_layer=0;
  int** active;
  int finalCall;
  struct SEG_ENCODING **segInfo;
  char* filename;
  PROP *from=NULL;
  double ****z=NULL;
  FILE* fp;
  

  /* useful references */
  eclOutInfo = &(massager->outputgrid[0].eclOutInfo);


  if (accumulateMode)
	{
	  filename = eclOutInfo->outputFilename;
	  fp = eclOutInfo->ecl_out_fp;
	  i_layer_st = i_layer;
	  i_layer_end = i_layer;
	}
  else
	{
	  filename = eclOutInfo->outputFilename_r;
	  fp = eclOutInfo->ecl_out_r_fp;
	  i_layer_st = 0;
	  i_layer_end = massager->nLayers-1;
	}

  if ((!accumulateMode || (i_layer == massager->nLayers-1)) && massager->verbosity>2) 
	{
	  printf("Writing ECLIPSE insitu file: %s\n",filename);
	}

  /* Eclipse dimensions */
  Nx = getNX(&eclOutInfo->eclOut);
  Ny = getNY(&eclOutInfo->eclOut);
  Nz = getNZ(&eclOutInfo->eclOut);

  int*** cnt_lookup = i3tensor(0,Nx,0,Ny,0,3);  // up to 4 possible duplications of x,y

  show_malloc_stats(); 



  for (i_layer = i_layer_st; i_layer<=i_layer_end; i_layer++)
	{
	  finalCall = ((i_layer == massager->nLayers-1)? 1 : 0);

	  //  i,j index table lookup to kriging solutions
	  // Segment info
	  segInfo = (massager->outputgrid[i_layer]).segInfo;    // common segInfo if common source grid
	  active = eclOutInfo->active;
	  cnt=0;
	  for (i=0;i<=Nx;i++)
		{
		  for (j=0;j<=Ny;j++)
			{
			  // actives are common to all grids
			  if (active[i][j]) 
				{
#if DUPLICATE_AT_FAULTS
				  if (massager->eclInfo.useECLIOsegments)
					{
					  for (i_dup=0; i_dup<segInfo[i][j].size; i_dup++)
						{
						  cnt_lookup[i][j][i_dup] = cnt;
						  cnt++;
						}
					}
				  else
					{
					  // presumably no duplicates: 
					  cnt_lookup[i][j][0] = cnt;
					  cnt++;
					}
#else
				  // no duplicates: 
				  cnt_lookup[i][j][0] = cnt;
				  cnt++;
#endif
				} 
			  else cnt_lookup[i][j][0] = -1;
			}
		}


	  // Geometry (zCorn) adjustment from thickness OR depth: perform on final pass
	  if (finalCall)
		{  
		  if (massager->thickness_index >= 0 || massager->d_index >= 0)
			{
			  int kB,iBB,jBB,kBB,c;
			  double dZ=0,Z,Zref,Zref_ext=0;
			  int massager_layer_indx=-1, massager_depth_layer_indx=-1;
			  int haveDepth = 0;

			  ecl_ref_layer = massager->outputgrid[i_layer].ECLIPSE_reference_depth_layer_number;

			  m = massager->thickness_index; 
			  if (massager->d_index >= 0) haveDepth=1;

			  // set up master index table for thicknesses, if available
			  int* layer_index_table = ivector(1,Nz);
			  for (kB=1;kB<=Nz; kB++) layer_index_table[kB] = -1;
			  if (massager->thickness_index >= 0)
				{
				  for (ii=0; ii<massager->nLayers; ii++)
					{
					  layer_index_table[massager->outputgrid[ii].ECLIPSE_layer_number] = ii;
					}
				}
			  if (massager->d_index >= 0)
				{
				  for (ii=0; ii<massager->nLayers; ii++)
					{
					  if (massager->outputgrid[ii].ECLIPSE_layer_number == ecl_ref_layer)
						{
						  massager_depth_layer_indx = ii;
						}
					}
				}

			  // Load corner points...

			  if (!accumulateMode && massager->eclInfo.writeInSitu)
				{
				  // A "realisation grid": restore zCorns from original grid: they get nuked each realisation
				  if (copyPropECL(&massager->outputgrid[0].eclOutInfo.eclOut, "ZCORN", 
								  &massager->outputgrid[0].ecl, "ZCORN", 0)) 
					{
					  fprintf(stderr, "writeECLIO_InSituOutput(): Could not copy ZCORN\n");
					  return EXIT_FAILURE; 
					}
				}

			  from = fetchProp(&eclOutInfo->eclOut,"ZCORN");
			  z = ((double ****)(from->p));  // should be eclOut->zCorn

			  // Preliminary check over zCorn values 
			  // to see if thicknesses are consistent...
#if 0
			  for (i=0,iB=1; i<=Nx; i++,iB++)
				{
				  for (j=0,jB=1; j<=Ny; j++,jB++) 
					{
					  if (iB<=Nx && jB<=Ny)
						{					  
						  for (kB=1; kB<=Nz ; kB++)
							{
							  if (iB>=1 && jB>=1 && iB<=Nx && jB<=Ny)
								dZ = z[iB][jB][kB][5] - z[iB][jB][kB][1];

							  // loop on adjacent columns
							  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
								{
								  for (iBB=iB; iBB>=iB-1; iBB--,c++)
									{
									  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny)
										{
										  int bug=0;
										  // consistent thicknesses?
										  bug = !APPROX(z[iBB][jBB][kB][c+4] - z[iBB][jBB][kB][c], dZ, 1e-6);
										  if (bug)
											{ 
											  printf("PRE ZCORN: discontinuous thicknesses at ecl_i,ecl_j = (%d,%d)\n",iB,jB);
#if DUPLICATE_AT_FAULTS
											  if (segInfo[i][j].size==1)
												{
												  printf("MISMATCHED THICKNESS IN CTS AREA %f %f\n",z[iBB][jBB][kB][c+4] - z[iBB][jBB][kB][c], dZ);
												}
#endif
											}
										}
									}
								}
							}
						} // active column
					}  // raster scan j
				} // raster scan i
#endif


			  // ...............Adjustment of grid..............................

			  // loop over grid edges (trace locations)
			  for (i=0,iB=1; i<=Nx; i++,iB++)
				{
				  for (j=0,jB=1; j<=Ny; j++,jB++) 
					{
					  if (active[i][j]) 
						{
						  for (i_dup=0; i_dup<segInfo[i][j].size; i_dup++)
							{
							  // think: thickness[i][j][massager_layer_indx] = krigMeans[massager_layer_indx][cnt_lookup[i][j][i_dup]][m];
							  // move corners at this trace relative to reference layer
							  // propagate x,y,z zCorns from reference layer

							  // seg number
							  // int seg = segInfo[i][j].values[i_dup];

							  // corner 'c' encoding for this (possibly duplicate) point
							  c_cmp = segInfo[i][j].corner_encoding[i_dup];

							  // Depth of reference surface if available
							  if (haveDepth && (massager_depth_layer_indx != -1)
								  && (cnt_lookup[iB-1][jB-1][i_dup] != -1))
								{
								  // Setting top of block. This may be the depth-sim surface if ref. layer
								  // depth of reference layer was kriged/simulated etc	
								  // Use this depth if it's available
								  Zref_ext = means[massager_depth_layer_indx][cnt_lookup[iB-1][jB-1][i_dup]][massager->d_index];
								}

							  // Going down from ref layer...
							  for (kB=ecl_ref_layer; kB<=Nz ; kB++)
								{
								  if (layer_index_table[kB] != -1)
									{
									  // propagate x,y,z zCorns) from reference layer						  
									  massager_layer_indx = layer_index_table[kB];
									  dZ = means[massager_layer_indx][cnt_lookup[i][j][i_dup]][m];
									}
								  else 
									{
									  // layer not massaged: use ecl thickness
									  // at i,j edges get from zCorn carefully. First case is most usual.
									  if (i<Nx && j<Ny)
										{
										  // most common (interior) case
										  dZ = z[iB][jB][kB][5] - z[iB][jB][kB][1];
										}
									  else if (i<Nx && j==Ny)
										{
										  // j-edge gridblock thickness
										  dZ = z[iB][jB-1][kB][7] - z[iB][jB-1][kB][3];
										}
									  else if (i==Nx && j<Ny)
										{
										  // i-edge gridblock thickness
										  dZ = z[iB-1][jB][kB][6] - z[iB-1][jB][kB][2];
										}
									  else if (i==Nx && j==Ny)
										{
										  // corner point. i-edge & j-edge gridblock thickness
										  dZ = z[iB-1][jB-1][kB][8] - z[iB-1][jB-1][kB][4];
										}
									}
								  // Adjust bottom zCorn; This loops around gridblocks 
								  // in order  2 1   with axis sense |j
								  //           4 3                   |__i
								  // 
								  if (dZ<0 && errCount<MAXWARNINGS) 
									{
									  printf("Thickness simulation error 1: check original eclipse grid for sane thicknesses!\n");
									  errCount++;
									  if (errCount==MAXWARNINGS) printf("Subsequent warnings suppressed...\n");
									}

								  // Adjust top of to match layer above
								  mightBeRefSurface = (haveDepth && (kB == ecl_ref_layer) && (massager_depth_layer_indx != -1)) ? 1 : 0;
								
								  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
									{
									  for (iBB=iB; iBB>=iB-1; iBB--,c++)
										{
										  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny && (c_cmp & (1<<c)))
											{
											  // top of block = bottom of block above
											  if (mightBeRefSurface  && (cnt_lookup[iB-1][jB-1][i_dup] != -1))
												{ 
												  // Setting top of block. This may be the depth-sim surface if ref. layer
												  // depth of reference layer was kriged/simulated etc	
												  // Use this depth if it's available
												  // Zref = krigMeans[massager_depth_layer_indx][cnt_lookup[iB-1][jB-1][i_dup]][massager->d_index];
												  // z[iBB][jBB][kB-1][c+4];
												  Zref = Zref_ext;
												  z[iBB][jBB][kB][c] = Zref;
												}
											  else if (kB-1 >= 1)
												{
												  // top of block = bottom of block above
												  // bottom of block above
												  Zref = z[iBB][jBB][kB-1][c+4];
												  z[iBB][jBB][kB][c] = Zref;
												}
											}
										}
									}

								  // set bottom
								  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
									{
									  for (iBB=iB; iBB>=iB-1; iBB--,c++)
										{
										  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny && (c_cmp & (1<<c)))
											{
											  // bottom = top + thickness
											  z[iBB][jBB][kB][c+4] = z[iBB][jBB][kB][c] + dZ;
											}
										}
									}
					  
								}
							  // Going UP from ref layer...
							  for (kB=ecl_ref_layer-1; kB>=1 ; kB--)
								{
								  // ditto, going up
								  if (layer_index_table[kB] != -1)
									{
									  // propogate x,y,z zCorns) from reference layer						  
									  massager_layer_indx = layer_index_table[kB];
									  dZ = means[massager_layer_indx][cnt_lookup[i][j][i_dup]][m];
									}
								  else
									{
									  // layer not massaged: use ecl thickness
									  // at i,j edges get from zCorn carefully. First case is most usual.
									  if (i<Nx && j<Ny)
										{
										  // most common (interior) case
										  dZ = z[iB][jB][kB][5] - z[iB][jB][kB][1];
										}
									  else if (i<Nx && j==Ny)
										{
										  // j-edge gridblock thickness
										  dZ = z[iB][jB-1][kB][7] - z[iB][jB-1][kB][3];
										}
									  else if (i==Nx && j<Ny)
										{
										  // i-edge gridblock thickness
										  dZ = z[iB-1][jB][kB][6] - z[iB-1][jB][kB][2];
										}
									  else if (i==Nx && j==Ny)
										{
										  // corner point. i-edge & j-edge gridblock thickness
										  dZ = z[iB-1][jB-1][kB][8] - z[iB-1][jB-1][kB][4];
										}
									}
								  if (dZ<0 && errCount<MAXWARNINGS) 
									{
									  printf("Thickness simulation error 2: check original eclipse grid for sane thicknesses!\n");
									  errCount++;
									  if (errCount==MAXWARNINGS) printf("Subsequent warnings suppressed...\n");
									}

								  // Adjust bottom to match layer below
								  if (kB+1 <= Nz) // always true
									{
									  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
										{
										  for (iBB=iB; iBB>=iB-1; iBB--,c++)
											{
											  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny && (c_cmp & (1<<c)))
												{
												  // bottom of layer = top of layer below
												  z[iBB][jBB][kB][c+4] = z[iBB][jBB][kB+1][c];
												}
											}
										}
									}
								  // Adjust top zCorn values......
								  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
									{
									  for (iBB=iB; iBB>=iB-1; iBB--,c++)
										{
										  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny && (c_cmp & (1<<c)))
											{
											  // top = bottom - thickness
											  z[iBB][jBB][kB][c] = z[iBB][jBB][kB][c+4] - dZ;
											}
										}
									}
					  
								} // loop over kB
							} // active column			  			  
						} // duplictes
					}  // raster scan j
				} // raster scan i

			  // CHECK! Post-distorted gridblocks are attached....
#if 0
			  for (i=0,iB=1; i<=Nx; i++,iB++)
				{
				  for (j=0,jB=1; j<=Ny; j++,jB++) 
					{
					  if (iB<=Nx && jB<=Ny)
						{					  
						  for (kB=1; kB<=Nz ; kB++)
							{
							  if (iB>=1 && jB>=1 && iB<=Nx && jB<=Ny)
								dZ = z[iB][jB][kB][5] - z[iB][jB][kB][1];

							  // loop on adjacent columns
							  for (jBB=jB,c=1; jBB>=jB-1; jBB--)
								{
								  for (iBB=iB; iBB>=iB-1; iBB--,c++)
									{
									  if (iBB>=1 && jBB>=1 && iBB<=Nx && jBB<=Ny)
										{
										  int bug=0;
										  // consistent thicknesses?
										  bug = !APPROX(z[iBB][jBB][kB][c+4] - z[iBB][jBB][kB][c], dZ, 1e-6);
										  if (bug)
											{ 
											  printf("POST ZCORN: discontinuous thicknesses at ecl_i,ecl_j = (%d,%d) <-> (%d,%d)\n",iB,jB,iBB,jBB);
#if DUPLICATE_AT_FAULTS
											  if (segInfo[i][j].size==1)
												{
												  printf("MISMATCHED THICKNESS IN CTS AREA %f %f\n",z[iBB][jBB][kB][c+4] - z[iBB][jBB][kB][c], dZ);
												}
#endif
											}
										}
									}
								}
							}
						} // active column
					}  // raster scan j
				} // raster scan i
#endif




			  free_ivector(layer_index_table,1,Nz);
			  layer_index_table=NULL;
			}
		  if (0)
			{
			  if (z==NULL)
				{
				  from = fetchProp(&eclOutInfo->eclOut,"ZCORN");
				  z = ((double ****)(from->p));  // should be eclOut->zCorn
				} 
			  checkXsection(filename,10,Ny,Nz,z);
			  checkYsection(filename,50,Nx,Nz,z);
			}
		}



	  /*************       Block centred properties   ***********/

	  // Treatment of duplicates: all properties except depth are assumed continuous across faults
	  // => use only cnt_lookup[i][j][0] nodes

	  // ptr to storage for contiguous cube array needed for eclio writing
	  float*** tmpProp;
	  float **workProp = NULL;
	  if (accumulateMode)
		{
		  // convert surfaces to properties...
		  ecl_layer = massager->outputgrid[i_layer].ECLIPSE_layer_number-1;


		  // loop over {property, STDDEV_property} pairs. Associated field is
		  // eclOutInfo->keyWords[eclOutInfo->nDuplicateKeys + mm]
		  for (mm=0,m_prop=0; (mm<2*massager->nPropertyNames); mm++)
			{
			  m = (int)(mm/2);
 
			  //if (m == massager->thickness_index || m == massager->d_index ) 
			  if (!strcmp(massager->propertyNames[m],"d") && ((mm&1) == 0))  // odd numbered props are stddevs
				{		  
				  // skip depth in rest of loop (not a body centred property)
				  // test above will accept STDDEV_d (which is the next mm)
				  continue;
				}

			  if (!workProp) workProp = matrix(0,Nx,0,Ny);

			  // copy corner properties to workspace array
			  for (i=0;i<=Nx;i++)
				{
				  for (j=0;j<=Ny;j++) 
					{
					  if (active[i][j]) 
						{
						  // massager->outputgrid[i_layer].eclOutInfo.property[m][i][j] = krigMeans[i_layer][cnt_lookup[i][j][0]][m];
						  // massager->outputgrid[i_layer].eclOutInfo.property[m+1][i][j] = krigStdDev[i_layer][cnt_lookup[i][j][0]][m+1];
						  if ((mm % 2)==0) //even mm's
							{
							  // a mean
							  // tmpProp[i][j][ecl_layer] = means[i_layer][cnt_lookup[i][j][0]][m];
							  workProp[i][j] = means[i_layer][cnt_lookup[i][j][0]][m];
							}
						  else
							{
							  // a stddev
							  // tmpProp[i][j][ecl_layer] = stddevs[i_layer][cnt_lookup[i][j][0]][m];
							  workProp[i][j] = stddevs[i_layer][cnt_lookup[i][j][0]][m];
							}
						}
					}
				}

			  // pointer to right page of big accumulating 4-cube
			  tmpProp = massager->eclInfo.property[m_prop];

			  // Basic Arithmetic averaging....
			  for (i=0;i<Nx;i++)
				{
				  for (j=0;j<Ny;j++) 
					{
					  // average 4 corners
					  for (iB=i,k=0,sum=0.0; iB<=i+1; iB++)
						{
						  for (jB=j; jB<=j+1; jB++)							
							{
							  if (active[iB][jB]) 
								{
								  sum += workProp[iB][jB];
								  k++;
								}
							}
						}
					  tmpProp[i][j][ecl_layer] = (k>0) ? sum/(float)k : D_MISSING;
					}
				}


			  /* Use property name as ECLIPSE keyword */
			  /* copy 3d array to a PROP */
			  // wait till last layer till fill prop and write
			  if (finalCall)
				{
				  name = eclOutInfo->keyWords[m_prop+1+eclOutInfo->nDuplicateKeys];
				  if (copyArrayToProp(&eclOutInfo->eclOut, (PROP *) NULL,
									  name,
									  (GENERIC_POINTER) tmpProp,
									  3, eclOutInfo->bounds, 'f', 1) )  
					/* last 1 makes sure the data fit (Block-centred) */ 
					{
					  fprintf(stderr,"Could not copy %s.\n",name);
					  return EXIT_FAILURE; 
					}
				}
			  m_prop++;
			}
		  if (workProp) 
			{
			  free_matrix(workProp,0,Nx,0,Ny); 
			  workProp=NULL;
			}

		  if (finalCall)
			{

			  /* 3D grid is finally populated. Write it out */
			  status =  writeRecordList(fp , 
										&eclOutInfo->eclOut, 
										eclOutInfo->nKeys, 
										eclOutInfo->keyWords);
			  if (status)  
				{
				  fprintf(stderr, "Error encountered in writeRecordList: status = %d.\n", status);
				  return EXIT_FAILURE; 
				}
			  /* and close I/O */
			  fclose(fp);

			  // debugPrintActnums(&eclOutInfo->eclOut,"act2.txt"); exit(0);
			  
			  if (massager->Nrealisations == 0)
				{ 
				  // don't need anything for later. Clean up now.
				  free_imatrix(eclOutInfo->bounds,1,3,1,2);
				  eclOutInfo->bounds = NULL;
				  if (eclOutInfo->nGridBlockKeys > 0)
					{
					  free_f4tensor(massager->eclInfo.property,0,eclOutInfo->nGridBlockKeys-1, 0,Nx-1, 0,Ny-1, 0,Nz-1);
					  massager->eclInfo.property=NULL;
					}
				  /* now get rid of eclOut data */
				  freeModel(&eclOutInfo->eclOut);
				}
			  if ((massager->verbosity>2)) 
				{
				  printf("Wrote ECLIPSE file %s\n",filename);
				}
			}
		} // ...accumulate mode
	  else  
		{
		  // Realisations. Convert surfaces to properties...
		  ecl_layer = massager->outputgrid[i_layer].ECLIPSE_layer_number-1;
		  for (m=0,m_prop=0; m<massager->nPropertyNames; m++)
			{
			  if (!strcmp(massager->propertyNames[m],"d"))
				{		  
				  // skip depth in rest of loop (not required)
				  continue;
				}

			  name = eclOutInfo->keyWords_r[m_prop+1+eclOutInfo->nDuplicateKeys];

			  printf("preparing prop %s, layer %d, for ECLIPSE file\n",massager->propertyNames[m],
					 massager->outputgrid[i_layer].ECLIPSE_layer_number);

			  if (!workProp) workProp = matrix(0,Nx,0,Ny);

			 
			  // copy corner properties to workspace array
			  for (i=0;i<=Nx;i++)
				{
				  for (j=0;j<=Ny;j++) 
					{
					  if (active[i][j]) 
						{
						  // a simulation
						  // tmpProp[i][j][ecl_layer] = means[i_layer][cnt_lookup[i][j][0]][m];
						  workProp[i][j] = means[i_layer][cnt_lookup[i][j][0]][m];							
						}
					}
				}	

			  // pointer to right page of big accumulating 4-cube
			  tmpProp = massager->eclInfo.property[m_prop];
		  
			  // Basic Arithmetic averaging....
			  for (i=0;i<Nx;i++)
				{
				  for (j=0;j<Ny;j++) 
					{
					  // average 4 corners
					  for (iB=i,k=0,sum=0.0; iB<=i+1; iB++)
						{
						  for (jB=j; jB<=j+1; jB++)							
							{
							  if (active[iB][jB]) 
								{
								  sum += workProp[iB][jB];
								  k++;
								}
							}
						}
					  tmpProp[i][j][ecl_layer] = (k>0) ? sum/(float)k : D_MISSING;
					}
				} 
			  
			  // On final layer, copy into eclipse
			  if (finalCall)
				{
				  printf("Copying props to ECLIPSE file\n");
				  /* Copy into eclipse */
				  /* copy 3d array to a PROP */
				  // wait till last layer till fill prop and write
				  if (copyArrayToProp(&eclOutInfo->eclOut, (PROP *) NULL,
									  name,
									  (GENERIC_POINTER) tmpProp,
									  3, eclOutInfo->bounds, 'f', 1) )  
					/* last 1 makes sure the data fit (Block-centred) */ 
					{
					  fprintf(stderr,"Could not copy %s.\n",name);
					  return EXIT_FAILURE; 
					}
				}
			  m_prop++;
			}
		  if (workProp) 
			{
			  free_matrix(workProp,0,Nx,0,Ny); 
			  workProp=NULL;
			}
		  
		  if (finalCall)
			{
			  /* 3D grid is finally populated. Write it out */
			  status =  writeRecordList(fp , 
										&eclOutInfo->eclOut, 
										eclOutInfo->nKeys_r, 
										eclOutInfo->keyWords_r);
			  if (status)  
				{
				  fprintf(stderr, "Error encountered in writeRecordList: status = %d.\n", status);
				  return EXIT_FAILURE; 
				}
		  
			  if ((massager->verbosity>2)) 
				{
				  printf("Wrote ECLIPSE file %s\n",filename);
				}
			}
		}
	}

  /* clean up */
  free_i3tensor(cnt_lookup,0,Nx-1,0,Ny-1,0,3);
  cnt_lookup = NULL;

  return EXIT_SUCCESS;
}

void checkXsection(char *suffix,int iXsection,int nch, int ndh,double ****zCorn)
{
  int i,j,iB,jB,kB;
  FILE *fpX;
  char *filename;

  printf("Writing const ecl_i Xsection at ecl_i=%d\n",iXsection);

  asprintf(&filename,"GridXsection-x.%s",suffix);
  fpX = fopen(filename,"w");
  // simple Xsection check: pipe 
  // stdout | awk '/XSECEND/{p=0}p==1{print $0}/XSECSTART/{p=1}' | gr -C --bitmap-size 1500x500 -B &
  iB = iXsection;

  // for drawing simplicity, draw like this   CABD:
  //        A --B
  //        |   |
  //        C---D

  for (j=0,jB=1; j<nch; j++,jB++) 
	{
	  for (kB=1; kB<=ndh ; kB++)
		{
		  fprintf(fpX,"%d %f\n",j,zCorn[iB][jB][kB][5]);
		  fprintf(fpX,"%d %f\n",j,zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%d %f\n",j+1,zCorn[iB][jB][kB][3]);
		  fprintf(fpX,"%d %f\n",j+1,zCorn[iB][jB][kB][7]);
		  fprintf(fpX,"%d %f\n",j,zCorn[iB][jB][kB][5]);
		}
	  // go across to right edge...
	  fprintf(fpX,"%d %f\n",j+1,zCorn[iB][jB][ndh][7]);
	}
  fclose(fpX);
}

void checkYsection(char *suffix,int jXsection,int nrh, int ndh,double ****zCorn)
{
  int i,iB,jB,kB;
  FILE *fpX;
  char *filename;

  printf("Writing const ecl_j Ysecton at ecl_j=%d\n",jXsection);
  asprintf(&filename,"GridXsection-y.%s",suffix);
  fpX = fopen(filename,"w");
  // simple Xsection check: pipe 
  // stdout | awk '/XSECEND/{p=0}p==1{print $0}/XSECSTART/{p=1}' | gr -C --bitmap-size 1500x500 -B &
  jB = jXsection;

  for (i=0,iB=1; i<nrh; i++,iB++) 
	{
	  for (kB=1; kB<=ndh ; kB++)
		{
		  fprintf(fpX,"%d %f\n",i,zCorn[iB][jB][kB][5]);
		  fprintf(fpX,"%d %f\n",i,zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%d %f\n",i+1,zCorn[iB][jB][kB][2]);
		  fprintf(fpX,"%d %f\n",i+1,zCorn[iB][jB][kB][6]);
		  fprintf(fpX,"%d %f\n",i,zCorn[iB][jB][kB][5]);
		}
	  // go across to right edge...
	  fprintf(fpX,"%d %f\n",i+1,zCorn[iB][jB][ndh][6]);
	}
  fclose(fpX);
}


// read massage-analysed data put out by InversionAnalyser --massage-analyse
int readMassageAnalysisData(struct MASSAGER *massager)
{
    char           *filename;
	double** data_ptr;
    FILE           *in;
	char *ptr;
	char	line_buffer[LINE_LENGTH];
	int i,j,k,line_no,nMatches=0,number_of_columns;
	int layerNo,found,skipLine,filterColumnNumber=-1;
	int nNonCovarianceColumns;
	int n_cov_entries=0,n_prop_entries=0,number_of_coord_columns=0,offset;
	double f;
	int** cov_column;
	int* cov_col_offset;
	int* cov_row_offset;
	fpos_t pos;

	massager->diagonalCovariancesOnly = (massager->Nrealisations==0)? 1 : 0;

	// binary extraction
	int is_binary=0;
	int machine_is_little_endian = 1;
	int short_scale = 32767;  // java standard MAX_SHORT
	int compressed_bytes_per_float = 4; // overwritten later

	if (massager->verbosity>1)
	  {
		printf("Reading Massage-analysis data from %s\n",massager->massage_analysis_file);
		printf("Diagonal covariance storage = %d\n",massager->diagonalCovariancesOnly);
	  }

    filename = massager->massage_analysis_file;
    if ((in = fopen(filename, "r")) == NULL) 
	  {
		warning("cannot open file %s",filename);

        return EXIT_FAILURE;
	  }
	// Parse:
	// title
	fgets(line_buffer,LINE_LENGTH,in);
	// binary test
	if (!strncmp(line_buffer,"Binary-compressed",17))
	  {
		is_binary = 1;
		machine_is_little_endian = !is_big_endian(); 
		compressed_bytes_per_float = 2;
		if (massager->verbosity>1) printf("Reading binary format\n");	
	  }
	else if (massager->verbosity>1) printf("Reading ascii format\n");	

	// number of columns
	if (fscanf(in,"%d\n",&number_of_columns) != 1)
	  {
		warning("Can't read no. of columns in file %s\n",filename);
		return EXIT_FAILURE;
	  }
	if (number_of_columns < (massager->nProperties+2))
	  {
		warning("Too few properties available in file %s\n",filename);
		return EXIT_FAILURE;
	  }

	/* 
	   Count columns titles and compute number of properties in file, plus number of
	   covariance entries 
	*/
	if (fgetpos(in,&pos)) warning("Can't position file pointer");
	for (i=0;i<number_of_columns;i++)
	  {
		if (fscanf(in,"%s",line_buffer)!=1)
		  {
			warning("Can't read column title %s in file %s\n",line_buffer,filename);
			return EXIT_FAILURE;
		  }
		// count entries that look like "C_xx"
		if (!strncmp(line_buffer,"C_",2)) n_cov_entries++;
	  }
	// since the file has the upper triangle of the covariance
	// we must have n_cov_entries = n_prop_entries*(n_prop_entries+1)/2
	// Compute n_prop_entries.
	for (n_prop_entries=0; (n_prop_entries*(n_prop_entries+1))/2 < n_cov_entries; n_prop_entries++);

	// number of columns for coordinates
	number_of_coord_columns = number_of_columns - n_prop_entries - n_cov_entries;

	// here's a brain-dead but simple way to remember the column numbers of the properties
	// in the file *in the natural ordering of the Summary-file*
	cov_column = imatrix(0,n_prop_entries-1,0,n_prop_entries-1);
	offset = number_of_columns - n_cov_entries;
	for (i=k=0; i<n_prop_entries; i++)
	  {
		for (j=i; j<n_prop_entries; j++) 
		  {
			cov_column[i][j] = offset+k;  // the zero-offset column number
			cov_column[j][i] = cov_column[i][j];  // transpose-proof
			k++;
		  }
	  }

	// rewind
	if (fsetpos(in,&pos)) warning("Can't reposition file pointer");

	// Set up check for a property name and layer number match.
	for (j=0; j<(massager->nProperties); j++) massager->properties[j].found=0;

	// column names
	offset = number_of_columns - n_cov_entries - n_prop_entries;
	for (i=0;i<number_of_columns;i++)
	  {
		if (fscanf(in,"%s",line_buffer)!=1)
		  {
			warning("Can't read column title %s in file %s\n",line_buffer,filename);
			return EXIT_FAILURE;
		  }

		// Check for a property name and layer number match. Compare only n characters to skip the _N suffixes if they are there
		for (j=0; j<(massager->nProperties); j++)
		  {
			massager->properties[j].applyFilter=0;
			if (!strncmp(line_buffer,massager->properties[j].name, strlen(massager->properties[j].name)))
			  {
				// name matches
				// extract layer number. We have prop_number, so point ahead and hope
				ptr = &line_buffer[1+strlen(massager->properties[j].name)];
				layerNo = atoi(ptr);
				// (1) check it's a valid number
				for (k=found=0; k<massager->nLayers; k++) 
				  {
					if (layerNo == massager->layerNumbers[k])
					  {
						found=1;
						break;
					  }
				  }
				// (2) check if it's the same layer number as the jth property
				if (layerNo == massager->properties[j].layerNumber)
				  {					
					// Yes. Make suitable entries
					massager->properties[j].inputColNumber = i;
					massager->properties[j].inputPropNumber = i - offset;
					massager->properties[j].found = 1;
					if (massager->inputgrid.filt.apply)
					  {
						if (!strcmp(massager->inputgrid.filt.propName, massager->properties[j].name))
						  {
							massager->properties[j].applyFilter=1;
							// hack for single column filter only. For multiple filters, we will need to apply
							// a more general algorithm
							filterColumnNumber = i;
						  }
					  }
					nMatches++;
					break;  // out of properties j loop
				  }
			  }
		  }

		// X word match ?
		if (strncmp(line_buffer, massager->inputgrid.xWord, strlen(massager->inputgrid.xWord))==0)
		  {
			// matches. Last two entries in massager->properties are x and y
			massager->properties[massager->nProperties].inputColNumber = i;
			massager->properties[massager->nProperties].name = massager->inputgrid.xWord;
			nMatches++;
		  }
		// Y word match ?
		if (strncmp(line_buffer, massager->inputgrid.yWord, strlen(massager->inputgrid.yWord))==0)
		  {
			// matches. Last two entries in massager->properties are x and y
			massager->properties[massager->nProperties + 1].inputColNumber = i;
			massager->properties[massager->nProperties + 1].name = massager->inputgrid.yWord;
			nMatches++;
		  }
	  }
	// Set up check for a property name and layer number match.
	for (j=0; j<(massager->nProperties); j++) 
	  {
		if (massager->properties[j].found==0)
		  {
			warning("Analyser file %s has no info on property %s in layer %d requested in file %s\n",
					filename,
					massager->properties[j].name,
					massager->properties[j].layerNumber,
					filename);
			return EXIT_FAILURE;
		  }
	  }
	if ((nMatches != massager->nProperties + 2))
	  {
		warning("Can't match all properties and X Y words in file %s\n",filename);
		return EXIT_FAILURE;
	  }

	// now count number of useful lines => number of points in grid

	massager->inputgrid.nPoints=0;
	if (!is_binary)
	  {
		int line_ok;
		fgetpos(in,&pos);
		for(;;)
		  {
			line_ok=1;
			skipLine=0;
			// check we can read number_of_columns on this line
			for (i=0;i<number_of_columns; i++)
			  {
				if (fscanf(in," %lf ",&f)!=1) line_ok=0;
				if (massager->inputgrid.filt.apply && (i==filterColumnNumber) && line_ok)
				  {
					if ((f<massager->inputgrid.filt.min) || 
						(f>massager->inputgrid.filt.max)) skipLine=1;
				  }
			  }
			if (line_ok) 
			  {
				if (!skipLine) massager->inputgrid.nPoints++;
			  }
			else break;
		  }
		printf("ASCII analyse-output file: expect %d lines\n",massager->inputgrid.nPoints);
		// rewind
		if (fsetpos(in,&pos)) warning("Can't reposition file pointer");
	  }
	else
	  {
		// compute num lines by byte count
		// skip dangling "\n"
		if (fgetc(in) != (int)'\n') warning("Can't clear end of header in %s",filename);
		fgetpos(in,&pos);
		long nbytes, n_bytes_per_line;
		// nbytes = n_cov_entries*compressed_bytes_per_float + n_prop_entries*4 + number_of_coord_columns*4 
		n_bytes_per_line = n_prop_entries*compressed_bytes_per_float // covar diag
		  + (n_cov_entries-n_prop_entries)*2                         // covar off diag (2 bytes shorts)
		  + n_prop_entries*4 + number_of_coord_columns*4;            // props + coords as floats
		// goto file end
		fseek(in,0,SEEK_END);
		nbytes = ftell(in);
		// rewind
		if (fsetpos(in,&pos)) warning("Can't reposition file pointer");
		nbytes -= ftell(in);
		
		massager->inputgrid.nPoints = nbytes / n_bytes_per_line;
		if ((nbytes % n_bytes_per_line) != 0)
		  {
			printf("Binary file looks corrupt. Guessing %d lines\n",massager->inputgrid.nPoints);
		  }
		printf("Binary analyse-output file: expecting %d trace-locations\n",massager->inputgrid.nPoints);
	  }

	// allocate memory. A property per layer and name, plus two for x,y
	for (j=0; j<(massager->nProperties + 2); j++)
	  {
		massager->properties[j].values = (double*)malloc((size_t) ((massager->inputgrid.nPoints) * sizeof(double)));
	  }

	// observations covariance upper-tri memory
	massager->property_covariance_tri = (float**)malloc((size_t)((massager->inputgrid.nPoints) * sizeof(float*)));
	massager->property_covariance_tri_size = massager->inputgrid.nPoints;
	for (i=0; i< massager->inputgrid.nPoints; i++)
	  {
		massager->property_covariance_tri[i] = getUpperTriMatrix(massager->nProperties,massager->diagonalCovariancesOnly);
	  }
	massager->property_NAN = cvector(0,massager->inputgrid.nPoints-1);
	memset(&massager->property_NAN[0],0,massager->inputgrid.nPoints);

	// point input grid at xy members
	massager->inputgrid.xValues = massager->properties[massager->nProperties].values;
	massager->inputgrid.yValues = massager->properties[massager->nProperties + 1].values;

	// ptrs to data for parsing convenience
	data_ptr = (double**)malloc((size_t) (number_of_columns * sizeof(double *)));
	// properties ptrs
	for (i=0; i<number_of_columns; i++)
	  {
		data_ptr[i]=NULL;
		for (j=0;j<(massager->nProperties+2); j++)
		  {

			if (i==(massager->properties[j].inputColNumber)) 
			  {
				data_ptr[i] = &(massager->properties[j].values[0]);
				break;
		 	  }			
		  }
	  }
	// covariance ptrs
	cov_row_offset = ivector(0,number_of_columns-1);
	cov_col_offset = ivector(0,number_of_columns-1);
	for (i=0; i<number_of_columns; i++)
	  {
		cov_row_offset[i] = cov_col_offset[i] = -1;  // flag to signal not-needed until shown otherwise
		for (j=0;j<(massager->nProperties); j++)
		  {
			for (k=0;k<(massager->nProperties); k++)
			  {
				// if the column number matches the covariance column for the property pair (j,k)
				// in *working* order (not necc. the Summary-file order), save lookup indices
				if (i == cov_column[massager->properties[j].inputPropNumber][massager->properties[k].inputPropNumber]) 
				  {
					cov_row_offset[i]=j;				
					cov_col_offset[i]=k;
				  }			
			  }
		  }
	  }
	// now read
	nNonCovarianceColumns = number_of_columns-n_cov_entries;

	int inputGridIndx=0;
	if (!is_binary)
	  {
		for (line_no=0; inputGridIndx<massager->inputgrid.nPoints; line_no++)
		  {
			// filtering flag
			skipLine=0;
			// properties
			for (i=0; i<nNonCovarianceColumns; i++)
			  {
				if ((nMatches = fscanf(in," %lf ",&f)) == 1)
				  {
					// apply filter
					if (massager->inputgrid.filt.apply && (i==filterColumnNumber))
					  {
						if ((f>=massager->inputgrid.filt.min) && 
							(f<=massager->inputgrid.filt.max))
						  {
							data_ptr[i][inputGridIndx] = f;
							if (isnan(f)) 
							  {
								warning("NANs in input massage-analyse means!!! (Smoother will zero-weight them.)");							
							  } 
						  }
						else 
						  {
							skipLine=1;
							// a bad trace. Skip this
							// printf("Skipped IG point using filter\n");
						  }
					  }
					else
					  {
						if (data_ptr[i]!=NULL) 
						  {
							data_ptr[i][inputGridIndx] = f;
							if (isnan(f)) 
							  {
								warning("NANs in input massage-analyse means!!! (Smoother will zero-weight them.)");
								massager->property_NAN[inputGridIndx]=1;							
							  }
						  }
					  }
				  }
				else 
				  {
					warning("Can't read number in file %s, file-line-number %d\n",filename,line_no+2+number_of_columns);
					return EXIT_FAILURE;
				  }
			  }
			// covariances. Last columns are always covariances
			for (;i<number_of_columns;i++)
			  {
				if (fscanf(in," %lf ",&f)==1 )
				  {
					// insert into table only if needed. Ensure Symmetry too
					if (!massager->diagonalCovariancesOnly || (cov_row_offset[i] == cov_col_offset[i]))
					  {
						if (!skipLine && (cov_row_offset[i]>=0) && (cov_col_offset[i]>=0))
						  {
							if (isnan((double)f)) 
							  {
								warning("NANs in input massage-analyse covariance!!! (Smoother will zero-weight them.)");
								massager->property_NAN[inputGridIndx]+=2;
							  }
							if (cov_row_offset[i] <= cov_col_offset[i])
							  set(massager->property_covariance_tri[inputGridIndx], f, cov_row_offset[i], cov_col_offset[i]);
							else
							  set(massager->property_covariance_tri[inputGridIndx], f, cov_col_offset[i], cov_row_offset[i]);
						  }
					  } 
				  }
				else
				  {
					warning("Can't read number in file %s, file-line-number %d\n",filename,line_no+2+number_of_columns);
					return EXIT_FAILURE;
				  }
			  }
			fscanf(in," \n");
			// print_matrix0(massager->property_covariance[inputGridIndx],massager->nProperties, massager->nProperties,"cov.txt");
			if (!skipLine) inputGridIndx++;
		  }
	  }
	else
	  {
		// Binary input

		float f_float;
		short f_short;

		for (line_no=0; inputGridIndx<massager->inputgrid.nPoints; line_no++)
		  {
			// filtering flag
			skipLine=0;
			// properties. Expect floats. Read all 4 bytes
			for (i=0; i<nNonCovarianceColumns; i++)
			  {
				f_float=0;
				nMatches = fread(&f_float,4,1,in); // read 4 bytes
				if (nMatches == 1)
				  {
					// swap bytes if necc. input is java big-endian
					if (machine_is_little_endian) swap_float_4(&f_float);
					// apply filter
					if (massager->inputgrid.filt.apply && (i==filterColumnNumber))
 					  {
						if ((f_float>=massager->inputgrid.filt.min) && 
							(f_float<=massager->inputgrid.filt.max))
						  {
							data_ptr[i][inputGridIndx] = f_float;
							if (isnan(f_float)) 
							  {
								warning("NANs in input massage-analyse means!!! (Smoother will zero-weight them.)");
								massager->property_NAN[inputGridIndx]=1;
							  }
						  }
						else 
						  {
							skipLine=1;
							// a bad trace. Skip this
							// printf("Skipped IG point using filter\n");
						  }
					  }
					else
					  {
						if (data_ptr[i]!=NULL) 
						  {
							data_ptr[i][inputGridIndx] = f_float;
							if (isnan(f_float)) 
							  {
								warning("NANs in input massage-analyse means!!! (Smoother will zero-weight them.)");
								massager->property_NAN[inputGridIndx]=1;

								{
								  int jjj;
								  for (jjj=0;jjj<massager->nProperties; jjj++)
									{
									  if (i==(massager->properties[jjj].inputColNumber))
										{
										  warning("NANs for layer %d",massager->properties[jjj].layerNumber);
										  break;
										}
									}  
								}

							  }
						  }
					  }
				  }
				else 
				  {
					if ((i==0) && feof(in))
					  {
						goto CLEAN_EXIT; // legal exits occur early with filters...
					  }
					else
					  {
						warning("Can't read binary float number in file %s, file-line-number %d\n",filename,line_no+2+number_of_columns);
						return EXIT_FAILURE;
					  }
				  }
			  }
			// covariances. Read only compressed_bytes_per_float bytes
			// Last columns are always covariances
			// diagonal is compressed_bytes_per_float big endian float
			float* cov_tri = massager->property_covariance_tri[inputGridIndx];
			for (;i<number_of_columns;i++)
			  {
				// read leading e.g. 2 bytes. other 2 bytes will be zero
				if (cov_row_offset[i] == cov_col_offset[i])
				  {
					// diagonals are compressed floats
					f_float=0;
					nMatches = fread(&f_float,compressed_bytes_per_float,1,in);
					// swap if required
					if (machine_is_little_endian) swap_float_4(&f_float);
					/*
					  if (f_float==0.0 && cov_row_offset[i]>=0.0)
					  {
					  printf("ZEROCOV: %f %f\n",massager->inputgrid.xValues[inputGridIndx],
					  massager->inputgrid.yValues[inputGridIndx]);
					  }
					*/
				  }
				else 
				  {
					// upper triangles are 2 byte signed shorts
					nMatches = fread(&f_short,2,1,in);
					// swap if required
					if (machine_is_little_endian) swap_short_2(&f_short);
					f_float = (float) f_short;
				  }  
				if (nMatches == 1)
				  {
					// printf("COVIN: %f \n",f_float);
					// insert into table only if needed. Ensure Symmetry too
					if (!skipLine && (cov_row_offset[i]>=0) && (cov_col_offset[i]>=0))
					  {
						if (!massager->diagonalCovariancesOnly || (cov_row_offset[i] == cov_col_offset[i]))
						  {
							if (isnan((double)f_float))
							  {
								warning("NANs in input massage-analyse covariance!!! (Smoother will zero-weight them.)"); 
								massager->property_NAN[inputGridIndx]+=2;
							  }
							if (cov_row_offset[i] <= cov_col_offset[i])
							  {
								set(cov_tri, f_float, cov_row_offset[i],cov_col_offset[i]);
							  }
							else
							  {
								set(cov_tri, f_float, cov_col_offset[i], cov_row_offset[i]);
							  }
						  }
					  }
				  }
				else
				  {
					warning("Can't read compressed float in file %s, file-line-number %d\n",filename,line_no+2+number_of_columns);
					return EXIT_FAILURE;
				  }
			  }
			// restore off-diagonal covariances
			if (!massager->diagonalCovariancesOnly)
			  {
				for (i=0; i<(massager->nProperties); i++)
				  {
					for (j=i+1; j<(massager->nProperties); j++)
					  {
						f_float = get(cov_tri,i,j) * sqrt(MAX(0.0,(get(cov_tri,i,i) * get(cov_tri,j,j))))/short_scale;
						set(cov_tri,f_float,i,j);
					  }
				  }
			  } 
			if (!skipLine) inputGridIndx++;
			// print_dmatrix0(cov,massager->nProperties, massager->nProperties,"cov.txt");
		  }
	  }

 CLEAN_EXIT:

	fclose(in);

	free(data_ptr);
	free_imatrix(cov_column,0,n_prop_entries-1,0,n_prop_entries-1); cov_column=NULL;
	free_ivector(cov_row_offset,0,number_of_columns-1); cov_row_offset=NULL;
	free_ivector(cov_col_offset,0,number_of_columns-1); cov_col_offset=NULL;

    return EXIT_SUCCESS;
}
  
// Add in systematic depth uncertainty (if this wasn't part of the massage-analyse
// covariance files (i.e. explicitly modelled in delivery)
int addSystematicDepthUncertainty(struct MASSAGER *massager)
{
  if (massager->systematicDepthUncertainty > 0)
	{
	  int i,j,k;
	  int *depthCols,nDepths;
	  double varDepth=0;
	  float** cov;
	  float* cov_tri;

	  // i is the input grid index. The covariance is blockdiagonal, each block
	  // indexed by i...

	  if (massager->diagonalCovariancesOnly)
		{
		  printf("Nonsensical to add reference depth uncertainty if using diagonal-only covariances! Exiting.\n");
		  exit(0);
		}

	  if (massager->verbosity>1)
		{
		  printf("Adding reference depth uncertainty of std. dev. = %.3f\n",
				 massager->systematicDepthUncertainty);
		}

	  depthCols = ivector(0,massager->nProperties);
	  nDepths=0;
	  for (j=0;j<(massager->nProperties);j++)
		{
		  if (!strcmp(massager->properties[j].name,"d")) depthCols[nDepths++]=j; // is a depth
		}
	  // variance component to add
	  varDepth = DSQR(massager->systematicDepthUncertainty);

	  for (i=0; i<massager->inputgrid.nPoints; i++)
		{
		  cov_tri = massager->property_covariance_tri[i];
#if 0
		  // Fragment to test pos. definiteness of input matrices
		  // print_matrix(cov,massager->nProperties,'d',NULL);
		  double max=0;
		  for (j=0;j<(massager->nProperties);j++) max=MAX(max,cov[j][j]);
		  for (j=0;j<(massager->nProperties);j++) cov[j][j] += 1e-6*max;
		  int info = clapack_wrapper_spotrf(CblasLower,massager->nProperties,&cov[0][0],massager->nProperties);

		  if (info>0) 
			{
			  // fail detected. 
			  warning ("Trace %d: Input Property Covariance not pos definite!\n",i);
			  exit(0);
			}
#endif
		  for (j=0;j<nDepths;j++)
			{
			  for (k=0;k<nDepths;k++)
				{
				  if (k>=j) 
					{
					  if (massager->diagonalCovariancesOnly && (depthCols[j] != depthCols[k])) continue;
					  set( cov_tri, varDepth + get(cov_tri,depthCols[j],depthCols[k]), depthCols[j], depthCols[k]);
					}
				}
			}
#if 0
		  // Fragment to test pos. definiteness of input matrices after depth ref. hack
		  // print_matrix(cov,massager->nProperties,'f',NULL);
		  double max=0;
		  for (j=0;j<(massager->nProperties);j++) max=MAX(max,cov[j][j]);
		  for (j=0;j<(massager->nProperties);j++) cov[j][j] += 1e-4*max;
		  int info = clapack_wrapper_spotrf(CblasLower,massager->nProperties,&cov[0][0],massager->nProperties);

		  if (info>0) 
			{
			  // fail detected. 
			  warning ("Trace %d: Input Property Covariance not pos definite!\n",i);
			  exit(0);
			}
#endif
		}
	  free_ivector(depthCols,0,massager->nProperties);
	}
  return EXIT_SUCCESS;
}


// free massage-analysed covar data put out by InversionAnalyser --massage-analyse
int freeMassageAnalysisCovariances(struct MASSAGER *massager)
{
  int i;

  if (massager->property_covariance_tri != NULL)
	{
	  // care with size here as defragging may have clawed back memory, and
	  // massager->property_covariance_tri_size < massager->inputgrid.nPoints
	  for (i=0; i< massager->property_covariance_tri_size; i++)
		{
		  if (massager->property_covariance_tri[i])
			{
			  freeUpperTriMatrix(massager->property_covariance_tri[i],massager->nProperties);
			  massager->property_covariance_tri[i] = NULL;
			} 
		}
	  free(massager->property_covariance_tri); 
	  massager->property_covariance_tri = NULL;
	}

  if (massager->property_NAN)
	{
	  free_cvector(massager->property_NAN,0,massager->inputgrid.nPoints-1);
	  massager->property_NAN = NULL;
	}

  return EXIT_SUCCESS;
}


// read massage-analysed data put out by InversionAnalyser --massage-analyse
int freeMassageAnalysisData(struct MASSAGER *massager)
{
  int i,j;

  // de allocate memory. A property per layer and name, plus two for x,y
  for (j=0; j<(massager->nProperties + 2); j++)
	{
	  free(massager->properties[j].values);
	  massager->properties[j].values = NULL;
	}

  // observations covariance memory
  freeMassageAnalysisCovariances(massager);  
 
  return EXIT_SUCCESS;
}



void allocateObservations(struct MASSAGER *massager, int layer_indx)
{
  struct WELL_OBSERVATIONS *wo = &massager->well_observations[layer_indx];

  if (wo->n > 0)
	{
	  wo->xValues = g_array_sized_new(FALSE, FALSE, sizeof (double), wo->n);
	  wo->yValues = g_array_sized_new(FALSE, FALSE, sizeof (double), wo->n);
	  wo->values = g_array_sized_new(FALSE, FALSE, sizeof (double), wo->n);
	  wo->propertyIndex = g_array_sized_new(FALSE, FALSE, sizeof (int), wo->n);
	  wo->property_name = g_array_sized_new(FALSE, FALSE, sizeof (char*), wo->n);
	}
}

void freeObservations(struct MASSAGER *massager, int layer_indx)
{
  int k;
  struct WELL_OBSERVATIONS *wo = &massager->well_observations[layer_indx];

  if (wo->n > 0)
	{
	  g_array_free(wo->xValues, TRUE);
	  wo->xValues=NULL;
	  g_array_free(wo->yValues, TRUE);
	  wo->yValues=NULL;
	  g_array_free(wo->values, TRUE);
	  wo->values=NULL;
	  g_array_free(wo->propertyIndex, TRUE); 
	  wo->propertyIndex=NULL;
	  g_array_free(wo->property_name, TRUE); 
	  wo->property_name=NULL;

	  for (k=0; k< massager->nPropertyNames; k++)
		{
		  if (wo->n_obs[k] > 0)
			{
			  free(wo->obs[k]); 
			  wo->obs[k] = NULL;
			}
		}
	}
}


// set up possible SU outputs
int setUpSUOutputs(struct MASSAGER *ma)
{
  if (ma->writeSUOut)
	{
	  
	  int m,j;
	  struct SU_OUT_INFO* suinfo;

	  // Allocate
	  ma->su_out_info = (struct SU_OUT_INFO*)malloc((size_t)(sizeof(struct SU_OUT_INFO))); // FIXME clean up!
	  suinfo = ma->su_out_info;
	  suinfo->nLayersPossibleForSUOutput = ma->nLayers;
	  suinfo->inTime = 0;
	  suinfo->getLastSurfacefromThickness = 0;
	  suinfo->getLastSurfacefromTbase = 0;

	  if (ma->t_index >= 0)
		{
		  // assume this means we want output in time
		  suinfo->inTime = 1;
		  if ((ma->verbosity > 1 )) printf("SU output in time\n");
		  if (ma->t_base_index >= 0)
			{
			  suinfo->getLastSurfacefromTbase = 1;
			}
		  else if ((ma->vp_eff_index >= 0) && (ma->thickness_index >= 0))
			{
			  // can compute bottom layers from vp_eff & thickness
			  suinfo->getLastSurfacefromThickness = 1;
			}
		  else
			{
			  suinfo->nLayersPossibleForSUOutput--;  // cannot show last layer
			  if (ma->verbosity > 1) 
				{
				  printf("SU output in time: cannot show last layer unless thickness and vp_eff massaged!\n");
				}
			}
		}
	  else if (ma->d_index >= 0)
		{
		  // assume this means we want output in depth
		  if (ma->verbosity > 1) printf("SU output in depth\n");
		  if (ma->thickness_index >= 0)
			{
			  // can compute bottom layers by adding thickness
			  suinfo->getLastSurfacefromThickness = 1;
			}
		  else
			{
			  suinfo->nLayersPossibleForSUOutput--;  // cannot show last layer
			  if (ma->verbosity > 1) printf("SU output in depth: cannot show last layer unless thickness massaged");
			}
		}
	  else
		{
		  printf("SU output not possible unless either depth (d) or time (t) is massaged\n");
		  return EXIT_FAILURE;
		}
	  
	  // count SU-depictable body properties. Must be same criterion as loop in SU output code
	  for (m=0,j=0; m < ma->nPropertyNames; m++) 
		{
		  if (indexIsADepictableSUBodyProperty(ma,m)) j++;
		}
	  suinfo->nBodyPropertiesPossibleForSUOutput = j;
	}
  return EXIT_SUCCESS;
}


// set up well observations property indexes
int setUpWellObservations(struct MASSAGER *massager)
{
  double x,y,val,denom;
  int i,j,jj,k,m;
  int convertNG_obs=0,newObsCount=0;
  struct WELL_OBSERVATIONS *wo;

  // First, convert all NG constraints to equivalents on thickness or net-sand if possible
  if ((massager->thickness_index>=0) && (massager->net_sand_index>=0) && (massager->NG_index>=0))
	{
	  printf("Converting NG well-observations to equivalent net-sand/thickness ones\n");
	  int matched=0;

	  for (i=0; i<massager->nLayers; i++)
		{
		  wo = &massager->well_observations[i];
		  newObsCount=0;
		  for (j=0,matched=0; j<wo->n; j++)
			{
			  // is this a NG observation?
			  if (!strcmp(g_array_index(wo->property_name,gchar*,j), massager->propertyNames[massager->NG_index]))
				{
				  // Yes. here are the x,y's
				  x = g_array_index (wo->xValues, gdouble, j);
				  y = g_array_index (wo->yValues, gdouble, j);
				  // Is there a thickness or net-sand at the same place?
				  for (jj=0; jj<wo->n; jj++)
					{
					  // is this a thickness observation?
					  if (!strcmp(g_array_index(wo->property_name,gchar*,jj), 
								  massager->propertyNames[massager->thickness_index]))
						{
						  if (APPROX(g_array_index (wo->xValues, gdouble, jj),x,WELL_XY_MATCH_EPS) && 
							  APPROX(g_array_index (wo->yValues, gdouble, jj),y,WELL_XY_MATCH_EPS))
							{
							  // and it's at the same place. 
							  // Create a new equivalent net-sand observation
							  // wo->property_name[j] = massager->propertyNames[massager->net_sand_index];
							  // wo->values[j] *= wo->values[jj];

							  val = g_array_index (wo->values, gdouble, j) * g_array_index (wo->values, gdouble, jj);
							  g_array_append_val(wo->values,val); 
							  g_array_append_val(wo->property_name, massager->propertyNames[massager->net_sand_index]); 
							  g_array_append_val(wo->xValues,x); 
							  g_array_append_val(wo->yValues,y); 
							  newObsCount++;
							  matched=1; 
							  break;
							}
						} 
					  else if (!strcmp(g_array_index(wo->property_name,gchar*,jj), massager->propertyNames[massager->net_sand_index]))
						{
						  // its a net-sand observation
						  if (APPROX(g_array_index (wo->xValues, gdouble, jj),x,WELL_XY_MATCH_EPS) && 
							  APPROX(g_array_index (wo->yValues, gdouble, jj),y,WELL_XY_MATCH_EPS))
							{
							  // and its at the same place. 
							  // Convert to thickness observation
							  // wo->property_name[j] = massager->propertyNames[massager->thickness_index];
							  // wo->values[j] = (wo->values[j]!=0)? wo->values[jj]/wo->values[j] : 0.0;
							  denom = g_array_index (wo->values, gdouble, j);
							  val = (denom!=0)? g_array_index (wo->values, gdouble, jj)/denom : 0.0;

							  g_array_append_val(wo->values,val); 
							  g_array_append_val(wo->property_name, massager->propertyNames[massager->thickness_index]); 
							  g_array_append_val(wo->xValues,x); 
							  g_array_append_val(wo->yValues,y); 
							  newObsCount++;		 
							  matched=1; 
							  break;
							}
						}
					}
				  if (!matched)
					{
					  printf("\nWarning: Unconverted NG well -observation in layer %d!\n",massager->layerNumbers[i]);
					}
				}
			}
		  // allow for extra observations
		  wo->n += newObsCount;
		}
	}


  for (i=0; i<massager->nLayers; i++)
	{
	  wo = &massager->well_observations[i];
	  wo->n_obs = ivector(0,massager->nPropertyNames-1);                  // FIXME free at end!!!
	  for (k=0; k<massager->nPropertyNames; k++) wo->n_obs[k]=0;
	 
	  for (j=0; j<wo->n; j++)
		{
		  /* flag properties as -1 initially */
		  // wo->propertyIndex[j] = -1;
		  k=-1;
		  g_array_append_val(wo->propertyIndex,k);

		  for (k=0; k< massager->nPropertyNames; k++)
			{
			  if (!strcmp(g_array_index (wo->property_name, gchar*, j), massager->propertyNames[k]))
				{
				  /* a match */
				  // wo->propertyIndex[j] = k;
				  g_array_insert_val(wo->propertyIndex,j,k);
				  (wo->n_obs[k])++;    // count # observations of property k in layer i
				}
			}
		  /* Now all observations of non-modelled properties have index -1 */
		}
	  // allocate convenience structures
	  wo->obs = (struct OBS**)malloc((size_t) (massager->nPropertyNames * sizeof(struct OBS*)));
	  // fill them
	  for (k=0; k< massager->nPropertyNames; k++)
		{
		  if (wo->n_obs[k] > 0)
			{
			  wo->obs[k] = (struct OBS*)malloc((size_t) (wo->n_obs[k] * sizeof(struct OBS)));
			  struct OBS* obs = wo->obs[k];
			  for (j=0,m=0; j<wo->n; j++)
				{
				  // if (wo->propertyIndex[j] == k)
				  if (g_array_index (wo->propertyIndex, gint , j) == k)
					{
					  // printf("Layer %d, added obs for %s\n",i+1,massager->propertyNames[k]);
					  obs[m].x =  g_array_index (wo->xValues, gdouble, j); // wo->xValues[j];
					  obs[m].y =  g_array_index (wo->yValues, gdouble, j); // wo->yValues[j];
					  obs[m].val =  g_array_index (wo->values, gdouble, j); // wo->values[j];
					  m++;
					}
				}
			}
		}
	}
  return EXIT_SUCCESS;
}


/* 
   Dump out file containing
   properties=A,B,....
   for use by viewing engines.
*/
int dumpBHPviewerPropertyKeyFile(struct MASSAGER *massager,char *pfilename, int includeUncertainties)
{
  int k,m;

  FILE *fp = NULL;

  if ((fp = fopen(pfilename,"w")) == NULL)
	{
	  printf("Cannot dump output file %s for SU properties\n",pfilename);
	  return EXIT_FAILURE;
	}
  fprintf(fp,"properties=");
  // means
  for (m=0,k=0; m< massager->nPropertyNames; m++) 
	{
	  if (indexIsADepictableSUBodyProperty(massager,m))
		{
		  // A depictable property
		  // Leave out t & t-base here (time-style) and d (depth style)
		  if (k>0) fprintf(fp,",");
		  fprintf(fp,"%s",massager->propertyNames[m]);
		  k++;
		}
	}
  if (includeUncertainties)
	{
	  // stddeviations
	  for (m=0; m< massager->nPropertyNames; m++) 
		{
		  // All stddevs are depictable
		  if (k>0) fprintf(fp,",");
		  fprintf(fp,"std_dev_%s",massager->propertyNames[m]);
		  k++;
		}
	}
  fprintf(fp,"\n");
  fclose(fp);
  
  return EXIT_SUCCESS;
}

/*
  check of whether index ought to be su - depictable.
  Returns  1 or 0.
  Exists to exclude t_base
*/
char indexIsADepictableSUBodyProperty(struct MASSAGER *ma, int m)
{
  //   if ((m != ma->thickness_index) && (m != ma->d_index) && (m != ma->t_index) && (m != ma->t_base_index))
  if ((m != ma->t_base_index))
	{
	  return '\1';
	}
  else 
	{
	  return '\0'; 
	}
}

int setUpUTMmapping(struct MASSAGER *massager)
{
  int i,j,k;
  double TOL = 1e-6;

  double DR[2][2];
  double DS[2][2];
  double det;
  struct PT2D pt;

  /* form DR = (r0-r1 | r0-r2) */
  DR[0][0] = massager->utmInfo.x[0] - massager->utmInfo.x[1];
  DR[0][1] = massager->utmInfo.x[0] - massager->utmInfo.x[2];
  DR[1][0] = massager->utmInfo.y[0] - massager->utmInfo.y[1];
  DR[1][1] = massager->utmInfo.y[0] - massager->utmInfo.y[2];
  /* form DS = (s0-s1 | s0-s2) */
  DS[0][0] = massager->utmInfo.inLine[0] - massager->utmInfo.inLine[1];
  DS[0][1] = massager->utmInfo.inLine[0] - massager->utmInfo.inLine[2];
  DS[1][0] = massager->utmInfo.xLine[0] - massager->utmInfo.xLine[1];
  DS[1][1] = massager->utmInfo.xLine[0] - massager->utmInfo.xLine[2];

  det = DS[0][0] * DS[1][1] - DS[1][0] * DS[0][1];
  if (fabs(det) < TOL)
	{
	  printf ("UTM coordinates in XML file are not unique, or collinear");
	  return EXIT_FAILURE;
	}
  /* form mapMatrix = DR.DS^{-1} */
  massager->utmInfo.mapMatrix[0][0] = (-DR[0][1]*DS[1][0] + DR[0][0]*DS[1][1])/det;
  massager->utmInfo.mapMatrix[0][1] = (DR[0][1]*DS[0][0] - DR[0][0]*DS[0][1])/det;
  massager->utmInfo.mapMatrix[1][0] = (-DR[1][1]*DS[1][0] + DR[1][0]*DS[1][1])/det;
  massager->utmInfo.mapMatrix[1][1] =  (DR[1][1]*DS[0][0] - DR[1][0]*DS[0][1])/det;

  /* set up offset */
  massager->utmInfo.mapOffset[0] = 0.0;
  massager->utmInfo.mapOffset[1] = 0.0;
  pt = getUTMSfromInlineXline(massager->utmInfo,massager->utmInfo.inLine[0],massager->utmInfo.xLine[0]);
  massager->utmInfo.mapOffset[0] = massager->utmInfo.x[0] - pt.x;
  massager->utmInfo.mapOffset[1] = massager->utmInfo.y[0] - pt.y;

  /* Do checks... */ 
  double scale=0.0;
  for (i=0;i<3;i++)
	{
	  pt = getUTMSfromInlineXline(massager->utmInfo,massager->utmInfo.inLine[i],massager->utmInfo.xLine[i]);
	  // printf ("UTM transformation check %g %g\n",pt.x - massager->utmInfo.x[i],pt.y - massager->utmInfo.y[i]);
	  scale = MAX(scale,MAX(fabs(pt.x),fabs(pt.y)));
	}
  if (scale==0.0) scale=1.0;
  scale *= TOL;
  for (i=0;i<3;i++)
	{
	  pt = getUTMSfromInlineXline(massager->utmInfo,massager->utmInfo.inLine[i],massager->utmInfo.xLine[i]);
	  if (((fabs(pt.x - massager->utmInfo.x[i])) > scale) || ((fabs(pt.y - massager->utmInfo.y[i])) > scale))
		{
		  printf ("UTM transformation is buggy, point %d\n",i);
		}  
	}


  return EXIT_SUCCESS;
}
 
struct PT2D getUTMSfromInlineXline(struct UTM_CONVERSION_INFO info,double inLine,double xLine)
{
  static struct PT2D pt;

  pt.x = info.mapMatrix[0][0] * inLine + info.mapMatrix[0][1] * xLine + info.mapOffset[0];
  pt.y = info.mapMatrix[1][0] * inLine + info.mapMatrix[1][1] * xLine + info.mapOffset[1];

  return pt;
}
