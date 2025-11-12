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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <glib.h>

#include "NR/nrutil.h"

/* Eclipse handling */ 
#include "ECLIO/eclio.h"  
#include "ECLIO/access.h"
#include "ECLIO/eclutil.h"
#include "ECLIO/blockops.h"
#include "ECLIO/set.h"

/* Neighbour searches */
#include "kdtree/nn.h"
#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"

/* Delaunay, bounding rectangles */
#include "cgeom/vdefs.h"
#include "cgeom/2Dhull.h"

#include "massager_types.h"
#include "clapack_cblas_wrappers.h"
#include "massager.h"
#include "eclipseHandler.h"
#include "util.h"
#include "errorHandler.h"

#define PINCH_DEFAULT 0.001
#define GAP_DEFAULT "GAP"
#define PINCHPV_DEFAULT 1.0E20
#define EPS 0.01
#define TOPBOT_DEFAULT "TOPBOT"
#define TMETH_DEFAULT "TOP"
#define DEFAULT_L_UNITS "METRES"

#define ACCEPTABLE_OFF_LATTICE_ECLIPSE_GRID_MAKING_FRACTION 0.05
#define ACCEPTABLE_EPSILON_COSINE_ECLIPSE_GRID_MAKING 0.05
#define MAX_ECLIPSE_GRID_MAKING_ASPECT_RATIO 20

/*
  To be called only if depth and thickness massaged. Gets min and max z of input data
*/
int getDepthRangeInputGrid(struct MASSAGER *massager, double *zMin, double *zMax)
{
  int i,j;
  int properties_indx_offset;
  double z;

  *zMin = DBL_MAX;
  *zMax = -DBL_MAX;

  if ((massager->d_index< 0) || (massager->thickness_index<0)) return EXIT_FAILURE;

  for (i=0; i< massager->nLayers ; i++)
	{ 
	  properties_indx_offset = i*massager->nPropertyNames;
	  for (j=0; j < massager->inputgrid.nPoints; j++) 
		{	
		  z = massager->properties[massager->d_index + properties_indx_offset].values[j];
		  *zMin = MIN(z, *zMin);
		  *zMax = MAX(z, *zMax);
		  if (i==(massager->nLayers-1))
			{
			  z += massager->properties[massager->thickness_index + properties_indx_offset].values[j];
			  *zMin = MIN(z, *zMin);
			  *zMax = MAX(z, *zMax);
			}
		}
	}
  printf("Input grid depth range: %g %g\n",*zMin,*zMax);

  return EXIT_SUCCESS; 
}


int makeEclipseGrid(struct MASSAGER *massager)
{
  int i,j,k,n;
  double *x,*y;
  double rBlockSq,r;
  double pt[2];
  struct NEAREST_REGULAR_GRID_INFO RGinfo;
  struct VORONOI_INFO vinfo;
  char isLine=0;

  // Eclipse variables
  struct ECL ecl;
  int nKeys;
  int nx,ny,nz;
  int transformToMap;
  double **axes=NULL;
  double u,v,X,Y;
  double z_line_bottom,z_line_top,delta_z;
  double ***c = NULL;
  PROP *pp = NULL;
  int ***vi = NULL;
  int **bnd = NULL;
  double ****z = NULL;
  int *k1 = NULL, *k2 = NULL;
  char **comp = NULL; char **join = NULL;

  n = massager->inputgrid.nPoints;

  // Eclipse grid will be based on input grid from seismic/Delivery
  // Points will be transformed to UTMs if UTM info is available
  if (massager->utmInfo.performInputGridUTMmapping) 
	{
	  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  // Eclipse points need transformation been done on the input grid
	  x = dvector(0,n-1);
	  y = dvector(0,n-1);
	  struct PT2D pt_utm;
	  for (j=0; j < n; j++) 
		{ 
		  pt_utm = getUTMSfromInlineXline(massager->utmInfo, 
										  massager->inputgrid.xValues[j],
										  massager->inputgrid.yValues[j]);
		  x[j] = pt_utm.x; 
		  y[j] = pt_utm.y;
		}
	}
  else
	{
	  // use original x, y from input grid
	  x = massager->inputgrid.xValues;
	  y = massager->inputgrid.yValues;
	}


  if (!massager->make_ecl_grid_options.MBR)
	{
	  /* Preferentially use Delaunay triangulation as basis for eclipse grid */
	  vinfo = getdiagram(0,1,n,'d',x,y);
	  // Then nearest regularly spaced grid
	  RGinfo = nearest_regular_lattice(&vinfo,0,n,'d',x,y,
									   ACCEPTABLE_EPSILON_COSINE_ECLIPSE_GRID_MAKING,
									   ACCEPTABLE_OFF_LATTICE_ECLIPSE_GRID_MAKING_FRACTION);
	  free_voronoi_info(&vinfo);
	  /* 
		 RGinfo now contains:
 
		 RGinfo.Nu, RGinfo.Nv        Number of Blocks in lattice
		 RGinfo.u_hat, RGinfo.v_hat  Unit vectors in the principal directions
		 RGinfo.u_min/max,RGinfo.v_min/max, are maximum extents in these directions
		 RGinfo.du/dv                Block widths in these directions
	  */
	  if (RGinfo.info == -2) 
		{
		  printf("Unable to match seismic grid to a regular areal 2D lattice: info status=%d\n",RGinfo.info);
		  printf("Fitting to 2D line model\n");
		  isLine=1;
		}
	  else if (RGinfo.info == -3) 
		{
		  printf("Unable to match seismic grid to a regular areal 2D lattice, or a regular line: info status=%d\n",RGinfo.info);
		  printf("Forcing equispaced 2D line model\n");
		  isLine=1;
		}
	  else if (RGinfo.info == -4) 
		{
		  printf("Unable to match seismic grid to a regular areal 2D lattice: info status=%d\n",RGinfo.info);
		  printf("Forcing equispaced 3D model\n");
		  // Random looking points. desperation. Use Minimum bounding rectangle & rough guess of cell sizes
		  if (getEclipseGriddingFromBoundingRectangle(n,x,y,&RGinfo) == EXIT_FAILURE)
			{
			  return EXIT_FAILURE;
			}
		}
	  else if (RGinfo.info == -10)
		{
		  printf("Memory alloc errors in makeEclipseGrid(). Giving up\n");
		  return EXIT_FAILURE;
		} 
	}
  else
	{
	  // Use minimum bounding rectangle as basis for grid
	  if (getEclipseGriddingFromBoundingRectangle(n,x,y,&RGinfo) == EXIT_FAILURE)
		{
		  return EXIT_FAILURE;
		}
	}
  
	// Alter spacing if requested
  if ((massager->make_ecl_grid_options.Nu!=-1) || (massager->make_ecl_grid_options.Nv!=-1))
	{
	  if (massager->make_ecl_grid_options.Nu != -1) RGinfo.Nu = massager->make_ecl_grid_options.Nu;
	  if (massager->make_ecl_grid_options.Nv != -1) RGinfo.Nv = massager->make_ecl_grid_options.Nv;
	
	  RGinfo.du = (RGinfo.u_max - RGinfo.u_min)/RGinfo.Nu;
	  RGinfo.dv = (RGinfo.v_max - RGinfo.v_min)/RGinfo.Nv;
	}



  /* Construct eclipse file with basic geometry */
#define NKEYS 9
  nKeys = NKEYS;
  char* keyWords[NKEYS+1] 
	= {NULL,"PINCH", "MAPAXES", "GRIDUNIT", "SPECGRID", "COORDSYS", "COORD","ZCORN", "ACTNUM","SEGMENTS"};
  // for now override massager->outputgrid[0].ECLIPSE_segments_keyword 

  

  // Dimensions
  if (isLine)
	{
	  // Nv will be 0, dv=0
	  nx = MAX(RGinfo.Nu,1);
	  ny = MAX(RGinfo.Nv,1);
	  // center the geometry on the line. Make blocks square.
	  RGinfo.dv = RGinfo.du;
	  RGinfo.v_min -= RGinfo.dv/2.0;
	  RGinfo.v_max += RGinfo.dv/2.0;
	}
  else
	{
	  nx = RGinfo.Nu;
	  ny = RGinfo.Nv;
	}
  
  if (1)
	{
	  // Show plan(s) in ascii file
	  int pass;
	  double X0,Y0;
	  FILE *outPlan= NULL;
	  char* filename;
	  for (pass=0;pass<=1;pass++)
		{
		  asprintf(&filename,"EclipseGriddingPlan%s.txt",(pass==0)?"":"Relative");
		  outPlan = fopen(filename, "w");
	  
		  fprintf(outPlan,"#m=0,S=2\n#Massage-analysis %s trace grid from massage-analysis "
				  "file & any geometric (UTM) transformation\n",(pass==0)?"absolute":"relative");
		  X0=x[0]; Y0=y[0];
		  for (i=0;i<n;i++) fprintf(outPlan,"%lg %lg\n",x[i]-X0,y[i]-Y0);
		  // Eclipse coord lines, as per loop below
		  fprintf(outPlan,"\n#m=0,S=1\n# Eclipse coord line coordinates\n");
		  for (j=0;j<=ny;j++) 
			{
			  v = j * RGinfo.dv + RGinfo.v_min;
			  for (i=0;i<=nx;i++) 
				{
				  u = i * RGinfo.du + RGinfo.u_min;
				  X = u * RGinfo.uhat[0] + v * RGinfo.vhat[0];
				  Y = u * RGinfo.uhat[1] + v * RGinfo.vhat[1];
				  fprintf(outPlan,"%lg %lg\n",X-X0,Y-Y0);
				}
			}
		  fclose(outPlan);
		  printf("Wrote Eclipse %s gridding plan to %s\n",((pass==0)?"":"relative"),filename); 
		  free(filename);
		}
	}



  nz = massager->nLayers;
  transformToMap = 1; 	/* force transform grid->map coords */

  initModel(&ecl);
  // PINCH
  s_pinch(&ecl, PINCH_DEFAULT, GAP_DEFAULT, PINCHPV_DEFAULT, TOPBOT_DEFAULT, TMETH_DEFAULT);

  // MAPAXES
  axes = dmatrix(1,3,1,2);
  axes[1][1] = 0;  // pt on y-axis
  axes[1][2] = 1;
  axes[2][1] = 0;  // origin
  axes[2][2] = 0;
  axes[3][1] = 1;  // pt on x-axis
  axes[3][2] = 0;
  s_mapAxes(&ecl, axes);

  // GRIDUNIT
  s_gridUnit(&ecl, DEFAULT_L_UNITS,"MAP");

  // SPECGRID
  s_specGrid(&ecl, nx, ny, nz,1,'F');

  // COORDSYS  
  k1 = ivector(1,1); k2 = ivector(1, 1);
  comp = cmatrix(1, 1, 0, TOKEN_LENGTH); 
  join = cmatrix(1, 1, 0, TOKEN_LENGTH);
  k1[1]=1;
  k2[1]=nz; 
  memset(join[1],'\0',TOKEN_LENGTH); 
  strcpy(comp[1],"INCOMP");
  s_coordSys(&ecl,1, k1, k2, comp, join);  

  // COORD
  // corner lines
  pp = addProp(&ecl, "COORD", 'd', 3, 0, nx, 0, ny, 1, 6);
  if (!pp || !pp->p) 
	{
	  fprintf(stderr,"makeEclipseGrid(): Could not init COORD\n");
	  exit(0);
	}
  c = (double ***)pp->p;

  // get rough depth limits from massaged input stats
  getDepthRangeInputGrid(massager,&z_line_top,&z_line_bottom);

  for (j=0;j<=ny;j++) 
	{
	  v = j * RGinfo.dv + RGinfo.v_min;
	  for (i=0;i<=nx;i++) 
		{
		  // for (k=1;k<=6;k++) {
		  // c[i][j][k] = {xt,yt,zt,xb,yb,zb} of line
		  u = i * RGinfo.du + RGinfo.u_min;
		  X = u * RGinfo.uhat[0] + v * RGinfo.vhat[0];
		  Y = u * RGinfo.uhat[1] + v * RGinfo.vhat[1];

		  c[i][j][1] = X;          // X
		  c[i][j][4] = c[i][j][1];
		  // y
		  c[i][j][2] = Y;
		  c[i][j][5] = c[i][j][2];
		  // z
		  c[i][j][3] = z_line_top;
		  c[i][j][6] = z_line_bottom;
		}
	}
  ecl.coordp = pp;
  ecl.coord = c;
  pp->set = 1;



  // ZCORN
  pp = addProp(&ecl, "ZCORN", 'd', 4, 1, nx, 1, ny, 1, nz, 1, 8);
  z = (double ****)pp->p;
  if (!z) 
	{
	  fprintf(stderr,"makeEclipseGrid(): Could not access ZCORN array data\n"); 
	  exit(0);
	}

  // Token equispaced placeholder depths
  delta_z = (z_line_bottom - z_line_top)/nz;
  for (i=1;i<=nx;i++) 
	{
	  for (j=1;j<=ny;j++) 
		{
		  for (k=1;k<=nz;k++) 
			{
			  double* zp = z[i][j][k];		
			  // use zCell function to get heights
			  //zCell(i,k,&zp[1],&zp[2],&zp[5],&zp[6],&m);
			  zp[1]=zp[2]= (k==1) ? z_line_top : z[i][j][k-1][5];
			  zp[5]=zp[6]= (zp[1] + delta_z);
			  // flat in y-direction
			  zp[3] = zp[1];
			  zp[4] = zp[2];
			  zp[7] = zp[5];
			  zp[8] = zp[6];
			}
		}
	}

  ecl.zCornp = pp;
  ecl.zCorn = z;
  pp->set = 1;

  // ACTNUM
  bnd = buildBndGrid(&ecl);
  pp = addPropN(&ecl,"ACTNUM",'i',3, bnd);
  if (!pp) 
	{
	  fprintf(stderr,"makeEclipseGrid(): Could not add ACTNUM for ACTNUM\n");exit(0);
	}
  vi = (int ***) (pp->p);
  if (0)
	{
	  // Set everything active
	  for (i=1;i<=nx;i++)
		for (j=1;j<=ny;j++) 
		  for (k=1;k<=nz;k++) vi[i][j][k]=1;
	}
  else
	{
	  // Stuff to switch off distant cells based on distance to original traces
	  int *neigh;
	  optkdNode *OptkdRoot;
	  double** Points;
	  
	  Points = dmatrix(0,n-1,0,1);
	  for (j=0; j < n; j++) 
		{
		  Points[j][0] = x[j];
		  Points[j][1] = y[j];
		}
		
	  // build kdtree
	  OptkdRoot = BuildOptTree(Points,n,2);
	  if (!OptkdRoot) 
		{ 
		  fprintf(stderr,"makeEclipseGrid(): Could not add refine actives. BuildOptTree() problems\n");
		}
	  else
		{
		  rBlockSq = SQR((RGinfo.du/2.0)) + SQR((RGinfo.dv/2.0));
		  // loop over eclipse block-centres
		  for (j=0;j<ny;j++) 
			{
			  v = (j+0.5) * RGinfo.dv + RGinfo.v_min;
			  for (i=0;i<nx;i++) 
				{
				  u = (i +0.5) * RGinfo.du + RGinfo.u_min;
				  pt[0] = u * RGinfo.uhat[0] + v * RGinfo.vhat[0]; // X
				  pt[1] = u * RGinfo.uhat[1] + v * RGinfo.vhat[1]; // Y
				  // get nearest neighbour of this point from kdtree
				  neigh = kdOptNNQuery(Points,2,&pt[0],1,EUCLIDEAN,OptkdRoot,1);
				  r = SQR(Points[neigh[1]][0] - pt[0]) + SQR(Points[neigh[1]][1] - pt[1]);
				  if (r>rBlockSq)
					{
					  // deem as inactive
					  for (k=1;k<=nz;k++) vi[i+1][j+1][k]=0;		
					}
				  else
					{
					  // deem as active
					  for (k=1;k<=nz;k++) vi[i+1][j+1][k]=1;	
					}
				  free(neigh);  // neighbourhood of jth observation
				  neigh=NULL;
				}
			}
		  KillOptTree(OptkdRoot);
		}
	  free_dmatrix(Points,0,n-1,0,1);
	}
  pp->set = 1;

  // SEGMENTS (default to 1, as no fault block determination possible)
  bnd = buildBndGrid(&ecl);
  pp = addPropN(&ecl,"SEGMENTS",'i',3, bnd);
  if (!pp) 
	{
	  fprintf(stderr,"makeEclipseGrid(): Could not add SEGMENTS\n");exit(0);
	}
  vi = (int ***) (pp->p);
  for (i=1;i<=nx;i++)
	for (j=1;j<=ny;j++) 
	  for (k=1;k<=nz;k++) vi[i][j][k]=1;
  pp->set = 1;


  // Tidying as per loadModel()
  setAuxPointers(&ecl);
  /* Additional operations */
  if (transformToMap) (void) rotateGridMap(&ecl);
  applyPinch(&ecl); /*Applies pinchout criteria if pinch->set*/
  if (applyPinch(&ecl))
	{
	  /*Applies pinchout criteria if pinch->set*/
	  fprintf(stderr,"makeEclipseGrid(): Error applying requested pinchout");
	}
	  
  /* Clean up */
  compressProp(&ecl); 

  

  /* now we write to the eclOut file */
  FILE *outEcl = NULL;
  char* filename = massager->outputgrid[0].filename;

  outEcl = fopen(filename, "w");
  if (!outEcl) 
	{
	  fprintf(stderr,"makeEclipseGrid():: Could not open output file %s, exiting\n", filename);
	}
  if (writeRecordList (outEcl, &ecl, nKeys , keyWords)) 
	{ 
	  fprintf(stderr, "makeEclipseGrid(): Error encountered in writeRecordList\n");
	}
  // close and flush, as read follows very soon
  char synced;
  synced = (fsync(fileno(outEcl)) == 0) ? 'y' : 'n';
  fclose(outEcl);

  printf("Wrote (synced=%c) 'receptacle' ECLIPSE file %s\n\n\n",synced,filename);

  if (massager->utmInfo.performInputGridUTMmapping)
	{
	  free_dvector(x,0,n-1); x=NULL;
	  free_dvector(y,0,n-1); y=NULL;
	}

  return EXIT_SUCCESS; 
}

int getEclipseGriddingFromBoundingRectangle(int n, double *x,double *y, struct NEAREST_REGULAR_GRID_INFO *RGinfo)
{
  int i,j,status;
  int size,*list;
  double **r=NULL;
  double u,v,block_size=0;
  double* corners=NULL;
  
  // copy xys into coords usable by twoDhull()

  r = dmatrix(0,n-1,0,1);
  for (j=0; j < n; j++) 
	{
	  r[j][0] = x[j];
	  r[j][1] = y[j];
	}
  
  list = twoDhull(n,(void**)r,'d',0,&size,1,&corners);

#if 0
  // Show plan in ascii file
  FILE *outMBR = NULL;
  char* filename = "EclipseMBRPlan.txt";
  outMBR = fopen(filename, "w");
  
  fprintf(outMBR,"#m=2,S=0\n");
  for (i=0;i<10;i+=2) fprintf(outMBR,"%lf %lf\n",corners[i%8],corners[(i+1)%8]);
  fprintf(outMBR,"\n");
  fprintf(outMBR,"#m=0,S=1\n# Points\n");
  for (i=0;i<n;i++) fprintf(outMBR,"%lg %lg\n",x[i],y[i]);
  fclose(outMBR);
  fprintf(stderr,"Wrote %s\n",filename);
#endif

  // corners now has smallest bounding rectangle, in anticlockwise order: adjacent 
  // indices give x,y's:   
  //       ^
  //       | vhat
  //   6,7  -------  4,5
  //       |       |
  //   0,1  -------  2,3    --> uhat
  u = hypot(corners[2]-corners[0],corners[3]-corners[1]);
  v = hypot(corners[6]-corners[0],corners[7]-corners[1]);

  RGinfo->u_min =0;
  RGinfo->u_max =0;
  RGinfo->v_min =0;
  RGinfo->v_max =0;

  if ((u>0) && (v>0))
	{
	  if (((v/u) < MAX_ECLIPSE_GRID_MAKING_ASPECT_RATIO) && ((u/v) < MAX_ECLIPSE_GRID_MAKING_ASPECT_RATIO))
		{
		  // not skinny: some disordered set....
		  block_size = sqrt(u*v/n);			
		}
	  else if ((u/v) >= MAX_ECLIPSE_GRID_MAKING_ASPECT_RATIO)
		{
		  // skinny in u direction
		  block_size = u/n;
		}
	  else if ((v/u) >= MAX_ECLIPSE_GRID_MAKING_ASPECT_RATIO)
		{
		  // skinny in v direction
		  block_size = v/n;
		}
	  RGinfo->Nu = MAX((int)(u/block_size + 0.5),1);
	  RGinfo->Nv = MAX((int)(v/block_size + 0.5),1);
	  RGinfo->du = u/RGinfo->Nu;
	  RGinfo->dv = v/RGinfo->Nv;
	  RGinfo->uhat[0] = (corners[2]-corners[0])/u;
	  RGinfo->uhat[1] = (corners[3]-corners[1])/u;
	  RGinfo->vhat[0] = (corners[6]-corners[0])/v;
	  RGinfo->vhat[1] = (corners[7]-corners[1])/v;
 	}
  else if ((v<=0) && (u>0))
	{
	  // skinny in u direction
	  block_size = u/n;
	
	  RGinfo->Nu = MAX((int)(u/block_size + 0.5),1);
	  RGinfo->Nv = 1;
	  RGinfo->du = u/RGinfo->Nu;
	  RGinfo->dv = RGinfo->du;
	  RGinfo->uhat[0] = (corners[2]-corners[0])/u;
	  RGinfo->uhat[1] = (corners[3]-corners[1])/u;
	  RGinfo->vhat[0] = -RGinfo->uhat[1];
	  RGinfo->vhat[1] = RGinfo->uhat[0];
	  RGinfo->v_min = -0.5*RGinfo->dv;
	  RGinfo->v_max =  0.5*RGinfo->dv;
	}
  else if ((u<=0) && (v>0))
	{
	  // skinny in v direction
	  block_size = v/n;

	  RGinfo->Nu = 1;
	  RGinfo->Nv = MAX((int)(v/block_size + 0.5),1);
	  RGinfo->dv = v/RGinfo->Nv;
	  RGinfo->du = RGinfo->dv;
	  RGinfo->vhat[0] = (corners[6]-corners[0])/v;
	  RGinfo->vhat[1] = (corners[7]-corners[1])/v;
	  RGinfo->uhat[1] = -RGinfo->vhat[0];
	  RGinfo->uhat[0] = RGinfo->vhat[1];
	  RGinfo->u_min = -0.5*RGinfo->du;
	  RGinfo->u_max =  0.5*RGinfo->du;
	}
  else
	{	  
	  printf("Minimum bounding rectangle is degenerate, & should not be gridding in "
			 "getEclipseGriddingFromBoundingRectangle(): i.e. bug. FIX.\n");
	  status = EXIT_FAILURE;
	}

  RGinfo->u_min += RGinfo->uhat[0] * corners[0] + RGinfo->uhat[1] * corners[1];
  RGinfo->u_max += RGinfo->uhat[0] * corners[2] + RGinfo->uhat[1] * corners[3];
  RGinfo->v_min += RGinfo->vhat[0] * corners[0] + RGinfo->vhat[1] * corners[1];
  RGinfo->v_max += RGinfo->vhat[0] * corners[6] + RGinfo->vhat[1] * corners[7];

  status = EXIT_SUCCESS;

	  
	  
 cleanup:

  free(corners);
  free(list);
  free_dmatrix(r,0,n-1,0,1);

  return status;
}


int writeEclipseRealisation(struct MASSAGER *massager,int real_no,float ***r_properties)
{
  int i,j,k;
  char* filename, *prefix=NULL;
  struct ECL_OUT_INFO *eclOutInfo;  

   eclOutInfo = &(massager->outputgrid[0].eclOutInfo);

   // Set output filename. Strip backwards to first suffix

   // filename = eclOutInfo->outputFilename;
   // for (i=strlen(filename); filename[i]!='.' && i>0; i--);
   // prefix = (i>0) ? strndup(filename,i) : strdup(filename);

   // asprintf(&filename,"%s.realisation.%d.ecl",prefix,real_no);
   // free(prefix);
   prefix = eclOutInfo->outputFilename;
   for (i=strlen(prefix); prefix[i]!='.' && i>0; i--);
   // asprintf replacement of strndup()
   asprintf(&filename,"%.*s.realisation.%d.ecl", (int)((i>0)?i:strlen(prefix)), prefix, real_no);

   eclOutInfo->outputFilename_r = filename;
   eclOutInfo->ecl_out_r_fp = fopen(filename,"w");
   if (eclOutInfo->ecl_out_r_fp == NULL)
	 {
	   printf("Cannot open output file %s for writing\n",filename);
	   return EXIT_FAILURE; 
	 } 

   /***** For realisations ****/

   /* Now we set up keywords and file I/O to the eclOut file */
   /* truncate the list for what we have available in the smashed model */
   /* I think 8 is the minimum we need */
   /* Write out  the kriged properties plus their std deviations */
   eclOutInfo->nDuplicateKeys = 9;

   // plus layers if present
   if (eclOutInfo->haveLayers)  eclOutInfo->nDuplicateKeys++;

   eclOutInfo->nGridBlockKeys_r = massager->nPropertyNames;
   if (massager->thickness_index >= 0)
	 {
	   // omit thickness: it's a geometry property
	   // eclOutInfo->nGridBlockKeys_r--;
	 }
   if (massager->d_index >= 0)
	 {
	   // also omit depth: it's a geometry property
	   eclOutInfo->nGridBlockKeys_r--;
	 }

   eclOutInfo->nKeys_r = eclOutInfo->nDuplicateKeys + eclOutInfo->nGridBlockKeys_r;

   eclOutInfo->keyWords_r = 
	 (char**)malloc((size_t)((eclOutInfo->nKeys_r+1) * sizeof(char*)));

   /* duplicated keys */
   /* duplicated keys */
   if (massager->verbosity>2)
	 {
	   fprintf(stdout,"Keys written to Realisations ECLIPSE grids:\n");
	   fprintf (stdout, 
				"Key Number   Key Value                       \n"
				"------------ --------------------------------\n");
	 }
   char* s;
   for (j=1;j<=eclOutInfo->nDuplicateKeys; j++)
	 {
	   eclOutInfo->keyWords_r[j] = eclOutInfo->keyWords[j]; 
	   if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,eclOutInfo->keyWords_r[j]);				  
	 }

   /* concatenate properties. We'll write the kriged means plus std deviations */

   for (k=0; k<massager->nPropertyNames; k++)
	 {
	   /* skip depth as it is visible in Petrel */
	   if (!strcmp(massager->propertyNames[k],"d")) continue;

	   /* property */
	   eclOutInfo->keyWords_r[j] = massager->propertyNames[k];

	   if (massager->verbosity>2) fprintf(stdout,"%12d %32s\n",j,eclOutInfo->keyWords_r[j]);
	   j++;
	 }

   // To store block-centred values, use:
   //  massager->eclInfo.property[p][i][j][k]
   // first index is property, then regular contiguous cubes i,j,k */

   // r_properties should behave like KrigMeans.....
   // ECLIO I/O
   // realisation data data is already fully 'accumulated'

   writeECLIO_InSituOutput(real_no,0,0,massager,r_properties,NULL);

   fclose(eclOutInfo->ecl_out_r_fp);
   free(filename);

   return EXIT_SUCCESS;
 }

 struct SEG_ENCODING **fetchSegLists(ECL *ecl,int *nNodes, int verbosity)
 {
   PROP *aP = NULL;
   PROP *actP= NULL;
   int ***a = NULL;
   int **a_act = (int **) NULL;
   int ***list= NULL;
   char routine[] = "fetchSegLists";
   char print;
   int i,ii,j,jj,k,kk,kp,found;
   int segCount,size,c,seg,status;
   struct INTCHAR *pr;

   int nx = ecl->specGrid.nx;
   int ny = ecl->specGrid.ny;
   int nz = ecl->specGrid.nz;
   int nrow = ecl->specGrid.nx+1;
   int ncol=ecl->specGrid.ny+1;

   int node_cnt=0;

   struct SEG_ENCODING **segInfo;

   // allocate segments pointer structure segInfo[i][j]{size,values[0..size-1],}

   /* allocate pointers to rows */
   segInfo=(struct SEG_ENCODING **) malloc((size_t)(nrow*sizeof(struct SEG_ENCODING*)));
   if (!segInfo) 
	 {
	   fprintf(stderr,("allocation failure 1 in fetchSegLists()")); 
	   return NULL;
	 };

   /* allocate rows and set pointers to them */
   segInfo[0]=(struct SEG_ENCODING *) malloc((size_t)((nrow*ncol)*sizeof(struct SEG_ENCODING)));
   if (!segInfo[0]) 
	 {
	   fprintf(stderr,("allocation failure 2 in fetchSegLists()")); 
	   return NULL;
	 }

   for (i=1; i<=ecl->specGrid.nx; i++) segInfo[i] = segInfo[i-1]+ncol;


   /* Get ActCoord */
   actP = prepareActCoord(ecl,0);
   if (!actP || !actP->set) return NULL;
   a_act = (int **)actP->p;
   if (!a_act) return NULL;




   /* Get segments */
   aP = fetchProp(ecl,"SEGMENTS");
   if (!aP || !aP->set || !aP->p) 
	 {
	   status = 3;
	   eclError(routine, status, "Must set SEGMENTS before calling\n");
	   return NULL;
	 }
   a = (int ***) aP->p;

   /* the actual calc */

   // show_malloc_stats(); 

   for (i=0;i<=nx;i++) 
	 {
	   for (j=0;j<=ny;j++) 
		 {
		   GPtrArray *segArray = g_ptr_array_new();
		   segCount=0;
		   // loop over adjoining block columns
		   for (jj=j,c=4; jj<=j+1; jj++)
			 {
			   for (ii=i; ii<=i+1; ii++,c--)
				 {
				   if (ii>=1 && ii<=nx && jj>=1 && jj<=ny)  // sensible segment indices
					 {
					   for (k=1;k<=nz;k++)                  // grid blocks
						 {
						   /* I think we should take the segment 
							  number of the non-negative labels
							  if a mixture is supplied, e.g. {-999,-999,2,2,2,2,}, 
							  assuming -ve means (don't know|too lazy)
						   */
						   if (a[ii][jj][k]>=0)             // useful segment label 
							 {
							   search_in_seg_array(segArray,a[ii][jj][k],c,&segCount); 
							 }
						 }
					 }
				 }
			 }

		   // count number of distinct segments

		   // we now have size entries
		   segInfo[i][j].size = segCount;
		   if (segCount>4)
			 {
			   printf("WARNING!!!: more than 4 segments meet at trace %d %d! Truncating to max of 4.\n",i,j);
			   segCount = 4;
			 }
		   if (a_act[i][j])   
			 {
			   if (segCount==0)
				 {
				   printf("FATAL ERROR: no sensible bounding segment labels at active trace indexed %d %d!\n",i,j);
				   // exit(0);
				 }
			   // count nodes, including duplicates
			   if (segCount > 0) node_cnt += segCount;
			 }
		   if (segCount > 0)
			 {
			   print = (char)((segCount>1 && verbosity>3) ? 1 : 0);
			   if (print)
				 {
				   printf("Duplicates at i,j= %d, %d, (segments,c): ",i,j);
				 }
			   segInfo[i][j].values = ivector(0, segCount-1);			  
			   segInfo[i][j].corner_encoding = cvector(0, segCount-1);	
			   for (k=0; k< MIN(segArray->len,4); k++)
				 { 
				   pr  = g_ptr_array_index(segArray,k);
				   segInfo[i][j].values[k] = pr->indx;
				   segInfo[i][j].corner_encoding[k] = pr->c; // Anded 'c' value, i.e. the sum (1<<c)
				   if (print) 
					 {
					   printf("(%d,%d) ",segInfo[i][j].values[k],(int)(segInfo[i][j].corner_encoding[k]));
					 }
				 }
			   if (print) printf("\n");
			 }
		   else 
			 {
			   segInfo[i][j].values = NULL;
			   segInfo[i][j].corner_encoding = NULL;
			 }
		   // clean up
		   for (k=0; k< segArray->len; k++) free(g_ptr_array_index(segArray,k));
		   g_ptr_array_free(segArray,TRUE);
		 }
	 }
   *nNodes = node_cnt;

   // show_malloc_stats(); 

   return segInfo;
 }


 void search_in_seg_array(GPtrArray *segArray,int seg, int c, int *segCount)
 {
   int i;
   int append, newSeg;
   struct INTCHAR *pr;
   char cn;

   // assume we have to add a new seg until contradicted
   newSeg=1;
   // linear search
   for (i=0; i< segArray->len; i++)
	 {
	   pr = g_ptr_array_index(segArray,i);
	   if (pr->indx == seg)
		 {
		   // segment matches: don't increment later
		   newSeg=0;
		   cn = (1<<(int)c);
		   // have we not seen 'c' yet?
		   if (!(pr->c & cn))
			 {
			   // Add it in-place.
			   pr->c = (char)(pr->c | cn);
			 }
		   break;
		 }
	 }
   if (newSeg) 
	 {
	   (*segCount)++;
	   // append entries
	   pr = (struct INTCHAR *)malloc((size_t)sizeof(struct INTCHAR));
	   pr->indx = seg;
	   pr->c = (1<<(int)c);
	   g_ptr_array_add(segArray,pr);
	 }
 }

void free_segInfo(struct MASSAGER *massager,struct SEG_ENCODING **segInfo,int Nx,int Ny)
 {
   int i,j;

   if (massager->eclInfo.useECLIOsegments)
	 {
	   for (i=0;i<=Nx;i++)
		 {
		   for (j=0;j<=Ny;j++)
			 {
			   if (segInfo[i][j].size>0)
				 {
				   free_ivector(segInfo[i][j].values,0,segInfo[i][j].size-1);
				   segInfo[i][j].values = NULL;			  
				   free_cvector(segInfo[i][j].corner_encoding,0,segInfo[i][j].size-1);
				   segInfo[i][j].corner_encoding=NULL;
				 }
			 }
		 }

	   free(segInfo[0]);
	   segInfo[0]=NULL;
	   free(segInfo);
	   segInfo = NULL;
	 }
 }

 void allocateStochEclipseMem(struct MASSAGER *massager)
 {
   int i,j;
   struct ECL_OUT_INFO *eclOutInfo;

   eclOutInfo = &(massager->outputgrid[0].eclOutInfo);
   eclOutInfo->r_properties = (float ****)malloc((size_t) (massager->Nrealisations * sizeof(float ***)));
   eclOutInfo->r_properties[0] = (float ***)malloc((size_t) (massager->Nrealisations * massager->nLayers * sizeof(float **)));
   for (i=1; i< massager->Nrealisations ; i++)
	 {
	   (eclOutInfo->r_properties)[i] = (eclOutInfo->r_properties)[i-1]+massager->nLayers;
	 }
   for (i=0; i< massager->Nrealisations ; i++)
	 {
	   for (j=0; j< massager->nLayers ; j++)
		 {
		   (eclOutInfo->r_properties)[i][j] = matrix(0,massager->outputgrid[j].nPoints-1,
													 0,massager->nPropertyNames-1);
		 }
	 }
 }

 void freeStochEclipseMem(struct MASSAGER *massager)
 {
   int i,j,Nx,Ny,Nz;
   struct ECL_OUT_INFO *eclOutInfo;
   PROP *from= NULL;
   ECL *ecl;

   eclOutInfo = &(massager->outputgrid[0].eclOutInfo);

   for (i=0; i< massager->Nrealisations ; i++)
	 {
	   for (j=0; j< massager->nLayers ; j++)
		 {
		   free_matrix(eclOutInfo->r_properties[i][j], 0,massager->outputgrid[j].nPoints-1,
					   0,massager->nPropertyNames-1);
		   eclOutInfo->r_properties[i][j]=NULL;
		 }
	 }
   free(eclOutInfo->r_properties[0]); eclOutInfo->r_properties[0] = NULL;
   free(eclOutInfo->r_properties); eclOutInfo->r_properties = NULL;

   from = fetchProp(&eclOutInfo->eclOut,"ZCORN");
   free_anyProp(from);

   ecl = &(massager->outputgrid[0].ecl);
   Nx = getNX(ecl);
   Ny = getNY(ecl);
   Nz = getNZ(ecl);
   eclOutInfo = &(massager->outputgrid[0].eclOutInfo);   
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

/*
  Scan over group of blocks centered on ecl_i,ecl_j to find one that encloses
  point (x,y) = testPt. Return eclipse segcol value if found, I_MISSING otherwise
*/
int getSegOfEclipseBlock(ECL *ecl,int **segcol, int ecl_i,int ecl_j,double z,double* testPt) 
{
  int i,j,ii,jj,inside=0;
  int range=2; 

  // range = range of indices i,j either side of ecl_i|j to examine. 2 seems 
  // a good compromise between efficiency and robustness
  int seg=I_MISSING;

  // inner 4 blocks have high hit rate. try these first
  for (i=MAX(0,ecl_i-1),ii=-1; i<=MIN(ecl->specGrid.nx-1,ecl_i+1); i++,ii++)
	{
	  for (j=MAX(0,ecl_j-1),jj=-1; j<= MIN(ecl->specGrid.ny-1, ecl_j+1); j++,jj++)
		{
		  inside = isInEclipseBlock(ecl,i,j,z,testPt);
		  if (inside) 
			{
			  seg = segcol[i+1][j+1];
			  // printf("SEGTEST %d %d\n",ii,jj);
			  return seg;
			}
		}
	}

  // if still no success, look in next shell
  for (i=MAX(0,ecl_i-range),ii=-range; i<=MIN(ecl->specGrid.nx-1,ecl_i+range); i++,ii++)
	{
	  for (j=MAX(0,ecl_j-range),jj=-range; j<= MIN(ecl->specGrid.ny-1, ecl_j+range); j++,jj++)
		{
		  // skip box inside
		  if (ii>=-1 && ii<=0 && jj>=-1 && jj<=0) continue;

		  inside = isInEclipseBlock(ecl,i,j,z,testPt);
		  if (inside) 
			{
			  seg = segcol[i+1][j+1];
			  // printf("SEGTEST %d %d\n",ii,jj);
			  return seg;
			}
		  // printf("SEGTEST %d %d %d\n",ii,jj,inside);
		}
	}
  return seg;
}


/*
  Test whether point testPt[2] (x,y) is in eclipse block with lower left corner indexed ecl_i,ecl_j, at this depth z.
  x,y's of quadrilateral corners are computed from coord values etc.
*/
int isInEclipseBlock(ECL *ecl,int ecl_i,int ecl_j,double z,double* testPt) 
{
  int i,j,k,dir,inside=0;
  double point[4];
  double pgon[4][2];
  double*** c;

  c = ecl->coord;
  
  // form square from i,j, corners, with z set from arg
  for(i = ecl_i,k=0; i<= ecl_i+1; i++)
	{
	  for (j=ecl_j; j<= ecl_j+1;j++)
		{ 
		  interpVert(z, &point[0], c[i][j]);
		  //  point[1,2,3]=x,y,z
		  k= (i==ecl_i) ? j-ecl_j : 3+ecl_j -j; // polygon in natural ordering (k will scan 0,1,3,2)
		  pgon[k][0] = point[1]; //x
		  pgon[k][1] = point[2]; //y
		}
	}
  /*
  for (k=0;k<4;k++)
	{
	  printf("%g %g\n",pgon[k][0] ,pgon[k][1]);
	}
  printf("%g %g\n\n",pgon[0][0] ,pgon[0][1]);
  */

  /* Used convex crossing algorithm for enclosure test */

  inside = CrossingsTest(pgon,4,testPt);

  return inside;
}

/* 
   Using pt_cnrs[0,1...11] = {x1,y1,z1, x2,y2,z2...}
   of 4 corners of a top face, compute coefficients to predict a new z
   from an x,y in r[0],r[1], by solving for bilinear_coeffs={A,B,C,D} 
   using the (vector) predictor
   z = A + Bu + Cv + D uv
   with u,v in [0,1] parameterising the top face
   Answer z returned in r[2].

   We solve for the A,B,C,D linearly, then the u,v of the new point (quadratic), then
   use the formula
   
*/
int doBilinearInterp(double* pt_cnrs,double* r,double *u,double *v)
{
 static double c[4][3]; 
 double det,minus_b;
 int i;

 // r = c[0] + u c[1] + v c[2] + u v c[3] 
 for (i=0;i<3;i++)
	{
	  c[0][i] = pt_cnrs[i];
	  c[1][i] = pt_cnrs[i+3] - pt_cnrs[i];
	  c[2][i] = pt_cnrs[i+6] - pt_cnrs[i];
	  c[3][i] = pt_cnrs[i+9] - pt_cnrs[i] -  c[1][i] -  c[2][i];
	}

 // Solve for u,v
 // determinant part
 det = (-4*(c[2][1]*c[3][0] - c[2][0]*c[3][1])*(c[0][1]*c[1][0] - c[0][0]*c[1][1] + c[1][1]*r[0] - c[1][0]*r[1]) + 
		 DSQR(c[1][1]*c[2][0] - c[1][0]*c[2][1] - c[0][1]*c[3][0] + c[0][0]*c[3][1] - c[3][1]*r[0] + c[3][0]*r[1]));

 if (det<0)
   {
	 // warning("Bilinear failed with -ve determinant");
	 return EXIT_FAILURE; 
   }

 // denominator
 *v = 2.0*(c[1][1]*c[3][0] - c[1][0]*c[3][1]);

 if ((*v)==0)
   {
	 warning("Bilinear failed with div. by zero");
	 // set for safety and return
	 *v = -1;
	 return EXIT_FAILURE; 
   }

 // else, must be OK
 det = sqrt(det);
 minus_b = -c[1][1]*c[2][0] + c[1][0]*c[2][1] - c[0][1]*c[3][0] + c[0][0]*c[3][1] - c[3][1]*r[0] + c[3][0]*r[1];
 // quadratic formula
 *u = (minus_b + det)/(*v);
 if ((*u)<0 || (*u)>1)
   {
	 // try other quadratic solution for luck...
	 *u = (minus_b - det)/(*v);
   }
 // use just x component to get v...
 (*v) = (r[0]-c[0][0]-(*u)*c[1][0])/(c[2][0]+(*u)*c[3][0]);
 
 // interpolate!
 r[2] = c[0][2]+ (*u) * c[1][2] + (*v) * c[2][2]+ (*u) * (*v) *c[3][2];


 return EXIT_SUCCESS; 
}

// debug func to draw grid top face and candidate point..
void debug_plotTestPt(double* pt_cnrs,double* r)
{
  int i,j;

  FILE *f;
  f = fopen("/tmp/junk.txt","w");
  for(i=0;i<5;i++)
	{
	  switch(i)
		{
		case 0:
		  j=0; break;
		case 1:
		  j=3; break;
		case 2:
		  j=9; break;
		case 3:
		  j=6; break;
		case 4:
		  j=0; break;
		}
	  fprintf(f,"%lg %lg\n",pt_cnrs[j],pt_cnrs[j+1]);
	}
  fprintf(f,"#m=0,S=2\n%lg %lg\n",r[0],r[1]);
  fclose(f);
  system("gr /tmp/junk.txt; ");
}
