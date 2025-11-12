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
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <glib.h>
extern int errno;


#include "ECLIO/eclio.h"  /* structure definition and a little miscellany */
#include "ECLIO/access.h" /* access high-level routines */
#include "ECLIO/eclutil.h" /* for model initialization only...my move to access */
#include "ECLIO/blockops.h" /* also added some lower-level routines, flagged below */
#include "ECLIO/copy.h"    /* ways to set records for moving data between models */
#include "ECLIO/read.h"    /* reads */


#include "commandline/commandline.h" 
#include "platform_specifics.h"

/****************************************************************************************

Simple code to print ascii-graphable files of eclipse ascii-grid
cross-sections.  x-values are taken as actual eclipse x,y's or as grid
integers ecl_i,ecl_j, and y-values are the zCorns of each block.  The
sequence of output numbers draws each grid trapezoid completely going
up columns, so piping straight to a line-graphing util makes a
suitable picture.

e.g.

% eclipseXsection x 5 eclipsefile.txt  | graph -TX

shows the xsection at i=5


*****************************************************************************************/



#define GRAY_0 0.2
#define GRAY_1 0.8
#define HEIGHT_INCHES 4
#define WIDTH_INCHES 6
#define USE_FILL_FUNC_OF_GRAPH 1

#define N_REQD_ECLIPSE_KEYWORDS 8
#define N_BONUS_ECLIPSE_KEYWORDS 2
static char* DESIRED_ECLIPSE_KEYWORDS[N_REQD_ECLIPSE_KEYWORDS+N_BONUS_ECLIPSE_KEYWORDS] 
     = {"PINCH","MAPAXES","GRIDUNIT","SPECGRID","COORDSYS","COORD","ZCORN","ACTNUM","SEGMENTS","NG"};

#define ARGS_USAGE "x|y|z indx eclipsefile.txt | stdout\n\nPrints x-section of eclipse file to stdout, at ecl_(i|j|k)=index, \
for x|y|z respectively\n\
For direct piped output with sensible colors, use, e.g. (plotutils):\n\
eclipseXsection x|y|z indx eclipsefile.txt | graph -TX -C -B --pen-colors '4=black:5=brown' --bitmap-size 1000x1000\n"

static  CLOption   helpOption = 
    { '\0', "help", 0, 0, "Prints x-section of eclipse file to stdout, at ecl_(i|j|k)=index, for x|y|z respectively. \
                           For direct piped output, use, e.g. (plotutils):\n \
				           graph -TX -C -B --pen-colors '4=black:5=brown' --bitmap-size 1000x1000\n"};


/* Options */
int useStdout=1; 
int drawRef=0;
int refN=0;
int printPS=0;
char* printPSfilename="GridXsection.ps";
float printPSsize=WIDTH_INCHES;
int setPrintPSsize=0;
int colourActives=0;
int skipInActives=0;
int showSegments=0;
int verbosity=0;

/* Handy globals */
int Nx,Ny,Nz;
int ***active=NULL;
int lattice=0;
int haveSegments=0; 
int minSeg=0,maxSeg=0;
int ***seg=NULL;
int **segcol=NULL;
int haveNG=0;
double*** NG=NULL;
ECL ecl;
double **** zCorn;
int xyzSupplied=0;
int iN,jN,kN;
float xN,yN,zN;
char* x_y_z=NULL;
int flipGridunits=0;

CLOption options[] =
  {
    {'v',NULL,1,  &verbosity, "verbosity",{{CL_INT_ARG,&verbosity,"level",NULL}}},
	{'a',NULL,0,  &colourActives, "colour grid by active flag"},
	{'\0',"skip-inactives",0,  &skipInActives, "don't draw inactive blocks"},
	{'s',NULL,0,  &showSegments, "colour grid by segments flag"},
	{'l',NULL,0,  &lattice, "display transverse cross sections using integer grid lattice, not natural x,y"},
	{'r',NULL,1,  &drawRef, "print top-surface only, for layer number reference-layer-number (to stdout)",{{CL_INT_ARG,&refN,"reference-layer-number",NULL}}},
	{'p',NULL,1,  &printPS, "print postscript image, with NG colouring if available",{{CL_STRING_ARG,&printPSfilename,"filename (or - for stdout)",
																					   NULL}}},
	{'\0',"point",1,  &xyzSupplied, "draw section nearest to supplied x,y,z point ",{{CL_STRING_ARG,&x_y_z,"x,y,z",NULL}}},
	{'\0',"pssize",1,  &setPrintPSsize, "postscript image width in inches",{{CL_FLOAT_ARG,&printPSsize,NULL,NULL}}},
	{'\0',"flip-gridunits",0,  &flipGridunits, "change input Eclipse grid from feet to meters"},
	CL_END_OPTIONS  /* defined as {0,NULL,..} to mark the end of the list   */
  };

void usage();
void printXsection(int iXsection,int nch, int ndh,double ****zCorn);
void printYsection(int jXsection,int nrh, int ndh,double ****zCorn);
void printLayerTop(int doX, int iXsection,int zIndx, int nrh, int nch,double ****zCorn);
void printPSCrossSection(ECL* ecl,int doY,int jXsection,int nrh, int ndh,double ****zCorn);
void printPlanSection(int kZsection,int Nx, int Ny,ECL *ecl,double ****zCorn);
void printPSPlanSection(int kZsection,int Nx, int Ny,ECL *ecl,double ****zCorn);

void getSegmentRange(ECL* ecl, int iS,int jS,int kS);
void printActives(int Nx, int Ny, int k,ECL *ecl);
int activateThinActnums(ECL *ecl); 
double getVertexXcoord(int i, int j, int k, int vrtx);
double getVertexYcoord(int i, int j, int k, int vrtx);
void getNearestBlockIndices(float X,float Y,float Z, int *i,int *j,int *k);
void drawSuppliedPoint(FILE* fp,double x,double y);
gint int_compare(gconstpointer a, gconstpointer b);
gboolean print_seg(gpointer key,gpointer value,gpointer data);

void usage()
{
  printf("Usage: eclipseXsection [-r n] x|y|z i_sec|j_sec eclipsefile.txt \n"); 
  printf("     -r draws reference layer top only (layer n)\n"); 
  printf("     -a colours by the active label\n"); 
  exit(1);
}

#ifndef SQR
#define SQR(a) ((a)*(a))
#endif

int main(int argc, char *argv[])
{ 

  CLContext   context;
  int i,j,k,ii,transformToMap,status=0,nKeys;
  int fnd=0;
  int doX=1,doZ=0; 
  int indx=0;
  float X,Y,Z;

  char * indx_ptr;
  char **keyWords = NULL;
  char msg[LINE_LENGTH];

  PROP *from=NULL;
  double ****z=NULL;
  PROP *aP = NULL; 

  char *filename;
  FILE *fp;


  clInit( &context, options,ARGS_USAGE, NULL );

  if (clParse( &context, argc, argv ))
    {
	  if (context.argc != 4) 
        {
		  fprintf(stderr, "\n");
		  fprintf(stderr, "Missing some command line arguments\n");
		  clPrintUsage(&context);		  
		  return EXIT_FAILURE;
        }
	  // x or y?
	  if (!strcmp(context.argv[1],"x")) doX=1;
	  else if (!strcmp(context.argv[1],"y")) doX=0;
	  else if (!strcmp(context.argv[1],"z")) doZ=1;
	  else clPrintUsage(&context);

	  // indx
	  indx = strtol(context.argv[2],&indx_ptr, 10);
	  if (errno == ERANGE || (indx_ptr==context.argv[2] && indx==0))
		{
		  fprintf(stderr,"Bad indx specification; require integer...\n"); 
		  usage();
		}
	  
	  // ecl file?
	  filename = context.argv[3];

	  // supplied point
	  if (xyzSupplied)
		{
		  char *delims=",";
		  char *token;
		  token = strtok(x_y_z,delims);
		  xN = atof(token);    /* convert first token to x */
		  token = strtok(NULL,delims);  /* get next token */
		  if (token == NULL ) {printf("specified point requires x,y,z as numbers\n");usage();} 
		  yN  = atof(token); 
		  token = strtok(NULL,delims);  /* get next token */
		  if (token == NULL ) {printf("specified point requires x,y,z as numbers\n");usage();} 
		  zN = atof(token);
		  fprintf(stderr,"Using point (%g,%g,%g)\n",xN,yN,zN);
		}
	}
  else return EXIT_FAILURE;


  

			  
  /* read Eclipse grids */
  transformToMap = 1; /* number of records to read and
						 force transform grid->map coords */
  if (verbosity>0) fprintf(stderr,"Loading Eclipse model from file %s\n",filename);

  // (1) simple read to check x,y's
  loadModelFromList(filename,&ecl,transformToMap,1,"SPECGRID");



  /* Eclipse dimensions */
  Nx = getNX(&ecl);
  Ny = getNY(&ecl);
  Nz = getNZ(&ecl); 
 

  if (doZ && indx>=Nz)
	{
	  fprintf(stderr,"Zindx value too high, must be in [0,%d]\n",Nz-1);
	  usage();
	}
  else
	{
	  if (doX && indx>=Nx)
		{
		  fprintf(stderr,"Xindx value too high, must be in [0,%d]\n",Nx-1);
		  usage();
		}
	  if (!doX && indx>=Ny)
		{
		  fprintf(stderr,"Yindx value too high, must be in [0,%d]\n",Ny-1);
		  usage();
		}
	  if (drawRef)
		{
		  if (refN<0 || refN>=Nz)
			{
			  fprintf(stderr,"ref layer too high, must be in [0,%d]\n",Nz-1);
			  usage();
			}
		}
	}

  // (2) Full read
  freeModel(&ecl);

  nKeys = fetchKeyList(filename, &keyWords, "--");
  /* List possible keywords */
  if (verbosity>0) 
	{
	  fprintf(stderr, 
			  "\nFound %d *probable* keys in file %s\n"
			  "Key Number   Key Value                        Message\n"
			  "------------ -------------------------------- ------- \n",
			  nKeys,filename);
	}
  for (j=1;j<=nKeys;j++) 
	{
	  checkKeyWord(keyWords[j],"(){}[]-_",msg,LINE_LENGTH);
	  if (verbosity>0) fprintf(stderr,"%12d %32s %s\n",j,keyWords[j], msg);
	}

  // reduce list to useful ones...
  // Mandatory words
  for (j=1,k=1;j<=nKeys;j++) 
	{ 
	  fnd=0;
	  for (i=0;i<N_REQD_ECLIPSE_KEYWORDS;i++)
		{
		  if (!strcmp(keyWords[j],DESIRED_ECLIPSE_KEYWORDS[i])) {
			fnd=1;break;
		  }
		}
	  if (fnd) keyWords[k++]=keyWords[j];
	}

  // Bonus words
  for (j=1;j<=nKeys;j++) 
	{
	  fnd=0;
	  for (ii=0,i=N_REQD_ECLIPSE_KEYWORDS;ii<N_BONUS_ECLIPSE_KEYWORDS; ii++,i++)
		{
		  // implied hack around Petrel's ornamented SEGMENTS word...
		  if (!strncmp(keyWords[j],DESIRED_ECLIPSE_KEYWORDS[i],
					   strlen(DESIRED_ECLIPSE_KEYWORDS[i]))) {
			fnd=1;break;
		  }
		}
	  if (fnd)
		{
		  if (ii==0)  haveSegments = 1;
		  if (ii==1)  haveNG = 1;
		  keyWords[k++]=keyWords[j];
		}
	}
  nKeys=k-1;

  if (showSegments && !haveSegments)
	{
	  fprintf(stderr,"SEGMENTS not available for colouring, exiting\n");
	  return EXIT_FAILURE;
	}

  /* After a bit of cleaning -- if you get a message here, program will probably crash 
	 when loading ... will provide a cushier exit later*/
  nKeys = editKeyList(nKeys, &keyWords);
  if (verbosity>0) 
	{
	  fprintf (stderr, 
		   "\nAfter editKeyList: %d *used* keys in file %s\n"
		   "Key Number   Key Value                        Message\n"
		   "------------ -------------------------------- ------- \n",
		   nKeys,filename);
	}
 for (j=1;j<=nKeys;j++) 
	{
	  checkKeyWord(keyWords[j],"(){}[]-_",msg,LINE_LENGTH);
	  if (verbosity>0) fprintf(stderr,"%12d %32s %s\n",j,keyWords[j], msg);
	}

  /* General method */
  /* load model */
  status = loadModelFromArraySansPinch(filename, &ecl, transformToMap, nKeys, keyWords);  
  if (status) 
	{
	  fprintf(stderr,"ERROR: Could not load model, exiting\n");
	  return EXIT_FAILURE;
	}
  
  // Actives
  aP = fetchProp(&ecl, "ACTNUM");
  if (!aP || !aP->p || !aP->set)
	{
	  fprintf(stderr,"Could not get ACTNUM\n");
	  return EXIT_FAILURE;
	}
  active = (int ***) aP->p;


  // Segments
  if (haveSegments)
	{
	  segcol = fetchSegCol(&ecl, 1, Nx, 1, Ny);
	}
  if (haveNG)
	{
	  from = fetchProp(&ecl,"NG");
	  NG = ((double ***)(from->p));
	}


  from = fetchProp(&ecl,"ZCORN");
  z = ((double ****)(from->p));  // should be eclOut->zCorn
  zCorn = z;

  // metricify if required
  if (flipGridunits) flipGridAndMapUnits(&ecl);

  // user supplied point
  if (xyzSupplied)
	{
	  getNearestBlockIndices(xN,yN,zN,&iN,&jN,&kN);
	  if (doZ) indx=kN-1;
	  else if (doX) indx = iN-1;
	  else indx = jN-1; // doY=1
	  
	}

  // Shading feedback
  if (colourActives) fprintf(stderr,"Colouring by actives\n");
  else if (showSegments && haveSegments) fprintf(stderr,"Colouring by segments\n");
  else if (haveNG) fprintf(stderr,"Colouring by NG\n");
  
  if (doZ)
	{
	  if (printPS)
		{
		  printPSPlanSection(indx,Nx,Ny,&ecl,z); 
		}
	  else printPlanSection(indx,Nx,Ny,&ecl,z); 
	}
  else
	{
	  if (drawRef)
		{
		  printLayerTop(doX,indx,refN,Nx,Ny,z);
		}
	  else 
		{
		  // main grids
		  if (doX) 
			{
			  if (printPS)
				{
				  printPSCrossSection(&ecl,0,indx,Ny,Nz,z);
				}
			  else 
				{
				  getSegmentRange(&ecl,indx+1,-1,-1);
				  printXsection(indx,Ny,Nz,z); 
				}
			}
		  else 
			{
			  if (printPS)
				{
				  printPSCrossSection(&ecl,1,indx,Nx,Nz,z);
				}
			  else 
				{
				  getSegmentRange(&ecl,-1,indx+1,-1);
				  printYsection(indx,Nx,Nz,z);
				}
			}
		}
	}

  fprintf(stderr, "Finished.\n");

  clDone( &context );
  return 0;
}

void printActives(int Nx, int Ny, int k,ECL *ecl)
{
  int i,j;
  int ***active=NULL;
  PROP *aP = NULL; 

  aP = fetchProp(ecl, "ACTNUM");
  if (!aP || !aP->p || !aP->set)
	{
	  fprintf(stderr,"Could not get ACTNUM\n");
	  exit(-1);
	}

  // printf("ACTIVES\n");

  active = (int ***) aP->p;
  for (i=1; i<=Nx; i++)
	{
	  for (j=1; j<=Ny; j++)
		{
		  printf("%d %d %d\n",i,j,active[i][j][k]);
		}
	}	  
}

void printPlanSection(int kZsection,int Nx, int Ny,ECL *ecl,double ****zCorn)
{
  int i,j,iB,jB,kB,m,segLast=0,segment,pass;
  float fill;
  FILE *fpX;
  char *filename;

  double ***xy;

  kB = kZsection+1;

  xy = fetchLayerTraces(ecl,kB,1, 0,Nx,0,Ny);

  getSegmentRange(ecl,-1,-1,kB);

  // asprintf(&filename,"GridPlanSection_i=%d.txt",kZsection);
  //  printf("Writing const ecl_i Xsection at ecl_i=%d to %s\n",kZsection,filename);

  fpX= stdout;

  // simple Xsection check: pipe 
  // stdout | awk '/XSECEND/{p=0}p==1{print $0}/XSECSTART/{p=1}' | gr -C --bitmap-size 1500x500 -B &

  // for drawing simplicity, draw like this   CABD:
  //        A --B
  //        |   |
  //        C---D 

  m=1;
  for (pass=0;pass<2;pass++)
	{
	  for (i=0; i<Nx; i++) 
		{
		  for (j=0; j<Ny; j++)
			{		  
			  if (skipInActives && active[i+1][j+1][kZsection+1]==0) continue;
			  if (colourActives)
				{
				  fill = 1.0;
				  m = (active[i+1][j+1][kZsection+1]) ? 2 : 1;
				}
			  else if (showSegments)
				{  
				  if (haveSegments)
					{
					  m=1;
					  segment = seg[i+1][j+1][kB];
					  fill = 0.1+0.9*(MAX(-1,segment)-minSeg)/((float)(maxSeg-minSeg));
					}
				} 
			  else if (haveNG) 
				{
				  // YELLOWS - BROWN 
				  m=5;
				  fill = MAX(0,NG[i+1][j+1][kB]);
				  fill = 1-fill;  // if plotting on brown scale
				} 
			  else
				{
				  fill=0.2;
				}
			  if (pass==1) 
				{
				  m=4;fill=-1;
				  // hilight nearest lines 
				  if (xyzSupplied)
					{
					  if (((i+1)==iN) || ((j+1)==jN)) m=3;
					}
				} // 4=black,5=brown

			  fprintf(fpX,"#m=%d,S=0,q=%f\n",m,fill); // fill from desired property
			  fprintf(fpX,"%lf %lf\n",xy[i][j][1], xy[i][j][2]);
			  fprintf(fpX,"%lf %lf\n",xy[i][j+1][1], xy[i][j+1][2]);
			  fprintf(fpX,"%lf %lf\n",xy[i+1][j+1][1], xy[i+1][j+1][2]);
			  fprintf(fpX,"%lf %lf\n",xy[i+1][j][1], xy[i+1][j][2]);
			  fprintf(fpX,"%lf %lf\n\n",xy[i][j][1], xy[i][j][2]);
			}
		}
	}
  drawSuppliedPoint(fpX,xN,yN);
}


void printPSPlanSection(int kZsection,int Nx, int Ny,ECL *ecl,double ****zCorn)
{
  int i,j,iB,jB,kB,ii,pass,segLast=0;
  int width,height;
  int segment;
  float ymin,ymax,xmin,xmax,col;
  PROP* aP=NULL;
  FILE *fpX;
  fpos_t pos;
  char *filename=NULL;


  double ***xy;

  xy = fetchLayerTraces(ecl, kZsection+1,1, 0,Nx,0,Ny);

  kB = kZsection+1;

  // for drawing simplicity, draw cells like this   CABD:
  //        A --B
  //        |   |
  //        C---D 

  getSegmentRange(ecl,-1,-1,kB);

  // find ranges
  ymin=xmin=DBL_MAX;
  ymax=xmax=-DBL_MIN;

  for (i=0; i<=Nx; i++) 
	{
	  for (j=0; j<=Ny; j++)
		{
		  ymin=MIN(ymin, xy[i][j][2]);
		  ymax=MAX(ymax, xy[i][j][2]);
		  xmin=MIN(xmin, xy[i][j][1]);
		  xmax=MAX(xmax, xy[i][j][1]);
		}
	}

  if(!strcmp(printPSfilename,"-")) fpX = stdout;
  else
	{
	  filename = printPSfilename;
	  fpX = fopen(filename,"w");
	}
  
  // drawing 
  width=(printPSsize/WIDTH_INCHES)*WIDTH_INCHES*72;
  height=(printPSsize/WIDTH_INCHES)*HEIGHT_INCHES*72;

  float xScale = width/(xmax-xmin);
  float yScale = height/(ymax-ymin);

  /* PostScript header */
  fprintf(fpX,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fpX,"%%%%Creator: eclipseXsection.c (James Gunning)\n");
  fgetpos(fpX,&pos);
  fprintf(fpX,"%%%%BoundingBox: %4d %4d %4d %4d\n",-1,1,width+1,height+1);
  fprintf(fpX,"                           \n"); // protection against fsetpos hack below
  fprintf(fpX,"%%%%EndComments\n");
  /* The +72 shifts the figure one inch from the lower left corner */
  //fprintf(fpX,"%d %d translate\n", 72, 72 );
  fprintf(fpX,".05 setlinewidth\n");
 
	
  fprintf(fpX,"0.0 setgray\n");

  for (pass=0;pass<=1;pass++)
	{
	  if (pass==1) 
		{
		  fprintf(fpX,"\n\n%%%%Grid boxes from here till Showpage\n\n");
		  fprintf(fpX,"1.0 setgray\n");
		  fprintf(fpX,"0 0 0 setrgbcolor\n");
		}
	  for (i=0; i< Nx  ; i++) 
		{
		  for (j=0; j<Ny ; j++)
			{	
			  if (skipInActives && active[i+1][j+1][kB]==0) continue;

			  fprintf(fpX,"newpath\n");

			  fprintf(fpX,"%lf %lf moveto\n",xScale * (xy[i][j][1] -xmin), yScale * (xy[i][j][2] - ymin));
			  fprintf(fpX,"%lf %lf lineto\n",xScale * (xy[i][j+1][1] -xmin), yScale * (xy[i][j+1][2] - ymin));
			  fprintf(fpX,"%lf %lf lineto\n",xScale * (xy[i+1][j+1][1] -xmin), yScale * (xy[i+1][j+1][2] - ymin));
			  fprintf(fpX,"%lf %lf lineto\n",xScale * (xy[i+1][j][1] -xmin), yScale * (xy[i+1][j][2] - ymin));
			  fprintf(fpX,"%lf %lf lineto\n",xScale * (xy[i][j][1] -xmin), yScale * (xy[i][j][2] - ymin));
			  fprintf(fpX,"closepath\n");
				  
			  if (pass==0) 
				{
				  if (colourActives) {
					// RED GREEN
					if (active[i+1][j+1][kB]) {
					  fprintf(fpX,"0 1 0 setrgbcolor\nfill\n");
					} else fprintf(fpX,"1 0 0 setrgbcolor\nfill\n");
				  } else if (showSegments) {
					// BLUE for -999's, yellows for others
					segment = seg[i+1][j+1][kB];
					col = 0.5+0.5*(MAX(-1,segment)-minSeg)/((float)(maxSeg-minSeg));
					if (segment>=0) {
					  fprintf(fpX,"%f %f 0 setrgbcolor\nfill\n",col,col);
					} else fprintf(fpX,"0 0 1 setrgbcolor\nfill\n");
				  } else if (haveNG) {
					// YELLOWS - BROWN 
					col = NG[i+1][j+1][kB];
					fprintf(fpX,"%f %f 0 setrgbcolor\nfill\n",col,col);
				  } else {
					col=0.9;
					fprintf(fpX,"%f %f 0 setrgbcolor\nfill\n",col,col);
				  }
				  // fprintf(fpX,"%f setgray\nfill\n",(haveNG)? NG[iB][jB][j] : GRAY_1);
				}
			  else fprintf(fpX,"stroke\n");
			}
		}
	}
  fprintf(fpX,"showpage\n%%%%EOF\n");
  
  /*
	fsetpos(fp,&pos);
	fprintf(fp,"%%%%BoundingBox: %4d %4d %4d %4d\n",(int)floor(bxmin-1),(int)floor(bymin-1),
	(int)ceil(bxmax+1),(int)ceil(bymax+1));
  */
  fclose(fpX);
  printf("Wrote %s\n",filename);
}



void printXsection(int iXsection,int nch, int ndh,double ****zCorn)
{
  int i,j,iB,jB,kB,segLast=0,segment;
  FILE *fpX;
  char *filename;

  asprintf(&filename,"GridXsection_i=%d.txt",iXsection);
  //  printf("Writing const ecl_i Xsection at ecl_i=%d to %s\n",iXsection,filename);

  fpX= stdout;

  // simple Xsection check: pipe 
  // stdout | awk '/XSECEND/{p=0}p==1{print $0}/XSECSTART/{p=1}' | gr -C --bitmap-size 1500x500 -B &
  iB = iXsection+1;

  // for drawing simplicity, draw like this   CABD:
  //        A --B
  //        |   |
  //        C---D

#if USE_FILL_FUNC_OF_GRAPH

  // Attempt to use "fill fraction" facility of plotutils' graph.
  // can use in color mode, m=1,2,3,4,5. Reserve 5 for lines

  int pass,color,linecol=5;
  float fill=1.0;
  linecol=4;  //black
  color=1;    //red
  if (!colourActives && !showSegments && haveNG) color=5; // brown
  for (j=0,jB=1; j<nch; j++,jB++) 
	{
	  for (kB=1; kB<=ndh ; kB++)
		{
		  if (skipInActives && active[iB][jB][kB]==0) continue;
		  if (colourActives) {
			color = active[iB][jB][kB] ? 2 : 1;
		  } else if (showSegments) {
			segment = seg[iB][jB][kB];
			fill = 0.1+0.9*(MAX(-1,segment)-minSeg)/((float)(maxSeg-minSeg));
		  } else if (haveNG) {
			fill = 1-NG[iB][jB][kB];
		  } else {
			color = 3; // light blue defaults
			fill=0.2;
		  }
		  for (pass=0;pass<=1;pass++)
			{		  
			  fprintf(fpX,"#m=%d,S=0,q=%f\n", pass ? linecol : color, pass ? -1.0 : fill); 
			  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,5),zCorn[iB][jB][kB][5]);
			  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,1),zCorn[iB][jB][kB][1]);
			  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,3),zCorn[iB][jB][kB][3]);
			  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,7),zCorn[iB][jB][kB][7]);
			  fprintf(fpX,"%f %lf\n\n",getVertexYcoord(iB,jB,kB,5),zCorn[iB][jB][kB][5]);
			}
		}
	}

#else

  for (j=0,jB=1; j<nch; j++,jB++) 
	{
	  if (haveSegments && (segcol[iB][jB] != segLast))
		{
		  fprintf(fpX,"#m=%d,S=0\n",segcol[iB][jB]);
		  segLast = segcol[iB][jB];
		}
	  for (kB=1; kB<=ndh ; kB++)
		{		  
		  if (skipInActives && active[iB][jB][kB]==0) continue;
		  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,5),zCorn[iB][jB][kB][5]);
		  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,1),zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,3),zCorn[iB][jB][kB][3]);
		  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,7),zCorn[iB][jB][kB][7]);
		  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,kB,5),zCorn[iB][jB][kB][5]);
		}
	  if (skipInActives && active[iB][jB][ndh]==0) continue;
	  // go across to right edge...
	  fprintf(fpX,"%f %lf\n",getVertexYcoord(iB,jB,ndh,7),zCorn[iB][jB][ndh][7]);
	}
#endif

  drawSuppliedPoint(fpX,(lattice ? (double)jN : yN), zN);
}

void printYsection(int jXsection,int nrh, int ndh,double ****zCorn)
{
  int i,iB,jB,kB,segLast=0;
  FILE *fpX;
  char *filename;

  asprintf(&filename,"GridXsection_j=%d.txt",jXsection);
  // printf("Writing const ecl_j Xsection at ecl_j=%d to %s\n",jXsection,filename);

  // fpX = fopen(filename,"w");
  fpX = stdout;

  // simple Xsection check: pipe 
  // stdout | awk '/XSECEND/{p=0}p==1{print $0}/XSECSTART/{p=1}' | gr -C --bitmap-size 1500x500 -B &
  jB = jXsection+1;


#if USE_FILL_FUNC_OF_GRAPH

  // Attempt to use "fill fraction" facility of plotutils' graph.
  // can use in color mode, m=1,2,3,4,5. Reserve 5 for lines

  int pass,color=1,linecol=4,segment;
  float fill=1.0;

  if (!colourActives && !showSegments && haveNG) color = 5;

  for (i=0,iB=1; i<nrh; i++,iB++) 
	{
	  for (kB=1; kB<=ndh ; kB++)
		{		  
		  if (skipInActives && active[iB][jB][kB]==0) continue;
		  if (colourActives) {
			color = active[iB][jB][kB] ? 2 : 1;
		  } else if (showSegments) {
			segment = seg[iB][jB][kB];
			fill = 0.1+0.9*(MAX(-1,segment)-minSeg)/((float)(maxSeg-minSeg));
		  } else if (haveNG) {
			fill = 1-NG[iB][jB][kB];
		  } else {
			color = 3;
			fill=0.2;
		  }
		  for (pass=0;pass<=1;pass++)
			{
			  fprintf(fpX,"#m=%d,S=0,q=%f\n", pass ? linecol : color, pass ? -1.0 : fill); 
			  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,5), zCorn[iB][jB][kB][5]);
			  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,1), zCorn[iB][jB][kB][1]);
			  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,2), zCorn[iB][jB][kB][2]);
			  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,6), zCorn[iB][jB][kB][6]);
			  fprintf(fpX,"%f %lf\n\n",getVertexXcoord(iB,jB,kB,5), zCorn[iB][jB][kB][5]);
			}
		}
	}

#else

  for (i=0,iB=1; i<nrh; i++,iB++) 
	{
	  if (haveSegments && (segcol[iB][jB] != segLast))
		{
		  fprintf(fpX,"#m=%d,S=0\n",segcol[iB][jB]);
		  segLast = segcol[iB][jB];
		}
	  for (kB=1; kB<=ndh ; kB++)
		{		  
		  if (skipInActives && active[iB][jB][kB]==0) continue;
		  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,5), zCorn[iB][jB][kB][5]);
		  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,1), zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,2), zCorn[iB][jB][kB][2]);
		  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,6), zCorn[iB][jB][kB][6]);
		  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,kB,5), zCorn[iB][jB][kB][5]);
		}
	  if (skipInActives && active[iB][jB][ndh]==0) continue;
	  // go across to right edge... 
	  fprintf(fpX,"%f %lf\n",getVertexXcoord(iB,jB,ndh,5), zCorn[iB][jB][ndh][6]);
	}
#endif
  // fclose(fpX);
  drawSuppliedPoint(fpX,(lattice ? (double)iN : xN), zN); 
}

/* pass ySec=1 if doing Ysexction, 0 otherwise */

void printPSCrossSection(ECL* ecl, int ySec,int indxXsection,int Nmax,int ndh,double ****zCorn)
{
  int i,iB=0,jB=0,kB,dk,ii,pass,segLast=0,segment;
  int width,height;
  float zmin,zmax;
  float col[3]={0,0,0};
  FILE *fpX;
  fpos_t pos;
  char *filename=NULL;
  double x,y;
  float zScale,xScale,z,xmin,xmax,ymin,ymax;


  if (ySec)
	{
	  jB = indxXsection+1;
	  getSegmentRange(ecl,-1,jB,-1);
	}
  else
	{
	  iB = indxXsection+1;
	  getSegmentRange(ecl,iB,-1,-1);
	}
  // find x,y ranges
  ymin=xmin=zmin=DBL_MAX;
  ymax=xmax=zmax=-DBL_MIN;


  for (ii=0; ii< Nmax ; ii++) 
	{
	  if (ySec) iB=(ii+1); else jB=(ii+1);
	  for (kB=1; kB<=ndh ; kB++)
		{		  
		  if (skipInActives && active[iB][jB][kB]==0) continue;

		  z = zCorn[iB][jB][kB][1]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
		  z = zCorn[iB][jB][kB][5]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
		  if (ySec) 
			{
			  z = zCorn[iB][jB][kB][2]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
			  z = zCorn[iB][jB][kB][6]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
			  x = getVertexXcoord(iB,jB,kB,5); xmin=MIN(x,xmin); xmax=MAX(x,xmax);
			  x = getVertexXcoord(iB,jB,kB,1); xmin=MIN(x,xmin); xmax=MAX(x,xmax);
			  x = getVertexXcoord(iB,jB,kB,2); xmin=MIN(x,xmin); xmax=MAX(x,xmax);
			  x = getVertexXcoord(iB,jB,kB,6); xmin=MIN(x,xmin); xmax=MAX(x,xmax);
			}
		  else
			{
			  z = zCorn[iB][jB][kB][3]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
			  z = zCorn[iB][jB][kB][7]; zmin=MIN(zmin, z); zmax=MAX(zmax,z);
			  y = getVertexYcoord(iB,jB,kB,5); ymin=MIN(y,ymin); ymax=MAX(y,ymax);
			  y = getVertexYcoord(iB,jB,kB,1); ymin=MIN(y,ymin); ymax=MAX(y,ymax);
			  y = getVertexYcoord(iB,jB,kB,3); ymin=MIN(y,ymin); ymax=MAX(y,ymax);
			  y = getVertexYcoord(iB,jB,kB,7); ymin=MIN(y,ymin); ymax=MAX(y,ymax); 
			}
		}
	} 
  if (verbosity>0)
	{
	  if (ySec)  fprintf(stderr,"Xmin/max=%g, %g\n",xmin,xmax);
	  else  fprintf(stderr,"Ymin/max=%g, %g\n",ymin,ymax);
	  fprintf(stderr,"Zmin/max=%g, %g\n",zmin,zmax);
	}
  // asprintf(&filename,"GridXsection_%c=%d.ps",(ySec)?'j':'c',indxXsection);
  if(!strcmp(printPSfilename,"-")) fpX = stdout;
  else
	{
	  filename = printPSfilename;
	  fpX = fopen(filename,"w");
	}

  // drawing 
  width=(printPSsize/WIDTH_INCHES)*WIDTH_INCHES*72;
  height=(printPSsize/WIDTH_INCHES)*HEIGHT_INCHES*72;

  zScale = (height/(zmax-zmin));
  if (ySec)
	{ 
	  xScale = ((float)width/(xmax-xmin));
	}
  else
	{
	  xScale = ((float)width/(ymax-ymin));
	}

  /* PostScript header */
  fprintf(fpX,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fpX,"%%%%Creator: eclipseXsection.c (James Gunning)\n");
  fgetpos(fpX,&pos);
  fprintf(fpX,"%%%%BoundingBox: %4d %4d %4d %4d\n",-1,1,width+1,height+1);
  fprintf(fpX,"                           \n"); // protection against fsetpos hack below
  fprintf(fpX,"%%%%EndComments\n");
  /* The +72 shifts the figure one inch from the lower left corner */
  //fprintf(fpX,"%d %d translate\n", 72, 72 );
  fprintf(fpX,".05 setlinewidth\n");
 
	
  fprintf(fpX,"0.0 setgray\n");

  for (pass=0;pass<=1;pass++)
	{
	  if (pass==1) 
		{
		  fprintf(fpX,"\n\n%%%%Grid boxes from here till Showpage\n\n");
		  fprintf(fpX,"1.0 setgray\n");
		  fprintf(fpX,"0 0 0 setrgbcolor\n");
		}
	  for (ii=0; ii< Nmax  ; ii++) 
		{
		  if (ySec) iB=(ii+1); else jB=(ii+1);
		  for (kB=1; kB<=ndh ; kB++)
			{	
			  if (skipInActives && active[iB][jB][kB]==0) continue;

			  // color
			  col[0]=col[1]=col[2]=0.0;
			  if (colourActives) {
				// RED GREEN
				if (active[iB][jB][kB]) {
				  col[1]=1;
				} else col[0]=1;
			  } else if (showSegments) {
				// BLUE for -999's, yellows for others
				segment = seg[iB][jB][kB];
				if (segment>=0) {
				  col[0] = col[1] = 0.5+0.5*(MAX(-1,segment)-minSeg)/((float)(maxSeg-minSeg));
				} else col[2]=1;
			  } else if (haveNG) {
				// YELLOWS - BROWN 
				col[0] = col[1] = NG[iB][jB][kB];
			  } else {
				col[0]=col[1]=0.9;
			  }

			  fprintf(fpX,"newpath\n");

			  if (ySec)
				{
				  // draw
				  fprintf(fpX,"%f %lf moveto\n",xScale*(getVertexXcoord(iB,jB,kB,5)-xmin), zScale*(-zCorn[iB][jB][kB][5]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexXcoord(iB,jB,kB,1)-xmin), zScale*(-zCorn[iB][jB][kB][1]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexXcoord(iB,jB,kB,2)-xmin), zScale*(-zCorn[iB][jB][kB][2]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexXcoord(iB,jB,kB,6)-xmin), zScale*(-zCorn[iB][jB][kB][6]+zmax));
				  fprintf(fpX,"closepath\n");
				  
				  if (pass==0) 
					{
					  fprintf(fpX,"%f %f %f setrgbcolor\nfill\n",col[0],col[1],col[2]);
					  // fprintf(fpX,"%f setgray\nfill\n",(haveNG)? NG[iB][jB][kB] : GRAY_1);
					}
				  else fprintf(fpX,"stroke\n");
				}
			  else
				{
				  fprintf(fpX,"%f %lf moveto\n",xScale*(getVertexYcoord(iB,jB,kB,5)-ymin), zScale*(-zCorn[iB][jB][kB][5]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexYcoord(iB,jB,kB,1)-ymin), zScale*(-zCorn[iB][jB][kB][1]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexYcoord(iB,jB,kB,3)-ymin), zScale*(-zCorn[iB][jB][kB][3]+zmax));
				  fprintf(fpX,"%f %lf lineto\n",xScale*(getVertexYcoord(iB,jB,kB,7)-ymin), zScale*(-zCorn[iB][jB][kB][7]+zmax));
				  fprintf(fpX,"closepath\n");
				  
				  if (pass==0) 
					{
					  fprintf(fpX,"%f %f %f setrgbcolor\nfill\n",col[0],col[1],col[2]);
					  // fprintf(fpX,"%f setgray\nfill\n",(haveNG)? NG[iB][jB][kB] : GRAY_1);
					}
				  else fprintf(fpX,"stroke\n");
				}
			}
		}
	}
  fprintf(fpX,"showpage\n%%%%EOF\n");
  
  /*
	fsetpos(fp,&pos);
	fprintf(fp,"%%%%BoundingBox: %4d %4d %4d %4d\n",(int)floor(bxmin-1),(int)floor(bymin-1),
	(int)ceil(bxmax+1),(int)ceil(bymax+1));
  */
  fclose(fpX);
  printf("Wrote %s\n",filename);
}

void printLayerTop(int doX, int iXsection,int zIndx, int nrh, int nch,double ****zCorn)
{
  int i,j,iB,jB,kB;
  FILE *fpX;

  fpX = stdout;

  //  printf("Writing const ecl_i Xsection at ecl_i=%d to %s\n",iXsection,filename);
  kB = zIndx;

  // Draw top of cells:
  if (doX)
	{
	  iB = iXsection+1;
	  for (j=0,jB=1; j<nch; j++,jB++) 
		{
		  fprintf(fpX,"%lf %f\n",getVertexYcoord(iB,jB,kB,1),zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%lf %f\n",getVertexYcoord(iB,jB,kB,3),zCorn[iB][jB][kB][3]);
		}
	}
  else
	{
	  jB = iXsection+1;
	  for (i=0,iB=1; i<nrh; i++,iB++) 
		{
		  fprintf(fpX,"%lf %f\n",getVertexXcoord(iB,jB,kB,1),zCorn[iB][jB][kB][1]);
		  fprintf(fpX,"%lf %f\n",getVertexXcoord(iB,jB,kB,2),zCorn[iB][jB][kB][2]);
		}
	}
}

void drawSuppliedPoint(FILE *fpX,double x,double y)
{
  if (xyzSupplied)
	{
	  // finally draw point as cross
#if USE_FILL_FUNC_OF_GRAPH
	  fprintf(fpX,"#m=3,S=5,q=0\n");
#else
	  fprintf(fpX,"#m=3,S=5\n");
#endif
	  fprintf(fpX,"%lf %lf\n",x, y);
	}
}

/*
  Checks segment range using unit-offset segment property (i.e.
  iS = [1...Nx]. Feed -1 as index to set plane
*/
void getSegmentRange(ECL* ecl, int iS,int jS,int kS)
{
  int i,j,k;

  PROP* aP;

  if (!haveSegments) return;
  
  aP = fetchProp(ecl,"SEGMENTS");
  if (!aP || !aP->set || !aP->p) {
	fprintf(stderr,"getSegmentRange(): SEGMENTS not available\n");
	exit(-1);
  }
  seg = (int ***) aP->p;
  minSeg=INT_MAX;
  maxSeg=-INT_MAX;

  GTree* seg_hash = g_tree_new(int_compare);
  gint* value_ptr;
  if (iS>0 && jS<0 && kS<0)
	{
	  for (j=1; j<=Ny; j++)
		{
		  for (k=1; k<=Nz; k++) 
			{
			  minSeg = MIN(minSeg, seg[iS][j][k]);
			  maxSeg = MAX(maxSeg, seg[iS][j][k]);
			  value_ptr = (gint*)&seg[iS][j][k];
			  if (!g_tree_lookup(seg_hash,value_ptr)) g_tree_insert(seg_hash,value_ptr,value_ptr);
			}
		}
	}
  else if (iS<0 && jS>0 && kS<0)
	{
	  for (i=1; i<=Nx; i++) 
		{
		  for (k=1; k<=Nz; k++)
			{
			  minSeg = MIN(minSeg, seg[i][jS][k]);
			  maxSeg = MAX(maxSeg, seg[i][jS][k]);
			  value_ptr = (gint*)&seg[i][jS][k];
			  if (!g_tree_lookup(seg_hash,value_ptr)) g_tree_insert(seg_hash,value_ptr,value_ptr);
			}
		}
	}
  else if (iS<0 && jS<0 && kS>0)
	{
	  for (i=1; i<=Nx; i++) 
		{
		  for (j=1; j<=Ny; j++)
			{
			  minSeg = MIN(minSeg, seg[i][j][kS]);
			  maxSeg = MAX(maxSeg, seg[i][j][kS]);
			  value_ptr = (gint*)&seg[i][j][kS];
			  if (!g_tree_lookup(seg_hash,value_ptr)) g_tree_insert(seg_hash,value_ptr,value_ptr);
			}
		}
	}
  else
	{
	  fprintf(stderr,"Segments getSegmentRange() args invalid!\n");
	}

  fprintf(stderr,"Segments range: [%d %d], coloured as [%d,%d]\n",minSeg,maxSeg,MAX(minSeg,-1),maxSeg);

  fprintf(stderr,"Used segments: ");
  g_tree_foreach(seg_hash,print_seg,NULL);
  fprintf(stderr,"\n");
  g_tree_destroy(seg_hash);

  minSeg = MAX(0,minSeg);
  minSeg--;
  maxSeg++;
}

gint int_compare(gconstpointer a, gconstpointer b)
{
  return (*(gint*)a < *(gint*)b) ? -1 : ((*(gint*)a > *(gint*)b) ? 1 : 0);
}

gboolean print_seg(gpointer key,gpointer value,gpointer data)
{
  fprintf(stderr,"%d ",*((gint *)key));
  return FALSE;
}

double getVertexXcoord(int i, int j, int k, int vrtx)
{
  int i_c,j_c;
  double pt[4];

  if (lattice) return (double)(i+((vrtx==2 || vrtx==4 || vrtx==6 || vrtx==8) ? 1: 0));

  //  c = ecl.coord; 
  i_c = (i-1) + ((vrtx&1) ? 0 : 1);
  j_c = (j-1) + ((vrtx==3 || vrtx==4 || vrtx==7 || vrtx==8) ? 1: 0);
  interpVert(zCorn[i][j][k][vrtx],&pt[0],(ecl.coord)[i_c][j_c]);
  //*x = pt[1];
  //*y = pt[1];
  return pt[1];
}

double getVertexYcoord(int i, int j, int k, int vrtx)
{
  int i_c,j_c;
  double pt[4];

  if (lattice) return (double)(j + ((vrtx==3 || vrtx==4 || vrtx==7 || vrtx==8) ? 1: 0));

  //  c = ecl.coord; 
  i_c = (i-1) + ((vrtx&1) ? 0 : 1);
  j_c = (j-1) + ((vrtx==3 || vrtx==4 || vrtx==7 || vrtx==8) ? 1: 0);
  interpVert(zCorn[i][j][k][vrtx],&pt[0],(ecl.coord)[i_c][j_c]);
  //*x = pt[1];
  //*y = pt[1];
  return pt[2];
}

void getNearestBlockIndices(float X,float Y,float Z, int *i,int *j,int *k)
{
  int iB,jB,kB,vrtx,i_c,j_c;
  double ***c,*zc;
  double ds,dmin,pt[4];

  dmin=DBL_MAX;

  for (iB=1; iB<=Nx; iB++)
	{
	  for (jB=1; jB<=Ny; jB++) 
		{
		  for (kB=1; kB<=Nz ; kB++)
			{
			  zc = zCorn[iB][jB][kB];
			  // examine vertex 1 (bottom corner) (smallest z, z goes down...)
			  // for (vrtx=1; vrtx<=4; vrtx++)
				{
				  vrtx = 1;

				  i_c = (iB-1) + ((vrtx&1) ? 0 : 1);
				  j_c = (jB-1) + ((vrtx==3 || vrtx==4 || vrtx==7 || vrtx==8) ? 1: 0);
				  interpVert(zc[vrtx],&pt[0],(ecl.coord)[i_c][j_c]);

				  ds = SQR(X - pt[1]) + SQR(Y - pt[2])+ SQR(Z -pt[3]); 
				  if (ds < dmin) 
					{
					  *i = iB;
					  *j = jB;
					  *k = kB;
					  dmin = ds;
					}
				}
			}
		}
	}
  fprintf(stderr,"Nearest block index %d %d %d\n",*i,*j,*k);
}
