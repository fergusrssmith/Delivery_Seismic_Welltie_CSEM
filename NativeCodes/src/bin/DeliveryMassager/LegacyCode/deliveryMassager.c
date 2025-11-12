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

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>


#include "commandline/commandline.h"
#include "std/std.h"
#include "std/error.h"

#include "NR/nrutil.h"
#include "kdtree/optkd.h"
#include "bayes/bayes.h"

#include "massager_types.h"
#include "deliveryMassager.h"
#include "inputgrid.h"
#include "xmlParser.h"
#include "solver.h"
#include "util.h"

#include "ECLIO/eclio.h"  /* structure definition and a little miscellany */
#include "ECLIO/access.h" /* access high-level routines */
#include "ECLIO/blockops.h" /* lower level routines */

#if 0
  // sparse library
  #ifdef USE_LASPACK
  #include "laspack/qmatrix.h"
  #include "laspack/vector.h"
  #include "laspack/errhandl.h"   
  #endif
#endif

#define ARGS_USAGE "inputfilename.xml"

// verbosity
int verbosity = 0;
// makes output grids by duplicating input ones 
int duplicate = 0;
// makes output grids = duplicate of input plus a jiggle factor, for grid 2,3....
int duplicate_random=0;
// force sparse usage
int force_sparse=0;
// Writes summary of output grid to ascii file
int ascii_out = 0;
// Forces ignoring of soft information in update
int ignore_soft = 0;
// To make detailed moving average trending
int useMAtrend=0;
// Number (if any) realisations
int Nrealisations=0;
/*
    The option array defining the structure of the command line options.
*/
CLOption options[] =
  {
    {'v',NULL,1,  &verbosity, "verbosity",{{CL_INT_ARG,&verbosity,NULL,NULL}}},
	{'d',NULL,0,  &duplicate, "duplicate input grid to make output grid",{{CL_INT_ARG,&duplicate,NULL,NULL}}},
	{'r',NULL,0,  &duplicate_random, "randomise duplicate grids",{{CL_INT_ARG,&duplicate_random,NULL,NULL}}},
	{'t',NULL,0,  &useMAtrend, "use moving average (variogram-driven) trend",{{CL_INT_ARG,&useMAtrend,NULL,NULL}}},
	{'s',NULL,0,  &force_sparse, "force use of sparse method even if direct possible",{{CL_INT_ARG,&force_sparse,NULL,NULL}}},
	{'a',NULL,0,  &ascii_out, "write summary of massaged output grid to ascii file",{{CL_INT_ARG,&ascii_out,NULL,NULL}}},
	{'\0',"IS",0, &ignore_soft, "ignore dense soft (delivery-inversion) information",{{CL_INT_ARG,&ascii_out,NULL,NULL}}},
    {'N',NULL,1,  &Nrealisations, "number of realisations",{{CL_INT_ARG,&Nrealisations,NULL,"number of realisations"}}},
    CL_END_OPTIONS  /* defined as {0,NULL,..} to mark the end of the list   */
  };


int main(int argc, char *argv[])
{ 
  // main-specified variables
  int i,j,len;
  CLContext   context;
  char* ptr;
  char* layer_number_list;
  char* property_name_list;
  char* XYcoordinates;
  char* variogram_list;
  char* variogram_description;

  // general structure we need to populate to get the parameters
  struct MASSAGER massager;


  

  /*
	Initialise the command line context variable.
	The third 'string' parameter is printed with the usage message
	to indicate what arguments, excluding the options and their arguments,
	are required on the command line.
  */
  
  clInit( &context, options,ARGS_USAGE, NULL );
        
  /*  Parse the command line and if successful, print out its contents.
	  If not successful the function prints out the usage message.
  */
  
  if (clParse( &context, argc, argv ))
    {
	  if (context.argc != 2) 
        {
		  fprintf(stderr, "\n");
		  fprintf(stderr, "Missing some command line arguments\n");
		  clPrintUsage(&context);		  
		  return EXIT_FAILURE;
        }

	  // options parsing
	  massager.verbosity = verbosity; 

	  massager.duplicate = duplicate;
	  massager.duplicate_random = duplicate_random;

	  massager.useDenseMethod = 1;   // this flags dense==possible
	  if (force_sparse) massager.useDenseMethod = 0;  // never use dense method

	  massager.writeAsciiOut = 0;
	  if (ascii_out) massager.writeAsciiOut = 1;  // never use dense method

	  massager.ignore_soft = 0;
	  if (ignore_soft) massager.ignore_soft = 1;  // ignore soft info....

	  massager.useMAtrend = 0;
	  if (useMAtrend) massager.useMAtrend = 1;  // use moving average trend

	  massager.Nrealisations=0;
	  if (Nrealisations) 
		{
		  massager.Nrealisations = Nrealisations;
		}

	  // Parse mandatory stuff 
	  // (1) master XML file
	  massager.massage_xml_file = context.argv[1];


	  /*  Clean up memory allocated by during the parse */
	  // clDone( &context );
    }
  else
    {
	  /*  Unsuccessful parse of the command line. */
	  return EXIT_FAILURE;
    }
  // parse XML file
  // read in pre-analysed data
  if (parseXMLfile(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
  
  // Set up outputs
  if (setUpOutputs(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  // perform well observation setups
  if (setUpWellObservations(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  // TODO output summary of XML file

  // read in pre-analysed data
  if (readMassageAnalysisData(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  // read in output grid
  if (readOutputGrids(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  // set up Massaging input grid
  if (setUpInputGrid(&massager.inputgrid) == EXIT_FAILURE) return EXIT_FAILURE;

  // set up Massaging interpolation
  // if (checkOutputGridisInterior(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
  if (setUpGridMapping(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  if (setUpObservationsForDenseSystem(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
  if (setUpFullyDensePriorModel(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  if (computePosteriorModel(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

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

// Set up mapping of input grid to output grid using kriging weights derived from
// variogram
int setUpGridMapping(struct MASSAGER *massager)
{
  int i,j;
  struct PT2D pt;

  if (massager->verbosity>1)
	printf("Setting up grid mapping (kriging)\n");

  // Decide on a method
  if ((massager->totalNOutputGridPoints <= MAX_N_DIRECT_METHOD)  // sensible # points
	  && (massager->useDenseMethod == 1)) // not forced off
	{
	  if (massager->verbosity>1)
		{
		  printf("Using dense methods...\n");
		}
	}
  else
	{
	  massager->useDenseMethod = 0;
	  if (massager->verbosity>1){
		printf("Using sparse methods...\n");
	  }
	}
  

  setUpMapping(massager);

  return EXIT_SUCCESS;
}

int setUpOutputs(struct MASSAGER *massager)
{
  int i;
  char filename[100];
  /* set up possible output files */
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


// Set up mapping of input grid to output grid using kriging weights derived from
// variogram
int setUpFullyDensePriorModel(struct MASSAGER *massager)
{
  if (massager->useDenseMethod)
	{
	  if (massager->verbosity>1)
		{
		  printf("Setting up fully dense prior model\n");
		}
	  return setUpDensePriorModel(massager);
	}
  return EXIT_SUCCESS;
}

// Compute posterior model using Bayes formulas

int computePosteriorModel(struct MASSAGER *massager)
{
  if (massager->useDenseMethod)
	{
	  if (massager->verbosity>1)
		{
		  printf("Computing dense posterior model\n");
		}
	  return computePosteriorModelFromFullDenseSystems(massager);
	}
  else
	{
	  if (massager->verbosity>1)
		{
		  printf("Solving model sequentially\n");
		}
	  // return solveSequentialSparseProblems(massager);
	  return solveScaledSmoothSparseProblems(massager);
	}
}



/*
  The observations require only the "y" vector and the
  "C_D" matrix, which is diagonal. We can set these up in
  canonical ordering independent of the chosen solution method
*/
int setUpObservationsForDenseSystem(struct MASSAGER *massager)
{
  int i,ii,j,jj,k,kk;
  int Ny;
  int i_y;
  int properties_indx_k;
  int properties_indx_kk;
  int j_y;
  int input_grid_indx_offset;
  double **C_D;
  double *y;
  double val;



  if (massager->useDenseMethod)
	{
	  if (massager->verbosity>1)
		{
		  printf("Setting up dense observations y and C_D\n");
		}

	  i_y=0;
	  Ny =  massager->solverInfo.nRowsA;
	  massager->solverInfo.y = dvector(0, Ny - 1);
	  y = massager->solverInfo.y;   // shortcut
	  
	  // Allocate
	  
	  massager->solverInfo.C_D = dmatrix(0, Ny - 1,0, Ny - 1);
	  C_D = massager->solverInfo.C_D;  // shortcut
	  // Zero C_D
	  memset(&C_D[0][0], 0, Ny * Ny * sizeof(double));


	  // loop over layers
	  for (i=0; i<  massager->nLayers ; i++)
		{
		  // (1) run over well observations first
		  // run over well observations......
		  for (j=0; j<massager->well_observations[i].n; j++)
			{
			  y[i_y] = massager->well_observations[i].values[j];
			  // scale hard observations by the value
			  if (fabs(massager->well_observations[i].values[j]) > 0)
				{
				  val = HARD_OBSERVATIONS_EPS * fabs(massager->well_observations[i].values[j]);
				}
			  else
				{
				  val = HARD_OBSERVATIONS_EPS;
				}
			  // diagonal entries on C
			  C_D[i_y][i_y] = val*val;
			  i_y++;
			}

		  // (2) Then the "soft" seismic data
		  // Loop over seismic grid (same for each layer)
		  input_grid_indx_offset = i_y;
		  for (j=0; j < massager->inputgrid.nPoints; j++)
			{
			  for (k=0; k< massager->nPropertyNames; k++)
				{
				  // need this to extract the right properties index from the unrolled array massager->properties
				  properties_indx_k = k + i*massager->nPropertyNames;
				  // observation
				  y[i_y] = massager->properties[properties_indx_k].values[j];

				  // covariance error of this observation. This is the covariance massager->property_covariance
				  // read in from "delivery", for the seismic trace j
				  j_y=0;
				  for (ii=0; ii<massager->nLayers ; ii++) 
					{
					  j_y += massager->well_observations[ii].n;
					  // double loop here is wasteful.... FIXME
					  for (jj=0; jj < massager->inputgrid.nPoints; jj++)
						{
						  for (kk=0; kk< massager->nPropertyNames; kk++)
							{
							  // populate C_D here
							  if (jj==j)
								{
								  properties_indx_kk = kk+ ii*massager->nPropertyNames;
								  val = massager->property_covariance[j][properties_indx_k][properties_indx_kk];
								  C_D[i_y][j_y] =  val;
								  // and symmetrise
								  C_D[j_y][i_y] = C_D[i_y][j_y];
								  if (i_y==j_y && C_D[j_y][i_y]==0) C_D[j_y][i_y] = VARIANCE_EPS;  // FIXME.. SENSIBLE????
								}
							  j_y++;
							}
						}
					}
				  i_y++;
				}
			}
		}
	}


  return EXIT_SUCCESS;
}


// Placeholder for reading output grid(s) using Chris' library
int readOutputGrids(struct MASSAGER *massager)
{
  int i,j,k,cnt;
  struct ECL ecl;
  int status,transformToMap,nKeys;
  int skipVoids,nrl,nrh,ncl,nch,ndh,active_size;
  double ***traces,x,y;
  int **active;

  // read and parse here

  if (massager->duplicate)
	{
	  double jiggle_scale=0.0;
	  for (i=0; i< massager->nLayers; i++)
		{
		  if (!massager->duplicate_random) 
			{
			  // Example of a simple case of 2 rectangular grids, mapping input to output FIXME
			  massager->outputgrid[i].isRectangular=1;
			  massager->outputgrid[i].nPoints = massager->inputgrid.nPoints;
			  massager->outputgrid[i].xValues = massager->inputgrid.xValues;
			  massager->outputgrid[i].yValues = massager->inputgrid.yValues;
			  massager->totalNOutputGridPoints=massager->outputgrid[i].nPoints; // sum over all output grids
			}
		  else
			{
			  jiggle_scale = 0.3;
			  printf("randomising output grid %d\n",i);

			  massager->outputgrid[i].isRectangular=0;
			  massager->outputgrid[i].nPoints = massager->inputgrid.nPoints;
			  massager->outputgrid[i].xValues =  dvector(0,massager->inputgrid.nPoints-1);
			  massager->outputgrid[i].yValues =  dvector(0,massager->inputgrid.nPoints-1);
			  for (j=0; j<massager->inputgrid.nPoints; j++)
				{
				  massager->outputgrid[i].xValues[j]= massager->inputgrid.xValues[j]+jiggle_scale*(0.5+(rand()/(double)RAND_MAX-0.5));
				  massager->outputgrid[i].yValues[j]= massager->inputgrid.yValues[j]-jiggle_scale*(0.5+(rand()/(double)RAND_MAX-0.5));
				}
			  massager->totalNOutputGridPoints += massager->outputgrid[i].nPoints; // sum over all output grids
			}
		}
	}
  else
	{
	  /* Map input grid to UTMs if useUTMconversion possible */
	  if (massager->utmInfo.useUTMconversion)
		{
		  struct PT2D pt;

		  if (setUpUTMmapping(massager) == EXIT_FAILURE) return EXIT_FAILURE;

		  /* loop over seismic inline/xline grid and convert to UTM's */
		  for (j=0; j < massager->inputgrid.nPoints; j++)
			{
			  pt = getUTMSfromInlineXline(massager->utmInfo, 
										 massager->inputgrid.xValues[j],
										 massager->inputgrid.yValues[j]);
			  if (massager->verbosity>3)
				{
				  printf("IG %g %g %15.10g %15.10g\n",massager->inputgrid.xValues[j],massager->inputgrid.yValues[j],pt.x,pt.y);
				}
			  massager->inputgrid.xValues[j] = pt.x; 
			  massager->inputgrid.yValues[j] = pt.y;
			}
		} 

 

	  for (i=0; i< massager->nLayers; i++)
		{
		  /* read Eclipse grids */
		  nKeys = 7; 
		  transformToMap = 1; /* number of records to read and
								 force transform grid->map coords */
		  if (massager->verbosity>1) printf("...loading model from file %s\n", massager->outputgrid[i].filename);
		  status = loadModel(massager->outputgrid[i].filename, &ecl, transformToMap, nKeys, 
							 "PINCH", "MAPAXES", "GRIDUNIT", "SPECGRID", 
							 "COORD", "ZCORN", "ACTNUM"); 
		  if (status) 
			{
			  printf("ERROR: Could not load model, exiting\n");
			  return EXIT_FAILURE;
			}
		  skipVoids = 1; /* don't include void blocks in gross count */
		  nrl = 0; 
		  nrh = getNX(&ecl); /* subranges n*l > 0 and n*h < nx are allowed */
		  ncl = 0; 
		  nch = getNY(&ecl); 
		  ndh = getNZ(&ecl);
		  traces = fetchTraces(&ecl, skipVoids, nrl, nrh, ncl, nch);
		  if (!traces)
			{ 
			  printf("ERROR: Returned null trace array\n");
			  return EXIT_FAILURE;
			}
		  printf("...fetching active trace indicators for model %s\n",ecl.name); 
		  active = fetchActCoord(&ecl, nrl, nrh, ncl, nch);
		  if (!active)
			{
			  printf("ERROR: Returned null active array\n");
			  return EXIT_FAILURE;
			}
		  printf("...fetched active array for model %s\n",ecl.name); 
		  
		  /* get grid size */
		  active_size = countActCoord(&ecl);
		  printf("Active grid size = %d\n",active_size);

		  /* alloc and load */
		  massager->outputgrid[i].isRectangular=0;
		  massager->outputgrid[i].nPoints = active_size;
		  massager->outputgrid[i].xValues =  dvector(0,active_size-1);
		  massager->outputgrid[i].yValues =  dvector(0,active_size-1);
		  massager->totalNOutputGridPoints += massager->outputgrid[i].nPoints; // sum over all output grids
 
		  cnt=0;
		  for (j=nrl;j<=nrh;j++)
			{
			  for (k=ncl;k<=nch;k++) 
				{
				  if (active[j][k]) 
					{
					  /* midpoints of traces: 6 values are x,y,z_top,x,y,z_bottom */
					  massager->outputgrid[i].xValues[cnt] = 0.5*(traces[j][k][1]+traces[j][k][4]);
					  massager->outputgrid[i].yValues[cnt] = 0.5*(traces[j][k][2]+traces[j][k][5]);
					  if (massager->verbosity>3)
						{
						  printf("GRID %15.10g %15.10g\n",massager->outputgrid[i].xValues[cnt],massager->outputgrid[i].yValues[cnt]);
						}
					  cnt++;
					}
				}
			}

		  /* free the ecl data  */
		  freeModel(&ecl);
 

		  if (massager->verbosity>1) printf("...finished loading model from file %s\n", massager->outputgrid[i].filename); 
		}
	}
  if (massager->verbosity>1)
	{
	  int i;
	  for (i=0; i< massager->nLayers; i++)
		{
		  printf("Total number of output grid points on layer %d = %ld\n",i,massager->totalNOutputGridPoints);
		}
	}

  return EXIT_SUCCESS;
} 


// read massage-analysed data put out by InversionAnalyser --massage-analyse
int readMassageAnalysisData(struct MASSAGER *massager)
{
    char           *filename;
	double** data_ptr;
    FILE           *in;
	char *ptr;
	char	line_buffer[LINE_LENGTH];
	int i,j,k,line_no,nMatches=0,number_of_columns,layerNo,found;
	int nNonCovarianceColumns;
	int n_cov_entries=0,n_prop_entries=0,offset;
	double f;
	int** cov_column;
	int* cov_col_offset;
	int* cov_row_offset;
	fpos_t pos;

	if (massager->verbosity>1)
	  printf("Reading Massage-analysis data\n");

    filename = massager->massage_analysis_file;
    if ((in = fopen(filename, "r")) == NULL) 
	  {
		warning("cannot open file %s",filename);

        return EXIT_FAILURE;
	  }
	// Parse:
	// title
	fgets(line_buffer,LINE_LENGTH,in);
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
	// here's a brain-dead but simple way to remember the column numbers of the properties
	// in the ile *in the natural ordering of the Summary-file*
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

	// now count number of useful lines
	massager->inputgrid.nPoints=0;
	fgetpos(in,&pos);
	{
	  int line_ok;
	  for(;;)
		{
		  line_ok=1;
		  // check we can read number_of_columns on this line
		  for (i=0;i<number_of_columns; i++)
			{
			  if (fscanf(in," %lf ",&f)!=1) line_ok=0;
			}
		  if (line_ok) massager->inputgrid.nPoints++;
		  else break;
		}
	}
	// printf("Found %d lines\n",massager->inputgrid.nPoints);
	// rewind
	if (fsetpos(in,&pos)) warning("Can't reposition file pointer");

	// allocate memory. A property per layer and name, plus two for x,y
	for (j=0; j<(massager->nProperties + 2); j++)
	  {
		massager->properties[j].values = (double*)malloc((size_t) ((massager->inputgrid.nPoints) * sizeof(double)));
	  }
	// observations covariance memory
	massager->property_covariance = d3tensor(0,(massager->inputgrid.nPoints)-1,0,(massager->nProperties)-1,0,(massager->nProperties)-1);

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
		cov_row_offset[i] = cov_col_offset[i] = -1;  // flag to assume not-needed until shown otherwise
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
	// read
	nNonCovarianceColumns = number_of_columns-n_cov_entries;

	for (line_no=0; line_no<massager->inputgrid.nPoints; line_no++)
	  {
		// properties
		for (i=0; i<nNonCovarianceColumns; i++)
		  {
			if ((nMatches = fscanf(in," %lf ",&f)) == 1)
			  {
				if (data_ptr[i]!=NULL) data_ptr[i][line_no] = f;
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
				if ((cov_row_offset[i]>=0) && (cov_col_offset[i]>=0))
				  {
					massager->property_covariance[line_no][cov_row_offset[i]][cov_col_offset[i]] = f;
					massager->property_covariance[line_no][cov_col_offset[i]][cov_row_offset[i]] = f;
				  }
			  }
			else
			  {
				warning("Can't read number in file %s, file-line-number %d\n",filename,line_no+2+number_of_columns);
				return EXIT_FAILURE;
			  }
		  }
		fscanf(in," \n");
	  }

	fclose(in);

	free(data_ptr);
	free_imatrix(cov_column,0,n_prop_entries-1,0,n_prop_entries-1);
	free_ivector(cov_row_offset,0,number_of_columns-1);
	free_ivector(cov_col_offset,0,number_of_columns-1);

    return EXIT_SUCCESS;
}

// set up well observations property indexes
int setUpWellObservations(struct MASSAGER *massager)
{
  int i,j,k,m;

  for (i=0; i<massager->nLayers; i++)
	{
	  massager->well_observations[i].n_obs = ivector(0,massager->nPropertyNames-1);                  // FIXME free at end!!!
	  for (k=0; k<massager->nPropertyNames; k++) massager->well_observations[i].n_obs[k]=0;

	  for (j=0; j<massager->well_observations[i].n; j++)
		{
		  for (k=0; k< massager->nPropertyNames; k++)
			{
			  if (!strcmp(massager->well_observations[i].property_name[j], massager->propertyNames[k]))
				{
				  massager->well_observations[i].propertyIndex[j] = k;
				  (massager->well_observations[i].n_obs[k])++;    // count # observations of property k in layer i
				}
			}
		}
	  // allocate convenience structures
	  massager->well_observations[i].obs = (struct OBS**)malloc((size_t) (massager->nPropertyNames * sizeof(struct OBS*)));
	  // fill them
	  for (k=0; k< massager->nPropertyNames; k++)
		{
		  if (massager->well_observations[i].n_obs[k] > 0)
			{
			  massager->well_observations[i].obs[k] = (struct OBS*)malloc((size_t) (massager->well_observations[i].n_obs[k] * sizeof(struct OBS)));
			  struct OBS* obs = massager->well_observations[i].obs[k];
			  for (j=0,m=0; j<massager->well_observations[i].n; j++)
				{
				  if (massager->well_observations[i].propertyIndex[j] == k)
					{
					  obs[m].x =  massager->well_observations[i].xValues[j];
					  obs[m].y =  massager->well_observations[i].yValues[j];
					  obs[m].val =  massager->well_observations[i].values[j];
					  m++;
					}
				}
			}
		}
	}
  return EXIT_SUCCESS;
}

int setUpUTMmapping(struct MASSAGER *massager)
{
  int i,j,k;


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
  if (fabs(det) < 1e-6)
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
  double eps = 1e-10;
  for (i=0;i<3;i++)
	{
	  pt = getUTMSfromInlineXline(massager->utmInfo,massager->utmInfo.inLine[i],massager->utmInfo.xLine[i]);
	  // printf ("UTM transformation check %g %g\n",pt.x - massager->utmInfo.x[i],pt.y - massager->utmInfo.y[i]);
	  if ((fabs(pt.x - massager->utmInfo.x[i]) > 1e-10) || fabs(pt.y - massager->utmInfo.y[i]))
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
