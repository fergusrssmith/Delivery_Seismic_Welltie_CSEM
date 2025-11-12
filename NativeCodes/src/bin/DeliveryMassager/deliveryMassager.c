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
// #include <errno.h>
#include <stdio.h>

#include <libxml/parser.h>


#include "ECLIO/eclio.h"  /* structure definition and a little miscellany */

#include "commandline/commandline.h" 
#include "errorHandler.h" 

#include "NR/nrutil.h"

#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"
#include "cgeom/vdefs.h"

#include "util.h" 
#include "massager_types.h"
#include "clapack_cblas_wrappers.h"
#include "massager.h"
#include "deliveryMassager.h" 
#include "xmlParser.h"
#include "solver.h" 
#include "coarseGridMaker.h"
#include "eclipseHandler.h"

#define ARGS_USAGE "massage-filename.xml"

// verbosity
int verbosity = 0;
// makes output grids by duplicating input ones 
int duplicate = 0;
// force sparse usage
int force_sparse=0;
// Writes summary of output grid to ascii file
int ascii_out = 0;
// Writes summary of output grid to su file
int su_out = 0;
// Writes summary of output grid to ECLIPSE file
int ecl_out = 0;
// Create eclipse grid from input seismic grid
int make_ecl_grid = 0;
char *make_ecl_options = NULL;
// Forces ignoring of soft information in update
int ignore_soft = 0; 
// Forces ignoring of well observation data in update
int ignore_well_data = 0;
// Number (if any) realisations
int Nrealisations=0;
// Random seed if realisations drawn
int Seed=-1;

// Acquire stats flag
int acquireStats=0;
int useStatsCutoff=0;
float statsCutoff=0;
// Overrride stats sorting property default name (net-sand)
int overrideSortPropName=0;
char* statsSortPropName = NULL;

// Decorations flags
int nDecorations=0;
int dumpPSdecoration=0; 
int dumpPSdecorationVal=0;
char dumpPSdecorationDir='\0';

// Depth uncertainty
int includeSystematicDepthUncertainty = 0;
float systematicDepthUncertainty = 0;

// make subgrid only
int makeSubgrid=0;

// make thin blocks active
int makeThinBlockActive = 0;
// patch segment labels 
int patchSegmentLabels = 0;
// alternative MAF file
int useAltMassageAnalyseFile = 0;
char* alt_massage_analysis_file = NULL;

// metricize
int flipGridunits=0;

// Pollution filter
int usePropertyMask=0;
char* maskPropName=NULL;
float propertyMaskMin=0;
float propertyMaskMax=0;

/*
    The option array defining the structure of the command line options.
*/
CLOption options[] =
  {
    {'v',NULL,1,  &verbosity, "verbosity",{{CL_INT_ARG,&verbosity,"level",NULL}}},
    {'\0',"seed",1,  &Seed, "random-seed",{{CL_INT_ARG,&Seed,NULL,NULL}}},
	{'d',NULL,0,  &duplicate, "duplicate input grid to make output grid. No ECLIPSE output possible."},
	{'a',NULL,0,  &ascii_out, "write summary of massaged output grid to ascii file"},
	{'\0',"su",0,  &su_out, "write summary of massaged output grid to (big-endian) SU file, Delivery-model-style. Requires -d."},
	{'\0',"ecl",0,  &ecl_out, "write summary of massaged output grid to eclipse file"},
	{'\0',"MAF",1,&useAltMassageAnalyseFile, 
	 "use supplied delivery_massage_analysis_filename, rather than named one in XML file",
	 {{CL_STRING_ARG,&alt_massage_analysis_file,"massage-analyse file",NULL}}},
	{'\0',"make-eclipse",1,  &make_ecl_grid, "make eclipse file from seismic grid geometry in massage-analyse file",
	 {{CL_STRING_ARG,&make_ecl_options,"options","comma separated selections from {d,MBR,Nu=?,Nv=?}, "
	   "e.g. MBR,Nu=2,Nv=4\n"
	   "\t\td=defaults,MBR=minumum bounding rectangle (might not be the rectangle aligned with seismic survey)\n"
	   "\t\tNu=,Nv= numbers of blocks in principal directions, overriding inferences from Delivery inversion spacing"}}},
	{'\0',"IS",0, &ignore_soft, "ignore dense soft (delivery-inversion) information (not yet functional)"},
	{'\0',"IW",0, &ignore_well_data, "ignore any hard well data in XML file"},
    {'N',NULL,1,  &Nrealisations, "number of realisations",{{CL_INT_ARG,&Nrealisations,"number",NULL}}},
	{'\0',"patch-ecl-segments",0,  &patchSegmentLabels, "homogenise segment labels vertically in eclipse input grid (use majority class)"},
	{'\0',"make-thin-ecl-blocks-active",0,  &makeThinBlockActive, " detect and make \"active\" any thin blocks in eclipse input grid"},
    {'\0',"systematic-depth-uncert",1, &includeSystematicDepthUncertainty, "specify reference depth standard deviation",
	 {{CL_FLOAT_ARG, &systematicDepthUncertainty,  "standard-deviation",NULL}}},
	{'\0',"stats",0,  &acquireStats, 
	 "show rough stats of volumetrics on net-sand, thickness if drawing realisations, sorted on net-sand",{{CL_INT_ARG,&ascii_out,NULL,NULL}}},
	{'\0',"stats-sorting-property",1,&overrideSortPropName, "non default (net-sand) property name to use in statistics CDFs",{{CL_STRING_ARG,&statsSortPropName,"property-name",NULL}}},
	{'\0',"stats-lower-cutoff",1,&useStatsCutoff, "lower cutoff for stats-sorting-property",{{CL_FLOAT_ARG,&statsCutoff,"value",NULL}}},
 	{'\0',"make-subgrid",0,  &makeSubgrid, "vertically decimate Eclipse grid in XML according to best match of depth+thickness\n\t\twell observation data in XML, to make coarse grid. No massaging performed."},
   {'\0',"decorate",1,&nDecorations, "decorate simulations (N per simulation) to match simulated NG",{{CL_INT_ARG,&nDecorations,"N",NULL}}},
    {'\0',"decorate-Xsection",2,&dumpPSdecoration, "print Xsection of decorations (x,y is direction of xsection, number is inline/xline or ecl index)",
	 {{CL_STRING_ARG, &dumpPSdecorationDir,"x|y number",NULL},
	  {CL_INT_ARG,&dumpPSdecorationVal,NULL,NULL}}},
    {'m',"property-mask",3,  &usePropertyMask, "Mask input traces to acceptable interval if required (grid reduction).",
	 {{CL_STRING_ARG, &maskPropName,"property-name",NULL},
	  {CL_FLOAT_ARG,&propertyMaskMin,"min",NULL}, 
	  {CL_FLOAT_ARG,&propertyMaskMax,"max",NULL}}},
	{'\0',"flip-gridunits",0,  &flipGridunits, "change input Eclipse gridunits to opposite units (feet <-> meters)"},
    CL_END_OPTIONS  /* defined as {0,NULL,..} to mark the end of the list   */
  };  
  



int main(int argc, char *argv[])
{ 
  // main-specified variables
  int len;
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

	  // massager initialisation (safe defaults)
	  init_massager(&massager);

	  // options parsing
	  massager.verbosity = verbosity; 

	  massager.duplicate = duplicate;

	  massager.makeThinBlockActive = makeThinBlockActive;

	  massager.flipGridunits = flipGridunits;

	  massager.patchSegmentLabels = patchSegmentLabels;

	  if (ascii_out) massager.writeAsciiOut = 1;
	  if (su_out) 
		{
		  if (!massager.duplicate)
			{
			  fprintf(stderr,"SU output requires -d (duplicated grid) option!\n"); 
			  return EXIT_FAILURE;
			}
		  massager.writeSUOut = 1;
		}

	  if (make_ecl_grid) 
		{
		  massager.make_ecl_grid = 1;
		  if (massager.duplicate) printf("--make-eclipse gets geometry from seismic, so -d is redundant (disabling -d)\n");
		  if (!ecl_out) printf("--make-eclipse effectively enables --ecl (--ecl enabled)\n");
		  massager.make_ecl_grid_options.make_ecl_options = make_ecl_options;
			
		  massager.duplicate = 0;
		  ecl_out = 1;
		}

	  // cannot write eclipse grids from duplicated files
	  if (ecl_out && !(massager.duplicate)) massager.eclInfo.writeECLOut = 1;

	  if (ignore_soft) 
		{
		  massager.ignore_soft = 1;  // ignore soft info....
		  printf("Beware!!: --IS sets all delivery covariances to the identity matrix, for all properties\n");
		}
	  
	  if (ignore_well_data) 
		{
		  massager.ignore_well_data = 1;  // ignore well data
		}

	  

	  if (useStatsCutoff)
		{
		  massager.statsCutoff = statsCutoff;
		}
	  if (Nrealisations) 
		{
		  massager.Nrealisations = Nrealisations;
		  massager.acquireStats = acquireStats;
		  asprintf(&(massager.statsSortPropName),"net-sand");
		}
	  if (overrideSortPropName)
		{
		  free(massager.statsSortPropName);
		  massager.statsSortPropName = statsSortPropName;
		}
	  if (Seed) 
		{
		  massager.seed = -(int)fabs((float)Seed);
		}

	  // decoration
	  if (nDecorations) 
		{
		  massager.nDecorations = nDecorations;
		  if (dumpPSdecoration)
			{
			  if (dumpPSdecorationDir!='x' && dumpPSdecorationDir!='y')
				{
				  fprintf(stderr, "decorate-Xsection must be in x or y direction\n");
				  return EXIT_FAILURE;
				}
			  massager.decoratorInfo.dumpPSdecoration=1;
			  massager.decoratorInfo.dumpPSdecorationVal = dumpPSdecorationVal;
			  massager.decoratorInfo.dumpPSdecorationDir = dumpPSdecorationDir;
			}
		}

	  // Depth uncert
	  if (includeSystematicDepthUncertainty)
		{
		  massager.systematicDepthUncertainty = systematicDepthUncertainty;
		} 

	  // masking
	  if (usePropertyMask)
		{
		  massager.inputgrid.filt.apply=1;
		  massager.inputgrid.filt.propName = maskPropName;
		  massager.inputgrid.filt.min = propertyMaskMin;
		  massager.inputgrid.filt.max = propertyMaskMax;
		}

	  // check make-eclipse options
	  if (parse_make_ecl_grid_options(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  if (useAltMassageAnalyseFile) 
		{
		  asprintf(&massager.massage_analysis_file,"%s",alt_massage_analysis_file);
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

  // Set up SU if requested
  if (setUpSUOutputs(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  // perform well observation setups
  if (setUpWellObservations(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  if (makeSubgrid)
	{
	  // Clunky use of massager framework as a subgrid generator only
	  if (makeOptimalCoarseLayeredECL(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  // cleanup
	  if (cleanUpSundryMemory(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}
  else
	{
	  // Usual massager workflow

	  // Set up outputs
	  if (setUpOutputs(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  // TODO output summary of XML file

	  // read in pre-analysed data
	  if (readMassageAnalysisData(&massager) == EXIT_FAILURE) return EXIT_FAILURE;


	  // add systematic depth uncertainty if required
	  if (addSystematicDepthUncertainty(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  // read in output grid
	  if (readOutputGrids(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

 	  // set up Massaging interpolation
	  // if (checkOutputGridisInterior(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
	  if (setUpGridMapping(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

	  if (solveScaledSmoothSparseProblems(&massager) == EXIT_FAILURE) return EXIT_FAILURE;

  
	  if (cleanUpMappingSetup(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
	  if (cleanUpGrids(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
	  if (cleanUpSundryMemory(&massager) == EXIT_FAILURE) return EXIT_FAILURE;
	}



  return EXIT_SUCCESS;
}
