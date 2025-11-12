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
#ifndef    MASSAGER_TYPES
#define    MASSAGER_TYPES

#define _GNU_SOURCE

#define DELIVERY_SOURCE

#include <stdio.h>
#include <glib.h>

#include "ECLIO/eclio.h"

/* MINIMUM list of useful keywords when reading ECLIO files */
#define N_DESIRED_ECLIPSE_KEYWORDS 10
static char* DESIRED_ECLIPSE_KEYWORDS[N_DESIRED_ECLIPSE_KEYWORDS] 
     = {"PINCH","MAPAXES","GRIDUNIT","SPECGRID","COORDSYS","COORD","ZCORN","ACTNUM","SEGMENTS","LAYERS"};

#define PETREL_HACKS 1 
#define MAX_N_DIRECT_METHOD 4000
#define HARD_OBSERVATIONS_EPS 0.0001      // C_D scaling entry for "well" observation
#define VARIANCE_EPS 1e-6                 // variance epsilon for zero properties
#define WELL_XY_MATCH_EPS 1e-6            // tolerance to regard well observations as co-located
#define DUPLICATE_AT_FAULTS 1
#define USE_IN_PLACE_SMOOTHING 1
#define MAXWARNINGS 10

#ifndef SQ
#define SQ(a) ((((a) == 0.0) ? 0.0 : (a)*(a)))
#endif

// 2D point
struct PT2D
{ 
  double x,y;
};

// index-char pair
struct INTCHAR
{
  int indx;
  char c;
};

// index pair
struct INTPAIR
{
  int indx;
  int val;
};

// short variable arrays
struct SEG_ENCODING
{
  short size;
  int* values;              // values[0...size-1] = list of distinct segment labels of segments meeting at corner (of 4 gridblocks)
  unsigned char* corner_encoding;
  // binary encoding of which corner 'c' values belong to segment
  // e.g. if values[0]=3 && corner_encoding[0] = 18 = 2^1+2^4, corners 1 and 4 belong to segment 3
};

// observation triple
struct OBS
{
  double x,y;        // position of observation 
  double val;        // value of observation
  double smoothVal;  // value of smooth mean surface value at this point (useful convenience)
};

// float+index pair
struct FLPAIR
{
  int indx;
  float val; 
};

// double+index pair
struct PAIR
{ 
  int indx;
  double val;
};

// PAIR_ARRAY 
struct PAIR_ARRAY
{
  int size;
  struct PAIR* values;
};

// FLPAIR_ARRAY
struct FLPAIR_ARRAY
{
  int size;
  struct FLPAIR* values;
};

struct UTM_CONVERSION_INFO
{
  int performInputGridUTMmapping;
  int performOutputGridUTMmapping;
  int utmTripleAvailable;
  double x[3];
  double y[3];
  double inLine[3];
  double xLine[3];
  double mapMatrix[2][2];
  double mapOffset[2];
};

// Property structure
struct PROPERTY
{
  char* name;
  double* values;
  int layerNumber;
  int normalised;
  int smooth_across_faults;
  int use_trend_smoothing_special_variogram;
  // convenience stuff
  int inputColNumber;            // which column in the Analyser summary file
  int inputPropNumber;           // which orderiing index in the Analyser summary file
  int found;                     // will be 1 if found in Analysis file
  int applyFilter;               // to exclude parts of the inversion if required
};

// Variogram structure
struct VARIOGRAM_2D
{
  int GSLIB_type;
  double sill,nugget,rangeX,rangeY;   // lengths presumably specified in same units as seismic grid etc.....
  double azimuth;   // azimuth here in degrees always

  /*
	utility variables for efficiency
  */
  double costh,sinth;
};

// 3D Variogram structure
struct VARIOGRAM_3D
{
  int GSLIB_type;
  double nugget,sill;
  double a_hmax,a_hmin,a_vert; 
  double ang1, ang2, ang3;

  /*
	utility variables for efficiency
  */
  double rotmat[3][3];
};

/* 
   Statistics may require constraints
*/
struct STATS_CONSTRAINTS
{
  int useStochasticDepth;
  char* sortingPropName;
  double Zmax,Zmin;
};

struct ECL_OUT_INFO
{
  ECL eclOut;
  char* outputFilename;   /* filename of output grid to dump massager ML output into */
  FILE* ecl_out_fp;       /* and associated file pointer */
  int haveLayers;
  int nKeys;
  int nDuplicateKeys;     /* the number of keywords duplicated from the ECLIO input file before properties are added */
  int nGridBlockKeys;     /* the number of gridblock keywords duplicated from the ECLIO input file before properties are added */
  int** bounds;
  float*** property;      /* storage for the post-massaged ML properties that will be written */
  float**** r_properties; /* storage for the post-massaged realisations properties that will be written */
  char **keyWords;
  int **active;           /* ptr to active gridblocks */

  /* stuff for ecl realisations output */
  char* outputFilename_r;   /* filename of output grid to dump massager ML output into */
  FILE* ecl_out_r_fp;       /* realisations file pointer */
  int nGridBlockKeys_r;     /* the number of gridblock keywords duplicated from the ECLIO input file before properties are added */
  int nKeys_r;
  char **keyWords_r;        /* realisation keywords */

  // stats constraints
  struct STATS_CONSTRAINTS statsConstr;
};

typedef struct ECL_GRID_OPTIONS
{
  char *make_ecl_options; // command string
  char use_default;
  int Nu,Nv;
  char MBR;
} ECL_GRID_OPTIONS;

struct SU_OUT_INFO
{
  int nLayersPossibleForSUOutput;
  int nBodyPropertiesPossibleForSUOutput;
  char inTime;
  char getLastSurfacefromThickness;
  char getLastSurfacefromTbase;
};

struct INPUTGRID_FILTER 
{
  int apply;
  char* propName;
  float min;
  float max;
};

struct INPUTGRID
{
  int nPoints;
  char* xWord;
  char* yWord;
  double* xValues;
  double* yValues;
  int *segnum;
  unsigned char *isUsedByInterpolation; 
  GHashTable *availableSegments;
  int utmXYHeaderWordsAvailable;
  struct INPUTGRID_FILTER filt;
  double traceArea;
  // convenience stuff
  int Nx,Ny;
  struct PT2D rFirst,rLast,dRx,dRy;
  // Stuff for use in kd trees
  double** Points;       // for use in KD tree
  optkdNode *OptkdRoot;  // in KD tree
};


struct OUTPUTGRID
{
  // specified filename
  char* filename;
  char* layer_name;
  char* ECLIPSE_segments_keyword;
  char* ECLIPSE_layers_keyword;
  int   ECLIPSE_layer_number;
  int   ECLIPSE_reference_depth_layer_number;
#if DUPLICATE_AT_FAULTS
  struct SEG_ENCODING **segInfo;
#endif
  // useful eclipse details
  ECL   ecl;                             /* Eclipse structure for input grid */
  int** ecl_active;                      /* matrix of active columns */
  int   ecl_Nx,ecl_Ny;                    /* matrix dimensions */

  struct ECL_OUT_INFO eclOutInfo;        /* Info to do with ECLIPSE file outputs */

  // generic info
  int nPoints;
  int isAMirror;
  int isRectangular;
  double* xValues; 
  double* yValues;
  double* zValues;
  float* areas;
  double xMax,xMin,yMax,yMin;
  // possible segmentation from eclipse files
  int *segnum;
  unsigned char *segnum_corner_encoding;
  // underlying Eclipse indexing
  int *ecl_i;
  int *ecl_j;
  int **segcol;
  // Stuff for use in kd trees
  double** Points;       // for use in KD tree
  optkdNode *OptkdRoot;  // in KD tree
};

struct WELL_OBSERVATIONS
{
  int n;                 // number of observations in this layer i = 0...n-1
  GArray* xValues;       // at position xValues[i], yValues[i],
  GArray* yValues;       // property_name[i] has value values[i]
  GArray* propertyIndex; // index of property in array of properties: propertyIndex[i] = property index k of ith observation
  GArray* values;
  GArray* property_name;
  // convenience stuff which we set up for kriging
  int* n_obs;            // array of the number of observations of property k: sum_k n_obs[k] = n
  struct OBS **obs;      // obs[k][j].x   is the (x,y position,value) of observation j of the kth property in this layer. j=0...n_obs[k]-1
};

// Internal truncated Gaussian simulation info
struct COLUMN_TGS_KRIG_WEIGHTS
{
  int N;  // N nbrs 
  double deviate_scale; // sqrt(kriging variance)
  // int* nbr_indices; // nbr indexes
  double** nbr_pointers; // nbr ptrs -> faster than indexes
  double* wts;      // wts[i] = weights j=0...N-1 for point i
};

// Internal truncated Gaussian simulation info
struct DECORATOR_INFO
{
  int Nz;                 // grid subdivision
  int decoration_layer;   // which massager layer to decorate

  double z_stretch;
  int Nattempts;
  double t_tree_timing; // 3D kdtree timing
 
  struct VARIOGRAM_3D variogram;
  
  // postscript outputs
  int dumpPSdecoration;
  int dumpPSdecorationVal;
  char dumpPSdecorationDir;

  // stuff for simulation
  int decoration_layer_indx;
  int NG_index,thickness_index,d_index;
  float** NGtarget;
  // float** NGerror;
  float** thickness;
  float** depth;
  long seed;             // Seed used for realisations

  struct COLUMN_TGS_KRIG_WEIGHTS **col_weights;
  int* threeDpath;
  double* simulation;
  short* facies_simulation;

  // Stuff for use in kd trees
  dynamicTree *dTree;
  double** Points;       // 3D grid for use in KD tree
  short*   flag;         // 3D flag for known hard data
  optkdNode *OptkdRoot;  // in KD tree
};

// Solver info
struct SOLVER_INFO
{
  // fancy lookup stuff
  struct FLPAIR_ARRAY** OGkrigInfoSoft;     // Arrays of kriging weights and indices for output grid, for soft info
  struct FLPAIR_ARRAY** WellkrigInfo;       // Arrays of kriging weights and indices for each grid, for well info

  float  ***smoothMeans;                  // Array for smooth mean estimate if used
  float  ***smoothPropertyCovariance;     // Smoothed version of "delivery" generated observations covariance: [i][j][k]
                                          // i = input grid location index
                                          // j,k = 0,1...(nProperties-1) property indices
  float  **smoothPropertyCovariance_tri;  // Upper triangular version of above
  int smoothPropertyCovariance_tri_size;  // Number of members of above that are still alloced after 'defragging'
};

// Eclipse handling
struct ECLIPSE_INFO
{
  int useECLIOsegments;              // use segments info in ECLIPSE files for fault blocks
  int commonEclipseGrid;             // is eclipse grid common to all layers
  int writeECLOut;                   // write massaged delivery files to eclipse file
  int writeInSitu;                   // write Ecl outputs in situ in one file
  float ****property;                // holds block-centred arrays for eclio inSitu output
};

// Surrogate for the data in a massager object
struct MASSAGER 
{
  struct INPUTGRID inputgrid;
  struct OUTPUTGRID *outputgrid;
  struct UTM_CONVERSION_INFO utmInfo;
  long totalNOutputGridPoints;
  int nPropertyNames;                // Number of property types, eg {NG, net-sand} = 2
  char** propertyNames;              // Property names
  int nProperties;                   // Number of property types * nLayers
  int nLayers;
  int *layerNumbers;
  struct PROPERTY *properties;       // array of properties we wish to krige
  float*** property_covariance;      // "delivery" generated observations covariance: [i][j][k]
                                     // i = input grid location index
                                     // j,k = 0,1...(nProperties-1) property indices
  float** property_covariance_tri;   // possible upper triangular equivalent of same
  int property_covariance_tri_size;  // Number of members of above that are still alloced after 'defragging'
  unsigned char* property_NAN;       // Flag to detect NAN's in traces
  int thickness_index;               // index of thickness if used.
  int net_sand_index;                // index of net-sand if used.
  int NG_index;                      // index of NG if used.
  int d_index;                       // index of d if used
  int t_index;                       // index of t if used.
  int t_base_index;                  // index of t_base if used.
  int vp_eff_index;                  // index of vp_eff if used.
  char* massage_xml_file;            // master XML file controlling everything
  char* massage_analysis_file;
  struct VARIOGRAM_2D *variograms;                       // list over output layers
  struct VARIOGRAM_2D trend_smoothing_variogram;         // trend smoother the same over all layers
  struct VARIOGRAM_2D trend_smoothing_special_variogram; // (optional) special trend smoother the same over all layers
  struct WELL_OBSERVATIONS *well_observations;           // list of well observations
  int verbosity;
  int duplicate;                     // make fake output grids by duplicating input ones
  int Nrealisations;
  int writeAsciiOut;                 // write massaged delivery files to summary ascii files 
  int writeSUOut;                    // write massaged delivery files to summary SU files
  int make_ecl_grid;                 // create eclipse grid from analyse geometry
  ECL_GRID_OPTIONS make_ecl_grid_options; // create eclipse grid options
  struct SU_OUT_INFO *su_out_info;
  int makeThinBlockActive;           // activate thin blocks in eclipse grids
  int flipGridunits;                      // convert input ecl to metric if requested
  int patchSegmentLabels;            // patch segment labels vertically in eclipse grid
  struct ECLIPSE_INFO eclInfo;       // Common information about eclipse file handling
  int ignore_soft;                   // forces ignoring of soft information
  int ignore_well_data;              // forces ignoring of well observations
  int diagonalCovariancesOnly;       // Flag to carry memory savings on covariances
  int acquireStats;                  // compute some useful stats as we go with realisations
  float statsCutoff;                 // lower bound for acquiring stats contributions (e.g. least accepatble net-sand)
  char* statsSortPropName;           // stats sorted on this property name
  int nDecorations;                  // do 3D simulation
  float systematicDepthUncertainty;
  struct SOLVER_INFO solverInfo;
  struct DECORATOR_INFO decoratorInfo; 
  char** ascii_out_filenames;
  char** ascii_out_realisations_filenames;
  FILE** ascii_out_fp;
  FILE** ascii_out_realisations_fp;
  long seed;                         // Seed used for realisations
};


#endif
