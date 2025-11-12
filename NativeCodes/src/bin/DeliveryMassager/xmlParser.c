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
/**
 * parse the XML <massaged_inversion> file and populate relevant
 * fields of the deliveryMassager structure
 */

#include <stdio.h>
#include <string.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
 
#include <glib.h>

#include "NR/nrutil.h"

#include "kdtree/optkd.h"
#include "kdtree/dynamic-tree.h"

#include "massager_types.h"
#include "xmlParser.h"
#include "variogram.h"
#include "decorator.h"
#include "clapack_cblas_wrappers.h"
#include "massager.h"

int parseXMLfile(struct MASSAGER *ma)
{
  int i,m;
  struct XML_PARSER_GLOBALS glbs;
  xmlDocPtr doc; /* the document tree */
  xmlNodePtr cur; 

  doc = xmlParseFile(ma->massage_xml_file);
  if (doc == NULL) 
	{
	  fprintf(stderr, "Failed to parse %s\n", ma->massage_xml_file);
	  return EXIT_FAILURE;
    }
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) 
	{
	  fprintf(stderr,"empty document\n");
	  xmlFreeDoc(doc);
	  return EXIT_FAILURE;
	}

  /* Initialisations */

  glbs.layer_indx=-1;
  glbs.in_layer=0;
  glbs.in_variogram=0;
  glbs.in_trend_variogram=0;
  glbs.in_special_trend_variogram=0;
  glbs.found_special_trend_variogram=0;
  glbs.in_observation=0;
  glbs.in_stats=0;
  glbs.in_decoration=0;
  ma->massage_analysis_file = NULL;
  ma->utmInfo.utmTripleAvailable = 0;
  ma->inputgrid.utmXYHeaderWordsAvailable = 1;
  glbs.in_utm_conversions=0;
  glbs.utm_point_indx=0;

  ma->eclInfo.useECLIOsegments = 0;
  ma->nLayers = 0;
  initialiseDecorator(&ma->decoratorInfo);

  if (parse(doc,cur,ma,&glbs) == EXIT_FAILURE)
	{
	  xmlFreeDoc(doc);
	  return EXIT_FAILURE;
	}
  xmlFreeDoc(doc);


  // Post-parse warnings and checks


  // Suitable coordinate warnings
  // parse() above will have detected utmTripleAvailable by now if possible
  if (((strcmp(ma->inputgrid.xWord,"gx") && strcmp(ma->inputgrid.xWord,"sx")))
	  || (strcmp(ma->inputgrid.yWord,"gy") && strcmp(ma->inputgrid.yWord,"sy")))
	{
	  ma->inputgrid.utmXYHeaderWordsAvailable = 0;
	  printf("Seismic contains no UTM header words {sx,sy,gx,gy}\n");
	  if (!ma->utmInfo.utmTripleAvailable)
		{
		  // no UTM triples either. Will be impossible to map to eclipse grid
		  if (!ma->duplicate)
			{
			  printf("XML has no UTM triples to enable conversions.\n");
			  printf("Impossible to map seismic to eclipse UTM grid. Please\n"
					 "provide UTM triples in xml or x,y seimic header words\n");
			  return EXIT_FAILURE;
			}
		  else
			{
			  printf("XML has no UTM triples to enable conversions.\n"
					 "Duplicating seismic grid, "
					 "so well data MUST be in seismic coordinates!\n");
			}
		}
	  else
		{
		  printf("Using UTM triples in XML to enable conversions: well data MUST be in UTM coordinates!\n");
		}
	}
  // global flags controlling subsequent mapping
  if (ma->utmInfo.utmTripleAvailable && !ma->inputgrid.utmXYHeaderWordsAvailable)
	{
	  ma->utmInfo.performInputGridUTMmapping = 1;
	}
  else
	{
	  ma->utmInfo.performInputGridUTMmapping = 0;
	}
  // output grid may need to be mapped only if duplicated & the input grid was mapped too
  if (ma->utmInfo.performInputGridUTMmapping && ma->duplicate)
	{
	  ma->utmInfo.performOutputGridUTMmapping = 1;
	}
  else  
	{
	  ma->utmInfo.performOutputGridUTMmapping = 0;
	}

  // Special variograms: test if supplied when requested
  for (m=0; m < ma->nPropertyNames; m++)
	{
	  if (ma->properties[m].use_trend_smoothing_special_variogram)
		{
		  if (!(glbs.found_special_trend_variogram))
			{
			  printf("Property %s requested extra trend_smoothing_special_variogram, but this was not "
					 " supplied in XML. Fix and re-run! \n",ma->propertyNames[m]);
			  return EXIT_FAILURE;
			}
		  printf("Special smoothing applied to property %s\n",ma->propertyNames[m]);
		}
	}


  // Eclipse sources
  if (ma->make_ecl_grid)
	{
	  if ((ma->thickness_index < 0) || (ma->d_index < 0))
		{
		  printf("Creating of ECLIPSE files only possible if both thickness and d are 'massaged'\n");
		  return EXIT_FAILURE;  
		}
	}

  if (isEclipseCommonSource(ma) == EXIT_FAILURE) return EXIT_FAILURE;

  if (!ma->duplicate && !ma->eclInfo.commonEclipseGrid)
	{
	   printf("Cannot yet handle case of multiple ECLIPSE input grids. Build a common one! Exiting\n");
	   return EXIT_FAILURE;
	}
  

  if (ma->nDecorations)
	{
	  if (ma->decoratorInfo.decoration_layer_indx==-1)
		{
		  printf("No decorations entry in XML file. Skipping decorations\n");
		  ma->nDecorations = 0;
		}
	  else
		{
		  // The ma->decoratorInfo.*_index variables have to be precomputed, 
		  // and correct for multiple layers

		  // initialisations to detect errors later on */
		  ma->decoratorInfo.thickness_index = -1;
		  ma->decoratorInfo.NG_index = -1;
		  ma->decoratorInfo.d_index = -1;
	  
		  int properties_indx_k;
		  for (m=0; m< ma->nPropertyNames; m++)
			{
			  properties_indx_k = m + ma->decoratorInfo.decoration_layer_indx*ma->nPropertyNames;
			  if (!strcmp(ma->propertyNames[m], "thickness"))
				ma->decoratorInfo.thickness_index = properties_indx_k;
			  if (!strcmp(ma->propertyNames[m], "NG"))
				ma->decoratorInfo.NG_index = properties_indx_k;
			  if (!strcmp(ma->propertyNames[m], "d"))
				ma->decoratorInfo.d_index = properties_indx_k;
			}
		  if (ma->decoratorInfo.thickness_index ==-1) 
			{
			  printf("Decoration requires property thickness to be massaged\n");
			  return EXIT_FAILURE;
			}
		  if (ma->decoratorInfo.NG_index ==-1) 
			{
			  printf("Decoration requires property NG to be massaged\n");
			  return EXIT_FAILURE;
			}
		  if (ma->decoratorInfo.d_index ==-1) 
			{
			  printf("Decoration requires property d to be massaged\n");
			  return EXIT_FAILURE;
			}
		}
	}

  return EXIT_SUCCESS;
}



// recursive parsing function

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define SWAPC(a,b) {cptr=(a);(a)=(b);(b)=cptr;}

int parse(xmlDocPtr doc, xmlNodePtr cur,struct MASSAGER *ma,struct XML_PARSER_GLOBALS *glbs) 
{
  int swap,pass;
  double val;
  xmlChar *ptr;
  xmlChar *ptr2;
  char *cptr;
  xmlXPathObjectPtr result;
  xmlXPathContextPtr context;
  struct VARIOGRAM_2D *variogram_ptr;

  // find number of layers at initialisation
  if (ma->nLayers==0)
	{
	  int i,j,k;

	  context = xmlXPathNewContext(doc);
	  /* keyword searches must be prepended with // */
	  /* layers */
	  //result = xmlXPathEvalExpression((const xmlChar *)"//layer", context);
	  result = xmlXPathEvalExpression((const xmlChar *)"/massaged_inversion/model_description/layer", context);
	  if (xmlXPathNodeSetIsEmpty(result->nodesetval)) 
		{
		  printf("No layers found in XML file.");
		  xmlXPathFreeContext(context);
		  return EXIT_FAILURE;
		}
	  ma->nLayers = result->nodesetval->nodeNr;
	  xmlXPathFreeObject(result); result = NULL;

	  /* properties. Full path here will prevent counting properties entries in well observations */
	  result = xmlXPathEvalExpression((const xmlChar *)"/massaged_inversion/model_description/property_name", context);
	  if (xmlXPathNodeSetIsEmpty(result->nodesetval))
		{
		  printf("No properties found in XML file. FIXME");
		  xmlXPathFreeContext(context);
		  return EXIT_FAILURE;
		}
	  ma->nPropertyNames = result->nodesetval->nodeNr;
	  ma->nProperties = ma->nPropertyNames * ma->nLayers;

	  // allocations
	  if (ma->writeAsciiOut)
		{
		  ma->ascii_out_filenames = 
			(char**)malloc((size_t) ((ma->nLayers) * sizeof(char*)));
		  ma->ascii_out_fp = 
			(FILE**)malloc((size_t) ((ma->nLayers) * sizeof(FILE*)));
		}
	  if (ma->Nrealisations>0)
		{ 
		  ma->ascii_out_realisations_filenames = 
			(char**)malloc((size_t) ((ma->nLayers) * sizeof(char*)));
		  ma->ascii_out_realisations_fp = 
			(FILE**)malloc((size_t) ((ma->nLayers) * sizeof(FILE*)));
		}
	  ma->layerNumbers = (int*)malloc((size_t) (ma->nLayers * sizeof(int)));
	  ma->propertyNames = (char **)malloc((size_t) ((ma->nPropertyNames) * sizeof(char *)));
	  // two extra properties for x,y's
	  ma->properties = (struct PROPERTY *)malloc((size_t) ((ma->nProperties+2) * sizeof(struct PROPERTY)));
	  ma->variograms = (struct VARIOGRAM_2D *)malloc((size_t) ((ma->nLayers) * sizeof(struct VARIOGRAM_2D)));
	  ma->outputgrid = (struct OUTPUTGRID *)malloc((size_t) ((ma->nLayers) * sizeof(struct OUTPUTGRID)));

 
	  // some sensible defaults
	  ma->eclInfo.writeInSitu = 0;
	  for (i=0; i < ma->nLayers; i++)
		{
		  ma->outputgrid[i].filename = NULL;
		  ma->outputgrid[i].layer_name = NULL; 
		  ma->outputgrid[i].ECLIPSE_segments_keyword = NULL;
		  ma->outputgrid[i].ECLIPSE_layer_number = -1;
		  ma->outputgrid[i].ECLIPSE_reference_depth_layer_number = INT_MIN;
		  ma->outputgrid[i].eclOutInfo.outputFilename = NULL;
		  ma->outputgrid[i].eclOutInfo.haveLayers = 0;
		  ma->outputgrid[i].eclOutInfo.statsConstr.useStochasticDepth = 0;
		  ma->outputgrid[i].eclOutInfo.statsConstr.Zmax = D_BIG;
		  ma->outputgrid[i].eclOutInfo.statsConstr.Zmin = -D_BIG;
		}
	  // ptrs to start of well observations
	  ma->well_observations = (struct WELL_OBSERVATIONS *)malloc((size_t) ((ma->nLayers) * sizeof(struct WELL_OBSERVATIONS)));

	  // loop over property names, since we have them at hand
	  // NB these will be permuted to put depth first, if specified
	  int I_to, I_from;
	  int *norms = ivector(0,ma->nPropertyNames-1);
	  int *smooth_across_faults = ivector(0,ma->nPropertyNames-1);
	  int *use_trend_smoothing_special_variogram = ivector(0,ma->nPropertyNames-1);

	  xmlChar *uri;
	  for (i=0; i < ma->nPropertyNames; i++) 
		{ 
		  ptr = xmlNodeListGetString(doc, result->nodesetval->nodeTab[i]->xmlChildrenNode, 1);
		  ma->propertyNames[i] = (char*)xmlStrdup(ptr);
		  xmlFree(ptr);
		  // normalisation attributes
		  uri = xmlGetProp(result->nodesetval->nodeTab[i], (xmlChar*)"normalised");
		  if (uri != NULL && !strcmp((char*)uri,"true"))
			{
			  norms[i] = 1;
			}
		  else norms[i]=0;
		  xmlFree(uri);

		  // krig-across-faults attributes
		  uri = xmlGetProp(result->nodesetval->nodeTab[i], (xmlChar*)"smooth_across_faults");
		  if (uri != NULL && !strcmp((char*)uri,"true"))
			{
			  smooth_across_faults[i] = 1;
			}
		  else smooth_across_faults[i] = 0;
		  xmlFree(uri);

		  // use_trend_smoothing_special_variogram  attributes
		  uri = xmlGetProp(result->nodesetval->nodeTab[i], (xmlChar*)"use_trend_smoothing_special_variogram");
		  if (uri != NULL && !strcmp((char*)uri,"true"))
			{
			  use_trend_smoothing_special_variogram[i] = 1;
			}
		  else use_trend_smoothing_special_variogram[i] = 0;
		  xmlFree(uri);

		  // series of checks
		  /* I think t_base is OK now
		  if (!strcmp(ma->propertyNames[i],"t_base"))
			{
			  printf("Massaging property t_base is not a suitable property to massage (belongs to only last layer)\n");
			  printf("(P.S. Final surface information reconstructable using vp_eff and thickness)\n");
			  return EXIT_FAILURE;
			}
		  */
		  char* name = ma->propertyNames[i];
		  int notNormed = (norms[i]==0) ? 1 : 0;
		  if (!strcmp(name,"phi_s") && notNormed)
			{
			  printf("Massaging property phi_s should be normalised! Switching normalisation ON\n");
			  norms[i] = 1;
			} 
		  if ((!strcmp(name,"NG") || !strcmp(name,"NG_hard")) && notNormed)
			{
			  printf("Massaging property %s should be normalised! Switching normalisation ON\n",name);
			  norms[i] = 1;
			}
		  if (!strncmp(name,"fluid_type_",11) && notNormed)
			{
			  printf("Massaging property fluid_type_* should be normalised! Switching normalisation ON\n");
			  norms[i] = 1;
			}
		  if ((!strcmp(name,"S_lsg") || !strcmp(name,"S_gas") || !strcmp(name,"S_oil")) && notNormed)
			{
			  printf("Massaging property %s should be normalised! Switching normalisation ON\n",name);
			  norms[i] = 1;
			}
		}
	  xmlXPathFreeObject(result); result = NULL;

	  // Re-sort to achieve these aims:
	  // (1) depth, "d", is always first, if specified
	  for (i=1; i < ma->nPropertyNames; i++) 
		{
		  if (!strcmp(ma->propertyNames[i],"d"))
			{
 			  SWAPC(ma->propertyNames[i],ma->propertyNames[0]);
			  SWAP(smooth_across_faults[i],smooth_across_faults[0]);
			  SWAP(norms[i],norms[0]);
			  SWAP(use_trend_smoothing_special_variogram[i],use_trend_smoothing_special_variogram[0]);
			  break;
			}
		}
	  // (2) and thickness comes *before* net-sand and both *before* NG (make NG last)
	  ma->thickness_index = -1;
	  ma->net_sand_index = -1;
	  ma->NG_index = -1;
	  ma->d_index = -1;
	  ma->t_index = -1;
	  ma->t_base_index = -1;
	  ma->vp_eff_index = -1;
	  fixLookupIndices(ma);
	  if ((ma->thickness_index>=0) && (ma->net_sand_index>=0))
		{
		  if (ma->thickness_index > ma->net_sand_index)
			{
			  SWAPC(ma->propertyNames[ma->thickness_index],ma->propertyNames[ma->net_sand_index]);
			  SWAP(smooth_across_faults[ma->thickness_index],smooth_across_faults[ma->net_sand_index]);
			  SWAP(norms[ma->thickness_index],norms[ma->net_sand_index]);
			  SWAP(use_trend_smoothing_special_variogram[ma->thickness_index],use_trend_smoothing_special_variogram[ma->net_sand_index]);
			}
		}
	  // make NG last
	  I_to = ma->nPropertyNames-1;
	  if (ma->NG_index>=0)
		{
		  SWAPC(ma->propertyNames[ma->NG_index],ma->propertyNames[I_to]);
		  SWAP(smooth_across_faults[ma->NG_index],smooth_across_faults[I_to]);
		  SWAP(norms[ma->NG_index],norms[I_to]);
		  SWAP(use_trend_smoothing_special_variogram[ma->NG_index],use_trend_smoothing_special_variogram[I_to]);
		}

	  // (3) Also, make vp_eff last,  or second last if NG is massaged
	  if (ma->NG_index >= 0) I_to--;
	  if (ma->vp_eff_index >= 0)
		{
		  // move vp_eff to location I_to (last,  or second last if NG is massaged)
		  I_from = ma->vp_eff_index;
		  SWAPC(ma->propertyNames[I_from],ma->propertyNames[I_to]);
		  SWAP(smooth_across_faults[I_from],smooth_across_faults[I_to]);
		  SWAP(norms[I_from],norms[I_to]);
		  SWAP(use_trend_smoothing_special_variogram[I_from],use_trend_smoothing_special_variogram[I_to]);
		}

	  // recompute *_index values, since swapping will have rendered them incorrect.
	  fixLookupIndices(ma);
	  
	  // User feedback
	  for (i=0; i < ma->nPropertyNames; i++) 
		{
		  printf("Massaging property: %s, normalised = %d\n",ma->propertyNames[i],norms[i]);
		}
	  if ((ma->thickness_index>=0) && (ma->net_sand_index>=0) && (ma->NG_index>=0))
		{
		  printf("\nWarning: NG simulation computed from net-sand and thickness: "
				 "NG well constraints are best applied directly on \nnet-sand and thickness. "
				 "Conversions to equivalent constraints will be attempted\n");
		}

	  // error check on stats
	  if (ma->acquireStats && (ma->Nrealisations>0)) 
		{
		  if (ma->duplicate && strcmp(ma->propertyNames[0],"d"))
			{
			  printf("WARNING!: No depth info to depth-truncate quantities for statistics. Stats will be over entire volume\n");
			  ma->acquireStats = 2;
			}
		  // check on sorting name
		  printf("Sorting stats by %s\n",ma->statsSortPropName);
		  int sortingNameOk = 0;
		  for (j=0; j < ma->nPropertyNames; j++)
			{ 
			  if (!strcmp(ma->statsSortPropName,ma->propertyNames[j])) sortingNameOk = 1;
			}
		  if (!sortingNameOk)
			{
			  printf("WARNING!: Specified statistics sorting property %s not available in XML! CDFs will be nonsense!\n",
					 ma->statsSortPropName);
			  exit(1); 
			}
		}
			  
			  
	  // and set up the properties for later parsing
	  for (i=k=0; i < ma->nLayers; i++) 
		{
		  for (j=0; j < ma->nPropertyNames; j++,k++)  
			{
			  ma->properties[k].name = ma->propertyNames[j];
			  ma->properties[k].normalised = norms[j];
			  ma->properties[k].smooth_across_faults = smooth_across_faults[j];
			  ma->properties[k].use_trend_smoothing_special_variogram = use_trend_smoothing_special_variogram[j];
			}
		}
	  free_ivector(norms,0,ma->nPropertyNames-1);
	  norms = NULL;
	  free_ivector(smooth_across_faults,0,ma->nPropertyNames-1);
	  smooth_across_faults =NULL;
	  free_ivector(use_trend_smoothing_special_variogram,0,ma->nPropertyNames-1);
	  use_trend_smoothing_special_variogram =NULL;

	  // clean up
	  xmlXPathFreeContext(context);
	}



  /* We don't care what the top level element name is */
  cur = cur->xmlChildrenNode;

  while (cur != NULL) 
	{
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"delivery_massage_analysis_filename")))
		{
		  // massage analysis filename 
		  if (ma->massage_analysis_file == NULL)
			{
			  ma->massage_analysis_file = (char*)xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			}
		}
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"x-coordinate_header_word")))
		{
		  // massage analysis x header word
		  ma->inputgrid.xWord = (char*)xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
		}
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"y-coordinate_header_word")))
		{
		  // massage analysis y header word
		  ma->inputgrid.yWord = (char*)xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
		}
	  /* layer entry */
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"layer")))
		{
		  /* increment layer number */
		  glbs->layer_indx++;

		  /* preset observations index */
		  glbs->observation_indx=-1;

		  {
			int j;
			/* find inversion layer number */

			xmlXPathContextPtr context;
			xmlXPathObjectPtr result;
			int m;

			context = xmlXPathNewContext(doc);
			// relative search
			context->node = cur;
			/* .//string is a relative search   http://www.w3.org/TR/xpath */
			result = xmlXPathEvalExpression((const xmlChar *)".//inversion_layer_number", context);
			if (xmlXPathNodeSetIsEmpty(result->nodesetval))
			  {
				printf("No inversion layer number found in XML file for layer %d\n",glbs->layer_indx);
			  }
			ptr = xmlNodeListGetString(doc, result->nodesetval->nodeTab[0]->xmlChildrenNode, 1);
			ma->layerNumbers[glbs->layer_indx] = atoi((char*)ptr);
			xmlFree(ptr);
			if (ma->layerNumbers[glbs->layer_indx] < 1)
			  {
				printf("Layer number not a positive number in XML file %s\n",ma->massage_xml_file);
				return EXIT_FAILURE;
			  }

			// and set up the layer Numbers for this layer
			for (j=0; j < ma->nPropertyNames; j++)  
			  {
				ma->properties[j+(glbs->layer_indx*ma->nPropertyNames)].layerNumber 
				  = ma->layerNumbers[glbs->layer_indx];
			  }
			xmlXPathFreeObject(result);result=NULL;


			/* count number of well observations for this layer */

			result = xmlXPathEvalExpression((const xmlChar *)".//well_observation", context);
			if (!xmlXPathNodeSetIsEmpty(result->nodesetval))
			  {
				printf("%d observations found in XML file for inversion layer %d\n",
					   result->nodesetval->nodeNr, ma->layerNumbers[glbs->layer_indx]);
			  }
			ma->well_observations[glbs->layer_indx].n = result->nodesetval->nodeNr;
			if (ma->ignore_well_data) 
			  {
				printf("Observations ignored...\n");
				ma->well_observations[glbs->layer_indx].n = 0;
			  }
			m = ma->well_observations[glbs->layer_indx].n;

			allocateObservations(ma,glbs->layer_indx);


			xmlXPathFreeContext(context);

			// recursively parse
			glbs->in_layer=1;
			if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
			glbs->in_layer=0;
		  
			xmlXPathFreeObject(result); result = NULL;
		  }
		}

	  /* UTM conversions */
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"utm-conversions")))
		{
		  /* count number of points */
		  context = xmlXPathNewContext(doc);
		  // relative search
		  context->node = cur;
		  result = xmlXPathEvalExpression((const xmlChar *)".//point", context);
		  if (!xmlXPathNodeSetIsEmpty(result->nodesetval))
			{
			  printf("%d utm points found in XML file\n",
					 result->nodesetval->nodeNr);
			}
		  if (result->nodesetval->nodeNr !=3)
			{
			  printf("please provide only 3 distinct non-collinear points for utm geometry conversion\n");
			  return EXIT_FAILURE;
			}
		  ma->utmInfo.utmTripleAvailable = 1;
		  xmlXPathFreeContext(context);
		  // recursively get bits
		  glbs->in_utm_conversions = 1;
		  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_utm_conversions = 0;
		}
	  if (glbs->in_utm_conversions)
		{
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"x")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  ma->utmInfo.x[glbs->utm_point_indx] = atof((char*)ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"y")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  ma->utmInfo.y[glbs->utm_point_indx] = atof((char*)ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"inline")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  ma->utmInfo.inLine[glbs->utm_point_indx] = atof((char*)ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"xline")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  ma->utmInfo.xLine[glbs->utm_point_indx] = atof((char*)ptr);
			  xmlFree(ptr);
			  glbs->utm_point_indx++;
			}
		}

	  /* Variogram stuff */

	  /* property variograms */
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"variogram")) && glbs->in_layer)
		{
		  // recursively get bits
		  glbs->in_variogram = 1;
		  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_variogram = 0;
		}
	  /* trend smoothing variogram */
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"trend_smoothing_variogram")))
		{
		  // recursively get bits
		  glbs->in_trend_variogram = 1;
		  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_trend_variogram = 0;
		}
	  /* special trend smoothing variogram */
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"trend_smoothing_special_variogram")))
		{
		  // recursively get bits
		  glbs->in_special_trend_variogram = 1;
		  glbs->found_special_trend_variogram = 1;
		  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_special_trend_variogram = 0;
		}
	  if (glbs->in_variogram  || glbs->in_trend_variogram || glbs->in_special_trend_variogram)
		{
		  if (glbs->in_variogram)
			{
			  variogram_ptr = &(ma->variograms[glbs->layer_indx]);
			}
		  else if (glbs->in_trend_variogram)
			{
			  variogram_ptr = &(ma->trend_smoothing_variogram);
			}
		  else 
			{
			  variogram_ptr = &(ma->trend_smoothing_special_variogram);
			}
		  /* variogram bits */
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"type")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  // (a) type
			  if (!strncmp((char*)ptr,"exp",3))
				{
				  variogram_ptr->GSLIB_type = 2;
				}
			  else if (!strncmp((char*)ptr,"sph",3))
				{
				  variogram_ptr->GSLIB_type = 1;
				}
			  else if (!strncmp((char*)ptr,"gau",3))
				{
				  variogram_ptr->GSLIB_type = 3;
				}
			  else 
				{
				  printf("variogram type must be exp,sph, or gau\n");
				  return EXIT_FAILURE;
				}
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"sill")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  variogram_ptr->sill = atof((char*)ptr);
			  xmlFree(ptr);
			  if (variogram_ptr->sill <= 0)
				{
				  printf("variogram sill must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"nugget")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  variogram_ptr->nugget = atof((char*)ptr);
			  xmlFree(ptr);
			  if (variogram_ptr->nugget < 0)
				{
				  printf("variogram nugget  must be >= 0 \n");
				  return EXIT_FAILURE;
				}
			  // Fix dodgy nugget on Gaussians if 0
			  if (variogram_ptr->GSLIB_type == 3) 
				{
				  if (variogram_ptr->nugget < variogram_ptr->sill * NUGGET_MIN)
					{
					  variogram_ptr->nugget = variogram_ptr->sill * NUGGET_MIN;
					  printf("Cleaned Gaussian nugget to %g\n",variogram_ptr->nugget);
					}
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"azimuth")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  variogram_ptr->azimuth = atof((char*)ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"rotated-x-range")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  variogram_ptr->rangeX = atof((char*)ptr);
			  xmlFree(ptr);
			  if (variogram_ptr->rangeX <= 0)
				{
				  printf("variogram rotated-x-range must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"rotated-y-range")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  variogram_ptr->rangeY = atof((char*)ptr);
			  xmlFree(ptr);
			  if (variogram_ptr->rangeY <= 0)
				{
				  printf("variogram rotated-y-range must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		}


	  if (glbs->in_layer)
		{
		  /* layer name */
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"name")))
			{
			  ptr = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			  ma->outputgrid[glbs->layer_indx].layer_name = (char*) xmlStrdup(ptr);
			  xmlFree(ptr);
			}

		  /* Well observations stuff */

		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"well_observation")) && (!ma->ignore_well_data))
			{
			  // recursively get bits
			  glbs->in_observation = 1;
			  glbs->observation_indx++;
			  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
			  glbs->in_observation = 0;
			}
		  if (glbs->in_observation == 1)
			{
			  /* well observation bits */
			  struct WELL_OBSERVATIONS *wo = &ma->well_observations[glbs->layer_indx];

			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"position_x")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  val = atof((char*)ptr);
				  g_array_insert_val(wo->xValues, glbs->observation_indx, val);

				  // wo->xValues[glbs->observation_indx] = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"position_y")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  val = atof((char*)ptr);
				  g_array_insert_val(wo->yValues, glbs->observation_indx,val);
				  // wo->yValues[glbs->observation_indx] = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"property_value")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  val = atof((char*)ptr);
				  g_array_insert_val(wo->values, glbs->observation_indx, val);
				  // wo->values[glbs->observation_indx] = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"property_name")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  // wo->property_name[glbs->observation_indx] = ptr;				  
				  ptr2 = xmlStrdup(ptr);
				  g_array_insert_val(wo->property_name, glbs->observation_indx, ptr2);
				  xmlFree(ptr);
				}
			}

		  if (!ma->duplicate)
			{
			  // ECLIPSE grid filename
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_grid_filename")))
				{
				  ptr = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
				  ma->outputgrid[glbs->layer_indx].filename = (char*)xmlStrdup(ptr);
				  xmlFree(ptr);
				}
			  // ECLIPSE grid segments if available filename
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_segments_keyword")))
				{
				  ma->eclInfo.useECLIOsegments = 1;
				  ma->outputgrid[glbs->layer_indx].ECLIPSE_segments_keyword = (char*)xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
				}
			  // ECLIPSE layer number if available
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_layer_number")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  ma->outputgrid[glbs->layer_indx].ECLIPSE_layer_number = atoi((char*)ptr);
				  xmlFree(ptr);
				}
			  // ECLIPSE ref-depth layer number if available
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_reference_depth_layer_number")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  ma->outputgrid[glbs->layer_indx].ECLIPSE_reference_depth_layer_number = atoi((char*)ptr);
				  xmlFree(ptr);
				}
			  // Output (ECLIPSE) grid filename if writing out
			  if (ma->eclInfo.writeECLOut)
				{
				  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_output_grid_filename")))
					{
					  ma->outputgrid[glbs->layer_indx].eclOutInfo.outputFilename = (char*)xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
					}
				}
			} 
		  // stats
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"stats")))
			{
			  // recursively get bits
			  glbs->in_stats = 1;
			  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
			  glbs->in_stats = 0;
			}
		  if (glbs->in_stats == 1)
			{
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"truncation_Zmax")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  ma->outputgrid[glbs->layer_indx].eclOutInfo.statsConstr.Zmax = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"truncation_Zmin")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  ma->outputgrid[glbs->layer_indx].eclOutInfo.statsConstr.Zmin = atof((char*)ptr);
				  xmlFree(ptr);
				}
			}

		  // decoration
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"decorations")))
			{
			  // recursively get bits
			  glbs->in_decoration = 1;
			  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
			  glbs->in_decoration = 0;
			}
		  if (glbs->in_decoration == 1)
			{
			  struct DECORATOR_INFO *decoratorInfo = &ma->decoratorInfo;
			  decoratorInfo->decoration_layer = ma->layerNumbers[glbs->layer_indx];
			  decoratorInfo->decoration_layer_indx = glbs->layer_indx;
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"Nz")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->Nz = atoi((char*)ptr);
				  if (decoratorInfo->Nz<1)
					{
					  printf("Must have decoration Nz>=1\n");
					  return EXIT_FAILURE;
					}
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"type")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  // (a) type
				  if (!strncmp((char*)ptr,"exp",3))
					{
					  decoratorInfo->variogram.GSLIB_type = 2;
					}
				  else if (!strncmp((char*)ptr,"sph",3))
					{
					  decoratorInfo->variogram.GSLIB_type = 1;
					}
				  else if (!strncmp((char*)ptr,"gau",3))
					{
					  decoratorInfo->variogram.GSLIB_type = 3;
					}
				  else 
					{
					  printf("decorator variogram type must be exp,sph, or gau\n");
					  return EXIT_FAILURE;
					}
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"sill")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.sill = atof((char*)ptr);
				  xmlFree(ptr);
				  if (decoratorInfo->variogram.sill <= 0)
					{
					  printf("decorator variogram sill must be positive\n");
					  return EXIT_FAILURE;
					}
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"nugget")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.nugget = atof((char*)ptr);
				  xmlFree(ptr);
				  if (decoratorInfo->variogram.nugget < 0)
					{
					  printf("decorator variogram nugget  must be >= 0 \n");
					  return EXIT_FAILURE;
					}
				  // Fix dodgy nugget on Gaussians if 0
				  if (decoratorInfo->variogram.GSLIB_type == 3) 
					{
					  if (decoratorInfo->variogram.nugget < decoratorInfo->variogram.sill * NUGGET_MIN)
						{
						  decoratorInfo->variogram.nugget = decoratorInfo->variogram.sill * NUGGET_MIN;
						  printf("Cleaned Gaussian nugget to %g\n",decoratorInfo->variogram.nugget);
						}
					}
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"a_hmax")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.a_hmax = atof((char*)ptr);
				  xmlFree(ptr);
				  if (decoratorInfo->variogram.a_hmax <= 0)
					{
					  printf("decorator variogram a_hmax must be positive\n");
					  return EXIT_FAILURE;
					}
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"a_hmin")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.a_hmin = atof((char*)ptr);
				  xmlFree(ptr);
				  if (decoratorInfo->variogram.a_hmin <= 0)
					{
					  printf("decorator variogram a_hmin must be positive\n");
					  return EXIT_FAILURE;
					}
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"a_vert")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.a_vert = atof((char*)ptr);
				  xmlFree(ptr);
				  if (decoratorInfo->variogram.a_vert <= 0)
					{
					  printf("decorator variogram a_vert must be positive\n");
					  return EXIT_FAILURE;
					}
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ang1")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.ang1 = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ang2")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.ang2 = atof((char*)ptr);
				  xmlFree(ptr);
				}
			  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ang3")))
				{
				  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
				  decoratorInfo->variogram.ang3 = atof((char*)ptr);
				  xmlFree(ptr);
				}
			}
		} // in layer


	  // recursion
	  if (parse(doc,cur,ma,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
	  cur = cur->next;
	}
  return EXIT_SUCCESS;
}

void fixLookupIndices(struct MASSAGER *ma)
{
  int i;
  for (i=0; i < ma->nPropertyNames; i++) 
	{
	  if (!strcmp(ma->propertyNames[i],"thickness")) ma->thickness_index=i;
	  if (!strcmp(ma->propertyNames[i],"net-sand")) ma->net_sand_index=i;
	  if (!strcmp(ma->propertyNames[i],"NG")) ma->NG_index=i;
	  if (!strcmp(ma->propertyNames[i],"d")) ma->d_index=i;
	  if (!strcmp(ma->propertyNames[i],"t")) ma->t_index=i;
	  if (!strcmp(ma->propertyNames[i],"t_base")) ma->t_base_index=i;
	  if (!strcmp(ma->propertyNames[i],"vp_eff")) ma->vp_eff_index=i;
	}
}

xmlChar* xmlNodeListGetStringWrapper(xmlDocPtr doc, xmlNodePtr cur)
{
  xmlChar *ptr;
  if (cur->xmlChildrenNode == NULL)
	{
	  printf("Parsing error for xml entry %s\n",cur->name );
	}
  ptr = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);

  return ptr;
}

#undef SWAP
#undef SWAPC
