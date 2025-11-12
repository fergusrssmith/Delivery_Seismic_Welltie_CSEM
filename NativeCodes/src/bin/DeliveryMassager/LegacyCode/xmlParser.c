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

#include "NR/nrutil.h"

#include "kdtree/optkd.h"
#include "bayes/bayes.h"
#include "massager_types.h"
#include "xmlParser.h"

int parseXMLfile(struct MASSAGER *massager)
{
  struct XML_PARSER_GLOBALS glbs;
  xmlDocPtr doc; /* the document tree */
  xmlNodePtr cur; 

  doc = xmlParseFile(massager->massage_xml_file);
  if (doc == NULL) 
	{
	  fprintf(stderr, "Failed to parse %s\n", massager->massage_xml_file);
	  return EXIT_FAILURE;
    }
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) 
	{
	  fprintf(stderr,"empty document\n");
	  xmlFreeDoc(doc);
	  return EXIT_FAILURE;
	}

  // initialise some useful flags
  glbs.layer_indx=-1;
  glbs.in_variogram=0;
  glbs.in_observation=0;
  massager->utmInfo.useUTMconversion = 0;
  glbs.in_utm_conversions=0;
  glbs.utm_point_indx=0;

  massager->nLayers = 0;
  if (massager->writeAsciiOut)
	{
	  massager->ascii_out_filenames = 
		(char**)malloc((size_t) ((massager->nLayers) * sizeof(char*)));
	  massager->ascii_out_fp = 
		(FILE**)malloc((size_t) ((massager->nLayers) * sizeof(FILE*)));
	}
  if (massager->Nrealisations>0)
	{
	  massager->ascii_out_realisations_filenames = 
		(char**)malloc((size_t) ((massager->nLayers) * sizeof(char*)));
	  massager->ascii_out_realisations_fp = 
		(FILE**)malloc((size_t) ((massager->nLayers) * sizeof(FILE*)));
	}
  if (parse(doc,cur,massager,&glbs) == EXIT_FAILURE)
	{
	  xmlFreeDoc(doc);
	  return EXIT_FAILURE;
	}
  xmlFreeDoc(doc);

  // Suitable coordinate warnings
  if (((strcmp(massager->inputgrid.xWord,"gx") && strcmp(massager->inputgrid.xWord,"sx")))
	  || (strcmp(massager->inputgrid.yWord,"gy") && strcmp(massager->inputgrid.yWord,"sy")))
	{
	  printf("Seismic contains no UTM header words {sx,sy,gx,gy}\n");
	  if (!massager->utmInfo.useUTMconversion)
		{
		  printf("And XML has no UTM triples to enable conversions. Well data MUST be in seismic coordinates!\n");
		}
	  else
		{
		  printf("Using UTM triples in XML to enable conversions: well data MUST be in UTM coordinates!\n");
		}
	}

  return EXIT_SUCCESS;
}



// recursive parsing function
int parse(xmlDocPtr doc, xmlNodePtr cur,struct MASSAGER *massager,struct XML_PARSER_GLOBALS *glbs) 
{
  xmlChar *ptr;
  xmlXPathObjectPtr result;
  xmlXPathContextPtr context;

  // find number of layers at initialisation
  if (massager->nLayers==0)
	{
	  int i,j,k;

	  context = xmlXPathNewContext(doc);
	  /* keyword searches must be prepended with // */
	  /* layers */
	  result = xmlXPathEvalExpression((const xmlChar *)"//layer", context);
	  if (xmlXPathNodeSetIsEmpty(result->nodesetval))
		{
		  printf("No layers found in XML file.");
		  xmlXPathFreeContext(context);
		  return EXIT_FAILURE;
		}
	  massager->nLayers = result->nodesetval->nodeNr;

	  /* properties. Full path here will prevent counting properties entries in well observations */
	  result = xmlXPathEvalExpression((const xmlChar *)"/massaged_inversion/model_description/property_name", context);
	  if (xmlXPathNodeSetIsEmpty(result->nodesetval))
		{
		  printf("No properties found in XML file. FIXME");
		  xmlXPathFreeContext(context);
		  return EXIT_FAILURE;
		}
	  massager->nPropertyNames = result->nodesetval->nodeNr;
	  massager->nProperties = massager->nPropertyNames * massager->nLayers;

	  // allocation
	  massager->layerNumbers = (int*)malloc((size_t) (massager->nLayers * sizeof(int)));
	  massager->propertyNames = (char **)malloc((size_t) ((massager->nPropertyNames) * sizeof(char *)));
	  // two extra properties for x,y's
	  massager->properties = (struct PROPERTY *)malloc((size_t) ((massager->nProperties+2) * sizeof(struct PROPERTY)));
	  massager->variograms = (struct VARIOGRAM_2D *)malloc((size_t) ((massager->nLayers) * sizeof(struct VARIOGRAM_2D)));
	  massager->outputgrid = (struct OUTPUTGRID *)malloc((size_t) ((massager->nLayers) * sizeof(struct OUTPUTGRID)));

	  // ptrs to start of well observations
	  massager->well_observations = (struct WELL_OBSERVATIONS *)malloc((size_t) ((massager->nLayers) * sizeof(struct WELL_OBSERVATIONS)));

	  // loop over property names, sice we have them at hand
	  int *norms = ivector(0,massager->nPropertyNames-1);

	  xmlChar *uri;
	  for (i=0; i < massager->nPropertyNames; i++) 
		{
		  ptr = xmlNodeListGetString(doc, result->nodesetval->nodeTab[i]->xmlChildrenNode, 1);
		  massager->propertyNames[i] = ptr;
		  // normalisation attributes
		  uri = xmlGetProp(result->nodesetval->nodeTab[i], "normalised");
		  if (uri !=NULL && !strcmp(uri,"true"))
			{
			  norms[i] = 1;
			}
		  else norms[i]=0;
		  xmlFree(uri);
		}
	  // and set up the properties for later parsing
	  for (i=k=0; i < massager->nLayers; i++) 
		{
		  for (j=0; j < massager->nPropertyNames; j++,k++)  
			{
			  massager->properties[k].name = massager->propertyNames[j];
			  massager->properties[k].normalised = norms[i];
			}
		}
	  free_ivector(norms,0,massager->nPropertyNames-1);

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
		  massager->massage_analysis_file = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
		}
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"x-coordinate_header_word")))
		{
		  // massage analysis x header word
		  massager->inputgrid.xWord = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
		}
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"y-coordinate_header_word")))
		{
		  // massage analysis y header word
		  massager->inputgrid.yWord = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
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
			massager->layerNumbers[glbs->layer_indx] = atoi(ptr);
			if (massager->layerNumbers[glbs->layer_indx] < 1)
			  {
				printf("Layer number not a positive number in XML file %s\n",massager->massage_xml_file);
				return EXIT_FAILURE;
			  }

			// and set up the layer Numbers for this layer
			for (j=0; j < massager->nPropertyNames; j++)  
			  {
				massager->properties[j+(glbs->layer_indx*massager->nPropertyNames)].layerNumber 
				  = massager->layerNumbers[glbs->layer_indx];
			  }



			/* count number of well observations for this layer */

			result = xmlXPathEvalExpression((const xmlChar *)".//well_observation", context);
			if (!xmlXPathNodeSetIsEmpty(result->nodesetval))
			  {
				printf("%d observations found in XML file for inversion layer %d\n",
					   result->nodesetval->nodeNr, massager->layerNumbers[glbs->layer_indx]);
			  }
			massager->well_observations[glbs->layer_indx].n = result->nodesetval->nodeNr;
			m = massager->well_observations[glbs->layer_indx].n;
			if (massager->well_observations[glbs->layer_indx].n > 0)
			  {
				massager->well_observations[glbs->layer_indx].xValues = 
				  dvector(0,m-1);
				massager->well_observations[glbs->layer_indx].yValues = 
				  dvector(0,m-1);
				massager->well_observations[glbs->layer_indx].values = 
				  dvector(0,m-1);
				massager->well_observations[glbs->layer_indx].propertyIndex = 
				  ivector(0,m-1);
				massager->well_observations[glbs->layer_indx].property_name = 
				  (char**)malloc((size_t) (m * sizeof(char*)));
			  }
			xmlXPathFreeContext(context);
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
		  massager->utmInfo.useUTMconversion = 1;
		  xmlXPathFreeContext(context);
		  // recursively get bits
		  glbs->in_utm_conversions = 1;
		  if (parse(doc,cur,massager,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_utm_conversions = 0;
		}
	  if (glbs->in_utm_conversions)
		{
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"x")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->utmInfo.x[glbs->utm_point_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"y")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->utmInfo.y[glbs->utm_point_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"inline")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->utmInfo.inLine[glbs->utm_point_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"xline")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->utmInfo.xLine[glbs->utm_point_indx] = atof(ptr);
			  xmlFree(ptr);
			  glbs->utm_point_indx++;
			}
		}

	  /* Variogram stuff */

	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"variogram")))
		{
		  // recursively get bits
		  glbs->in_variogram = 1;
		  if (parse(doc,cur,massager,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_variogram = 0;
		}
	  if (glbs->in_variogram == 1)
		{
		  /* variogram bits */
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"type")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  // (a) type
			  if (!strncmp(ptr,"exp",3))
				{
				  massager->variograms[glbs->layer_indx].GSLIB_type = 2;
				}
			  else if (!strncmp(ptr,"sph",3))
				{
				  massager->variograms[glbs->layer_indx].GSLIB_type = 1;
				}
			  else if (!strncmp(ptr,"gau",3))
				{
				  massager->variograms[glbs->layer_indx].GSLIB_type = 3;
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
			  massager->variograms[glbs->layer_indx].sill = atof(ptr);
			  xmlFree(ptr);
			  if (massager->variograms[glbs->layer_indx].sill <= 0)
				{
				  printf("variogram sill must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"nugget")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->variograms[glbs->layer_indx].nugget = atof(ptr);
			  xmlFree(ptr);
			  if (massager->variograms[glbs->layer_indx].nugget < 0)
				{
				  printf("variogram nugget  must be >= 0 \n");
				  return EXIT_FAILURE;
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"azimuth")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->variograms[glbs->layer_indx].azimuth = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"rotated-x-range")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  ptr = xmlNodeListGetStringWrapper(doc, cur);
			  massager->variograms[glbs->layer_indx].rangeX = atof(ptr);
			  xmlFree(ptr);
			  if (massager->variograms[glbs->layer_indx].rangeX <= 0)
				{
				  printf("variogram rotated-x-range must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"rotated-y-range")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->variograms[glbs->layer_indx].rangeY = atof(ptr);
			  xmlFree(ptr);
			  if (massager->variograms[glbs->layer_indx].rangeY <= 0)
				{
				  printf("variogram rotated-y-range must be positive\n");
				  return EXIT_FAILURE;
				}
			}
		}

	  /* Well observations stuff */

	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"well_observation")))
		{
		  // recursively get bits
		  glbs->in_observation = 1;
		  glbs->observation_indx++;
		  if (parse(doc,cur,massager,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
		  glbs->in_observation = 0;
		}
	  if (glbs->in_observation == 1)
		{
		  /* well observation bits */
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"position_x")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->well_observations[glbs->layer_indx].xValues[glbs->observation_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"position_y")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->well_observations[glbs->layer_indx].yValues[glbs->observation_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"property_value")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->well_observations[glbs->layer_indx].values[glbs->observation_indx] = atof(ptr);
			  xmlFree(ptr);
			}
		  if ((!xmlStrcmp(cur->name,(const xmlChar *)"property_name")))
			{
			  if (!(ptr = xmlNodeListGetStringWrapper(doc, cur))) return EXIT_FAILURE;
			  massager->well_observations[glbs->layer_indx].property_name[glbs->observation_indx] = ptr;
			}
		}


	  // Output (ECLIPSE) grid filename
	  if ((!xmlStrcmp(cur->name,(const xmlChar *)"ECLIPSE_grid_filename")))
		{
		  massager->outputgrid[glbs->layer_indx].filename = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
		}

	  // recursion
	  if (parse(doc,cur,massager,glbs) == EXIT_FAILURE) return EXIT_FAILURE;
	  cur = cur->next;
	}
  return EXIT_SUCCESS;
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
