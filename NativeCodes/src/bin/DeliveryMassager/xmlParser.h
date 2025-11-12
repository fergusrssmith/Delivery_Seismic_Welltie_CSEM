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

#ifndef _xmlParser_h_
#define _xmlParser_h_

struct XML_PARSER_GLOBALS 
{
  int layer_indx;
  char in_layer;
  char in_trend_variogram;
  char in_special_trend_variogram;
  char found_special_trend_variogram;
  char in_variogram;
  int observation_indx;
  char in_observation;
  char in_utm_conversions;
  int utm_point_indx;
  char in_stats;
  char in_decoration;
};
int parseXMLfile(struct MASSAGER*);
int parse(xmlDocPtr, xmlNodePtr, struct MASSAGER*, struct XML_PARSER_GLOBALS*);
xmlChar* xmlNodeListGetStringWrapper(xmlDocPtr doc, xmlNodePtr cur);
void fixLookupIndices(struct MASSAGER *ma);
#endif
