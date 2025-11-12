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

#include "std/std.h"
#include "std/error.h"
#include "kdtree/optkd.h"

#include "bayes/bayes.h"
#include "massager_types.h"
#include "inputgrid.h"

int setUpInputGrid(struct INPUTGRID *ig)
{
  double norm,dotP;

  // set up input grid parameters
  if (0)
	{
	  // Obsolete stuff assuming input grid is rastered and complete

	  if (ig->nPoints<2)
		{
		  warning("Too few points in input grid");
		  return EXIT_FAILURE;
		}
	  ig->rFirst.x = ig->xValues[0];
	  ig->rFirst.y = ig->yValues[0];
	  ig->rLast.x = ig->xValues[ig->nPoints-1];
	  ig->rLast.y = ig->yValues[ig->nPoints-1];
	  // square of norm of vector joining first two points
	  norm = (ig->xValues[1] -ig->xValues[0])*(ig->xValues[1] -ig->xValues[0])
	    +(ig->yValues[1] -ig->yValues[0])*(ig->yValues[1] -ig->yValues[0]);
	  ig->dRx.x = (ig->xValues[1] - ig->rFirst.x)/norm;
	  ig->dRx.y = (ig->yValues[1] - ig->rFirst.y)/norm;
	  dotP = (ig->rLast.x - ig->rFirst.x)*ig->dRx.x + (ig->rLast.y - ig->rFirst.y)*ig->dRx.y;
	  ig->Nx = 1+rint(dotP);  // number of x values along raster
	  if (ig->Nx<1)
		{
		  warning("input grid is nonrectangular");
		}
	  if (ig->Nx < (ig->nPoints))
		{
		  // looks like a full rectangular grid. OK
		  norm = (ig->xValues[ig->Nx] - ig->xValues[0])*(ig->xValues[ig->Nx] - ig->xValues[0])
	        +(ig->yValues[ig->Nx] - ig->yValues[0])*(ig->yValues[ig->Nx] - ig->yValues[0]);
		  ig->dRy.x = (ig->xValues[ig->Nx] - ig->rFirst.x)/norm;
		  ig->dRy.y = (ig->yValues[ig->Nx] - ig->rFirst.y)/norm;
		  dotP = (ig->rLast.x - ig->rFirst.x)*ig->dRy.x + (ig->rLast.y - ig->rFirst.y)*ig->dRy.y;
		  ig->Ny = rint(dotP);  // number of y values along raster
		}
	  else
		{
		  ig->Ny=1; 
		  ig->dRy.x=0;
		  ig->dRy.y=0;
		}
	  // now we can get bounding indices of any x,y point with the getBoundingIndices method
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
