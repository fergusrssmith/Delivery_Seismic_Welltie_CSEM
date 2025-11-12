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

#include "kdtree/optkd.h"
#include "bayes/bayes.h"
#include "massager_types.h"
#include "variogram.h"

/*
  Set up efficiency variables for extensive use of corr2D
  with a fixed variogram.
*/
void setUpcorr2D(struct VARIOGRAM_2D *variogram)
{
  double tmp;
  tmp = M_PI*variogram->azimuth/180.0;
  variogram->costh = cos(tmp);
  variogram->sinth = sin(tmp); 
  /* maybe a lookup table later for exp() etc....*/
}

/*
  Return the 2D covariance of these two points, as specified
  by the variogram. Initialise with setUpcorr2D().

  Conventions are as per GSLIB: x = easting, y=Northing,
  azimuth = clockwise.

  FIXME: complete for all standard variogram types.......
*/
double corr2D(double x1,double y1,double x2,double y2,struct VARIOGRAM_2D *variogram)
{
  double dr,dx,dy,dxp,dyp;
  double nugget=0.0;

  dx = x1-x2;
  dy = y1-y2;

  dxp = (variogram->costh*dx - variogram->sinth*dy) / variogram->rangeX;
  dyp = (variogram->sinth*dx + variogram->costh*dy) / variogram->rangeY;

  if (dx==0 && dy==0) nugget = variogram->nugget;

  switch (variogram->GSLIB_type)
	{
	case 1:
	  {
		/* spherical */
		{
		  double x;
		  x = sqrt(dxp*dxp + dyp*dyp);
		  return nugget+variogram->sill * ((x<=1.0) ? (1.0-1.5*x+0.5*x*x) : 0);
		}
		break;
	  }
	case 2:
	  {
		/* exponential */
		return nugget+variogram->sill * exp(-3.0*sqrt(dxp*dxp + dyp*dyp));
		break;
	  }
	case 3:
	  {
		/* Gaussian */
		return nugget+variogram->sill * exp(-3.0*(dxp*dxp + dyp*dyp));
		break;
	  }
	}
  // never gets here
  return 0;
}
