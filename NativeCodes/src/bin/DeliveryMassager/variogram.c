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
#include "kdtree/dynamic-tree.h"
#include "massager_types.h"
#include "variogram.h"


static double DEG2RAD = 0.01745329251994;   /* pi/180 */
static double EPSLON = 1e-20;   /* eps */
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

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
  if (variogram->GSLIB_type==3) // Gaussian
	{
	  // safeguard with small nugget
	  if (variogram->nugget/variogram->sill < NUGGET_MIN)
		{
		  variogram->nugget = variogram->sill * NUGGET_MIN;
		}
	}
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

  if ((x1==x2) && (y1==y2)) nugget = variogram->nugget;

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

/*
  Return the 2D covariance of these two points, as specified
  by the variogram, but assuming exponential type. Safer for generic
  interpolation weights. Coordinate transformations/rotations etc same as for
  corr2D(), i.e.  Initialise with setUpcorr2D().

  Conventions are as per GSLIB: x = easting, y=Northing,
  azimuth = clockwise.

  FIXME: complete for all standard variogram types.......
*/
double corr2Dexp(double x1,double y1,double x2,double y2,struct VARIOGRAM_2D *variogram)
{
  double dr,dx,dy,dxp,dyp;
  double nugget=0.0;

  dx = x1-x2;
  dy = y1-y2;

  dxp = (variogram->costh*dx - variogram->sinth*dy) / variogram->rangeX;
  dyp = (variogram->sinth*dx + variogram->costh*dy) / variogram->rangeY;

  if ((x1==x2) && (y1==y2)) nugget = variogram->nugget;

  /* exponential */
  return nugget+variogram->sill * exp(-3.0*sqrt(dxp*dxp + dyp*dyp));
}




/*
  Set up efficiency variables for extensive use of corr3D
  with a fixed variogram.
*/
void setUpcorr3D(struct VARIOGRAM_3D *variogram)
{
  int i,j;
  double min,max,sum,aa1,aa2,aa3;
  double afac1,afac2,sina,sinb,sint,cosa,cosb,cost;
  double anis1,anis2,alpha,beta,theta;

  aa1 = fabs(variogram->a_hmax);
  aa2 = fabs(variogram->a_hmin);
  aa3 = fabs(variogram->a_vert);
  max = MAX(MAX(aa1,aa2),aa3);

  /* compute equivalents to GSLIB's anis1 anis2 parameters */
  anis1 = aa2/max;
  anis2 = aa3/max;

  /* now continue on with usual GSLIB stuff.... */
  //   variogram->ang1             Azimuth angle for principal direction
  //   variogram->ang2             Dip angle for principal direction
  //   variogram->ang3             Third rotation angle
  //   anis1            First anisotropy ratio
  //   anis2            Second anisotropy ratio
  //
  // Converts the input angles to three angles which make more
  //  mathematical sense:
  //
  //         alpha   angle between the major axis of anisotropy and the
  //                 E-W axis. Note: Counter clockwise is positive.
  //         beta    angle between major axis and the horizontal plane.
  //                 (The dip of the ellipsoid measured positive down)
  //         theta   Angle of rotation of minor axis about the major axis
  //                 of the ellipsoid.
  //
  if (variogram->ang1>=0.0 && variogram->ang1<270.0)
	{
	  alpha = (90.0   - variogram->ang1) * DEG2RAD;
	}
  else
	{
	  alpha = (450.0  - variogram->ang1) * DEG2RAD;
	}
  beta  = -1.0 * variogram->ang2 * DEG2RAD;
  theta =        variogram->ang3 * DEG2RAD;
  //
  // Get the required sines and cosines:
  //
  sina  = sin(alpha);
  sinb  = sin(beta);
  sint  = sin(theta);
  cosa  = cos(alpha);
  cosb  = cos(beta);
  cost  = cos(theta);
  //
  // Construct the rotation matrix in the required memory:
  //
  afac1 = 1.0 / MAX(anis1,EPSLON);
  afac2 = 1.0 / MAX(anis2,EPSLON);
  variogram->rotmat[0][0] =       (cosb * cosa);
  variogram->rotmat[0][1] =       (cosb * sina);
  variogram->rotmat[0][2] =       (-sinb);
  variogram->rotmat[1][0] = afac1*(-cost*sina + sint*sinb*cosa);
  variogram->rotmat[1][1] = afac1*(cost*cosa + sint*sinb*sina);
  variogram->rotmat[1][2] = afac1*( sint * cosb);
  variogram->rotmat[2][0] = afac2*(sint*sina + cost*sinb*cosa);
  variogram->rotmat[2][1] = afac2*(-sint*cosa + cost*sinb*sina);
  variogram->rotmat[2][2] = afac2*(cost * cosb);

  // Then, the quadratic distance is ||rotmat.dr||^2 = dr^T.rotmat^T.rotmat.dr/aa1^2
  // as per the GSLIB function sqdist.
  for (i=0;i<3;i++) {
	for (j=0;j<3;j++) {
	  variogram->rotmat[i][j]/=aa1;
	}
  }
}

/*
  Return the 3D covariance of these two points, as specified
  by the variogram. Initialise with setUpcorr3D().

  Conventions are as per GSLIB: x = Easting, y=Northing, z=depth
  azimuth = clockwise.

  FIXME: complete for all standard variogram types.......
*/
double corr3D(double x1,double y1,double z1,double x2,double y2,double z2,struct VARIOGRAM_3D *variogram)
{
  int i,j;
  double dr[3],da[3],ds;
  double nugget=0.0;

  dr[0] = x2-x1;
  dr[1] = y2-y1;
  dr[2] = z2-z1;

  // quadratic distance is ||rotmat.dr||^2 = dr^T.rotmat^T.rotmat.dr/aa1^2
  for (i=0,ds=0.0;i<3;i++)
	{
	  for (j=0,da[i]=0;j<3;j++) da[i]+=variogram->rotmat[i][j]*dr[j];
	  ds += da[i]*da[i];
	}
  
  if ((x1==x2) && (y1==y2) && (z1==z2)) nugget = variogram->nugget;

  switch (variogram->GSLIB_type)
	{
	case 1:
	  {
		/* spherical */
		{
		  double x;
		  x = sqrt(ds);
		  return nugget+variogram->sill * ((x<=1.0) ? (1.0-1.5*x+0.5*x*x) : 0);
		}
		break;
	  }
	case 2:
	  {
		/* exponential */
		return nugget+variogram->sill * exp(-3.0*sqrt(ds));
		break;
	  }
	case 3:
	  {
		/* Gaussian */
		return nugget+variogram->sill * exp(-3.0*(ds));
		break;
	  }
	}
  // never gets here
  return 0;
}
