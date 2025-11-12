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
/* 
   Minimum acceptable nugget for Gaussian variograms
*/
#define NUGGET_MIN 1e-6

void setUpcorr2D(struct VARIOGRAM_2D *variogram);
double corr2D(double x1,double y1,double x2,double y2,struct VARIOGRAM_2D *variogram);
double corr2Dexp(double x1,double y1,double x2,double y2,struct VARIOGRAM_2D *variogram);
void setUpcorr3D(struct VARIOGRAM_3D *variogram);
double corr3D(double x1,double y1,double z1,double x2,double y2,double z2,struct VARIOGRAM_3D *variogram);

