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



package au.csiro.JamesGunning.waveletextraction;

import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.lib.util.*;


/**
   Data class for storing well trajectory information.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class UTM
{
  /**
	 number of points in survey
  */
  int N;
  /**
	 measured depth
  */
  double[] Zm;
  /**
	 true vertical depth
  */
  double[] Z;
  /**
	 X
  */
  double[] X;
  /**
	 Y
  */
  double[] Y;

  /**
	 First and last elements of Z[]
  */
  double zMax,zMin;

  /** 
	  Useful spline for Z->X interpolation
  */
  Spline xSpline;

  /** 
	  Useful spline for Z->Y interpolation
  */
  Spline ySpline;

  /**
	 flag we set to speedup interpolations etc
  */
  boolean wellIsVertical=false;
  /**
	 Call this to set up useful natural splines xSpline, ySpline once the UTM data are in place
  */
  void setUpXYSplines(ErrorHandler EH, CheckshotAndMarkerData chkShotandMarkerData)
  {
	int i;
	double x,y,z,zm;

	// check survey XY interpolation possible over this measured depth range
	double ZmMin = (chkShotandMarkerData.Zchm)[0];
	double ZmMax = chkShotandMarkerData.Zchm[chkShotandMarkerData.Nch-1];

	if (ZmMin<Zm[0])
	  {
		EH.println("        WARNING: Survey data does not start high enough for depth %10.6g",ZmMin);
		EH.println("        Adding extra survey point by linear extrapolation: Beware!");

		x = ((X[1]-X[0])/(Zm[1]-Zm[0])) * (ZmMin - Zm[0]) + X[0];
		y = ((Y[1]-Y[0])/(Zm[1]-Zm[0])) * (ZmMin - Zm[0]) + Y[0];
		z = ((Z[1]-Z[0])/(Zm[1]-Zm[0])) * (ZmMin - Zm[0]) + Z[0];

		EH.println("        New survey point [X Y Z Zm]: %10.6g %10.6g %10.6g %10.6g",
				   x,y,z,ZmMin);
		N++;
		X = ArrayUtils.prepend(x,X);
		Y = ArrayUtils.prepend(y,Y);
		Z = ArrayUtils.prepend(z,Z);
		Zm = ArrayUtils.prepend(ZmMin,Zm);
	  }

	if (ZmMax>Zm[N-1])
	  {
		EH.println("        WARNING: Survey data does not go deep enough for depth %10.6g",ZmMax);
		EH.println("        Adding extra survey point by linear extrapolation: Beware!");

		x = ((X[N-2]-X[N-1])/(Zm[N-2]-Zm[N-1])) * (ZmMax - Zm[N-1]) + X[N-1];
		y = ((Y[N-2]-Y[N-1])/(Zm[N-2]-Zm[N-1])) * (ZmMax - Zm[N-1]) + Y[N-1];
		z = ((Z[N-2]-Z[N-1])/(Zm[N-2]-Zm[N-1])) * (ZmMax - Zm[N-1]) + Z[N-1];

		EH.println("        New survey point point [X Y Z Zm]: %10.6g %10.6g %10.6g %10.6g",
				   x,y,z,ZmMax);
		N++;
		X = ArrayUtils.append(x,X);
		Y = ArrayUtils.append(y,Y);
		Z = ArrayUtils.append(z,Z);
		Zm = ArrayUtils.append(ZmMax,Zm);
	  }

	xSpline = new Spline(Z, X, N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
	ySpline = new Spline(Z, Y, N, Double.MAX_VALUE, Double.MAX_VALUE,EH);

	
	zMin = Z[0];
	zMax = Z[N-1];
  }

  /**
	 Simple checks for monotonicity of survey etc
  */
  void checkSurveyDataForSanity(ErrorHandler EH)
  {
	int i;
	for (i=1;i<N;i++)
	  {
		if (Zm[i] <= Zm[i-1])
		  {
			EH.die_err("         WARNING: Survey data has out-of-order or repeated point! Zm="+Zm[i]);
		  }
	  }
	// setup verticality flag
	wellIsVertical=true;
	for (i=1;i<X.length;i++)
	  {
		if ((X[i] != X[0]) || (Y[i] != Y[0])) {wellIsVertical=false; break;}
	  }
	if (wellIsVertical) EH.println("                         Well is vertical.");
  }
}
