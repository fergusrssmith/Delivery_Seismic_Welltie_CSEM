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

import au.csiro.JamesGunning.lib.ErrorHandler;
/** 
	Useful info for making UTM conversions from seismic
*/
public class UTMConversionInfo
{
  double[] x = new double[3];
  double[] y = new double[3];
  double[] inLine = new double[3];
  double[] xLine = new double[3];
  double[][] fwdMatrix = new double[2][2];
  double[][] revMatrix = new double[2][2];
  double[] mapOffset = new double[2];

  UTMxyValues pt = new UTMxyValues();

  boolean setUpUTMmapping(ErrorHandler EH)
  {
	int i,j,k;
	double[][] DR = new double[2][2];
	double[][] DS = new double[2][2];
	double det;
	UTMxyValues pt;

	/* form DR = (r0-r1 | r0-r2) */
	DR[0][0] = x[0] - x[1];
	DR[0][1] = x[0] - x[2];
	DR[1][0] = y[0] - y[1];
	DR[1][1] = y[0] - y[2];

	/* form DS = (s0-s1 | s0-s2) */
	DS[0][0] = inLine[0] - inLine[1];
	DS[0][1] = inLine[0] - inLine[2];
	DS[1][0] = xLine[0] - xLine[1];
	DS[1][1] = xLine[0] - xLine[2];

	det = DS[0][0] * DS[1][1] - DS[1][0] * DS[0][1];
	if (Math.abs(det) < 1e-6)
	  {
		EH.warn("UTM coordinates in XML file are not unique, or collinear");
		return false;
	  }

	/* form fwdMatrix = DR.DS^{-1} */
	fwdMatrix[0][0] = (-DR[0][1]*DS[1][0] + DR[0][0]*DS[1][1])/det;
	fwdMatrix[0][1] = (DR[0][1]*DS[0][0] - DR[0][0]*DS[0][1])/det;
	fwdMatrix[1][0] = (-DR[1][1]*DS[1][0] + DR[1][0]*DS[1][1])/det;
	fwdMatrix[1][1] =  (DR[1][1]*DS[0][0] - DR[1][0]*DS[0][1])/det;
	
	/* set up offset */
	mapOffset[0] = 0.0;
	mapOffset[1] = 0.0;
	pt = getUTMSfromInlineXline(inLine[0],xLine[0]);
	mapOffset[0] = x[0] - pt.x;
	mapOffset[1] = y[0] - pt.y;
	
	// reverse mapping 
	det = fwdMatrix[0][0] * fwdMatrix[1][1] - fwdMatrix[0][1] * fwdMatrix[1][0];
	revMatrix[0][0] = fwdMatrix[1][1]/det;
	revMatrix[0][1] = -fwdMatrix[0][1]/det;
	revMatrix[1][0] = -fwdMatrix[1][0]/det;
	revMatrix[1][1] = fwdMatrix[0][0]/det;


	// checks
	double eps = 1e-10;
	for (i=0;i<3;i++)
	  {
		pt = getUTMSfromInlineXline(inLine[i],xLine[i]);
		// System.out.println("UTM transformation check "+(pt.x - x[i])+" "+(pt.y - y[i]));
		if ((Math.abs(pt.x - x[i]) > 1e-10) || (Math.abs(pt.y - y[i])> 1e-10))
		  {
			EH.warn("UTM fwd transformation is buggy, point "+i);
		  }
		pt = getInlineXlineFromUTMS(x[i],y[i]);
		// ySystem.out.println("UTM transformation check "+(pt.x - inLine[i])+" "+(pt.y - xLine[i]));
		if ((Math.abs(pt.x - inLine[i]) > 1e-10) || (Math.abs(pt.y - xLine[i])> 1e-10))
		  {
			EH.warn("UTM reverse transformation is buggy, point "+i);
		  }
	  }
	return true;
  }

  /**
	 Compute the x,y value of these inline, xline coordinates and return as a
	 UTMxyValues object
  */
  UTMxyValues getUTMSfromInlineXline(double inLine,double xLine)
  {
	pt.x = fwdMatrix[0][0] * inLine + fwdMatrix[0][1] * xLine + mapOffset[0];
	pt.y = fwdMatrix[1][0] * inLine + fwdMatrix[1][1] * xLine + mapOffset[1];
	
	return pt;
  }

  /**
	 Compute the (full-precision) inline, xline coordinates from x,y and return as a
	 UTMxyValues object
  */
  UTMxyValues getInlineXlineFromUTMS(double x,double y)
  {
	double dx,dy;

	dx = x - mapOffset[0];
	dy = y - mapOffset[1];
	pt.x = revMatrix[0][0] * dx + revMatrix[0][1] * dy;
	pt.y = revMatrix[1][0] * dx + revMatrix[1][1] * dy;
	
	// System.out.println("inline = " +pt.x+ " xline = "+pt.y);
	return pt;
  }

  public class UTMxyValues
  {
	double x;
	double y;
  }
}
