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
import au.csiro.JamesGunning.lib.Numerical.*;
/**
   Class for performaing Z<->Zm interpolation based on well survey
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class WellSurveyDepthConverter
{
  private ErrorHandler EH;
  private UTM utmData;
  private Spline zSpline;

  /**
	 Constructor. Set up cubic splines from the UTM data for each coordinate
  */
  public WellSurveyDepthConverter(UTM utmData,ErrorHandler EH)
  {
	this.EH = EH;
	this.utmData = utmData;
	initialiseDepthConverter();
  }

  private void initialiseDepthConverter()
  {
	if (true)
	  {
		// set up a simple spline to get Z values for any Zm, so we can map the upscaled logs
		// to true depth
		zSpline = new Spline(utmData.Zm, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
	  }
	else
	  {
		// here's a chunk of code to get Zm values if the survey has only X,Y,Z
		int i,j,nSeg;
		double t,dt,x0,x1,y0,y1,z0,z1,length,old_length;
		double[] xParametric, Zm;
		Spline xSpline;
		Spline ySpline;
		Spline zSpline;
		
		double EPS=1.0e-7; // FIXME set sane somewhere
		// use i=0,1,... as fake x-data
		xParametric = new double[utmData.N];
		for (i=0; i< utmData.N ; i++) xParametric[i]=(double)i;

		// create parametric splines
		xSpline = new Spline(xParametric, utmData.X, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
		ySpline = new Spline(xParametric, utmData.Y, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
		zSpline = new Spline(xParametric, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
		// now compute sequence of arc-lengths vs true depth pairs using the
		// formulae for arc lengths of cubic
		length=0;
		Zm = new double[utmData.N];
		for(i=0; i< utmData.N-1; i++)
		  {
			old_length=Double.MAX_VALUE;
			// get segment length by extrapolating sum of segments as discretisation bisects
			for (nSeg=1,dt=1;;dt/=2.0,nSeg*=2)
			  {
				length=0.0;
				x0=xSpline.splint(xParametric[i],true);
				y0=ySpline.splint(xParametric[i],true);
				z0=zSpline.splint(xParametric[i],true);
				t=xParametric[i]+dt;
				// sum of straight line segments
				for(j=1;j<=nSeg;j++,t+=dt)
				  {
					x1=xSpline.splint(t,true);
					y1=ySpline.splint(t,true);
					z1=zSpline.splint(t,true);
					length += Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
					x0=x1;
					y0=y1;
					z0=z1;
				  }
				if (Math.abs(length-old_length)<EPS*length) break;
				old_length=length;
			  }
			// ok the length is now right
			Zm[i]=length + ((i>0) ? Zm[i-1] : 0);
			//System.out.println("length = "+Zm[i]);
		  }
		// and we need to set up this naturalspline.......
		zSpline = new Spline(Zm, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);
	  }
  }

  /**
	 Depth converter. 
	 Given an array of measured depths contained within the survet range, return an array of the
	 equivalent true vertical depths, interpolated out of the spline zSpline built by the
	 Constructor
	 @param measuredDepths Array of depths which *must* fall in the range Zm[0]...Zm[N-1]
	 of the UTM object fed to the constructor
  */
  public double[] getTrueVerticalDepths(double[] measuredDepths)
  {
	double[] trueDepths = measuredDepths.clone();
	for (int i=0; i<trueDepths.length ; i++)
	  {	  
		trueDepths[i] = zSpline.splint(measuredDepths[i],true);
		//System.out.println(	trueDepths[i]+" "+measuredDepths[i]);
	  }
	return trueDepths;
  }
}
