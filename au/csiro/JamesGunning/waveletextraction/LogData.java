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
/**
   Data class for storing general log information. Not a very general
   design, but should be efficient.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

import cern.colt.list.*;


import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;


public class LogData
{
  /**
	 measured depth
  */
  double[] Zm;
  // unit to be determined from LAS file
  WellLogDefinitions.MD_unit Zm_unit = WellLogDefinitions.MD_unit.UNKNOWN; 

  /**
	 the raw log  DT (p-wave velocity) data
  */
  double[] DT;

  // unit to be determined from LAS file
  WellLogDefinitions.DT_unit DT_unit = WellLogDefinitions.DT_unit.UNKNOWN;

  // Density. Unit not used - just store string for now
  String RHO_unit="?";


  /**
	 the raw log DT (shear-wave velocity) data
  */
  double[] DTshear=null;
  /**
	 Density data
  */
  double[] rho;
  /**
	 Optional gamma log
  */
  double[] gamma;
  /**
	 Greatest useful depth of log data (determined after reading)
  */
  double ZmMax;
  /**
	 Least useful depth of log data (determined after reading)
  */
  double ZmMin;

  boolean haveGamma = true;
  
  /**
	 Is true if there exists data just above longest good interval found by a call to
	 setupBracketOfInterest().
  */
  boolean dataExistsBelowLongestGoodInterval;
  /**
	 Is false if there exists data just below longest good interval found by a call to
	 setupBracketOfInterest().
  */
  boolean dataExistsAboveLongestGoodInterval;
  /**
	 index corresponding to ZmMin
  */
  int iStartLongestGoodInterval;
  /**
	 index corresponding to ZmMax
  */
  int iEndLongestGoodInterval;
  /**
	 Greatest useful time of log data (from convenient time-to-depth map)
  */
  double tMax;
  /**
	 Least useful time of log data (from convenient time-to-depth map)
  */
  double tMin;
  /** 
	  Survey data. Set of X,Y,TVD,MD points
  */
  UTM survey;

  /**
	 Minimum amount of consecutive decent log data when testing
  */
  public static final int MIN_CONSECUTIVE_DECENT_LOG_DATA_COUNT=100;
  
  private double zStartGoodData, zEndGoodData;

  private ErrorHandler _EH;
  /**
	 Constructor. Create room for survey by default
  */
  public LogData(ErrorHandler EH)
  {
	survey = new UTM();
	this._EH = EH;
  }
  /**
	 Method to clean up log data and restrict it to named interval.
	 Checks for bad data in this interval as well.
	 Assumes ZmMin < ZmMax
  */
  void cleanUpDataAndCheckIntegrity(double ZmMin, double ZmMax)
  {
	int i, iStart, iEnd;

	iStart = Math.max(Model.getStraddlingLowerIndex(ZmMin,Zm)+1,0);
	iEnd = Math.min(Model.getStraddlingLowerIndex(ZmMax,Zm),Zm.length);
	for (i=iStart;i<iEnd;i++)
	  {
		if (DT[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad DT data at depth " + Zm[i] + ". Can't proceed");
		  }
		if (DTshear[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad shear data at depth " + Zm[i] + ". Can't proceed");
		  }
		if (rho[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad density data at depth " + Zm[i] + ". Can't proceed");
		  }
	  }

	// If we got to here, things are OK: reduce data.
	Zm = ArrayUtils.subselect(iStart,iEnd,Zm);
	DT = ArrayUtils.subselect(iStart,iEnd,DT);
	DTshear = ArrayUtils.subselect(iStart,iEnd,DTshear);
	rho = ArrayUtils.subselect(iStart,iEnd,rho);
	if (gamma != null && gamma.length != 0) 
	  {
		gamma = ArrayUtils.subselect(iStart,iEnd,gamma);
	  }
  }

  /**
	 Method to clean up log data and restrict it to interval
	 spanning named point.
	 Checks for bad data in this interval as well.
	 Assumes Zm_ref
  */
  void cleanUpDataPlaceholder(double Zm_ref)
  {
	int i, iStart, iEnd;


	// _EH.println("FIXME 2");
	// FIXME 2 restore these to below...
	// boolean allowInterpolation,boolean eliminateSpikes,double maxDepthInterpolation)

	setupBracketOfInterest(Zm[0],Zm[Zm.length-1],false,false,0);

	iStart = Math.max(Model.getStraddlingLowerIndex(ZmMin,Zm)+1,0);
	iEnd = Math.min(Model.getStraddlingLowerIndex(ZmMax,Zm),Zm.length);
	for (i=iStart;i<iEnd;i++)
	  {
		if (DT[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad DT data at depth " + Zm[i] + ". Can't proceed");
		  }
		if (DTshear[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad shear data at depth " + Zm[i] + ". Can't proceed");
		  }
		if (rho[i]<FundamentalConstants.BAD_NUMBER)
		  {
			throw new RuntimeException("Bad density data at depth " + Zm[i] + ". Can't proceed");
		  }
	  }

	// If we got to here, things are OK: reduce data.
	Zm = ArrayUtils.subselect(iStart,iEnd,Zm);
	DT = ArrayUtils.subselect(iStart,iEnd,DT);
	DTshear = ArrayUtils.subselect(iStart,iEnd,DTshear);
	rho = ArrayUtils.subselect(iStart,iEnd,rho);
	if (gamma != null && gamma.length != 0) 
	  {
		gamma = ArrayUtils.subselect(iStart,iEnd,gamma);
	  }
  }

  /**
	 Finds a bracket of the longest continuous 
	 logged interval of good data in the provided depth interval. The
	 fields ZmMax, ZmMin are  set to correspond to this as well. 
	 If allowInterpolation==true, holes in the log data of length < maxDepthInterpolation
	 will be filled by interpolation.

	 Returns false if nothing useful 
	 is found.
  */
  boolean setupBracketOfInterest(double zLow,double zHigh,
								 boolean allowInterpolation,boolean eliminateSpikes,
								 double maxDepthInterpolation)
  {
	int i,j,iEnd,iStart=0,iEndInt=0,length,bestLength=0;
	double ZmEnd,ZmStart;

	// initialise
	iStartLongestGoodInterval = 0;
	dataExistsAboveLongestGoodInterval = false;
	dataExistsBelowLongestGoodInterval = false;

	i = Math.min(Math.max(Model.getStraddlingLowerIndex(zLow,Zm)+1,0),Zm.length-1);
	iEnd = Math.min(Math.max(Model.getStraddlingLowerIndex(zHigh,Zm),0),Zm.length-1);

	// loop to get biggest intact data sequence
	boolean isNotFirstInterval = false;
	ZmEnd=Zm[i];
	iEndInt = i;
	// Define threshold of usable data. Here, things must be +ve, not just >-999
	// could beFundamentalConstants.BAD_NUMBER (usually -999)
	for (;i<=iEnd;i++)
	  {
		// search for start of interval
		for (;i<=iEnd;i++)
		  {
			if (DT[i]>FundamentalConstants.POS_THRESHOLD && DTshear[i]>FundamentalConstants.POS_THRESHOLD
				&& rho[i] >FundamentalConstants.POS_THRESHOLD)
			  {
				// good data point. We are at the start of an interval
				ZmStart = Zm[i];
				if (allowInterpolation && isNotFirstInterval
					&& ((ZmStart - ZmEnd) <= maxDepthInterpolation))
				  {
					// test to see if we patch up the break just encountered
					// form linear interpolations of three logs over missing data interval
					_EH.println("         WARNING: interpolating bad log data at depth "
								+Zm[iEndInt+1]+" to "+Zm[i-1]);

					double rhoSlope,rhoIntercept;
					double DTSlope,DTIntercept;
					double DTshearSlope,DTshearIntercept;
					double denom = ((i-iEndInt) != 0) ? (i-iEndInt) : 1.0;
					DTSlope = (DT[i]-DT[iEndInt])/denom;
					DTIntercept = DT[i] - DTSlope*i;
					DTshearSlope = (DTshear[i]-DTshear[iEndInt])/denom;
					DTshearIntercept = DTshear[i] - DTshearSlope*i;
					rhoSlope = (rho[i]-rho[iEndInt])/denom;
					rhoIntercept = rho[i] - rhoSlope*i;
					// now hard-overwrite logs in this gap with interpolations
					for (j=iEndInt+1; j<i; j++)
					  {
						DT[j] = DTIntercept + DTSlope*j;
						DTshear[j] = DTshearIntercept + DTshearSlope*j;
						rho[j] = rhoIntercept + rhoSlope*j;
					  }
					// don't update iStart
				  }
				else
				  {
					// we have a genuine new start
					iStart=i;
				  }
				isNotFirstInterval=true;
				break;
			  }
		  }

		// search for end of interval
		for (;i<=iEnd;i++)
		  {
			if ((DT[i]<=FundamentalConstants.POS_THRESHOLD) || 
				(DTshear[i]<=FundamentalConstants.POS_THRESHOLD)
				|| (rho[i]<=FundamentalConstants.POS_THRESHOLD)
				|| (i==iEnd))
			  {
				if (i<iEnd)
				  {
					iEndInt = i-1;
				  }
				else
				  {
					iEndInt = i;
				  }
				length = iEndInt-iStart;
				if (length > bestLength)
				  {
					bestLength = length;
					iStartLongestGoodInterval = iStart;
				  }
				ZmEnd = Zm[i];
				break;
			  }
		  }
	  }
	if (bestLength==0) 
	  {
		return false;
	  }
	else
	  {
		// set commensurate minimum depth
		_EH.println("         Longest log-suite interval: Zm = "+Zm[iStartLongestGoodInterval]+" ending at Zm = "+Zm[iEndInt-1]);
		ZmMin = Zm[iStartLongestGoodInterval];
		iEndLongestGoodInterval = iStartLongestGoodInterval + bestLength;
		ZmMax = Zm[iEndLongestGoodInterval];

		if (iStartLongestGoodInterval>0 && Zm[iStartLongestGoodInterval-1]<=zLow)
		  {
			dataExistsAboveLongestGoodInterval = true;
		  }
		if (iEndLongestGoodInterval<(Zm.length-1) && Zm[iEndLongestGoodInterval+1]>=zHigh)
		  {
			dataExistsBelowLongestGoodInterval = true;
		  }
		// despike if required
		if (eliminateSpikes) eliminateLogSpikes();

 		return true;
	  }
  }

  void eliminateLogSpikes()
  {
	int i,j,iStart=-1,iEnd=-1;
	int NSMOOTH=30;
	double DT_s,s;
	double[] wts= new double[2*NSMOOTH+1];
	double[] DT_cpy;
	double EPS=0.15;

	DT_cpy = DT.clone();

	// smoothing weights
	for (j=-NSMOOTH,s=0; j<=NSMOOTH; j++)
	  {
		wts[j+NSMOOTH]= 1.0-Math.abs(j/((double)NSMOOTH));
		s += wts[j+NSMOOTH];
	  }
	for (j=-NSMOOTH; j<=NSMOOTH; j++)
	  {
		wts[j+NSMOOTH] /=s;
	  }
	for (i=iStartLongestGoodInterval; i<=iEndLongestGoodInterval; i++)
	  {
		iStart=iEnd=-1;
		for (DT_s=0.0,s=0.0,j=Math.max(iStartLongestGoodInterval, i-NSMOOTH); 
			 j<=Math.min(i+NSMOOTH, iEndLongestGoodInterval); j++)
		  {
			DT_s += wts[(i-j)+NSMOOTH] * DT_cpy[j];
			s +=  wts[(i-j)+NSMOOTH];
		  }
		DT_s /=s;
		if (Math.abs(DT_cpy[i] - DT_s) > EPS*DT_s)
		  {
			iStart = (i>iStartLongestGoodInterval)? i-1 : i;
			//DT[i] = DT_s;
			//System.out.println("SMOOTH HERE "+Zm[i]);
		  }
		if (iStart>=0)
		  {
			// search fwd to end of transgression
			i++;
			for (; i<=iEndLongestGoodInterval; i++)
			  {
				for (DT_s=0.0,s=0.0,j=Math.max(iStartLongestGoodInterval, i-NSMOOTH); 
					 j<=Math.min(i+NSMOOTH, iEndLongestGoodInterval); j++)
				  {
					DT_s += wts[(i-j)+NSMOOTH] * DT_cpy[j];
					s +=  wts[(i-j)+NSMOOTH];
				  }
				DT_s /=s;
				if (Math.abs(DT_cpy[i] - DT_s) < EPS*DT_s)
				  {
					iEnd = (i<iEndLongestGoodInterval)? i+1 : i;
					break;
				  }
			  }
		  }
		if (iStart>=0 && iEnd>=0)
		  {
			// found a finite interval with transgressions
			for (i=iStart+1; i<iEnd; i++)
			  {
				DT[i] = (i-iStart)*(DT_cpy[iEnd]-DT_cpy[iStart])/(iEnd-iStart) + DT_cpy[iStart];
				DTshear[i] = (i-iStart)*(DTshear[iEnd]-DTshear[iStart])/(iEnd-iStart) + DTshear[iStart];
				rho[i] = (i-iStart)*(rho[iEnd]-rho[iStart])/(iEnd-iStart) + rho[iStart];
			  }
		  }
	  }
	if (false)
	  {
		for (i=iStartLongestGoodInterval; i<=iEndLongestGoodInterval; i++)
		  {
			System.out.println("LOG: "+Zm[i]+" "+DT_cpy[i]+" "+DT[i]);
		  }
		//System.exit(0);
	  }
	// 
  }

  /**
	 Returns first index of longest-good-data interval just found in
	 setupBracketOfInterest().
	 Use/call immediately after setupBracketOfInterest()
	 for number to be sensible.
  */
  int getStartIndexOfLongestDataIntervalJustFound()
  {
	return iStartLongestGoodInterval;
  }

  /**
	 Returns last good index of longest-good-data interval just found in
	 setupBracketOfInterest().
	 Use/call immediately after setupBracketOfInterest()
	 for number to be sensible.
  */
  int getLargestIndexOfLongestDataIntervalJustFound()
  {
	return iEndLongestGoodInterval;
  }


  /**
	 Returns largest index of data above seedDepth point (in Zm) which brackets a continuous 
	 logged interval of good data below  the seedDepth. The
	 field ZmMax is set to correspond to this as well. 
	 Returns -1 for a bad start (zSeed has bad log data)
  */
  int setBiggestIndexOfGoodDataAfterThisPoint(double zSeed)
  {
	int i;

	i = Math.min(Math.max(Model.getStraddlingLowerIndex(zSeed,Zm)+1,0),Zm.length-1);
	// starting point check
	if (DT[i]<FundamentalConstants.POS_THRESHOLD || DTshear[i]<FundamentalConstants.POS_THRESHOLD
			|| rho[i]<FundamentalConstants.POS_THRESHOLD)
	  {
		// bad start
		return -1;
	  }
	else
	  {
		i++;
		for (;i<Zm.length;i++)
		  {
			if (DT[i]<FundamentalConstants.POS_THRESHOLD || DTshear[i]<FundamentalConstants.POS_THRESHOLD
				|| rho[i]<FundamentalConstants.POS_THRESHOLD)
			  {
				break;
			  }
		  }
		i = Math.max(0,i-1);
	  }


	// set commensurate greatest depth
	ZmMax = Zm[i];

	return i;
  }


  /**
	 Returns suggested data inteval of complete data after 
	 seedDepth point (in Zm) which brackets a continuous 
	 logged interval of good data below  the seedDepth.
  */
  boolean searchForGoodDataIntervalAfterThisPoint(double zSeed)
  {
	int i;
	int iStart,iEnd,cnt=0;

	i = Math.min(Math.max(Model.getStraddlingLowerIndex(zSeed,Zm)+1,0),Zm.length-1);
	if (i == Zm.length-1) return false;

	for (;i<Zm.length;i++)
	  {
		if (!(DT[i]<FundamentalConstants.POS_THRESHOLD || DTshear[i]<FundamentalConstants.POS_THRESHOLD
			  || rho[i]<FundamentalConstants.POS_THRESHOLD))
		  {
			break;
		  }
	  }
	iStart=i;
	zStartGoodData = Zm[i];
	// we are at a point with good data. Start counting
	for (;i<Zm.length;i++)
	  {
		cnt++;
		if ((DT[i]<FundamentalConstants.POS_THRESHOLD || DTshear[i]<FundamentalConstants.POS_THRESHOLD
			  || rho[i]<FundamentalConstants.POS_THRESHOLD))
		  {
			break;
		  }
	  }
	iEnd = i-1;
	zEndGoodData = Zm[i-1];

	if (cnt > MIN_CONSECUTIVE_DECENT_LOG_DATA_COUNT) 
	  {
		return true;
	  }
	else
	  {
		return false;
	  }

	// Start of interval
	// System.out.println("Start of log-suite interval at Zm = "+Zm[iStart]+" ending at Zm = "+Zm[iEnd-1]);
  }

  /**
	 After searchForGoodDataIntervalAfterThisPoint() is called, this will return the
	 Zm of the start of the interval of good log data
  */
  double getGoodDataIntervalAfterThisPointZmin()
  {
	return zStartGoodData;
  }
  /**
	 After searchForGoodDataIntervalAfterThisPoint() is called, this will return the
	 Zm of the end of the interval of good log data
  */
  double getGoodDataIntervalAfterThisPointZmax()
  {
	return zEndGoodData;
  }

  double getSonicLogScaling()
  {
	double DT_scaling;

	if (DT_unit == WellLogDefinitions.DT_unit.S_PER_M)	  
	  DT_scaling = 1.0;	  
	else
	  DT_scaling = FundamentalConstants.SONIC_LOG_TIME_SCALING;

	return DT_scaling;
  }
 
  String getDTandRhoUnits()
  { 
	String s;

	if (DT_unit == WellLogDefinitions.DT_unit.US_PER_FT)
	  {
		s="us/ft";
	  }
	else if (DT_unit == WellLogDefinitions.DT_unit.US_PER_M)
	  {
		s="us/m";		
	  }
	else if (DT_unit == WellLogDefinitions.DT_unit.S_PER_M)
	  {
		s="s/m";
	  }
	else 
	  {
		s="?";
	  }
	// String for unit_DT unit_DTS unit_rho
	return s+" "+s+" "+RHO_unit;
  }

  /**
	 Sort log data into increasing order if reverse order detected (common in LAS files)
  */
  void reverseDepthOrderOfData()
  {
	int i,j,size;
	double tmp;
	size = Zm.length-1;

	for(i=0,j=size; i<j; i++,j--) 
	  {
		tmp = Zm[i];
		Zm[i] = Zm[j];
		Zm[j] = tmp;
	  }

	for(i=0,j=size; i<j; i++,j--) 
	  {
		tmp = DT[i];
		DT[i] = DT[j];
		DT[j] = tmp;
	  }

	if (DTshear!=null)
	  {
		for(i=0,j=DTshear.length-1; i<j; i++,j--) 
		  {
			tmp = DTshear[i];
			DTshear[i] = DTshear[j];
			DTshear[j] = tmp;
		  }
	  }

	if (gamma!=null)
	  {
		for(i=0,j=gamma.length-1; i<j; i++,j--) 
		  {
			tmp = gamma[i];
			gamma[i] = gamma[j];
			gamma[j] = tmp;
		  }
	  }

	for(i=0,j=size; i<j; i++,j--) 
	  {
		tmp = rho[i];
		rho[i] = rho[j];
		rho[j] = tmp;
	  }
  }

  /**
	 Sort log data into increasing order if reverse order detected (common in LAS files)
  */
  void muteData(double MDu,double MDl)
  {
	int i,n;
	n = Zm.length;
	double tmp;
	for(i=0;i<n;i++)
	  {
		tmp = Zm[i];
		if((Zm[i] < MDu) || (Zm[i] > MDl)) DT[i] = FundamentalConstants.BAD_NUMBER;
	  }
  }
}
