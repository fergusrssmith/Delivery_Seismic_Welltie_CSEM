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

import java.util.*;
import java.io.*;

import cern.colt.list.*;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.stat.segmentation.*;

/**
   Data class for storing checkshot and marker information
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class CheckshotAndMarkerData
{
  // ptr to underyling extractor
  private WaveletExtractor _WE;

  private ErrorHandler EH;

  // well index in _WE.
  int well_indx;

  /**
	 Set to true if checkshot data has MD
	 depending how this is set, various conversions are applied in the right places etc
  */
  boolean USE_MD_IN_CHECKSHOTS = true;
  /**
	 Set to true if marker data has MD
	 depending how this is set, various conversions are applied in the right places etc
  */
  boolean USE_MD_IN_MARKERS = true;

  /**
	 number of points in checkshot arrays
  */
  int Nch=0;
  /**
	 Checkshot measured-depth
  */
  double[] Zchm;
  /**
	measured checkshot one-way time
  */
  double[] Tch;
  /**
	 Estimated (prior) std deviations in checkshot times
  */
  double[] sigmaTch;

 /**
	 number of points in marker arrays
  */
  int Nmk=0;
  /**
	 measured marker depth
  */
  double[] Zmkm;
  /**
	measured marker one-way time
  */
  double[] Tmk;
  /**
	 Estimated (prior) std deviations in marker times
  */
  double[] sigmaTmk;

  TimeToDepthMap.TimeToDepthTriples[] t2d;

  /**
	 Rough estimate of mean velocity in natural units
  */
  double meanNaturalVelocity = 0;
  
  /**
	 Detail for auto-generated checkshots if used
  */
  AutoCheckshotInfo autoCheckshotInfo = null;

  public CheckshotAndMarkerData(WaveletExtractor WE, int well_indx)
  {
	this._WE = WE;
	this.EH = WE.EH;
	this.well_indx = well_indx;
  }

  /**
	 Method to setUp marker and checkshot data as the same fundamental
	 kinds of constraints. call after the basic checkshots and markers 
	 information has been set up. This creates and initialises the arrays
	 Z_mean, T_mean, sigma_T.
  */
  void setUpCheckShotAndMarkerData(UTM utmData)
  {
	int i,j,klo,khi,k;
	double zm,wt,t;
	Spline zSpline;

	// complete description is union of markers and checkshots
	t2d = new TimeToDepthMap.TimeToDepthTriples[Nmk+Nch];

	// spline to convert MD's in checkshots/markers etc to TVD's
	zSpline = new Spline(utmData.Zm, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,EH);

	// convert measured checkshot depths to true depths
	for (i=0; i<Nch; i++)
	  {
		t2d[i] = new TimeToDepthMap.TimeToDepthTriples();
		t2d[i].Z_mean = zSpline.splint(Zchm[i],true);
		t2d[i].T_mean = Tch[i];
		t2d[i].sigma_T = sigmaTch[i];
	  }
	
	// loop over markers to convert them also
	for (i=0,j=Nch;i<Nmk;i++,j++)
	  {
		// depth
		zm = Zmkm[i];
		// find true time by adding deviation back to linearly interpolated
		// checkshot data

		// bounding index in checkshot array.
		if (zm <Zchm[0]) EH.die_err("Marker data of depth "+zm+" above checkshot range");
		if (zm > Zchm[Zchm.length-1]) EH.die_err("Marker data of depth "+zm+" below checkshot range");
		// bisection to find bracket
		klo = 0;
		khi = Zchm.length-1;
		while (khi-klo > 1) 
		  {
			k=(khi+klo) >> 1;
			if (Zchm[k] > zm) khi=k;
			else klo=k;
		  }
		// now Zchm[klo] <= zm < Zchm[khi]  && khi=klo+1
		wt = (Zchm[khi] - zm)/(Zchm[khi] -  Zchm[klo]);
		// linearly interpolated value
		t = wt * Tch[klo]+(1-wt)*Tch[khi];
		// true value after correction from marker
		t += Tmk[i];

		// Add these to the complete descriptions
		t2d[j]=new TimeToDepthMap.TimeToDepthTriples();
		// true depth 
		t2d[j].Z_mean = zSpline.splint(zm,false);
		t2d[j].T_mean = t;
		t2d[j].sigma_T = sigmaTmk[i];
	  }
	// now ensure triples are depth sorted
	Arrays.sort(t2d);	
  }

  /**
	 Useful method to check checkshot data is sane. Checks data ordering, and
	 anything else you might want to insert here
  */
  void checkCheckshotAndMarkerDataForSanity(LogData logData, GeneralWellInfo WI)
  {
	int i,iStart,iEnd, iCheck,indx;
	boolean timesInMilliseconds = false;
	double testTime, Z_new, Zmin,Zmax,z,T_new, sigmaT_new;
	double tStartCheck,tEndCheck;
	double tStartCheckOffset,tEndCheckOffset;
	int synchIndx;
	Property systematicTProp=null;

	// Checkshot

	// (a) Apply some simple checks
	if ((Nch < 2)) //  && !WI.checkshotAutoMake)
	  {
		EH.die_err("Must be at least two checkshot points for extraction to be possible");
	  }
	if (_WE.verbosity>5) EH.println("Combined checkshot and marker data:");
	if (_WE.TIME_UNITS_AUTO_DETECTION)
	  {
		// test if mean checkshot time looks like ms
		for (i=0,testTime=0; i<Nch; i++) testTime += Tch[i];
		testTime/= Nch;
		if (testTime > FundamentalConstants.TIME_UNITS_AUTO_DETECTION_TIME_MAX)
		  {
			timesInMilliseconds=true;
		  }
	  }
	for (i=0; i<Nch; i++)
	  {
		// times may be in ms. May need to be forced into seconds....
        if (_WE.TIME_UNITS_AUTO_DETECTION)
          {
			// if (Tch[i] > FundamentalConstants.TIME_UNITS_AUTO_DETECTION_TIME_MAX)
			if (timesInMilliseconds)
			  {
				Tch[i] *= 1.0e-3;
				sigmaTch[i] *= 1.0e-3;
			  }
		  }
		if (_WE.verbosity > 5) EH.println("%10.6f  %10.6f",Zchm[i],Tch[i]);
		if (i>0)
		  {
			if (Zchm[i] <= Zchm[i-1]) EH.die_err("Checkshot data out of depth order, row "+i);
			if (Zchm[i] == Zchm[i-1]) EH.die_err("Checkshot data has depth repetitions, row "+i);
			if (Tch[i] <= Tch[i-1]) EH.die_err("Checkshot data out of time order, row "+i);
			if (Tch[i] == Tch[i-1]) EH.die_err("Checkshot data has time repeats, row "+i);
			// order -of-magnitude for velocities ??
		  }
		if (Zchm[i] < 0) EH.die_err("Checkshot data has -ve depths, row "+i);
		if (Tch[i] < 0) EH.die_err("Checkshot data has -ve times, row "+i);
		if (sigmaTch[i] < 0) EH.die_err("Checkshot data has -ve time uncertainties, row "+i);
	  }

	// (b) time extent of checkshot
	testTime = Tch[0] + FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * sigmaTch[0];
	if (_WE.m.registrationErrorsAreCommon)
	  {
		synchIndx = _WE.wellInfo[well_indx].synchronizeWellIndex;
		systematicTProp = _WE.seismicInfo[synchIndx][0].systematicTProp;
	  }
	// If we allow use of logs right up to _WE.tStart, we may have to 
	// compute reflectivities above _WE.start, so be safe
	tStartCheck = (_WE.assumeNoReflectionsOutsideLogs)? (_WE.tStart - _WE.tCodaMax) : _WE.tStart;
	if ((testTime > tStartCheck)) //  && !WI.checkshotAutoMake)
	  {
		EH.println("         WARNING! Checkshot has no safe point before start of estimation interval ["+
				   _WE.tStart+","+_WE.tEnd+"]");
		EH.println("         Adding a safe new checkshot point by linear extrapolation: Beware!");

		// Add a new checkshot at 4 or so standard deviations above the interval start
		// with extra margin of a few dt.
		// Make it have zero std. deviation, so it can't wander into the test interval and wreak
		// havoc
		sigmaT_new = 0.0;
		T_new = tStartCheck - 2*_WE.dt-FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * sigmaT_new;
		Z_new = (Zchm[1]-Zchm[0])*(T_new-Tch[0])/(Tch[1]-Tch[0]) + Zchm[0];
		EH.println("         New checkshot point [MD, t, sigmaT]: %.6g %.4g %.4g",
				   Z_new, T_new, sigmaT_new);
		Nch++;
		Zchm = ArrayUtils.prepend(Z_new,Zchm);
		Tch = ArrayUtils.prepend(T_new,Tch);
		sigmaTch = ArrayUtils.prepend(sigmaT_new,sigmaTch);

		// save for use in autoSegmentModelselect
		if (autoCheckshotInfo != null)
		  {
			autoCheckshotInfo.safeUpperCheckshot_Tch = T_new;
			autoCheckshotInfo.safeUpperCheckshot_Zchm = Z_new;
		  }
	  }

	// end of interval
	testTime = Tch[Nch-1] - FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * sigmaTch[Nch-1];

	// If we allow use of logs right up to _WE.tEnd, we may 
	// have to compute reflectivities past _WE.tEnd,
	// so be safe
	tEndCheck = (_WE.assumeNoReflectionsOutsideLogs) ? (_WE.tEnd - _WE.tPrecursorMin) : _WE.tEnd;
	if ((testTime < tEndCheck)) //  && !WI.checkshotAutoMake)
	  {
		EH.println("         WARNING! Checkshot has no safe point after end of estimation interval ["+
				   _WE.tStart+","+_WE.tEnd+"]");
		EH.println("         Adding a safe new checkshot point by linear extrapolation: Beware!");

		// Add a new checkshot at 4 or so standard deviations after tEnd
		// Make it have zero std. deviation, so it can't wander into the test interval and wreak
		// havoc
		// OLD: sigmaT_new = sigmaTch[Nch-1];
		sigmaT_new = 0.0;
		T_new = tEndCheck +  2*_WE.dt + FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * sigmaT_new;
		Z_new = (Zchm[Nch-2]-Zchm[Nch-1])*(T_new-Tch[Nch-1])/(Tch[Nch-2]-Tch[Nch-1]) + Zchm[Nch-1];
		if (true)
		  {
			// extend to maximum of well depth ?
			if (logData.Zm[logData.Zm.length-1] > Z_new)
			  {
				EH.println("         New checkshot point extended to end of log data");
				Z_new = logData.Zm[logData.Zm.length-1];
			  }
			// commensurate linearly extrapolated point...
			T_new = ((Tch[Nch-2]-Tch[Nch-1])/(Zchm[Nch-2]-Zchm[Nch-1]))*(Z_new-Zchm[Nch-1]) + Tch[Nch-1];
		  }
		if (T_new < Tch[Nch-1])
		  {
			// suggested extrapolation point is out of order. Bleat.
			EH.println("         Suggested new checkshot point is [MD, t, sigmaT]: %.6g %.4g %.4g",
					   Z_new, T_new, sigmaT_new);
			EH.println("         But this falls within existing checkshot! Recommend ");
			EH.println("         you deepen last checkshot point and/or reduce its time error!");
			EH.die();
		  }

		EH.println("         New checkshot point [MD, t, sigmaT]: %.6g %.4g %.4g",
				   Z_new,T_new, sigmaT_new);
		Nch++;
		Zchm = ArrayUtils.append(Z_new,Zchm);
		Tch = ArrayUtils.append(T_new,Tch);
		sigmaTch = ArrayUtils.append(sigmaT_new,sigmaTch);

		// save for use in autoSegmentModelselect
		if (autoCheckshotInfo != null)
		  {
			autoCheckshotInfo.safeLowerCheckshot_Tch = T_new;
			autoCheckshotInfo.safeLowerCheckshot_Zchm = Z_new;
		  }
	  }

	// (c) Construct estimates of range of log data necessary under the 
	// usual stochastic model for checkshots
	tStartCheckOffset = tStartCheck;
	tEndCheckOffset = tEndCheck;

	// (1) Smallest Zm
	for (iStart=0; iStart<Nch; iStart++)
	  {
		if (Tch[iStart] > tStartCheckOffset) 
		  {
			break;
		  }
	  }
	iStart--;
	// iStart now indexes the checkshot just prior to the estimation interval [tStartCheckOffset,tEndCheckOffset]
	// Here's an estimate of the starting depth which logs ought to be present for, using
	// linear interpolation over the spanning interval
	Zmin = ((Zchm[iStart+1]-Zchm[iStart])/(Tch[iStart+1]-Tch[iStart]))*(tStartCheckOffset - Tch[iStart])+ Zchm[iStart];

	// (2) Largest Zm
	for (iEnd=Nch-1; iEnd>0; iEnd--)
	  {
		if (Tch[iEnd] < tEndCheckOffset) 
		  {
			break;
		  }
	  }
	// iEnd now indexes the checkshot just prior to end of the estimation interval [tStart,tEnd]
	// Here's an estimate of the final depth which logs ought to be present for, using
	// linear interpolation over the spanning interval
	Zmax  = ((Zchm[iEnd+1]-Zchm[iEnd])/(Tch[iEnd+1]-Tch[iEnd]))*(tEndCheckOffset - Tch[iEnd])+ Zchm[iEnd];
	EH.println("         Logs need to span from measured depth %-10.6g to %-10.6g",Zmin,Zmax);


	// Suggested depths from span of useful data in logs
	boolean foundLogData = logData.setupBracketOfInterest(Zmin,Zmax,_WE.allowLogInterpolation,
														  _WE.eliminateSpikes,
														  _WE.logInterpolationMaxDist);

	if (!foundLogData)
	  {
		EH.println("         WARNING! Logs not useful in estimation time interval [%-8.4g, %-8.4g]%n"
				  +"                  (estimated depths =  [%-8.4f %-8.4f])",
				   _WE.tStart, _WE.tEnd, Zmin, Zmax);
		EH.println("                  Try another interval!");
		EH.die_err("");
	  }
	
	// now get time corresponding to start of good log data interval
	logData.tMin = computeApproximateTimeDepthPair(logData.ZmMin);
	EH.println("         Useful log data interval found starting at %-10.6g, which has time approx %-10.6g",
			   logData.ZmMin, logData.tMin);

	if (!logData.dataExistsAboveLongestGoodInterval)
	  {
		EH.println("         WARNING! Logs not complete before start of needed interval [%-9.6f,%-9.6f].",
				   tStartCheck, tEndCheck);
		EH.println("                  Recommend a narrower estimation interval!");
		if (!_WE.ALLOW_INCOMPLETE_LOGS)
		  {
			EH.die_err("");
		  }
	  }
	//	iCheck = logData.getLargestIndexOfLongestDataIntervalJustFound();


	// compute corresponding approximate time from checkshot...
	logData.tMax = computeApproximateTimeDepthPair(logData.ZmMax);
	EH.println("         Useful log data interval ending at %-10.6g,  which has time approx %-10.6g",
			   logData.ZmMax, logData.tMax);
	if (!logData.dataExistsBelowLongestGoodInterval)
	  {
		EH.println("         WARNING! Logs not complete to end of estimation interval [%.6g, %.6g].",
				   tStartCheck, tEndCheck);
		EH.println("                  Recommend a narrower estimation interval!");
		if (!_WE.ALLOW_INCOMPLETE_LOGS)
		  {
			EH.die_err("");
		  }
	  }


	// Markers. time unit autodetection not so easy here, so assume units are the
	// same as we noticed in the checkshot. Ordering doesn't especially matter here
	for (i=0; i<Nmk; i++)
	  {
		// times may be in ms. May need to be forced into seconds....
        if (_WE.TIME_UNITS_AUTO_DETECTION)
          {
			if (timesInMilliseconds)
			  {
				Tmk[i] *= 1.0e-3;
				sigmaTmk[i] *= 1.0e-3;     
			  }
		  }
		if (Zmkm[i] < 0) EH.die_err("Marker data has -ve depths, row "+i);
		if (sigmaTmk[i] < 0) EH.die_err("Marker data has -ve time uncertainties, row "+i);
	  }

	// And form velocity estimate (for units work)
	double dZ,dT;
	dZ = (Zchm[Nch-1] - Zchm[0]);
	dT = 0.5*(Tch[Nch-1] - Tch[0]);  // since times are two-way
	meanNaturalVelocity = dZ/((dT!=0)? dT : 1);
  }

  /**
	 Useful method for getting an approximate time-depth pair from checkshot. Returns checkshot time
	 corresponding to depth z
  */
  double computeApproximateTimeDepthPair(double z)
  {
	int i;
	double T_new;
	// bracket this in checkshot...
	for (i=0; i<Zchm.length-1; i++) if (Zchm[i] > z) break; 
	// commensurate linearly interpolated point...
	if (i == 0)
	  {
		EH.warn("         WARNING: extrapolating checkshot up to obtain time estimate for out-of-range depth");
		i=1;
	  }
	if (i == Zchm.length-1)
	  {
		EH.warn("         WARNING: extrapolating checkshot down to obtain time estimate for out-of-range depth");
	  }
	T_new = ((Tch[i-1]-Tch[i])/(Zchm[i-1]-Zchm[i]))*(z-Zchm[i-1]) + Tch[i-1];
	
	return T_new;
  }


  /**
	 Generate checkshot data from full segmentation of sonic velocity log, and 
	 prepare for use in model selection.

	 Not fully debugged and properly automatic yet. See FIXME notes and issues below.
  */
  void autoGenerateCheckshots(LogData logData, GeneralWellInfo wellInfo)
  {

	SegmentedLogData sd;

	if (logData.Zm_unit == WellLogDefinitions.MD_unit.UNKNOWN || 
		logData.DT_unit == WellLogDefinitions.DT_unit.UNKNOWN )
	  {
		EH.die_err("Measured depth and DT units not available from log file "+wellInfo.logFilename+
				   " autoGenerateCheckshots() not possible. Giving up!");
	  }

	// Segment logs by slowness (1/velocity)

	// simple checks
	if (Zchm[0] >logData.Zm[0])
	  {
		EH.println("WARNING! autoGenerateCheckshots(): Well log data starts before first checkshot point");
	  }
	if (Zchm[Zchm.length-1] 
		< logData.Zm[logData.Zm.length-1])
	  {
		EH.println("WARNING! autoGenerateCheckshots(): Well log extends below last checkshot point");
	  }
	// FIXME HERE: get suitable replacement
	logData.cleanUpDataPlaceholder(Zchm[0]);

	if (_WE.verbosity>0) EH.println("Segmenting velocity log data ...well "+wellInfo.logFilename);
	sd = new SegmentedLogData(wellInfo,EH,_WE.verbosity);

	// Now segment this slowness data
	// EH.println("FIXME 4: min seg size here is 40 Compute from 1/4 lambda, plus log spacing, plus mean velocity???");
	DataSegmentor segmentor = new DataSegmentor(logData.DT,0.0f,40,-1,false);

	autoCheckshotInfo = new AutoCheckshotInfo(logData,wellInfo);
	// set permanent reference points from input files:
	autoCheckshotInfo.Zch_ref = Zchm[0];
	autoCheckshotInfo.Tch_ref = Tch[0];
	autoCheckshotInfo.sigTch_ref = sigmaTch[0];

	// times may be in ms. May need to be forced into seconds....
	if (_WE.TIME_UNITS_AUTO_DETECTION)
	  {
		if (autoCheckshotInfo.Tch_ref > FundamentalConstants.TIME_UNITS_AUTO_DETECTION_TIME_MAX)
		  {
			autoCheckshotInfo.Tch_ref *= 1.0e-3;
			autoCheckshotInfo.sigTch_ref *= 1.0e-3;
		  }
	  }

	autoCheckshotInfo.segIndexValuesArray = segmentor.dynamicProgrammingMLSegmentWithDecimationIfRequired(true);

	// set initial segmentation by most-numerous changepoints BIC choice (last model)
	autoGenSetCheckshotModel(autoCheckshotInfo.segIndexValuesArray.length-1);
  }

	/**
	   Not fully debugged and properly automatic yet. See FIXME notes and issues below.
	*/
  void autoGenSetCheckshotModel(int model_indx)
  {
	int i,j,k,n,jmin,jmax=0,cnt,kRef=-1;
	double Z,dT_tot,dt_Factor=1.0;
	boolean badData=false;
	IntArrayList segIndexValues;
	LogData logData;
	GeneralWellInfo wellInfo;

	segIndexValues = autoCheckshotInfo.segIndexValuesArray[model_indx];
	logData = autoCheckshotInfo.logData;  // useful ref.
	wellInfo = autoCheckshotInfo.wellInfo;  // useful ref.
	if (segIndexValues==null) 
	  {
		// data size probably too big: use simpler approx method
		EH.die_err("autoGenerateCheckshots(): segmentation failed...");
	  }


	// FIXME Velocity units : assuming logs in us/ft or us/m
	double Zch_ref, Tch_ref, sigTch_ref;

	if (logData.Zm_unit == WellLogDefinitions.MD_unit.M )
	  {
		if (logData.DT_unit == WellLogDefinitions.DT_unit.US_PER_FT)
		  {
			dt_Factor = 1.0/(FundamentalConstants.METRES_PER_FOOT
							 * FundamentalConstants.SONIC_LOG_TIME_SCALING);
		  }
		else if (logData.DT_unit == WellLogDefinitions.DT_unit.US_PER_M)
		  {
			dt_Factor = 1.0/FundamentalConstants.SONIC_LOG_TIME_SCALING;
		  }
		else if (logData.DT_unit == WellLogDefinitions.DT_unit.S_PER_M)
		  {
			dt_Factor = 1;			
		  }
		else
		  {
			EH.die_err("autoGenerateCheckshots(): units uninterpretable from LAS file, giving up.");
		  }
	  }
	else if (logData.Zm_unit == WellLogDefinitions.MD_unit.FT)
	  {
		if (logData.DT_unit == WellLogDefinitions.DT_unit.US_PER_FT)
		  {
			dt_Factor = 1.0/FundamentalConstants.SONIC_LOG_TIME_SCALING;	  
		  }
		else if (logData.DT_unit == WellLogDefinitions.DT_unit.US_PER_M)
		  {
			dt_Factor = FundamentalConstants.METRES_PER_FOOT
			          /FundamentalConstants.SONIC_LOG_TIME_SCALING;	  
		  }
		else if (logData.DT_unit == WellLogDefinitions.DT_unit.S_PER_M)
		  {
			dt_Factor = FundamentalConstants.METRES_PER_FOOT;
		  }
		else
		  {
			EH.die_err("autoGenerateCheckshots(): units uninterpretable from LAS file, giving up.");
		  }
	  }
	else
	  {
		EH.die_err("autoGenerateCheckshots(): units uninterpretable from LAS file, giving up.");
	  }

		  
	// reference point from checkshot first data point
	Zch_ref = autoCheckshotInfo.Zch_ref;
	Tch_ref = autoCheckshotInfo.Tch_ref;
	sigTch_ref = autoCheckshotInfo.sigTch_ref;

	// ref point must lie in segmented sonic interval, or things are hopeless...
	jmin = segIndexValues.getQuick(0);
	Z = (logData.Zm[jmin]+logData.Zm[jmin+1])/2.0 ;
	if (Zch_ref < Z)
	  {
		EH.die_err("Cannot do auto-make of checkshot with reference point above sonic log interval: Zm = "+ Zch_ref);
	  }
	jmin = segIndexValues.getQuick(segIndexValues.size()-1);
	Z = (logData.Zm[jmin]+logData.Zm[jmin+1])/2.0 ;
	if (Zch_ref > Z)
	  {
		EH.die_err("Cannot do auto-make of checkshot with reference point below sonic log interval: Zm = "+ Zch_ref);
	  }

	// From here onwards, ref point is within a segment...

	// Now make new checkshots and markers.
	// Make this as long as number of changepoints + 1, for the reference point
	n = segIndexValues.size();

	// Check if ref point and segment fall within same logging interval. 
	int indxRefKnot = -1;
	for (i=0;i<n;i++)
	  {
		j = segIndexValues.getQuick(i);
		if ((logData.Zm[j] <= Zch_ref) && (Zch_ref < logData.Zm[j+1]))
		  {
			// ref point and segment fall within same logging interval. We will amalgamate
			// the knot here
			indxRefKnot = i;
			break;
		  }
	  }

	Nch = (indxRefKnot==-1)?  n+1 : n;
	int Ninternal = Nch;

	// extra points above and below if needed
	if (autoCheckshotInfo.safeUpperCheckshot_Tch != Double.MAX_VALUE) Nch++;
	if (autoCheckshotInfo.safeLowerCheckshot_Tch != Double.MAX_VALUE) Nch++;


	// Use no markers
	Nmk = 0;

	Zchm = new double[Nch];
	Tch = new double[Nch];
	sigmaTch = new double[Nch];

	boolean refPointInLastSegment=false;
	boolean foundRefPoint=false;

	k=0;
	// possible safe upper checkshot point
	if (autoCheckshotInfo.safeUpperCheckshot_Tch != Double.MAX_VALUE)
	  {
		Zchm[k] = autoCheckshotInfo.safeUpperCheckshot_Zchm;
		Tch[k] = autoCheckshotInfo.safeUpperCheckshot_Tch;
		sigmaTch[k] = 0;
		k++;
	  }

	// Now internal members
	Tch[k] = 0;
	sigmaTch[k] = sigTch_ref;
	// FIXME
	// EH.println("FIXME 3 UNITS BELOW MAKE SENSIBLE see bit with dt_Factor");
	for (i=0;i<n-1;)
	  {
		// seg start
		jmin = segIndexValues.getQuick(i);
		if (refPointInLastSegment) 
		  {
			// second half of intersecting segment
			jmin = jmax; 
			refPointInLastSegment = false; 
		  }
		Zchm[k] = (logData.Zm[jmin]+logData.Zm[jmin+1])/2.0;  // MD

		// seg end
		jmax = segIndexValues.getQuick(i+1);
		Zchm[k+1] = (logData.Zm[jmax]+logData.Zm[jmax+1])/2.0; // MD

		if (!foundRefPoint)
		  {
			if ((Zchm[k] <= Zch_ref) && (Zch_ref < Zchm[k+1]))
			  {
				// Segment covers ref. point. Linear search for ref point in segment
				{
				  for (j=jmin; j<= jmax+1; j++)
					{
					  if ((logData.Zm[j] <= Zch_ref) && (Zch_ref < logData.Zm[j+1]))
						{
						  // have bracketed ref point with log data.
						  // as crude approximation, treat as if
						  // ref point falls precisely in the middle of the 2 samples
						  jmax = j;
						  refPointInLastSegment = foundRefPoint = true;
						  kRef = k+1;
						  break;
						}
					}
				  if (!refPointInLastSegment) EH.die_err("Failed bracketing search in autoGenerateCheckshots()! BUG! FIXME!");
				}
			  }
		  }
		// Sum time in segment j=jmin+1....jmax
		dT_tot = 0;
		for (cnt=0,j = jmin+1; j <= jmax; j++,cnt++)
		  {
			if (logData.DT[j] <= FundamentalConstants.BAD_NUMBER)
			  {
				badData = true;
				break;
			  }
			dT_tot += logData.DT[j]*(logData.Zm[j+1] - logData.Zm[j]);
		  }
		// dT_tot *= (logData.Zm[jmax] - logData.Zm[jmin]);
		if (badData)
		  {
			EH.die_err("Cannot do auto make of checkshot with missing sonic log data, depth "+ logData.Zm[j]);
		  }
		// Here dT_tot has units [DT sonic-units] x [welllog meas.depth units], e.g. (uS/ft)*m
		// Calculation here assumes logs in uS/ft [or S/m], need to convert to S
		// points at end of bracket		
		
		Tch[k+1] = Tch[k] + 2.0 * dt_Factor * dT_tot; // since times are two-way
		sigmaTch[k+1] = sigTch_ref;
		// advance if segment spanned at least one log measurement
		if (cnt>0) k++;

		// goto next segment unless examining second half of segment intersecting ref point
		if (!refPointInLastSegment) i++;
	  }
	if (kRef==-1) EH.die_err("Failed to set reference time in autoGenerateCheckshots()! BUG! FIXME!");


	// Now re-reference with respect to reference point
	dT_tot = Tch_ref - Tch[kRef];
	// skip fixed safe upper checkshot pt if required
	k = (autoCheckshotInfo.safeUpperCheckshot_Tch == Double.MAX_VALUE) ? 0 : 1;
	for (j=0;j<Ninternal;j++,k++)
	  {
		// Tch[k] <- (Tch[k] - Tch[kRef]) + Tch_ref;
		Tch[k] += dT_tot;
	  }

	if (autoCheckshotInfo.safeLowerCheckshot_Tch != Double.MAX_VALUE)
	  {
		k = Nch-1;
		// possible safe upper checkshot point
		Zchm[k] = autoCheckshotInfo.safeLowerCheckshot_Zchm;
		Tch[k] = autoCheckshotInfo.safeLowerCheckshot_Tch;
		sigmaTch[k] = 0;
	  }


	// meanNaturalVelocity = ???
	double dZ,dT;
	dZ = (Zchm[Nch-1] - Zchm[0]);
	dT = 0.5*(Tch[Nch-1] - Tch[0]);  // since times are two-way
	meanNaturalVelocity = dZ/((dT!=0)? dT : 1);

	PrintWriter out;
	String filename = new String(Util.fileBasename(wellInfo.checkshotFilename)+"_automake-test."+model_indx+".txt");
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }
	out.printf("Full optimally-segmented checkshot data from reference point (line %d)%n3%nMD%nT%nsigmaT%n",(kRef+1));	
	for (k=0;k<Nch;k++)
	  {
		out.printf("%10.8g %10.8g %10.8g%n",Zchm[k],Tch[k],sigmaTch[k]);
	  }
	out.close();
	EH.println("Wrote test autogenerated checkshot "+filename);
  }




  /**
	 A useful checking method for generating fake checkshot data from the upscaled well log.
	 Development purposes only.
	 Here we assume all the segmented data is in place, including the true depths.
	 The velocity is taken as constant up to the first logged data
  */
  void GenerateFakeCheckshotdata(SegmentedLogData segLogData, GeneralWellInfo wellInfo)
  {
	int i,j,N;
	double[] times;
	double z,zm,t;
	double ZmMin,ZMin;
	double ZSPACING=500;          // 500ft spaced checkshot
	Spline zSpline, tSpline;

	N = segLogData.Z.length;
	times = new double[N];
	// downward compute the times, taking all the slab above the
	// upscaled chunk to have the same velocity as the top layer
	times[0] = segLogData.Z[0]/segLogData.Vp[0];
	for (i=1;i<N;i++)
	  {
		times[i] = (segLogData.Z[i]-segLogData.Z[i-1])/segLogData.Vp[i-1] + times[i-1];
	  }
	// now make some useful splines
	// will need to do conversions   Zm  -> Z -> t
	zSpline = new Spline(segLogData.Zm, segLogData.Z,segLogData.Z.length, Double.MAX_VALUE, Double.MAX_VALUE,EH);
	tSpline = new Spline(segLogData.Z,times, segLogData.Z.length, Double.MAX_VALUE, Double.MAX_VALUE,EH);

	// now, fake checkshot at DeltaZ spacing to hole bottom
	// Write as a pseudo GEOEAS file
	String cfilename = Util.fileBasename(wellInfo.checkshotFilename);
	EH.println("Printing checkshots to : "+cfilename);
	PrintWriter out;
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(cfilename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + cfilename);
	  }
	out.println("Fake checkshot data%n3%nMD%nT%nsigmaT");	

	zm=0;
	ZmMin = segLogData.Zm[0];
	ZMin = segLogData.Z[0];
	for (zm=0; zm < segLogData.Zm[segLogData.Zm.length - 1] ; zm += ZSPACING)
	  {
		z = (zm>=ZmMin) ? zSpline.splint(zm,true) : zm*ZMin/ZmMin;  // true depth from measured depth
		t = (z>=ZMin) ? tSpline.splint(z,true) : z*times[0]/ZMin;   // time from true depth	
		out.println(zm+" "+t+" 0.002"); // checkshot data pretend there's 2ms uncertainty
	  }
	// make sure there's a checkshot at the end of the log data, so all markers will be bracketed
	z = zSpline.splint(segLogData.Zm[segLogData.Zm.length - 1],true);  // true depth from measured depth
	t = tSpline.splint(z,true);   // time from true depth	
	out.println(zm+" "+t+" 0.002"); // checkshot data pretend there's 2ms uncertainty
	out.close();
  } 

  /**
	 A useful checking method for generating fake marker data from the Checkshot data
	 and the upscaled well log.
	 Development purposes only.
  */
  void GenerateFakeMarkerdata(SegmentedLogData segLogData, GeneralWellInfo wellInfo)
  {
	int i,j,k,N,index;

	double zm,wt,t,dt,contrast;
	double NMARKERS=5;          // generate markers at 5 largest impedance-contrast points
	SortableNumbers[] cons;

	N = segLogData.Z.length;
	// compute and sort the impedance contrasts
	cons = new SortableNumbers[N-1];
	cons[0] = new SortableNumbers(0,0);
	for (i=1;i<N-1;i++)
	  {
		contrast = (segLogData.rho[i]*segLogData.Vp[i] - segLogData.rho[i-1]*segLogData.Vp[i-1])/
		  (segLogData.rho[i]*segLogData.Vp[i] + segLogData.rho[i-1]*segLogData.Vp[i-1]);
		cons[i]= new SortableNumbers(Math.abs(contrast),i);
	  }
	// sort into increasing order
	Arrays.sort(cons);
	
	// now, fake markers at most important reflections
	// Write as a pseudo GEOEAS file
	String mfilename = wellInfo.markersFilename;
	EH.println("Printing markers to : "+mfilename);
	PrintWriter out;
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(Util.fileBasename(mfilename))));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + mfilename);
	  }

	out.println("Fake marker data%n4%nMD%nDeltaT%nsigmaDeltaT");	
	for (i=0,j=N-2;i<NMARKERS;i++,j--)
	  {
		// true depth
		index = cons[j].index;
		zm = segLogData.Zm[index];

		// may need these...this doesn't work yet, but we don't need it here
		// bounding index in checkshot array
		//for (k=-1;zm>Zchm[k+1];k++) {}
		// now Zchm[k] < zm < Zchm[k+1]
		//wt = (Zchm[k+1] - zm)/(Zchm[k+1] -  Zchm[k]);
		// linearly interpolated value
		//t = wt * Tch[k]+(1-wt)*Tch[k+1];

		// generate some kind of deviate for the marker wrt this time
		// take 4ms uncertainty here
		dt = (Math.random()-0.5)*0.004;
		out.println(zm+" "+dt+" 0.004"); // marker data
	  }
	out.close();
  }

  /**
	 Class containing all the useful information required
	 in autoGenerated checkshots
  */
  public static class AutoCheckshotInfo
  {
	// link to underlying log data
	LogData logData;
	GeneralWellInfo wellInfo;
	// permanent reference points from checkshot file first line
	double Zch_ref, Tch_ref, sigTch_ref;

	// Safe upper and lower checkshots of sigmaT=0 if needed
	// the MAX_VALUE flags "unused"
	double safeUpperCheckshot_Zchm;
	double safeUpperCheckshot_Tch = Double.MAX_VALUE;
	double safeLowerCheckshot_Zchm;
	double safeLowerCheckshot_Tch = Double.MAX_VALUE;

	// Array of segments
	IntArrayList[] segIndexValuesArray=null;

	AutoCheckshotInfo(LogData l,GeneralWellInfo wI)
	{
	  this.logData = l;
	  this.wellInfo = wI;
	}
  }
}
