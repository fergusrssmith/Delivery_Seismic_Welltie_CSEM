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
 
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import static  au.csiro.JamesGunning.lib.util.MathUtil.sqr;

import cern.colt.list.*;
import cern.jet.stat.*;
import cern.jet.random.*;
import cern.jet.random.engine.MersenneTwister;

import java.util.*;
import java.io.*;
/**
   Basic model variables for wavelet extraction.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class Model
{
  /**
	 ptr to underlying wavelet extractor object and thus to all its fields
  */
  WaveletExtractor _WE;
  int nWells;
  int nStacks=0;

  Petrophysics petro = null;
  TimeLapseInfo timelapseinfo=null;
  
  // flag is switched on if multi-stack registration errors are common or there is only one stack
  boolean registrationErrorsAreCommon = false;

  /**
	 Switch this on to use log-data consistent t<->d conversion for
	 times intermediate between the major checkshot/marker points. If
	 off, linear interpoaltion is used. The latter is probably near
	 enough if the markers are in good places.
  */
  boolean useFineSegmentedInterpolation = false;

  /**
	 flag to indicate status of current time ordering map
  */
  boolean timesAreDecentlyOrdered = true;

  /**
	 If verbosity requires angle dumps, do this only once.
   */
  private boolean haveDumpedAngles = false;

  TimeToDepthMap[] timeToDepthMaps;
  /**
	 Array of noise levels, one per stack
  */
  Property[] seismic_noise;

  Wavelet wavelet;
  /**
	 For debugging purposes, sometimes the X sections show
	 discontinuities, which screw up the minimiser and wreak
	 havoc. These are bugs. To help detect these, switch this on, find
	 "debuggingXSections" in the source, and do stuff in the chunks of code dependent on this...
	 Switch off when debugging finished.
  */
  boolean debuggingXSections = false;

  boolean debug; 

  int debug_cntr=0;

  double[] debugArray=null;     // count of number of segment tops found in current interval, if debug=true

  /**
	 switch for amplitude caching
  */
  public static boolean USE_CACHED_REFLECTIVITY = true;

  /**
	 List of field titles output in properties.su file, as read by BHPcommand
  */
  String[] propertyNames;

  /**
	 the reflectivity B factor (coefficient of theta^2) in the Zoeppritz eqn.
	 elemnet 0 is the shared one if in common, or an extra one for each stack is created if
	 free
  */
  Property[] B_AVO;

  // Master property
  Property deltaThomsen_shale; 

  // work variables
  double[] farWavStretch;
  double[] farWavAmp;

  /**
	 Factor to model wavelet stretch
  */
  Property[] farStackWaveletTimeStretch;
  /**
	 Useful hashmap of these
  */
  ArrayList farStackWaveletTimeStretchList = null;

  /**
	 Fixed variables to model wavelet stretch: these retain the values 
	 in the XML, since this property goes through a logit map which
	 alters the means/std deviations
  */
  double farStackWaveletTimeStretchMean=1.05;
  double farStackWaveletTimeStretchSigma=0;

  /**
	 XML parsed amplitude params for the far wavelets
  */
  double farStackWaveletAmplitudeMean=1;
  double farStackWaveletAmplitudeSigma=0;


  /**
	 extra factor to model wavelet stretch
  */
  Property[] farStackWaveletAmplitude;

  /** 
	 A remapper for the constrained stretch factor if needed
  */

  PropertyMapper farStackWaveletTimeStretchMapper;


  private StackAngleLoggingData[] stackAngleLoggingData;
  private int criticalReflectionErrorLoggingCount = 0;

  /**
	 Print critical reflection errors to stdout a maximum of this many times
  */
  public static final byte MAX_CRITICAL_REFL_LOGGING_COUNT = 20;

  /*
	 Saturation handling entities
  */
  SaturationInfo[] satInfo=null;

  public Model(WaveletExtractor WE)
  {
	int i;

	this._WE = WE;

	petro = new Petrophysics();  // default member hasAnisotropy = false

	debug=false;     // LEAVE OFF UNLESS DEBUGGING!
	if (debug)       _WE.EH.println("TURN OFF SEG DEBUGGING FOR PRODUCTION!");
  }

  /**
	 Set up stretch parameters as soon as possible, i.e when nStacks is known
  */
  void initialiseStretchParameters()
  {
	int i;

	if (_WE.freeAVO)
	  {
		B_AVO = new Property[nStacks];
		for (i=0; i<nStacks; i++)
		  {
			B_AVO[i] = new Property(1.0,false,true);
			// set mean to default to 1.0
			B_AVO[i].setMean(1.0);
		  }
	  }
	else
	  {
		B_AVO = new Property[1];
		B_AVO[0] = new Property(1.0,false,true);
		// set mean to default to 1.0
		B_AVO[0].setMean(1.0);
	  }

	farStackWaveletTimeStretch = new Property[nStacks];
	farStackWaveletTimeStretchList = new ArrayList();
	farStackWaveletAmplitude = new Property[nStacks];
	for (i=0; i<(nStacks-1); i++)
	  {
		farStackWaveletTimeStretch[i] = new Property(1.0,false,false);
		farStackWaveletTimeStretch[i].setMean(1.0);
		farStackWaveletTimeStretchList.add(farStackWaveletTimeStretch[i]);
		
		farStackWaveletAmplitude[i] = new Property(1.0,false,false);
		farStackWaveletAmplitude[i].setMean(1.0);
	  }
	
	if (nStacks > 1)
	  {
		farWavAmp = new double[_WE.m.nStacks];
		farWavStretch = new double[_WE.m.nStacks];
	  }
  }


  /** 
	  Set up common Properties shared by all wavelet models. Call only after
	  number of stacks is known and a fair guess for a characteristic amplitude
	  is known
  */
  void setUpCommonProperties(double seismicAmplitudeGuess)
  {
	seismic_noise = new Property[nStacks];
	for (int i=0; i<nStacks; i++)
	  {
		double noiseGuess = Math.abs(seismicAmplitudeGuess);
		seismic_noise[i] = new 
		  Property(noiseGuess, false,true);
	  }
  }


  
  /**
	 Determine if all wells are "synchronised", and if so, what the
	 indx of the master well is. If not return -1
  */
  int getCommonTregSynchIndx()
  {
	int j;
	boolean hasCommonTreg;
	GeneralWellInfo _WI;
	int synchIndx = 0;

	// determine if all wells synchronised
	// First, find a potential "master"
	for (j=0; j<_WE.m.nWells; j++)
	  {
		if (_WE.wellInfo[j].synchronizeWellName == null) {synchIndx = j; break;}
	  }
	// synchIndx is first indx of well that is a "master". All others *must* point to it
	// assume all registration effects synchronised unless we find a contradiction
	hasCommonTreg = true;
	for (j=0; j<_WE.m.nWells; j++)
	  {
		if (j==synchIndx) continue;
		_WI = _WE.wellInfo[j];
		if (_WI.synchronizeWellName == null) // different well & looks free
		  {
			hasCommonTreg = false;
			break;
		  }
		if (!_WI.synchronizeWellName.equals(_WE.wellInfo[synchIndx].well_name)) hasCommonTreg = false;
	  }
	if (hasCommonTreg) return synchIndx;
	else return -1;
  }

  /**
	 Set current model variables for specified wavelet. 
	 This ultimately will point all the pointers in Model to
	 the relevant fields in the current model.
  */
  void setSharedVariables(ModelVariables theseModelVariables)
  {
	// set wavelet
	setWavelet(theseModelVariables);

	// set anisotropy mode
	petro.hasAnisotropy = theseModelVariables.hasAnisotropy;
  }

  /**
	 Set current model variables for specified wavelet. 
	 This ultimately will point all the pointers in Model to
	 the relevant fields in the current model.
  */
  void setWavelet(ModelVariables theseModelVariables)
  {
	wavelet =  theseModelVariables.wavelet;
  }



  /**
	 Updates time to depth maps based on current values in the 
	 Properties fields
  */
  double updateTimeToDepthMaps()
  {
	double timesMisOrderingPenalty=0.0;
	timesAreDecentlyOrdered=true;
	
	for (int k=0;k<nWells;k++)
	  {
	    timesMisOrderingPenalty += timeToDepthMaps[k].updateTimesAndReturnMisorderingPenalty();
		if (! timeToDepthMaps[k].areCurrentTimesDecentlyOrdered()) timesAreDecentlyOrdered=false;
	  }
	return timesMisOrderingPenalty;
  }


  /** 
	  Recomputes log properties under saturation shifts 
  */
  void recomputeLogProperties(boolean after)
  {
	int i;   
	for (i=0;i<_WE.m.nWells;i++)
	  {
		_WE.sLogData[i].recomputeSegProperties(_WE.m,_WE.dumpPandSreflectivity,_WE.doShuey,after);
	  }
  }

  /**
	 Reflectivity copy
  */
  
  void makeReflectivityCopy()
  {
	int i,k;
	ConvolutionalModel.BaseVintageReflectivities _BVR = _WE.convModel.BVR;

	// FIXME copy here, not clone
	_BVR.Reff = ArrayUtils.deepCopy(_WE.convModel.Reff);

  }

 
  /**
	 Compute the sequence of reflectivities at the segmented
	 boundaries, for all wells, for all stacks, on the time axis. Much
	 interpolation involved, so operates best if range checks etc are
	 enforced with calls to the methods checklogDataRanges() and
	 checkTimeToDepthRanges() below.
  */
  void computeReflectivities(boolean useCache, boolean computePandSreflectivities)
  {
	int i,j,k,indx,iStart,iEnd,iSlogStart=0,iSlogEnd=0,jEnd;
	int firstSegmentedLogIndex=0;
	double alpha=1, sum, dZ, Vp, zLeft, zRight, tLeft;
	double A_neg1,A_0,A_1,A_2,tSegj,wt,tmp1,tmp2;
	double Vdenom, ref_coeff=0.0, dR,sqrt_Vdenom;
	TimeToDepthMap t2dmap=null;
	int Zsz;
	double tStartUsefulSegments, tEndUsefulSegments;
	double[] segTimes=null;     // mapping of layer tops to times, using current model (from segmented data)
	double[] segDepths=null;    // layer tops in true depth (from segmented data)
	double[] segVelocity=null;  // layer velocities (from segmented data)
	double[] _R_near=null;      // ptr to near ref coeff for current log
	double[] R1_AVO = null;     // ptr to AVO sin^2(theta) term
	double[] Vsq_AVO = null;    // save of segments vp_1^2
	double[] R=null;            // working space ptr
	double[] R_P=null;          // working space ptr for extra P reflectivity
	double[] R_S=null;          // working space ptr for extra S reflectivity
	double[] R_Rho=null;        // working space ptr for extra Density reflectivity
	double[] deltaZs=null;      // ptr to segmented log deltaVs
	double[] deltaRho=null;     // ptr to segmented log deltaRho
	double[] _Rps_AVO=null;     // ptr to segmented Rps_AVO term1 
	double[] _Rps_AVO2=null;    // ptr to segmented Rps_AVO term 2 

	SeismicStackInfo _SSI;
	// Shuey
	double[] R2_Shuey_AVO = null;
	// Zoeppritz pointers if used
	double[] _Vp = null,_Vs = null,_rho = null;

	double sth_sq=0;

	boolean stackAnglePrescribed = false;
	boolean isPSdata = false;

	boolean intervalHasCompleteSegmentedData=false;
	boolean intervalHasSomeData=false;
	boolean seekingFirstInterval;

	// Logging
	boolean dumpAngles=false;
	if (!haveDumpedAngles && _WE.dumpAnglesToFiles)
	  {
		dumpAngles = true;
		if (!haveDumpedAngles)
		  {
			stackAngleLoggingData = new StackAngleLoggingData[nWells];
			for (k=0;k<nWells;k++) stackAngleLoggingData[k] = new StackAngleLoggingData(nStacks);
		  }
	  }
	// If we use Cache, all the in-situ values are accepted as valid. Take care!!!!
	// _WE.EH.println("RCACHE "+(useCache?1:0));
	if (!(useCache && USE_CACHED_REFLECTIVITY))
	  {
		for (k=0;k<nWells;k++)
		  {
			t2dmap = timeToDepthMaps[k];
			SegmentedLogData _sLogData = _WE.sLogData[k];

			segTimes = _sLogData.T;
			segDepths = _sLogData.Z;
			segVelocity = _sLogData.Vp;
			if (debug) debugArray = new double[segDepths.length];
			Zsz = segDepths.length-1;


			// tStartUsefulSegments = _WE.tStart;
			// tEndUsefulSegments = _WE.tEnd;
			tStartUsefulSegments = _WE.convModel.tStartLegalRange;
			tEndUsefulSegments = _WE.convModel.tEndLegalRange;

			// (1) find data points in timeToDepth map that bracket tStart & tEnd
			iStart = getStraddlingLowerIndex(tStartUsefulSegments,t2dmap.t);
			if (iStart<0)
			  {
				_WE.EH.warn("tStart out of range!");
			  }
			iEnd = 1 + getStraddlingLowerIndex(_WE.tEnd,t2dmap.t);

			if (iEnd >= t2dmap.t.length)
			  {
				// we might be out-of-range here, or, the bracketing procedure was needlessly conservative
				// Check the latter to see if we can get away without a rigorous double bracket:
				// (i.e. is setting iEnd = t2dmap.t.length-1 acceptable?)
				if (segTimes[segTimes.length-1] > _WE.tEnd)
				  {
					// yes we can get away with it
					iEnd = t2dmap.t.length-1;
				  }
				else
				  {
					// do it anyway, but complain
					iEnd = t2dmap.t.length-1;
					_WE.EH.warn("iEnd out of range!");
				  }
			  }

			// (2) loop over adjacent-pairs in the time depth map
			// interpolating the set of segmented-interval reflections
			seekingFirstInterval = true;
			for (i=iStart ; i<iEnd ; i++)
			  {
				// Major bracket left & right edge depth				
				zLeft = t2dmap.Z[i];
				zRight= t2dmap.Z[i+1];

				// now bracket this range in the segmented logs
				// look out for missing data...
				intervalHasCompleteSegmentedData = false;
				intervalHasSomeData = true;
				if (seekingFirstInterval) 
				  {
					// first bracket, or the last one didn't work
					iSlogStart = getStraddlingLowerIndex(zLeft,segDepths);
					iSlogEnd = getStraddlingLowerIndex(zRight,segDepths);
					intervalHasCompleteSegmentedData = false;
					// there are six possibilities here
					if (iSlogStart<0 && iSlogEnd<0)
					  {					
						// no data at all
						intervalHasSomeData=false;
					  }
					else if (iSlogStart<0 && iSlogEnd>=0 && iSlogEnd<Zsz)
					  {
						// some data
						seekingFirstInterval = false;
						iSlogStart=0;
					  }
					else if (iSlogStart<0 && iSlogEnd>=Zsz)
					  {
						// some data
						seekingFirstInterval = false;
						iSlogEnd = Zsz;
					  }
					else if (iSlogStart>=0 && iSlogStart<Zsz && iSlogEnd>=0 && iSlogEnd<Zsz)
					  {
						// some data & it's complete
						seekingFirstInterval = false;
						intervalHasCompleteSegmentedData = true;
					  }
					else if (iSlogStart>=0 && iSlogStart<Zsz && iSlogEnd>=Zsz)
					  {
						// some data
						seekingFirstInterval = false;
						iSlogEnd = Zsz;
					  }
					else if (iSlogStart>=Zsz && iSlogEnd>=Zsz)
					  {
						// no data
						intervalHasSomeData=false;
					  }
					// we might be out-of-range here, or, the bracketing procedure was needlessly conservative
					// Check the latter to see if we can get away without a rigorous double bracket:
					// (i.e. is setting iSlogStart=0 acceptable?)

					// By here, iSlogStart && iSlogEnd are at safe values

					if (!seekingFirstInterval) firstSegmentedLogIndex = iSlogStart+1;
				  }
				else
				  {
					// subsequent bracket
					iSlogEnd = getStraddlingLowerIndex(zRight,segDepths);
					intervalHasCompleteSegmentedData = false;
					// there are six possibilities here
					if (iSlogStart<0 && iSlogEnd<0)
					  {					
						// no data at all	
						intervalHasSomeData=false;
					  }
					else if (iSlogStart<0 && iSlogEnd>=0 && iSlogEnd<Zsz)
					  {
						// some data
						iSlogStart=0;
					  }
					else if (iSlogStart<0 && iSlogEnd>=Zsz)
					  {
						// some data
						iSlogEnd = Zsz;
					  }
					else if (iSlogStart>=0 && iSlogStart<Zsz && iSlogEnd>=0 && iSlogEnd<Zsz)
					  {
						// some data & it's complete
						intervalHasCompleteSegmentedData = true;
					  }
					else if (iSlogStart>=0 && iSlogStart<Zsz && iSlogEnd>=Zsz)
					  {
						// some data
						iSlogEnd = Zsz;
					  }
					else if (iSlogStart>=Zsz && iSlogEnd>=Zsz)
					  {
						// no data
						intervalHasSomeData=false;
					  }
				  }
				if (debug) debugArray[i]=iSlogEnd -iSlogStart ;
				if (intervalHasSomeData)
				  {
					if (useFineSegmentedInterpolation && intervalHasCompleteSegmentedData)
					  {
						// Can only use useFineSegmentedInterpolation if we have complete data
						// first, get well log-adjustment velocity factor alpha
						dZ = segDepths[iSlogStart+1] - zLeft;
						Vp= segVelocity[iSlogStart];
						sum = dZ/Vp;
						for (j=iSlogStart+1;j<iSlogEnd;j++)
						  {
							dZ = segDepths[j+1] - segDepths[j];
							Vp = segVelocity[j];
							sum += dZ/Vp;
						  }
						dZ = zRight - segDepths[iSlogEnd];
						Vp = segVelocity[iSlogEnd];
						sum += dZ/Vp;
						alpha = (t2dmap.t[i+1] - t2dmap.t[i])/sum;			
						// ok now loop over subsections, computing times
						dZ = segDepths[iSlogStart+1] - zLeft;
						Vp = segVelocity[iSlogStart];
						sum = dZ/Vp;
						segTimes[iSlogStart+1] = t2dmap.t[i] + alpha*sum;
						for (j=iSlogStart+1; j<iSlogEnd; j++)
						  {
							dZ = segDepths[j+1] - segDepths[j];
							Vp= segVelocity[j];
							sum += dZ/Vp;
							segTimes[j+1] = t2dmap.t[i] + alpha*sum;
							//System.out.println(segDepths[j]+" "+Vp);
						  }
					  }
					else
					  {
						// Loop over subinterval of checkshot data, computing intermediate 
						// times at segmentation boundaries by simple linear interpolation
						// This ignores the velocity variations on internal segments
						tmp1 =  (t2dmap.t[i+1]-t2dmap.t[i])/(zRight-zLeft);
						for (j=iSlogStart+1; j<=iSlogEnd; j++)
						  {
							segTimes[j] =  t2dmap.t[i] + tmp1*(segDepths[j] - zLeft);
						  }
					  }
				  }
				if (!seekingFirstInterval)
				  {
					// increment bracketing indices
					iSlogStart = iSlogEnd;
				  }			
			  }



			// Now compute sequence of reflection coefficients for each stack

			_R_near =_sLogData.R_near;
			R1_AVO = _sLogData.R1_AVO;
			Vsq_AVO = _sLogData.Vsq_AVO;
			if (_WE.doShuey)
			  {
				R2_Shuey_AVO = _sLogData.R2_Shuey_AVO;
			  }
			else if (_WE.doIsoZoeppritz)
			  {
				_Vp = _sLogData.Vp;
				_Vs = _sLogData.Vs;
				_rho = _sLogData.rho;
			  }
			if (_WE.wellInfo[k].needPSreflectivities)
			  {
				_Vp = _sLogData.Vp;
				_Rps_AVO = _sLogData.Rps_AVO;
				_Rps_AVO2 = _sLogData.Rps_AVO2;
				/*
				  // for ***Scat*** below
				  _Vs = _sLogData.Vs;
				  _rho = _sLogData.rho;
				*/
			  }
			if (computePandSreflectivities)
			  {
				deltaZs = _sLogData.deltaZs;
				deltaRho= _sLogData.deltaRho;
				R_P = _WE.convModel.Reff_P[k];
				R_S = _WE.convModel.Reff_S[k];
				R_Rho = _WE.convModel.Reff_Rho[k];
				Arrays.fill(R_P,0.0);
				Arrays.fill(R_S,0.0);
				Arrays.fill(R_Rho,0.0);
			  }


			for (i=0; i<nStacks; i++)
			  {
				_SSI = _WE.seismicInfo[k][i];
				stackAnglePrescribed = _SSI.isStackAnglePrescribed();
				isPSdata = _SSI.isPSdata;
				Vdenom = _SSI.getVdenom();
				sqrt_Vdenom = _SSI.getsqrtVdenom();				  
				R = _WE.convModel.Reff[k][i];
				Arrays.fill(R,0.0);

				double coeff,B;
				B = _WE.freeAVO ? B_AVO[i].value : B_AVO[0].value;
				if (_WE.doShuey || _WE.doIsoZoeppritz)
				  {
					coeff = B;
				  }
				else
				  {
					if (stackAnglePrescribed)
					  {
						coeff = B * MathUtil.sqr(_SSI.theta_specified);
					  }
					else
					  {
						coeff = B/Vdenom;
					  }
				  }
				// loop over segmented layer tops, where reflections occur
				jEnd = iSlogEnd;
				if (_WE.doIsoZoeppritz) jEnd = Math.min(jEnd, _rho.length-1); // One short item in _rho,_Vp arrays etc.
				// jEnd = Math.min(jEnd, _rho.length-1);  // **Scat** needed for test below
				for (j=firstSegmentedLogIndex; j<= jEnd ; j++)
				  {
					//				if (debug) debugArray[j] = segTimes[j];
					tSegj = segTimes[j];

					if (tSegj >= _WE.convModel.tStartLegalRange 
						&& tSegj <_WE.convModel.tEndLegalRange)
					  {
						// if (debug) debugArray[j]+=1.0;
						// Master reflection Zoeppritz equation
						if (_WE.doShuey || _WE.doIsoZoeppritz)
						  {
							// sin^2(theta_1)
							sth_sq = stackAnglePrescribed ? _SSI.sinsq_theta_specified : Vsq_AVO[j]/Vdenom;
										  
							if (sth_sq>=1)
							  {
								logCriticalReflections(0,sth_sq,i,tSegj);
								ref_coeff = 0;
							  }
							else  
							  {
								if (!isPSdata) 
								  {
									if (_WE.doShuey)
									  {
										// Shuey equation:   NB. tan^2 sin^2 = sin^4/(1-sin^2)
										ref_coeff = _R_near[j] + 
										  coeff * (R1_AVO[j] * sth_sq + R2_Shuey_AVO[j]*sth_sq*sth_sq/(1.0-sth_sq));
									  }
									else
									  {									
										// Zoeppritz
										ref_coeff = GeophysicsUtils.zoeppritz_Rpp_Rps(false,sth_sq,_Vp[j-1],_Vp[j],_Vs[j-1],_Vs[j],_rho[j-1],_rho[j]);
										if (ref_coeff == -Double.MAX_VALUE)
										  {
											logCriticalReflections(1,sth_sq,i,tSegj);
											ref_coeff = 0;
										  }
									  }
								  }
								else
								  {
									// Haven't done Shuey RPS yet. Do Zoeppritz for both...
									// Zoeppritz
									ref_coeff = GeophysicsUtils.zoeppritz_Rpp_Rps(true,sth_sq,_Vp[j-1],_Vp[j],_Vs[j-1],_Vs[j],_rho[j-1],_rho[j]);
									if (ref_coeff == -Double.MAX_VALUE)
									  {
										logCriticalReflections(1,sth_sq,i,tSegj);
										ref_coeff = 0;
									  }
								  }
							  }
						  }						
						else
						  {
							// Original linearized O(theta^2) approximations
							if (!isPSdata) 
							  {
								// PP coefficient
								dR = coeff * R1_AVO[j];
								if (!stackAnglePrescribed) dR *= Vsq_AVO[j];
								ref_coeff = _R_near[j] + dR;
							  }
							else
							  {
								// PS coefficient, linearised, approximating theta1 \approx theta
								if (!stackAnglePrescribed) 
								  {
									double theta1 = Math.sqrt(Vsq_AVO[j]/Vdenom);  // units-corrected, correct theta1 value
									ref_coeff = theta1 * (_Rps_AVO[j] + sqr(theta1) *_Rps_AVO2[j]); // Linearised Aki-Richards, assuming theta1 ~ theta
								  }
								else 
								  {
									ref_coeff = _SSI.theta_specified * (_Rps_AVO[j] + sqr(_SSI.theta_specified) *_Rps_AVO2[j]);
								  }
							  }
						  }

						if (petro.hasAnisotropy) // _WE.doAnisotropicModel
						  {
							// watch out for unitsConversionFactor attached to DeltaR_anis (like Vsq_AVO too)
							ref_coeff += deltaThomsen_shale.value * _sLogData.DeltaR_anis[j]/Vdenom;
						  }

						if (dumpAngles)
						  {
							// Show angles 
							if (!(_WE.doShuey || _WE.doIsoZoeppritz))
							  {
								// sth_sq = Vsq_AVO[j]/Vdenom; // sin^2(theta_1) (not computed above)
								// sin^2(theta_1) (not computed above)
								sth_sq = stackAnglePrescribed ? _SSI.sinsq_theta_specified : Vsq_AVO[j]/Vdenom;
							  }
							if (i==0) stackAngleLoggingData[k].times.add(tSegj);
							stackAngleLoggingData[k].theta[i].add(Math.toDegrees(Math.asin(Math.sqrt(sth_sq))));
						  }


						// Distribute this reflection, which occurs at segTimes[j], onto the
						// regular-spaced arrays Reff

						indx = (int)((tSegj - _WE.convModel.tStartReff)/_WE.dt);

						// left bracketing time on array of equispaced times for Synthetic
						tLeft = (indx*_WE.dt+_WE.convModel.tStartReff);
						wt = (tSegj-tLeft)/_WE.dt;
						if (true)
						  {
							// Use four-point Lagrange formula. 
							// Compute the four coefficients A_{-1},A_0,A_1,A_2
							// as per p 879 of Abramowitz & Stegun....// N.B. 1/6 = 0.166666..
							A_neg1 = -wt*(wt-1.0)*((tmp1 = wt-2.0))*0.166666666666;
							A_0 = ((tmp2 = wt*wt-1.0))*tmp1*0.5;
							A_1 = -wt*(wt+1.0)*tmp1*0.5;
							A_2 = wt*tmp2*0.166666666666;
							// Lagrange 4-pt  interpolation
							R[indx+2] += A_2 * ref_coeff;
							R[indx+1] += A_1 * ref_coeff;
							R[indx] += A_0 * ref_coeff;
							R[indx-1] += A_neg1 * ref_coeff;

							// extra reflectivities if desired
							if (computePandSreflectivities && (i==0))
							  {
								// Will only be called if _WE.dumpPandSreflectivity is true
								// (1) P: R_near is the normal incidence P-P guy
								ref_coeff = _R_near[j];
								// Lagrange 4-pt  interpolation
								R_P[indx+2] += A_2 * ref_coeff;
								R_P[indx+1] += A_1 * ref_coeff;
								R_P[indx] += A_0 * ref_coeff;
								R_P[indx-1] += A_neg1 * ref_coeff;

								// (2) S: deltaZs is the S-S impedance delta
								ref_coeff = deltaZs[j];
								// Lagrange 4-pt  interpolation
								R_S[indx+2] += A_2 * ref_coeff;
								R_S[indx+1] += A_1 * ref_coeff;
								R_S[indx] += A_0 * ref_coeff;
								R_S[indx-1] += A_neg1 * ref_coeff;

								// (3) Rho: deltaRho is the Density reflection coefficient
								ref_coeff = deltaRho[j];
								// Lagrange 4-pt  interpolation
								R_Rho[indx+2] += A_2 * ref_coeff;
								R_Rho[indx+1] += A_1 * ref_coeff;
								R_Rho[indx] += A_0 * ref_coeff;
								R_Rho[indx-1] += A_neg1 * ref_coeff;
							  }
						  }
						else
						  {
							// sinc formula for testing: very demanding but differentiable
							for (int kk =-indx;;kk++)
							  {
								if (indx+kk >= R.length) break;
								double sinc_wt=0;
								double arg = (kk-wt)*Math.PI;
								sinc_wt = (arg!=0) ? Math.sin(arg)/arg : 1.0;
								R[indx+kk] += sinc_wt * ref_coeff;
							  }
						  }
					  }
				  }
				//Util.printfileZO("R.txt",R,R.length); // FIXME REMOVE
				//Util.printfileZO("Rnear.txt",_R_near,_R_near.length); // FIXME REMOVE
			  }
		  }
	  }
	if (dumpAngles)
	  {
		// mark this so subsequent calls don't dump angles
		writeOutAngles();
		haveDumpedAngles = true;
		// GC: don't needthis anymore
		stackAngleLoggingData = null;
	  }
  }

 void computeReflectivityDifference()
  {
	int i,j,k,N;
	double[] _R,_R0;
	ConvolutionalModel.BaseVintageReflectivities _BVR = _WE.convModel.BVR;

	for (i=0;i<nWells;i++)
	  {
		for (j=0; j<nStacks; j++)
		  {			  
			_R0 = _BVR.Reff[i][j];
			_R = _WE.convModel.Reff[i][j];
			N = _R.length;
			for(k=0;k<N;k++) _R[k] -= _R0[k];
		  }
	  }
  }


  void logCriticalReflections(int kind,double sth_sq,int i_stack,double t)
  {
	if (criticalReflectionErrorLoggingCount < MAX_CRITICAL_REFL_LOGGING_COUNT)
	  {
		switch (kind)
		  {
		  case 0:
			// Dix
			_WE.EH.warn("Dix equation, stack "+(i_stack+1)+
					" time "+t+" yields impossible angle (sin[theta]="+Math.sqrt(sth_sq)+
					"): check stack params! Setting refl coeff to 0! ");
			break;
		  case 1:
			// Zoeppritz
			_WE.EH.warn("Zoeppritz equation, stack "+(i_stack+1)+
						" time "+t+" yields critical reflection: check stack params! Setting refl coeff to 0! ");
			break;
		  }
		_WE.EH.warn("Suggestions: (1) --dump-angles to see angles. (2) Check stack specs, "
					+"(3) Check log data for fast layers - try also --eliminate-log-spikes ");
		criticalReflectionErrorLoggingCount++;
	  } 
	else if (criticalReflectionErrorLoggingCount == MAX_CRITICAL_REFL_LOGGING_COUNT)
	  {
		_WE.EH.warn("Critical reflections: suppressing further errors...");
		criticalReflectionErrorLoggingCount++;
	  }
  }


  /**
	 File dump of stack angles 
  */
  void writeOutAngles()
  {
	int i,j,k;
	PrintWriter out;
	StackAngleLoggingData _stackAngleLoggingData;
	// open ascii summary file
	for (k=0;k<nWells;k++)
	  {
		String filename = Util.fileBasename(_WE.wellInfo[k].logFilename).replaceAll("(?i)\\.las|\\.txt","");
		String fullFilename = "StackAngles."+filename+".txt";
		_WE.EH.println("Writing stack angles to "+fullFilename);
		try 
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(fullFilename)));
			out.printf("Initial segment incident angles for well %s%n%d%n",filename,(1+_WE.m.nStacks));
			out.printf("Segment-top-time%n");
			for (j=0;j<_WE.m.nStacks;j++)
			  {
				out.printf("stack_%d-incident_angle(degrees)%n",(j+1));
			  }
			_stackAngleLoggingData = stackAngleLoggingData[k];
			// skip last point
			int length  = _stackAngleLoggingData.times.size()-1;

			for (i=0;i<length;i++)
			  {
				out.printf("%10.6g ",_stackAngleLoggingData.times.get(i));
				for (j=0;j<_WE.m.nStacks;j++)
				  {
					out.printf("%10.6g ",_stackAngleLoggingData.theta[j].get(i));
				  }
				out.println();
			  }
			out.close();
			_WE.EH.println("Wrote data to: " + fullFilename);
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + fullFilename);
		  }
	  }
  }

  /**
	 Useful function in creating synthetic minicubes: add a little scaled noise to reflectivity
	 data before convolving
  */
  void addNoiseToReflectivity()
  {
	int i,j,k;
	double var=0,noise;
	double[] R;
	DoubleArrayList dList;
	if (Math.abs(_WE.syntheticNoiseLevel) > 0)
	  {

		// get average reflectivity variance over all wells, stacks
		for (k=0;k<nWells;k++)
		  {
			for (i=0; i<nStacks; i++)
			  {
				dList = new DoubleArrayList(_WE.convModel.Reff[k][i]);
				var += Descriptive.sampleVariance(dList,0);
				dList=null;
			  }
		  }
		var/=(nStacks*nWells);
		noise = Math.abs(_WE.syntheticNoiseLevel)*Math.sqrt(var);
		AbstractDistribution Ngen = new Normal(0,1,new MersenneTwister(_WE.useSeed ? _WE.seed : MersenneTwister.DEFAULT_SEED));
		for (k=0;k<nWells;k++)
		  {
			for (i=0; i<nStacks; i++)
			  {
				R = _WE.convModel.Reff[k][i];
				for(j=0;j<R.length;j++)
				  {
					R[j] += noise * Ngen.nextDouble();
				  }
			  }
		  }
	  }
  }



  /**
	 Run over all wells to check if time to depth ranges and range of log data is adequate
  */
  boolean checkTimeToDepthRanges()
  {
	int k,Nwarn,iStart,iEnd;
	TimeToDepthMap t2dmap;
	boolean isOk=true;
	for (k=0;k<nWells;k++)
	  {
		t2dmap = timeToDepthMaps[k];

		// Does checkshot and marker data embrace estimation interval [tStart,tEnd] ?

		Nwarn = t2dmap.t.length-1;
		iStart = getStraddlingLowerIndex(_WE.tStart,t2dmap.t);
		if (iStart<0)
		  {
			_WE.EH.warn("Checkshot and Marker data for well "+
							  _WE.wellInfo[k].logFilename+
							  " do not extend to early enough times to "+
							  "perform calculations on the time interval ["+_WE.tStart+","+_WE.tEnd+"]");
			isOk =false;
		  }
		if (iStart==Nwarn) 
		  {
			_WE.EH.warn("Checkshot and Marker data for well "+
					_WE.wellInfo[k].logFilename+
					" do not extend to late enough times to "+
					"perform calculations on the time interval on ["+_WE.tStart+","+_WE.tEnd+"]");
			isOk =false;
		  }

		iEnd   = getStraddlingLowerIndex(_WE.tEnd,t2dmap.t);
		if (iEnd<0)
		  {
			_WE.EH.warn("Checkshot and Marker data for well "+
					_WE.wellInfo[k].logFilename+
					" do not extend to early enough times to "+
					"estimate wavelet on ["+_WE.tStart+","+_WE.tEnd+"]");
			isOk =false;
		  }
		if (iEnd==Nwarn)
		  {
			_WE.EH.warn("Checkshot and Marker data for well "+
					_WE.wellInfo[k].logFilename+
					" do not extend to late enough times to "+
					"estimate wavelet on ["+_WE.tStart+","+_WE.tEnd+"]");
			isOk =false;
		  }
	  }
	return isOk;
  }

  /**
	 Run over all wells to check if range of log data is adequate. Should be
	 called after checkTimeToDepthRanges() to verify sensibility of current data and
	 mapping configuration
  */
  boolean checklogDataRanges()
  {
	int i,k,indx,indxRec,Nwarn,iStart,iEnd;
	TimeToDepthMap t2dmap;
	double z, depths[];
	boolean isOk=true;

	for (k=0;k<nWells;k++)
	  {
		depths = _WE.sLogData[k].Z;
		Nwarn = depths.length-1;
		t2dmap = timeToDepthMaps[k];
		// we assume her that the next two calls don't fail.
		// Should be ok if checkTimeToDepthRanges() has been called before
		iStart = getStraddlingLowerIndex(_WE.tStart,t2dmap.t);
		iEnd   = 1+getStraddlingLowerIndex(_WE.tEnd,t2dmap.t);
		// loop over endpoint segments in time depth map
		for (i=0;i<=1;i++)
		  {
			// Major bracket left & right edge depth				
			z = t2dmap.Z[(i==0) ? iStart : iEnd];
			indx = getStraddlingLowerIndex(z,depths);
			if (indx < 0) 
			  {
				indxRec = getStraddlingLowerIndex(depths[0],t2dmap.Z);
				_WE.EH.warn("Log data in well "+
							_WE.wellInfo[k].logFilename+
							" do not extend to early enough depths to "+
							"estimate wavelet on ["+_WE.tStart+","+_WE.tEnd+"]. "+
							"Minimim time is about " + t2dmap.t[indxRec+1]);
				isOk=false;
			  }
			if (indx == Nwarn) 
			  {
				indxRec = getStraddlingLowerIndex(depths[Nwarn],t2dmap.Z);
				_WE.EH.warn("Log data in well "+
							_WE.wellInfo[k].logFilename+
							" do not extend to large enough depths to "+
							"estimate wavelet on ["+_WE.tStart+","+_WE.tEnd+"]. "+
							"Maximum time is about " + t2dmap.t[indxRec]);
				isOk=false;
			  }
		  }
	  }
	return isOk;
  }

  /**
	 Write checkshot interval velocities to file
  */
  void computeCheckshotIntervalVelocities()
  {
	int i,j,k;
	double vP;
	TimeToDepthMap t2dmap=null;

	for (k=0;k<nWells;k++)
	  {
		t2dmap = timeToDepthMaps[k];

		for (i=0; i<t2dmap.t.length-1; i++)
		  {
			vP = (t2dmap.Z[i+1]-t2dmap.Z[i])/(t2dmap.t[i+1]-t2dmap.t[i]);
			//_WE.EH.println("CSIV: "+i+" "+vP);
		  }
	  }
  }

 /**
	Assemble list of property titles in layer-model like String list.
	Make sure this is consistent with
	writeOutStochasticProperties(). This gets called once, not long
	after parsing all the files, so propertyNames can then be used as
	a reference in writing out the properties to the SU files.
  */
  String[] getPropertiesStringList()
  {
	int i,j,synchIndx;
	ArrayList strings = new ArrayList();

	// (2) checkshot and marker Vp's, t's
	for (j=0;j<nWells;j++)
	  {
		strings.add(_WE.wellInfo[j].well_name+"_checkshot_and_marker_interval_top_time");
		strings.add(_WE.wellInfo[j].well_name+"_checkshot_and_marker_interval_Vp");
	  }

	// (3) Well log errors ....??

	// (4)  time registration error for different stacks
	for (i=0; i<_WE.m.nWells; i++)
	  {
		// positioning errors may come from a master...get these
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;

		for (j=0; j<_WE.m.nStacks; j++)
		  {
			Property systematicTProp = _WE.seismicInfo[synchIndx][_WE.synchroniseStackTregistrationErrors ? 0 : j].systematicTProp;
			if (j==0 || !_WE.synchroniseStackTregistrationErrors)
			  {
				if (systematicTProp.sigma > 0)
				  {
					strings.add("t_shift_"+_WE.seismicInfo[i][j].filename);
				  }
			  }
		  }
	  }

	// (5) positioning error dx dy for each well
	for (i=0; i<_WE.m.nWells; i++)
	  {
		Property systematicXProp;
		Property systematicYProp;
		// positioning errors may come from a master...get these
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;
		systematicXProp = _WE.wellInfo[synchIndx].systematicXProp;
		systematicYProp = _WE.wellInfo[synchIndx].systematicYProp;

		if (systematicXProp.sigma > 0)
		  {
			strings.add("X_shift_"+_WE.wellInfo[i].well_name);
		  }
		if (systematicYProp.sigma > 0)
		  {
			strings.add("Y_shift_"+_WE.wellInfo[i].well_name);
		  }
	  }

	// (6) AVO's if modelled
	if (_WE.m.B_AVO[0].sigma>0)
	  {
		if (_WE.freeAVO)
		  {
			for (i=0; i<nStacks; i++)
			  {
				strings.add("B_AVO_stack_"+(i+1));
			  }
		  }
		else
		  {
			strings.add("B_AVO_stack_1");
		  }
	  }

	// 6.5 Anisotropy
	if (petro.hasAnisotropy) // _WE.doAnisotropicModel)
	  {
		if (_WE.m.deltaThomsen_shale.sigma > 0)
		  {
			strings.add("Thomsen_shale_delta");
		  }
	  }


	// (7) far stack wavelet scaling if modelled
	if (nStacks > 1 && _WE.useStretchedModel)
	  {
		for (j=1;j<nStacks;j++)
		  {
			if (_WE.m.farStackWaveletTimeStretchSigma >0)
			  {
				strings.add("stack_"+(j+1)+"_wavelet_t_stretch");
			  }
			if (_WE.m.farStackWaveletAmplitudeSigma >0)
			  {
				strings.add("stack_"+(j+1)+"_wavelet_amplitude");
			  }
		  }
	  }

	// (8) Noise levels:
	for (j=0;j<nStacks;j++)
	  {
		strings.add("stack_"+(j+1)+"_noise_rms");
	  }


	// save name in this class variable (it's used in writeOutStochasticProperties())
	propertyNames = (String[])(strings.toArray(new String[0]));

	return propertyNames;
  }


  /**
	 Write checkshot interval velocities to file
  */
  void writeOutStochasticProperties(int realisationNumber)
  {
	int i,ii,j,jMax=0,k,ns,nBlock,nT;
	double vP;
	TimeToDepthMap t2dmap=null;
	int synchIndx;
	Property systematicXProp;
	Property systematicYProp;


	if (true) return;	

	// Compute maximum Block size, which occurs for well jMax.
	// Let's make the size of the blocks equal to the largest number of checkshot *intervals*
	// by looking at all the checkshot+marker sets
	nBlock = 0;
	for (j=0;j<nWells;j++)
	  {
		if ((timeToDepthMaps[j].t.length-1) > nBlock)
		  {
			jMax = j;
			nBlock = timeToDepthMaps[j].t.length - 1;
		  }
	  }


	// compute trace size

	ns = propertyNames.length * nBlock;

	// and the last block must consist of (nBlock+1) values, as a surrogate for times
	ns += (nBlock+1);     


	// set header fields.
	_WE.tr.set("f1",0f);
	_WE.tr.set("ns",(short)ns);
	// make the header field go from 1...N
	_WE.tr.set("mark",(short)(realisationNumber+1));

	// now start writing stuff
	k=0;
	// (1,2) t2d times and velocities values
	for (j=0;j<nWells;j++)
	  {
		t2dmap = timeToDepthMaps[j];
		nT = t2dmap.t.length;
		// times (omit last time)
		for (i=0; i<(nT-1); i++)
		  {
			_WE.tr.set_data_value((float)t2dmap.t[i],k++);
		  }
		// pad if necessary
		for (i=0; i < (nBlock - (nT-1))  ; i++) _WE.tr.set_data_value(0f, k++);

		// velocities
		for (i=0; i<(nT-1); i++)
		  {
			vP = (t2dmap.Z[i+1]-t2dmap.Z[i])/(t2dmap.t[i+1]-t2dmap.t[i]);
			_WE.tr.set_data_value((float)vP,k++);
		  }
		// pad if necessary (extra one needed here!)
		for (i=0; i < (nBlock - (nT-1))  ; i++) _WE.tr.set_data_value(0f, k++);
	  }

	// (4)  time registration error for different stacks
	for (i=0; i<_WE.m.nWells; i++)
	  {
		// positioning errors may come from a master...get these
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;

		for (j=0; j<_WE.m.nStacks; j++)
		  {
			int stackIndx = _WE.synchroniseStackTregistrationErrors ? 0 : j;
			Property systematicTProp = _WE.seismicInfo[synchIndx][stackIndx].systematicTProp;
			if (systematicTProp.sigma > 0)
			  {
				for (ii=0; ii < nBlock ; ii++)
				  {
					_WE.tr.set_data_value((float)systematicTProp.value, k++);
				  }
			  }
		  }
	  }

	// (5) positioning error dx dy for each well
	for (i=0; i<_WE.m.nWells; i++)
	  {
		// positioning errors may come from a master...get these
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;
		systematicXProp = _WE.wellInfo[synchIndx].systematicXProp;
		systematicYProp = _WE.wellInfo[synchIndx].systematicYProp;

		if (systematicXProp.sigma > 0)
		  {
			for (ii=0; ii < nBlock ; ii++)
			  {
				_WE.tr.set_data_value((float)systematicXProp.value, k++);
			  }
		  }
		if (systematicYProp.sigma > 0)
		  {
			for (ii=0; ii < nBlock ; ii++)
			  {
				_WE.tr.set_data_value((float)systematicYProp.value, k++);
			  }	
		  }
	  }

	// (6) AVO's if modelled
	if (_WE.m.B_AVO[0].sigma>0)
	  {
		if (!_WE.freeAVO)
		  {
			for (i=0; i < nBlock ; i++)
			  {
				_WE.tr.set_data_value((float)_WE.m.B_AVO[0].value, k++);
			  }
		  }
		else
		  {
			for (j=0; j<nStacks; j++)
			  {
				for (i=0; i < nBlock ; i++)
				  {
					_WE.tr.set_data_value((float)_WE.m.B_AVO[j].value, k++);
				  }
			  }
		  }
	  }
	// free stack scale
	if (nStacks>1 && _WE.useStretchedModel)
	  {
		for (j=0; j<nStacks; j++)
		  {
			if (_WE.m.farStackWaveletTimeStretchSigma >0)
			  {
				double stretch;

				if (_WE.m.farStackWaveletTimeStretchMapper!=null) 
				  {
					stretch = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(_WE.m.farStackWaveletTimeStretch[j].value);
				  }
				else
				  {
					stretch = _WE.m.farStackWaveletTimeStretch[j].value;
				  }
				for (i=0; i < nBlock ; i++) _WE.tr.set_data_value((float)stretch, k++);
			  }
			if (_WE.m.farStackWaveletAmplitudeSigma >0)
			  {
				for (i=0; i < nBlock ; i++)
				  {	
					_WE.tr.set_data_value((float)_WE.m.farStackWaveletAmplitude[j].value, k++);
				  }
			  }
		  }
	  }

	// 6(b) anisotropy
	if (petro.hasAnisotropy) // _WE.doAnisotropicModel)
	  {
		if (_WE.m.deltaThomsen_shale.sigma >0)
		  {
			for (i=0; i < nBlock ; i++)
			  {	
				_WE.tr.set_data_value((float)_WE.m.deltaThomsen_shale.value, k++);
			  }
		  }
	  }


	// (7) Noise levels: duplicate across all the blocks
	for (j=0;j<nStacks;j++)
	  {
		for (ii=0; ii < nBlock ; ii++)
		  {
			_WE.tr.set_data_value((float)_WE.m.seismic_noise[j].value, k++);
		  }
	  }

	// The properties file must have a final block of (nBlock+1) values
	// we choose this to be the times of the longest checkshot
	t2dmap = timeToDepthMaps[jMax];
	for (ii=0; ii <= nBlock ; ii++)
	  {
		_WE.tr.set_data_value((float)t2dmap.t[ii], k++);
	  }
	  
	// Finally, write it
	_WE.suwriter.fputtr(_WE.out_strm,_WE.tr);
  }


  /**
	 Get index p such that xvalues[p] <= x < xValues[p+1]. xValues
	 must be already sorted into increasing order.  Returns -1 if
	 below range, xValues.length-1 if above range. Bisection used.
  */
  static int getStraddlingLowerIndex(double x, double xValues[])
  {
	int k,klo,khi;

	// bounding index
	if (x < xValues[0]) return -1;
	if (x > xValues[xValues.length-1]) return xValues.length-1;

	// bisection to find bracket
	klo = 0;
	khi = xValues.length-1;
	while (khi-klo > 1) 
	  {
		k=(khi+klo) >> 1;
		if (xValues[k] > x) khi=k;
		else klo=k;
	  }
	return klo;
  }

  static class StackAngleLoggingData
  {
	DoubleArrayList times;
	DoubleArrayList[] theta;

	StackAngleLoggingData(int nStacks)
	{
	  times = new DoubleArrayList();
	  theta = new DoubleArrayList[nStacks];
	  for (int i=0; i<nStacks;i++) theta[i] = new DoubleArrayList();
	}
  }
}
