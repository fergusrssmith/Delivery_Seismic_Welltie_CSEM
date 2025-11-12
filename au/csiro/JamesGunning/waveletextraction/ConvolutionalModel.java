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
import cern.colt.matrix.impl.*;

import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.Numerical.djbFFT.DJBFFT;
import au.csiro.JamesGunning.lib.util.*;

/**
   Basic convolution model variables and methods 
   for wavelet extraction.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class ConvolutionalModel
{
  /**
	 Flag to indicate Z,->t mapping. <i> Tests indicate that Splines
	 for the ZT mapping can cause serious flapping of the map if the
	 sequence of points is not especially smooth (as can happen under
	 possible configurations of the stochastic model): the usual
	 spline-kicking phenomenon. Linear mapping (which is used if this
	 is off) tends to be more robust under these circumstances.  </i>
  */
  public static final boolean useSplinesForZTmapping = false;

  /**
	 Flag to indicate rigorous exclusion of window edge regions
	 polluted by padding etc.  <i> Tests indicate that this is
	 necessary to prevent discontinuities in the objective function
	 caused by reflections wandering in and out of the
	 (_WE.tStart,_WE.tEnd) inversion interval.  These discontinuities
	 are large and seriously distress the minimiser.  </i>
	 <p>

	Excluding edge regions means also that the number of data points in the likelihood will vary
	as the wavelet span is changed, and the Bayesian machinery doesn't cope very well with
	this. For this reason, we set the range of data in the likelihood to be the <i>smallest</i>
	which will be seen by any of the wavelets
  */
  public static final boolean RIGOROUSLY_EXCLUDE_EDGE_REGIONS = true;

  /**
	 Set this to be true if you want the true (first trace) and synthetic
	 (second trace) SU files to include reflectivity in the third trace
	 as well
  */
  public final static boolean WRITE_REFLECTIVITIES_TO_TRUE_PLUS_SYNTH_FILES = true;
  /**
	 ptr to underlying wavelet extractor object and thus to all its fields
  */
  WaveletExtractor _WE;
  /**
	 ptr to underlying model
  */
  Model _m;

  BaseVintageReflectivities BVR=null;

  /**
	 Allowable start time of this wavelet extraction, for all wavelets
  */
  double tStartGlobal;

  /**
	 Allowable end time of this wavelet extraction, for all wavelets
  */
  double tEndGlobal;

  
  /**
	 Arrays of (optional) normal incidence reflectivities on equal spacing for each stack
	 Reff_P[i][k] = reflection for ith well, sample time k.
	 (Should not differ with stack I think)
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_P;
  /**
	 Arrays of (optional) shear reflectivities on equal spacing for each stack
	 (Should not differ with stack I think)
	 Reff_S[i][k] = shear reflection for ith well, sample time k.
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_S;
  /**
	 Arrays of (optional) density reflectivities on equal spacing for each stack
	 (Should not differ with stack I think)
	 Reff_S[i][k] = density reflection for ith well, sample time k.
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_Rho;

  /**
	 Arrays of reflectivities on equal spacing for each stack
	 Reff[i][j][k] = reflection for ith well, jth stack, sample time k
  */
  double[][][] Reff;

  /**
	 Arrays of synthetic seismic, 
	 S_synth[i][j][k] = reflection for ith well, jth stack, sample time k
  */
  double[][][] S_synth;
  /**
	 Number of distinct time samples used in forming synthetics
  */
  int nZ;
  /**
	 Snapshot arrays of true seismic, at current position of well, variable configuration etc 
	 S_true[i][j][k] = seismic for ith well, jth stack, sample time k
  */
  double[][][] S_true;
  /**
	 Snapshot arrays of errors (true seismic-synthetic), at current position of well, variable configuration etc 
	 errors[i][j][k] = error for ith well, jth stack, sample time k
  */
  double[][][] errors;

  /**
	 Snapshot arrays of index ranges of errors, at current position of well, variable configuration etc 
	 kStarts[i][j] = start index for ith well.
  */
  int[] kStarts;
  /**
	 Snapshot arrays of index ranges of errors, at current position of well, variable configuration etc 
	 kEnds[i][j] = end index for ith well.
  */
  int[] kEnds;

  /**
	 Snapshot cache arrays of amplitudes
	 cachedAmplitude[i][j][k] = seismic (data) amplitude for ith well, jth stack, sample time k
  */
  double[][][] cachedAmplitude;

  /**
	 Allow code to save trajectory values
  */
  public static final boolean SAVE_TRAJECTORY_VALUES = true;
  /**
	 Snapshot arrays of X position:
	 trajectoryX[i][k] = X coordinate for ith well, sample time k
  */
  double[][] trajectoryX;
  /**
	 Snapshot arrays of Y position:
	 trajectoryY[i][k] = Y coordinate for ith well, sample time k
  */
  double[][] trajectoryY;

  /**
	 Snapshot of current SS errors; a vector indexed by stack number
  */
  double[] SSerrors;
  
  
  /**
	 Starting time of range in which true reflections can occur
  */
  double tStartLegalRange;
 /**
	 End time of range in which true reflections can occur
  */
  double tEndLegalRange;
  /**
	 Starting time of reflectivity sequence
  */
  double tStartReff=0;

  /**
	 Final time of reflectivity sequence
  */
  double tEndReff=0;

  /**
	 Effective number of points in "data", for each stack
  */
  double[] dimData;

  /**
	 Last chisq computed in computeSyntheticNegLogLikelihood()
  */
  double chisq_last;


  private static final double TWOLOG2PI = 3.67575134;


  /**
	 Switch for miniCubeExceptionHandling
  */
  boolean softMiniCubeExceptionHandling = false;


  /**
	 switch for amplitude caching
  */
  public static boolean USE_CACHED_AMPLITUDES = true;

  /**
	 Penalty term to penalise minicube exceptions where the well
	 wanders outside the legal seismic minicube region
  */
  public static final double MINICUBE_EXCEPTION_PENALTY = 1000;

 

  public ConvolutionalModel(WaveletExtractor WE)
  {
	this._WE = WE;
	this._m = WE.m;

	// default values
	tStartGlobal = _WE.tStart;
	tEndGlobal = _WE.tEnd;

	if (!_WE.useJNI) DJBFFT.disableNativeFFT();
	if (_WE.verbosity>2)
	  {
		if (DJBFFT.haveNativeFFT()) 
		  {
			_WE.EH.println("Using JNI native FFT");
		  }
		else
		  {
			_WE.EH.println("Using non-optimised Java FFT");
		  }
	  }
  }
  /**
	 Set up required objects for convolutional modelling
  */
  void initialiseConvolutionalModel()
  {
	int i,N;
	double tStartEffective, tEndEffective;

	// Start synthetic seismic at integral multiple of dt just outside window.
	// possible extra time range (1 extra sample) needed for 4 pt Lagrange interpolation

	// If we assume zero reflectivity outside logs, we ought to start computing
	// reflectivity and synthetics outside the time window also:
	tStartEffective = _WE.tStart;
	tEndEffective = _WE.tEnd;
	if (_WE.assumeNoReflectionsOutsideLogs)
	  {
		tStartEffective = Double.MAX_VALUE;
		for (i=0;i<_WE.m.nWells;i++)
		  {
			tStartEffective = Math.min(tStartEffective,_WE.logData[i].tMin - _WE.tCodaMax);
		  }
		if (tStartEffective > _WE.tStart)
		  {
			// the start of window is already safe. Leave it as specified
			tStartEffective = _WE.tStart;
		  }
		else
		  {
			// start a bit earlier
			tStartEffective = Math.max(tStartEffective , _WE.tStart - _WE.tCodaMax);
		  }
		tEndEffective = -Double.MAX_VALUE;
		for (i=0;i<_WE.m.nWells;i++)
		  {
			tEndEffective = Math.max(tEndEffective,_WE.logData[i].tMax - _WE.tPrecursorMin);
		  }
		if (tEndEffective < _WE.tEnd)
		  {
			// the end of the window is already safe. Leave it as specified
			tEndEffective = _WE.tEnd;
		  }
		else
		  {
			// end a bit later
			tEndEffective = Math.min(tEndEffective , _WE.tEnd - _WE.tPrecursorMin);
		  }
	  }

	tStartLegalRange = (Math.floor(tStartEffective/_WE.dt))*_WE.dt;
	tStartReff = tStartLegalRange - _WE.dt;
	tEndLegalRange = tEndEffective;
	tEndReff = tEndLegalRange+ 2*_WE.dt;
	// for reflection coefficients, use this many samples....
	N = (int)Math.ceil((tEndReff-tStartReff)/_WE.dt);

	if (N<=4) _WE.EH.die_err("Too little data for noise estimation: N observations <=4");

	// set up stack arrays
	Reff = new double[_m.nWells][_m.nStacks][N];
	// synthetic seismic
	S_synth = new double[_m.nWells][_m.nStacks][N];

	// Set up true seismic snapshot arrays
	// true seismic: snapshot at current position also used as cached amplitudes
	S_true = new double[_m.nWells][_m.nStacks][N];
	// vectors of errors
	errors = new double[_m.nWells][_m.nStacks][N];
	// index bounds
	kStarts = new int[_m.nWells];
	kEnds = new int[_m.nWells];
	// well trajectory
	if (SAVE_TRAJECTORY_VALUES)
	  {
		trajectoryY = new double[_m.nWells][N];
		trajectoryX = new double[_m.nWells][N];
	  }
	// extra reflectivity storage
	if (_WE.dumpPandSreflectivity)
	  {
		Reff_P   = new double[_m.nWells][N];
		Reff_S   = new double[_m.nWells][N];
		Reff_Rho = new double[_m.nWells][N];
	  }

	// PS Keep nZ fixed after this!!
	nZ = N;

	// vector of SS errors for stacks
	SSerrors = new double[_m.nStacks];
	dimData = new double[_m.nStacks];
  }


  /**
	 compute a set of synthetic seismics from the reflection
	 coefficients precomputed in the Reff field. Computed using
	 current model wavelet, current anisotropy etc,
	 for all stacks and all wells.
  */
  void computeSyntheticSeismic()
  {
	int i,k,n;
	int pwStart,pwNumber;
	int wavPaddingLength, wFinalIndex, wLength;
	double[] R;         // working space ptr
	double[] padded_R;
	double[] padded_wavelet;
	double[] wavelet_coeffs;
	double[] padded_synth;
	Wavelet wavelet;

	// get wavelet from current model wavelet
	wavelet =_WE.m.wavelet;
	wavPaddingLength = getPaddinglength((!_WE.useStretchedModel) ? wavelet.values[0] : wavelet.values[1]);

	// Padding length of data in FFT's 
	// Protect for longer wavelets than data sequences
	n = Math.max(getPaddinglength(Reff[0][0]), wavPaddingLength);

	padded_R= new double[n];
	padded_wavelet= new double[n];
	padded_synth = new double[2*n];
	for (k=0;k<_m.nWells;k++)
	  {
		for (i=0; i<_m.nStacks; i++)
		  {
			// reflectivity sequence
			R = Reff[k][i];

			// generate convolution
			// Pad out R
			Arrays.fill(padded_R,0.0);
			Arrays.fill(padded_wavelet,0.0);
			Arrays.fill(padded_synth,0.0);
			Arrays.fill(S_synth[k][i],0.0);
			System.arraycopy(R,0,padded_R,0,R.length);

			// Pad out convolution
			if ((i==0) || !_WE.useStretchedModel)
			  {
				// normal wavelet
				wavelet_coeffs = wavelet.values[0];
				// +ve time part of wavelet
				wFinalIndex = wavelet.finalFineSample;
				wLength = wavelet.getLength();
			  }
			else
			  {
				// stretched wavelet
				wavelet_coeffs = wavelet.values[i];
				// +ve time part of wavelet
				wFinalIndex = wavelet.finalStretchedFineSample;
				wLength = wavelet.getStretchedLength();
			  }

			System.arraycopy(wavelet_coeffs,0,padded_wavelet,0,wFinalIndex+1);
			// -ve time part
			pwStart = wFinalIndex + 1 + n - wLength;
			pwNumber = wLength - 1 - wFinalIndex;
			System.arraycopy(wavelet_coeffs,wFinalIndex+1,
							 padded_wavelet, pwStart, pwNumber);
			if (_WE.m.debug) // Xsection debugging. Off for production
			  {
				Util.printfileZO("wavelet."+k+".txt",padded_wavelet,n);
			  }
			//Util.printfileZO("w.txt",padded_wavelet,padded_wavelet.length);

			// Get response
			DJBFFT.convolve(padded_R,n,padded_wavelet,padded_synth);
			// copy valid part of data to synthetics
			System.arraycopy(padded_synth,0,S_synth[k][i],0,Math.min(nZ,n));
			if (_WE.m.debug) // Xsection debugging. Off for production
			  {
				Util.printfileZO("Synth.well_"+k+".stack_"+i+"."+_WE.m.debug_cntr+".txt",
								 padded_synth,n);
			  }
		  }
	  }
  }

  /**
	 Switch on/off exception handling and X,Y limiting in computeSyntheticNegLogLikelihood()
  */
  void switchOnSoftMinicubeExceptionHandling(boolean on)
  {
	softMiniCubeExceptionHandling = on;
  }


// TIMING   long tStartTime; 
// TIMING   double tTime=0; 
// TIMING   double estimatedTime=0; 
// TIMING   double estimatedInnerTime=0; 
// TIMING   double tInnerTime=0; 

  /**
	 Method to compute seismic mismatch (-log(Likelihood)) after synthetic is computed
  */
  double computeSyntheticNegLogLikelihood(boolean workInLogNoiseParameterisation, boolean useCachedAmplitude)
  {
	int i,j,k;
	int inLineDenom=0,xLineDenom=0;	
	int iLeft=0,jLeft=0;
	int tIndx;
	boolean doFull3D=false;
	double tIndxDouble;
	double X,Y,Z,t,tRegOffset;
	double inLine,xLine;
	double amp=0,tmp1,tmp2,chisq,chiSqSum;
	double tAbove;
	double tSeis,wt, wt_first=0,A_neg1=0,A_0=1,A_1=0,A_2=0;
	double xWt,yWt;
	double seis00=0,seis01=0,seis10=0,seis11=0;
	double[] Ssynth =null;
	double[] Strue  =null;
	double[] err    =null;
	double[] wellTrajectoryX = null;
	double[] wellTrajectoryY = null;

	// boolean cacheAmpIsTheSameAsCurrentAmp=true; // FIXME DELETE LATER

	double	noiseStdDev, logNoiseVar;
	// DenseDoubleMatrix3D trueSeis; // JGDENSEHACK
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D trueSeis=null;
	UTM utmData=null;
	LinearInterpolation zInterp=null;
	Spline zSpline;
	Spline xSpline=null;
	Spline ySpline=null;
	UTM survey=null;
	double zMin=0,zMax=0;
	SeismicStackInfo seisInfo=null;
	UTMConversionInfo.UTMxyValues seisCoordsFromUTM=null;
 	Wavelet wavelet;
	double tStart,tEnd,samplingWeight;
	int kStart;
	int indxMax=0;
	int inlineMaxMinusOne=0, xlineMaxMinusOne=0;
	double miniCubeXYExceptionPenalty = 0.0;
	double miniCubeTExceptionPenalty = 0.0;

	int synchIndx;
	Property systematicXProp;
	Property systematicYProp;
	Property systematicTProp;

	// get wavelet from current model wavelet
	wavelet =_WE.m.wavelet;
	// sampling weight to compensate for oversampling
	samplingWeight = Math.sqrt(_WE.dt/wavelet.dtCoarse);

	// setup
	Arrays.fill(SSerrors,0);
	Arrays.fill(dimData,0.0);

	// sum over wells and stacks

// TIMING 	estimatedTime=0; 
// TIMING 	estimatedInnerTime=0; 
// TIMING 	long startInnerTime=0,tStartTotalTime;

	chisq = 0.0;
	for (i=0;i<_m.nWells;i++)
	  {
		// Time range might possibly vary with well
		if (RIGOROUSLY_EXCLUDE_EDGE_REGIONS)
		  {
			// Global region:
			tStart = tStartGlobal;
			tEnd = tEndGlobal;
		  }
		else
		  {
			tStart = _WE.tStart;
			tEnd = _WE.tEnd ;
		  }
		if (_WE.ALLOW_INCOMPLETE_LOGS)
		  {
			if (_WE.assumeNoReflectionsOutsideLogs)
			  {
				// Assume no significant reflections coming from the 
				// coda/precursor segments just outside the logs
				tStart = Math.max(tStart,_WE.logData[i].tMin);
				tEnd = Math.min(tEnd,_WE.logData[i].tMax);
			  }
			else
			  {
				// make sure we are safely inside the maximum log data range if we
				// can't assume absence of reflections outside of it
				tStart = Math.max(tStart,_WE.logData[i].tMin + _WE.multiModelSetUp.maxCodaMargin);
				tEnd = Math.min(tEnd,_WE.logData[i].tMax - _WE.multiModelSetUp.maxPrecursorMargin);
			  }
		  }




		// starting index wrt Reff array. Has to be internal in window
		kStart = (int)Math.ceil((tStart - tStartReff)/_WE.dt);

		// Well positioning errors may come from a master...get these.
		// If no master, it will come from itself (synchronizeWellIndex == i)
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;
		systematicXProp = _WE.wellInfo[synchIndx].systematicXProp;
		systematicYProp = _WE.wellInfo[synchIndx].systematicYProp;

		// trajectories
		if (SAVE_TRAJECTORY_VALUES)
		  {
			wellTrajectoryX = trajectoryX[i];
			wellTrajectoryY = trajectoryY[i];
		  }

		// set up splines for below if required
		TimeToDepthMap t2dMap  = _WE.m.timeToDepthMaps[i];
		// System.out.println("CACHEUSE " + ((useCachedAmplitude && USE_CACHED_AMPLITUDES)?"1":"0"));

		if (!(useCachedAmplitude && USE_CACHED_AMPLITUDES))
		  {
			if (useSplinesForZTmapping)
			  {
				zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
									 Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
			  }
			else
			  {
				zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);
			  }

			// the following ought to be done once and for all......
			utmData = _WE.logData[i].survey;
			xSpline = utmData.xSpline;
			ySpline = utmData.ySpline;
			zMin = utmData.zMin;
			zMax = utmData.zMax;
		  }

		for (j=0; j<_m.nStacks; j++)
		  {
   			// Registration errors may come from a master...get this.
			systematicTProp = _WE.seismicInfo[synchIndx][_WE.synchroniseStackTregistrationErrors ? 0 : j].systematicTProp;
			tRegOffset = systematicTProp.value;

			// synthetic starting at times tStartReff in steps of _WE.dt
			// is S_synth[i][j]
			Ssynth = S_synth[i][j];
			// true seismic snapshot variables
			Strue = S_true[i][j];
			// errors
			err = errors[i][j];
			// Arrays.fill(err,0.0);
			
			// noise level for this stack
			noiseStdDev = _WE.m.seismic_noise[j].value;

			// only need all this apparatus if true amplitude must be extracted afresh from minicube
			if (!(useCachedAmplitude && USE_CACHED_AMPLITUDES))
			  {
				// actual seismic data is in
				seisInfo = (_WE.seismicInfo[i][j]);
				trueSeis = seisInfo.data;
				indxMax = trueSeis.columns() - 3;
				inlineMaxMinusOne = seisInfo.inline_max - 1;
				xlineMaxMinusOne = seisInfo.xline_max - 1;
				inLineDenom = ((seisInfo.inlineSpacing!=0)? seisInfo.inlineSpacing : 1 );
				xLineDenom = ((seisInfo.xlineSpacing!=0)? seisInfo.xlineSpacing : 1 );
				doFull3D = (!seisInfo.is1DdataAndUseSnapping && !seisInfo.is2DdataAndUseSnapping);
				// get x,y positions within minicube by a spline of x(t) y(t)
				// generated by current timetodepth map
				// zSpline or zInterp maps t->Z
			  }

			// starting indx and (true) time
			// start of window region
			// NB k will start close to zero, but may not be exactly zero!
			t = tStartReff + _WE.dt * (double)kStart;
			chiSqSum = 0.0;
			for (k=kStart; t<=tEnd ; k++, t+=_WE.dt)
			  {
				if (useCachedAmplitude && USE_CACHED_AMPLITUDES)
				  {
					amp = Strue[k];
				  }
				else
				  {
// TIMING 					tStartTotalTime = System.nanoTime(); 

					// Extract the amplitude afresh from the minicube by computed the well position etc etc
					// get true depth at this time
					if (useSplinesForZTmapping) 
					  {
						Z = zSpline.splint(t,k!=kStart);
					  }
					else
					  {
						Z = zInterp.interpolate(t,k!=kStart);
					  }
					// and therefore current well position
					if (Z<=zMin)
					  {
						X = xSpline.splint(zMin,false);
						Y = ySpline.splint(zMin,false);
					  }
					else if (Z>=zMax)
					  {
						X = xSpline.splint(zMax,false);
						Y = ySpline.splint(zMax,false);
					  }
					else
					  {
						// we are here 99% of the time...
						if (utmData.wellIsVertical)
						  {
							X = utmData.X[0];
							Y = utmData.Y[0];
						  }
						else
						  {
							X = xSpline.splint(Z,k!=kStart);
							Y = ySpline.splint(Z,k!=kStart);
						  }
					  }

					// But there may be a systematic error in this...
					X += systematicXProp.value;
					Y += systematicYProp.value;

					// convert to inlines/xlines. Three styles possible
					if (seisInfo.is1DdataAndUseSnapping)
					  {
						// Data set is single trace only
						inLine = seisInfo.inline_min; // SNAPPED (seisCoordsFromUTM.x);
						xLine = seisInfo.xline_min;   // SNAPPED (seisCoordsFromUTM.y);

					  }
					else if (seisInfo.is2DdataAndUseSnapping)
					  {
						// Data set is 2D line. Use snapping to project onto inline/x axis: xline=xline_min
						seisCoordsFromUTM = seisInfo.utmInfo.getInlineXlineFromUTMS(X,Y);
						if (seisInfo.fastDirection == 'i')
						  {
							inLine = seisCoordsFromUTM.x;
							xLine = seisInfo.xline_min;    // SNAPPED (seisCoordsFromUTM.y);
						  }
						else
						  {
							inLine = seisInfo.inline_min;  // SNAPPED (seisCoordsFromUTM.x);
							xLine = seisCoordsFromUTM.y;
						  }
					  }
					else if (!_WE.useUTMconversion)
					  {
						// use 3D header words
						seisInfo.computeIndices(X,Y);
						inLine = seisInfo.getInLine();
						xLine = seisInfo.getXLine();
					  }
					else
					  {
						// use global UTM's
						seisCoordsFromUTM = _WE.utmInfo.getInlineXlineFromUTMS(X,Y);
						inLine = seisCoordsFromUTM.x;
						xLine = seisCoordsFromUTM.y;
					  }
					// save trajectory if needed - same for all stacks
					if (SAVE_TRAJECTORY_VALUES && j==0)
					  {
						wellTrajectoryX[k] = inLine;
						wellTrajectoryY[k] = xLine;
					  }



					// then get amplitude at this (inLine,xLine,t) from minicube using bilinear interpolation
					// If 1D, everything is safe.Skip this.
					if (!seisInfo.is1DdataAndUseSnapping)
					  {
						if (seisInfo.is2DdataAndUseSnapping)
						  {
							if (softMiniCubeExceptionHandling)
							  {
								if (seisInfo.fastDirection == 'i')
								  {
									// edge-limited safe handling. xline is always safe
									// Put position on edge of 2Dline, but add in a penalty
									if (inLine < seisInfo.inline_min)
									  {
										miniCubeXYExceptionPenalty = seisInfo.inline_min - inLine;
										inLine = seisInfo.inline_min;
									  }
									if (inLine > inlineMaxMinusOne)
									  {
										miniCubeXYExceptionPenalty = inLine - seisInfo.inline_max;
										inLine = inlineMaxMinusOne;
									  }
								  }
								else
								  {
									// edge-limited safe handling. inline is always safe
									// Put position on edge of 2Dline, but add in a penalty
									if (xLine < seisInfo.xline_min)
									  {
										miniCubeXYExceptionPenalty = seisInfo.xline_min - xLine;
										xLine = seisInfo.xline_min;
									  }
									if (xLine > xlineMaxMinusOne)
									  {
										miniCubeXYExceptionPenalty = xLine - seisInfo.xline_max;
										xLine = xlineMaxMinusOne;
									  }
								  }
							  }
							else
							  {					
								// catastrophic handling
								if (seisInfo.fastDirection == 'i')
								  {
									if (inLine < seisInfo.inline_min)
									  {
										throw new RuntimeException("Well position "+inLine+" < than inLine range of 2D line: "+seisInfo.inline_min);
									  }
									if (inLine > inlineMaxMinusOne)
									  {
										throw new RuntimeException("Well position "+inLine+" > than inLine range of 2D line: "+seisInfo.inline_max);
									  }
								  }
								else
								  {
									if (xLine < seisInfo.xline_min)
									  {	
										throw new RuntimeException("Well position "+xLine+" < than xLine range of 2D line: "+seisInfo.xline_min); 
									  }
									if (xLine > xlineMaxMinusOne)
									  {
										throw new RuntimeException("Well position "+xLine+" > than xLine range of 2D line: "+seisInfo.xline_max); 
									  }
								  }
							  }
						  }
						else
						  {
							// usual full 3D interpolation. Check inLine, xLine
							if (softMiniCubeExceptionHandling)
							  {
								// edge-limited safe handling
								// Put position on edge of minicube, but add in a penalty
								if (inLine < seisInfo.inline_min)
								  {
									miniCubeXYExceptionPenalty = seisInfo.inline_min - inLine;
									inLine = seisInfo.inline_min;
								  }
								if (inLine > inlineMaxMinusOne)
								  {
									miniCubeXYExceptionPenalty = inLine - seisInfo.inline_max;
									inLine = inlineMaxMinusOne;
								  }
								if (xLine < seisInfo.xline_min)
								  {
									miniCubeXYExceptionPenalty = seisInfo.xline_min - xLine;
									xLine = seisInfo.xline_min;
								  }
								if (xLine > xlineMaxMinusOne)
								  {
									miniCubeXYExceptionPenalty = xLine - seisInfo.xline_max;
									xLine = xlineMaxMinusOne;
								  }
							  }
							else
							  {					
								// catastrophic handling
								if (inLine < seisInfo.inline_min)
								  {
									throw new RuntimeException("Well position "+inLine+" < than inLine range of minicube: "+seisInfo.inline_min);
								  }
								if (inLine > inlineMaxMinusOne)
								  {
									throw new RuntimeException("Well position "+inLine+" > than inLine range of minicube: "+seisInfo.inline_max);
								  }
								if (xLine < seisInfo.xline_min)
								  {
									throw new RuntimeException("Well position "+xLine+" < than xLine range of minicube: "+seisInfo.xline_min);
								  }
								if (xLine > xlineMaxMinusOne)
								  {
									throw new RuntimeException("Well position "+xLine+" > than xLine range of minicube: "+seisInfo.xline_max);
								  }
							  }
						  }
					  }


					// Perform interpolation by 4 pt Lagrange on the 4 bordering traces to get 4 neighbours at the
					// correct time, then bilinear in the plane

				
					// Compute seismic time: may differ by a registration error from true time
					tSeis = t + tRegOffset;
					// And pull out seismic data value at this time....
					// time index wrt true Seismic cube
					tIndxDouble = (tSeis - seisInfo.tStartData)/_WE.dt;
					// tIndx = (int)Math.floor(tIndxDouble);
					tIndx = (int)tIndxDouble;
					//_WE.EH.println("OFFSETS "+	tIndxDouble+" "+tIndx);				
					if (softMiniCubeExceptionHandling)
					  {
						// edge-limited safe handling
						// Put position on edge of minicube, but add in a penalty (which is *smooth* in t)
						if (tIndx<1)
						  {
							// clamp it at the top of the seismic
							miniCubeTExceptionPenalty += (1 - tIndxDouble);
							tIndx = 1;
						  }
						if (tIndx > indxMax)
						  {
							// clamp it at the bottom of the seismic
							miniCubeTExceptionPenalty += (tIndxDouble - indxMax);
							tIndx = indxMax;
						  }
					  }
					else
					  {
						// fussy handling
						if (tIndx<1)
						  {
							throw new RuntimeException("time "+ tSeis+" required in seismic above minicube");
						  }
						if (tIndx > indxMax)
						  {
							throw new RuntimeException("time "+ tSeis+" required in seismic below minicube");
						  }
					  }

					// Compute the four coefficients A_{-1},A_0,A_1,A_2
					// as per p 879 of Abramowitz & Stegun....
					// NB 1/6 = 0.16666666666666666666;

					// TIMING startInnerTime = tStartTotalTime; // System.nanoTime(); 

					tAbove = tIndx * _WE.dt + seisInfo.tStartData;
					wt = (tSeis - tAbove)/_WE.dt;
					// DON'T DO THIS: wt can flip between 0 and 0.999999999999
					// because of rounding, so take care!!!!! Recompute weights each time!
					// if ((k == kStart)) // tSeis goes up in steps of dt, so weights are always same as first...
					  {
						A_neg1 = -wt*(wt-1.0)*((tmp1 = wt-2.0))*0.166666666666666;
						A_0 = ((tmp2 = wt*wt-1.0))*tmp1*0.5;
						A_1 = -wt*(wt+1.0)*tmp1*0.5;
						A_2 = wt*tmp2*0.166666666666666;
						wt_first=wt;
					  }
					// lower left corner indices of 4 traces surrounding point (inLine,xLine,t)....
					// iLeft = (int)Math.floor((inLine - seisInfo.inline_min)/inLineDenom);
					// jLeft = (int)Math.floor((xLine - seisInfo.xline_min)/xLineDenom);
					// faster: use floor behaviour of (int) for positive values...
					iLeft = (int)((inLine - seisInfo.inline_min)/inLineDenom);
					jLeft = (int)((xLine - seisInfo.xline_min)/xLineDenom);

					if (doFull3D)
					  {
						// full 3D case: most common case
						// Four corner values required for Bilinear from Lagrange
						// bilinear weights
						xWt = (inLine - (seisInfo.inline_min + iLeft*seisInfo.inlineSpacing))/(double)inLineDenom;
						yWt = (xLine - (seisInfo.xline_min + jLeft*seisInfo.xlineSpacing))/(double)xLineDenom;
						try
						  {
							if (true)
							  {
								amp =  trueSeis.interpolate3DQuick(iLeft,jLeft,tIndx,wt_first,xWt,yWt);
							  }
							else
							  {
								seis00 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
								  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
								iLeft++;
								
								seis10 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
								  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
								jLeft++;
								
								seis11 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
								  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
								iLeft--;
								
								seis01 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
								  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
								jLeft--;

								// Bilinear formula for amplitude
								amp = (1-xWt)*(1-yWt)*seis00 + xWt*(1-yWt)*seis10 + xWt*yWt*seis11 + (1-xWt)*yWt*seis01;
							  }
						  }
						catch (ArrayIndexOutOfBoundsException e)
						  {
							_WE.EH.println("i ="+iLeft+" j="+jLeft+" k="+tIndx+" cols ="+trueSeis.columns());
						  }
					  }
					else if (seisInfo.is1DdataAndUseSnapping)
					  {
						// 1D case. Linear interp. 
						// Snapping projects everything onto the 00 point, so iLeft=jLeft=0, xWt=yWt=0
						seis00 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
						  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
						amp = seis00;
					  }
					else if (seisInfo.is2DdataAndUseSnapping)
					  {
						if (seisInfo.fastDirection == 'i')
						  {
							// 2D case. Snapping and linear interp. There are only inlines. xLine == seisInfo.xline_min always
							// Snapping projects everything onto the inline axis, so jLeft=0, yWt=0
							seis00 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
							  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
							iLeft++;
							seis10 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
							  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
							iLeft--;
							// Linear interp formula for amplitude
							xWt = (inLine - (seisInfo.inline_min + iLeft*seisInfo.inlineSpacing))/(double)inLineDenom;
							amp = (1-xWt)*seis00 + xWt*seis10;
						  }
						else
						  {
							// 2D case. Snapping and linear interp. There are only xlines. inLine == seisInfo.inline_min always
							// Snapping projects everything onto the xline axis, so iLeft=0, xWt=0
							seis00 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
							  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
							jLeft++;
							seis01 =  A_neg1 * trueSeis.getQuick(iLeft,jLeft,tIndx-1) + A_0* trueSeis.getQuick(iLeft,jLeft,tIndx)
							  + A_1 * trueSeis.getQuick(iLeft,jLeft,tIndx+1) + A_2* trueSeis.getQuick(iLeft,jLeft,tIndx+2);
							jLeft--;
							// Linear interp formula for amplitude
							xWt = (xLine - (seisInfo.xline_min + jLeft*seisInfo.xlineSpacing))/(double)xLineDenom;
							amp = (1-xWt)*seis00 + xWt*seis01;
						  }
					  }
// TIMING 					xWt = System.nanoTime();  	
// TIMING 					estimatedTime += (xWt - tStartTotalTime); 					
// TIMING 					estimatedInnerTime += (xWt- startInnerTime);  


					// Save useful amplitude for later. Used in both caching and I/O
					Strue[k] = amp;
				  }
				// OK Now we've got the amplitude !!

				// compute & store the normalised error, downweighted by the sampling weight
				err[k] = samplingWeight*(Ssynth[k]-amp)/noiseStdDev;

				// if we're ignoring seismic, turn all this off
				if (_WE.ignoreSeismic) err[k]=0.0;

				// now add up the chisq of the mismatch for this stack
				chiSqSum += err[k]*err[k];
				// SSerrors[j] += err[k]*err[k];


				// Caching test: comment out above line to apply
				/*
				  if (USE_CACHED_AMPLITUDES) 
				    {
				      if (amp!=Strue[k]) 
					    {
						  cacheAmpIsTheSameAsCurrentAmp = false;
    					}
					  Strue[k] = amp;
				    }
				*/
			  } // sum over samples (k)
			// remember array bounds
			if (j==0)
			  {
				kStarts[i] = kStart;
				kEnds[i] = k-1;
			  }
			// Cache tests
			/*
			  if (USE_CACHED_AMPLITUDES) // FIXME DELETE LATER
			  {
				if (useCachedAmplitude && !cacheAmpIsTheSameAsCurrentAmp)
				  {
					// bad prediction !!!
					_WE.EH.println("BAD CACHE "+(cacheAmpIsTheSameAsCurrentAmp?1:0)+(useCachedAmplitude?1:0)+
								   ((miniCubeXYExceptionPenalty!=0)?1:0)+((miniCubeTExceptionPenalty!=0)?1:0)
								   +" " +maxCacheErr); // FIXME DELETE LATER
				  }
				else
				  {
					// good
					_WE.EH.println("CACHE "+(cacheAmpIsTheSameAsCurrentAmp?1:0)+(useCachedAmplitude?1:0)+
								   ((miniCubeXYExceptionPenalty!=0)?1:0)+((miniCubeTExceptionPenalty!=0)?1:0)
								   +" " +maxCacheErr); // FIXME DELETE LATER
				  }
			  }
			*/

			chisq += chiSqSum;               // total running sum
			SSerrors[j] += chiSqSum;         // total running sum for jth stack

			// Make SSerrors the raw (unnormalised) error measurement
			// SSerrors[j] *= noiseStdDev*noiseStdDev;
			// number of "effective" data points used: summed over all wells
			// total running data count for jth stack
			dimData[j] += (samplingWeight*samplingWeight*(k-kStart));


			// And finally the piece from the determinant of Cnoise
			// the dimensionality of the "data" is taken to be dimData at this point
			if (!_WE.ignoreSeismic) 
			  {
				//_WE.EH.printf("ChsqSeis (SS)="+chisq+" ");
				// The determinant is ndata * log(sigma^2)
				// Take ndata = (samplingWeight*samplingWeight*((k-kStart)
				// introduce fixed normalisation constant here to prevent this term swamping the SSE term
				tmp1 = noiseStdDev/_WE.multiModelSetUp.typicalRMSSeismicAmplitude;
				logNoiseVar = Math.log(tmp1*tmp1);
				double nDataEff = samplingWeight*samplingWeight*((k-kStart));
				chisq += nDataEff * logNoiseVar;
				// normalisation may vary with dimension of data
				chisq += nDataEff * TWOLOG2PI;  
				//_WE.EH.println(" (+det) "+chisq);
			  }
		  } // sum over stacks
	  }     // sum over wells

// TIMING 	tTime += 1.0e-6*estimatedTime; 
// TIMING 	tInnerTime += 1.0e-6*estimatedInnerTime; 

	// Make SSerrors the raw (unnormalised) error measurement
	for (j=0; j<_m.nStacks; j++)
	  {
		// noise level for this stack
		noiseStdDev = _WE.m.seismic_noise[j].value;
		SSerrors[j] *= noiseStdDev*noiseStdDev;
	  }

	// Add in out-of-bounds penalties (usually 0)
	chisq += (miniCubeXYExceptionPenalty>0)? miniCubeXYExceptionPenalty*miniCubeXYExceptionPenalty : 0;
	chisq += (miniCubeTExceptionPenalty>0)? miniCubeTExceptionPenalty*miniCubeTExceptionPenalty : 0;

	if (miniCubeTExceptionPenalty>0 )
	  {
		_WE.EH.println("warning: minicube t penalty "+ miniCubeTExceptionPenalty);
	  }
	if (miniCubeXYExceptionPenalty>0 )
	  {
		_WE.EH.println("warning: minicube xy penalty "+ miniCubeXYExceptionPenalty);
	  }

	chisq/=2.0;
	
	// JEFFREY's term: one per stack-noise parameter
	// Jeffrey's prior term is 1/|sigma|, so the neg-log
	// equivalent is log(sigma) = 0.5 log(sigma^2)
	if (!_WE.fixNoise)
	  {
		if (_WE.ignoreSeismic) 
		  {
			// desire noise to have a Gaussian posterior N(1,1)
			// we should have _WE.m.seismic_noise[j].mean = 1
			for (j=0; j<_m.nStacks; j++)
			  {
				// Broad Gaussian form exactly as per completeModeCalculation() function
				// mean is taken to be modelVariables.initialGuess[i], before working in log space
				if (workInLogNoiseParameterisation)
				  {
					wt = Math.log(_WE.m.seismic_noise[j].value) - Math.log(_WE.m.seismic_noise[j].mean);
				  }
				else
				  {
					wt = _WE.m.seismic_noise[j].value - _WE.m.seismic_noise[j].mean;
				  }
				chisq += 0.5*wt*wt;
			  }
		  }
		else
		  {
			if (ModeFinder.ABSORB_JEFFREYS_PRIOR_INTO_LIKELIHOOD)
			  {
				if (!workInLogNoiseParameterisation)
				  {
					// If working in log-noise space, the Jacobian has removed this term
					for (j=0; j<_m.nStacks; j++)
					  {
						logNoiseVar = Math.log(_WE.m.seismic_noise[j].value * _WE.m.seismic_noise[j].value);
						chisq += 0.5*logNoiseVar;
					  }
				  }

				// Weak prior term on log noise: only operates if -IS is switched on
				if (false && workInLogNoiseParameterisation) 
				  {
					for (j=0; j<_m.nStacks; j++)
					  {
						// Broad Gaussian form exactly as per completeModeCalculation() function
						// mean is taken to be modelVariables.initialGuess[i], before working in log space
						wt = Math.log(_WE.m.seismic_noise[j].value) - Math.log(_WE.m.seismic_noise[j].mean);
						chisq += 0.5*wt*wt*ModeFinder.PRIOR_INVERSE_VARIANCE_FOR_LOG_NOISE;
					  }
				  }
			  }
		  }
	  }


	chisq_last = chisq;

	//if (_WE.verbosity>5) _WE.EH.println("seismic chisq = "+chisq);

	return chisq;
  }

  // Instrumentation for timings...
  void writeOptimiserTimings()
  {
// TIMING 	_WE.EH.println("TOTAL TIMES: "+ tTime*.001+ " inner: "+tInnerTime*0.001); 
  }

  /**
	 Here we compute (as per the linear estimate) the multi-stack noise covariance: see the
	 appendix in the WaveletExtractor paper
  */
  double[][] computeMultiStackNoiseCovarianceMatrix()
  {
	int i,j,jp,k;
	double sum, noiseStdDev_j, noiseStdDev_jp;
	double[] SSerrorsCovSqrtDiag;
	double[][] err_i;
	double[][] SSerrorsCov;

	SSerrorsCov = new double[_m.nStacks][_m.nStacks];
	SSerrorsCovSqrtDiag = new double[_m.nStacks];
	for (j=0; j<_m.nStacks; j++)
	  {
		noiseStdDev_j = _WE.m.seismic_noise[j].value;
		for (jp=0; jp<_m.nStacks; jp++)
		  {
			noiseStdDev_jp = _WE.m.seismic_noise[jp].value;
			sum=0;
			for (i=0;i<_m.nWells;i++)
			  {
				err_i = errors[i];
				// err_i = S_synth[i];  // can be used to check correlation of fwd models b/w stacks - often about 0.85ish
				for (k=kStarts[i]; k<=kEnds[i]; k++)
				  {
					// err[k] = samplingWeight*(Ssynth[k]-amp)/noiseStdDev;
					sum += err_i[j][k] * err_i[jp][k];
				  }
			  }
			sum *= (noiseStdDev_j * noiseStdDev_jp);
			// Assume here dimData[j] = dimData[jp]. Pretty sure this is always true
			SSerrorsCov[j][jp] = sum/dimData[j]; 
		  } // stack jp
	  }     // stack j

	// Convert to corr. coeff form
	for (j=0; j<_m.nStacks; j++) SSerrorsCovSqrtDiag[j] = Math.sqrt(SSerrorsCov[j][j]);

	for (j=0; j<_m.nStacks; j++)
	  {
		SSerrorsCov[j][j] = 1.0;
		for (jp=j+1; jp<_m.nStacks; jp++)
		  {
			SSerrorsCov[j][jp] /= (SSerrorsCovSqrtDiag[j]*SSerrorsCovSqrtDiag[jp]);
			SSerrorsCov[jp][j] = SSerrorsCov[j][jp];
		  }
	  }
	if (_WE.verbosity>3) 
	  {
		_WE.EH.printf("Noise level correlation coefficient matrix:%n");	
		for (j=0; j<_m.nStacks; j++)
		  {
			for (jp=0; jp<_m.nStacks; jp++)
			  {
				_WE.EH.printf("%-10.6g ",SSerrorsCov[j][jp]);
			  }
			_WE.EH.printf("%n");
		  }
	  }
	SSerrorsCovSqrtDiag=null;

	return SSerrorsCov;
  }



  /**
	 Compute rms of true seismic at wells from most recent cached value, ignoring zeros.
	 Averaged over all wells and stacks.
  */
  double[] getTrueSeisRMSatWells()
  {
	int i,j,k,cnt=0;
	double[] _data;
	double sum=0.0;
	double[] rms = new double[_WE.m.nStacks];

	for (j=0; j<_WE.m.nStacks; j++)
	  {
		sum=0.0;
		for (i=0;i<_WE.m.nWells;i++)
		  {
			// Ignore zeros as they are almost certainly from the padding at missing log data
			_data = _WE.convModel.S_true[i][j];
			for (k=0;k<_data.length;k++) 
			  {
				if (_data[k]!=0.0) 
				  {
					sum+=_data[k]*_data[k];
					cnt++;
				  }
			  }
		  }
		rms[j] = sum/((cnt>0)?cnt:1.0);
		rms[j] = Math.sqrt(rms[j]);
	  }

	return rms;
  }

  /**
	 <p>
	 Write out 
	 SU versions of current synthetics. Designed to be called after
	 computeSyntheticNegLogLikelihood() has been called with sane
	 variables in place. The files written will have true and synthetic 
	 traces written to adjacent traces: </p>
	 <p>
	 Traces appear in this order:
	 </p>
	 (1) true seismic <br>
	 (2) synthetic seismic <br>
	 (3) (maybe) reflectivity<br>
	 (4) (maybe) X values <br>
	 (5) (maybe) Y values <br>
	 (6) (maybe) direct (normal incidence) P reflectivity
	 (7) (maybe) direct S reflectivity
	 <p>
	 There are distinct files for each
	 stack & each well.
  */
  void writeOutTrueAndSyntheticSUset(String prefix)
  {
	int i,j,k,n,synchIndx;
	double MLoffset_val=0.0,x=0,y=0;
	double[] S;
	String filename=null;
	FileOutputStream out=null;
	Property systematicTProp;

	if (_WE.dumpPandSreflectivity)
	  {
		if (_WE.useSaturationModel) _WE.m.recomputeLogProperties(true);
		// compute P and S reflectivities
		// FIXME WRONG! for saturations
		_WE.m.computeReflectivities(false,true);
	  }

	n = S_true[0][0].length;


	// set requisite fields
	// first, clear all other stuff that might be distracting & wrong
	_WE.tr.zeroHeader();

	// default f1: overridden with cube offset later if possible
	//_WE.tr.set("f1",(float)tStartReff);
	_WE.tr.set("delrt",(short)(Math.rint(1e3*tStartReff))); // in ms. Best approx to tStartReff
	_WE.tr.set("ns",(short)n);	
	_WE.tr.set("dt",_WE.m.wavelet.dt);

	// make traces for each well
	for (i=0;i<_m.nWells;i++)
	  {
		synchIndx = _WE.wellInfo[i].synchronizeWellIndex;
		

		// a few useful header words perhaps
		if (SAVE_TRAJECTORY_VALUES)
		  {
			// find first non-zero inline crossline
			S = trajectoryX[i];
			for (k=0;k<nZ;k++)
			  {
				if ((x=S[k]) != 0.0)
				  {
					y = trajectoryY[i][k];
					break;
				  }
			  }
			_WE.tr.set(_WE.seismicInlineWord,(int)Math.rint(x));
			_WE.tr.set(_WE.seismicXlineWord,(int)Math.rint(y));
		  }

		// and for each stack
		for (j=0; j<_m.nStacks; j++)
		  {
			// suitable file for data
			filename = (prefix+"."+(Util.fileBasename(_WE.wellInfo[i].logFilename)).replaceAll("(?i)\\.las|\\.txt","")
						+".stack."+(j+1)+".su");
			_WE.EH.println("Writing synthetic and true seismic to "+filename);

	
			systematicTProp = _WE.seismicInfo[synchIndx][_WE.synchroniseStackTregistrationErrors ? 0 : j].systematicTProp;
			MLoffset_val = systematicTProp.value;


			_WE.tr.set("f1",(float)(tStartReff + MLoffset_val));

			try 
			  {
				out = new FileOutputStream(filename);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" + filename);
			  }

			// write out traces. True first, then synthetic
			// if (writeReflectivities), reflectivity third
			// NB k will start close to zero, but may not be exactly zero!
			// set a header field to help filter these
			_WE.tr.set("mark",(short)0);
			// true data
			S = S_true[i][j];
			for (k=0; k<n; k++)
			  {
				_WE.tr.set_data_value((float)(S[k]) ,k);
			  }
			_WE.suwriter.fputtr(out,_WE.tr);
			_WE.tr.set("mark",(short)1);
			// synthetic data
			S = S_synth[i][j];
			for (k=0;k<nZ;k++)
			  {
				_WE.tr.set_data_value((float)(S[k]),k);
			  }
			_WE.suwriter.fputtr(out,_WE.tr);

			if (WRITE_REFLECTIVITIES_TO_TRUE_PLUS_SYNTH_FILES)
			  {
				_WE.tr.set("mark",(short)2);
				// Reflection coefficients
				S = Reff[i][j];
				for (k=0;k<nZ;k++)
				  {
					_WE.tr.set_data_value((float)S[k],k);
				  }
				_WE.suwriter.fputtr(out,_WE.tr);
			  }
			if (SAVE_TRAJECTORY_VALUES)
			  {
				// X values
				_WE.tr.set("mark",(short)3);
				S = trajectoryX[i];

				for (k=0;k<nZ;k++)
				  {
					_WE.tr.set_data_value((float)S[k],k);
				  }
				// pre and postpad region with endpointvalues for simplicity of plotting
				preAndPostPadTrace(_WE.tr,nZ);
				_WE.suwriter.fputtr(out,_WE.tr);


				_WE.tr.zero_data_array();  // don't propogate any wierd padding...
				// Y values
				_WE.tr.set("mark",(short)4);
				S = trajectoryY[i];
				for (k=0;k<nZ;k++)
				  {
					_WE.tr.set_data_value((float)S[k],k);
				  }
				// pre and postpad region with endpointvalues for simplicity of plotting
				preAndPostPadTrace(_WE.tr,nZ);
				_WE.suwriter.fputtr(out,_WE.tr);
				_WE.tr.zero_data_array();   // don't propogate any  wierd padding...
			  }
			if (WRITE_REFLECTIVITIES_TO_TRUE_PLUS_SYNTH_FILES)
			  {
				if (_WE.dumpPandSreflectivity)
				  {
					// direct P reflectivity
					_WE.tr.set("mark",(short)5);
					S = Reff_P[i];
					for (k=0;k<nZ;k++) _WE.tr.set_data_value((float)S[k],k);
					_WE.suwriter.fputtr(out,_WE.tr);

					// direct S reflectivity
					_WE.tr.set("mark",(short)6);
					S = Reff_S[i];
					for (k=0;k<nZ;k++) _WE.tr.set_data_value((float)S[k],k);
					_WE.suwriter.fputtr(out,_WE.tr);

					// direct Density reflectivity
					_WE.tr.set("mark",(short)7);
					S = Reff_Rho[i];
					for (k=0;k<nZ;k++) _WE.tr.set_data_value((float)S[k],k);
					_WE.suwriter.fputtr(out,_WE.tr);
				  }
			  }

			// set this back to default 0
			_WE.tr.set("mark",(short)0);
			try
			  {
				_WE.suwriter.flush(out); // must do this!
				out.close();
				_WE.EH.println("Wrote su data to: " + filename);
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close" + filename);
			  }
		  }
	  }
  }

  /**
	 Simple routine to pre and post pad trace values 
	 with values from first and last nontrivial values
	 e.g. [0,0,0,.....0,1,5,2,3,0,0,.....0] goes to
	 [1,1,1,.....1,1,5,2,3,3,3,.....3]
  */
  private void preAndPostPadTrace(SuSegyFormat tr, int n)
  {
	int i;
	double y=0;

	// leading padding
	for (i=0;i<n;i++)
	  {
		if (tr.get_data_value(i)!=0)
		  {
			y = tr.get_data_value(i--);
			break;
		  }
	  }
	if (i==n) return; // (must be all zeros)
	for (;i>=0;i--) tr.set_data_value((float)y,i);
	for (i=n-1;i>=0;i--)
	  {
		if (tr.get_data_value(i)!=0)
		  {
			y = tr.get_data_value(i++);
			break;
		  }
	  }
	for (;i<n;i++) tr.set_data_value((float)y,i);
  }



  /**
	 Strawman: write out 
	 SU versions of current synthetics. spanning transverse range
	 ample enough range to cover well range.
	 @param makeTilted: switch on to produce a "tilted" synthetic minicube, with correct values
	 at the well location. Otherwise the minicube data set is flat (same for all traces)
	 
  */
  void writeOutFullSyntheticSUset(boolean makeTilted)
  {
	int i,j,k,n;
	int inLineCntr,xLineCntr,inLine,xLine,xLineMin,inLineMin;
	Wavelet wavelet;
	boolean MAKE_SU_DATA_SET_RIGHT_HANDED = true;// makes a right-handed data set

	double dt, DT_EPSILON=1e-4;                  // FIXME promote
	boolean ADD_PADDING_TO_MINICUBE = true;      // FIXME promote
	int TRANSVERSE_DISTANCE_TRACE_SPACING=20;    // FIXME promote
	double DISTORT_TIME = 0.004;                 // max ms distortion of seismic minicube between traces 
                                                 // from well 0.002 is mild, 0.02 is severe
	double aDistort=0, bDistort=0;

	double slope, largest_sigma_t, Zmin, Zmax, ZNext;
	double t,z,Xmin, Xmax, Ymin, Ymax;
	double val;
	double synth[]; // useful ptr
	UTM utmData;
	Spline zSpline;
	LinearInterpolation zInterp;
	FileOutputStream out;
	String filename;

	// get wavelet from current model wavlet
	wavelet =_WE.m.wavelet;


	n = S_synth[0][0].length;
	// set requisite fields
	_WE.tr.set("ns",(short)n);	
	_WE.tr.set("f1",(float)tStartReff);
	_WE.tr.set("dt",wavelet.dt);

	// times for amplitude spline
	double times[] = new double[n];
	for (j=0; j<n; j++) times[j]=tStartReff+j*_WE.dt;
	
	_WE.EH.printf("Writing synthetic seismic minicubes, noise level ="+Math.abs(_WE.syntheticNoiseLevel));
	_WE.EH.println(", using tracf as the (inline) X coordinate, and fldr as the (xline) Y coordinate");

	// long sTime=System.currentTimeMillis();

	// make minicubes traces for each well
	for (i=0;i<_m.nWells;i++)
	  {
		// determine transverse range of well log in window
		TimeToDepthMap t2dMap  = _WE.m.timeToDepthMaps[i];

		if (useSplinesForZTmapping)
		  {
			zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
								 Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		  }
		else
		  {
			zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);
		  }
		// gradient by FD's. Try to get natural dt
		dt = Math.abs(DT_EPSILON*(Math.max(1.0,t2dMap.t[t2dMap.t.length-1] - t2dMap.t[0])));
		DoubleArrayList list = new DoubleArrayList(t2dMap.sigma_t);
		largest_sigma_t = cern.jet.stat.Descriptive.max(list);

		// shallowest depth// get true depth at this time
		if (useSplinesForZTmapping) 
		  {
			Zmin = zSpline.splint(_WE.tStart,false);
			ZNext = zSpline.splint(_WE.tStart+dt,true);
		  }
		else
		  {
			Zmin = zInterp.interpolate(_WE.tStart,false);
			ZNext = zInterp.interpolate(_WE.tStart+dt,true);
		  }
		slope = Math.abs((ZNext - Zmin)/dt);
		Zmin -= FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * largest_sigma_t*slope;

		// deepest depth
		if (useSplinesForZTmapping) 
		  {
			Zmax = zSpline.splint(_WE.tEnd,false);
			ZNext = zSpline.splint(_WE.tEnd+dt,true);
		  }
		else
		  {
			Zmax = zInterp.interpolate(_WE.tEnd,false);
			ZNext = zInterp.interpolate(_WE.tEnd+dt,true);
		  }
		slope = Math.abs((ZNext - Zmax)/dt);
		Zmax += FundamentalConstants.T_UNCERTAINTY_N_STD_DEVS_FACTOR * largest_sigma_t*slope;

		// loop over UTM well trajectory to get largest x/y range associated with this true depth range
		// spline to convert MD's in checkshots/markers etc to TVD's
		utmData = _WE.logData[i].survey;
		zSpline = new Spline(utmData.Zm, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		// find where to start. Avoid reverse spline as well may not be monotonic
		for (j=0; j<utmData.Z.length; j++)
		  {
			if (utmData.Z[j]>Zmin)
			  {
				j--;
				break;
			  }
		  }
		// OK find biggest cube;
		Xmin = Ymin = Double.MAX_VALUE;
		Xmax = Ymax = -Double.MAX_VALUE;
		for (; j<utmData.Z.length ; j++)
		  {
			Xmin = Math.min(Xmin,utmData.X[j]);
			Ymin = Math.min(Ymin,utmData.Y[j]);
			Xmax = Math.max(Xmax,utmData.X[j]);
			Ymax = Math.max(Ymax,utmData.Y[j]);
			if (utmData.Z[j] > Zmax) break;
		  }
		if (_WE.verbosity>0) _WE.EH.println("Bounding box, well "+(i+1)+
									", Xmin="+Xmin+" Xmax="+Xmax+" Ymin="+Ymin+" Ymax="+Ymax);

		// OK now tStart, tEnd, Xmin, Xmax, Ymin, Ymax define box of data we need
		// Write out a seismic cube with the synthetic data stretched approximately
		// across this cube.

		// Spatial distortion to get interesting (but valid) seismic minicubes
		// Create data via
		//    S(x,y,t) = synth(t - aDistort*(x-x(t)) + bDistort*(y-y(t)))
		// where x(t) from x(z_MD(t)) etc. this guarantees S(x(t),y(t),t) = synth(t)
		// Get necessary splines & suitable distorion scale
		if (useSplinesForZTmapping) 
		  {
			zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
								 Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		  }
		else
		  {
			zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);
		  }

		Spline xLineSpline = new Spline(utmData.Z, utmData.X, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		Spline inLineSpline = new Spline(utmData.Z, utmData.Y, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		if (makeTilted)
		  {
			aDistort = DISTORT_TIME/TRANSVERSE_DISTANCE_TRACE_SPACING;
			bDistort = DISTORT_TIME/TRANSVERSE_DISTANCE_TRACE_SPACING;
		  }

		double xLineWell, inLineWell;
		

		for (j=0; j<_m.nStacks; j++)
		  {
			// suitable file for data
			filename = ("synth_seis."+Util.fileBasename(_WE.wellInfo[i].logFilename)+".stack."+(j+1)+".su");
			try 
			  {
				out = new FileOutputStream(filename);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" + filename);
			  }

			// write out mini-cube
			// synthetic data
			synth = S_synth[i][j];
			// amplitude spline
			Spline aSpline = new Spline(times, synth, n, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);

			double xLineRange,inLineRange;
			if (MAKE_SU_DATA_SET_RIGHT_HANDED)
			  {
				xLineRange = (Ymax-Math.floor(Ymin))/TRANSVERSE_DISTANCE_TRACE_SPACING;
				inLineRange = (Xmax-Math.floor(Xmin))/TRANSVERSE_DISTANCE_TRACE_SPACING;
				xLineMin = (int)Math.floor(Ymin);
				inLineMin = (int)Math.floor(Xmin);
			  }
			else
			  {
				xLineRange = (Xmax-Math.floor(Xmin))/TRANSVERSE_DISTANCE_TRACE_SPACING;
				inLineRange = (Ymax-Math.floor(Ymin))/TRANSVERSE_DISTANCE_TRACE_SPACING;
				xLineMin = (int)Math.floor(Xmin);
				inLineMin = (int)Math.floor(Ymin);
			  }
			if (ADD_PADDING_TO_MINICUBE)
			  {
				xLineMin -= TRANSVERSE_DISTANCE_TRACE_SPACING;
				inLineMin -= TRANSVERSE_DISTANCE_TRACE_SPACING;
				inLineRange +=2.01; // must have 3 traces if well is dead straight, min of 2 otherwise
				xLineRange +=2.01;
			  }

			for (xLine = xLineMin, xLineCntr=0; xLineCntr < (int)Math.ceil(xLineRange);
				 xLine += TRANSVERSE_DISTANCE_TRACE_SPACING, xLineCntr++)
			  {
				for (inLine=inLineMin,inLineCntr=0; inLineCntr<(int)Math.ceil(inLineRange);
					 inLine += TRANSVERSE_DISTANCE_TRACE_SPACING,inLineCntr++)
				  {
					// xLine direction
					//_WE.tr.set("fldr",xLine);
					_WE.tr.set(_WE.seismicXlineWord,xLine);
					// inline direction
					//_WE.tr.set("tracf",inLine);
					_WE.tr.set(_WE.seismicInlineWord,inLine);
					for (k=0;k<n;k++)
					  {
						// time of sample
						t = tStartReff+k*_WE.dt;
						// System.out.println(" xLine "+xLine+" inLine "+inLine+" t "+t);
						// well location at this time via t2d and UTM splines t->z->xLine,inLine 
						// True vertical depth
						if (useSplinesForZTmapping) 
						  {
							z = zSpline.splint(t,k!=0);
						  }
						else
						  {
							z = zInterp.interpolate(t,k!=0);
						  }						
						// => well location
						xLineWell = xLineSpline.splint(z,true);
						inLineWell = inLineSpline.splint(z,true);
						// equivalent distorted time, truncated to within range
						t += aDistort*(xLine - xLineWell) + bDistort*(inLine - inLineWell);
						if (t<times[0]) t=times[0];
						if (t>times[n-1]) t=times[n-1];
						// amplitude at equivalent time
						val = aSpline.splint(t,true);
						_WE.tr.set_data_value((float)val,k);
					  }
					_WE.suwriter.fputtr(out,_WE.tr);
				  }
			  }
			try
			  {
				_WE.suwriter.flush(out); // must do this!
				out.close();
				_WE.EH.println("Wrote su data to: " + filename);
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close" + filename);
			  }
		  }
	  }
  }

  /**
	 Write out 
	 SU versions of current synthetics, one trace per well	 
  */
  void writeOutAVOSyntheticSUset()
  {
	int i,j,k,n;

	Wavelet wavelet;

	double dt, DT_EPSILON=1e-4;                     // FIXME promote
	double AVO_GATHER_MAX_ANGLE = 50*Math.PI/180.0; // FIXME promote
	int N_AVO_GATHER_TRACES = 10;                   // FIXME promote




	double t,z,largest_sigma_t,Zmin,ZNext,Zmax;
	double val,theta,Xst;
	double synth[]; // useful ptr
	UTM utmData;
	Spline zSpline;
	LinearInterpolation zInterp;
	FileOutputStream out=null;
	String filename=null;


	// make required synthetic wavelet

	if (_WE.external_wavelet_filename==null)
	  {
		_WE.m.wavelet = new Wavelet(_WE, _WE.syntheticFreq, _WE.syntheticSampleTime, !_WE.BBsynthetic);
	  }
	else
	  {
		_WE.m.wavelet = new Wavelet(_WE,_WE.external_wavelet_filename);
		if (_WE.m.wavelet.dt != _WE.syntheticSampleTime)
		  {
			_WE.EH.println("Warning: dt of wavelet different to specified dt. Using wavelet dt!");
		  }
	  }
	//	_WE.m.wavelet = new Wavelet(_WE, _WE.syntheticFreq, _WE.syntheticSampleTime, true);
	wavelet =_WE.m.wavelet;


	n = S_synth[0][0].length;
	// set requisite fields
	_WE.tr.set("ns",(short)n);	
	_WE.tr.set("f1",(float)tStartReff);
    _WE.tr.set("delrt",(short)(Math.rint(1e3*tStartReff))); // in ms. Best approx to tStartReff
	_WE.tr.set("dt",wavelet.dt);

	// times for amplitude splineg
	double times[] = new double[n];
	for (j=0; j<n; j++) times[j]=tStartReff+j*_WE.dt;

	
	_WE.EH.println("Writing synthetic AVO gather traces");

	// do this for only one stack per well:
	_WE.m.nStacks = 1;
	_WE.lastWaveletFilenames = new String[_WE.m.nStacks];

	for (i=0; i<_m.nWells; i++)
	  {
		// suitable file for data
		filename = ("synth_seis."+Util.fileBasename(_WE.wellInfo[i].logFilename)+".AVO-gather.su");
		try 
		  {
			out = new FileOutputStream(filename);
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open" + filename);
		  }

		// determine transverse range of well log in window
		TimeToDepthMap t2dMap  = _WE.m.timeToDepthMaps[i];

		if (useSplinesForZTmapping)
		  {
			zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
								 Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		  }
		else
		  {
			zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);
		  }
		// gradient by FD's. Try to get natural dt
		dt = Math.abs(DT_EPSILON*(Math.max(1.0,t2dMap.t[t2dMap.t.length-1] - t2dMap.t[0])));

		DoubleArrayList list = new DoubleArrayList(t2dMap.sigma_t);
		largest_sigma_t = cern.jet.stat.Descriptive.max(list);

		// shallowest depth// get true depth at this time
		if (useSplinesForZTmapping) 
		  {
			Zmin = zSpline.splint(_WE.tStart,false);
			ZNext = zSpline.splint(_WE.tStart+dt,true);
		  }
		else
		  {
			Zmin = zInterp.interpolate(_WE.tStart,false);
			ZNext = zInterp.interpolate(_WE.tStart+dt,true);
		  }

		// deepest depth
		if (useSplinesForZTmapping) 
		  {
			Zmax = zSpline.splint(_WE.tEnd,false);
			ZNext = zSpline.splint(_WE.tEnd+dt,true);
		  }
		else
		  {
			Zmax = zInterp.interpolate(_WE.tEnd,false);
			ZNext = zInterp.interpolate(_WE.tEnd+dt,true);
		  }

		// loop over UTM well trajectory to get largest x/y range associated with this true depth range
		// spline to convert MD's in checkshots/markers etc to TVD's
		utmData = _WE.logData[i].survey;
		zSpline = new Spline(utmData.Zm, utmData.Z, utmData.N, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		if (useSplinesForZTmapping) 
		  {
			zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
								 Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
		  }
		else
		  {
			zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);
		  }
				
		// write out trace
		// synthetic data
		synth = S_synth[i][0];

		double dtheta = AVO_GATHER_MAX_ANGLE/(double)N_AVO_GATHER_TRACES;

		// start gather with -B value
		_WE.m.B_AVO[0].value = -1;

		_WE.tr.set("f2",(float)((1-N_AVO_GATHER_TRACES)*dtheta*180.0/Math.PI));
		_WE.tr.set("d2",(float)(dtheta*180.0/Math.PI));


		for (j=-N_AVO_GATHER_TRACES; j<=N_AVO_GATHER_TRACES; j++)
		  {
			// flip to +ve B for +ve angles
			if (j==0) _WE.m.B_AVO[0].value = 1;

			// make a pseudo stack (Xmax=Min) at angle theta 
			theta = j*dtheta;
			Xst = theta * _WE.seismicInfo[i][0].T_event * _WE.seismicInfo[i][0].V_stack;
			_WE.seismicInfo[i][0].offsetXmax = _WE.seismicInfo[i][0].offsetXmin = Xst;
			if (_WE.seismicInfo[i][0].isStackAnglePrescribed()) 
			  {
				_WE.seismicInfo[i][0].theta_specified = Math.abs(theta);
			  }
			_WE.seismicInfo[i][0].computeStackParams();

			// recompute reflectivities and synthetic
			if (_WE.useSaturationModel) _WE.m.recomputeLogProperties(false);

			_WE.m.computeReflectivities(false,false);

			if (_WE.useSaturationModel) 
			  {
				// copy reflectivities
				_WE.m.makeReflectivityCopy();
				// compute effect of sat change,
				_WE.m.recomputeLogProperties(true);
				// new reflectivities
				_WE.m.computeReflectivities(false,false);
				// then difference
				_WE.m.computeReflectivityDifference();
			  }
			_WE.convModel.computeSyntheticSeismic();

			// amplitude spline
			Spline aSpline = new Spline(times, synth, n, Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);

			for (k=0;k<n;k++)
			  {
				// time of sample
				t = tStartReff+k*_WE.dt;
				// System.out.println(" xLine "+xLine+" inLine "+inLine+" t "+t);
				// well location at this time via t2d and UTM splines t->z->xLine,inLine 
				// True vertical depth
				if (useSplinesForZTmapping) 
				  {
					z = zSpline.splint(t,k!=0);
				  }
				else
				  {
					z = zInterp.interpolate(t,k!=0);
				  }						
				if (t<times[0]) t=times[0];
				if (t>times[n-1]) t=times[n-1];
				val = aSpline.splint(t,true);
				_WE.tr.set_data_value((float)val,k);
			  }

			// let f2 be the angle
			_WE.tr.set("f2",(float)(j*dtheta*180.0/Math.PI));

			_WE.suwriter.fputtr(out,_WE.tr);
		  }
	  }
	try
	  {
		_WE.suwriter.flush(out); // must do this!
		out.close();
		_WE.EH.println("Wrote su data to: " + filename);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close" + filename);
	  }
  }


  /**
	 get suitable padding length for data sequence R (power of two above length of R)
  */
  private int getPaddinglength(double R[])
  {
	int n,ng;

	n=R.length;
	ng=0;
	while ((n >>= 1) > 0 ) ng++;
	// In general R.length not a power of 2, so find the right n = 2^ng  bigger than R.length
	if (R.length != (1<<ng))
	  {
		// R.length not a power of 2
		n=(1<<(ng+1));
	  }
	else
	  n = R.length; // it was already a power of 2
	return n;
  }

  /**
	 Debugging routine for Z-t spline
  */
  private void debugZtSpline(TimeToDepthMap t2dMap)
  {
	int k;
	double t,Z;
	String filename = "junk.txt";

	_WE.EH.println("Printing spline to : "+filename);
	PrintWriter out;
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(Util.fileBasename(filename))));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }
	Spline zSpline = new Spline(t2dMap.t, t2dMap.Z, t2dMap.t.length, 
										Double.MAX_VALUE, Double.MAX_VALUE,_WE.EH);
	// coarse data
	for (k=0; k<t2dMap.t.length; k++)
	  {
		out.println(t2dMap.t[k]+" "+t2dMap.Z[k]);
	  }
	out.println("");

	// starting indx and time
	k = (int)Math.ceil((_WE.tStart - tStartReff)/_WE.dt);
	t = tStartReff + _WE.dt * (double)k;
	for (; t<=_WE.tEnd; t+=_WE.dt)
	  {
		// get true depth at this time
		Z = zSpline.splint(t,false);
		out.println(t+" "+Z);
	  }
	out.close();
  }

  /**
	 data class holding reflectivity arrays for time-lapse problems
	 typically working space for reflectivity at vintage 0
  */
  public static class BaseVintageReflectivities
  {
  /**
	 Arrays of (optional) normal incidence reflectivities on equal spacing for each stack
	 Reff_P[i][k] = reflection for ith well, sample time k.
	 (Should not differ with stack I think)
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_P;
  /**
	 Arrays of (optional) shear reflectivities on equal spacing for each stack
	 (Should not differ with stack I think)
	 Reff_S[i][k] = shear reflection for ith well, sample time k.
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_S;
  /**
	 Arrays of (optional) density reflectivities on equal spacing for each stack
	 (Should not differ with stack I think)
	 Reff_S[i][k] = density reflection for ith well, sample time k.
	 Available only with --dump-PS-reflectivity
  */
  double[][] Reff_Rho;

  /**
	 Arrays of reflectivities on equal spacing for each stack
	 Reff[i][j][k] = reflection for ith well, jth stack, sample time k
  */
  double[][][] Reff;
  }
}
