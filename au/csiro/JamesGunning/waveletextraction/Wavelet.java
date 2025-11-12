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

import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;
import static au.csiro.JamesGunning.lib.util.MathUtil.*;




/**
   Basic set up to describe a wavelet
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class Wavelet
{
 /**
	 This constant guarantees samples of a Ricker-2 wavelet generated
	 random signal will pick up 95% of the equivalent degrees-of-freedom (read energy)
	 if sampled at a time rate dt = SAMPLE_TIME_MAGIC_CONSTANT/f_peak
	 where f_peak is the peak spectral frequency of the Ricker-2 spectrum
  */
  public static final double SAMPLE_TIME_MAGIC_CONSTANT = 0.253;
  /**
	 This constant denotes the fraction of the seismic power spectrum
	 that sets the Nyquist sample rate of the knot points that define
	 the wavelet, as related to the above.
  */
  public static final double BAND_EDGE_FRACTION = 0.95;
  /**
	 Ricker-2 wavelets have this ratio of the mean frequency (first
	 moment of the power spectrum) to the peak frequency
  */
  private static final double MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER=1.124;

 /**
	Calling Extractor
 */
  private WaveletExtractor _WE;

/**
	Calling Extractor's error handler
 */
  private ErrorHandler _EH;

  /** 
	  coarse index of start of wavelet. Sample here is zero. Only coarse samples
	  <i>between</i> firstCoarseSample and finalCoarseSample are free (nonzero)
   */
  int firstCoarseSample;
  /** 
	  coarse index of end of wavelet. Sample here is zero. Only coarse samples
	  between firstCoarseSample and finalCoarseSample are free (nonzero)
   */
  int finalCoarseSample;
  
  /** 
	  index of start of wavelet
   */
  int firstFineSample;
  /** 
	  index of end of wavelet
   */
  int finalFineSample;
  /** 
	  precursor time span of wavelet (secs)
   */
  double precursorDuration;
  /** 
	  coda time span of wavelet (secs)
   */
  double codaDuration;

  /**
	 Sampling rate of coarse samples
  */
  double dtCoarse;
  /**
	 time values on coarse scale
  */
  double tCoarse[];

  /**
	 Sequence of samples at coarse scale defined by Nyquist rate
  */
  double valuesCoarse[];
  /**
	 property managers for the coarse coefficients
  */
  Property[] valuesCoarseProp;

  /**
	 wavelet sampling time in microseconds
	 (is a short header field)
  */
  short dt;
  /**
	 wavelet sampling time in seconds. Always consistemt with dt
  */
  double delta_t;

  /**
	 Sequence of samples at standard seismic scale defined by seismic dt (2ms,4ms etc),
	 obtained by Splines from coarse values.
	 values[0][0..] is basic wavelet, values[1] is 
	 first stack stretched wavelet etc
  */
  double values[][];

  /**
	 Commencement time of wavelet in seconds
	 (same as su's f1: a <b>float</b> in the headers)
  */
  double[] f1;

  /** 
	  Spline object to connect coarse and fine samples
  */
  Spline tSpline;

  /**
	 phase spectrum as computed by computePhaseSpectrum() 
	 is accessible here after calling the method computePhaseSpectrum().
	 phase[i] = atan(phaseIm[i]/phaseRe[i]);
  */
  double[] phase;
  double[] phaseRe;
  double[] phaseIm;


  /**
	 ratio of second peak to primary peak which is used in certain branching decisions
	 related to peak detection
  */
  public static final double SECOND_PEAK_CRITICAL_RATIO = 0.8;

  /**
	 Variable indicating status of second peak which is usable and
	 sensble directly after a getPeakTime() call. Will be true if
	 second peak has amplitude A_2: 
	 |A_2| < SECOND_PEAK_CRITICAL_RATIO * | A_peak |
  */
  boolean secondPeakIsStrong = false;

  /**
	 The frequency spacing of the phase spectrum phase[] is also available here
	 after a call to computePhaseSpectrum();
  */
  double phaseVectorFrequencySpacing;

  /**
	 Epsilon used in (one method) of determining the band edge 
	 in the estimateWaveletIndependentDT()
	 method. It's the ratio of the FFT at the bandedge to the peak FFT value
  */
  public static final double BAND_EDGE_EPSILON=0.01; 


  // Fields used in peak detection
  private double tTarget;
  private boolean seekNearestPeak=false;
  double peakTime;

  // ----------------------- Stretched wavelet fields (start) ----------------------------

  /**
	 The far stack wavelet may be derived from the basic wavelet by stretching.
	 Here we define some basic fields that may be of use in such a model.
	 The may not be used in most simple cases
  */

  boolean useStretchedModel = false;

  /**
	Maximum stretch allowed for embedded stretched wavelet if used
  */
  public static double MAX_WAVELET_STRETCH = 2.0;
  /** 
	  index of start of Stretched wavelet
   */
  int firstStretchedFineSample;
  /** 
	  index of end of Stretched wavelet
   */
  int finalStretchedFineSample;
  /** 
	  precursor time span of Stretched wavelet (secs)
   */
  double precursorStretchedDuration;
  /** 
	  coda time span of Stretched wavelet (secs)
   */
  double codaStretchedDuration;

  /**
	 number of wavelets embedded in this representation via stretching
  */
  int nWavelets=1;

 // ----------------------- Stretched wavelet fields (end) ----------------------------


  /**
	 Make new wavelet object with parameters computed from the
	 bandedge.  If makeStandardRickerType==true, will create a full
	 Wavelet object with a standard Ricker-2 type wavelet.
  */
  public Wavelet(WaveletExtractor WE, double fBandEdge, double suggested_sampling_dt, 
				 boolean makeStandardRickerType)
  {
	this._WE = WE;
	this._EH = _WE.EH;

	if (makeStandardRickerType)
	  {
		makeRickerWavelet(fBandEdge,suggested_sampling_dt);
	  }
	else
	  {
		// bband
		makeBBWavelet(fBandEdge,suggested_sampling_dt);
	  }
  }
  
  /**
	 Make new wavelet object with coarse samples extending (inclusive) from
	 firstCoarseSample ...finalCoarseSample, spaced at coarse_scale_dt.
	 The fine scale wavelet at spacing fine_scale_dt is setup to be computed from this using
	 cubic splines.<b>
	 Must have finalCoarseSample > firstCoarseSample+1

	 If useStretchedModel == true, constructor will create a model
	 with a second, stretched internal representation, mapped from the
	 primary variables, which models a dispersed wavelet of different
	 amplitude.
  */
  public Wavelet(WaveletExtractor WE, int firstCoarseSample, int finalCoarseSample, 
				 double fine_scale_dt, double coarse_scale_dt)
  {
	int i,j;

	this._WE = WE;
	this._EH = _WE.EH;
	this.useStretchedModel = WE.useStretchedModel;

	dtCoarse = coarse_scale_dt;
	dt=(short)(fine_scale_dt*1e6);
	delta_t = dt*1e-6;
	this.firstCoarseSample = firstCoarseSample;
	this.finalCoarseSample= finalCoarseSample;


	if (finalCoarseSample <= firstCoarseSample+1)
	  {
		// silly wavelet
		_EH.die_err("Wavelet must have positive span");
	  }
	// x values for abscissa in spline
	tCoarse = new double[finalCoarseSample-firstCoarseSample+1];
	for (i=firstCoarseSample,j=0; i<=finalCoarseSample; i++,j++) 
	  {
		tCoarse[j] = i* dtCoarse;
	  }

	// set up coarse values. First and last are initialised to zero
	valuesCoarse = new double[finalCoarseSample-firstCoarseSample+1];
	// the properties array need only point at internal points
	valuesCoarseProp = new Property[finalCoarseSample-firstCoarseSample-1];
	for (i=firstCoarseSample+1,j=0; i<finalCoarseSample; i++,j++) 
	  {
		valuesCoarseProp[j] = new Property(0,false,false);
	  }


	if (WE.useWaveletPrior)
	  {
		if ((finalCoarseSample - firstCoarseSample -1 ) != _WE.waveletPrior.amp_means.size())
		  {
			_WE.EH.die_err("Wavelet spans must be specified the same as "+
						   "the prior-wavelet derivation if using prior wavelets");
		  }
		if (Math.abs(tCoarse[1] - _WE.waveletPrior.f1_prior)/delta_t > 1e-3)
		  {
			_WE.EH.die_err("tCoarse[1] "+tCoarse[1]+" _WE.waveletPrior.f1_prior"+_WE.waveletPrior.f1_prior);
			_WE.EH.die_err("Wavelet spans must be specified so as to have the same precursor/coda"+
						   " as the prior-wavelet derivation wavelet: "+tCoarse[1] +" "+_WE.waveletPrior.f1_prior);
		  }

		// OK, alter the prior:
		for (i=firstCoarseSample+1,j=0; i<finalCoarseSample; i++,j++) 
		  {
			valuesCoarseProp[j].setMean(_WE.waveletPrior.amp_means.get(j));
		  }
	  }




	// Associated fine scale factors;
	finalFineSample = (int)Math.floor(finalCoarseSample*dtCoarse/delta_t);
	codaDuration = finalFineSample*delta_t;

	firstFineSample = (int)Math.ceil(firstCoarseSample*dtCoarse/delta_t);
	precursorDuration = -firstFineSample*delta_t;

	// Multi wavelet entries
	values = new double[WE.m.nStacks][];
	f1 = new double[WE.m.nStacks];




	// common wavelet
	values[0] = new double[finalFineSample-firstFineSample+1];
	f1[0] = firstFineSample*delta_t;
  
	// Stretched wavelet entries
	setUpStretching();
	
	// phase spectrum stuff

	int n = FFT.getNearestPowerOfTwo(values[0].length);
	phaseVectorFrequencySpacing = 1.0/((double)n * delta_t);
	phase = new double[n/2];
	phaseRe = new double[n/2];
	phaseIm = new double[n/2];


	// Stretched wavelet entries
// 	if (WE.useStretchedModel)
// 	  {
// 		// Associated fine scale factors;
// 		finalStretchedFineSample = (int)Math.floor(finalCoarseSample*dtCoarse*MAX_WAVELET_STRETCH/delta_t);
// 		codaStretchedDuration = finalStretchedFineSample*delta_t;

// 		firstStretchedFineSample = (int)Math.ceil(firstCoarseSample*dtCoarse*MAX_WAVELET_STRETCH/delta_t);
// 		precursorStretchedDuration = -firstStretchedFineSample*delta_t;
		
// 		stretchedValues = new double[finalStretchedFineSample - firstStretchedFineSample + 1];
// 		f1Stretched = firstStretchedFineSample*delta_t;
// 	  }

  }
  
  /**
	 Read in and set up a wavelet from the named wavelet file. This is only
	 used in conjunction with synthetics.
  */
  public Wavelet(WaveletExtractor WE, String wavelet_filename)
  {
	int i,j,ns;
	FileInputStream in;
	Wavelet wavelet=null;
	double[] _values;

	this._WE = WE;
	this._EH = _WE.EH;

	// open the wavelet file.
	try
	  {
		in = new FileInputStream(wavelet_filename);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + wavelet_filename);
	  }
	if (WE.sureader.fgettr(in,WE.trW)<1)
	  throw new RuntimeException("No wavelet data in " + wavelet_filename);
	// simple tests on critical fields
	if ((ns = (int) WE.trW.get_short("ns"))<2)
	  throw new RuntimeException("Wavelet is unacceptably short " + wavelet_filename);
	if (WE.trW.get_short("dt")==0)
	  throw new RuntimeException("Wavelet has weird dt=0 in " + wavelet_filename);

	dt = WE.trW.get_short("dt");
	delta_t = dt*1e-6;

	// For want of a better choice, set the coarse samples to that from the wavelet file
	dtCoarse = ((double)dt)/1.0e6;

	f1 = new double[WE.m.nStacks];

	if (WE.trW.get_short("delrt") != 0)
	  {
		// looks like delrt field defines window start (ms)
		f1[0] = (double)WE.trW.get_short("delrt")*1e-3; // to secs
	  }
	else
	  {
		f1[0] = (double) WE.tr.get_float("f1");
	  }

	// f1/dtCoarse should be an integer, but round to be sure fp error doesn't kill us
	firstCoarseSample = (int)Math.rint(f1[0]/dtCoarse);
	finalCoarseSample = ns - 1 + firstCoarseSample;
	tCoarse = new double[finalCoarseSample-firstCoarseSample+1];
	valuesCoarse = new double[finalCoarseSample-firstCoarseSample+1];

	// the properties array need only point at internal points
	valuesCoarseProp = new Property[finalCoarseSample-firstCoarseSample-1];
	for (i=firstCoarseSample+1,j=0; i<finalCoarseSample; i++,j++) 
	  {
		valuesCoarseProp[j] = new Property(0,false,false);
	  }

	// Associated fine scale factors. Make sure these are the same as the minicubes;

	finalFineSample = (int)Math.rint(finalCoarseSample*dtCoarse/delta_t);
	codaDuration = finalFineSample*delta_t;

	firstFineSample = (int)Math.rint(firstCoarseSample*dtCoarse/delta_t);
	precursorDuration = -firstFineSample*delta_t;

	// Multi wavelet entries
	values = new double[WE.m.nStacks][];
	values[0] = new double[finalFineSample-firstFineSample+1];

	// Get coarse wavelet samples.
	// In trace, wavelet is in ordinary form.
	for (j=firstCoarseSample,i=0; i<ns; i++,j++)
	  {
		tCoarse[i] = j*dtCoarse;
		valuesCoarse[i]  = (double) WE.trW.get_data_value(i);
	  }

	// make Spline interpolator
	tSpline = new Spline(tCoarse,valuesCoarse,tCoarse.length,0.0,0.0,_EH);

	
	// make fine wavelet by cubic interpolation
	// complete calculation of fine wavelet by cubic interpolation
	// coda
	_values = values[0];

	for (i=0; i<=finalFineSample; i++)
	  {
		_values[i] = tSpline.splint(i*delta_t,true);
	  }
	// precursor: wraparound bit
	for (j=finalFineSample-firstFineSample,i=-1; i>=firstFineSample; i--,j--)
	  {
		_values[j] = tSpline.splint(i*delta_t,true);
	  }


	setUpStretching();
	if (_WE.useStretchedModel)
	  {
		for (i=0; i < (_WE.m.nStacks-1); i++)
		  {
			_WE.m.farWavStretch[i] = _WE.m.farStackWaveletTimeStretchMean;
			_WE.m.farWavAmp[i] = _WE.m.farStackWaveletAmplitudeMean;
		  }
		computeFarWaveletValues();
	  }


	// close file
	try
	  {
		in.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + wavelet_filename);
	  }
  }


  private void setUpStretching()
  {
	int i;

	if (_WE.useStretchedModel)
	  {
		// Associated fine scale factors;
		finalStretchedFineSample = (int)Math.floor(finalCoarseSample*dtCoarse*MAX_WAVELET_STRETCH/delta_t);
		codaStretchedDuration = finalStretchedFineSample*delta_t;
		
		firstStretchedFineSample = (int)Math.ceil(firstCoarseSample*dtCoarse*MAX_WAVELET_STRETCH/delta_t);
		precursorStretchedDuration = -firstStretchedFineSample*delta_t;
		for (i=1; i<_WE.m.nStacks; i++)
		  {
			values[i] = new double[finalStretchedFineSample - firstStretchedFineSample + 1];
			f1[i] = firstStretchedFineSample*delta_t;
		  }
	  }
  }


  /*
	Invoke this every time the wavelet coarse parameters are changed. It maps the
	coarse Parameters onto the coarse grid, recomputes the spline, and then 
	recomputes the fine values.

	If the wavelet carries a stretched representation as well
	(useStretchedModel == true), the scale and time stretch in
	_WE.m.farWavStretch
	and
	_WE.m.farWavAmp
	are
	applied to the stretched wavelets. The calling routine *must* guarantee
	that 1 <= _WE.m.farWavStretch[k] <= MAX_WAVELET_STRETCH.

	If useStretchedModel == false, the arguments are merely ignored.
  */
  void updateWavelet()
  {
	int i,j,k;
	double t;
	double[] _values;

	// map Property values onto coarse samples
	// endpoints are left deliberately at zero
	for (i=0;i<valuesCoarseProp.length;i++) 
	  {
		valuesCoarse[i+1] = valuesCoarseProp[i].value;
	  }
	// make Spline interpolator
	tSpline = new Spline(tCoarse,valuesCoarse,tCoarse.length,0.0,0.0,_EH);

	// complete calculation of fine wavelet by cubic interpolation
	// (1) first (near) stack..
	_values = values[0];
	// coda
	for (i=0; i<=finalFineSample; i++)
	  {
		_values[i] = tSpline.splint(i*delta_t,true);
	  }
	// precursor: wraparound bit
	for (j=finalFineSample-firstFineSample,i=-1; i>=firstFineSample; i--,j--)
	  {
		_values[j] = tSpline.splint(i*delta_t,true);
	  }

	// stretched wavelet
	if (useStretchedModel) computeFarWaveletValues();
  }

  private void computeFarWaveletValues()
  {
	int i,j,k;
	double t,stretchTimeFac,stretchAmp;
	double[] _values;

	for (k=0; k < (_WE.m.nStacks-1); k++)
	  {
		_values = values[k+1];
		
		stretchTimeFac = Math.max(1,Math.min(_WE.m.farWavStretch[k],MAX_WAVELET_STRETCH));
		stretchAmp = _WE.m.farWavAmp[k];
		// complete calculation of fine wavelet by cubic interpolation
		// coda
		for (i=0; i<=finalStretchedFineSample; i++)
		  {
			t = i*delta_t/stretchTimeFac;
			_values[i] = (t < tCoarse[tCoarse.length-1]) ? stretchAmp * tSpline.splint(t,true) : 0.0;
		  }
		// precursor: wraparound bit
		for (j=finalStretchedFineSample-firstStretchedFineSample,i=-1; i>=firstStretchedFineSample; i--,j--)
		  {
			t = i*delta_t/stretchTimeFac;
			_values[j] = (t > tCoarse[0]) ? stretchAmp * tSpline.splint(t,true) : 0.0;
		  }
	  }
  }


  /**
	Time shift in-situ wavelet by the delay
  */
  void embedTimeShift(double delay)
  {
	int i;
	double t,tmin,tmax;
	boolean debug=false;

	// endpoints are left deliberately at zero
	tmin=tCoarse[0];
	tmax=tCoarse[tCoarse.length-1];

	// System.out.println("Wavelet shift embedding delay = "+delay);
	if (debug)
	  {
		for (t=tmin; t<=tmax; t+=dtCoarse/10.0)
		  {
			System.out.println("WC "+t+" "+tSpline.splint(t,true));
		  }
		System.out.println("WC");
	  }

	// only impose shift on internal coarse properties (ends fixed at zero)
	for (i=1;i<valuesCoarseProp.length;i++) 
	  {
		// preset to zero for safety
		valuesCoarseProp[i-1].value = 0;
		t = tCoarse[i] + delay;
		if ((t >= tmin) && (t <= tmax))
		  {
			valuesCoarseProp[i-1].value = tSpline.splint(t,false);
		  }
	  }

	updateWavelet();

	if (debug)
	  {
		for (t=tmin; t<=tmax; t+=dtCoarse/10.0)
		  {
			System.out.println("WC "+t+" "+tSpline.splint(t,true));
		  }
		System.out.println("WC");
	  }
  }

  /**
	 Get mean wavelet time (a "proxy" for phase) computed as the sample mean of
	 time, weighted by the current wavelet coefficients squared. To be used
	 only after update wavelet is called. If use8thPowerWeights==true,
	 we use |coefficients|^8 as the weight. This emphasizes the largest peak. 
  */
  double getMeanWaveletTime(boolean use8thPowerWeights)
  {
	int i,j;
	double sum,norm,wt,tMean;

	sum=0.0;
	norm=0.0;

	for (i=firstCoarseSample,j=0; i<=finalCoarseSample; i++,j++) 
	  {
		wt = sqr(valuesCoarse[j]*valuesCoarse[j]);
		if (use8thPowerWeights) wt = sqr(wt*wt);
		norm += wt;
		sum += tCoarse[j]*wt;
	  }

	tMean = (norm>0) ? (sum/norm) : 0.0;


	return tMean;
  }

  /**
	 Experimental: 
	 Get wavelet peak time (a "proxy" for phase) computed from parabolic fit of
	 5 points centred on largest sample. To be used
	 only after update wavelet is called
  */
  double getPeakTime5ptfit(boolean debug)
  {
	int i,j,jMax=0,iMax=0;
	double max,den,tp=0;

	max=0;
	for (i=firstCoarseSample,j=0; i<=finalCoarseSample; i++,j++) 
	  {
		if (Math.abs(valuesCoarse[j]) > max)
		  {
			jMax=j;
			iMax=i;
			max = Math.abs(valuesCoarse[j]);
		  }
	  }
	if (jMax<(valuesCoarse.length-2) && jMax>1)
	  {
		double a,b;
		a = (2*valuesCoarse[jMax-2] -valuesCoarse[jMax-1]
			 -2*valuesCoarse[jMax] -valuesCoarse[jMax+1] +2*valuesCoarse[jMax+2])/(14.0*dtCoarse*dtCoarse);
		b =  (-2*valuesCoarse[jMax-2] -valuesCoarse[jMax-1]
			  +valuesCoarse[jMax+1] +2*valuesCoarse[jMax+2])/(10.0*dtCoarse);
		tp = ((a!=0) ? -b/(2.0*a) : 0) + tCoarse[jMax];
		if (debug)
		  {
			System.out.println("Jmax = "+jMax);
			double c = (-3* valuesCoarse[jMax-2] + 12* valuesCoarse[jMax -1]+ 
						17 * valuesCoarse[jMax]+ 12 * valuesCoarse[jMax+ 1]- 3* valuesCoarse[jMax+ 2])/35.0;
			for (i=firstCoarseSample; i<=finalCoarseSample; i++) 
			  {
				System.out.println("WC "+i*dtCoarse+" "+tSpline.splint(i*dtCoarse,true));
			  }
			System.out.println("WC");
			double t;
			for (t=tCoarse[jMax-2];t<=tCoarse[jMax+2];t+=dtCoarse/5.0)
			  {
				double tt =  t - tCoarse[jMax];
				System.out.println("WC "+ t + " "+(a*tt*t+b*tt+c));
			  }
		  }
	  }

	return tp;
  }

  /**
	 Query for status of secondary peak in wavelet. Call after a
	 getPeakTime() call. Returns true if |second peak| is larger than
	 SECOND_PEAK_CRITICAL_RATIO times the amplitude of the largest
	 peak
  */
  boolean isSecondPeakStrong()
  {
	return secondPeakIsStrong;
  }

  /**
	 Get wavelet peak time computed from cubic spline
	 and analytical maximum/minimum. Use one that is nearest to
	 the specified time t.
  */
  double getNearestPeakTime(double t, boolean constrainPeakArrivalSign, boolean peakArrivalPositive)
  {
	seekNearestPeak = true;
	tTarget = t;
	double tp = getPeakTime(false, constrainPeakArrivalSign, peakArrivalPositive);
	seekNearestPeak = false;

	return tp;
  }
  /**
	 Get wavelet peak time computed from cubic spline and analytical
	 maximum/minimum.  If useWaveletMoments = false (internal, and the
	 default), uses strict peak-hunting, otherwise fake moment-based
	 method is used. For the default method, if
	 restrictToSpecifiedSign==true, will look only for peaks of sign
	 given by "seekPositive", otherwise "seekPositive" is ignored
  */
  double getPeakTime(boolean debug, boolean restrictToSpecifiedSign, boolean seekPositive)
  {
	double tp=0,tmp;
	boolean useWaveletMoments = false;

	if (useWaveletMoments) 
	  {
		// Use 8th-power moment to estimate peak position (always differentable)
		tp = getMeanWaveletTime(true);

		// Also, have no knowledge of second peak, so set it to weak by default
		secondPeakIsStrong = false;

		peakTime = tp;
		return tp;
	  }
	else
	  {
		// use peak-detection algorithm. Usually good, but has jumps.

		secondPeakIsStrong = false; // assume this until shown otherwise

		// get peak time of spline
		tp = tSpline.getXatMinOrMax(restrictToSpecifiedSign,seekPositive);

		if (tp==tCoarse[0] || tp == tCoarse[tCoarse.length-1]) tp=0;

		// extra code to allow smooth transition between two competing peaks
		double tp_comp, yp, xpeak2, ypeak, ypeak2, theta;
		double[] y = tSpline.ya;
		double[] xa = tSpline.xa;
		double[] y2a = tSpline.y2a;
		int i,j,jMax=0, jMax2=0;
		int n = tSpline.n;
		double max=0, max2=0;
		// find largest |sample of specified flavour
		if (!restrictToSpecifiedSign)
		  {
			// find largest |sample|
			for (j=0; j<n; j++) 
			  {
				if (Math.abs(y[j]) > max)
				  {
					jMax=j;
					max = Math.abs(y[j]);
				  }
			  }
		  }
		else
		  {
			if (seekPositive)
			  {
				max = y[0];
				for (j=0; j<n; j++) 
				  {
					if (y[j] > max)
					  {
						jMax=j;
						max = y[j];
					  }
				  }
			  }
			else
			  {
				max = y[0];
				for (j=0; j<n; j++) 
				  {
					if (y[j] < max)
					  {
						jMax=j;
						max = y[j];
					  }
				  }
			  }
		  }

		// If it's internal, see if secondary peak is also internal.
		if (jMax<=(n-2) && jMax>=1)
		  {
			boolean isPeak,isTrough;
			// do sequential hunt. Are guaranteed y=0 at ends
			if (!restrictToSpecifiedSign)
			  {
				// Unsigned search
				for (j=1; j<(n-1); j++) 
				  {
					isPeak = ((y[j] > y[j-1]) && (y[j] > y[j+1]));
					isTrough = ((y[j] < y[j-1]) && (y[j] < y[j+1]));
					if ((isPeak || isTrough) && (j!=jMax) && (Math.abs(y[j]) > max2))
					  {
						jMax2=j;
						max2 = Math.abs(y[j]);
					  }
				  }
			  }
			else
			  {
				// Signed search
				for (j=1; j<(n-1); j++) 
				  {
					isPeak = ((y[j] > y[j-1]) && (y[j] > y[j+1]));
					isTrough = ((y[j] < y[j-1]) && (y[j] < y[j+1]));
					if ((isPeak && seekPositive && (y[j]>max2) || (isTrough && !seekPositive && y[j]<max2)) && (j!=jMax))
					  {
						jMax2=j;
						max2 = y[j];
					  }
				  }
			  }

			// if second-biggest sample was not adjacent, and internal, we possibly have competing peaks
			if ((jMax2<=(n-2)) && (jMax2>=1)           // internal
				&& (Math.abs(jMax2-jMax)>=2))          // not adjacent
			  {
				// Value of first peak
				ypeak = tSpline.splint(tp,true);

				// get nature of second peak
				int klo,khi;
				boolean maximise;
				double h,h2,h4,b,a,x,det,A;
				double min, xMin, xMax;

				maximise = ((y[jMax2] >= y[jMax2+1]) && (y[jMax2] >= y[jMax2-1]));

				// loop over two intervals on either side of max
				max = -Double.MAX_VALUE;
				min = Double.MAX_VALUE;
				xMin = xMax = xa[jMax2]; // overwritten almost certainly

				for (klo=jMax2-1;klo<=jMax2;klo++)
				  {
					khi = klo+1;	
					if (klo<0 || khi>(n-1)) continue;

					h = xa[khi]-xa[klo];
					if (h == 0.0) 
					  {
						throw new RuntimeException("Bad x in getXatMinOrMax()");
					  }
					h2 = h*h;
					h4 = h2*h2;

					// in this interval, the cubic spline has a parabolic derivative. Find its zeros
					tmp = y2a[klo] - y2a[khi];
					det = 36.0 * h4 * y2a[khi] * y2a[khi] + 
					  12.0 * h2 * tmp *(-6.0*y[klo] + 6.0*y[khi] + h2*(y2a[klo] + 2.0*y2a[khi]));
					if (det>=0 && tmp != 0)
					  {
						det = Math.sqrt(det);
						// loop over two roots of quadratic: j is the +/- sign.
						for (j=-1;j<=1;j+=2)
						  {
							A = j*((-j)*6.0*h2*y2a[khi] + det)/(6*h2*tmp);
							if (A>=0 && A<=1) // root is within interval
							  {
								x = xa[khi] - h*A;
								a = (xa[khi]-x)/h;
								b = (x-xa[klo])/h;
								yp = a*y[klo]+b*y[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
								if (yp>max)
								  {
									max = yp; xMax=x;
								  }
								if (yp<min)
								  {
									min = yp; xMin = x;
								  }
							  }
						  }
					  }
				  }
				xpeak2 = maximise ? xMax : xMin;
				ypeak2 = maximise ? max : min;


				// Now the two peaks are at ypeak2 and ypeak, with times xpeak2 and tp respectively.
				// Treat this in on of two different styles:
				if (seekNearestPeak)
				  {
					// (a) Snap to peak nearest requested target time
					if (Math.abs(tp - tTarget) <= Math.abs(xpeak2 - tTarget))
					  {
						peakTime = tp;
						return tp;
					  }
					else
					  {
						peakTime = xpeak2;
						return xpeak2;
					  }
				  }
				else
				  {
					// (b) Adjudicate between competing peaks
					theta = Math.abs(ypeak2/((ypeak!=0)? ypeak : 1));
					if (theta>1.0)
					  {
						// ypeak2 is actually bigger than ypeak: swap times
						tmp = xpeak2; 
						xpeak2 = tp; 
						tp = tmp;
						theta = 1.0/theta;
					  }
					double theta_c = SECOND_PEAK_CRITICAL_RATIO; // critical threshold for competition
					// Test if peaks seem to be competing and are the same sign.
					// If so, strike a compromise
					if (theta > theta_c)
					  {
						// mark second peak as strong
						secondPeakIsStrong = true;
						// is it the same sign?
						if (ypeak2*ypeak > 0 )
						  {
							// peaks starting to compete. Strike compromise
							tp_comp = ((1.0-theta_c)*tp + (theta-theta_c)*xpeak2)/(1.0+theta-2*theta_c);
							tp = tp_comp; 
						  }
					  }
				  }
			  }
		  }

		peakTime = tp;
		return tp;
	  }
  }


  /**
	 Get wavelet absolute peak
  */
  double getCoarseAbsPeak()
  {
	int j;
	double max;

	max=0;
	for (j=0; j<valuesCoarse.length; j++) 
	  {
		if (Math.abs(valuesCoarse[j]) > max)
		  {
			max = Math.abs(valuesCoarse[j]);
		  }
	  }
	return max;
  }

  // useful debugging function. Splines can do tricky things!
  void debugPrintFineWavelet()
  {
	int i,n;
	double t,tmin,tmax;

	tmin=tCoarse[0];
	tmax=tCoarse[tCoarse.length-1];
	n =(int)Math.ceil((tmax - tmin)/(dtCoarse/10.0));
	double[] y = new double[n];
	double[] times = new double[n];
	for (i=0,t=tmin; (t<=tmax) && (i<n); i++,t+=dtCoarse/10.0)
	  {
		times[i] = t;
		y[i] = tSpline.splint(t,i>0);
	  }
	 Util.printfileZO("fine-wavelet.txt",times,y,n);
  }

  // useful debugging function. Splines can do tricky things!
  void debugPrintCoarseWavelet()
  {
	Util.printfileZO("coarse-wavelet.txt",tCoarse,valuesCoarse,tCoarse.length);
  }


  short getFinalIndex()
  {
	return (short)finalFineSample;
  }
  /**
	   precursor, or lead time of wavelet in seconds, 
	   (i.e. extent of influence for times > t_reflection)
	   Valid only after dt field has been set.
  */
  double getDTprecursor()
  {
	return  (-firstFineSample*delta_t);
  }
  /**
	 coda length of wavelet in seconds;  (i.e. extent of influence for times < t_reflection)
	 Valid only after dt field has been set.
  */
  double getDTcoda()
  {
	return (double)(finalFineSample)*delta_t;
  }

  /**
	 number of values in the values[] array
   */
  short getLength()
  {
	return (short)(finalFineSample-firstFineSample+1);
  }
  /**
	 number of values in the values[1..][] array
   */
  short getStretchedLength()
  {
	return (short)(finalStretchedFineSample-firstStretchedFineSample+1);
  }

 /**
	 Write out 
	 SU version of current wavelet, on fine scale.
	 If addExtraPadding==true, extra padding will be added so overall window
	 spans tMin to tMax, otherwise these parameters are ignored.
  */

  void writeOutSUWavelet(String fileName, boolean addSUSuffix,boolean addExtraPadding, 
						 double tMin, double tMax,boolean append,int mark_value,boolean haveTseismicErrors)
  {
	int j,k,kk,kStart,n;
	int nPrecursorPad=0,nCodaPad=0;
	double val,tStart,f1_tmp;
	double[] S;
	double[] _values;
	String filename=null;
	FileOutputStream out=null;
	int synchIndx;

	// determine whether "offset" style wavelets can be written for this problem (needs all wells synchronized)
	// if not, synchIndx=-1
	synchIndx = _WE.m.getCommonTregSynchIndx();

	// set requisite fields
	f1_tmp = f1[0];
	_WE.trW.set("f1",(float)f1_tmp);
	_WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp))); // in ms. Best approx to f1
	_WE.trW.set("ns",getLength());	
	_WE.trW.set("dt",dt);
	_WE.trW.set("mark",(short)mark_value);

	if (addExtraPadding)
	  {
		nPrecursorPad = Math.max(0,(int)Math.rint((-precursorDuration-tMin)/delta_t));
		nCodaPad = Math.max(0,(int)Math.rint((tMax-codaDuration)/delta_t));
		if (nPrecursorPad>0) 
		  {
			f1_tmp = f1[0] - nPrecursorPad*delta_t;
			_WE.trW.set("f1",(float)f1_tmp);
			_WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp)));
		  }
		_WE.trW.set("ns",(short)(getLength()+nPrecursorPad+nCodaPad));
	  }
	// suitable file for data
	if (addSUSuffix) 
	  {
		filename = Util.fileBasename(fileName).replaceAll(".su","").concat(".su");
	  }
	else filename = new String(Util.fileBasename(fileName));

	_WE.EH.println("Writing ML wavelet to "+filename);
	_WE.lastWaveletFilenames[0] = filename; // remember !

	try 
	  {
		out = new FileOutputStream(filename,append);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open" + filename);
	  }
	_values = values[0];


	// precursor zero padding if required
	for (k=0; k<nPrecursorPad; k++) _WE.trW.set_data_value(0f,k);
	// precursor (values[] is in wrap-around order)
	for (kk=finalFineSample+1; kk<_values.length; kk++,k++)
	  {
		_WE.trW.set_data_value((float)_values[kk],k);
	  }
	// coda
	for (kk=0; kk<=finalFineSample; kk++,k++)
	  {
		_WE.trW.set_data_value((float)_values[kk],k);
	  }
	// coda zero padding
	for (kk=0; kk<nCodaPad; kk++,k++) _WE.trW.set_data_value(0f,k);
	// write it
	_WE.suwriter.fputtr(out,_WE.trW);
	try
	  {
		_WE.suwriter.flush(out); // must do this!
		out.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close" + filename);
	  }

	// Write "offset" wavelets if possible & desirable
	if ((synchIndx != -1) && haveTseismicErrors)
	  {
		// offset wavelets for each stack
		for (j=0; j<_WE.m.nStacks; j++)
		  {
			if (useStretchedModel && j>0) break; // only near stack if stretching is happening

			String oFilename = Util.fileBasename(filename).replaceAll(".su",".stack."+(j+1)+".offset.su");
			_WE.EH.println("Writing ML \"offset\" wavelet to "+oFilename);
			try 
			  {
				out = new FileOutputStream(oFilename,append);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" + oFilename);
			  }
			// put time-shifts in header words
			f1_tmp = f1[j] + _WE.seismicInfo[synchIndx][_WE.synchroniseStackTregistrationErrors ? 0 : j].systematicTProp.value;
			_WE.trW.set("f1",(float)f1_tmp);
			// Zero delrt, so folks know it is a funky number in f1
			// _WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp))); // in ms. Best approx to f1
			_WE.trW.set("delrt",(short)0);
			// write it
			_WE.suwriter.fputtr(out,_WE.trW);
			try
			  {
				_WE.suwriter.flush(out); // must do this!
				out.close();
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close" + filename);
			  }
		  }
	  }

	


	// Far wavelet if stretch used
	if (useStretchedModel)
	  {
		for (j=1; j<_WE.m.nStacks; j++)
		  {
			// set requisite fields
			f1_tmp = f1[j];
			_WE.trW.set("f1",(float)f1_tmp);
			_WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp))); // in ms. Best approx to f1
			_WE.trW.set("ns",getStretchedLength());	
			_WE.trW.set("dt",dt);
			_WE.trW.set("mark",(short)mark_value);

			if (addExtraPadding)
			  {
				nPrecursorPad = Math.max(0,(int)Math.rint((-precursorStretchedDuration-tMin)/delta_t));
				nCodaPad = Math.max(0,(int)Math.rint((tMax-codaStretchedDuration)/delta_t));
				if (nPrecursorPad>0) 
				  {
					f1_tmp = f1[j] - nPrecursorPad*delta_t;
					_WE.trW.set("f1",(float)f1_tmp);
					_WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp))); // in ms. Best approx to f1
				  }
				_WE.trW.set("ns",(short)(getStretchedLength()+nPrecursorPad+nCodaPad));			
			  }
			// suitable file for data
			filename = fileName.replaceAll(".su","").concat(".stack."+(j+1)+".stretched.su");

			_WE.EH.println("Writing ML wavelet to "+filename);
			_WE.lastWaveletFilenames[j] = filename;

			try 
			  {
				out = new FileOutputStream(filename,append);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" + filename);
			  }
			// precursor zero padding if required
			for (k=0; k<nPrecursorPad; k++) _WE.trW.set_data_value(0f,k);
			// precursor (values[] is in wrap-around order)
			_values = values[j];
			for (kk=finalStretchedFineSample+1; kk<_values.length; kk++,k++)
			  {
				_WE.trW.set_data_value((float)_values[kk],k);
			  }
			// coda
			for (kk=0; kk<=finalStretchedFineSample; kk++,k++)
			  {
				_WE.trW.set_data_value((float)_values[kk],k);
			  }
			// coda zero padding
			for (kk=0; kk<nCodaPad; kk++,k++) _WE.trW.set_data_value(0f,k);
			// write it
			_WE.suwriter.fputtr(out,_WE.trW);

			try
			  {
				_WE.suwriter.flush(out); // must do this!
				out.close();
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close" + filename);
			  }

			// Write "offset" wavelets if possible
			if ((synchIndx != -1) && haveTseismicErrors)
			  {
				// offset wavelets for each stack
				String oFilename = Util.fileBasename(filename).replaceAll(".su",".offset.su");
				_WE.EH.println("Writing ML \"offset\" wavelet to "+oFilename);
				try 
				  {
					out = new FileOutputStream(oFilename,append);
				  } 
				catch (IOException e) 
				  {
					throw new RuntimeException("couldn't open" + oFilename);
				  }
				// put time-shifts in header words
				f1_tmp = f1[j] + _WE.seismicInfo[synchIndx][_WE.synchroniseStackTregistrationErrors ? 0 : j].systematicTProp.value;
				_WE.trW.set("f1",(float)f1_tmp);
				// Zero delrt, so folks know it is a funky number in f1
				// _WE.trW.set("delrt",(short)(Math.rint(1e3*f1_tmp))); // in ms. Best approx to f1
				_WE.trW.set("delrt",(short)0);
				// write it
				_WE.suwriter.fputtr(out,_WE.trW);
				try
				  {
					_WE.suwriter.flush(out); // must do this!
					out.close();
				  }
				catch (IOException e) 
				  {
					throw new RuntimeException("couldn't close" + filename);
				  }
			  }
		  }
	  }
  }



  // Mock broadband wavelet, for experiments. Two types, hand-tuned from 
  // experiments (private JG: see WaveletFilters.nb)
  private void makeBBWavelet(double f_at_band_edge,double suggested_sampling_dt)
  {
	int N,i,j,ii,k;
	double t,norm,tmax;
	double[] _values;

	// splines on t-scale with bandedge of 1.
	// mildly broadband
	double[] tsp1= {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60};
	double[] ysp1 ={1., 0.0776222, -0.160865, -0.0279797, -0.0651356, -0.0623286, 
				   -0.0504867, -0.0451806, -0.0416235, -0.0354069, -0.02853, -0.0231756, 
				   -0.018456, -0.0134745, -0.00900451, -0.00564087};
	// very broadband
	double[] tsp2 = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68,
					 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 
					 132, 136, 140, 144, 148};
	double[] ysp2 = {1., 0.0171062, -0.0545585, 0.00321552, -0.0324364, -0.0176984, 
					 -0.0187235, -0.0214439, -0.0192669, -0.0186839, -0.019309, 
					 -0.0184234, -0.0177171, -0.0176836, -0.0170576, -0.0163361, 
					 -0.0160016, -0.01542, -0.0146963, -0.0142016, -0.0136133, -0.0128953, 
					 -0.012315, -0.0117159, -0.0110148, -0.0103972, -0.00979922, 
					 -0.00912863, -0.00850938, -0.00792913, -0.0073031, -0.00670956, 
					 -0.00616406, -0.00559539, -0.00504861, -0.00455279, -0.00405187, 
					 -0.00356722};
	
	double[] tsp,ysp;
	if (true)
	  {
		// mild
		tsp=tsp1; ysp=ysp1;
	  }
	else
	  {
		// very
		tsp=tsp2; ysp=ysp2;
	  }
	Spline sp = new Spline(tsp,ysp,tsp.length, Double.MAX_VALUE, Double.MAX_VALUE,_EH);

	N = tsp.length-1;
	tmax = tsp[N-1];

	_EH.println("Making BB wavelet: edge freq. = %10.6g",f_at_band_edge); 
	dtCoarse = 1/(6.3*f_at_band_edge); // specific to spline form above

	// compute suitable finer-scale canonical sampling rate
	
	// band edge Nyquist rate
	delta_t = 0.5/f_at_band_edge;

	// find canonical sampling time just less than this (1ms,2ms,4ms etc)
	if (delta_t<FundamentalConstants.BASE_SCALING_TIME)
	  {
		throw new RuntimeException("Desired synthetic seismic has unreasonably high bandwidth");
	  }
	for (t=FundamentalConstants.BASE_SCALING_TIME ;; t*=2.0)
	  {
		if (t>delta_t) break;
		if (t>1000*FundamentalConstants.BASE_SCALING_TIME) 	
		  {
			throw new RuntimeException("Desired synthetic seismic has unreasonably low bandwidth");
		  }
	  }
	// desired sampling time is half of last attempt, since it failed
	delta_t = 0.5*t;
	// compare to users suggested time
	if (delta_t < 0.99 * suggested_sampling_dt)
	  {
		_EH.printf("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is too long "
						+"for bandwidth desired. I suggest "+(1000*delta_t)+"ms.%n"
						+"I will use " + (1000*delta_t)+"ms despite what you suggest.%n");
		
	  }
	else if ((delta_t > 1.01 * suggested_sampling_dt) &&  (delta_t < 4 * suggested_sampling_dt))
	  {
		_EH.printf("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is needlessly "
						+"short for bandwidth desired. I suggest "+(1000*delta_t)+"ms%n"
						+"I'll do what you say anyway.%n");
		// override sensible value with users idea
		delta_t = suggested_sampling_dt;
	  }
	else if (delta_t > 4 * suggested_sampling_dt)
	  {
		throw new RuntimeException("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is "
								   +"absurdly short for bandwidth desired. I suggest "+(1000*delta_t)+"ms.");
	  }
	else
	  {
		// a good suggestion
		delta_t = suggested_sampling_dt;
	  }
	// sample time in mu_s
	dt = (short)(1e6 * delta_t);

	// wavelet is dead (to ~10^-3) by t = 6 in spline, so
	finalCoarseSample = N;
	firstCoarseSample = -N;
	// x values for abscissa in spline
	tCoarse = new double[2*N+1];
	valuesCoarse = (double[])tCoarse.clone();

	// store Coarse wavelet in ordinary way............................
	for (j=0,i=-N,norm=0; i<=N; i++,j++)
	  {
		ii = (i<0) ? -i : i;
		tCoarse[j] = (((i<0)?-1:1) * tsp[ii]) * dtCoarse;
		valuesCoarse[j] = ysp[ii];
		norm += valuesCoarse[j]*valuesCoarse[j];
	  }
	norm=Math.sqrt(norm);
	// make Spline interpolator
	tSpline = new Spline(tCoarse,valuesCoarse,tCoarse.length,0.0,0.0,_EH);

	// make fine wavelet by cubic interpolation
	finalFineSample = (int)(dtCoarse*tsp[N]/delta_t);
	firstFineSample = -finalFineSample;

	f1 = new double[_WE.m.nStacks];
	f1[0] = firstFineSample*delta_t;

	// Multi wavelet entries
	values = new double[_WE.m.nStacks][];
	values[0] = new double[2*finalFineSample+1];
	// System.out.println("");
	_values = values[0];
	for (i=0;i<=finalFineSample;i++)
	  {
		_values[i] = tSpline.splint(i*delta_t,true);
	  }
	// wraparound bit
	for (j=2*finalFineSample,i=1;i<=finalFineSample;i++,j--)
	  {
		_values[j] = _values[i];
	  }



	// stretched wavelet stuff

	setUpStretching();
	if (_WE.useStretchedModel)
	  {
		for (k=0; k < (_WE.m.nStacks-1); k++)
		  {
			_WE.m.farWavStretch[k] = _WE.m.farStackWaveletTimeStretchMean;
			_WE.m.farWavAmp[k] = _WE.m.farStackWaveletAmplitudeMean;
		  }
		computeFarWaveletValues();
	  }
  }

  /** 
	  Make Ricker wavelet for near wavelet.
	  Calculations here assume spectral form W(f) ~ f^2 exp(-f^2 lambda^2 / 2)
	  or, in time,   W(t) ~ (1-(2 pi t/lambda)^2) * exp(-(2 pi t)^2 / 2 lambda^2).
	  From frequency spectrum, lambda is defined from bandedge 
	  by lambda = 3/freq_at_bandedge.
	  From suggested freq_at_bandedge, we compute lambda, a sensible sampling rate dt, and
	  make the wavelet etc.
	  Bandedge relates to Peak frequency by fBE/f_peak=3/sqrt(2)
  */
  private void makeRickerWavelet(double f_at_band_edge,double suggested_sampling_dt)
  {
	int i,j,k;
	double spacing,lambda,f_peak;
	double x,ls,t;
	double norm;
	double[] _values;

	_EH.printf("%nCreating synthetic wavelet%n");

	// Sampling rate from magic Rule
	f_peak = f_at_band_edge/2.12;    // standard ratio of f_peak to fBandedge. 14% power point on tail of spectrum

	_EH.println("Ricker peak freq. = %10.6g",f_peak); // FIXME may want to remove
	dtCoarse = SAMPLE_TIME_MAGIC_CONSTANT/f_peak;

	lambda = 3.0/f_at_band_edge;

	// compute suitable finer-scale canonical sampling rate
	
	// band edge Nyquist rate
	delta_t = 0.5/f_at_band_edge;

	// find canonical sampling time just less than this (1ms,2ms,4ms etc)
	if (delta_t<FundamentalConstants.BASE_SCALING_TIME)
	  {
		throw new RuntimeException("Desired synthetic seismic has unreasonably high bandwidth");
	  }
	for (t=FundamentalConstants.BASE_SCALING_TIME ;; t*=2.0)
	  {
		if (t>delta_t) break;
		if (t>1000*FundamentalConstants.BASE_SCALING_TIME) 	
		  {
			throw new RuntimeException("Desired synthetic seismic has unreasonably low bandwidth");
		  }
	  }
	// desired sampling time is half of last attempt, since it failed
	delta_t = 0.5*t;
	// compare to users suggested time
	if (delta_t < 0.99 * suggested_sampling_dt)
	  {
		_EH.printf("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is too long "
						+"for bandwidth desired. I suggest "+(1000*delta_t)+"ms.%n"
						+"I will use " + (1000*delta_t)+"ms despite what you suggest.%n");
		
	  }
	else if ((delta_t > 1.01 * suggested_sampling_dt) &&  (delta_t < 4 * suggested_sampling_dt))
	  {
		_EH.printf("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is needlessly "
						+"short for bandwidth desired. I suggest "+(1000*delta_t)+"ms%n"
						+"I'll do what you say anyway.%n");
		// override sensible value with users idea
		delta_t = suggested_sampling_dt;
	  }
	else if (delta_t > 4 * suggested_sampling_dt)
	  {
		throw new RuntimeException("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is "
								   +"absurdly short for bandwidth desired. I suggest "+(1000*delta_t)+"ms.");
	  }
	else
	  {
		// a good suggestion
		delta_t = suggested_sampling_dt;
	  }
	// sample time in mu_s
	dt = (short)(1e6 * delta_t);

	// wavelet is dead (to ~10^-3) by t = 0.7 * lambda, so
	finalCoarseSample = (int)Math.ceil(0.7*lambda/dtCoarse);
	firstCoarseSample = -finalCoarseSample;
	// x values for abscissa in spline
	tCoarse = new double[2*finalCoarseSample+1];
	valuesCoarse = new double[2*finalCoarseSample+1];

	// compute coarse wavelet
	spacing = lambda/(2.0*Math.PI*dtCoarse);
	ls = spacing * spacing;

	// store Coarse wavelet in ordinary way............................
	//   System.out.println("wavelets:");
	for (j=0,i=firstCoarseSample,norm=0;i<=finalCoarseSample;i++,j++)
	  {
		x = i * i;
		tCoarse[j] = i*dtCoarse;
		valuesCoarse[j] = (1.0-x/ls)*Math.exp(-0.5*x/ls);
		norm += valuesCoarse[j]*valuesCoarse[j];
	  }
	norm=Math.sqrt(norm);
	/*
	System.out.println("Wavelet coeffs: (t, amp)");
	for (i=0;i<tCoarse.length;i++)
	  {
		// Omit normalisation, since this doesn't matter that much, and
		// it's useful if we know the central coefficient is exactly 1
		// Normalisation:
		//valuesCoarse[i]/=norm;
		System.out.println(tCoarse[i]+" "+valuesCoarse[i]);
	  }
	*/
	// make Spline interpolator
	tSpline = new Spline(tCoarse,valuesCoarse,tCoarse.length,0.0,0.0,_EH);

	// make fine wavelet by cubic interpolation
	finalFineSample = (int)(finalCoarseSample*dtCoarse/delta_t);
	firstFineSample=-finalFineSample;

	f1 = new double[_WE.m.nStacks];
	f1[0] = firstFineSample*delta_t;

	// Multi wavelet entries
	values = new double[_WE.m.nStacks][];
	values[0] = new double[2*finalFineSample+1];
	// System.out.println("");
	_values = values[0];
	for (i=0;i<=finalFineSample;i++)
	  {
		_values[i] = tSpline.splint(i*delta_t,true);
	  }
	// wraparound bit
	for (j=2*finalFineSample,i=1;i<=finalFineSample;i++,j--)
	  {
		_values[j] = _values[i];
	  }



	// stretched wavelet stuff

	setUpStretching();
	if (_WE.useStretchedModel)
	  {
		for (k=0; k < (_WE.m.nStacks-1); k++)
		  {
			_WE.m.farWavStretch[k] = _WE.m.farStackWaveletTimeStretchMean;
			_WE.m.farWavAmp[k] = _WE.m.farStackWaveletAmplitudeMean;
		  }
		computeFarWaveletValues();
	  }
  }


  /**
	 a method to estimate suitable wavelet independent sampling rate from
	 sample data
  */
  public static double estimateWaveletIndependentDT(double[] data,double dt_data, 
													int verbosity, ErrorHandler EH)
  {
	// Take FFT of data to estimate best sampling rate
	int i,m,n;
	double mOn2;
	double[] padded_data;
	double delta_t_indep_samples,f_peak;
	double i_edge,f_edge;
	

	m=n=data.length;

	// FFT power of 2 stuff
	int ng=0;
	while ((n >>= 1) > 0 ) ng++;
	// In general data.length not a power of 2, so find the right n = 2^ng  bigger than data.length
	if (data.length != (1<<ng))
	  {
		// data.length not a power of 2
		n=(1<<(ng+1));
	  }
	else
	  n = data.length; // it was already a power of 2
	// now pad out data
	padded_data = new double[n];
	// +ve time part of data
	System.arraycopy(data,0,padded_data,0,m);
	// Welch window data to taper nicely before FFT
	mOn2=m/2.0;
	for (i=0;i<m;i++)padded_data[i]*=(1.0-sqr((i-mOn2)/mOn2));

	// get FFT
	FFT.realftZO(padded_data,n,1);

	// Method to extract spectral peak from first moment (more robust)
	int i_max=0;
	double spec, max=0.0,Esum=0.0;
	double fft[] = new double[n/2];
	double freqs[] = new double[n/2];
	double i_peak, num=0,den=0;
	double freq_den = (2.0*(double)n * dt_data);
	// Util.printfileZO("Padded_FFT.txt",padded_data,padded_data.length);
	for (i=2;i<=n-2;i+=2)
	  {
		spec = sqr(padded_data[i])+sqr(padded_data[i+1]);
		if (spec>=max)
		  {
			max = spec;
			i_max=i;   // peak occurs at this index
		  }
		fft[i/2]=spec;
		freqs[i/2]=(double)i/freq_den;        // freqs[j] = 2 j /(2 n dt) = j/(n dt), j=1...n/2-1
		num += spec*i;  // first moment
		den += spec;    // normalisation
		Esum += spec;
	  }
	// So, from padded_data[0...n-1],
	// the DFT frequencies are f_j = j/(n dt_data), j=-n/2...n/2  (NR, equation 12.1.5)
	// The peak-frequency f_j occurs at j = i_max/2, since after the FFT,
	// padded_data[] alternates Re, Im parts:  
	// padded_data[2,3] = {Re(f_1),Im(f_1)}, padded_data[4,5] = {Re(f_2),Im(f_2)}, etc. 
	
	// Alternatively, define i_peak from first moment plus magic Ricker rule
	// (less noisy approximation to i_max above: should be comparable)
	i_peak = (num/den)/MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER;

	// Chunk to print out spectrum. It's typically quite noisy
	// OLD Util.printfileZO("FFT_data_spectrum.txt",freqs,fft,Math.min(4*i_max,n/4));
	Util.printfileZO("FFT_data_spectrum.txt",freqs,fft,n/2);

	fft=null;   // GC
	freqs=null; // GC

	f_peak = (i_peak/2.0)/((double)n * dt_data); // OLD (wrong) i_peak/((double)n * dt_data);
	EH.println("Inferred (1st moment-based) data peak frequency=%1g Hz",f_peak);

	// Where the tail cutoff has been placed if you graph "FFT_data_spectrum.txt"
	// This is the Nyquist freq associated with delta_t_indep_samples below

	// Estimate of band edge is f_peak/(2*SAMPLE_TIME_MAGIC_CONSTANT)

	EH.println("Inferred band-edge Ricker-rule cutoff estimated (from peak) at effective frequency %1g ",
			   f_peak/(2.0*SAMPLE_TIME_MAGIC_CONSTANT));

	// Here's the 95% point doing a trapezoid-based integral
	double sum=0;
	i_edge = i_peak/(2.0*SAMPLE_TIME_MAGIC_CONSTANT);  // make sure this has a sane default value

	Esum *= BAND_EDGE_FRACTION; // 95% Energy point
	for (i=2;i<=n-2;i+=2)
	  {
		spec = sqr(padded_data[i])+sqr(padded_data[i+1]);
		sum += spec;
		if (sum > Esum)
		  {
			i_edge = (spec>0) ? ((i-2)+ (Esum - (sum-spec))/spec) : (i-2);
			break;
		  }
	  }
	double f_edge_trapez = i_edge/(2.0*(double)n * dt_data);
	EH.println("Inferred band-edge cutoff from spectral integral=%1g Hz",f_edge_trapez);

	// But let's use the sampling time corresponding to the maximum of the frequencies
	// f_peak/SAMPLE_TIME_MAGIC_CONSTANT, f_edge_trapez

	// Set error-sampling time of all synthetics
	delta_t_indep_samples = Math.min(SAMPLE_TIME_MAGIC_CONSTANT/f_peak, 1.0/(2.0*f_edge_trapez));


	padded_data=null;

	return delta_t_indep_samples;
  }

  /**
	 Compute phase spectrum of current wavelet. Print to a file if
	 filename != null.
	 Also, print duplicate copy over range [phaseConstraintIndxStart,phaseConstraintIndxEnd]
	 if these values are non-zero
  */
  void computePhaseSpectrum(String filename, int phaseConstraintIndxStart, int phaseConstraintIndxEnd)
  {
	int i,j,k,kk,m,n;
	//double mOn2;
	double[] padded_data;
	double[] _values;

	_values = values[0];
	m = _values.length;
	n = FFT.getNearestPowerOfTwo(_values.length);

	// now pad out data
	padded_data = new double[n];
	// precursor (_values[] is in wrap-around order)
	for (k=n-1,kk=m-1; kk>=(m-1+firstFineSample); kk--,k--)
	  {
		padded_data[k] = _values[kk];
	  }
	// coda
	for (k=0,kk=0; kk<=finalFineSample; kk++,k++)
	  {
		padded_data[k] = _values[kk];
	  }
	// Welch window data to taper nicely before FFT
	// mOn2 = m/2.0;
	// for (i=0;i<m;i++) padded_data[i]*=(1.0-sqr((i-mOn2)/mOn2));
	// get FFT

	// Util.printfileZO("Padded_wavelet.txt",padded_data,padded_data.length);

	FFT.realftZO(padded_data,n,1);

	double freqs[] = new double[n/2];
	if (phase==null || (phase.length!=n)) phase = new double[n/2];

	//for (i=0,j=0; i<=n/2; i+=2,j++)
	for (i=0,j=0; i<=n-2; i+=2,j++)
	  {
		double spec = sqr(padded_data[i])+sqr(padded_data[i+1]);
		// fft[i/2=j]=spec;
		freqs[j] = phaseVectorFrequencySpacing*(double)j;
		phase[j] = Math.toDegrees(Math.atan2(padded_data[i+1],padded_data[i]));
		phaseRe[j] = padded_data[i];
		phaseIm[j] = padded_data[i+1];
	  }
	// Chunk to print out phase spectrum. It's typically quite noisy
	if ((phaseConstraintIndxStart!=0) && (phaseConstraintIndxEnd!=0))
	  {
		//for (i=0;i<Math.min(phase.length,phaseConstraintIndxStart);i++) phase[i]=0;
		//for (i=phase.length-1;i>Math.max(0,phaseConstraintIndxEnd);i--) phase[i]=0;
	  }
	if (filename != null)
	  {
		//OLD Util.printfileZO(filename,freqs,phase,n/4);
		Util.printfileZO(filename,freqs,phase,n/2);
	  }
	freqs=null; // GC
  }

  void computeSpectrum(String filename)
  {
	int i,j,k,kk,m,n,Nm;
	double[] padded_data;
	double[] _values;
	double S;
	
	_values = values[0];
	m = _values.length;
	n = FFT.getNearestPowerOfTwo(_values.length);

	// now pad out data
	padded_data = new double[n];
	// precursor (_values[] is in wrap-around order)
	for (k=n-1,kk=m-1; kk>=(m-1+firstFineSample); kk--,k--) padded_data[k] = _values[kk];
	// coda
	for (k=0,kk=0; kk<=finalFineSample; kk++,k++) padded_data[k] = _values[kk];
	// get FFT

	// Util.printfileZO("Padded_wavelet.txt",padded_data,padded_data.length);

	FFT.realftZO(padded_data,n,1);

	Nm = n/2-1;

	double freqs[] = new double[n/2];
	double spec[] = new double[n/2];

	// compute pwr up & including  Nyquist, i=0...Nm
	// i=0
	freqs[0] = 0;
	spec[0] = sqr(padded_data[0]);
	// i=1...Nm-1
	for (i=1; i<Nm; i++)
	  {
		S = sqr(padded_data[2*i]) + sqr(padded_data[2*i+1]);	  
		freqs[i] = phaseVectorFrequencySpacing*(double)i;
		spec[i] = S;
	  }
	// i=Nm
	freqs[Nm] = phaseVectorFrequencySpacing*(double)Nm;
	spec[Nm] = sqr(padded_data[1]);
  

	// Ampltude power exponent:  W_(f) ~ f^p. Estimate p
	double Aslope[] = new double[n/2];
	for (i=0; i<(Nm-1); i++)
	  {
		if(i==0)
		  {
			Aslope[i] = ((log(spec[i+1])-log(spec[i]))/phaseVectorFrequencySpacing)
			  /log(freqs[i]+.1);
			Aslope[i]=0; // not reliable
		  }
		else
		  {
			Aslope[i] = ((log(spec[i+1])-log(spec[i-1]))/(2*phaseVectorFrequencySpacing))
			  /log(freqs[i]);
		  }
	  }
	Util.printfileZO("Pwr_slope.txt",freqs,Aslope,Nm-1);
	
	// Chunk to print out pwr spectrum. It's typically quite noisy
	if (filename != null)
	  {
		Util.printfileZO(filename,freqs,spec,Nm);
	  }
  	freqs=null; // GC
	spec=null;
  }
}
 
