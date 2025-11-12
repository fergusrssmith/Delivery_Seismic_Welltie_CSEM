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



package au.csiro.JamesGunning.seismicinverter;

import java.util.*;
import java.io.*;
import cern.colt.list.*;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.util.*;





/**
   Class to attempt to find spike positions by deconvolution and threshold detection.
*/	 
public final class EventDetector
{
  Model _m;
  ModelSampler _MS;                          // reference to _m.MS

  /*
	Recursive event enumeration members
  */
  private int[] pos,best_pos;
  private double[] tE;
  int Np,Ne,iE,iEbest=0;
  double bestNegLogPost;

  /*
	  If this is set to 
	  * 0, layer tops can pile up on the same event. making the 
	  layer combinations more numerous. This is good if the event detection finds fewer events than there
	  are layers.

	  * If 1, the next layer time will fall 
	  no closer than the next detected event. Good in that there are fewer combinations, but
	  we may miss a pinchout...
	*/
  public static int INCREMENT_FOR_EVENTS = 0;

  /**
	 maximum number of combinations to examine. If things are worse than this, 
	 the class will abandon and return null.
  */
  public static int MAX_COMBINATIONS = 10000;

  /**
	 Relative noise level used in deconvolution padding of wavelet spectrum
  */
  public static double DECON_EPS=0.1;

  /*
	Debugging members
  */
  static boolean DEBUG_DUMP_SYNTH_SEIS = false;
  static boolean DEBUG = false;
  static String SYNTH_NEAR_FILENAME = "event_detection_synth.su";

  private SuWriter suwriter;
  private SuSegyFormat tr_near_synth;
  private FileOutputStream near_synth_out_strm;

  /**
	 Constructor.
	 nonTimeProperties must be an ArrayList of PropertyMapper objects, pointing to
	 suitable layer properties to optimise over
  */ 

  public EventDetector(Model m)
  {
	this._m = m;
	this._MS = _m.MS;
  }

  /**
	 Returns an array of best model times t[0..._m.n_active_layers]
	 as determined by deconvolution and thresholding. If the number
	 of combinations exceeds MAX_COMBINATIONS, we give up and return null.
  */
 public double[] detectLayerTimes()
  {
	int i,ii,j,k,N,Ndata;
	int[] peakIndexes;


	//	Debug-dumper test

	if (_m._SI.verbosity>2)
	  {
		_m._SI.EH.printf("Starting event detection...%n");
	  }
	// set all parameters at their prior+isopach posterior means
	for (k=0;k<_MS.mu.size();k++)
	  {
		_MS.property_sample[k] = _MS.mu.getQuick(k);
		_MS.setValuesFromIndex(k,_MS.property_sample[k]);
	  }
	// order the times
	_MS.orderTimes();
	// generate effective impedance etc for all layers
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK


	// retrieve the near-stack seismic
	SuSegyFormat tr = _m._SI.tr_near;
	StackData stack_data = _m.SM.nearstack_data;

	double t_start_comparison = _m.MPD.MetaInfo[1].t.value - stack_data.getTabove();
	double t_end_comparison = _m.MPD.MetaInfo[_m.n_layers-1].t.value + stack_data.getTbelow();
	// Data. Traces may be delayed if in eg deep water, and delrt!=0 (delrt in ms)
	double t_true_start = ((double)tr.get_short("delrt"))*1e-3;
	double delta_t_true = ((double)tr.get_short("dt"))*1e-6;
	int ns_true = (int)tr.get_short("ns");

	int jstart = (int)Math.ceil((t_start_comparison - t_true_start)/delta_t_true);
	int jend = (int)Math.floor((t_end_comparison - t_true_start)/delta_t_true);

	Ndata = jend-jstart+1;
	N = FFT.getNearestPowerOfTwo(Ndata);
	if (N == Ndata) N *=2;
	double[] S = new double[N];
	// copy data
	for (i=0,j=jstart; j<=jend; j++,i++) S[i] = tr.get_data_value(j);

	// now detect events in S[]
	OneDSeismicConvolutionalModel.Wavelet wavelet = _m.SM.near_wavelet;

	// Set up padded wavelet just as in generateSyntheticTrace_kernel()
	double[] padded_wavelet = new double[N];
	int wl = wavelet.getLength();
	int wfi = wavelet.getFinalIndex();
	// +ve time part of wavelet
	System.arraycopy(wavelet.W,0,padded_wavelet,0,wfi+1);
	// -ve time part
	System.arraycopy(wavelet.W,wfi+1,
					 padded_wavelet, wfi + 1 + N - wl,
					 wl - 1 - wfi);

	// inject wavelet into complex array (imag bits will allocate to 0 )
	double[] complex_wavelet = new double[2*N];
	double[] cmplx_wavelet_fft_amp = new double[N];

	for (i=0;i<N;i++) complex_wavelet[2*i]=padded_wavelet[i];
	padded_wavelet=null;

	// Now get FFT of wavelet
	FFT.four1ZO(complex_wavelet,N, 1);

	if (DEBUG) Util.printfileZO("cmplx_wav_FFT.txt",complex_wavelet,2*N);

	int n=N;
	// zero pad data
    double[] complex_data = new double[2*N];
	double[] deconvolved_data = new double[2*N];

	// form data sequence with Welch window
	if (true)
	  {
		double wt;
		double Non2 = Ndata/2.0;
		for (i=0;i<Ndata;i++) 
		  {
			wt = 1.0 - MathUtil.sqr((i - Non2)/Non2);
			complex_data[2*i] = wt*S[i];
		  }
	  }
	else
	  {
		for (i=0; i<Ndata; i++) 
		  {
			complex_data[2*i] = S[i];
		  }
	  }
	if (DEBUG) Util.printfileZO("true.txt",S,N);
	if (DEBUG) Util.printfileZO("cmplx_seis.txt",complex_data,2*N);
	// Now get FFT of data
	FFT.four1ZO(complex_data,N, 1);
	if (DEBUG) Util.printfileZO("cmplx_seis_fft.txt",complex_data,2*n);

	/* now deconvolve, with bandlimiting */
	double D,Dmax=0;
	for (i=0,ii=1;i<N;i++,ii+=2) 
	  {
		D = complex_wavelet[ii-1]*complex_wavelet[ii-1]+complex_wavelet[ii]*complex_wavelet[ii];
		Dmax = Math.max(Dmax,D);
	  }
	for (i=0,ii=1;i<N;i++,ii+=2) 
	  {
		D = complex_wavelet[ii-1]*complex_wavelet[ii-1]+complex_wavelet[ii]*complex_wavelet[ii];
		if (D==0.0)
		  {
			throw new IllegalArgumentException("Deconvolving at response zero in EventDetector");
		  }
		if (D < DECON_EPS*Dmax) D=DECON_EPS*Dmax;
		if (true) 
		  {
			// deconvolve 
			// +ve frequencies
			deconvolved_data[ii-1] = (complex_wavelet[ii-1] * complex_data[ii-1] + 
									  complex_wavelet[ii] * complex_data[ii])/(D*N);// Re
			deconvolved_data[ii] =  (complex_wavelet[ii-1] * complex_data[ii] - 
									 complex_wavelet[ii] * complex_data[ii-1])/(D*N);// Im
		  }
		else
		  {
			// examine only effect of windowing and band-limiting
			deconvolved_data[ii-1] = complex_data[ii-1]/N;
			deconvolved_data[ii] = complex_data[ii]/N;
		  }
		// -ve ones. Should be conjugate
		if (i>0)
		  {
			deconvolved_data[2*N - ii + 1] = deconvolved_data[ii-1]; // Re
			deconvolved_data[2*N - ii + 2] = -deconvolved_data[ii];  // Im
		  }
		else 
		  {
			// might as well set the DC term to zero for this application
			deconvolved_data[0] = deconvolved_data[1] = 0.0;
			// exact copy....
			//deconvolved_data[0] = complex_data[1]/N;
			//deconvolved_data[1] = 0.0;
		  }
	  }
	if (DEBUG) Util.printfileZO("deconFFT.txt",deconvolved_data,2*N);
	// and invert FFT to get spikes
	FFT.four1ZO(deconvolved_data,N, -1);
	// map deconvolved image into its first half, since every second imaginary entry is zero anyway
	for (i=0,ii=0;i<N;i++,ii+=2) 
	  {
		deconvolved_data[i]=deconvolved_data[ii];
	  }
	if (DEBUG) Util.printfileZO("decon.txt",deconvolved_data,N);


	/* 
	   Do peak detection on deconvolved data to form amplitude-sorted list of
	   possible peak locations
	*/
	peakIndexes = formListOfLikelyLayerLocations(deconvolved_data,Ndata);
	// if no detectable peaks, give up
	if (peakIndexes == null) return null;

	/* 
	   Now loop over list of possible partitions of internal layer times placed at these events.
	   We must place the Np internal reflections at all possible
	   internal choices of Np events from the list of Ne events. This is most simply
	   done recursively
	*/
	Ne = peakIndexes.length;
	Np =_m.n_active_layers-1;
	pos = new int[Np];
	best_pos = new int[Np];

	
	

	// care with combinatorial size
	// default is use deplicate mode. Check this first.
	if (INCREMENT_FOR_EVENTS==0)
	  {
		// Number of combinations here is Binomial(Ne+Np-1, Np)
		if (cern.jet.math.Arithmetic.binomial(Ne + Np - 1, Np) > MAX_COMBINATIONS)
		  {
			// first, try switching to non-duplicate mode (less combinations)
			INCREMENT_FOR_EVENTS = 1;
		  }
	  }
	if (INCREMENT_FOR_EVENTS==1)
	  {
		// This is a standard combinations problem. Goes as Binomial(Ne,Np)
		if (cern.jet.math.Arithmetic.binomial(Ne, Np) > MAX_COMBINATIONS)
		  {
			// too big. Give up by returning null
			return null;
		  }
		if (Np>Ne)
		  {
			// cannot partition with incrementing
			if (_m._SI.verbosity>2)
			  {
				_m._SI.EH.println("Event detector warning: fewer detected events than layers, skipping event detection");
			  }
			return null;
		  }
	  }

	if (Ne < Np && _m._SI.verbosity>2)
	  {
		_m._SI.EH.println("Event detector warning: fewer detected events than layers");
	  }

	// Compute associated event times for these Ne events.
	// Natural order is in descending
	// amplitude size, so these are not in time order
	tE = new double[Ne];
	for (i = 0; i < Ne; i++) 
	  {
		tE[i] = (peakIndexes[i]+jstart)*delta_t_true + t_true_start;
	  }
	// sort into ascending order
	Arrays.sort(tE);

	// Set up recursive enumeration loop
	bestNegLogPost = Double.MAX_VALUE;
	iE=0;
	if (DEBUG_DUMP_SYNTH_SEIS) initialiseDebugSyntheticSeismicOutput();

	// Now recurse
	// loop over combinations in array
	ArrayList posIndexes = Partitions.createPartitionList(Ne,Np,(INCREMENT_FOR_EVENTS==0));
	// impossible case
	if (posIndexes==null) return null;


	for (i = 0; i < posIndexes.size(); i++) 
	  {
		evaluateCombo((int[])posIndexes.get(i));
	  }
	posIndexes = null;


	// On completion, restore best times
	setTimes(best_pos);
	if (_m._SI.verbosity>2)
	  {
		_m._SI.EH.print("Event detection best times for combo "+iEbest+": ");
		for (i=0; i<Np;i++) _m._SI.EH.printf("%7.4f ",tE[best_pos[i]]);
		_m._SI.EH.println("");
	  }
	if (DEBUG_DUMP_SYNTH_SEIS) finalizeDebugSyntheticSeismicOutput();


	// last call to setTimes will have restored optimal times in
	// _MS.property_sample, so just return this
	if (false)
	  {
		double[] optTimes = new double[_m.n_active_layers+1];
		System.arraycopy(_MS.property_sample,0, optTimes,0,_m.n_active_layers+1);	
		return optTimes;
	  }
	else
	  {
		return _MS.property_sample;
	  }
  }


  private void evaluateCombo(int[] position)
  {
	int i;

	// combination constructed.
	// evaluate raw likelihood from associated combination of event times

	if (_m._SI.verbosity>4) 
	  {
		_m._SI.EH.printf("Combo %5d: ",iE);
		for (i=0; i<Np;i++) _m._SI.EH.printf("%3d ",position[i]);
		for (i=0; i<Np;i++) _m._SI.EH.printf("%7.4f ",tE[position[i]]);
	  }
	setTimes(position);

	// evaluate likelihood and store best params
	double val = _MS.debugComputeNegLogPosterior();
		
	if (_m._SI.verbosity>4) 
	  {
		_m._SI.EH.printf("%7.4f%n",val);
	  }
	if (val < bestNegLogPost)
	  {
		// store state
		bestNegLogPost = val;
		System.arraycopy(position,0,best_pos,0,Np);
		iEbest = iE;
	  }
	if (DEBUG_DUMP_SYNTH_SEIS)
	  {
		debugDumpSynthetics(iE);
	  }
	iE++;
  }




  /*
	function to set internal layer times according to the sequence of locations
	tE[position[0]], tE[position[1]], ....
  */

  private void setTimes(int[] position)
  {
	int i;
	// set all parameters at their prior+isopach posterior means
	// set top and base times
	if (!true)
	  {
		// from prior
		for (i=0; i<=(_m.n_active_layers); i+=_m.n_active_layers)
		  {
			_MS.property_sample[i] = _MS.mu.getQuick(i);
			_MS.setValuesFromIndex(i,_MS.property_sample[i]);
		  }
	  }
	else
	  {
		// from event first and last
		_MS.property_sample[0] = Math.min(tE[position[0]], _MS.mu.getQuick(0));
		_MS.setValuesFromIndex(0,_MS.property_sample[0]);
		_MS.property_sample[_m.n_active_layers] = Math.max(tE[position[Np-1]], _MS.mu.getQuick(_m.n_active_layers));
		_MS.setValuesFromIndex(_m.n_active_layers,_MS.property_sample[_m.n_active_layers]);
	  }
	// set internal layer times as per this event combination
	for (i=1;i<=Np;i++)
	  {
		_MS.property_sample[i] = tE[position[i-1]];
		_MS.setValuesFromIndex(i,_MS.property_sample[i]);
	  }
	// order the times
	_MS.orderTimes();
	// generate effective impedance etc for all layers
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
  }


  // Now form list of indexes of peaks.....
  private int[] formListOfLikelyLayerLocations(double[] S,int n)
  {
	/*
	  Ok, here, the deconvolved seismic is in the first
	  half of S. Search it for ( peaks and return the
	  best candidates in locations[]
	*/
	int i,j,k;
	int[] peak_locations;
	double[] peak_abs_values;

	SortUtils[] peaks = new SortUtils[n];

	k=0;
	for (i=0;i<n;i++) 
	  {
		// go uphill till a peak found
		while ((S[i]<S[i+1]) && (i<n)) {i++;}
		if (i!=0 && i<n)  // no peak allowed if we start downhill
		  {
			peaks[k]=new SortUtils();
			peaks[k].x = Math.abs(S[i]);
			peaks[k].indx = i;
			k++;
		  }
		// go downhill
		while ((S[i] >= S[i+1]) && (i<n)) {i++;}
		if (i<n)
		  {
			peaks[k]=new SortUtils();
			peaks[k].x = Math.abs(S[i]);
			peaks[k].indx = i;
			k++;
		  }
	  }
	// Ok there are now k-1 peaks. Sort in ascending order.
	if (k>2) 
	  {
		Arrays.sort(peaks,1,k);
		// Store the biggest ones (from end of list)
		IntArrayList locations = new IntArrayList();
		
		for (i=k-1,j=0; i>=1 ; i--,j++) 
		  {
			locations.add(peaks[i].indx);
			// System.out.println(peaks[i].indx + " " + peaks[i].x);
		  }
		locations.trimToSize();
		return locations.elements();
	  }
	else
	  { 
		return null;
	  }
  }

  /** 
	  debug method to dump synthetic seismic of realisation.  the
	  combo number is written into the
	  _m._SI.MF.realisationNumberHeaderField ("mark"?) header word
  */
  private void debugDumpSynthetics(int combo_number)
  {
	// might as well inherit inline/xline fields from model
	_m.SM.setInlineXlineHeaderFields(_m._SI.tr_model,tr_near_synth);
	// dump out synthetic
	// set requisite fields
	tr_near_synth.set("ns",(short)_m.SM.ns);
	tr_near_synth.set("f1",(float)_m.SM.t_start);
	tr_near_synth.set("dt",_m.SM.near_wavelet.dt);
	tr_near_synth.set("delrt",(short)Math.rint((_m.SM.t_start*1e3)));
	// and samples  from near synthetic
	for (int j=0;j<_m.SM.ns;j++) 
	  {
		tr_near_synth.set_data_value((float)_m.SM.padded_synth_near[j],j); 
	  }	
	// have to update this field whether accepted or not
	tr_near_synth.set(_m._SI.MF.realisationNumberHeaderField,(short)combo_number);	
	suwriter.fputtr(near_synth_out_strm, tr_near_synth);
  }


 private void initialiseDebugSyntheticSeismicOutput()
  {
	suwriter = new SuWriter();  // set up synthetic seismic output
	tr_near_synth = new SuSegyFormat();		
	// Nears. Open stream.
	try 
	  {
		near_synth_out_strm = (FileOutputStream) new 
		  FileOutputStream(SYNTH_NEAR_FILENAME);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open" 
								   + SYNTH_NEAR_FILENAME);
	  }
  }

  private void finalizeDebugSyntheticSeismicOutput()
  {
	if (near_synth_out_strm!=null)
	  {
		try
		  {
			suwriter.flush(near_synth_out_strm); // must do this!
			near_synth_out_strm.close();
			if (_m._SI.verbosity>2) 
			  {
				_m._SI.EH.println("%nWrote event-detection debug synthetic data to " + SYNTH_NEAR_FILENAME);
			  }
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close " + SYNTH_NEAR_FILENAME);
		  }
	  }
  }
}
