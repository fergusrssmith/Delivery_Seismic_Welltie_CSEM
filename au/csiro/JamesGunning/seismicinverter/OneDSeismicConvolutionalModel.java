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

import java.io.*;
import java.util.Arrays;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.MathUtil;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.Numerical.djbFFT.DJBFFT;
import au.csiro.JamesGunning.sulib.*;

import static au.csiro.JamesGunning.lib.util.MathUtil.*;

/**
   Simple 1d, no-multiples convolutional model of seismic response
   using Zoepritz reflection coefficients to second order
   in incidence angle, first order in property-changes.
*/

//public class  OneDSeismicConvolutionalModel implements SeismicModel
public final class OneDSeismicConvolutionalModel
{
  /**
	 reference to current model, from which the seismic forward model can
	 access the current layer properties via _m.LP[i]
  */
  Model _m;

  /**
	 Some of the calls here are reusable in contexts other than that of the
	 SeismicInverter. to use them, initialise the OneDSeismicConvolutionalModel
	 object with hasInverterContext = false
  */
  private boolean hasInverterContext = false;

  private int n_layers;

  private ErrorHandler _EH;

  StackData nearstack_data = null;
  StackData farstack_data  = null;
  boolean doSNchecks = false;  // assum only one check unless noise level varies

  
  Wavelet near_wavelet=null;
  Wavelet far_wavelet=null;

  String  wavelet_filename;

  /**
	 AVO constants that define deviation of reflectivity
	 from theoretical response A_avo=B_avo=1. These
	 may vary for each trace, or may be set to constants
	 See section "Synthetic seismic" in documentation
  */  
  boolean A_B_avo_synchronised = false;
  protected ModelPriorDescriptions.PriorDetails A_avo;
  protected ModelPriorDescriptions.PriorDetails sigma_A_avo;
  protected ModelPriorDescriptions.PriorDetails B_avo;
  protected ModelPriorDescriptions.PriorDetails sigma_B_avo;

  /**
	 And the LayerProperties variables which will be
	 used in the model sampler to represent that actual AVO
	 coefficient variables
  */
  public LayerProperties.Property A;
  public LayerProperties.Property B;

  private int Nstacks;
  int Ncoarse_tot;        // Total of Ncoarse[i_stack] summed over all stacks - updated on the fly

  /** 
	  Useful constants for trace-based information about
	  <b>current</b> synthetic
  */
  int ns;
  double delta_t;                // time sampling interval in seconds in synthetics
  double delta_t_indep_samples;  // time sampling interval in seconds for seismic-error acquisition (slower than above)
  double delta_t_indep_nearest;  // sampling time nearest to delta_t_indep_samples which is a multiple of delta_t
  double delta_t_ratio_correction_factor;      // =sqrt(delta_t_indep_nearest/delta_t_indep_samples)
  double t_start;                // start of time interval over which to compute synthetics
  double t_end;                  // end of time interval over which to compute synthetics
  boolean useExclusionZone = false; // use exclusion zone
  double t_start_exclusion[];    // utility variables for exclusion zone, if used
  double t_end_exclusion[];      // as above.
  double[] padded_synth_near;
  double[] padded_synth_far;
  double[] seismic_error_vector;  // error at each coarse sample time, concatenated over all stacks
  double hiddenSeismicError;      // true synthetic seismic error computed whether seismic is ignored (-IS) or not
  double logDetCnoiseAnd2PiPieces;// Determinant corresponding to seismic noise likelihood

  double[] R_near;                // suite of reflection coefficients from layer to layer above, in order
  double[] R_far;                 // suite of reflection coefficients from layer to layer above, in order
  double[] t_fc;                  // suite of internal fluid contact reflection times, in order
  double[] R_near_fc;             // suite of internal fluid contact reflection coefficients, in order
  double[] R_far_fc;              // suite of internal fluid contact reflection coefficients, in order
  double[] layer_top_times;       // internal variable for layer top times
  /**
	 Bunch of constants relating to making of synthetics with Ricker wavelets
  */
  /**
	 Default tail Ricker frequency in Hz
  */
  double f_at_band_edge = 40.0;
  /**
	 User suggested wavelet sampling dt (s)
  */
  double suggested_sampling_dt = 0.004;

  /**
	 This constant guarantees samples of a Ricker-2 wavelet generated
	 random signal will pick up 95% of the equivalent degrees-of-freedom (read energy)
	 if sampled at a time rate dt=SAMPLE_TIME_MAGIC_CONSTANT/f_peak
	 where f_peak is the peak spectral frequency of the Ricker-2 spectrum
  */
  private static final double SAMPLE_TIME_MAGIC_CONSTANT = 0.253;

  /**
	 Ricker-2 wavelets have this ratio of the mean frequency (first
	 moment of the power spectrum) to the peak frequency
  */
  public static final double MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER=1.124;

  private boolean debugDumpFiles = false;     // FIXME remove with all references later on
  long debugCounter = 0 ;                     // FIXME remove?
  private boolean incrementDebugCounter=true; // FIXME remove later on


  public OneDSeismicConvolutionalModel(boolean hasInverterContext, Model m)
  {
	this._m = m;
	this.hasInverterContext = hasInverterContext;
	// error handler
	_EH = m._SI.EH;
	if (!m._SI.useJNI) DJBFFT.disableNativeFFT();
	if (_m._SI.verbosity>2)
	  {
		if (DJBFFT.haveNativeFFT()) 
		  {
			_EH.println("Using JNI native FFT");
		  }
		else
		  {
			_EH.println("Using non-optimised Java FFT");
		  }
	  }
  }
 
  /**
	 To use some methods outside the context of the Model and
	 SeismicInverter classes, use this constructor
  */
  public OneDSeismicConvolutionalModel(InversionAnalyser IA)
  {
	this.hasInverterContext=false;
	_EH = IA.EH;
	if (!IA.useJNI) 
	  {
		DJBFFT.disableNativeFFT();	
		_EH.warn("Using non-optimised Java FFT");
	  }
  }

  public void setNearstackData(StackData sd)
  {
	this.nearstack_data = sd;
  }

  public void setFarstackData(StackData sd)
  {
	this.farstack_data = sd;
  }
  
  public StackData  getNearstackData()
  {
	return nearstack_data;
  }

  public StackData  getFarstackData()
  {
	return farstack_data;
  }


  /** 
	  Generic function to set up calculations required for seismic model
  */
  public boolean initialiseSeismicModel(int n_layers)
  {
	if (hasInverterContext)
	  {
		if (_m._SI.makeSyntheticAndWavelet)
		  {
			setSyntheticParameters((double)_m._SI.syntheticFreq, 
								   (double)_m._SI.syntheticSampleTime);
			makeRickerWavelet();
		  }
		else
		  {
			// read in wavelets
			readInWavelets();
			// check SN at initialisation (done at least once)
			if (!isSNratioOK())
			  {
				_m._SI.EH.println("Signal to noise ratio is unrealistically large. "+
						   "Check specified noise level in XML");
				return false;
			  }
		  }
	  }

	R_near = new double[n_layers];
	R_far  = new double[n_layers];

	// contacts Ref's
	t_fc = new double[n_layers];
	//	d_fc  = new double[n_layers];
	R_near_fc = new double[n_layers];  // FIXME check size!!!!!
	R_far_fc  = new double[n_layers];

	layer_top_times = new double[n_layers+1]; // includes t_base
	this.n_layers = n_layers;

	return true;
  }

  /**
	 precompute of log noises for efficiency
  */
  void precomputeLogNoise()
  {
	if (nearstack_data != null && nearstack_data.noise_rms != null)
	  {
		nearstack_data.computeLogNoise();
	  }
	if (farstack_data != null && farstack_data.noise_rms != null)
	  {
		farstack_data.computeLogNoise();
	  }
  }

  /** 
	  Check if current SN ratio seems credible
  */
  boolean isSNratioOK()
  {
	double ratio;
	if (nearstack_data != null && nearstack_data.noise_rms != null)
	  {
		ratio = 10.0*nearstack_data.noise_rms.value/near_wavelet.max_abs;
		if (_m._SI.verbosity>2 && ratio<0.01)
		  {
			_m._SI.EH.println("WARNING: noise-to-signal ratio (very approx. "+
						   1.0/ratio+") is very aggressive! Inverter may fail!");
		  }
		if (ratio < SeismicInverter.MIN_ACCEPTABLE_NS_RATIO)
		  {
			return false;
		  }
		nearstack_data.computeLogNoise();
	  }
	if (farstack_data != null && farstack_data.noise_rms != null)
	  {
		ratio = 10.0*farstack_data.noise_rms.value/far_wavelet.max_abs;
		if (_m._SI.verbosity>2 && ratio<0.01)
		  {
			_m._SI.EH.println("WARNING: noise-to-signal ratio (very approx. "+
						   1.0/ratio+") is very aggressive! Inverter may fail!");
		  }
		if (ratio < SeismicInverter.MIN_ACCEPTABLE_NS_RATIO)
		  {
			return false;
		  }
		nearstack_data.computeLogNoise();
	  }
	return true;
  }

  /**
	 Set time range over which to compute synthetic seismic
  */
  public void setSyntheticAcquisitionRange(double tStart, double tEnd)
  {
	t_start = tStart;
	t_end = tEnd;
	if (t_start >= t_end)
	  throw new RuntimeException("Have set t_start = "+ tStart+
								 " for synthetic seismic to more than t_end = "+ tEnd); 
	// (This should never happen)
  }

  void generateSyntheticTrace()
  {
	int j;
	double t_end_max   = -Double.MAX_VALUE;
	double t_start_min = Double.MAX_VALUE/10;
	Wavelet wavelet=null;
	StackData stack_data=null;

	// set up layer top times, with t_base
	for (j=0; j<n_layers; j++) layer_top_times[j] = _m.LP[j].t.value;
	layer_top_times[n_layers] = _m.LP[n_layers-1].t_base.value;

	// First figure out the time interval required for a meaningful synthetic
  	// Make the interval the range of times spanned by the trace prior mean-times,
	// from top of layer 2 to top of last layer, extended by largest wavelet codas or precursors
	// Further, round this to microseconds so the delrt field can be precise
	for (int stack_cntr=0;stack_cntr<=1;stack_cntr++)
	  {
		if (stack_cntr==0)
		  {
			wavelet = near_wavelet;	
		  }
		else if (stack_cntr==1)
		  {
			wavelet = far_wavelet;
		  }
		// there may be only a near stack
		if (wavelet==null) continue;
		// start at first reflector minus maximal DTcoda over all wavelets
		// double t_try = _m.MPD.MetaInfo[1].t.value - wavelet.getDTcoda();
		// Safety: see
		double t_try = _m.MPD.min_interesting_reflector_time_mean - wavelet.getDTcoda();
		t_start_min = (t_try<t_start_min)? t_try : t_start_min;
		// end at max(DTprecursor) past last reflector
		// t_try = _m.MPD.MetaInfo[_m.n_layers-1].t.value + wavelet.getDTprecursor();
		t_try = _m.MPD.max_interesting_reflector_time_mean + wavelet.getDTprecursor();
		t_end_max = Math.max(t_end_max, t_try);
	  }
    // set range from this, but rounding so the delrt fields can be precise
	if (_m._SI.makeSyntheticAndWavelet || _m._SI.dumpSyntheticSeismic)
	  {
		setSyntheticAcquisitionRange(_m._SI.tSyntheticStart,_m._SI.tSyntheticEnd);
	  }
	else 
	  {
		setSyntheticAcquisitionRange(MathUtil.fCeil(t_start_min* 1.0e3)/1.0e3, t_end_max);		
	  }
	// QC the range if dumpSyntheticSeismic is set: it has to be over a larger range than needed
	// for the seismic error acquistion
	if (_m._SI.dumpSyntheticSeismic)
	  {
		if (_m._SI.tSyntheticStart > MathUtil.fCeil(t_start_min* 1.0e3)/1.0e3)
		  {
			_EH.die_err("Synthetic range ("+_m._SI.tSyntheticStart+","+_m._SI.tSyntheticEnd+
						") supplied starts at too late a time for valid seismic error acquisition in model");
		  }
		if (_m._SI.tSyntheticEnd < t_end_max)
		  {
			_EH.die_err("Synthetic range ("+_m._SI.tSyntheticStart+","+_m._SI.tSyntheticEnd+
						") supplied ends at too early a time for valid seismic error acquisition in model");
		  }
	  }
	// now get the actual synthetic seismic
	generateSyntheticTrace_kernel(null);
  }

  /**
	 The kernel of the synthetic seismic generator. Usable without a sane Model and 
	 SeismicInverter context, but easily wrapped when that context exists 
	 (see  generateSyntheticTrace() ). In this case the Wavelet argument will be null. 

	 Without the inverter context, requires a Wavelet
	 object to be sensibly setup. The internal fields layer_top_times[]
	 and R_near[] have to be sensibly filled before calling this kernel.
	 note the R_near[i] is an array of successive reflection coefficients, which,
	 in the presence of pinchouts, may not necessarily mean the reflection coefficient
	 from layer i-1 to layer i.
  */

  void generateSyntheticTrace_kernel(Wavelet wavelet)
  {
	int i,j,k,stack_cntr,R_cntr,R_fc_cntr;
	double V_denom=0.0, ref_coeff=0, av_theta_sq, wt;
	double vp_eff_1,vp_eff_2,vs_eff_1,vs_eff_2,rho_eff_1,rho_eff_2;
	double t, t_next;
	double frac_drho;
	double frac_dvp;
	double frac_dvs;
	double vs_on_vp;

	double[] R;
	double[] padded_R;
	double[] padded_wavelet;
	// Wavelet wavelet=null;
	StackData stack_data=null;
	boolean isPSdata=false;
	boolean handleContacts=false;
	double th_info;

	// An interface between layers is notated _1 = "above", _2 = "below"
	int i_1,i_2,indx_2=0,indx_1;
	LayerProperties LP_1,LP_2;
	// A fluid contact interface uses these
	int i_fc_1,i_fc_2;
	LayerProperties.InternalFluidContactInfo _FCI;
	LayerProperties.InternalPermeableFluidContactInfo _PFCI;

	// zero these arrays so no old coefficients get into the lists if there are pinchouts
	if (hasInverterContext) 
	  {
		Arrays.fill(R_near,0);
		Arrays.fill(R_far,0);
		Nstacks = 2; 
		if (_m._SI.doSimpleFluidUnitsAndContacts) 
		  {
			handleContacts = true;
			Arrays.fill(t_fc,0);
			Arrays.fill(R_near_fc,0);
			Arrays.fill(R_far_fc,0);
		  }
	  }
	else
	  {
		Nstacks = 1;
	  }
	// loop over near, far stacks
	for (stack_cntr=0; stack_cntr < Nstacks; stack_cntr++)
	  {
		// in inverter context, there is a different wavelet for each stack
		if (hasInverterContext) 
		  {
			switch (stack_cntr)
			  {
			  case 0:
				wavelet = near_wavelet;
				stack_data = nearstack_data;
				break;
			  case 1:
				wavelet = far_wavelet;
				stack_data = farstack_data;
				break;
			  }
			// there may be only a near stack
			if (stack_data==null || wavelet==null) continue;
		  }

		// sampling time same as wavelet
		delta_t = ((double)wavelet.dt/1.0e6);
		// number of samples in total synthetic interval
		// ns = 1 + (int)((t_end - t_start)/delta_t);
		ns = 1 + MathUtil.fCeil_int((t_end - t_start)/delta_t);

		// reflection coeff array
		R = new double[ns];
		// debugCounter++;  // TEST REMOVE!
		// useful factor for angles. this is what distinguishes stacks
  		if (hasInverterContext) 
		  {
			// Dix equation
			V_denom = pow4(stack_data.getVstack()) * sqr(stack_data.getTime())/stack_data.getOffsetMeanSquare();
			// stackAnglePrescribed = stack_data.isStackAnglePrescribed();
			if (_m._SI.doSimpleFluidUnitsAndContacts) handleContacts = true;
			isPSdata = stack_data.isPSdata;
			if (!isPSdata)
			  {
				th_info = stack_data.sq_theta_specified;
			  }
			else
			  {
				th_info = stack_data.theta_specified;
			  }			
		  }
		else
		  {
			th_info = -Double.MAX_VALUE; // correct for Analyser context: uses stored, not computed refl coeff values
		  }
		if (false)
		  {
			// Tried-and tested legacy version!

			// run down time sequence to pick up layer changes & compute refl coeffs
			i_1=0;
			i_2=1;
			// occasionally layer times end up above where we start (t_start). Trap this.
			for (j=0; j<n_layers; j++)
			  {
				if (t_start >= layer_top_times[j])
				  {
					i_1=j;
					i_2=j+1;
				  }
				else
				  {
					break;
				  }
			  }
			if (hasInverterContext) LP_1 =  _m.LP[i_1];

			// first *layer* reflection can not be earlier than indx=1. Start here
			R_cntr=1;
			// Where contact occurs can be anywhere. Starting at zero maximises the number we can stuff into n layers
			R_fc_cntr=0;     

			// loop over sample times
			for (j=0, t=t_start ; t<=t_end && i_2<n_layers && j<ns-1 ; t+=delta_t, j++)
			  {
				if ((t >= layer_top_times[i_2]) && (i_2 < n_layers))
				  {
					// we just jumped into a new layer.....				
					// check for pinch-outs. Loop until we see a genuine new layer
					do
					  {
						i_2++;
					  }
					while (i_2<n_layers && t >= layer_top_times[i_2]);
					i_2--;
					if (hasInverterContext) LP_2 =  _m.LP[i_2];

					// OK, at this point, i_2 & i_1 delimit a range of possible layer-top times that
					// fit in the time interval [t-delta_t. t].
					// Nearly always there is only one layer top time in this interval, but there may be several.
					// When there are several, they may be (a) ordered - a whole bunch of very thin layers, so
					// the reflection coefficient must be summed over all such transitions, or 
					// (b) they may describe a pinchout, and there will be contrasts between 
					// non-adjacent layer numbers.
					// The bookkeeping for this has to be careful
					i=i_1;
					k=i+1;
					while (k<=i_2)
					  {
						if (hasInverterContext)
						  {
							t_next = (k+1 > n_layers-1) ? _m.LP[k].t_base.value : layer_top_times[k+1];
						  }
						else
						  {
							t_next = (k+1 > n_layers-1) ? t_end : layer_top_times[k+1];
						  }
						if (t_next > layer_top_times[k])
						  {
							// no pinch out. Layer boundaries follow in time						
							indx_1 = i;
							indx_2 = k;
							if (hasInverterContext)
							  {
								// Compute the actual reflection coefficient
								LP_1 =  _m.LP[indx_1];
								LP_2 =  _m.LP[indx_2];
								ref_coeff = getReflCoeff(isPSdata,th_info,V_denom,LP_1.rho_eff,LP_1.vp_eff,LP_1.vs_eff,LP_2.rho_eff,LP_2.vp_eff,LP_2.vs_eff);

								// in case we fall outside the model, put no reflections (should never happen in a decent model)
								if (t<layer_top_times[0] || t>_m.LP[n_layers-1].t_base.value) ref_coeff = 0.0;
								// store it for outputting when we need to
								if (stack_cntr==0) R_near[R_cntr] = ref_coeff;
								if (stack_cntr==1) R_far[R_cntr] = ref_coeff;
								R_cntr++;
							  }
							else
							  {
								// Use the one already set up for us outside this call
								ref_coeff = R_near[R_cntr++];
							  }
							// distribute this across the 4 nearest times suitably
							// (Lagrange 4pt interpolation)
							wt = (t - layer_top_times[indx_2])/delta_t;

							if (wt<0 || wt>1.0) 
							  {
								_EH.warn("t="+t+", layer_top_times[indx_2]="+layer_top_times[indx_2]+", wt = "+wt+
										 " i_1="+i_1+", i_2="+i_2+", cntr = "+debugCounter);
								if (hasInverterContext) dumpDetail();
								throw new RuntimeException("Interpolation wrong! fixme!"); // (This should never happen)
							  }
							if (j>1) 
							  {
								LagrangeExtirpolate(j, R, wt, ref_coeff);
							  }
							else
							  {
								R[j] += ref_coeff;
							  }
							i=k;
							k=i+1;
						  }
						else
						  {
							// a pinch out. Look for the next real layer
							k++;
						  }
					  }

					// and start looking for the next transition
					i_1=i_2;
					i_2++;
				  }
			  }
		  }

		if (true)
		  {
			// New version for fluid contacts
			boolean haveReflection;
			R_cntr=1;
			R_fc_cntr=0;

			if (hasInverterContext)
			  {
				// for (indx_1=0; indx_1<n_layers-1; indx_1++)
				int indx_1_max = ( handleContacts ? n_layers : (n_layers-1));
				for (indx_1=0; indx_1 < indx_1_max; indx_1++)
				  {
					if (layer_top_times[indx_1+1] > layer_top_times[indx_1])
					  {
						// OK, layer indx_1 has non-zero thickness
						// save properties
						LP_1 =  _m.LP[indx_1];
						rho_eff_1 = LP_1.rho_eff; vp_eff_1 = LP_1.vp_eff; vs_eff_1 = LP_1.vs_eff;

						// Do contacts internal reflections (imples inverter context)
						if (handleContacts)
						  {
							_FCI = _m.LP[indx_1].fcInfo;
							if (!_FCI.is_permeable)
							  {
								// Pure shale layer; keep track of hidden t_fc times
								for (i_fc_1=0; i_fc_1<_FCI.n_internal_sublayers-1; i_fc_1++)
								  {
									if (_FCI.t_contact[i_fc_1+1] > _FCI.t_contact[i_fc_1])
									  {
										// OK, sublayer i_fc_1 has non-zero thickness
										// Seek internal reflection
										haveReflection = false;
										for (i_fc_2=i_fc_1+1; i_fc_2<_FCI.n_internal_sublayers; i_fc_2++)
										  {
											if (_FCI.t_contact[i_fc_2+1] > _FCI.t_contact[i_fc_2])
											  {
												// layer i_fc_2 is first one after i_fc_1 with non-zero thickness. Can support a reflection
												haveReflection = true;
												break;
											  }
										  }
										if (haveReflection)
										  {
											if (stack_cntr==0) t_fc[R_fc_cntr] = _FCI.t_contact[i_fc_2];
											R_fc_cntr++;
										  }
									  }
								  }
							  } 
							else if (_FCI.n_internal_sublayers>1)
							  {
								_PFCI = (LayerProperties.InternalPermeableFluidContactInfo) _FCI;
								// permeable layers (_PFCI != null && _PFCI.n_internal_sublayers>1)
								for (i_fc_1=0; i_fc_1<_PFCI.n_internal_sublayers-1; i_fc_1++)
								  {
									if (_PFCI.t_contact[i_fc_1+1] > _PFCI.t_contact[i_fc_1])
									  {
										// OK, sublayer i_fc_1 has non-zero thickness
										// Seek internal reflection
										haveReflection = false;
										for (i_fc_2=i_fc_1+1; i_fc_2<_PFCI.n_internal_sublayers; i_fc_2++)
										  {
											if (_PFCI.t_contact[i_fc_2+1] > _PFCI.t_contact[i_fc_2])
											  {
												// layer i_fc_2 is first one after i_fc_1 with non-zero thickness. Can support a reflection
												haveReflection = true;
												break;
											  }
										  }
										if (haveReflection)
										  {
											// Compute the actual internal reflection coefficient
											ref_coeff = getReflCoeff(isPSdata,th_info,V_denom,_PFCI.rho_eff[i_fc_1],_PFCI.vp_eff[i_fc_1],
																	 _PFCI.vs_eff[i_fc_1],_PFCI.rho_eff[i_fc_2],
																	 _PFCI.vp_eff[i_fc_2],_PFCI.vs_eff[i_fc_2]);
											// Save sublayer 2 properties since they're going to be needed as "1" when we jump into next layer
											rho_eff_1 = _PFCI.rho_eff[i_fc_2]; vp_eff_1 = _PFCI.vp_eff[i_fc_2]; vs_eff_1 = _PFCI.vs_eff[i_fc_2];

											// store ref_coeff for outputting when we need to
											if (stack_cntr==0) 
											  {
												t_fc[R_fc_cntr] = _PFCI.t_contact[i_fc_2];
												R_near_fc[R_fc_cntr] = ref_coeff;
											  }
											if (stack_cntr==1) R_far_fc[R_fc_cntr] = ref_coeff;
											R_fc_cntr++;

											// Index of time sample in R array just below reflection
											j = 1 + MathUtil.fFloor_int((_PFCI.t_contact[i_fc_2] - t_start)/delta_t);
											t = t_start + j*delta_t;
											// ensure it's within block (takes care of rounding errors)
											t = Math.min(Math.max(t,_PFCI.t_contact[i_fc_2]), _PFCI.t_contact[i_fc_2]+delta_t);
											if (j<(ns-1) && t>=t_start && t<=t_end)
											  {
												// distribute this across the four nearest times (lagrange 4pt interpolation)
												wt = (t - _PFCI.t_contact[i_fc_2])/delta_t;
												wt = Math.max(Math.min(1.0,wt),0);
												if (false)
												  {
													// guaranteed not to happen
													if (wt<0 || wt>1.0) 
													  {
														_EH.warn("Contact reflections 2: t="+t+", _PFCI.t_contact[i_fc_2]="+_PFCI.t_contact[i_fc_2]+", wt = "+wt);
														dumpDetail();
														throw new RuntimeException("Interpolation wrong! fixme!"); // (This should never happen)
													  }
												  }
												if (j>1) LagrangeExtirpolate(j, R, wt, ref_coeff); 
												else R[j] += ref_coeff;
											  }
										  }
										i_fc_1 = i_fc_2-1; // ensures next full layer examined is i_fc_1 since we hit i_fc_1++ on loop.
									  }
								  }
								// Save last layer properties before exiting, as interface to next real layer
							  }
						  }
					  
						// Seek reflection coeff at jump into next layer
						haveReflection = false;
						for (indx_2=indx_1+1; indx_2<n_layers; indx_2++)
						  {
							// base time of next layer
							t_next = (indx_2 == n_layers-1) ? _m.LP[indx_2].t_base.value : layer_top_times[indx_2+1];
							if (t_next > layer_top_times[indx_2])
							  {
								// layer indx_2 is first one after indx_1 with non-zero thickness. Can support a reflection
								haveReflection=true;
								break;
							  }
						  }
						if (haveReflection)
						  {
							// Compute the actual reflection coefficient
							LP_2 =  _m.LP[indx_2];
							rho_eff_2 = LP_2.rho_eff; vp_eff_2 = LP_2.vp_eff; vs_eff_2 = LP_2.vs_eff;

							// Under contacts model, check if layer 2 has internal sublayers, so we get the right _2 properties
							if (handleContacts)
							  {
								_FCI = LP_2.fcInfo;
								if (_FCI.is_permeable && (_FCI.n_internal_sublayers > 1)) // skip shales
								  {
									_PFCI = (LayerProperties.InternalPermeableFluidContactInfo)_FCI;
									for (i_fc_1=0; i_fc_1<_PFCI.n_internal_sublayers; i_fc_1++)
									  {
										if (_PFCI.t_contact[i_fc_1+1] > _PFCI.t_contact[i_fc_1])
										  {
											// OK, sublayer i_fc_1 has non-zero thickness
											break;
										  }
									  }
									// overwrite lower properties with those of uppermost nonzero thickness sublayer just found
									rho_eff_2 = _PFCI.rho_eff[i_fc_1]; vp_eff_2 = _PFCI.vp_eff[i_fc_1]; vs_eff_2 = _PFCI.vs_eff[i_fc_1];
								  }
							  }

							ref_coeff = getReflCoeff(isPSdata,th_info,V_denom,rho_eff_1,vp_eff_1,vs_eff_1,rho_eff_2,vp_eff_2,vs_eff_2);

							// store it for outputting when we need to
							if (stack_cntr==0) R_near[R_cntr] = ref_coeff;
							if (stack_cntr==1) R_far[R_cntr] = ref_coeff;
							R_cntr++;



							// HERE do cross layer reflection between indx_1 & indx_2
							// Index of time sample in R array just below reflection
							j = 1 + MathUtil.fFloor_int((layer_top_times[indx_2] - t_start)/delta_t);
							t = t_start + j*delta_t;
							// ensure it's within block (takes care of rounding errors)
							t = Math.min(Math.max(t,layer_top_times[indx_2]), layer_top_times[indx_2]+delta_t);
							if (j<(ns-1) && t>=t_start && t<=t_end)
							  {
								// distribute this across the four nearest times (lagrange 4pt interpolation)
								wt = (t - layer_top_times[indx_2])/delta_t;
								wt = Math.max(Math.min(1.0,wt),0);
								if (false)
								  {
									// guaranteed not to happenif (wt<0 || wt>1.0) 
									{
									  _EH.warn("Contact reflections 2: t="+t+", layer_top_times[indx_2]="+layer_top_times[indx_2]+", wt = "+wt);
									  dumpDetail();
									  throw new RuntimeException("Interpolation wrong! fixme!"); // (This should never happen)
									}
								  }
								if (j>1) 
								  {
									LagrangeExtirpolate(j, R, wt, ref_coeff);
								  }
								else
								  {
									R[j] += ref_coeff;
								  }
							  }

						  }

						indx_1 = (indx_2-1); // ensures next full layer examined is indx_1 since we hit indx_1++ on loop.
					  }
				  }
			  }   // hasInverter context
			else
			  {
				// Call from InversionAnalyser!! No LayerProperties available etc
				// stuff the known fields in at the apposite times...
				// 2-pass fill: pass = (1) layer reflections, (0) fluid contact reflections
				int pass;
				double t_ref=0;
				boolean foundFirstNonZeroThicknessLayer = false;

				for (pass=1; pass>=0; pass--)
				  {
					// layers only need to go to indx_1=n_layers-2, contacts to indx_1=n_layers-1
					for (indx_1=0; indx_1<n_layers-pass; indx_1++)
					  {
						if (pass==1)
						  {
							// ***layer-jump*** reflections from Analyser
							// layer-jump test as above: thick layer - to - thick layer
							haveReflection = false;
							if (layer_top_times[indx_1+1] > layer_top_times[indx_1])
							  {
								if (!foundFirstNonZeroThicknessLayer) foundFirstNonZeroThicknessLayer = true;
								// OK, layer indx_1 has non-zero thickness
								// Seek reflection coeff at jump into next non-zero thickness layer
								for (indx_2=indx_1+1; indx_2<n_layers; indx_2++)
								  {
									// base time of next layer : note use of t_end!!! LP not available
									t_next = (indx_2 == n_layers-1) ? t_end : layer_top_times[indx_2+1];
									if (t_next > layer_top_times[indx_2])
									  {
										// layer indx_2 is first one after indx_1 with non-zero thickness. Can support a reflection
										haveReflection=true;
										t_ref = layer_top_times[indx_2];
										// use this ref_coeff
										ref_coeff = R_near[R_cntr++];
										break;
									  }
								  }
							  }
							else
							  {
								// Haven't yet found first non-zero thickness layer
								// special care needed..move ahead one
								if (!foundFirstNonZeroThicknessLayer) indx_2 = indx_1 + 1;
							  }
						  }
						else
						  {
							// ***fluid contact*** reflections from Analyser
							t_ref = t_fc[R_fc_cntr];
							// use this ref_coeff
							ref_coeff = R_near_fc[R_fc_cntr];
							haveReflection = (ref_coeff!=0.0);
							R_fc_cntr++;
						  }

						if (haveReflection)
						  {
							// Stuff in reflection...
							// Index of time sample in R array just below reflection
							j = 1 + MathUtil.fFloor_int((t_ref - t_start)/delta_t);
							t = t_start + j*delta_t;
							// ensure it's within block (takes care of rounding errors)
							t = Math.min(Math.max(t,t_ref), t_ref+delta_t);
							if (j<(ns-1) && t>=t_start && t<=t_end)
							  {
								// distribute this across the four nearest times (lagrange 4pt interpolation)
								wt = (t - t_ref)/delta_t;
								wt = Math.max(Math.min(1.0,wt),0);
								if (j>1) 
								  {
									LagrangeExtirpolate(j, R, wt, ref_coeff);
								  }
								else
								  {
									R[j] += ref_coeff;
								  }
							  }
						  }
						if (pass==1) indx_1 = (indx_2-1); // ensures next full layer examined is indx_1 since we hit indx_1++ on loop.
					  }
				  }
			  }
		  }




		// now convolve this reflectivity sequence
		int n;
		// In general R.length not a power of 2, so find the smallest n = 2^q bigger than R.length
		n = FFT.getNearestPowerOfTwo(R.length);

		// now pad out R
		padded_R = new double[n];
		padded_wavelet = new double[n];
		switch (stack_cntr) 
		  {
		  case 0:
			padded_synth_near = new double[2*n];
			break;
		  case 1:
			padded_synth_far = new double[2*n];
			break;
		  }
		System.arraycopy(R,0,padded_R,0,R.length);
		int wl = wavelet.getLength();
		int wfi = wavelet.getFinalIndex();
		// +ve time part of wavelet
		System.arraycopy(wavelet.W,0,padded_wavelet,0,wfi+1);
		// -ve time part
		System.arraycopy(wavelet.W,wfi+1,
						 padded_wavelet,wfi + 1 + n - wl,
						 wl - 1 - wfi);

		if (debugDumpFiles && stack_cntr==0) 
		  {
			Util.printfileZO("Wpad_"+debugCounter+".txt",padded_wavelet,padded_wavelet.length);
			Util.printfileZO("Rpad"+debugCounter+".txt",padded_R,padded_R.length);
		  }

		switch (stack_cntr) 
		  {
		  case 0:
			//FFT.convlvZO_clean(padded_R,n,padded_wavelet,1,padded_synth_near);
			DJBFFT.convolve(padded_R,n,padded_wavelet,padded_synth_near);
			break;
		  case 1:
			//FFT.convlvZO_clean(padded_R,n,padded_wavelet,1,padded_synth_far);
			DJBFFT.convolve(padded_R,n,padded_wavelet,padded_synth_far);
			break;
		  }

		if (debugDumpFiles && stack_cntr==0)  
		  {
			Util.printfileZO("Synth"+debugCounter+".txt",padded_synth_near,padded_synth_near.length);
		  }
	  }
	// clean up to help GC
	R=null;
	padded_R=null;
	padded_wavelet=null;
  } 


  double getReflCoeff(boolean isPSdata,double th_info,double V_denom, double rho_1,double vp_1,double vs_1,double rho_2,double vp_2,double vs_2)
  { 
	double frac_drho,frac_dvp,frac_dvs,beta,vpBar,vsBar;
	double R,th;

	// Dix equation
	// av_theta_sq = vp_1^2/V_denom;
	vpBar = 0.5*(vp_2 + vp_1);
	vsBar = 0.5*(vs_2 + vs_1);
	frac_drho = (rho_2 - rho_1)/((rho_2 + rho_1)/2.0);
	frac_dvp = (vp_2 - vp_1)/vpBar; 
	frac_dvs = (vs_2 - vs_1)/vsBar;
	beta = vsBar/vpBar;
							
	if (!isPSdata)
	  {
		// PP reflection
		if (th_info>=0) 
		  {  // usually prescribed-angle mode
			th = B.value * th_info;
		  }
		else
		  {  // ray-trace mode
			th = B.value * (vp_1 * vp_1 /V_denom);
		  }
		
		R = (A.value * (frac_drho + frac_dvp)/2.0 + th * (frac_dvp/2.0 - 2.0*(beta * beta)*(frac_drho + 2.0*frac_dvs)));
	  }
	else
	  {
		// PS reflectivity (See Appendix 7.1 in Wavelet Extraction writeup)
		// Omit A, B fudge factors for PS. th_info is just theta here
		if (th_info>=0) 
		  { // usually prescribed-angle mode
			th = th_info;
		  }
		else
		  { 
			// ray-trace mode
			th = Math.sqrt(vp_1 * vp_1 /V_denom);	
		  }
		// usually prescribed-angle mode
		R = th * (-0.5*(frac_drho + 2.0 * beta * (frac_drho + 2.0*frac_dvs)) + 
				  sqr(th) * ((1.0 + beta*(8 + 9.0*beta))*frac_drho  +  beta*(16.0 + 24.0*beta)*frac_dvs)/12.0);
	  }
	return R;
  }

  /**
	 extirpolate the reflection coefficient ref_coeff on to the 
	 array R, given the ref_coeff occurs between j-1 and j at j-wt
  */
  void LagrangeExtirpolate(int j, double[] R, double wt, double ref_coeff)
  {
	// Use four-point Lagrange formula. 
	// Compute the four coefficients A_{-1},A_0,A_1,A_2
	// as per p 879 of Abramowitz & Stegun....
	double A_neg1,A_0,A_1,A_2,tmp1,tmp2;
	
	A_neg1=-wt*(wt-1.0)*((tmp1=wt-2.0))/6.0;
	A_0=((tmp2=wt*wt-1.0))*tmp1/2.0;
	A_1=-wt*(wt+1.0)*tmp1/2.0;
	A_2=wt*tmp2/6.0;

	// Lagrange 4-pt  interpolation
	R[j-2] += A_2 * ref_coeff;
	R[j-1] += A_1 * ref_coeff;
	R[j] += A_0 * ref_coeff;
	R[j+1] += A_neg1 * ref_coeff;
  }

  private void dumpDetail()
  {
	for (int i=0;i<n_layers;i++)
	  {
		System.out.println("Layer "+i+", t = "+_m.LP[i].t.value);
	  }
	System.out.println("Layer "+(n_layers)+", tbase = "+_m.LP[n_layers-1].t_base.value);
  }
  
  /**
	 Compute synthetics and acquire seismic mismatch error. Returns
	 chisq/2 for the mismatch (i.e. the negative log likelihood)
  */
  public double computeSyntheticAndSeismicError()
  {
	int Nstacks=2;             // So far, two stacks, (nears and fars)
	int stack_cntr,i,j;
	int j_start,j_end;
	int ns_true;               // num of samples in current trace from seismic stack
	int[] Ncoarse;             // number of coarse samples in each stack used to accumulate error
	// int Ncoarse_tot;        // Total "" summed over all stacks
	double i_cont;
	double t;                  // time in comparison loop
	double[] t_start_comparison; // time at which valid comparison can be started
	double[] t_end_comparison;   // time at which to end comparisons
	double t_true_start;       // start of recording time of true seismic trace
	double delta_t_true;       // sampling time of true seismic
	double wt;                 // interpolation weight for true seismic
	double amp;                // true seismic amplitude
	Wavelet wavelet=null;
	SuSegyFormat tr=null;
	StackData stack_data=null;
	double[] S=null;           // reference to synthetic seismic for appropriate stack
	double noise;              // will be rms noise level for the appropriate stack
	double chisq=0.0;
	
	hiddenSeismicError=0.0;  

	// first, generate synthetics
	generateSyntheticTrace();

	// number of fine samples in total synthetic interval
	// ns = 1 + (int)((t_end - t_start)/delta_t);
	ns = 1 + MathUtil.fCeil_int((t_end - t_start)/delta_t);
	
	Ncoarse=new int[Nstacks];
	t_start_comparison = new double[Nstacks];
	t_end_comparison = new double[Nstacks];
	// Count total number of error-mismatch data points, loop over near, far stacks
	Ncoarse_tot=0;
	logDetCnoiseAnd2PiPieces=0;
	for (stack_cntr=0;stack_cntr<=1;stack_cntr++)
	  {
		if (stack_cntr==0)
		  {
			wavelet = near_wavelet;
			stack_data = nearstack_data;
		  }
		else if (stack_cntr==1)
		  {
			wavelet = far_wavelet;
			stack_data = farstack_data;
		  }
		// there may be only a near stack
		if (wavelet==null) continue;
		// Allowable range over which synthetic seismic is to be *used*
		// NEW: override with user choice (default is DTcoda)
		t_start_comparison[stack_cntr] = _m.MPD.MetaInfo[1].t.value - stack_data.getTabove();
		
		// NEW: override with user choice (default is DTprecursor)
		t_end_comparison[stack_cntr] = _m.MPD.MetaInfo[_m.n_layers-1].t.value + stack_data.getTbelow();

		//_m._SI.EH.print("time comparison range "+ t_start_comparison[stack_cntr] + " - " +	t_end_comparison[stack_cntr]);
		// Number of coarse samples we will use
		Ncoarse[stack_cntr] = (int)Math.rint((t_end_comparison[stack_cntr] - t_start_comparison[stack_cntr])
											 /delta_t_indep_nearest) + 1;
		Ncoarse_tot += Ncoarse[stack_cntr];
		if (Ncoarse_tot<=0) 
		  {
			throw new RuntimeException("Synthetic seismic region is too narrow to get meaningful seismic. Must broaden the model!");
		  }

		// Note: uncorrelated noise model!
		// noise = stack_data.getNoise();
		// last term in next eqn should always be Math.log(stack_data.noise_rms.value);
		logDetCnoiseAnd2PiPieces += Ncoarse[stack_cntr] * stack_data.log_noise_rms; 
		logDetCnoiseAnd2PiPieces += Ncoarse[stack_cntr] * MathUtil.HALF_LOG_TWOPI;

		// switch off seismic likelihood if desired
		if (_m._SI.ignoreSeismic)
		  {
			// I think we need this...
			logDetCnoiseAnd2PiPieces = 0;
		  }
		else if (false) // I mean here something like the below: if (useExclusionZone && isTimeinAnExclusionZone(t))
		  {
			// FIXME! - need to put in adjustment for hole, but efficiently....
		  }
	  }




	// now loop over near, far stacks, accumulating error over valid interval
	// Fill also the array of seismic error values
	seismic_error_vector= new double[Ncoarse_tot];


	int k=0;
	for (stack_cntr=0;stack_cntr<=1;stack_cntr++)
	  {
		if (stack_cntr==0)
		  {
			wavelet = near_wavelet;
			tr = _m._SI.tr_near;
			S = padded_synth_near;
			stack_data = nearstack_data;
		  }
		else if (stack_cntr==1)
		  {
			wavelet = far_wavelet;
			tr = _m._SI.tr_far;
			S = padded_synth_far;
			stack_data = farstack_data;
		  }

		// there may be only a near stack
		if (wavelet==null) continue;

		// Data. Traces may be delayed if in eg deep water, and delrt!=0 (delrt in ms)
		t_true_start = ((double)tr.get_short("delrt"))*1e-3;
		delta_t_true = ((double)tr.get_short("dt"))*1e-6;
		ns_true = (int)tr.get_short("ns");

		// noise level
		noise = stack_data.getNoise();

		// Range of indices of padded_synth that are valid (inclusive)
		// OLD		j_start = (int)Math.rint(wavelet.getDTprecursor()/delta_t);
		j_start = (int)Math.rint((t_start_comparison[stack_cntr]-t_start)/delta_t);
		// this will be exact by costruction of delta_t_indep_nearest
		int delta_j =  (int)Math.rint(delta_t_indep_nearest/delta_t);
		j_end = j_start + (Ncoarse[stack_cntr]-1) * delta_j;

		// SAFETY HACK FIXME SHOULD NEVER BE NEEDED!!! FIX SO NEXT LINE IS NOT NEEDED!
		j_start=Math.max(0,j_start); 


		// Run over samples that can legitimately be compared
		// Protect against i being too big or small
		if (t_start_comparison[stack_cntr] - t_true_start < 0)
		  throw new RuntimeException("True seismic traces do not start at early enough times"+
									 " in computeSyntheticAndSeismicError().");
		if (MathUtil.fFloor((t_end_comparison[stack_cntr] - t_true_start)/delta_t_true) + 2 >= ns_true)
		  {
			throw new RuntimeException("True seismic traces do not contain late enough times"+
									   " in computeSyntheticAndSeismicError().");
		  }


		// OLD 
		// for (j=j_start, t=t_start_comparison[stack_cntr] ; j<=j_end ; t+=delta_t_indep_nearest, j+=delta_j)
		// NEW
		double t_start_compare  = j_start*delta_t + t_start;
		for (j=j_start, t=t_start_compare; j<=j_end ; t+=delta_t_indep_nearest, j+=delta_j)
		  {
			i_cont = (t - t_true_start)/delta_t_true;
			// i = (int)MathUtil.fFloor(i_cont);
			// Are guaranteed i_cont>0 from tests above, so
			i = (int)i_cont;
			wt = i_cont - i; 
			if (true)
			  {
				// linearly interpolate true seismic
				// amp = tr.get_data_value(i)*(1.0-wt) + tr.get_data_value(i+1)*wt;
				// fast access using prestored data
				amp = stack_data.currentTrace[i]*(1.0-wt) + stack_data.currentTrace[i+1]*wt;
			  }
			else
			  {
				// Lagrange 4pt interpolation: p 879 of Abramowitz....
				double A_neg1,A_0,A_1,A_2,tmp1,tmp2; 
				A_neg1=-wt*(wt-1.0)*((tmp1=wt-2.0))/6.0;
				A_0=((tmp2=wt*wt-1.0))*tmp1/2.0;
				A_1=-wt*(wt+1.0)*tmp1/2.0;
				A_2=wt*tmp2/6.0;
				amp =  tr.get_data_value(i-1)*A_neg1 + tr.get_data_value(i)*A_0 
				  + tr.get_data_value(i+1)*A_1 +  tr.get_data_value(i+2)*A_2;
			  }

			double synth_amp = S[j];                         // synth seismic at this sample time
			seismic_error_vector[k]=(amp-synth_amp)/noise;   // remember error vector

			// Stru_vector[k] = amp; _m._SI.EH.println("REMOVE!!");
			// Ssyn_vector[k] = synth_amp; _m._SI.EH.println("REMOVE!!");


			// to get samples on the synthetic grid, we probably oversampled just a little
			// undo this effect (delta_t_indep_nearest is usually a little smaller 
			// than delta_t_indep_samples):
			seismic_error_vector[k]  *= delta_t_ratio_correction_factor;
			// switch off seismic likelihood if desired
			if (_m._SI.ignoreSeismic)
			  {
				hiddenSeismicError += sqr(seismic_error_vector[k]);
				seismic_error_vector[k]=0;
			  }
			else if (useExclusionZone && isTimeinAnExclusionZone(t))
			  {
				hiddenSeismicError += sqr(seismic_error_vector[k]);
				seismic_error_vector[k]=0;
				// Remove seismic in synthetic
				S[j] = 0;
			  }
			chisq += sqr(seismic_error_vector[k]);


			k++;
		  }

		// stuff to clean synthetic seismics if desired
		// (i.e. set synthetic to zero outside error-acquisition range)
		if (_m._SI.truncateSynthetics)
		  {
			// precursor
			Arrays.fill(S,0,j_start,0.0);
			// coda
			Arrays.fill(S,j_end+1,S.length,0.0);
		  }

		// Util.printfileZO("Strue.txt",Stru_vector,Stru_vector.length); _m._SI.EH.println("REMOVE!!");
		// Util.printfileZO("Ssyn.txt",Ssyn_vector,Ssyn_vector.length); _m._SI.EH.println("REMOVE!!");
		// System.exit(0); _m._SI.EH.println("REMOVE!!");

		if (debugDumpFiles && stack_cntr==0)
		  {
			//_m._SI.EH.println("j-range is " + j_start + ","+j_end);
			// Util.printfileZO("DAmp.txt",amp_debug,amp_debug.length);
			//_m._SI.EH.println("Dumped Damp.txt, chisq = "+chisq+ ", tst,te="+ 
			// t_start_comparison[stack_cntr]+" " +t_end_comparison[stack_cntr] );
			Util.printfileZO("error"+debugCounter+".txt",seismic_error_vector,seismic_error_vector.length);
			if (incrementDebugCounter) debugCounter++;
		  }
	  }

	chisq/=2.0;

	if (_m._SI.ignoreSeismic) 
	  {
		hiddenSeismicError/=2.0;
	  }
	else
	  {
		hiddenSeismicError = chisq;
	  }


	// clean up to help GC
	Ncoarse = null;
	t_start_comparison = null;
	t_end_comparison = null;
	S=null; 

	if (_m._SI.verbosity>5) _m._SI.EH.println("Seismic negLog= %10.4g ",chisq);

 	return chisq; // make it a true negative log(likelihood)
  }

  /**
	 Returns error vector of model mismatches at all relevant data points,
	 scaled by their respective error standard deviations
  */
  double[] getErrorVector()
  {
	return seismic_error_vector;
  }



  /** 
	  method to dump synthetic seismic of realisation.
	  If newState==true, the current synthetic (from the proposal) will be dumped, otherwise the
	  synthetic of the last accepted proposal (the current chain state) will be dumped. The
	  realisation number is written into the header field realisationNumberHeaderField
  */
  public void dumpSynthetics(boolean newState, int realisation_number, boolean write_state)
  {
	if (_m._SI.dumpSyntheticSeismic || _m._SI.makeSyntheticAndWavelet || _m._SI.dumpOptimalSeismic) 
	  {	
		// nears
		if (near_wavelet!=null && nearstack_data!=null)
		  {
			if (newState)
			  {
				if (_m._SI.makeSyntheticAndWavelet)
				  {
					// entirely new synthetics may do well to inherit 
					// all the fields from the model prior traces
				    _m._SI.tr_near_synth.copyHeader(_m._SI.tr_model);
				  }
				else if (_m._SI.dumpSyntheticSeismic)
				  {
					// might as well inherit inline/xline fields from model
					setInlineXlineHeaderFields(_m._SI.tr_model,_m._SI.tr_near_synth);
				  }
				// dump out synthetic
				// set requisite fields
				_m._SI.tr_near_synth.set("ns",(short)ns);
				_m._SI.tr_near_synth.set("f1",(float)t_start);
				_m._SI.tr_near_synth.set("dt",near_wavelet.dt);
				_m._SI.tr_near_synth.set("delrt",(short)Math.rint((t_start*1e3)));
				// and samples 
				for (int j=0;j<ns;j++) 
				  {
					_m._SI.tr_near_synth.set_data_value((float)padded_synth_near[j],j); 
				  }	
			  }
			// have to update this field whether accepted or not
			_m._SI.tr_near_synth.set(_m._SI.MF.realisationNumberHeaderField,(short)realisation_number);	
			// only write state if insisted on
			if (write_state) 
			  {
				_m._SI.suwriter.fputtr(_m._SI.near_synth_out_strm,_m._SI.tr_near_synth);
			  }
		  }
		// fars
		if (far_wavelet!=null && farstack_data!=null)
		  {
			if (newState)
			  {
				if (_m._SI.makeSyntheticAndWavelet)
				  {
					// entirely new synthetics may do well to inherit 
					// all the fields from the model prior traces
				    _m._SI.tr_far_synth.copyHeader(_m._SI.tr_model);
				  }
				else if (_m._SI.dumpSyntheticSeismic)
				  {
					// might as well inherit inline/xline fields from model
					setInlineXlineHeaderFields(_m._SI.tr_model,_m._SI.tr_far_synth);
				  }
				_m._SI.tr_far_synth.set("ns",(short)ns);	
				_m._SI.tr_far_synth.set("f1",(float)t_start);
				_m._SI.tr_far_synth.set("dt",far_wavelet.dt);
				_m._SI.tr_far_synth.set("delrt",(short)Math.rint(t_start*1e3));
				// and samples
				for (int j=0;j<ns;j++) 
				  {
					_m._SI.tr_far_synth.set_data_value((float)padded_synth_far[j],j); 
				  }
			  }
			_m._SI.tr_far_synth.set(_m._SI.MF.realisationNumberHeaderField,(short)realisation_number);	
			// only write state if insisted on
			if (write_state) _m._SI.suwriter.fputtr(_m._SI.far_synth_out_strm,_m._SI.tr_far_synth);
		  }
	  }
  }

  // Copy inline/xline words from src su trace to dest
  void setInlineXlineHeaderFields(SuSegyFormat tr_src, SuSegyFormat tr_dest)
  {
	tr_dest.set("ep",tr_src.get_int("ep"));
	tr_dest.set("cdp",tr_src.get_int("cdp"));
	tr_dest.set("tracl",tr_src.get_int("tracl"));
	tr_dest.set("fldr",tr_src.get_int("fldr"));
  }

  /**
	 Useful function to get estimate of length scale of wavelet.
	 Calculations here assume spectral form W(f) ~ f^2 exp(-f^2 lambda^2 / 2)
	 or, in time,   W(t) ~ (1-(2 pi t/lambda)^2) * exp(-(2 pi t)^2 / 2 lambda^2).
	 From frequency spectrum, lambda is defined from bandedge 
	 by lambda = 3/freq_at_bandedge.
  */
  double getRickerCorrelationLengthLambda()
  {
	return 3.0/f_at_band_edge;
  }

  /**
	 Useful function to get estimate of band edge of wavelet.
  */
  double getWaveletBandEdgeFrequency()
  {
	return f_at_band_edge;
  }

  void setSyntheticParameters(double tail_f, double dt)
  {
	if (tail_f<0)
	  throw new RuntimeException("Synthetic seismics must have +ve frequency");
	this.f_at_band_edge = tail_f;
	this.suggested_sampling_dt = dt;
  }

  /** 
	  Make Ricker wavelet for near wavelet.
	  Calculations here assume spectral form W(f) ~ f^2 exp(-f^2 lambda^2 / 2)
	  or, in time,   W(t) ~ (1-(2 pi t/lambda)^2) * exp(-(2 pi t)^2 / 2 lambda^2).
	  From frequency spectrum, lambda is defined from bandedge 
	  by lambda = 3/freq_at_bandedge.
	  From suggested freq_at_bandedge, we compute lambda, a sensible sampling rate dt, and
	  make the wavelet etc.
  */
  private void makeRickerWavelet()
  {
	int i,j,N;
	double spacing,dt,lambda;
	double x,ls,t;
	double norm;

	_m._SI.EH.print("%nCreating synthetic wavelet%n");

	lambda = 3.0/f_at_band_edge;
	// corresponding Nyquist rate
	dt = 1.0/(2.0*f_at_band_edge);

	// find canonical sampling time just less than this (1ms,2ms,4ms etc)
	if (dt<_m._SI.FSC.BASE_SCALING_TIME)
	  {
		throw new RuntimeException("Desired synthetic seismic has unreasonably high bandwidth");
	  }
	for (t=_m._SI.FSC.BASE_SCALING_TIME ;; t*=2.0)
	  {
		if (t>dt) break;
		if (t>1000*_m._SI.FSC.BASE_SCALING_TIME) 	
		  {
			throw new RuntimeException("Desired synthetic seismic has unreasonably low bandwidth");
		  }
	  }
	// desired sampling time is half of last attempt, since it failed
	dt = 0.5*t;
	// compare to users suggested time
	if (dt < 0.99 * suggested_sampling_dt)
	  {
		_m._SI.EH.print("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is too long "
						+"for bandwidth desired. I suggest "+(1000*dt)+"ms.%n"
						+"I will use " + (1000*dt)+"ms despite what you suggest.%n");
		
	  }
	else if ((dt > 1.01 * suggested_sampling_dt) &&  (dt < 4 * suggested_sampling_dt))
	  {
		_m._SI.EH.print("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is needlessly "
						+"short for bandwidth desired. I suggest "+(1000*dt)+"ms%n"
						+"I'll do what you say anyway.%n");
		// override sensible value with users idea
		dt = suggested_sampling_dt;
	  }
	else if (dt > 4 * suggested_sampling_dt)
	  {
		throw new RuntimeException("Suggested sampling time " + (1000*suggested_sampling_dt)+"ms is "
								   +"absurdly short for bandwidth desired. I suggest "+(1000*dt)+"ms.%n");
	  }
	else
	  {
		// a good suggestion
		dt = suggested_sampling_dt;
	  }

	// wavelet is dead (to 3*10^-6) by t = 0.9 * lambda, so
	N = MathUtil.fCeil_int(0.9*lambda/dt);
	spacing = lambda/(2.0*Math.PI*dt);

	
	near_wavelet = new Wavelet((short)(2 * N - 1), (short)(N-1));
	far_wavelet = near_wavelet;

	ls = sqr(spacing);

	for (i=0;i<N;i++)
	  {
		x=sqr(i);
		near_wavelet.W[i]=(1.0-x/ls)*Math.exp(-0.5*x/ls);
	  }
	// wraparound bit
	for (j=2*N-2,i=1;i<N;i++,j--)
	  {
		near_wavelet.W[j]=near_wavelet.W[i];
	  }
	for (norm=0.0,i=0;i<near_wavelet.getLength();i++)
	  {
		norm+=sqr(near_wavelet.W[i]);
	  }
	norm=Math.sqrt(norm);
	for (i=0;i<near_wavelet.getLength();i++)
	  {
		near_wavelet.W[i]/=norm;
	  }
	near_wavelet.f1 = -(N-1.0)*dt;
	// sample time in mu_s
	near_wavelet.dt = (short)(1e6 * dt);

	computeBestErrorSamplingTime(near_wavelet);

	writeOutWavelet(near_wavelet);
  }

 /**
	 Write a wavelet for all output stacks
  */
  private void writeOutWavelet(Wavelet wavelet)
  {
	writeOutWavelet(wavelet,null);
  } 

 /**
	 Write a Ricker wavelet to both near and far wavelet_filename files
  */
  private void writeOutWavelet(Wavelet wavelet,String suffix)
  {
	int i,j;
	SuSegyFormat tr = new SuSegyFormat();
	FileOutputStream near_out, far_out=null;
	_m._SI.suwriter = new SuWriter();

	// open the wavelet file. Assume the same wavelet for nears and fars
	int Nstacks=2;
	String[] wavelet_filenames = new String[Nstacks];
	try 
	  {
		wavelet_filenames[0] = nearstack_data.wavelet_filename + ((suffix!=null) ? suffix : "");
		near_out = (FileOutputStream) new FileOutputStream(wavelet_filenames[0]);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open" + wavelet_filenames[0]);
	  }
	if (farstack_data!=null)
	  {
		try
		  {
			wavelet_filenames[1] = farstack_data.wavelet_filename + ((suffix!=null) ? suffix : "");
			far_out = (FileOutputStream) new FileOutputStream(wavelet_filenames[1]);
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open" + wavelet_filenames[1]);
		  }
	  }
	// set requisite fields

	tr.set("ns",wavelet.getLength());	
	tr.set("f1",(float)wavelet.f1);
	tr.set("dt",wavelet.dt);
	tr.set("delrt",(short)Math.rint(wavelet.f1*1e3));

	// and wavelet samples
	int N,Npre,Ncoda;
	N = wavelet.getLength();
	Npre = (int)Math.rint(-1e6*wavelet.f1/wavelet.dt);
	Ncoda = N-Npre;
	// time half t<0
	for (j=Ncoda,i=0; i<Npre; i++,j++)
	  {
		tr.set_data_value((float)wavelet.W[j],i);
	  }    
	// time half t>=0
	for (j=0;j<Ncoda;j++,i++)
	  {
		tr.set_data_value((float)wavelet.W[j],i); 
	  }

	_m._SI.suwriter.fputtr(near_out,tr);
	if (farstack_data!=null) _m._SI.suwriter.fputtr(far_out,tr);

	// near wavelet
	try
	  {
		_m._SI.suwriter.flush(near_out); // must do this!
		near_out.close();
		if (_m._SI.verbosity>0) 
		  _m._SI.EH.println("%nWrote synthetic wavelet to: " + wavelet_filenames[0]);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close" + wavelet_filenames[0]);
	  }
	// far wavelet
	if (farstack_data!=null)
	  {
		try
		  {
			_m._SI.suwriter.flush(far_out); // must do this!
			far_out.close();
			if (_m._SI.verbosity>0) 
			  _m._SI.EH.println("%nWrote synthetic wavelet to: " + wavelet_filenames[1]);
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close" + wavelet_filenames[1]);
		  }
	  }
	tr=null;
	far_out=null;
	near_out=null;
	_m._SI.suwriter=null;
  }

  /**
	 Read in a wavelet from both near and far wavelet_filename files
  */
  private void readInWavelets()
  {
	int stack_cntr;
	SuReader sureader = new SuReader();
	SuSegyFormat tr = new SuSegyFormat();
	StackData stack_data=null;
	Wavelet wavelet=null;

	for (stack_cntr=0;stack_cntr<=1;stack_cntr++)
	  {
		if (stack_cntr==0)
		  {
			stack_data = nearstack_data;
		  }
		else if (stack_cntr==1)
		  {
			stack_data = farstack_data;
		  }
		// there may not be far-stack data
		if (stack_cntr==1 && stack_data==null) continue;
		
		switch (stack_cntr)
		  {
		  case 0:
			wavelet = near_wavelet = getWavelet(sureader, tr ,stack_data.wavelet_filename);
			computeBestErrorSamplingTime(near_wavelet);
			break;
		  case 1:
			wavelet = far_wavelet = getWavelet(sureader, tr ,stack_data.wavelet_filename);
			break;
		  }
		// verify if tAbove or tBelow entries are sensible wrt the wavelet coda and precursor
		// The rules are:
		// 0 <= tAbove <= wavelet_coda
		// 0 <= tBelow <= wavelet_precursor
		if (!Double.isNaN(stack_data.getTabove()))  // t-above specified in XML
		  {
			if (stack_data.getTabove() > wavelet.getDTcoda())
			  {
				throw new RuntimeException("t_above value for stack cannot be greater than wavelet coda = "+
										   (1000*wavelet.getDTcoda())+" ms");
			  }
		  }
		else
		  {
			// set to DTCoda value (sensible default)
			stack_data.setTabove(wavelet.getDTcoda());
		  }
		if (!Double.isNaN(stack_data.getTbelow())) // t-below specified in XML
		  {
			if (stack_data.getTbelow() > wavelet.getDTprecursor())
			  {
				throw new RuntimeException("t_below value for stack cannot be greater than wavelet precursor = "+
										   (1000*wavelet.getDTprecursor())+ " ms");
			  }
		  }
		else
		  {
			// set to DTprecursor value (sensible default)
			stack_data.setTbelow(wavelet.getDTprecursor());
		  }

	  }
  }
  
  /**
	 Read in and set up a wavelet from the named wavelet file, given
	 already setup SuReader and SuSegyFormat fields
  */
  Wavelet getWavelet(SuReader sureader, SuSegyFormat tr, String wavelet_filename)
  {
	double delta_t,f1;
	int i,j,final_indx;
	FileInputStream in;
	Wavelet wavelet=null;
	// open the wavelet file. Assume the same wavelet for nears and fars
	try
	  {
		in = (FileInputStream) new FileInputStream(wavelet_filename);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + wavelet_filename);
	  }
	if (sureader.fgettr(in,tr)<1)
	  throw new RuntimeException("No wavelet data in " + wavelet_filename);
	// simple tests on critical fields
	if ((ns = (int) tr.get_short("ns"))<2)
	  throw new RuntimeException("Wavelet is unacceptably short " + wavelet_filename);
	if ((short)tr.get_short("dt")==0)
	  throw new RuntimeException("Wavelet has weird dt=0 in " + wavelet_filename);
	delta_t = (double)((short)tr.get_short("dt")/1.0e6);

	if ((f1 = (double)tr.get_short("delrt")) != 0) // delrt in ms
	  {
		f1 *= 1e-3;
	  }
	else 
	  {
		f1 = tr.get_float("f1");
	  }


	// f1/delta_t should be an integer, but round to be sure fp error doesn't kill us
	final_indx = ns -1 + (int)Math.rint(f1/delta_t);

	wavelet = new Wavelet((short)ns,(short)final_indx);

	wavelet.dt = (short)tr.get_short("dt");
	wavelet.f1 = f1;

	// Get wavelet samples.
	// In trace, wavelet is in ordinary form.
	// Copy to wrap-around form
	// time half t<0
	for (j=ns-final_indx-2,i=ns-1;j>=0;i--,j--)
	  {
		wavelet.W[i] = (double) tr.get_data_value(j);
	  }
	// time half t>=0
	for (j=ns-final_indx-1,i=0 ; j<ns ; j++,i++)
	  {
		wavelet.W[i] = (double) tr.get_data_value(j);
	  }
		

	if (hasInverterContext && _m._SI.waveletTrimEpsilon > 0)
	  {
		// trim wavelets
		int imax=0,imin=0,final_indx_new,ns_new;
		double u,eps = _m._SI.waveletTrimEpsilon;
		for(i=0,u=0;i<wavelet.W.length;i++) u = Math.max(u, wavelet.W[i]);
		wavelet.max_abs = u;
		u *= eps;
		// coda
		for (j=ns-1; j>=0; j--) if (Math.abs(tr.get_data_value(j)) > u) {imax=j; break;}
		// precursor
		for (j=0; j<ns;j++) if (Math.abs(tr.get_data_value(j)) > u) {imin=j; break;}

		final_indx_new = imax - (ns - final_indx - 1);
		ns_new = imax-imin+1;
		// overwrite wavelet
		wavelet = new Wavelet((short)ns_new,(short)final_indx_new);
		wavelet.dt = (short)tr.get_short("dt");
		f1 = -((ns-final_indx-1) - imin)*delta_t;
		wavelet.f1 = f1;
		for (j=imin,i=final_indx_new+1; j<(ns-final_indx-1); i++,j++)
		  {
			wavelet.W[i] = (double) tr.get_data_value(j);
		  }
		// time half t>=0
		for (j = ns-final_indx-1,i=0 ; i<=final_indx_new ; j++,i++)
		  {
			wavelet.W[i] = (double) tr.get_data_value(j);
		  }	
		if (_m._SI.verbosity>0) 
		  {
			_EH.printf("Wavelets trimmed with eps = %.2g%n",eps);
		  }			
	  }
	//writeOutWavelet(wavelet,".short.su");

	try
	  {
		in.close();
		if (hasInverterContext)
		  {
			if (_m._SI.verbosity>0) 
			  {
				_EH.println("Read wavelet from " + wavelet_filename);
			  }
		  }
		else
		  {
			_EH.warn("Read wavelet from " + wavelet_filename);
		  }
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + wavelet_filename);
	  }
	return wavelet;
  }





  private void computeBestErrorSamplingTime(Wavelet wavelet)
  {
	// estimate best sampling time for independent samples
	// for a signal generated by delta-functions convolved with
	// this wavelet
	double[] W, padded_wavelet;

	W=wavelet.W ;

	// Take FFT of wavelet to estimate best sampling rate
	int n=W.length;
	int ng=0;
	while ((n >>= 1) > 0 ) ng++;
	// In general W.length not a power of 2, so find the right n = 2^ng  bigger than W.length
	if (W.length != (1<<ng))
	  {
		// W.length not a power of 2
		n=(1<<(ng+1));
	  }
	else
	  n = W.length; // it was already a power of 2
	// now pad out W
	padded_wavelet = new double[n];
	int wl = wavelet.getLength();
	int wfi = wavelet.getFinalIndex();
	// +ve time part of wavelet
	System.arraycopy(wavelet.W,0,padded_wavelet,0,wfi+1);
	// -ve time part
	System.arraycopy(wavelet.W,wfi+1,
						 padded_wavelet,wfi + 1 + n - wl,
						 wl - 1 - wfi);


	FFT.realftZO(padded_wavelet,n,1);

	int i,i_dft,i_max=0;
	double i_peak,spec, f_peak, max=0.0;
	// Ordinary sampling time
	delta_t = ((double)wavelet.dt*1.0e-6);

	// detect peak wavelet amplitude for SN comparisons	
	for (i=0,max=0.0;i<W.length;i++) max = Math.max(max,Math.abs(W[i]));
	wavelet.max_abs = max;

	max=0.0;
	if (!true)
	  {
		// set the peak frequency from the absolute peak
		for (i=2;i<=n/2;i+=2)
		  {
			spec = sqr(padded_wavelet[i])+sqr(padded_wavelet[i+1]);
			if (spec>=max)
			  {
				max = spec;
				i_max=i;
			  }
		  }
		// take peak at literal peak
		i_peak = i_max/2.0;
	  }
	else
	  {
		// set the peak frequency from the fft-weighted index (first moment)
		// this is the "mean" frequency. If the spectrum is hump shaped, this will be a good guess
		double fft[] = new double[n/2];
		double freqs[] = new double[n/2];
		double num=0,den=0, freq_den=((double)n * delta_t);
		// OLD for (i=2,i_dft=1;i<=n/2;i+=2,i_dft++)
		for (i=2,i_dft=1; i<=n-2; i+=2,i_dft++)
		  {
			spec = sqr(padded_wavelet[i])+sqr(padded_wavelet[i+1]);
			if (spec>=max)
			  {
				max = spec;
				i_max=i_dft;
			  }
			fft[i_dft]=spec;
			freqs[i_dft]=(double)i_dft/freq_den;
			num+=spec*i_dft;  // first moment
			den+=spec;    // normalisation
		  }
		i_peak = (num/den)/MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER;
		_m._SI.EH.println("Wavelet peak detected at index "+i_max/2);
		_m._SI.EH.println("Wavelet mean-peak detected at effective index %4g",i_peak/2.0);
		if (_m._SI.verbosity>4) 
		  {
			Util.printfileZO("FFT_wavelet_spectrum.txt",freqs,fft,fft.length);
		  }
		fft=null;   // GC
		freqs=null; // GC
	  }

	
	f_peak = i_peak/((double)n * delta_t);
	if (_m._SI.verbosity>0)
	  {
		_m._SI.EH.println("Wavelet peak frequency=%3g Hz",f_peak);
	  }

	// Set error-sampling time of all synthetics
	delta_t_indep_samples = SAMPLE_TIME_MAGIC_CONSTANT/f_peak;

	// delta_t_indep_samples = delta_t; _m._SI.EH.println("REMOVE!!"); 

	// This is just totally unacceptable: data must be badly sampled
	if (delta_t_indep_samples < 0.3*delta_t)
	  throw new RuntimeException("wavelet spectrum is wierd! Independent-samples"
								 +" sampling time is less than sampling-time");

	delta_t_indep_nearest = delta_t * MathUtil.fFloor(delta_t_indep_samples/delta_t);
	if (delta_t_indep_samples < delta_t)
	  {
		// This is bad, so let's winge a lot, but keep going
		_m._SI.EH.println("WARNING! Wavelet spectrum has enough high frequency content%n"+
						  "that independent-samples sampling time is less than sampling-time");
		_m._SI.EH.println("Desired error-acquisition dt = %5g ms, will use %5g ms",
						  1000*delta_t_indep_samples, 1000*delta_t);
		// snap it to the dt and live with it
		delta_t_indep_nearest = delta_t;
	  }
	else
	  {
		if (_m._SI.verbosity>0)
		  {
			_m._SI.EH.println("Inferred error-acquisition dt = %5g ms, using corrected subsampling at %5g ms",
							  1000*delta_t_indep_samples, 1000*delta_t_indep_nearest);
		  }
	  }

	


	delta_t_ratio_correction_factor = Math.sqrt((delta_t_indep_nearest/delta_t_indep_samples));
	// seismic_error_factor = delta_t_indep_nearest/delta_t_indep_samples;
  }


  /**
	 Convenience method for computing a set of "exclusion time windows" based
	 on a nominated start layer number and "end" layer number. Assume 
	 named layer numbers are legal, i.e.
	 1 <= layer_number <= n_layers
  */
  void computeExclusionZones()
  {
	int i;
	double SIGMA_MARGIN = 3.0;  // safety margin number of std deviations
	double t, t_marg_coda, t_marg_pre, sig_t;
	boolean thisLayerUseExclusionZone;
	StackedModelInfo.SeismicExclusionLayers[] seismicExclusionLayers;

	seismicExclusionLayers = _m.stackedModelInfo.seismicExclusionLayers;

	// assume no exclusion unless sensible zone can be defined from inputs
	useExclusionZone = false;
	if (seismicExclusionLayers != null)
	  {
		t_start_exclusion = new double[seismicExclusionLayers.length];
		t_end_exclusion =  new double[seismicExclusionLayers.length];

	
		t_marg_coda =  near_wavelet.getDTcoda();
		if (far_wavelet != null) t_marg_coda = Math.max(t_marg_coda, far_wavelet.getDTcoda());
		t_marg_pre =  near_wavelet.getDTprecursor();
		if (far_wavelet != null) t_marg_pre = Math.max(t_marg_pre, far_wavelet.getDTprecursor());

		for (i=0; i< seismicExclusionLayers.length; i++)
		  {
			thisLayerUseExclusionZone = false;

			// start at first reflector minus maximal DTcoda over all wavelets
			if (seismicExclusionLayers[i].seismic_exclusion_start_layer<=0)
			  {
				// make this inactive
				t_start_exclusion[i] = Double.MAX_VALUE;
			  }
			else
			  {	
				t = _m.MPD.MetaInfo[seismicExclusionLayers[i].seismic_exclusion_start_layer-1].t.value - t_marg_coda;
				sig_t = _m.MPD.MetaInfo[seismicExclusionLayers[i].seismic_exclusion_start_layer-1].sigma_t.value;
				t -= SIGMA_MARGIN * sig_t;
		
				t_start_exclusion[i] = t;

				thisLayerUseExclusionZone = true;
				useExclusionZone = true;
			  }

			// end time
			if (seismicExclusionLayers[i].seismic_exclusion_end_layer > _m.n_layers)
			  {
				// make inactive
				t_end_exclusion[i] = -Double.MAX_VALUE;
			  }
			else
			  {		
				if (seismicExclusionLayers[i].seismic_exclusion_end_layer < _m.n_layers)
				  {
					// use top time of next layer
					t = _m.MPD.MetaInfo[seismicExclusionLayers[i].seismic_exclusion_end_layer].t.value + t_marg_pre;
					sig_t = _m.MPD.MetaInfo[seismicExclusionLayers[i].seismic_exclusion_end_layer].sigma_t.value;
				  }
				else
				  {
					// use base time of last layer
					t = _m.MPD.t_base.value + t_marg_pre;
					sig_t = _m.MPD.sigma_t_base.value;
				  }
				t += SIGMA_MARGIN * sig_t;
				t_end_exclusion[i] = t;

				thisLayerUseExclusionZone = true;
				useExclusionZone  = true;
			  }
			if (thisLayerUseExclusionZone)
			  {
				if (_m._SI.verbosity>2) 
				  {
					if (t_start_exclusion[i] < t_end_exclusion[i])
					  {
						_EH.printf("Warning: seismic exclusion zone %d: time window: [%6g, %6g]%n",(i+1),
								   t_start_exclusion[i],t_end_exclusion[i]);
					  }
					else
					  {
						_EH.printf("Warning: seismic exclusion zone %d: time window inoperative%n",i+1);
					  }
				  }
			  }
		  }
	  }
  }

  // return true if the time t falls in any of the exclusion zones, used
  // in the Stacked-Model facility
  private boolean isTimeinAnExclusionZone(double t)
  {
	int i;
	for (i=0; i<t_start_exclusion.length; i++)
	  {
		if ((t > t_start_exclusion[i]) && (t < t_end_exclusion[i])) 
		  {
			return true;
		  }
	  }
	return false;
  }


  public static final class Wavelet
  {
	/** 
		length of wavelet
	*/
	private short size;

	/** 
		wavelet coefficients
	*/
	double[] W;
	  
	/**
	   sample index of final positive time sample
	*/
	short final_index;

	/**
		 wavelet sampling time in microseconds
		 (is a short header field)
	  */
	short dt;
	  
	/**
		 Commencement time of wavelet in seconds
		 (same as su's f1: a <b>float</b> in the headers)
	  */
	double f1;
	/**
	   Maximum absolute sample
	*/
	double max_abs;

	public Wavelet(short size, short final_index)
	{
	  if (size<1)
		throw new RuntimeException("Wavelets must be created with length >= 1");
	  if (final_index >= size)
		throw new RuntimeException("Wavelets must be created with final index < size");
	  W = new double[size];
	  this.size = (short) size;
	  this.final_index = final_index;
	}

	short getLength()
	{
	  return size;
	}
	short getFinalIndex()
	{
	  return final_index;
	}
	/**
	   precursor, or lead time of wavelet in seconds, 
	   (i.e. extent of influence for times > t_reflection)
	   Valid only after dt field has been set.
	*/
	double getDTprecursor()
	{
	  return (double) (size-1-final_index)*dt*1e-6;
	}
	/**
	   coda length of wavelet in seconds;  (i.e. extent of influence for times < t_reflection)
	   Valid only after dt field has been set.
	*/
	double getDTcoda()
	{
	  return (double)(final_index)*dt*1e-6;
	}
  }
}
