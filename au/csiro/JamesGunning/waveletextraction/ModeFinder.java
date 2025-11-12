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

import java.io.*;
import java.util.*;

import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import cern.colt.matrix.*;
import cern.colt.list.*;
import cern.jet.math.*;


import au.csiro.JamesGunning.lib.Numerical.ModifiedCholeskyLeftDecomposition;
import au.csiro.JamesGunning.lib.Minpack.optimization.*;
import au.csiro.JamesGunning.lib.util.*;

/**
   Optimising class for finding modes for a given wavelet
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/


public class ModeFinder implements Uncmin_methods 
{
  // Link to calling extractor
  private WaveletExtractor _WE;

  ModelVariables modelVariables;

  /**
	 Skip BFGS and optimise using Gauss-Newton
  */
  public static boolean DO_GAUSS_NEWTON_ONLY = true;

  /**
	 be selective about remappings (omit noise)
  */
  public static boolean DO_PARTIAL_LINEAR_REMAPPING = true;

  /**
	 Set this to use parametrisation (additive squares) which guarantees
	 monotonicity in the t2d map
  */
  boolean reparametriseToEnforceT2DMonoticity = !true;
  /**
	 Absorb Jeffrey's prior into likelihood and use weak Gaussian prior
	 onn log(noise)
  */
  public static final boolean ABSORB_JEFFREYS_PRIOR_INTO_LIKELIHOOD = true;

  /**
	 Enable use of approximate Bayesian information criterion when comparing modes
	 of different dimensionality (instead of formal Bayes factors from Laplace method) 
  */
  public static final boolean ENABLE_BIC_USE = true;

  /**
	 Epsilon to be used in finite differences. Forward
	 difference suggest this ought to be O(sqrt(10*machine epsilon)), where
	 the 10 comes from likely propogated roundoff in the numerical model [d=O(100) parameters
	 with roundoff diffusion going like sqrt(d)]. Postscript: experiments have shown that
	 O(10^-5)-O(10-6) seem to be most robust. Change from this value WITH EXTREME CAUTION.
  */
  public static final double DESIGN_MATRIX_EPS=1e-6;

  /**
	 Mode searcher will retry this many Newton refinements
	 of a mode if faulty Newton convergence is detected. If still
	 faulty, we give up.
  */
	static public int MAX_NUM_FAULTY_NEWTON_REFINEMENTS = 30;

  /**
	 The line search and subsequently the Newton method sometimes 
	 terminates at the bottom of a "cliff-face" causes by
	 schizophrenic peak detection. Handle such a termination by adding
	 large penalties.
  */
  static public boolean HANDLE_AMBIGUOUS_PEAK_DETECTION_TERMINATION_USING_PENALTY = true;
 /**
	 The line search and subsequently the Newton method sometimes
	 terminates at the bottom of a "cliff-face" causes by
	 schizophrenic peak detection. This is the penalty we add
	 large penalties.
  */
  static public double AMBIGUOUS_PEAK_DETECTION_TERMINATION_PENALTY = 0;

  /**
	 Maximum scaled step size of parameter vector in Gauss-Newton method.
  */
  public static double MAX_STEP_SIZE=20.0;

  /**
	 Maximum number of Newton steps used in Newton methods of computing the posterior 
	 covariance and final solution
  */
  public static int MAX_NEWTON_STEPS = 50;

  // Variable for Newton optimisations
  /**
	 Epsilon for Newton method terminations. Usually f=O(d), and we stop optimising when
	 |delta f | < 0.0001 |f|. This is the 0.0001
  */
  public static double NEWTON_EPS = 0.0001;
  /**
   Alternative epsilon for Newton method terminations. We can stop optimising when
   |delta f | < 0.005, which is within statistically reasonable accuracy. 
   This is the 0.005
  */
  public static double ABS_NEWTON_EPS = 0.005;
  /**
	 There is an epsilon for termination of Newton methods associated with floating point
	 rounding too
  */
  public static double EPS_FLOATINGPOINT=1e-8;
  /**
	 Newton steps with line search step less than
	 FRAC_SUSPICIOUS_NEWTON_STEP of the full Newton step and yielding improvements
	 less than EPS_SUSPICIOUS_NEWTON_STEP (absolute) will fire off a steepest
	 descent line search as well.
  */
  public static double EPS_SUSPICIOUS_NEWTON_STEP = 0.05;
  public static double FRAC_SUSPICIOUS_NEWTON_STEP = 0.15;
  /**
	Exponential reparametrisation of noise requires safeguarding against over/underflow using this...
  */
  public static double LOG_NOISE_SAFETY_MARGIN = 10;
  /**
	Exponential reparametrisation of noise requires safeguarding against over/underflow using this...
	This limits noise < e^LOG_NOISE_UPPER_SAFETY_MARGIN * guessed value, which is reasonable
  */
  public static double LOG_NOISE_UPPER_SAFETY_MARGIN = 3;
  
  /**
	 Number of steps in Xsections when they are printed
  */
  public static int N_XSECTION_STEPS = 20;

  /**
	 Threshold to print dodgyXsection files if the cross section
	 through the max. likelihood (ML) point appears to have a point of
	 lower -log(posterior) than the ML point. The candidate point must
	 be DODGY_XSECTION_THRESHOLD less than the ML point to initiate a
	 dodgyXsection printout.
  */
  public static double DODGY_XSECTION_THRESHOLD = 0.01;

  /**
	 Threshold to print dodgyXsection files if the cross section
	 through the max. likelihood (ML) point appears to have a discontinuous Xsection
	 through the ML point. The candidate point must
	 be DODGY_XSECTION_JUMP_THRESHOLD different than the ML point to initiate a
	 dodgyXsection printout at either side of the ML point.
  */
  public static double DODGY_XSECTION_JUMP_THRESHOLD = 0.5;

  /**
	 The log(noise) variables have a weak Gaussian prior with this variance
  */
  public static double PRIOR_INVERSE_VARIANCE_FOR_LOG_NOISE = 1/1000.0;
  /**
	 Occasionally the optimiser thinks it can produce zero noise. For stability we have
	 to threshold this.
  */
  public static double MINIMUM_TOLERABLE_SN_RATIO = 1e-6;


  /**
	 This allows the stringing together of the other (non synthetic related)
	 constraints in the Bayesian updating formulae
  */
  boolean absorbOtherConstraintsInLikelihood = true;
  /**
	 After an optimisation step, if the problem is degenerate (noise ->0), this will be set to true
  */
  boolean problemHasDegenerateNoise = false;

  boolean doingMultiStartSolutions=false;
  
  private int optimiserPass=0; // keeps track of which optimisation peass we are in
  private double fNewNewton; // value at new point
  private double fNewNewtonOld; // value at old point
  private double delta_fNewNewton; // change from old point to new point
  private DoubleMatrix2D HessianPostMinimiser;

  /* Peak arrival times in line search */
  private double tPeakStartLineSearch;
  private double tPeakEndLineSearch;



  private int it_count=0;
  private int multistart_ctr=0;
  private boolean abortNow=false;

  private double negLogSeisLikelihood;       // these are set by each f_to_minimise call
  private double negLogPriorLikelihood;      // these are set by each f_to_minimise call


  // Control variables for solver.
  // see dfault_f77 for example of how this is done in Uncmin....

  double[] f_best = new double[2];             // f[1] has the function minimum at exit
  int[] info = new int[2];                     // status information about the code
  double[] fscale = new double[2];
  int[] method = new int[2];
  int[] is_exp = new int[2];
  int[] msg = new int[2];
  int[] n_good_digits = new int[2];
  int[] max_its = new int[2];
  int[] is_gradient_supplied = new int[2];
  int[] is_hessian_supplied = new int[2];
  double[] trust_radius = new double[2];
  double[] gradient_tol = new double[2];
  double[] max_step_size = new double[2];
  double[] step_tol = new double[2];

  // workspace variables
  private double[][] Hessian_wk;
  private double[] Hessian_diag_wk;
  private double[] fullParameters_best;
  private double[] grad_sub_best;
  private double[] typical_values;
  private double[] fullParameters;

  private double[] noiseLevelsLinear;


  // Typical scales, normalised for nosie transformations
  private double[] typScale;

  // possible communication of line search directions
  private DoubleMatrix1D ls_fullParameters;

  // Debugging. Switch both off (production)
  private boolean loadFaultyNewtonParams = !true;
  private boolean dumpFaultyNewtonParams = !true;
  
  /** 
	  Constructor
  */
  public ModeFinder(WaveletExtractor WE, ModelVariables modelVariables)
  {
	this._WE = WE;
	this.modelVariables = modelVariables;
  }


  /**
	 Find optimal solution, starting from prior mean guess of all parameters
  */
  
  public void findOptimalSolution()
  {
	int i,n_tot;
	double fMin;
	boolean bestPeakArrivalPositive=false;
	ConvolutionalModel _convModel;
	ModelVariables.ModelVariablesPosteriorSnapshot mVSnapshot=null;

	// shortcut
	_convModel = _WE.convModel;

	// to suboptimise for wavelet-coefficients first, 
	// before global optimisation using depth parameters & other free parameters
	boolean applySubsetoptimisationFirst = true;    

	// optimiser setup
	// Tolerances...
	double epsilon  = 1.12e-16; // default 1.12e-16;
	// FIXME can the gradient tolerance can be weakened??????
	gradient_tol[1] = Math.pow(epsilon,1.0/3.0);
	max_step_size[1] = 0;     // 0.0 default
	step_tol[1]=Math.sqrt(epsilon);
	trust_radius[1]=-1;   // -1 is default

	// Flags for methods/outputs...
	method[1]=1;          // 1 = line-search
	// NB. UNCMIN can be a bit wierd under BFGS if scaling is dodgy. I believe scaling issues are
	// fixed so BFGS (Secant) is safe
	is_exp[1]=(_WE.useSecantOptimisation ? 1 : 0);  
	msg[1]=8;             // msg=8 seems to suppress most of the noise... :-). msg=0,1 is very noisy
	if (_WE.verbosity>=5) msg[1]=0;
	n_good_digits[1]=-1;  // -1 gives defaults
	max_its[1]=1000;
	is_gradient_supplied[1]=0;
	is_hessian_supplied[1]=0;


	problemHasDegenerateNoise = false;


	if (applySubsetoptimisationFirst)
	  {
		// make submodel consisting only of wavelet coefficients 
		modelVariables.initialiseModel(true);
	  }
	else
	  {
		// single optimise on full model only
		modelVariables.initialiseModel(false);
	  }

	// Set size of system
	n_tot = modelVariables.nParams;
	noiseLevelsLinear = new double[_WE.m.nStacks];
	allocateWorkspace(n_tot);

	// initial starting values and scales

	// FIXME get these right. Note unit offsetting
	for (i=0; i<n_tot; i++)
	  {
		if (reparametriseToEnforceT2DMonoticity && 
			i>modelVariables.firstTDparamIndex && i<=modelVariables.lastTDparamIndex)
		  {
			// loop over subrange of parameters pertaining to t2d conversion.
			// additive squares guarantees monotonicity
			for (i=modelVariables.firstTDparamIndex+1; i<=modelVariables.lastTDparamIndex; i++)
			  {
				fullParameters[i+1] = Math.sqrt(modelVariables.initialGuess[i] 
												- modelVariables.initialGuess[i-1]);
				typical_values[i+1] = fullParameters[i+1];
			  }
			i--; // so we don't double increment
		  }
		else
		  {
			// ordinary properties. We optimise using the difference from the initial guess
			fullParameters[i+1] = 0.0;
			// Set typical scaling. take care here as it affects convergence
			if (i<modelVariables.nGaussian)
			  {
				typical_values[i+1] = modelVariables.typicalScale[i];
			  }
			else
			  {
				// This is used to set initial Hessian.
				typical_values[i+1] = Math.log(modelVariables.typicalScale[i]);
			  }
		  }
	  }

	fscale[1]=(double)n_tot;  // assumes chisq \sim n_tot if properly scaled
	
	it_count=0;
	optimiserPass=1;
	
	// Now do it !
	check_f_to_minimize(fullParameters,false);
	if (applySubsetoptimisationFirst)
	  {
		// Get initial solution with fixed T2D, AVO etc. This initial solution will 
		// likely be decent if there is limited T2D freedom, registration effects etc, but possibly very poor 
		// if a strong phase constraint applies and we're not using a smart method. 
		if (!_WE.constrainPhase || (_WE.constrainPhase && ModelVariables.CONSTRAIN_PHASE_VIA_PHASE_LINE))
		  {
			if (_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign)  
			  {
				// make first solution the "natural" solution, unhampered by time shifts
				// Achieve this by switching off this part of the likelihood
				_WE.peakArrivalStdDev *= 1e20;
			  }
			// Let's take the wavelet-coefficients problem to be fully linear in the coefficients, and
			// use the Bayesian least squares problem. This analytic solution
			// is faster than optif9_f77() for the approximate first pass
			completeModeCalculation(true,true);
			// Solution is now in-place in *unmapped* in-situ mode variable mu_post
			// Extract mode variables to raw Parameters
			reverseMapModeParametersToRawParameters(fullParameters_best);
			// Map raw parameters to in-situ variables
			mapParameters(fullParameters_best,true);
			// and thence from in-situ variables to offset parameters
			reverseMapParameters(fullParameters_best,true);
			fMin = f_best[1] = f_to_minimize(fullParameters_best);

			if (_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign)  _WE.peakArrivalStdDev /= 1e20; // restore tolerance
		  }
		else  // use  Uncmin for subset optimisation
		  {
			if (_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign)  
			  {
				// make first solution the "natural" solution, unhampered by time shifts
				// Achieve this by switching off this part of the likelihood
				_WE.peakArrivalStdDev *= 1e20;
			  }
			Uncmin_f77.optif9_f77(n_tot, fullParameters, this,
								  typical_values,fscale,method,is_exp,msg,n_good_digits,
								  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
								  gradient_tol,max_step_size,step_tol,
								  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
			// if a premature abort, get out quickly
			if (abortNow) 
			  {
				_WE.EH.println("UNCMIN: Preliminary Abort(). FIXME");
				return;
			  }
			if (info[1]>=4)
			  {
				_WE.EH.println("UNCMIN: Preliminary subset optimisation did not converge!!!");
			  }
			// Ensure MAP parameters are in place
			fMin = f_to_minimize(fullParameters_best);
			if (_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign)  _WE.peakArrivalStdDev /= 1e20; // restore tolerance
		  }
		if (_WE.printMinimiserXsections)
		  {
			_WE.peakArrivalStdDev *= 1e20;
			debugPrintMinimiserXsections("postMinimiserXsection.TDfixed."
										 +_WE.multiModelSetUp.modeNumber+".",
										 fullParameters_best);
			_WE.peakArrivalStdDev /= 1e20;
		  }
		if (!_WE.ignoreSeismic)
		  {
			// update noise estimates if not ignoring seismic
			for (i=0; i<_WE.m.nStacks; i++) 
			  {
				noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);
				if (noiseLevelsLinear[i] < (MINIMUM_TOLERABLE_SN_RATIO*
											_WE.multiModelSetUp.typicalRMSSeismicAmplitude))
				  {
					_WE.EH.warn("Warning: noise estimate underflow - thresholding for safety");
					noiseLevelsLinear[i] = MINIMUM_TOLERABLE_SN_RATIO
					  *_WE.multiModelSetUp.typicalRMSSeismicAmplitude;
				  }	
				if (!modelVariables.fixNoise)
				  {
					_WE.m.seismic_noise[i].setValue(_WE.ignoreSeismic ? 1.0 : noiseLevelsLinear[i]);
				  }
			  }
		  }
	  }
	else
	  {
		Uncmin_f77.optif9_f77(n_tot, fullParameters, this,
							  typical_values,fscale,method,is_exp,msg,n_good_digits,
							  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
							  gradient_tol,max_step_size,step_tol,
							  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
		// if a premature abort, get out quickly
		if (abortNow) 
		  {
			_WE.EH.println("UNCMIN: Preliminary Abort(). FIXME"); return;
		  }
		if (info[1]>=4)
		  {
			_WE.EH.println("UNCMIN: Solution did not converge!!!");
		  }
		// Ensure true parameters are in place
		f_best[1] = f_to_minimize(fullParameters_best);
	  }

	// and compute commensurate noise level estimates
	for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);

	
	if (_WE.verbosity>3) 
	  { 
		_WE.EH.println("Wavelet span: (time span: %-10.5g %-10.5g). Coarse samples span %2d %2d, "
					   +"internal degrees of freedom: %2d", -_WE.m.wavelet.precursorDuration,
					   _WE.m.wavelet.codaDuration, _WE.m.wavelet.firstCoarseSample,
					   _WE.m.wavelet.finalCoarseSample,
					   _WE.m.wavelet.valuesCoarseProp.length);
		_WE.EH.println("Minimiser its = "+it_count+", fmin = "+ f_best[1]+ ", solution:");
		for (i=0; i<modelVariables.nParams; i++) 
		  {
			_WE.EH.printf("%10.8g ",((Property)modelVariables.mapFromIndex.get(i)).value);
		  }
		_WE.EH.println("");
		if (modelVariables.fixNoise)
		  {
			_WE.EH.printf("Noise level(s) (linear estimate)");
			for (i=0; i<_WE.m.nStacks; i++) 
			  {
				_WE.EH.printf("%10.8g ",noiseLevelsLinear[i]);
			  }
			_WE.EH.printf("%nFixed (guessed) noise level: %10.8g%n",_WE.m.seismic_noise[0].value);
		  }
		else
		  {
			_WE.EH.printf("Noise level(s) (linear estimate)");
			for (i=0; i<_WE.m.nStacks; i++) 
			  {
				_WE.EH.printf("%10.8g ",noiseLevelsLinear[i]);
			  }
			_WE.EH.printf("%nMinimiser Noise level(s) %n");
			for (i=0; i<_WE.m.nStacks; i++) _WE.EH.printf("%10.8g ",_WE.m.seismic_noise[i].value);
			_WE.EH.printf("%n");
		  }
		_WE.EH.printf("Ndata (effective) points in seismic (stack 1) = %5.2g%n",_convModel.dimData[0]);
	  }

	


	if (applySubsetoptimisationFirst)
	  {
		// Second pass with true prior
		optimiserPass=2;

		// this will overwrite in-situ fields with default guesses
		modelVariables.initialiseModel(false);
		//.......Set size of system......
		n_tot = modelVariables.nParams;
		allocateWorkspace(n_tot);
		// initial starting values and scales
		// FIXME get these right. Note unit offsetting
		modelVariables.recomputeInitialGuessFromInSituFields();
		initialiseParametersAndInitialValues(fullParameters);
		
		if (DO_GAUSS_NEWTON_ONLY)
		  {
			// Skip following BFGS optimisation. But ensure fullParameters_best is properly setup			

			// Fix to include multi-start

			copy(fullParameters,fullParameters_best);
			fMin = completeModeCalculation(false,false);     // STD LINEARISED VERSION
			// Then remap from in-situ variables to offset parameters
			reverseMapParameters(fullParameters_best,true);

			if (_WE.constrainPeakArrival && !_WE.constrainPeakArrivalSign)
			  {
				// agnostic peak solution requested, but did it solve OK?
				if (modelVariables.modeHasPeakAmbiguity)
				  {
					// No: it came out ambiguous. Try again, using both styles of peak sign
					if (_WE.verbosity>3) 
					  {
						_WE.EH.println("Ambiguous peak-arrival constraint: attempting sign-specific arrival-time constraint");
					  }

					double fMin_h = Double.MAX_VALUE;
					double[] fullParameters_best_h = fullParameters_best.clone();
					boolean peakArrivalPositive;
					bestPeakArrivalPositive=false;

					modelVariables.constrainPeakArrivalSign = true;
					for (peakArrivalPositive = false;; peakArrivalPositive=true)
					  {
						// Loop over both styles of peak detection if it's required and user 
						// is agnostic, otherwise
						// use the default and break out after the first pass through this loop
						modelVariables.peakArrivalPositive = peakArrivalPositive;

						if (_WE.verbosity>3) 
						  {
							_WE.EH.println("Attempting sign-specific arrival-time constraint with sign "+(peakArrivalPositive?"+ve":"-ve"));
						  }
						// initial starting values and scales again
						initialiseParametersAndInitialValues(fullParameters);
						copy(fullParameters,fullParameters_best);

						fMin = completeModeCalculation(false,false);     // STD LINEARISED VERSION
						// Then remap from in-situ variables to offset parameters
						reverseMapParameters(fullParameters_best,true);
						// Compute commensurate noise level estimates
						for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);

						// Keep track of best solution
						if (fMin<fMin_h)
						  {
							fMin_h = fMin;
							copy(fullParameters_best,fullParameters_best_h);
							bestPeakArrivalPositive = peakArrivalPositive;
						  }
						if (false)
						  {
							// write out multistart wavelets - debugging
							_WE.m.wavelet.writeOutSUWavelet("MLwavelet.Multistart.base."+(peakArrivalPositive?"p":"n"),true,true,
															_WE.multiModelSetUp.tMinLongestWavelet,
															_WE.multiModelSetUp.tMaxLongestWavelet,false,0,false);
						  }
						// on -ve peak, take snapshot in case we have to restore it
						if (!peakArrivalPositive) mVSnapshot = modelVariables.takePosteriorSnapshot();
						// escape on second pass
						if (peakArrivalPositive) break;
					  }

					if (!bestPeakArrivalPositive)
					  {
						// have to restore -ve peak solution
						// restore to condition of best solution
						modelVariables.peakArrivalPositive = bestPeakArrivalPositive;
						// restore best solution
						copy(fullParameters_best_h,fullParameters_best);
						// this will restore correct covariances, mu_post values etc etc
						mVSnapshot.restoreFromSnapshot(modelVariables);
						// evaluate...
						fMin = f_to_minimize(fullParameters_best);
						// Compute commensurate noise level estimates
						for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);
						// fMin_h = completeModeCalculation(false,true);  
						f_best[1] = fMin_h;
					  }
					if (_WE.verbosity>3) 
					  {
						_WE.EH.println("Best solution of sign-specific peak-arrival constraining = "+fMin_h);
					  }	
				  }
			  }
			// Update noise level estimates
			for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);
		  }
		else
		  {
			fscale[1]=(double)n_tot;  // assumes chisq \sim n_tot if properly scaled
	
			it_count=0;
			msg[1]=8; // Needs to be reset. msg=8 seems to suppress most of the noise... :-). msg=0,1 is very noisy
			if (_WE.verbosity>=5) msg[1]=0;
			is_gradient_supplied[1]=0; // Needs to be reset.
			is_hessian_supplied[1]=0;  // Needs to be reset.
		
			if (_WE.printMinimiserXsections)
			  {
				debugPrintMinimiserXsections("preMinimiserXsection."+_WE.multiModelSetUp.modeNumber+".",
											 fullParameters);
			  }

			Uncmin_f77.optif9_f77(n_tot, fullParameters, this,
								  typical_values,fscale,method,is_exp,msg,n_good_digits,
								  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
								  gradient_tol,max_step_size,step_tol,
								  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
			// if a premature abort, get out quickly
			if (abortNow) {_WE.EH.println("UNCMIN: Preliminary Abort(). FIXME");return; }
			if (info[1]>=4)
			  {
				_WE.EH.println("UNCMIN: Solution did not converge!!!");
			  }
			// Ensure MAP parameters are in place
			fMin = f_to_minimize(fullParameters_best);
			// and compute commensurate noise level estimates
			for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);

		  }



		// multiple random-point minimisations
		if (_WE.multiStart && !reparametriseToEnforceT2DMonoticity)
		  {
			doingMultiStartSolutions = true;
			// Random start points come from hypercube in those parameters likely to cause strong nonlinearity
			// which I suspect are positioning parameters.
			int dim=0,Np;
			
			// Accumulate list of parameters to multistart from
			IntArrayList indxList = new IntArrayList(0);
			// t-registration errors
			if (modelVariables.haveTseismicErrors)
			  {
				for (i=modelVariables.firstTregistrationParamIndex; i<= modelVariables.lastTregistrationParamIndex; i++)
				  {
					indxList.add(i);
				  }
			  }
			// XY minicube errors
			if (modelVariables.haveXYseismicErrors)
			  {
				for (i=modelVariables.firstXYpositioningParamIndex; i<=modelVariables.lastXYpositioningParamIndex; i++)
				  indxList.add(i);
			  }
			indxList.trimToSize();

			int[] hyperIndices = indxList.elements();
			dim = hyperIndices.length;

			if (dim!=0)
			  {
				float[][] offsets;
				Np = (int)Math.max(1.0, Math.floor(Math.pow(WaveletExtractor.MAX_GLOBAL_OPT_STARTING_POINTS, 1.0/dim)));
				// should never need more than 5...
				Np = Math.min(Np,5); 
				// Np=50; // for basins plots

				offsets = au.csiro.JamesGunning.lib.Numerical.Hypercube.createHypercubeList(dim,Np);

				_WE.EH.println("Multistart solutions: (numbers of starts = "+offsets.length+") (number of params = "+dim+")");
				int j;
				double scaleInit, fMin_h;
				double[] fullParameters_best_h = fullParameters_best.clone();

				boolean peakArrivalPositive;
			    int multistart_ctr_glbl=0,multistart_ctr_best=0;

				// Take snapshot of best solution so far
				mVSnapshot = modelVariables.takePosteriorSnapshot();
				fMin_h = fMin;

				for (peakArrivalPositive = false;; peakArrivalPositive=true)
				  {
					// Loop over both styles of peak detection if it's required and user 
					// is agnostic, otherwise
					// use the default and break out after the first pass through this loop
					if (_WE.constrainPeakArrival && !_WE.constrainPeakArrivalSign)
					  {
						// User wants peak arrival alignment, but is unsure of sign
						modelVariables.peakArrivalPositive = peakArrivalPositive;
						modelVariables.constrainPeakArrivalSign = true;
						_WE.EH.println("Peak detection sign positive = "+ peakArrivalPositive);
					  }
					// loop over points in hypercube
					for (multistart_ctr=0; multistart_ctr<offsets.length; multistart_ctr++,multistart_ctr_glbl++)
					  {
						// Set start point according to current hypercube point
						for (i=0,j=0; i<n_tot; i++)
						  {
							if ((j<dim) && (i == hyperIndices[j]))
							  {
								// a hypercube param. Offsets are in [0..1], scale to [-2...2]
								// start points are now the initial guess +/- [-2...2] * typical_values
								scaleInit = 4 * (offsets[multistart_ctr][j]-0.5);
								j++; // wait for next hypercube index
							  }
							else scaleInit = 0; // don't scatter other dimensions

							if (DO_PARTIAL_LINEAR_REMAPPING)
							  {
								if (i>=modelVariables.nGaussian)
								  {
									// log-noise parameterisation
									fullParameters[i+1] = Math.log(Math.abs(modelVariables.initialGuess[i]));
								  }
								else
								  {
									if (i>=modelVariables.firstTregistrationParamIndex && i<= modelVariables.lastTregistrationParamIndex)
									  {
										// the typical *scale* for rego errors is too narrow to get well
										// separated start points. But the prior uncertainty sigma_{t,reg} is often useful,
										// and a quarter-loop is also a good idea.
										// Since C_{ii}=sigma_{t,reg}^2, and
										// Quarter loop of the time interval 0.25/_WE.multiModelSetUp.f_peak
										double tReg = Math.min(0.25/_WE.multiModelSetUp.f_peak,
															   Math.sqrt(modelVariables.C_prior.getQuick(i,i)));
										fullParameters[i+1] = scaleInit * tReg;
									  }
									else
									  {
										// ordinary properties. Shifted by their initial guess
										fullParameters[i+1] = scaleInit * typical_values[i+1];	
									  }
								  }
							  }
							else
							  {
								// ordinary properties. Shifted by their initial guess
								fullParameters[i+1] = scaleInit * typical_values[i+1];
							  }
						  }

						// Fire off local minimiser
						it_count=0;


						if (DO_GAUSS_NEWTON_ONLY)
						  {
							// Ensure fullParameters_best is properly setup						
							copy(fullParameters,fullParameters_best);

							fMin = completeModeCalculation(false,false);     // STD LINEARISED VERSION
							// Then remap from in-situ variables to offset parameters
							reverseMapParameters(fullParameters_best,true);
						  }
						else
						  {
							msg[1]=8; // Needs to be reset. msg=8 seems to suppress most of the noise... :-). msg=0,1 is very noisy
							if (_WE.verbosity>=5) msg[1]=0;
							is_gradient_supplied[1]=0; // Needs to be reset.
							is_hessian_supplied[1]=0;  // Needs to be reset.
		
							_WE.EH.printf("Multistart init vals= ");
							for (i=0; i<dim; i++) _WE.EH.printf("%10.8g ", fullParameters[hyperIndices[i]+1]);
							_WE.EH.println("");

							Uncmin_f77.optif9_f77(n_tot, fullParameters, this,
												  typical_values,fscale,method,is_exp,msg,n_good_digits,
												  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
												  gradient_tol,max_step_size,step_tol,
												  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
							// if a premature abort, get out quickly
							if (abortNow) {_WE.EH.println("UNCMIN: Preliminary Abort(). FIXME");return;}
							if (info[1]>=4)
							  {
								_WE.EH.println("UNCMIN: Solution did not converge!!!");
							  }
							// Ensure MAP parameters are in place
							fMin = f_to_minimize(fullParameters_best);
						  }

						if (_WE.verbosity>3) 
						  {
							_WE.EH.printf("Multistart solution: its=%d val=%10.8g%n",it_count,fMin);						
						  }

						// Compute commensurate noise level estimates
						for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]); 

						// Keep track of best solution
						if (fMin<fMin_h)
						  {
							fMin_h = fMin;
							copy(fullParameters_best,fullParameters_best_h);
							bestPeakArrivalPositive = peakArrivalPositive;
							mVSnapshot = modelVariables.takePosteriorSnapshot();
							multistart_ctr_best = multistart_ctr_glbl;
						  }

						if (false)
						  {
							// write out multistart wavelets
							_WE.m.wavelet.writeOutSUWavelet("MLwavelet.Multistart."+(modelVariables.peakArrivalPositive?"p.":"n.")
															+multistart_ctr_glbl,true,true,
															_WE.multiModelSetUp.tMinLongestWavelet,
															_WE.multiModelSetUp.tMaxLongestWavelet,false,0,false);
							_WE.m.wavelet.getPeakTime(false,modelVariables.constrainPeakArrivalSign,modelVariables.peakArrivalPositive);
						  }

					  } // loop over Hypercube points

					if (!_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign)
					  {
						// Don't care about peak arrival, or we do care, but have specified the sign
						// => only one pass needed
						break;
					  }
					// break out on second pass
					if (peakArrivalPositive) break;
				  } // loop over peak-arrival sign

				// make sure modelVariables is left in state corresponding to best solution
				if (_WE.constrainPeakArrival && !_WE.constrainPeakArrivalSign)
				  {
					modelVariables.peakArrivalPositive = bestPeakArrivalPositive;
				  }

				multistart_ctr = 0; // reset for debugging graphs

				// restore best solution
				copy(fullParameters_best_h,fullParameters_best);
				// this will restore correct covariances, mu_post values etc etc
				mVSnapshot.restoreFromSnapshot(modelVariables);
				// evaluate...
				fMin = f_to_minimize(fullParameters_best);
				f_best[1] = fMin;
				// f_best[1] = fMin_h = completeModeCalculation(false,true);  
				// Compute commensurate noise level estimates
				for (i=0; i<_WE.m.nStacks; i++) noiseLevelsLinear[i] = Math.sqrt(_convModel.SSerrors[i]/_convModel.dimData[i]);
				if (_WE.verbosity>3) 
				  {
					_WE.EH.println("Best multistart solution was attempt "+multistart_ctr_best+" with minimum "+fMin_h);
				  }

				fullParameters_best_h = null;
				offsets=null;
			  }
			doingMultiStartSolutions = false;
		  }

		// Construct Hessian at MAP point by finite diffs
		if (_WE.USE_FULL_NUMERICAL_HESSIAN) 
		  {
			constructHessianByFDs(fullParameters_best,Hessian_wk);
		  }

		if (_WE.printMinimiserXsections)
		  {
			debugPrintMinimiserXsections("postMinimiserXsection."+_WE.multiModelSetUp.modeNumber+".",
											 fullParameters_best);
		  }

		if (_WE.verbosity>3) 
		  {
			_WE.EH.println("Minimiser its = "+it_count+", fmin = "+ f_best[1]+ ", solution:");
			for (i=0; i<modelVariables.mapFromIndex.size(); i++) 
			  {
				_WE.EH.printf("%10.8g ",((Property)modelVariables.mapFromIndex.get(i)).value);
			  }
			_WE.EH.println("");
			if (!modelVariables.fixNoise)
			  {
				_WE.EH.printf("Noise level(s) (linear estimate)");
				for (i=0; i<_WE.m.nStacks; i++) 
				  {
					_WE.EH.printf("%10.8g ",noiseLevelsLinear[i]);
				  }
				_WE.EH.printf("%nMinimiser noise level(s) ");
				for (i=0; i<_WE.m.nStacks; i++) _WE.EH.printf("%10.8g ",_WE.m.seismic_noise[i].value);
				_WE.EH.printf("%n");
			  }
			else
			  {
				// Classical formula for linear regression
				_WE.EH.printf("Noise level(s) (linear estimate)");
				for (i=0; i<_WE.m.nStacks; i++) 
				  {
					_WE.EH.printf("%10.8g ",noiseLevelsLinear[i]);
				  }
				_WE.EH.printf("%nFixed (guessed) noise level: %10.8g%n",_WE.m.seismic_noise[0].value);
				_WE.EH.printf("Ndata (effective) points = %5.2g%n",_convModel.dimData[0]);
			  }
		  }
	  }
	// flag degenerate noise estimate
	for (i=0; i<_WE.m.nStacks; i++) 
	  {
		if (noiseLevelsLinear[i] <= MINIMUM_TOLERABLE_SN_RATIO *_WE.multiModelSetUp.typicalRMSSeismicAmplitude)
		  {
			problemHasDegenerateNoise = true;
		  }
	  }

	// complete setup of mode (covariance etc)
	if (DO_GAUSS_NEWTON_ONLY)
	  {
		// mode likely to be already available
		if (_WE.multiStart && !reparametriseToEnforceT2DMonoticity)
		  {
			// must recompute if multistart used; but only once
			if (!_WE.USE_FULL_NUMERICAL_HESSIAN)
			  {
				completeModeCalculation(false,true);     // STD LINEARISED VERSION
			  }
			else
			  {
				completeModeCalculationFromBFGS(f_best[1]);   // TO BE DEBUGGED
			  }
		  }
	  }
	else
	  {
		if (!_WE.USE_FULL_NUMERICAL_HESSIAN)
		  {
			// full iterative version
			completeModeCalculation(false,false);         // STD LINEARISED VERSION
		  }
		else
		  {
			completeModeCalculationFromBFGS(f_best[1]);   // TO BE DEBUGGED
		  }
	  }

	// GC
	// Hessian_wk = null; require this from now
	Hessian_diag_wk = null;
	grad_sub_best = null;
	noiseLevelsLinear = null;
  }

  /**
	 simple wrapper for matching vector copies
  */
  void copy(double[] from,double[] to)
  { 
	System.arraycopy(from,0,to,0,from.length);
  }

  /**
	 Construct Hessian by finite differences. Standard formula 7.21 of Nocedal and
	 Wright. The epsilon seems big because (1) the LogDet term in the prior
	 causes cancellation issues if ep ~ 10^{-NDIGITS/3), as it ought
	 to be. (2) Also, there are discretisation issues which wreck the
	 first derivative (piecewise continuity only). So I reckon NDIGITS
	 = about 6
  */
  private void constructHessianByFDs(double[] fullParametersBest, double[][] hessGuess)
  {
	int i,j,nT;
	double val,f_MAP,hessEntry,ep_i,ep_j;
	double EPSILON_FD=1e-2;
	double[] params = new double[fullParametersBest.length];
	double[] f_perm;


	nT = modelVariables.nParams;


	double[] untranslatedParams = fullParametersBest.clone();

	// ensure in-situ state is set....
	mapParameters(fullParametersBest,true); 
	// set up raw parameters vector
	for (i=0;i<(modelVariables.nParams);i++) 
	  {    
		untranslatedParams[i+1] = ((Property)modelVariables.mapFromIndex.get(i)).value;
		// last entries in untranslatedParams are log(noise) values
		if (i>=modelVariables.nGaussian) untranslatedParams[i+1] = Math.log(Math.abs(untranslatedParams[i+1]));
	  }	

	// Make copy into working vector
	copy(untranslatedParams,params);


	HessianPostMinimiser = new DenseDoubleMatrix2D(nT,nT);

	if (true)
	  {
		f_MAP = f_to_minimize_raw(params);
		f_perm = new double[nT];
		double[] ep = new double[nT];
		if (false)
		  {
			// hessGuess comes in as L factor of Hessian after BFGS finishes.
			// Compute diag elements ep[i] from (L.L^T)_{ii}
			for (i=1;i<=nT;i++) 
			  {
				for (j=1;j<=i;j++) ep[i-1] += hessGuess[i][j] * hessGuess[i][j];
			  }
		  }
		// initial run down diagonal with rough guess of ep to get suitable refined scaling ep[i]
		for (i=0;i<nT;i++)
		  {
			ep[i] = modelVariables.typicalScale[i] * EPSILON_FD ;
			params[i+1] = untranslatedParams[i+1] + ep[i];
			f_perm[i] = f_to_minimize_raw(params);
			params[i+1] = untranslatedParams[i+1];
		  }

		for (i=0;i<nT;i++)
		  {
			params[i+1] = untranslatedParams[i+1] + 2*ep[i];
			val = f_to_minimize_raw(params);
			hessEntry = (val - 2*f_perm[i] + f_MAP)/(ep[i]*ep[i]);
			// _WE.EH.println("HESS-app "+i+" "+hessEntry);
			// recompute suitable ep from Hessian estimate
			// ep[i] = (hessEntry>0) ? EPSILON_FD/Math.sqrt(hessEntry) : ep[i];
			ep[i] = (hessEntry>0) ? 0.3/Math.sqrt(hessEntry) : ep[i];   // FIXME!!!!!!!!!!!!!
			// restore
			params[i+1] = untranslatedParams[i+1];
		  }
		// now fully compute Hessian
		for (i=0;i<nT;i++)
		  {
			params[i+1] = untranslatedParams[i+1] + ep[i];
			f_perm[i] = f_to_minimize_raw(params);
			params[i+1] = untranslatedParams[i+1];
		  }

		for (i=0;i<nT;i++)
		  {
			params[i+1] = untranslatedParams[i+1] + ep[i];
			for (j=0;j<=i;j++)
			  {
				if (j<i) 
				  {
					params[j+1] = untranslatedParams[j+1] + ep[j];
				  }
				else params[j+1] = untranslatedParams[j+1] + 2*ep[j];
				val = f_to_minimize_raw(params);
				hessEntry = (val - f_perm[i] - f_perm[j] + f_MAP)/(ep[i]*ep[j]);
				//_WE.EH.println("HESS "+i+" "+j+" "+Math.max(Math.abs(val - f_MAP),Math.abs(f_perm[j]-f_MAP))/val);
				// if (i==j) _WE.EH.println("HESS "+i+" "+hessEntry+" "+params[j+1]+" "+untranslatedParams[j+1]+" "+(f_perm[i] - f_MAP));
				HessianPostMinimiser.setQuick(i,j,hessEntry);
				HessianPostMinimiser.setQuick(j,i,hessEntry);
				// restore
				params[j+1] = untranslatedParams[j+1];
			  }
			params[i+1] = untranslatedParams[i+1];
		  }
		// restore MAP value
		f_to_minimize_raw(params);
		// Util.print(ep,nT);
		ep=null;
	 	//Util.debugPrintMatrixToFile("HessianFD."+_WE.multiModelSetUp.modeNumber+".txt",HessianPostMinimiser);
	  }
	else
	  {
		// old code that is known to work, but slower
		double numerator;
		f_MAP = f_to_minimize_test(params,untranslatedParams);
		for (i=0;i<nT;i++)
		  {
			ep_i = modelVariables.typicalScale[i] * EPSILON_FD; 
			for (j=0;j<=i;j++)
			  {
				ep_j = modelVariables.typicalScale[j] * EPSILON_FD;
				/* term 1 */
				params[i+1] = untranslatedParams[i+1]+ep_i+((i == j)?ep_j:0.0);
				params[j+1] = untranslatedParams[j+1]+ep_j+((i == j)?ep_i:0.0);
				numerator =  f_to_minimize_test(params,untranslatedParams);
				/* term 2 */
				params[i+1] = untranslatedParams[i+1]+ep_i-((i == j)?ep_j:0.0);
				params[j+1] = untranslatedParams[j+1]-ep_j+((i == j)?ep_i:0.0);
				numerator -=  f_to_minimize_test(params,untranslatedParams);
				/* term 3 */
				params[i+1] = untranslatedParams[i+1]-ep_i+((i == j)?ep_j:0.0);
				params[j+1] = untranslatedParams[j+1]+ep_j-((i == j)?ep_i:0.0);
				numerator -=  f_to_minimize_test(params,untranslatedParams);
				/* term 4 */
				params[i+1] = untranslatedParams[i+1]-ep_i-((i == j)?ep_j:0.0);
				params[j+1] = untranslatedParams[j+1]-ep_j-((i == j)?ep_i:0.0);
				numerator +=  f_to_minimize_test(params,untranslatedParams);
				/* form matrix elements. Use symmetry */
				hessEntry = numerator/(4.0*ep_i*ep_j);  /* see Gelman 9.4 */
				HessianPostMinimiser.setQuick(i,j,hessEntry);
				HessianPostMinimiser.setQuick(j,i,hessEntry);
				/* restore params before changing rows/cols */
				params[j+1] = untranslatedParams[j+1];
			  }
			params[i+1] = untranslatedParams[i+1];
		  }
		// restore MAP value
		f_to_minimize_test(untranslatedParams,untranslatedParams);
		Util.debugPrintMatrixToFile("HessianFD."+_WE.multiModelSetUp.modeNumber+".txt",HessianPostMinimiser);
	  }
  }

  // test function for FD Hessian
  private double f_to_minimize_test(double[] x,double[] xopt)
  {
	int i,n=x.length-1;
	double s=0;
	if (true)
	  {
		for (i=0;i<n;i++)
		  {
			s += (x[i+1]-xopt[i+1])*(x[i+1]-xopt[i+1])/(modelVariables.typicalScale[i]*modelVariables.typicalScale[i]);
		  }
		for (i=0;i<n-1;i++)
		  {
			s += (x[i+1]-xopt[i+1])*(x[i+2]-xopt[i+2])/(modelVariables.typicalScale[i]*modelVariables.typicalScale[i+1]);
		  }
	  }
	else
	  {
		for (i=8;i<10;i++)
		  {
			s += (x[i+1]-xopt[i+1])*(x[i+1]-xopt[i+1])/(modelVariables.typicalScale[i]*modelVariables.typicalScale[i]);
		  }
	  }
	return s/2.0;
  }



  /**
	 The minimiser etc will sometimes drive the parameters to wierd values. We cope
	 with this in various ways in the minimiser, but we would like the base case starting values to be sensible.
	 Here is where we check this
  */
  void check_f_to_minimize(double[] parameters,boolean raw)
  {
	// switch on checking, and call to see if there's an exception
	_WE.convModel.switchOnSoftMinicubeExceptionHandling(false);
	if (raw)
	  {
		f_to_minimize_raw(parameters);
	  }
	else
	  {
		f_to_minimize(parameters);
	  }
	// switch off checking
	_WE.convModel.switchOnSoftMinicubeExceptionHandling(true);
  }

  /**
	 Minimisation function for set of remapped parameters
	 Here parameters[] is the unit-offset suite of remapped parameters we wish to optimise
  */
  public double f_to_minimize(double parameters[])
  {
	double chisq;

	// set in situ parameters with remapping
	chisq = mapParameters(parameters,true);
	chisq += f_to_minimize_kernel();

	return chisq;
  }

  /**
	 Minimisation function for set of raw parameters
	 raw_parameters[] (still unit-offset). Returns negative
	 log of posterior of model.
  */
  public double f_to_minimize_raw(double raw_parameters[])
  {
	double chisq;

	// set in situ parameters with remapping
	chisq = mapParameters(raw_parameters,false);
	chisq += f_to_minimize_kernel();
 
	return chisq;
  }


  /**
	 Raw kernel to evaluate minimisation function, assuming 
	 all parameters have been mapped in place
  */

  public double f_to_minimize_kernel()
  {
	int i;
	double stretch, chisq=0.0;
	double timesMisOrderingPenalty;

	// update time-depth map
	timesMisOrderingPenalty =_WE.m.updateTimeToDepthMaps();

	if (!_WE.m.timesAreDecentlyOrdered)
	  {
		_WE.EH.warn("warning: bad time ordering");
	  }

	if (_WE.useStretchedModel)
	  {
		// remap wav. stretch parameters
		for (i=0; i< (_WE.m.nStacks-1); i++)
		  {
			if (_WE.m.farStackWaveletTimeStretchMapper!=null) 
			  {
				_WE.m.farWavStretch[i] = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(_WE.m.farStackWaveletTimeStretch[i].value);
			  }
			else
			  {
				_WE.m.farWavStretch[i] = _WE.m.farStackWaveletTimeStretch[i].value;
			  }
			
			_WE.m.farWavAmp[i] = _WE.m.farStackWaveletAmplitude[i].value;
		  }
	  }
	// update Wavelet
	_WE.m.wavelet.updateWavelet();
			
	// compute basic reflectivity arrays
	if (_WE.useSaturationModel) _WE.m.recomputeLogProperties(false); 
	computeFastReflectivities();

	if (_WE.useSaturationModel) 
	  {
		// copy reflectivities
		_WE.m.makeReflectivityCopy();
		// compute effect of sat change,
		_WE.m.recomputeLogProperties(true);
		// new reflectivities
		computeFastReflectivities();
		// then difference
		_WE.m.computeReflectivityDifference();
	  }
		
	// compute synth seismic 
	_WE.convModel.computeSyntheticSeismic();
		
	// Get its associated likelihood
	chisq += computeFastSyntheticNegLogLikelihood(true);


	// Plus non-synthetic constraints
	chisq += modelVariables.computeAuxiliaryConstraints();

	// prior contributions
	chisq += modelVariables.computePriorNegLogLikelihood(false);

	if (timesMisOrderingPenalty!=0)
	  {
		// absorb this also to try and get things back into order
		chisq += timesMisOrderingPenalty;
	  }

	// keep track of number of function calls
	it_count++;

	if (_WE.verbosity>5) 
	  {
		// call this to get chisqVelMismatch	chiSqPhase
		modelVariables.computeAuxiliaryConstraints();
		// summarise......
		_WE.EH.println("chisq contribs (timeMisorderingPenalty, synthNLL, velNLL, phaseNLL, peakNLL, prior) "+
					   timesMisOrderingPenalty+ " "+computeFastSyntheticNegLogLikelihood(true)
					   + " " +modelVariables.chisqVelMismatch+"  "+modelVariables.chiSqPhase + "  " + 
					   + modelVariables.chisqPeakMismatch
					   + " " +modelVariables.computePriorNegLogLikelihood(false));
	  }

	return chisq;
  }


  
  /**
	 Fast computation of reflectivities, using caching
  */
  void computeFastReflectivities()
  {
	boolean useCachedReflectivity=false;
	
	if (Model.USE_CACHED_REFLECTIVITY)
	  {
		useCachedReflectivity = canUseCache(modelVariables.reflectivityCachingTestProperties, 
											modelVariables.reflectivityCachingOldValues);
	  }
	// FIXME add saturations to caching etc
	if (_WE.useSaturationModel) useCachedReflectivity = false;

	_WE.m.computeReflectivities(useCachedReflectivity,false);
  }




  /**
	 Fast link to synthetic negative log likelihood, using possible amplitude caching
  */
  double computeFastSyntheticNegLogLikelihood(boolean workInLogNoiseParameterisation)
  {
	boolean useCachedAmplitude=false;
	if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
	  {
		useCachedAmplitude = canUseCache(modelVariables.amplitudeCachingTestProperties, modelVariables.amplitudeCachingOldValues);
	  }
	return _WE.convModel.computeSyntheticNegLogLikelihood(workInLogNoiseParameterisation,useCachedAmplitude);
  }


  /**
	 Method to find whether a particular cache can be used
  */
  private boolean canUseCache(ArrayList<Property> cacheCheckArray, double[] cacheCheckOldValues)
  {
	int i;
	boolean cacheIsCorrect;
	double val;
	// (1) check if cache checking properties are the same as the old ones

	cacheIsCorrect = true;      // assume so till proven otherwise
	for (i=0; i<cacheCheckOldValues.length; i++)
	  {
		val = (cacheCheckArray.get(i)).value;
		if (val != cacheCheckOldValues[i])
		  {
			cacheIsCorrect = false;
		  }
		// store new ones as we go
		cacheCheckOldValues[i] = val;
	  }

	// never skip the first call
	if (it_count==0)
	  {
		cacheIsCorrect = false;
	  }

	return cacheIsCorrect;
  }

  /** 
	  Resets it_count so functions are *fully* evaluated on next call
  */
  void resetCache()
  {
	it_count = 0;
  }

  /**
    Initialise suite of parameters
  */
  void initialiseParametersAndInitialValues(double[] parameters)
  {
	int i;

	for (i=0; i<modelVariables.nParams; i++)
	  {
		if (reparametriseToEnforceT2DMonoticity && 
			i>modelVariables.firstTDparamIndex && i<=modelVariables.lastTDparamIndex)
		  {
			// loop over subrange of parameters pertaining to t2d conversion.
			// additive squares guarantees monotonicity
			for (i=modelVariables.firstTDparamIndex+1; i<=modelVariables.lastTDparamIndex; i++)
			  {
				parameters[i+1] = Math.sqrt(modelVariables.initialGuess[i] 
												- modelVariables.initialGuess[i-1]);
				typical_values[i+1] = parameters[i+1];
			  }
			i--; // so we don't double increment
		  }
		else
		  {
			if (DO_PARTIAL_LINEAR_REMAPPING)
			  {
				if (i>=modelVariables.nGaussian)
				  {
					// log-noise parameterisation
					parameters[i+1] = Math.log(Math.abs(modelVariables.initialGuess[i]));
				  }
				else
				  {
					// ordinary properties. Shifted by their initial guess
					parameters[i+1] = 0.0;
				  }
			  }
			else
			  {
				// ordinary properties. Shifted by their initial guess
				parameters[i+1] = 0.0;
			  }
			// Set typical scaling. take care here as it affects convergence
			if (i<modelVariables.nGaussian)
			  {
				typical_values[i+1] = modelVariables.typicalScale[i];
			  }
			else
			  {
				// This is used to set initial Hessian.
				typical_values[i+1] = Math.log(modelVariables.typicalScale[i]);
			  }
		  }
	  }
  }

  /**
	 Set the properties in various objects to the values required by the current 
	 parameter set. parameters[] is of length modelVariables.nParams+1.
	 The parameters are added to the initialGuess field to produce the
	 true in-situ properties
  */
  private double mapParameters(double[] parameters, boolean offsetRemapParameters)
  {
	int i;
	double t,penalty=0;
	Property prop;

	if (offsetRemapParameters)
	  {
		for (i=0; i<modelVariables.nParams; i++)
		  {
			// note parameters is unit offset
			if (reparametriseToEnforceT2DMonoticity && 
				i>modelVariables.firstTDparamIndex && i<=modelVariables.lastTDparamIndex)
			  {
				// ADDITIVE SQUARES PARAMETERISATION TO GUARANTEE MONOTONICITY
				t = ((Property)modelVariables.mapFromIndex.get(modelVariables.firstTDparamIndex)).value;
				// loop over subrange of parameters pertaining to t2d conversion.
				// additive squares guarantees monotonicity
				for (i=modelVariables.firstTDparamIndex+1; i<=modelVariables.lastTDparamIndex; i++)
				  {
					prop = (Property)modelVariables.mapFromIndex.get(i);
					t += parameters[i+1]*parameters[i+1];
					prop.setValue(t);
				  }
				i--; // so we don't double increment
			  }
			else
			  {
				if (DO_PARTIAL_LINEAR_REMAPPING)
				  {
					if (i<modelVariables.nGaussian)
					  {
						// linear mapping for ordinary parameters
						prop = (Property)modelVariables.mapFromIndex.get(i);
						prop.setValue(parameters[i+1] + modelVariables.initialGuess[i]);
					  }
					else
					  {
						// Noise remapping in log-noise space
						prop = (Property)modelVariables.mapFromIndex.get(i);
						prop.setValue(safeNoiseValue(parameters[i+1],typical_values[i+1]));
					  }
				  }
				else
				  {
					// LINEAR REMAPPING RELATIVE TO INITIAL GUESS
					prop = (Property)modelVariables.mapFromIndex.get(i);
					if (i<modelVariables.nGaussian)
					  {
						prop.setValue(parameters[i+1] + modelVariables.initialGuess[i]);
					  }
					else
					  {
						prop.setValue(safeNoiseValue(parameters[i+1] + modelVariables.initialGuess[i],
														typical_values[i+1]));
					  }
				  }
			  }
		  }
	  }
	else
	  {
		// simple raw mapping, no shifts, no positivity guarantees
		// only noise may be remapped
		for (i=0; i<modelVariables.nParams; i++)
		  {
			prop = (Property)modelVariables.mapFromIndex.get(i);
			if (i>=modelVariables.nGaussian)
			  {	
				// care here to prevent over/underflow
				prop.setValue(safeNoiseValue(parameters[i+1],typical_values[i+1]));
			  }
			else
			  {
				prop.setValue(parameters[i+1]);
			  }
		  }
	  }
	penalty += checkEndTimes();
	return penalty;
  }



  /**
	 Under exponential transformation, noise cannot thrash about too much, or we are in trouble.
	 Regulate it here
  */
  private double safeNoiseValue(double logNoise, double logNoiseScale)
  {
	double sigma;

	if (logNoise < (logNoiseScale - LOG_NOISE_SAFETY_MARGIN))
	  {
		sigma = (logNoiseScale - LOG_NOISE_SAFETY_MARGIN);
	  }
	else if (logNoise > (logNoiseScale + LOG_NOISE_SAFETY_MARGIN))
	  {
	    sigma = (logNoiseScale + LOG_NOISE_UPPER_SAFETY_MARGIN);
	  }
	else
	  {
		sigma = logNoise;
	  }
	return Math.exp(sigma);
  }


  /**
	 Set the parameters[] vector according to the in situ *mode* variables mu_post 
  */
  private void reverseMapModeParametersToRawParameters(double[] parameters)
  {
	int i;

	for (i=0; i<modelVariables.nParams; i++)
	  {
		// mu_post is always in log-noise space
		if (i<modelVariables.nGaussian)
		  {
			parameters[i+1] = modelVariables.mu_post.getQuick(i);
		  }
		else
		  {
			// log-noise to actual noise
			parameters[i+1] = Math.exp(modelVariables.mu_post.getQuick(i));
		  }
	  }
  }

  /**
	 Set the parameters[] vector according to the in situ variables. 
	 An exact inverse of mapParameters()
  */
  private void reverseMapParameters(double[] parameters, boolean offsetRemapParameters)
  {
	int i;
	Property prop,prop2;

	if (offsetRemapParameters)
	  {
		for (i=0; i<modelVariables.nParams; i++)
		  {
			// note parameters is unit offset
			if (reparametriseToEnforceT2DMonoticity && 
				i>modelVariables.firstTDparamIndex && i<=modelVariables.lastTDparamIndex)
			  {
				// ADDITIVE SQUARES PARAMETERISATION TO GUARANTEE MONOTONICITY
				for (i=modelVariables.firstTDparamIndex+1; i<=modelVariables.lastTDparamIndex; i++)
				  {
					prop = (Property)modelVariables.mapFromIndex.get(i);
					prop2 = (Property)modelVariables.mapFromIndex.get(i-1);
					parameters[i+1] = Math.sqrt(Math.abs(prop.value - prop2.value));
				  }
				i--; // so we don't double increment
			  }
			else
			  {
				prop = (Property)modelVariables.mapFromIndex.get(i);
				if (DO_PARTIAL_LINEAR_REMAPPING)
				  {
					if (i<modelVariables.nGaussian)
					  {
						// linear mapping for ordinary parameters
						parameters[i+1] = prop.value - modelVariables.initialGuess[i];
					  }
					else
					  {
						// remapping for noise...
						parameters[i+1] = Math.log(Math.abs(prop.value));
					  }
				  }
				else
				  {
					// LINEAR REMAPPING RELATIVE TO INITIAL GUESS
					// prop.setValue(parameters[i+1] + modelVariables.initialGuess[i]);
					parameters[i+1] = prop.value - modelVariables.initialGuess[i];
					if (i>=modelVariables.nGaussian)
					  {
						parameters[i+1] = Math.log(prop.value) - modelVariables.initialGuess[i];
					  }
				  }
			  }
		  }
	  }
	else
	  {
		// simple raw mapping, no shifts, no positivity guarantees
		for (i=0; i<modelVariables.nParams; i++)
		  {
			prop = (Property)modelVariables.mapFromIndex.get(i);
			// prop.setValue(parameters[i+1]);
			parameters[i+1] = prop.value;
			if (i>=modelVariables.nGaussian)
			  {
				parameters[i+1] = Math.log(parameters[i+1]);
			  }
		  }
	  }
  }


  
  public double checkEndTimes()
  {
	Property prop;
	int i;
	double t,penalty=0;

	// Now apply some checks on out-of control time maps.
	// (1) Ensure they span the desired time range ***at least***
	// We should pick _WE.(tStart,tEnd) well within the time span of the
	// knots, but if the knot-sequence wanders inside this range (the prior
	// should prevent this, but it can't be guaranteed) we're in big trouble.
	// This hack might help: we add a huge penalty to bad knot configurations....
	for (i=0; i<_WE.m.nWells; i++)
	  {
		TimeToDepthMap t2dmap = _WE.m.timeToDepthMaps[i];
		// old ->	modelVariables.mapFromIndex.get(modelVariables.firstTDparamIndex);		
		prop = t2dmap.tProp[0];
		if (prop!=null)
		  {
			t = prop.value;
			if (t > _WE.tStart)
			  {
				// first knot has crept inside [tStart,tEnd]
				penalty += (t - _WE.tStart);
				prop.setValue(_WE.tStart);
			  }
		  }
		//	prop = (Property)modelVariables.mapFromIndex.get(modelVariables.lastTDparamIndex);
		prop = t2dmap.tProp[t2dmap.tProp.length-1];
		if (prop!=null)
		  {
			t = prop.value;
			if (t < _WE.tEnd)
			  {
				// last knot has crept inside [tStart,tEnd]
				penalty += (_WE.tEnd - t);
				prop.setValue(_WE.tEnd);
			  }
		  }
	  }
	// if (penalty!=0) _WE.EH.println("Iteration "+it_count+" penalty = "+penalty);
	return penalty;
  }

  /**
	 Bayesian update computed direct from BFGS Hessian marix.
  */
  void completeModeCalculationFromBFGS(double negLogPost)
  {
	int i,j,nG,nN;
	double val,sum;
	ModifiedCholeskyLeftDecomposition cholDec;
	DoubleMatrix2D Htmp;  // FIXME could remove

	// Number of Gaussian parameters
	nG = modelVariables.nGaussian;
	// Number of noise parameters
	nN = modelVariables.nParams - modelVariables.nGaussian;

	// embed the MAP solution
	modelVariables.mu_post = new DenseDoubleMatrix1D(nG+nN);
	for (i=0; i<modelVariables.nParams; i++) 
	  {    
		val = ((Property)modelVariables.mapFromIndex.get(i)).value;
		// Set local Gaussian model vector
		if (i<nG) 
		  {
			modelVariables.mu_post.setQuick(i,val);
		  }
		else
		  {// a noise parameter with log transformation
			val = Math.log(Math.abs(val));
			modelVariables.mu_post.setQuick(i,val);
		  }
	  }

	
  
	// get nearest safely PD Cholesky factor from F.D. Hessian 
	// use supersafe prescaling
	if (true)  Htmp = HessianPostMinimiser.like().assign(HessianPostMinimiser); // FIXME REMOVE

	cholDec = new ModifiedCholeskyLeftDecomposition(HessianPostMinimiser,true,false);

	//Util.debugPrintMatrixToFile("HessianFDafter."+_WE.multiModelSetUp.modeNumber+".txt",HessianPostMinimiser);

	modelVariables.L_post = cholDec.getL();
	modelVariables.L_Permutation = cholDec.getPermutation(); 

	// recompute safe Hessian, and store diagonal elements
	DoubleMatrix2D Cinv = HessianPostMinimiser.like();
	SeqBlas.seqBlas.dgemm(false,true,1.0,modelVariables.L_post, modelVariables.L_post,0,Cinv);
	modelVariables.CtildeInvDiag = new double[nG+nN];
	for (i=0; i<(nG+nN); i++)
	  {
		modelVariables.CtildeInvDiag[modelVariables.L_Permutation[i]] = Cinv.getQuick(i,i);
	  }
	//Util.print(modelVariables.CtildeInvDiag,(nG+nN));

	Cinv = null;

	// We need the Cholesky inverse. Do this in-place
	ModifiedCholeskyLeftDecomposition.getCholeskyLeftInverse(modelVariables.L_post);

	if (false)
	  {
		// Debug
		DoubleMatrix2D Ct = HessianPostMinimiser.like();
		DoubleMatrix2D result = HessianPostMinimiser.like();
		Util.debugPrintMatrixToFile("Linv."+_WE.multiModelSetUp.modeNumber+".txt",modelVariables.L_post);
		SeqBlas.seqBlas.dgemm(true,false,1.0,modelVariables.L_post, modelVariables.L_post,0,result);
		for (i=0; i<(nG+nN); i++)
		  {
			for (j=0; j<(nG+nN); j++)
			  {
				Ct.setQuick(modelVariables.L_Permutation[i],modelVariables.L_Permutation[j],result.getQuick(i,j));
			  }
		  }
		SeqBlas.seqBlas.dgemm(false,false,1.0,Htmp,Ct,0,result);
		Util.debugPrintMatrixToFile("HessianInv."+_WE.multiModelSetUp.modeNumber+".txt",Ct);
		Util.debugPrintMatrixToFile("HessianInvCheck."+_WE.multiModelSetUp.modeNumber+".txt",result);
	  }


	// Compute log(noise) variances while we're here
	modelVariables.CtildeDiag = new double[nG+nN];
	for (i=0; i<(nG+nN); i++)
	  {
		for(j=i,sum=0.0; j<(nG+nN); j++) sum += MathUtil.sqr(modelVariables.L_post.getQuick(j,i));
		// insert posterior variance in correct (de-permuted) index in Ctilde diagonal position
		modelVariables.CtildeDiag[modelVariables.L_Permutation[i]] = sum;
	  }
	//Util.print(modelVariables.CtildeDiag,(nG+nN));

	// Update other terms required for sampling
	// 1/2(Determinant of C_eff)
	// since C_eff = P.L_post^T.L_post.P^T and |P|=1
	double halfLogDetCeff = 0.0;

	for (i=0; i<(nG+nN); i++)
	  {
		halfLogDetCeff += Math.log(modelVariables.L_post.getQuick(i,i));
	  }

	modelVariables.halfLogDetCeff = halfLogDetCeff;  // FIXME don't really need this variable

	// Mode Laplace factor
	// negLogPost = f_to_minimize_raw(newParams);
	modelVariables.logModeWeight = halfLogDetCeff - negLogPost;

	// Here we could use the approximate Bayes Information Criterion BIC = 0.5*n_p*log(n_d)
	// approach. But this depends only on the likelihood...  
	// But since BIC is only asymptotic, the Laplace estimator already acquired is probably better......
	// The number of parameters differs only in the wavelet coefficients, 
	// so that's all we'll include here
	int nParams = _WE.m.wavelet.valuesCoarseProp.length;
	if (_WE.useBIC && ENABLE_BIC_USE)
	  {
		if (_WE.ignoreSeismic)
		  {
			// if we ignore Seismic, it is as if we make no observations,and thus 
			// there should be no BIC (all models are equiprobable)
			modelVariables.logModeWeight = 0.0;
		  }
		else
		  {
			double nBIC=0.0;
			for (i=0;i<_WE.m.nStacks;i++) 
			  {
				// Since n_eff/N = dt/dtCoarse and
				// dimData = sqrt(n_eff/N) * N,  pull in the extra sqrt factor
				nBIC += _WE.convModel.dimData[i] * Math.sqrt(_WE.dt/_WE.m.wavelet.dtCoarse);
			  }
			// The standard BIC rule	
			modelVariables.logModeWeight = -computeFastSyntheticNegLogLikelihood(true)
			  - 0.5*nParams*Math.log(nBIC);
		  }
	  }

	modelVariables.noiseCorrCoeffMatrix = _WE.convModel.computeMultiStackNoiseCovarianceMatrix();

	if (_WE.verbosity>2)
	  {
		// FIXME this hack applies only for first stack......

		_WE.EH.println("(1)MODE (2)No.-params (3)Laplace-log(weight) (4)-log(posterior) (5)SSE(stack1) "+
					   "(6)-log(likelihood) (7)-log(prior) (8)noise(stack-1)-rms (9)halfLogDetCeff (10)halflogDetCp "+
					   "(11)halflogDetCp-halfLogDetCeff");
		_WE.EH.println("MODE: %3d %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g",nParams,
					   modelVariables.logModeWeight,negLogPost,
					   _WE.convModel.SSerrors[0],
					   computeFastSyntheticNegLogLikelihood(true),
					   modelVariables.computePriorNegLogLikelihood(false),
					   _WE.m.seismic_noise[0].value,
					   modelVariables.halfLogDetCeff,
					   modelVariables.logDetC_prior_only/2.0,
					   modelVariables.logDetC_prior_only/2.0 -modelVariables.halfLogDetCeff);
	  }
  }

  /**
	 Sensitivity matrix calculation, with noise parameters updated in
	 log-Gaussian approximation (so log(sigma_i) has Gaussian prior,
	 where the sigma_i are noise standard deviations).

	 Bayesian update etc done here.  The update is done in the
	 reparametrised space, so no issues with truncation etc ought to
	 occur.
  */
  double completeModeCalculation(boolean runQuietly, boolean onePassOnly)
  {
	Blas blas = SeqBlas.seqBlas;
	int i,j,k,k_b,l,nG,nN,nMeasurements,nSeisMeasurements,row,col,column;
	int i_trial;
	double epsScale,denom,delta,sum,val=0,f_0,f_prime=0,negLogPost=0,stp_size=0;
	double[] untranslatedParams;
	double[] newParams;
	double[][][] baseErrors;
	double[] baseAuxiliaryErrors=null;
	double[] tilde_dy_sq_sums=null;
	DoubleMatrix2D X,Xseismic;
	int[] blockXstarts=null, blockXsizes=null;
	DoubleMatrix2D C_prefactor;
	DoubleMatrix2D QtDY=null,QtDY_X=null;
	DoubleMatrix1D rhs,rhs1,rhs2=null;
	DoubleMatrix1D m_full;
	DoubleMatrix1D mG;          // model parameters for Gaussian components (first part of full vector)
	DoubleMatrix1D mN=null;     // current log(noise) parameters
	DoubleMatrix1D m_full_old;  // vector of previous step's m
	DoubleMatrix1D mGwork;      // workspace
	DoubleMatrix1D invCp_noise=null;
	DoubleMatrix1D muPrior_noise=null;
	DoubleMatrix1D tilde_dy,tilde_dy_seismic;
	DoubleMatrix1D delta_M;          // combined updates {delta_m,delta_mn}
	DoubleMatrix1D delta_m;          // updates on Gaussian parameters
	DoubleMatrix1D delta_mn=null;    // updates on Noise    parameters
	DoubleMatrix2D Cinv,Cinv11,Cinv12=null,Cinv21=null,Cinv22=null;
	ModifiedCholeskyLeftDecomposition cholDec=null;

	double factor;
	double[] err;

	// Use of finite difference gradient vector
	boolean computeFDgrad = true;
	DoubleMatrix1D gradF=null;

	// Do. X.X^T blockwise. Big holes correspond to missing log data for large windows
	boolean USE_BLOCK_XXT_PRODUCT = true;

 	// here we have found a minimum, and seek the local sensitivity matrix

	// assume noise is decent until proven otherwise
	problemHasDegenerateNoise = false;

	// Number of Gaussian parameters
	nG = modelVariables.nGaussian;
	// Number of noise parameters
	nN = modelVariables.nParams - modelVariables.nGaussian;

	// maximum number of distinct data measurements
	nSeisMeasurements = _WE.m.nWells * _WE.m.nStacks * _WE.convModel.nZ;
	nMeasurements = nSeisMeasurements;
	if (absorbOtherConstraintsInLikelihood)
	  {
		nMeasurements += modelVariables.dimConstraints;
	  }


	X = new DenseDoubleMatrix2D(nMeasurements,nG);
	// view of the design matrix proper to the seismic only
	Xseismic = X.viewPart(0,0,nSeisMeasurements,nG);
	if (USE_BLOCK_XXT_PRODUCT)
	  {
		blockXstarts = new int[_WE.m.nWells * _WE.m.nStacks + 1];
		blockXsizes = blockXstarts.clone();
	  }
	tilde_dy = new DenseDoubleMatrix1D(nMeasurements);
	tilde_dy_seismic = tilde_dy.viewPart(0,nSeisMeasurements);	
	C_prefactor = new DenseDoubleMatrix2D(nG,nG);
	modelVariables.CtildeInvDiag = new double[nG+nN];
	if (nN>0) 
	  {
		// Q^transpose DY involves only seismic mismatches/measurements
		// this is hat{Q} in the writeup
		QtDY = new DenseDoubleMatrix2D(nN, nSeisMeasurements);
	  }
	m_full = new DenseDoubleMatrix1D(nG+nN);
	// m_full_old = m_full.like();
	mG = m_full.viewPart(0,nG);
	if (nN>0) mN = m_full.viewPart(nG,nN);
	mGwork = new DenseDoubleMatrix1D(nG);
	// update space
	delta_M = new DenseDoubleMatrix1D(nG+nN);
	delta_m = delta_M.viewPart(0,nG);
	if (nN>0) delta_mn= delta_M.viewPart(nG,nN);
	// initialise Cpn to be diagonal
	if (nN>0) 
	  {
		invCp_noise = new DenseDoubleMatrix1D(nN).assign(PRIOR_INVERSE_VARIANCE_FOR_LOG_NOISE);
		muPrior_noise = new DenseDoubleMatrix1D(nN).assign(Math.log(modelVariables.initialGuess[nG]));
	  }

	Cinv = new DenseDoubleMatrix2D(nG+nN, nG+nN);
	// selection views on the blocks of Cinv for convenience of manipulation
	Cinv11 = Cinv.viewPart(0,0,nG,nG);
	if (nN>0) 
	  {
		Cinv21 = Cinv.viewPart(nG,0,nN,nG);
		Cinv12 = Cinv.viewPart(0,nG,nG,nN);
		Cinv22 = Cinv.viewPart(nG,nG,nN,nN);
	  }
	// and the RHS
	rhs = new DenseDoubleMatrix1D(nG+nN);
	rhs1 = rhs.viewPart(0,nG);
	if (nN>0) rhs2 = rhs.viewPart(nG,nN);

	// gradients
	if (computeFDgrad)
	  {
		gradF = new DenseDoubleMatrix1D(nG+nN);
	  }

	// ensure in-situ state is set....
	mapParameters(fullParameters_best,true); 
	// set up raw parameters vector
	// & typical scale for fwd differences
	untranslatedParams = new double[modelVariables.nParams+1];
	newParams = untranslatedParams.clone();
	typScale = modelVariables.typicalScale.clone();

	/*
	  if (modelVariables.nParams>11) 
	  {
		((Property)modelVariables.mapFromIndex.get(11)).setValue(0.2);
		_WE.EH.println("REMOVE THIS CHUNK!!!!!!!!");
	  }
	*/

	for (i=0; i<modelVariables.nParams; i++) 
	  {    
		untranslatedParams[i+1] = ((Property)modelVariables.mapFromIndex.get(i)).value;

		// Set local Gaussian model vector
		if (i<nG) 
		  {
			m_full.setQuick(i,untranslatedParams[i+1]);

			typScale[i] = modelVariables.typicalScale[i];
		  }
		else
		  {// a noise parameter with log transformation
			val = Math.log(Math.abs(untranslatedParams[i+1]));
			m_full.setQuick(i,val);
			untranslatedParams[i+1] = val;
			// we expect log noise to get narrowed on log scale by approx sqrt(nMeasurements)...
			typScale[i] = Math.log(modelVariables.typicalScale[i])/Math.sqrt(Math.max(1.0,nMeasurements));
		  }
	  }


	if (loadFaultyNewtonParams && !onePassOnly)
	  {
		// Debugging......
		_WE.EH.println("FIXME REMOVE THIS CHUNK");
		readDebugParametersState("debugparams.bin",untranslatedParams);
		for (i=0; i<modelVariables.nParams; i++) 
		  {
			if (i<nG) 
			  {
				// fullParameters_best[i+1] = untranslatedParams[i+1];
				m_full.setQuick(i,untranslatedParams[i+1]);
			  }
			else
			  {
				m_full.setQuick(i,untranslatedParams[i+1]);
				//y fullParameters_best[i+1] = Math.exp(Math.abs(untranslatedParams[i+1]));
			  }
		  }
		// mapParameters(fullParameters_best,true); 
	  }



	if (nN>0) 
	  {
		tilde_dy_sq_sums = new double[nN];
	  }
	modelVariables.CtildeDiag = new double[nG+nN];

	int maxNewtonSteps = MAX_NEWTON_STEPS;

	modelVariables.useNearestPeakArrivalConstraint = false;

	for (i_trial=0; i_trial<MAX_NUM_FAULTY_NEWTON_REFINEMENTS; i_trial++)
	  {
		modelVariables.modeHasPeakAmbiguity = false;

		for (int NewtonIts=0; NewtonIts<maxNewtonSteps; NewtonIts++)
		  {
			// (A) Get local sensitivity matrix X and delta_y by finite diffs.
			// call at expansion points m_0
			f_0 = val = f_to_minimize_raw(untranslatedParams);
			if (NewtonIts==0) fNewNewtonOld = f_0;

			// feedback
			if ((_WE.verbosity > 3) && (NewtonIts==0) && !runQuietly)
			  {
				_WE.EH.println("Pre-update point ");
				for (i=1;i<=nG;i++) _WE.EH.printf("%10.8g ",untranslatedParams[i]);
				// and the noise levels...
				for (j=0;j<nN;j++,i++) _WE.EH.printf("%10.8g ",untranslatedParams[i]);
				_WE.EH.println("");
				if (nN==0) 
				  {
					_WE.EH.printf("Noise levels fixed at (stacks(s)): ");
					for (j=0;j<_WE.m.nStacks;j++) _WE.EH.printf("%10.8g ",_WE.m.seismic_noise[j].value);
					_WE.EH.println("");
				  }
				/*
				  _WE.EH.println("%nPrior Std deviations....");
				  for (i=0; i<modelVariables.C_prior.rows(); i++) 
				  {
				  _WE.EH.printf("%10.8g ",Math.sqrt(modelVariables.C_prior.getQuick(i,i)));
				  }
				*/
				if (_WE.verbosity==4)
				  {
					_WE.EH.println("f_update "+NewtonIts+" = "+val);
				  }
				else
				  {
					_WE.EH.println("f_update Newton_its= "+NewtonIts+" func_eval_count= "+it_count+" f = "+val);
				  }
			  }
			baseErrors = ArrayUtils.deepCopy(_WE.convModel.errors);

			// Util.printfileZO("seis.txt",_WE.convModel.S_true[0][0],_WE.convModel.S_true[0][0].length); // FIXME DELETE 

			// Form residual vector dy = (y-f(m,mn))/noise
			for (j=0,row=0; j<_WE.m.nWells; j++) {
			  for (k=0; k<_WE.m.nStacks; k++) {
				for (l=0; l<_WE.convModel.nZ; l++,row++) {
				  tilde_dy.setQuick(row, -baseErrors[j][k][l]);
				}
			  }
			  // Util.debugPrintMatrixToFile("tilde_dy."+NewtonIts+".Acc-well."+j+".txt",tilde_dy); // FIXME DELETE
			  // Util.printfileZO("Ssynth."+NewtonIts+".well."+j+".txt",_WE.convModel.S_synth[j][0],_WE.convModel.nZ); // FIXME DELETE
			  // Util.printfileZO("Reff."+NewtonIts+".well."+j+".txt",_WE.convModel.Reff[j][0],_WE.convModel.nZ); // FIXME DELETE
			}


			// Plus any non-synthetic constraints
			if (absorbOtherConstraintsInLikelihood)
			  {
				baseAuxiliaryErrors = modelVariables.constraintErrors.clone();
				for (j=0; j<modelVariables.dimConstraints; j++,row++) 
				  {
					tilde_dy.setQuick(row,-modelVariables.constraintErrors[j]); // here tilde_dy = (y_aux - f_aux(m))/sigma 
				  }
			  }
			// Util.debugPrintMatrixToFile("tilde_dy."+NewtonIts+".txt",tilde_dy);
			// Util.debugPrintMatrixToFile("tilde_dy_seismic."+NewtonIts+".txt",tilde_dy_seismic);

			// get sensitivity matrix X and delta_y, for Gaussian parameters
			for (i=1,column=0; i<=nG; i++,column++)
			  {
				val = untranslatedParams[i]; // oldvalue
				// epsScale =  Math.max(typScale[i-1],Math.abs(val)); // suitable scale of perturbation
				epsScale =  typScale[i-1];                        // fixes issues with large offsets......
				denom = DESIGN_MATRIX_EPS * epsScale;             //  = newvalue - oldvalue
				untranslatedParams[i] = val + denom;              //  = newvalue

				// call function with perturbation
				f_prime = f_to_minimize_raw(untranslatedParams);

				// populate design matrix
				// next three loops are  over X's rows, or observations, for a fixed column
				// term for Finite diff denominator & standard error removal
				factor = 1.0/denom;
				row = 0;
				for (j=0,k_b=0; j<_WE.m.nWells; j++)
				  {
					for (k=0; k<_WE.m.nStacks; k++)
					  {
						err = _WE.convModel.errors[j][k]; 
						if (USE_BLOCK_XXT_PRODUCT)
						  {
							row+=_WE.convModel.kStarts[j];
							blockXstarts[k_b] = row;
							for (l=_WE.convModel.kStarts[j]; l<=_WE.convModel.kEnds[j]; l++,row++)
							  {
								delta = factor*(err[l] - baseErrors[j][k][l]);
								X.setQuick(row,column,delta);
							  }
							blockXsizes[k_b] = row - blockXstarts[k_b];
							k_b++;
							row += (_WE.convModel.nZ -1 - _WE.convModel.kEnds[j]);
						  }
						else
						  {
							for (l=0; l<_WE.convModel.nZ; l++,row++)
							  {
								delta = factor*(err[l] - baseErrors[j][k][l]);
								X.setQuick(row,column,delta);
							  }	
						  }
					  }
				  }
				// absorb remaining "observations" in next block
				if (USE_BLOCK_XXT_PRODUCT) blockXstarts[k_b]=row;

				// Plus any non-synthetic constraints
				if (absorbOtherConstraintsInLikelihood)
				  {
					// populate design matrix
					// loops over X's rows, or observations, for a fixed column
					err = modelVariables.constraintErrors; 
					for (j=0; j<modelVariables.dimConstraints; j++,row++) 
					  {
						delta = factor*(err[j] - baseAuxiliaryErrors[j]);
						X.setQuick(row,column,delta);
					  }
				  }
				// end of block
				if (USE_BLOCK_XXT_PRODUCT) blockXsizes[k_b] = row - blockXstarts[k_b];

				if (computeFDgrad)
				  {
					// store -grad for convenience
					gradF.setQuick(column,-factor*(f_prime - f_0));
				  }
				//Util.print(untranslatedParams,untranslatedParams.length);
				// restore parameter i at base value
				untranslatedParams[i] = val;
			  }

			if (computeFDgrad)
			  {
				// Continue on using same FD's for noise derivatives, if required
				// get sensitivity matrix X and delta_y, for Gaussian parameters
				for (; i<=nG+nN; i++,column++)
				  {
					val = untranslatedParams[i];                      // oldvalue
					epsScale = typScale[i-1];                         // NOT!!! Math.max(typScale[i-1],Math.abs(val)); // suitable scale of perturbation
					denom = DESIGN_MATRIX_EPS * epsScale;             //  = newvalue - oldvalue
					untranslatedParams[i] = val + denom;              //  = newvalue

					// call function with perturbation
					f_prime = f_to_minimize_raw(untranslatedParams);
			
					// Get function FD gradient terms
					factor = 1.0/denom;
					// store -grad for convenience
					gradF.setQuick(column,-factor*(f_prime - f_0));

					// restore parameter i at base value
					untranslatedParams[i] = val;
				  }
			  }


			// Util.debugPrintMatrixToFile("X."+NewtonIts+".txt",X); // FIXME remove
			// Util.debugPrintMatrixToFile("Xseismic."+NewtonIts+".txt",Xseismic); // FIXME remove
			if (false)
			  {
				// debugging test for pure linear model
				DoubleMatrix1D testRHS = tilde_dy.like();
				blas.dgemv(false,1.0,X,mG,0.0,testRHS);
				// add update to existing point
				testRHS.assign(tilde_dy, Functions.plus);
				Util.debugPrintMatrixToFile("testVec."+NewtonIts+".txt",testRHS); // FIXME remove
				Util.debugPrintMatrixToFile("X."+NewtonIts+".txt",X); // FIXME remove
			  }
			// OK now we have the design matrix X for the Gaussian parameters and
			// the tilde_dy vector

			// (B) Form RHS of equation for Gaussian model parameters update
		
			if (nN>0) 
			  {
				// Current log(noise) vector is in mN:

				// Set up QtDY
				// No contributions are needed here from auxiliary constraints
				for (j=0,col=0; j<_WE.m.nWells; j++)
				  {
					for (k=0; k<nN; k++)
					  {
						for (l=0; l<_WE.convModel.nZ; l++,col++)
						  {
							QtDY.setQuick(k,col,-2*tilde_dy.getQuick(col));
						  }
					  }
				  }
				// Util.debugPrintMatrixToFile("QtDY.txt",QtDY);		  

				// Update Cinv22
				// Form sums over (tilde_dy)^2 for each noise component (stack)
				// I don't think any contributions are needed here from auxiliary constraints
				Arrays.fill(tilde_dy_sq_sums, 0.0);
				for (j=0,row=0; j<_WE.m.nWells; j++)
				  {
					for (k=0; k<nN; k++)
					  {
						for (l=0; l<_WE.convModel.nZ; l++,row++)
						  {
							tilde_dy_sq_sums[k] += tilde_dy.getQuick(row) * tilde_dy.getQuick(row); 
						  }
					  }
				  }

				// clear accumulators
				if ((NewtonIts>0) || (i_trial>0)) Cinv22.assign(0.0);
				
				// Write this Phat information straight into Cinv22
				for (i=0;i<nN;i++) Cinv22.setQuick(i, i, 2.0*tilde_dy_sq_sums[i]);
				// At this point Cinv22 = Phat

				// rhs2 is initialised as -e_st - almost certainly the same for all stacks i
				for (i=0;i<nN;i++) 
				  {
					rhs2.setQuick(i, _WE.ignoreSeismic ? 0 : -_WE.convModel.dimData[i]);
				  }
				// Form rhs2 -= (1/2) QtDY.tilde_dy 
				blas.dgemv(false,-0.5,QtDY,tilde_dy_seismic,1.0,rhs2);
				if (true)
				  {
					blas.dgemv(false,-0.5,QtDY,tilde_dy_seismic,0.0,rhs2);
					for (i=0;i<nN;i++) 
					  {
						rhs2.setQuick(i, rhs2.getQuick(i)+( _WE.ignoreSeismic ? 0 : -_WE.convModel.dimData[i]));
					  }
				  }



				// And then rhs2 -= Cpn^{-1}.(m_0^{(n)}-mu^{(n)})
				for (i=0;i<nN;i++) 
				  {
					rhs2.setQuick(i,rhs2.getQuick(i) - 
								  (mN.getQuick(i)-muPrior_noise.getQuick(i))*invCp_noise.getQuick(i));
				  } 
			  }
			// form X.tilde_dy term
			// Util.debugPrintMatrixToFile("tilde_dy."+NewtonIts+".txt",tilde_dy); // FIXME DELETE

			blas.dgemv(true,1.0,X,tilde_dy,0.0,rhs1);
			// first make sure mGwork is set for -ve the Gaussian parameters
			for (i=0;i<nG;i++) mGwork.setQuick(i,-untranslatedParams[i+1]);
			// add mu_prior
			mGwork.assign(modelVariables.mu_prior, Functions.plus);
			// then mult by CP^{-1}
			if (true)
			  {
				blas.dgemv(false,1, modelVariables.Cinv_prior, mGwork,1, rhs1);
			  }
			else
			  {
				DoubleMatrix1D tmpRhs = mGwork.like();
				blas.dgemv(false,1, modelVariables.Cinv_prior, mGwork,0,tmpRhs);
				rhs1.assign(tmpRhs, Functions.plus);
			  }
			// .........OK now the RHS for the delta_m equation is in rhs={rhs1,rhs2}....................

			// Util.debugPrintMatrixToFile("rhs."+NewtonIts+".txt",rhs);
			// Util.debugPrintMatrixToFile("gradF."+NewtonIts+".txt",gradF);
			// Util.debugPrintMatrixToFile("X."+NewtonIts+".txt",X);

			if (computeFDgrad)
			  {
				// use numerical FD rhs:
				rhs.assign(gradF);
			  }

			// Form the Coefficient matrix
			// (1) Form C_prefactor = X^T.sigma.X + Cprior^-1
			// (the sigma bit is already absorbed into X)
			// This is the inverse covariance of the updated model
			C_prefactor.assign(modelVariables.Cinv_prior);
			if (USE_BLOCK_XXT_PRODUCT)
			  {
				// X^T.X = sum_blocks_i X_i^T.X_i
				for (i=0; i<blockXstarts.length; i++)
				  {
					DoubleMatrix2D Xblock = X.viewPart(blockXstarts[i],0,blockXsizes[i],nG);
					blas.dgemm(true,false,1.0,Xblock,Xblock,1.0,C_prefactor);					  
				  }
			  }
			else
			  {
				blas.dgemm(true,false,1.0,X,X,1.0,C_prefactor);
			  }

			// Take copy of this while it's available
			Cinv11.assign(C_prefactor);

			if (nN>0)
			  {
				// Form -Q^T DY X piece
				QtDY_X = QtDY.zMult(Xseismic,null,-1,0,false,false);
				// Take copy of this while it's there
				Cinv21.assign(QtDY_X);
				// copy to Cinv12
				Cinv12.assign(Cinv21.viewDice());
				// form Cinv22 inverse and add the Noise prior while we're here
				for (i=0;i<nN;i++)
				  {
					Cinv22.setQuick(i, i, Cinv22.getQuick(i,i) + invCp_noise.getQuick(i));
				  }
				// Cinv22 is now actually (C_pn^{-1}+Phat)
			  }
			// Store inverse diagonal elements
			for (i=0;i<nG+nN;i++) modelVariables.CtildeInvDiag[i] = Cinv.getQuick(i,i);

			// ..........Ok now the coefficient matrix Cinv is in place ...............
			// Util.debugPrintMatrixToFile("Cinv."+NewtonIts+".txt",Cinv);

			// (C) now solve the Cinv.{delta_m,delta_mn} = rhs system
			// Prescale to look after wierd scalings that can occur (ft, secs)
			cholDec  = new ModifiedCholeskyLeftDecomposition(Cinv,true,false);
			cholDec.solve(delta_M,rhs);


			// Add update to existing point. 

			// Perform step-size limitation. Scale wrt typical modelVariables.typicalScale,
			// since delta_m=0 is a possibility
			for (i=0,stp_size=0.0; i<nG; i++) 
			  {
				stp_size += MathUtil.sqr(delta_m.getQuick(i)/typScale[i]);
			  }
			if (nN>0)
			  {
				// typScale = log(modelVariables.typicalScale[i]) for noise params
				for (; i<nG+nN; i++) 
				  {
					stp_size += MathUtil.sqr(delta_m.getQuick(i)/((typScale[i]!=0) ? typScale[i] : 1.0));
				  }
			  }
			stp_size = Math.sqrt(stp_size);
			if (stp_size > MAX_STEP_SIZE)
			  {
				// reduce step. scale all components by length, so overall step is MAX_STEP_SIZE
				delta_m.assign(cern.jet.math.Mult.mult(MAX_STEP_SIZE/stp_size));
				if (nN>0) delta_mn.assign(Mult.mult(MAX_STEP_SIZE/stp_size));
			  }

			mG.assign(delta_m, Functions.plus);
			if (nN>0) mN.assign(delta_mn, Functions.plus);


			// embed & de-transform these into newParams
			for (i=0;i<nG;i++) newParams[i+1] = mG.getQuick(i);
			// update noise...
			if (nN>0) 
			  {
				for (j=0; j<nN; j++,i++) 
				  {
					newParams[i+1] = mN.getQuick(j);
					if (Math.exp(newParams[i+1]) < (MINIMUM_TOLERABLE_SN_RATIO*
													_WE.multiModelSetUp.typicalRMSSeismicAmplitude))
					  {
						_WE.EH.warn("Warning: noise underflow - thresholding for safety");
						newParams[i+1] = Math.log(MINIMUM_TOLERABLE_SN_RATIO
												  *_WE.multiModelSetUp.typicalRMSSeismicAmplitude);
						problemHasDegenerateNoise = true;
					  }
				  }
			  }

			if (dumpFaultyNewtonParams)
			  {
				_WE.EH.println("FIXME REMOVE THIS CHUNK");
				writeDebugParametersState("debugparams.bin",untranslatedParams);
			  }


			// do line search along joining line for best point, leaving best fit embedded in newParams...
			getBestPointAlongJoiningLine(untranslatedParams,newParams,f_0,rhs,cholDec);

			

			// To produce no line_search update on pass (testing purposes)
			// set 2-> desired number. With (optimiserPass>2), is harmless (never applies)
			if (optimiserPass>2)
			  {
				System.arraycopy(untranslatedParams,1,newParams,1,untranslatedParams.length-1);
				_WE.EH.warn("COMMENT OUT LINE ABOVE ME!!!!");
			  }

			// embed it
			negLogPost = f_to_minimize_raw(newParams);

			// provide feedback
			if (_WE.verbosity>3 && !runQuietly)
			  {
				_WE.EH.println("Newton update "+NewtonIts);
				if (_WE.verbosity==4) 
				  {
					_WE.EH.println("%nf_update %d = %10.8g%n", NewtonIts, negLogPost);
				  }
				else
				  {
					_WE.EH.printf("f_update Newton_its= %d func_eval_count= %d f= %10.8g%n", NewtonIts, it_count, negLogPost);
				  }
				for (i=1; i<=nG+nN; i++) _WE.EH.printf("%10.8g ",newParams[i]);
				_WE.EH.println("");
			  }

			

			System.arraycopy(newParams,1,untranslatedParams,1,untranslatedParams.length-1);

			// embed parameters vector in current views mG,mN,m_full
			for (i=0;i<nG+nN;i++) m_full.setQuick(i,untranslatedParams[i+1]);

			// Rough termination criterion
			// set Newton variables
			fNewNewton = negLogPost;
			delta_fNewNewton = fNewNewton - fNewNewtonOld;
			fNewNewtonOld = fNewNewton;

			if ( onePassOnly ||
				 // (((NewtonIts>0) && (2.0*Math.abs(delta_fNewNewton) < NEWTON_EPS * (Math.abs(fNewNewton) + EPS_FLOATINGPOINT)))))
				 (((NewtonIts>=0) && (Math.abs(delta_fNewNewton) < ABS_NEWTON_EPS))))
			  {
				// Accept Newton termination only if last line search did not involve a 
				// discontinuous peak-arrival jump, and if we're using strict peak-hunting
				if (_WE.constrainPeakArrival 
					&& (Math.abs(tPeakStartLineSearch - tPeakEndLineSearch) > _WE.dt))
				  {
					// Looks bad.  We're trapped in a discontinuous minima.
					if (HANDLE_AMBIGUOUS_PEAK_DETECTION_TERMINATION_USING_PENALTY)
					  {
						// add penalty to give mode posterior and bail out
						negLogPost += AMBIGUOUS_PEAK_DETECTION_TERMINATION_PENALTY;
						modelVariables.modeHasPeakAmbiguity = true;
						if (_WE.verbosity>3 && !runQuietly) 
						  {
							_WE.EH.println("Newton terminated with discontinuous peak-arrival behaviour, after "+
										   (NewtonIts+1)+" iterations, "+ it_count+ " function evaluations");
						  }
						break;
					  }
					else
					  {
						if (!modelVariables.useNearestPeakArrivalConstraint)
						  {
							// We are using strict (arbitrated) peak detection, and the peak is still ambiguous
							// on the last line search. Flip to "nearest peak" detection
							_WE.EH.println("WARNING: unsatisfactory Newton termination criterion: peak detection discontinuity");  
							// Start Newton afresh
							NewtonIts = -1;
							modelVariables.useNearestPeakArrivalConstraint = true;
						  }
					  }
				  }
				else
				  {				
					if (_WE.verbosity>3 && !runQuietly) 
					  {
						_WE.EH.println("Newton terminated after "+(NewtonIts+1)+" iterations, "+ it_count+ " function evaluations");
					  }
					break;
				  }
			  } 
			else
			  {
				// not converged. If we haven't done this by maxNewtonSteps something's awry
				if ((NewtonIts==(maxNewtonSteps-1)) && (maxNewtonSteps>1) && (_WE.verbosity>3))
				  {
					_WE.EH.println("WARNING: unsatisfactory Newton termination criterion");  
					// Probably the mode's not worth trusting, but continue anyway
				  }
			  }
		  }
		
		

		// make sure line search is embedded
		for (i=0;i<nG+nN;i++) mG.setQuick(i,untranslatedParams[i+1]);

		// store salient quantities for this mode
		// Cholesky L factor 
		modelVariables.L_post = cholDec.getL();

		//Util.debugPrintMatrixToFile("Lpost.txt",modelVariables.L_post); // FIXME DELETE

		// (& we need its inverse remember)
		cholDec.getCholeskyLeftInverse(modelVariables.L_post);
		// Util.debugPrintMatrixToFile("LpostInv.txt",modelVariables.L_post); // FIXME DELETE
		modelVariables.L_Permutation = cholDec.getPermutation(); 
		modelVariables.mu_post = m_full;

		// Compute log(noise) variances while we're here
		for (i=0; i<(nG+nN); i++)
		  {
			for(j=i,sum=0.0; j<(nG+nN); j++) sum += MathUtil.sqr(modelVariables.L_post.getQuick(j,i));
			// insert posterior variance in correct (de-permuted) index in Ctilde diagonal position
			modelVariables.CtildeDiag[modelVariables.L_Permutation[i]] = sum;
		  }

		// Check Newton convergence by collocation
		if (onePassOnly)
		  {
			// this is only a rough guess. Go in blithely
			break;
		  }
		else if (_WE.constrainPeakArrival && !_WE.constrainPeakArrivalSign && modelVariables.modeHasPeakAmbiguity)
		  {
			// Discontinuous solution due to peak ambiguity. Outside this function we will try
			// again, with constrainPeakArrivalSign=true, so we might as well give up now.
			break;
		  }
		else
		  {
			// final point should be embedded by modeHasFaultyNewtonConvergence()...
			if (modeHasFaultyNewtonConvergence(i_trial))
			  {
				if (false)
				  {
					scaleVectorByCurrentEstimatedUncertainty(ls_fullParameters, false);
					scaleVectorByCurrentEstimatedUncertainty(gradF, true);

					double norm = SeqBlas.seqBlas.ddot(gradF,ls_fullParameters);
					double gradNorm = SeqBlas.seqBlas.dnrm2(gradF); 
					double sdNorm = SeqBlas.seqBlas.dnrm2(ls_fullParameters); 
					_WE.EH.println("FIXME COMMENT OUT : Projection of gradient in "
								   +"suggested FAULTY-Newton direction  = "+norm/(gradNorm*sdNorm));
				  }
				// We will have to keep iterating...
				if (_WE.verbosity>3) _WE.EH.println("FAULTY NEWTON.. "+i_trial);
				
				if (dumpFaultyNewtonParams)
				  {
					System.exit(0);
				  }

				// embed improved point coming from modeHasFaultyNewtonConvergence() in untranslatedParams
				for (i=0;i<nG+nN;i++) untranslatedParams[i+1] = modelVariables.mu_post.getQuick(i);
			  }
			else
			  {
				// if (_WE.verbosity>3) _WE.EH.println("NEWTON NOT FAULTY.. "+i_trial);
				// looks OK. Break out out Newton
				break;
			  }
		  }
	  }


	// Update other terms required for sampling
	// 1/2(Determinant of C_eff)
	// since C_eff = P.L_post^T.L_post.P^T and |P|=1
	double halfLogDetCeff = 0.0;

	for (i=0; i<(nG+nN); i++)
	  {
		halfLogDetCeff += Math.log(modelVariables.L_post.getQuick(i,i));
	  }

	modelVariables.halfLogDetCeff = halfLogDetCeff;  // FIXME don't really need this variable

	// Mode Laplace factor
	modelVariables.logModeWeight = halfLogDetCeff - negLogPost;

	// Here we could use the approximate Bayes Information Criterion BIC = 0.5*n_p*log(n_d)
	// approach. But this depends only on the likelihood...  
	// But since BIC is only asymptotic, the Laplace estimator already acquired is probably better......
	// The number of parameters differs only in the wavelet coefficients, 
	// so that's all we'll include here
	int nParams = _WE.m.wavelet.valuesCoarseProp.length;
	if (_WE.useBIC && ENABLE_BIC_USE)
	  {
		if (_WE.ignoreSeismic)
		  {
			// if we ignore Seismic, it is as if we make no observations,and thus 
			// there should be no BIC (all models are equiprobable)
			modelVariables.logModeWeight = 0.0;
		  }
		else
		  {
			double nBIC=0;
			for (i=0;i<_WE.m.nStacks;i++) 
			  {
				// Since n_eff/N = dt/dtCoarse and
				// dimData = sqrt(n_eff/N) * N,  pull in the extra sqrt factor
				nBIC += _WE.convModel.dimData[i] * Math.sqrt(_WE.dt/_WE.m.wavelet.dtCoarse);
			  }
			// The standard BIC rule	
			modelVariables.logModeWeight = -computeFastSyntheticNegLogLikelihood(true)
			  - 0.5*nParams*Math.log(nBIC);
			_WE.EH.println("BIC: "+modelVariables.logModeWeight);
		  }
	  }
	// Get noise covariance
	modelVariables.noiseCorrCoeffMatrix = _WE.convModel.computeMultiStackNoiseCovarianceMatrix();


	if (_WE.verbosity>2 && !runQuietly)
	  {
		// FIXME this hack applies only for first stack......

		_WE.EH.println("(1)MODE (2)No.-params (3)Laplace-log(weight) (4)-log(posterior) (5)SSE(stack1) "+
					   "(6)-log(likelihood) (7)-log(prior) (8)noise(stack-1)-rms (9)halfLogDetCeff (10)halflogDetCp "+
					   "(11)halflogDetCp-halfLogDetCeff");
		_WE.EH.println("MODE: %3d %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g",nParams,
					   modelVariables.logModeWeight,negLogPost,
					   _WE.convModel.SSerrors[0],
					   computeFastSyntheticNegLogLikelihood(true),
					   modelVariables.computePriorNegLogLikelihood(false),
					   _WE.m.seismic_noise[0].value,
					   modelVariables.halfLogDetCeff,
					   modelVariables.logDetC_prior_only/2.0,
					   modelVariables.logDetC_prior_only/2.0 -modelVariables.halfLogDetCeff);
	  }

	return negLogPost;
  }


  private void writeDebugParametersState(String filename, double[] data)
  {
	try 
	  {
		// Create an output stream to the file.
		FileOutputStream file_output = new FileOutputStream(filename);
		// Wrap the FileOutputStream with a DataOutputStream
		DataOutputStream data_out = new DataOutputStream (file_output);

		// Write the data to the file in an integer/double pair
		for (int i=0; i < data.length; i++) {
		  data_out.writeDouble(data[i]);
		}
		// Close file when finished with it..
		file_output.close();
	  }
    catch (IOException e) 
	  {
		System.out.println ("IO exception = " + e );
	  }
	System.out.println ("Wrote debug.bin");
  }

  private void readDebugParametersState(String filename, double[] data)
  {
	try 
	  {
		// Wrap the FileInputStream with a DataInputStream
		FileInputStream file_input = new FileInputStream(filename);
		DataInputStream data_in    = new DataInputStream(file_input);
		
		int i=0;
		while (i<data.length) 
		  {
			try 
			  {
				data[i++] = data_in.readDouble();
			  }
			catch (EOFException eof) 
			  {
				System.out.println ("End of File");
				break;
			  }
		  }
		data_in.close ();
	  } 
	catch  (IOException e) 
	  {
		System.out.println ( "IO Exception =: " + e );
	  }
  }



  private void scaleVectorByCurrentEstimatedUncertainty(DoubleMatrix1D params, boolean isGradient)
  {
	int i;
	double scale;
	for (i=0; i<params.size(); i++)
	  {
		scale = 1.0/Math.sqrt(modelVariables.CtildeInvDiag[i]);
		if (isGradient) params.setQuick(i,params.getQuick(i) * scale);
		else params.setQuick(i, params.getQuick(i)/scale);
	  }
  }

  /**
	 Line search algorithm (rough and ready...)
	 Note f77ish unit offset convention.
	 The vector gradF is the gradient vector in the *descent* direction at xOld.
	 f0 is the function value (already known) at xOld.
  */
  private void getBestPointAlongJoiningLine(double[] xOld, double[] xNew, double f0, DoubleMatrix1D gradF,
											ModifiedCholeskyLeftDecomposition cholDec)
  {
	int i,j,n,i_restarts,pass=0,firstPass;
	double frac,FRAC_DECREMENT=0.1; // FIXME PROMOTE
	double negLogPost=f0, negLogPostBest=Double.MAX_VALUE;
	double[] xTry = xOld.clone();
	double[] xBest = xOld.clone();
	double fracBest=1.0;

	boolean searchImprovedOnThisPass = false;
	boolean usingSD = false;

	n = xTry.length;
	double fracStart = 1.0;
	double fracEnd = 0.0;
	double fracDecrement = FRAC_DECREMENT;

	if (true)
	  {
		// Multi-pass line search by collocation. More expensive than cubics etc, but
		// leads to better solutions and then fewer Newton iterations in outer loop:
		// cost is negligible compared to FD gradient evaluation.
		for (pass=0,firstPass=0; pass<=1; pass++,firstPass++)
		  {
			for (frac=fracStart; frac>=fracEnd; frac-=fracDecrement)
			  {
				for (i=1;i<n;i++)
				  {
					xTry[i] = frac*xNew[i]+(1.0-frac)*xOld[i];
				  }
				negLogPost = f_to_minimize_raw(xTry);

				if (_WE.constrainPeakArrival && !_WE.ignoreSeismic && firstPass==0)
				  {			
					// remember location of main peak for start and end of coarsest-scale line search
					// wavelet.peakTime is set by the peak search buried within f_to_minimize_raw()
					if (frac==fracStart)
					  {
						// Newton step
						tPeakEndLineSearch = modelVariables.wavelet.peakTime;
					  }
					else if (frac <= (fracEnd+0.5*fracDecrement))
					  {
						tPeakStartLineSearch = modelVariables.wavelet.peakTime;
						// _WE.EH.println("FIXME PEAK GAP %9.6g %9.6g",tPeakStartLineSearch, tPeakEndLineSearch);
					  }	
				  }

				if (_WE.verbosity>3)
				  {
					if (frac<=1.0)
					  {
						// _WE.EH.println((usingSD?"Steepest desc. ":"")+ "Line-search target point: f_post="+negLogPost);
						_WE.EH.printf("Mapped variables (frac=%4.3f): ",frac);
						for (i=1; i<n; i++) 
						  {
							_WE.EH.printf("%9.6g ",((Property)modelVariables.mapFromIndex.get(i-1)).value);
						  }
						_WE.EH.println("");
					  }
					_WE.EH.println("Line-search %2d: %4.3f %10.8g %d",pass,frac, negLogPost, (it_count-1));
					/*
					  _WE.EH.println("Line-search FIXME %2d: %4.3f %10.8g %d %10.8g",
					  pass,negLogPost, (it_count-1),
					  _WE.m.wavelet.getPeakTime(false,modelVariables.constrainPeakArrivalSign,modelVariables.peakArrivalPositive));
					*/
				  }
				if (negLogPost < negLogPostBest)
				  {
					negLogPostBest = negLogPost;
					fracBest = frac;
					System.arraycopy(xTry,1,xBest,1,n-1);
					// if some improvement....
					if (Math.abs(frac)>1e-5) 
					  {
						searchImprovedOnThisPass = true;
					  }
					else
					  {
						searchImprovedOnThisPass = false;
					  }
				  }
			  }
			if (usingSD) 
			  {
				// Last attack will be SD, if used. Can break here.
				break;
			  }

			if (pass==0 && searchImprovedOnThisPass && Math.abs(negLogPostBest - negLogPost) < EPS_SUSPICIOUS_NEWTON_STEP 
				&& fracBest <= 4*FRAC_SUSPICIOUS_NEWTON_STEP)
			  {
				// Newton "kind of" worked, but not very well on [0,1] interval. 
				// we got a narrow improvement, and nowhere near the full Newton step
				// Try steepest descent  for luck (as well)...
				searchImprovedOnThisPass = false;
				usingSD = true;
				if (_WE.verbosity>3)
				  {
					_WE.EH.println("Line-search: suspiciously bad Newton step - trying steepest descent (as well)...");
				  }
			  }

			// _WE.EH.println("Line-search end pass: negLogPostBest = "+negLogPostBest+", fracBest = "+fracBest);
			if (searchImprovedOnThisPass) 
			  {
				// break on improvement
				break;
			  }
			else
			  {
				if (!usingSD && fracStart >= 0.5)
				  {
					// Finer Newton search: reduce forward search by a decade
					fracStart /=10.0;
					fracEnd = 0;
					fracDecrement/=10.0;
					pass--;
				  }
				else
				  {
					if (_WE.verbosity>3)
					  {
						_WE.EH.println("Line-search: switching to steepest descent:");
					  }
					pass--;
					usingSD = true;

					// switch to steepest descent
					// (1) Bracketing phase
					double alpha_0 = 0.0;
					double deltaFtarget = 1.0;
					DoubleMatrix1D gradDirection = gradF;

					if (false)
					  {
						// Gradient check
						gradDirection = gradF.copy();
						gradF.assign(0);
						int column;
						double f_0 = f_to_minimize_raw(xOld);
						double val,denom,factor,f_prime,epsScale;
						double maxEps=0,minEps=1e50;
						for (i=1,column=0; i<=modelVariables.nGaussian; i++,column++)
						  {
							val = xOld[i]; // oldvalue
							// epsScale = Math.max(typScale[i-1],Math.abs(val)); // suitable scale of perturbation
							epsScale = typScale[i-1]; // suitable scale of perturbation
							denom = DESIGN_MATRIX_EPS * epsScale;             //  = newvalue - oldvalue

							double scale = 1.0/Math.sqrt(modelVariables.CtildeInvDiag[column]);
							minEps=Math.min(denom/scale, minEps);
							maxEps=Math.max(denom/scale, maxEps);

							xOld[i] = val + denom;              //  = newvalue

							// call function with perturbation
							f_prime = f_to_minimize_raw(xOld);

							factor = 1.0/denom;
							// store -grad for convenience
							gradF.setQuick(column,-factor*(f_prime - f0));
							// restore parameter i at base value

							xOld[i] = val;
						  }
						// Continue on using same FD's for noise derivatives, if required
						// get sensitivity matrix X and delta_y, for Gaussian parameters
						for (; i<=modelVariables.nParams ; i++,column++)
						  {
							val = xOld[i];                      // oldvalue
							// epsScale = Math.max(typScale[i-1],Math.abs(val)); // suitable scale of perturbation
							epsScale = typScale[i-1]; // suitable scale of perturbation
							denom = DESIGN_MATRIX_EPS * epsScale;             //  = newvalue - oldvalue
							double scale = 1.0/Math.sqrt(modelVariables.CtildeInvDiag[column]);
							minEps=Math.min(denom/scale, minEps);
							maxEps=Math.max(denom/scale, maxEps);

							xOld[i] = val + denom;              //  = newvalue
						  
							// call function with perturbation
							f_prime = f_to_minimize_raw(xOld);
						  
							// Get function FD gradient terms
							factor = 1.0/denom;
							// store -grad for convenience
							gradF.setQuick(column,-factor*(f_prime - f_0));
						  
							// restore parameter i at base value
							xOld[i] = val;
						  }
						// _WE.EH.println("EPS SCALE MAX/MIN = "+minEps+" "+maxEps);
					  }
					// Scale transformations
					if (false)
					  {
						_WE.EH.println("GRAD RESCALED...CHECK???");
						// Scale transformations
						// Rescale gradient to direction implied by natural scales
						// implied by the diagonal of C^{-1}
						gradDirection = gradF.like();
						for (i=0; i<n-1; i++)
						  {
							double s = Math.sqrt(Math.max(0.0,modelVariables.CtildeInvDiag[i]));
							s = (s>0)? 1.0/s : 1.0;
							gradDirection.setQuick(i,s*gradF.getQuick(i));
						  }
						// Make gradDirection a unit vector
						double norm = SeqBlas.seqBlas.dnrm2(gradDirection); 
						gradDirection.assign(Mult.div(norm));
						// gradDirection is now unit vector in its original direction
						norm = SeqBlas.seqBlas.ddot(gradF,gradDirection);
						// scale by projection of gradF in gradDirection direction
						gradDirection.assign(Mult.mult(norm));
					  }
					if (false)
					  {
						// project gradient in direction of current largest eigendirection
						// at this point remember L is cholesky factor of inverse
						DoubleMatrix1D eDir = new DenseDoubleMatrix1D(getLargestEigenDirection(cholDec));
						double norm = SeqBlas.seqBlas.ddot(eDir,gradF);
						_WE.EH.println("SD-direction.slowest-eVector = "+ norm);
						gradDirection = eDir.copy();
						gradDirection.assign(Mult.mult(norm));
					  }
					for (i=0;i<(n-1);i++)
					  {
						alpha_0 = Math.max(alpha_0,Math.abs(gradDirection.getQuick(i)));
					  }
					if (alpha_0<=0) break; // seems pointless to continue.

					// Protection in case of discontinuities giving enormous gradients...
					double alpha_0_min=Double.MAX_VALUE,tmp,EPS_STEP_SIZE=0.01;
					for (i=0;i<(n-1);i++)
					  {
						tmp = Math.abs(gradDirection.getQuick(i));
						if (tmp>0) alpha_0_min = Math.min(alpha_0_min, modelVariables.typicalScale[i]/tmp);
					  }
					alpha_0_min *= EPS_STEP_SIZE;


					// start at 1/10th of gradient-inspired step size:
					alpha_0 = 0.1*deltaFtarget/(alpha_0*alpha_0);
					// but no smaller than natural length scales would imply as reasonable
					alpha_0 = Math.max(alpha_0, alpha_0_min);
					// no more than 20 doublings (safe upper bound)
					// if we need to restart, no more than 5 times (this is alpha_0 *= 1e-5)
					for (j=0,i_restarts=0; (j<=20) && (i_restarts<=5); alpha_0 *= 2.0, j++)
					  {
						for (i=1;i<n;i++)
						  {
							xNew[i] = xOld[i] + alpha_0*gradDirection.getQuick(i-1);
						  }
						negLogPost = f_to_minimize_raw(xNew);
						// _WE.EH.println("GRAD SCALING j="+j+" f="+negLogPost);
						if ((j==0) && (negLogPost > f0)) 
						  {
							// first step was too aggressive. Back off 
							// alpha by 0.1 (0.05 anticipates alpha_0*=2.0)
							alpha_0 *= 0.05;
							j=-1;
							i_restarts++;
							continue;
						  }
						if (negLogPost > f0) 
						  {
							// interval is now bracketed
							// _WE.EH.println("Line-search scaling: new-current="+(negLogPost-negLogPostBest));
							break;
						  }
					  }
					// Simple test of nearest point
					fracStart = 1;
					fracEnd = 0;
					fracDecrement=.1;
					for (i=1;i<n;i++)
					  {
						xNew[i] = xOld[i] + fracDecrement*alpha_0*gradDirection.getQuick(i-1);
					  }
					negLogPost = f_to_minimize_raw(xNew);
					if (negLogPost >= f0)
					  {
						fracStart = 0.1;
						fracDecrement = 0.01;
						if (_WE.verbosity>3)
						  {
							_WE.EH.println("Steepest descent line-search: contracting range...");
						  }
					  }
					else
					  {
						for (i=1;i<n;i++) xNew[i] = xOld[i] + alpha_0*gradDirection.getQuick(i-1);
						negLogPost = f_to_minimize_raw(xNew);
					  }
					// (2) Now resume loop with xNew set at this downhill point.
				  }
			  }
		  }

		if (_WE.verbosity>3)
		  {
			_WE.EH.println("Line-search: fixed at frac= %4.3g, f=%10.8g",fracBest, negLogPostBest);
		  }
		// copy best into xNew
		System.arraycopy(xBest,1,xNew,1,n-1);
		xTry = null;
		xBest = null;
	  }
	else
	  {
		// Traditional parabolic and cubic line-searches
		// Cheap, but leads to poorer solutions and more outer loop iterations
		double ALPHA = 1.0e-4, TOLX = 1.0e-7;
		double a,frac2=0.0,fracin,b,disc,f2=0.0;
		double rhs1,rhs2,slope,temp,test,tmpfrac;
		boolean check = false;

		negLogPostBest = f_to_minimize_raw(xOld);

		slope=0.0;
		for (i=1;i<n;i++) slope += (-gradF.getQuick(i-1))*(xNew[i]-xOld[i]);
		if (slope >= 0.0) _WE.EH.println("Roundoff problem in getBestPointAlongJoiningLine()");
		test=0.0;
		for (i=1;i<n;i++) 
		  {
		  temp=Math.abs((xNew[i]-xOld[i]))/Math.max(Math.abs(xOld[i]),1.0);
		  if (temp > test) test=temp;
		}
		fracin=TOLX/test;
		frac=1.0;
		for (;;) 
		  {
			for (i=1;i<n;i++) xTry[i]=xOld[i]+frac*(xNew[i]-xOld[i]);
			negLogPost=f_to_minimize_raw(xTry);
			if (_WE.verbosity>3)
			  {
				if (frac==1.0)
				  {
					_WE.EH.println("Line-search target point: f_post="+negLogPost);
					for (i=1; i<n; i++) 
					  {
						_WE.EH.printf("%9.6g ",((Property)modelVariables.mapFromIndex.get(i-1)).value);
					  }
					_WE.EH.println("");
				  }
				_WE.EH.println("Line-search %2d: %4.3f %10.8g %d",pass,frac,negLogPost,it_count-1);
			  }
			if (frac < fracin) 
			  {
				for (i=0;i<n;i++) xTry[i]=xOld[i];
				check=true;
				break;
			  } 
			else if ((negLogPost <= negLogPostBest+ALPHA*frac*slope))
			  {
				// normal exit
				break;
			  }
			else 
			  {
				if (frac == 1.0)
				  {
					tmpfrac = -slope/(2.0*(negLogPost-negLogPostBest-slope));
				  }
				else 
				  {
					rhs1=negLogPost-negLogPostBest-frac*slope;
					rhs2=f2-negLogPostBest-frac2*slope;
					a=(rhs1/(frac*frac)-rhs2/(frac2*frac2))/(frac-frac2);
					b=(-frac2*rhs1/(frac*frac)+frac*rhs2/(frac2*frac2))/(frac-frac2);
					if (a == 0.0) 
					  {
						tmpfrac = -slope/(2.0*b);
					  }
					else 
					  {
						disc=b*b-3.0*a*slope;
						if (disc < 0.0) tmpfrac=0.5*frac;
						else if (b <= 0.0) tmpfrac=(-b+Math.sqrt(disc))/(3.0*a);
						else tmpfrac=-slope/(b+Math.sqrt(disc));
					  }
					if (tmpfrac>0.5*frac) tmpfrac=0.5*frac;
				  }
			  }
			frac2=frac;
			f2 = negLogPost;
			frac = Math.max(tmpfrac,0.01*frac);    // minimum step of 0.01
		  }
		if (check)
		  {
			// xNew too close to xOld. Reset to xOld
			System.arraycopy(xOld,1,xNew,1,n-1);
			frac = 0.0;
			negLogPost = negLogPostBest;
		  }
		// breaks exit here
		System.arraycopy(xTry,1,xNew,1,n-1);
		fracBest = frac;
		negLogPostBest = negLogPost; 

		if (_WE.verbosity>3)
		  {
			_WE.EH.println("Line-search: fixed at frac= %4.3g, f=%10.8g",fracBest,negLogPostBest);
		  }
		xTry = null;
		xBest = null;

		return;
	  }
  }

  /* 
	 Calls here will be with Linverse = Cholesky factor L of permuted form of Cpost: P.L.L^T.P^T = C^{-1} for some P
  */
  private double[] getLargestEigenDirection(ModifiedCholeskyLeftDecomposition cholDec)
  {
	int i,j,k,N,M;
	double[] diagScale;
	double sum;
	DoubleMatrix2D L_post=null, Linverse;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D Cperm;


	N = modelVariables.nParams;
	double[] eValues = new double[N];

	_WE.EH.println("getLargestEigenDirection(): rebuild with section uncommented, linked to jlapack files etc!");
	System.exit(0); 
	
	// Commented out of code to remove jlapack dependencies (this routine never called anyway)

	// SECTION HERE SHOULD BE UNCOMMENTED..........
	// Rebuild linking against  f2jutil.jar, blas.jar, lapack.jar, xerbla.jar
	// Commented out of code to remove jlapack dependencies (this routine never called anyway)
	/*
	
	Linverse = cholDec.getL();

	// form permuted posterior C
	Cperm = 
	  new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(N,N);
	if (!true)
	  {
		// FIXME !!!!! EIGENVALUES BELOW SHOULD SELECT SMALLEST, NOT LARGEST
		// SEE (below) SMALLEST
		L_post = Linverse;
		SeqBlas.seqBlas.dgemm(false,true,1,L_post,L_post,0,Cperm);
	  }
	else 
	  {
		L_post = Linverse.copy();
		ModifiedCholeskyLeftDecomposition.getCholeskyLeftInverse(L_post);
		SeqBlas.seqBlas.dgemm(true,false,1,L_post,L_post,0,Cperm);
	  }

	// scale it
	diagScale = diagonalScaleMatrix(Cperm);

	// Get M smallest eigenvalues (since we have Cinverse, not C)
	M=1;
	double[] eValues = new double[N];
	double[] eVectors = new double[N*M];
	double[] work = new double[8*N];
	int[] iwork = new int[5*N];
	int[] ifail = new int[N];
	org.netlib.util.intW m = new org.netlib.util.intW(0);  
	org.netlib.util.intW info = new org.netlib.util.intW(0); 
	org.netlib.lapack.Dsyevx.dsyevx("v","i","u",
									N,Cperm.getElements(),0,N, 
									0.0,0.0,
									N-M+1,N,    // SMALLEST 1...M, OR LARGEST N-M+1...N
									0.0,
									m,
									eValues, 0,
									eVectors,0,N, 
									work,0,work.length,
									iwork,0, 
									ifail, 0,
									info);
	for (i=0;i<N;i++) 
	  {
		// rescale to get true direction
		eVectors[i] = diagScale[i]*eVectors[i];
	  }
	int[] L_Permutation = cholDec.getPermutation(); 
	for (i=0,sum=0.0;i<N;i++) 
	  {
		// undo permutation
		j = L_Permutation[i];
		eValues[j] = eVectors[i];
		sum += eValues[j]*eValues[j];
	  }

	sum = Math.sqrt(sum);

	for (i=0;i<N;i++) 
	  {
		// rescale to unit vector
		eValues[i] /= sum;
	  }


	// GC
	diagScale = null;
	Cperm = null;
	work = null;
	eVectors = null;
	iwork=null;
	ifail=null;
	m=null;
	info=null;


	*/
	// ............... SECTION HERE SHOULD BE UNCOMMENTED

	return eValues;
  }

  /**
	 Print to files cross sections of -log(posterior)
	 after minimiser terminates
  */
  void debugPrintMinimiserXsections(String file_prefix,double[] centreParameters)
  {
	int i,ii,j,nG;
	double val,newval,scale;
	double[] fvals,xvals;
	Property prop;
	String filename = null;

	// Number of Gaussian parameters
	nG = modelVariables.nGaussian;
	double[] untranslatedParams = centreParameters.clone();

	// ensure in-situ state is set....
	mapParameters(centreParameters,true); 
	// set up raw parameters vector
	for (i=0;i<(modelVariables.nParams);i++) 
	  {    
		untranslatedParams[i+1] = ((Property)modelVariables.mapFromIndex.get(i)).value;
		// last entries in untranslatedParams are log(noise) values
		if (i>=nG) untranslatedParams[i+1] = Math.log(Math.abs(untranslatedParams[i+1]));
	  }	
	fvals = new double[2*N_XSECTION_STEPS+1];
	xvals = new double[2*N_XSECTION_STEPS+1];

	
	for (i=0; i<modelVariables.nParams; i++)
	  {
		filename = file_prefix+i+".txt";

		// debugging stuff. Leave infrastructure
		if (_WE.m.debuggingXSections)
		  {
			// Debugging infrastructure. Find a problematic Xsection, and insert here :-)
			// Turn on debuggingXSections
			if (filename.equals("preMinimiserXsection.0.R_near22.txt"))
			  {
				_WE.m.debug = true;
			  }
			else
			  {
				_WE.m.debug = false;
			  }
		  }

		val = untranslatedParams[i+1];
		scale = Math.abs(typical_values[i+1]);
		if (i>=nG)
		  { // A noise parameter
			if (!_WE.ignoreSeismic) scale = 0.2*val; // seems to be about right
			else scale = 0.2;                        // since mean will be 1
		  }
		prop = (Property)modelVariables.mapFromIndex.get(i);
		// time-to-depths: typical scale (prior sigma)
		if (i>=modelVariables.firstTDparamIndex && i<=modelVariables.lastTDparamIndex) 
		  {
			scale = 0.1*modelVariables.typicalScale[i];
		  }
		for (j=0,ii=-N_XSECTION_STEPS; ii<=N_XSECTION_STEPS; ii++,j++)
		  {
			newval = val + scale *((double)ii/(double)N_XSECTION_STEPS);
			untranslatedParams[i+1] = newval;
			// call function with perturbation
			fvals[j] = f_to_minimize_raw(untranslatedParams);
			xvals[j] = prop.value; // mapped parameter	(value set in f_to_minimize_raw()) (true in-situ variables...)
			if (_WE.m.debug) _WE.m.debug_cntr++;
		  }
		// restore value
		untranslatedParams[i+1] = val;
		f_to_minimize_raw(untranslatedParams);
		Util.printfileZO(filename,xvals,fvals,fvals.length);
	  }
  }





  /**
	 Print to files cross sections of -log(posterior)
  */
  void debugPrintModeXsections(String file_prefix,boolean printOnlyIfDodgy,boolean printOnlyIfDiscontinuousAtCentre)
  {
	int i,ii,j,k,nG;
	double val,newval,centreVal,scale;
	double[] currentParams;
	double[] fvals,xvals;
	Property prop;
	boolean isSatProp = false, isFarStretchProp=false;
	SaturationInfo _sInfo=null;

	// Number of Gaussian parameters
	nG = modelVariables.C_prior.rows();

	currentParams = new double[modelVariables.nParams+1];
	for (i=0;i<modelVariables.nParams;i++) currentParams[i+1] = modelVariables.mu_post.getQuick(i);

	fvals = new double[2*N_XSECTION_STEPS+1];
	xvals = new double[2*N_XSECTION_STEPS+1];
	
	for (i=0; i<modelVariables.nParams; i++)
	  {
		val = modelVariables.mu_post.getQuick(i);
		prop = (Property)modelVariables.mapFromIndex.get(i);
		isFarStretchProp=false;
		if (_WE.m.farStackWaveletTimeStretchList.contains(prop) && _WE.m.farStackWaveletTimeStretchMapper!=null) 
		  {
			isFarStretchProp = true;
		  }
		if (_WE.useSaturationModel)
		  {
			// lookup for saturations. Only variable ones appear, so they will have a mapper
			isSatProp = (_sInfo = modelVariables.getParentSatInfo(prop)) != null;
		  }	

		// go +/=- 2 "std deviations" (they aren't strictly this, its
		// really +/- 2 lengths l defined by l = (Ctilde^{-1})_{ii}
		scale = 2.0/Math.sqrt(modelVariables.CtildeInvDiag[i]);
		for (j=0,ii=-N_XSECTION_STEPS; ii<=N_XSECTION_STEPS; ii++,j++)
		  {
			newval = val + scale *((double)ii/(double)N_XSECTION_STEPS);
			currentParams[i+1] = newval;
			// call function with perturbation
			fvals[j] = f_to_minimize_raw(currentParams);
			xvals[j] = prop.value; // mapped parameter	(value set in f_to_minimize_raw()) (true in-situ variables...)
			// some properties may have been nonlinearly remapped. Retrieve true values here...
			if (isFarStretchProp)
			  {
				xvals[j] = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(prop.value);
			  }
			else if (isSatProp)
			  {
				xvals[j] = _sInfo.satPropMapper.getTrueValue(prop.value);
			  }		  
		  }
		// restore value
		currentParams[i+1] = val;
		centreVal = f_to_minimize_raw(currentParams);
		if (printOnlyIfDiscontinuousAtCentre)
		  {
			if ((Math.abs(fvals[N_XSECTION_STEPS-1] - centreVal) > DODGY_XSECTION_JUMP_THRESHOLD)
				||(Math.abs(fvals[N_XSECTION_STEPS+1] - centreVal) > DODGY_XSECTION_JUMP_THRESHOLD))
			{
			  Util.printfileZO(file_prefix+i+".jumpy_optimum.txt",xvals,fvals,fvals.length);
			}
		  }
		else if (printOnlyIfDodgy)
		  {
			// ML value doesn't seem to be right
			DoubleArrayList list = new DoubleArrayList(fvals);
			double max_entry = cern.jet.stat.Descriptive.max(list);

			if (max_entry < centreVal -DODGY_XSECTION_THRESHOLD)
			  {
				// XSection is dodgy
				Util.printfileZO(file_prefix+i+".dodgy_optimum.txt",xvals,fvals,fvals.length);
			  }
			// or scaling seems really off
			val = max_entry - centreVal;
			if (val< 1.5 || val > 2.5 )
			  {
				// XSection is dodgy
				Util.printfileZO(file_prefix+i+".dodgy_scaling.txt",xvals,fvals,fvals.length);
				if (!true)
				  {
					// curvature profile
					ii=N_XSECTION_STEPS;
					double[] curv = new double[N_XSECTION_STEPS];
					for (j=1;j<=N_XSECTION_STEPS;j++)
					  {
						curv[j-1] = (fvals[ii+j]+fvals[ii-j]-2*fvals[ii])/MathUtil.sqr(xvals[ii+j]-xvals[ii]);
					  }
					
					Util.printfileZO(file_prefix+i+".curvature.txt",curv,curv.length);
				  }
			  }
		  }
		else
		  {
			Util.printfileZO(file_prefix+i+".txt",xvals,fvals,fvals.length);
		  }
	  }
  }


  /**
	 Print to files cross sections of -log(posterior) along principal axis of Gaussian
	 approximation
  */
  void debugPrintModePrincipalXsections(String file_prefix,boolean printOnlyIfDodgy)
  {
	int i,ii,j,jj,k,indx,N,M,nG;
	double val,scale,lambda;
	double[] currentParams;
	double[] fvals,xvals;
	double centreVal=0;
	double fWorst=Double.MAX_VALUE;
	int dirWorst=-1;
	DoubleMatrix2D L_post;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D Cperm;
	EigenvalueDecomposition eDecomp;
	DoubleMatrix2D Lambda,uTransp;
	DoubleMatrix1D eM,eWork,mu_post;	
	int[] L_Permutation;
	double[] diagScale=null;

	// Number of Gaussian parameters
	nG = modelVariables.nGaussian;
	N = modelVariables.nParams;
	currentParams = new double[N+1];
	for (i=0;i<N;i++) currentParams[i+1] = modelVariables.mu_post.getQuick(i);


	fvals = new double[2*N_XSECTION_STEPS+1];
	xvals = new double[2*N_XSECTION_STEPS+1];
	
	// Require eigensystem of modelVariables.L_post...
	// ensure system is pure lower triangular
	L_post = modelVariables.L_post;
	// form symmetric permuted covariance Cperm = L_post^T.L_post
	Cperm = 
	  new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(N,N);
	SeqBlas.seqBlas.dgemm(true,false,1,L_post,L_post,0,Cperm);

	double[] eValues = new double[N];

	if (true)
	  {
		// now decompose
		eDecomp = new EigenvalueDecomposition(Cperm);
		Lambda = eDecomp.getD();
		uTransp = eDecomp.getV();
		for (i=0;i<N;i++) eValues[i] = Lambda.getQuick(i,i);
		// print all eigendirections
		M=N;
	  }
	else
	  {
		// Netlib (jlapack) routines for selected eigenvalues....UNCOMMENTED to remove jlapack dependencies

		/* FROM HERE.....

		// scale it
		diagScale = diagonalScaleMatrix(Cperm);
	
		// now decompose
		uTransp = Cperm.like();
		
		// scan M largest eigenvalues
		M=N;
		double[] eVectors = new double[N*M];
		double[] work = new double[8*N];
		int[] iwork = new int[5*N];
		int[] ifail = new int[N];
		org.netlib.util.intW m = new org.netlib.util.intW(0);  
		org.netlib.util.intW info = new org.netlib.util.intW(0); 
		org.netlib.lapack.Dsyevx.dsyevx("v","i","u",
										N,Cperm.getElements(),0,N, 
										0.0,0.0,
										N-M+1,N,
										0.0,
										m,
										eValues, 0,
										eVectors,0,N, 
										work,0,work.length,
										iwork,0, 
										ifail, 0,
										info);
		uTransp = Cperm.like().assign(0);
		for (i=0,k=0;i<M;i++) {
		  for (j=0;j<N;j++,k++) 
			{
			  uTransp.setQuick(j,i,eVectors[k]);
			}
		}
		*/
		// ..... TO HERE
	  }

	eM = modelVariables.mu_post.like();
	eWork = modelVariables.mu_post.like();
	mu_post = modelVariables.mu_post;
	L_Permutation = modelVariables.L_Permutation;

	for (i=0; i<M; i++)
	  {
		// eM = perturbations in principal directions i space
		eM.assign(0); 
		scale = Math.sqrt(Math.max(eValues[i],0.0));
		for (j=0,ii=-N_XSECTION_STEPS; ii<=N_XSECTION_STEPS; ii++,j++)
		  {
			lambda = 2.0*((double)ii/(double)N_XSECTION_STEPS);
			if (true)
			  {
				// unit vector in ith eigen-direction
				eM.setQuick(i,lambda*scale);
				uTransp.zMult(eM,eWork,1,0,false);
				if (diagScale != null)
				  {
					for (k=0; k<N; k++) eWork.setQuick(k, eWork.getQuick(k) * diagScale[k]);
				  }
				// now eWork  = uTransp.eM
			  }
			else
			  {
				// test: form eWork = L_post^T.eM
				eM.setQuick(i,lambda);
				eWork.assign(eM);
				SeqBlas.seqBlas.dtrmv(false, true, false, L_post,eWork);
			  }
			// now form P.eWork and write into currentParams[]
			for (jj=0;jj<N;jj++) 
			  {
				// undo permutation & add back mean
				indx = L_Permutation[jj];
				val = eWork.getQuick(jj) + mu_post.getQuick(indx);
				// Write it into in-situ field for the permuted index
				// prop = (Property)modelVariables.mapFromIndex.get(indx);
				// prop.setValue(val);  (has index indx)
				currentParams[indx+1] = val;
			  }
			// call function with perturbation
			fvals[j] = f_to_minimize_raw(currentParams);
			xvals[j] = lambda;
			if (ii==0) centreVal = fvals[j];
			if (fvals[j] < fWorst && ii!=0)
			  {
				fWorst = fvals[j];
				dirWorst = i;
			  }
		  }
		// centreVal = f_to_minimize_raw(currentParams);
		if (printOnlyIfDodgy)
		  {
			DoubleArrayList list = new DoubleArrayList(fvals);
			double min_entry = cern.jet.stat.Descriptive.min(list);

			if (min_entry < centreVal - DODGY_XSECTION_THRESHOLD)
			  {
				// XSection is dodgy
				Util.printfileZO(file_prefix+i+".txt",xvals,fvals,fvals.length);
			  }
		  }
		else
		  {
			Util.printfileZO(file_prefix+i+".txt",xvals,fvals,fvals.length);
		  }
	  }
	if (dirWorst>=0)
	  {
		_WE.EH.println("Xsection with worst mislocated minimum = "+dirWorst);
	  }
  }



  /**
	 Rescale matrix to unit-diagonal form. Return scaling vector
  */
  private double[] diagonalScaleMatrix(DoubleMatrix2D C_symm)
  {
	int i,j,N;
	N = C_symm.rows();

	double[] diagScale = new double[N];
	for (i=0;i<N;i++) 
	  {
		if (Math.abs(C_symm.getQuick(i,i)) > 0.0)
		  {
			diagScale[i] = Math.sqrt(Math.abs(C_symm.getQuick(i,i)));
		  }
		else
		  {
			diagScale[i] = 1.0;
		  }
	  }
	for (i=0;i<N;i++) {
	  for (j=0;j<N;j++) {
		C_symm.setQuick(i,j,C_symm.getQuick(i,j)/(diagScale[i]*diagScale[j]));
	  }
	}
	return diagScale;
  }


  /** 
	  A method to examine the posterior surface near the MAP (maximum aposteriori
	  likelihood) points. The Newton solver is sometimes fooled, and gets
	  stuck at little "ledges", close to a steep downhill slope, or in the upper
	  bowl of an asymmetric "W" shaped double minima.
	  To detect this we will check the posterior surface at points where the time values are sampled
	  at -1,-0.5,+0.5,+1 std deviations away from the MAP point, taken one at a time (like a very
	  simple Box Benkhen scheme). Points of MAP probablity lower than some Newton threshold less
	  than the MAP point are taken as evidence of the above pathology, and flagged.
	  If so, we return true, to indicate wierdness of the surface (probable false convergence)
  */
  private boolean modeHasFaultyNewtonConvergence(int i_trial)
  {
	int i,j,ii;
	double val,newval,MAPnegLogPost,scale;
	double[] currentParams;
	double[] bestParams;
	double fval,fbest;
	boolean looksFaulty = false;


	int NSTEPS=4; // left and right of MAP point

	currentParams = new double[modelVariables.nParams+1];
	for (i=0;i<modelVariables.nParams;i++) currentParams[i+1] = modelVariables.mu_post.getQuick(i);
	bestParams = currentParams.clone();


	// MAP point 
	MAPnegLogPost = f_to_minimize_raw(currentParams);
	fbest =	MAPnegLogPost - DODGY_XSECTION_THRESHOLD;

	for (i=0; i<modelVariables.nParams; i++)
	  {
		val = modelVariables.mu_post.getQuick(i);
		// prop = (Property)modelVariables.mapFromIndex.get(i);
		// go +/=- 2 std deviations
		scale = 2.0/Math.sqrt(modelVariables.CtildeInvDiag[i]);
		for (j=0,ii=-NSTEPS; ii<=NSTEPS; ii++,j++)
		  {
			if (ii==0) continue;
			newval = val + scale *((double)ii/(double)NSTEPS);
			currentParams[i+1] = newval;
			// call function with perturbation
			fval = f_to_minimize_raw(currentParams);
			// examine point
			if (fval < fbest)
			  {
				// XSection is dodgy
				looksFaulty = true;
				fbest = fval;
				copy(currentParams,bestParams);
			  }
		  }
		// restore value
		currentParams[i+1] = val;
		f_to_minimize_raw(currentParams);
	  }
	if (looksFaulty) 
	  {
		if (false)
		  {
			debugPrintModeXsections("FaultyNewtonXsection."+_WE.multiModelSetUp.modeNumber+"."+
									multistart_ctr+"."+i_trial+".",false,true);
		  }
		if (false)
		  {
			// print out Xsection in this direction......
			NSTEPS=30;
			double[] fvals = new double[2*NSTEPS+1];
			double[] xvals = new double[2*NSTEPS+1];
			ls_fullParameters = new DenseDoubleMatrix1D(modelVariables.nParams);

			for (j=0,ii=-NSTEPS; ii<=NSTEPS; ii++,j++)
			  {
				for (i=0; i<modelVariables.nParams; i++)
				  {
					val = modelVariables.mu_post.getQuick(i);	
					ls_fullParameters.setQuick(i,(bestParams[i+1] - val));
					newval = val + (bestParams[i+1] - val) * 2.0*((double)ii/(double)NSTEPS);
					currentParams[i+1] = newval;
				  }
				// call function with perturbation
				fvals[j] = f_to_minimize_raw(currentParams);
				xvals[j] = 2.0*((double)ii/(double)NSTEPS);
			  }
			Util.printfileZO("FaultyNewtonXsection."+_WE.multiModelSetUp.modeNumber+"."+
							 multistart_ctr+"."+i_trial+".txt",xvals,fvals,fvals.length);
		  }
		// embed new values in mu_post
		for (i=0;i<modelVariables.nParams;i++) modelVariables.mu_post.setQuick(i,bestParams[i+1]);
	  }
	// leave this function with best parameters in-situ
	f_to_minimize_raw(bestParams);

	return looksFaulty;
  }


  /**
	 Workspace set up
  */
  private void allocateWorkspace(int n_tot)
  {
	// needed workspace
	Hessian_wk = new double[n_tot+1][n_tot+1];
	Hessian_diag_wk = new double[n_tot+1];
	fullParameters_best = new double[n_tot+1];
	grad_sub_best = new double[n_tot+1];
	typical_values = new double[n_tot+1];
	fullParameters = new double[n_tot+1];
  }



  /**
   The solver throws execution to here everytime it gets a major new
   point. This is where we can monitor  things.
  */
  public void updateNewPoint(double fnew, double[] xnew)
  {
	if (_WE.verbosity>4) 
	  {
		_WE.EH.printf("its= "+it_count+" f= "+ fnew+ " params= ");
		for (int i=1; i<xnew.length; i++) 
		  {
			_WE.EH.printf("%9.6g ",xnew[i]);
		  }
		_WE.EH.println("");
	  }
  }

 public boolean abortNow()
  {
	return abortNow;
  }

  public boolean wasAborted()
  {
	return abortNow;
  }

 public void setFlags(boolean setOrReset)
  {
  }

  /**
	 Returns component of likelihood due to negLogPriorLikelihood
	 to be used only after a successful f_to_minimize() call
  */
  public double getNegLogPriorLikelihood()
  {
	return negLogPriorLikelihood;
  }

  /**
	 Returns component of likelihood due to negLogSeisLikelihood
	 to be used only after a successful f_to_minimize() call
  */
  public double getNegLogSeisLikelihood()
  {
	return negLogSeisLikelihood;
  }


  // null methods required since we can't figure these out analytically
  public void gradient(double times[], double g[])
  {};
  public void hessian(double times[], double h[][]) 
  {};
}
