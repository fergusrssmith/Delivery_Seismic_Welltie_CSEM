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



package au.csiro.JamesGunning.seismicinverter4D;

import java.util.*;


import cern.colt.matrix.*;
import cern.colt.matrix.linalg.*;

import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;

import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.util.MathUtil;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;


/**
   A class to draw samples of the model using the MCMC technique, with a Metropolis algorithm
   and smaller scaled local hops. Mainly an extension of the ModelSampler class.

   Implements a Metropolis sampling algorithm described in Section 11.5 of Gelman et al,
   Bayesian Data Analysis (Chapman and Hall). The Metropolis sampler alternates
   mode-jumping (different fluid combos and or different minima in the likelihood surface
   for the same fluids) with small hops from from a MultiGaussian distro centred at the current
   point and with a covariance set at the estimated mode covariance scaled by a magic small 
   number (about 2.4/srqt(d)), as described in Gelman.
*/
public final class RandomWalkSampler
{
  /**
	 pointer to underlying model
  */
  Model _m;
  /**
	 and ModelSampler
  */
  ModelSampler _MS;
  /**
	 and ErrorHandler
  */
  ErrorHandler _EH;

  /**
	 The model-jumping random walk sampler will propose model-jumps with probability
	 equal to this constant. We suggest 50% prob. of trying a model jump, and 50%
	 doing a random hop within a model
  */
  public final static double MODEL_JUMPING_PROPOSAL_RATE=0.5;

  /**
	 When making a model-jump, we will propose a new mode in proportion to its estimated
	 posterior weight if this is true. If false, hybrid weights are used
  */
  public final static boolean USE_POSTERIOR_WEIGHTS_IN_MODE_JUMP_PROPOSALS = true;
/**
	 When making a model-jump, we can force this to flip into another mode
	 by switching this on. It seesm a good idea to promote mixing. The biasing
	 is compensated for in the acceptance step
  */
  public final static boolean MODE_JUMPS_MUST_LEAVE_OLD_MODE = true; // Best left TRUE

  /**
	 Sampling efficiency according to Gelman is about 0.3/d for a d-dimensional
	 Gaussian distribution sampled by a hopping method with proposal scale covariance
	 equal to about 2.4/sqrt(d) times the true posterior covariance. This is the
	 efficiency constant.
  */
  public static final double GELMAN_EFFICIENCY_CONSTANT = 0.3;

  /**
	 Read blurb for GELMAN_EFFICIENCY_CONSTANT. This is the
	 Numerator constant from which we compute the initial scale.
	 I.e. scale factor = GELMAN_NUMERATOR_CONSTANT/sqrt(d).
  */
  public static final double GELMAN_NUMERATOR_CONSTANT = 2.4;

  /**
	 Boolean switch that controls whether the random walker 
	 attached to a mode is allowed to jump into a different pinchout
	 configuration
  */
  public static final boolean forbidLayerCodeChangingJumps = false; 
 // FIXME: would like to apply this but we sometimes seem to get stuck in Voronoi problems 

  private boolean modeHopIsPossible=true;  // useful to track mode hopping permissibility

  /**
	 This is the numerator constant from which we compute the initial
	 scale for MALA style samplers. See G.O.Roberts, "Optimal scaling
	 for various Metropolis Hastings algorithms",
	 Stat. Sci. 16(4),351-367
  */
  public static final double MALA_NUMERATOR_CONSTANT = 1.31;

  /**
	 Switching this on allows for occasional big steps in properties that may be
	 thresholded, and thus susceptible to slow diffusion outside the truncation range
	 N.B. this seems to break the 3 layer wedge model when it's on....
  */
  public static final boolean USE_RANDOM_BIGSTEPS_IF_REQUIRED = true;

  /**
	 This constant defines the number of posterior standard deviations
	 from the posterior mean a truncated-property boundary must be in
	 order for the big-jump random sampler not to be activated.
	 Eg, if <NG>+BIG_JUMP_ACTIVATION_MARGIN * sigma_NG > 1.0, we
	 will need to activate random big-jumps for the NG property.
  */
  public static final double BIG_JUMP_ACTIVATION_MARGIN = 3.0;  

  private int mc_cntr=0;
  private int fs_old;
  private double negLogLikSeisIsopachDepthVS, negLogLikPriorOnly;
  private LocalMode last_accepted_mode;
  private double negLogQtheta,negLogDetCphi;

  boolean useMALAsampler;
  private double q_MALA=0.0;




  // Big step flag for properties that wander into truncated ranges
  private boolean allowRandomBigSteps;
  // Big step probability for within-mode hpos
  private double randomBigStepsProbability=0.25;

  private long debug_cnt=0;
  /**
	 Draw samples of the model using the MCMC technique, with a Metropolis algorithm
	 and smaller scaled local hops
  */

  RandomWalkSampler(Model m)
  {
	this._m = m;
	this._MS = m.MS;
	this._EH = m._SI.EH;
	useMALAsampler = false; // don't seem to get many improvements out of MALA
	allowRandomBigSteps = USE_RANDOM_BIGSTEPS_IF_REQUIRED;
  }

  void sampleFromLocalModesUsingRandomWalkSampler()
  {
	if (false) 
	  {
		// Importance sampling across ensemble. Not correct yet
		// FIXME this needs testing+fixing (low priority)
		sampleFromLocalModesUsingEnsembleRandomWalkSampler();
	  }
	else
	  {
		// Andrieu style model jumping and RWS
		sampleFromLocalModesUsingModelJumpingRandomWalkSampler();		
	  }
  }

  // Ensemble-looping sampler methods

  /**
	 Here we run chains in parallel for all modes, and sample the
	 actual mode state by the current cdf computed from the likelihoods
  */

  private void sampleFromLocalModesUsingEnsembleRandomWalkSampler()
  {
	int i,j,k,realisation_number=0;
	int realisation_label_number=0, decimation_constant=1;
	int n_modes;
	boolean newState = false, write_state=false;
	boolean allowTraceOutput=true;
	double log_r,log_r_old=1.0;
	int[] N_fs;
	int[] N_jumps_accepted;
	double[] running_mode_cdf;
	double[] rel_log_lik;
	boolean[] accept_flag;
	double negLogLikMin;
	double val;

	// variables for mode moments analysis.......
	boolean getMoments=false;
	DenseDoubleMatrix1D[]  exper_S=null;
	DenseDoubleMatrix1D   _exper_S=null;    // useful ptr	
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D[]  exper_SS=null;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D   _exper_SS=null;   // useful ptr
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D[]  exper_samples=null;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D    exper_negLogLiks=null;
	// ..........to here

	if (_m._SI.verbosity>3) 
	  {
		_EH.println("Sample mode fluid-state | layers present, code | full-logLik prior-log-lik seis+iso+depthSClogLik | "+
						  "| log-lik-ratio hop-success chosen-mode%n");
	  }

	// set up vector of last states for each mode
	n_modes = _MS.localModes.size();
	for (k=0;k<n_modes;k++)
	  {
		_MS.mode = _MS.localModes.get(k);
		// set up initial state of MCMC chain
		// Note! use the mean vector here, as this will always generate a valid initial state
		// that has an acceptable pinchout configuration
		_MS.mode.last_accepted_state =  _MS.mode.QP_optimum.toArray();
		_MS.mode.last_accepted_R_pp = new double[_m.SM.Nstacks][_m.n_layers];
		if (_m._SI.doSimpleFluidUnitsAndContacts) 
		  {
			_MS.mode.last_accepted_R_pp_fc = new double[_m.SM.Nstacks][_m.n_layers];
			_MS.mode.last_accepted_t_fc = new double[_m.Nvintages][_m.n_layers];
		  }
	  }

	// output will be decimated according to Gelman's rule-of-thumb (assumes optimal Hopping sampling)
	if (_m._SI.decimateRandomWalkSamplerOutput) 
	  {
		decimation_constant = (int)(((double)_MS.max_mode_size)/GELMAN_EFFICIENCY_CONSTANT);
	  }

	// If we expressly suppress realisations, forbid output
	if (!_m._SI.dumpRealisations) allowTraceOutput = false;

	// writing may be suppressed if this is a calibration run
	if (_MS.CM != null && _MS.CM.RafteryStyleDiagnostics)
	  {
		if (_MS.CM.isRafteryCalibrationRun())
		  {
			allowTraceOutput = false;
		  }
		else
		  {
			// on second pass, will need to thin out, or on first pass, if we aren't
			// doing Raftery-auto-decimation, this will be a factor of unity
			decimation_constant *= _MS.CM.getRafteryWorstCaseThinningRecommendation(); 
		  }
	  }

	// loop over number of samples in MCMC chain
	int Nsamples=_m._SI.n_samples_per_trace_undecimated * decimation_constant;
	// but this may perhaps be overridden
	if (_m._SI.haveSpecifiedPreDecimationNsamples) 
	  {
		Nsamples=_m._SI.n_samples_per_trace_undecimated;
	  }
	
	if (getMoments)
	  {
		exper_S = new DenseDoubleMatrix1D[n_modes];
		exper_SS = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D[n_modes];
		for (k=0;k<n_modes;k++)
		  {
			_MS.mode = _MS.localModes.get(k);
			exper_S[k] = new DenseDoubleMatrix1D(_MS.mode.last_accepted_state.length);
			exper_SS[k] = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(_MS.mode.last_accepted_state.length, _MS.mode.last_accepted_state.length);
		  }
		// actual samples stuff
		exper_samples = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D[n_modes];
		for (k=0;k<n_modes;k++)
		  {
			_MS.mode = _MS.localModes.get(k);
			exper_samples[k] = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(Nsamples, _MS.mode.last_accepted_state.length);
		  }
		exper_negLogLiks = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(n_modes,Nsamples);		
	  }

	// place to stuff temporary samples as candidates
	_MS.property_sample = new double[_MS.max_mode_size];



	// local cdf acquired across all modes
	running_mode_cdf = new double[n_modes];
	// relative log-likelihoods for feedback
	rel_log_lik = new double[n_modes];
	accept_flag = new boolean[n_modes];

	// running estimates of acceptance rates
	N_jumps_accepted = new int[n_modes];
	Arrays.fill(N_jumps_accepted,0);
	N_fs = new int[_m.FS.possibleFluidStates.size()];

	// Now generate samples
	for (mc_cntr=1;mc_cntr<=2*Nsamples;mc_cntr++)
	  {
		// flag to signify new transition in MCMC chain
		newState = false;

		// Compute sample number
		realisation_number = mc_cntr-Nsamples;		

		// (1) sampling over m - marginal

		// Loop over all modes. Here we are sampling from the marginal for all the m_k, with
		// the model-label k summed out. The marginal is just the product of likelihood-posteriors
		// over all models. We should get a faux Brownian motion running for all modes
		// Each mode is then sampling a posterior that honours prior, isopach, plus seismic.
		for (k=0;k<n_modes;k++)
		  {
			if (_m._SI.verbosity>3) _EH.printf("%4d %3d ",mc_cntr, k);

			_MS.setMode(_MS.localModes.get(k));

			// Make Hop
			formHoppingProposalDraw();

			// Get likelihood (prior + isopach + seismic added here)
			log_r = computeHoppingLikelihood(true);
			
			// Apply Metropolis rule, using this mode's last state
			// (last accepted:  prior + isopach + seismic)
			log_r_old = _MS.mode.lastAcceptedNegLogAllLikelihoods + 
			            _MS.mode.lastAcceptedNegLogLikPriorOnly;
			rel_log_lik[k] = log_r - log_r_old;

			// Test to see if we accept the new state using the Metropolis rule
			if (_m._SI.verbosity>3) 
			  {
				_EH.printf("%14.8g ",rel_log_lik[k]); // the log(likelihood ratio)
			  }
			// Apply Metropolis rule.
			// Make comparison with care to avoid overflow in wierd cases. 
			// Also always accept the first candidate to start the chain
			if ((log_r <= log_r_old) || _MS.Ugen.nextDouble() <= Math.exp(-(log_r - log_r_old)) || (mc_cntr==1))
			  {
				accept_flag[k]=true;
				if (_m._SI.verbosity>3) _EH.print("Ah  -%n");
				fs_old = _MS.fs;
				newState = true;
				// remember state
				System.arraycopy(_MS.property_sample,0,_MS.mode.last_accepted_state,
								 0,_MS.mode.last_accepted_state.length);
				ArrayUtils.copy(_m.SM.R_pp,_MS.mode.last_accepted_R_pp);
				//ArrayUtils.copy(_m.SM.R_far,_MS.mode.last_accepted_R_far);
				if (_m._SI.doSimpleFluidUnitsAndContacts) 
				  {
					ArrayUtils.copy(_m.SM.R_pp_fc,_MS.mode.last_accepted_R_pp_fc);
					ArrayUtils.copy(_m.SM.t_fc,_MS.mode.last_accepted_t_fc);
				  }
				_MS.mode.lastAcceptedNegLogAllLikelihoods = negLogLikSeisIsopachDepthVS;
				_MS.mode.lastAcceptedNegLogLikPriorOnly = negLogLikPriorOnly;

				if (realisation_number>0) N_jumps_accepted[k]++;
			  }
			else
			  {
				accept_flag[k]=false;
				if (_m._SI.verbosity>3) _EH.print("Rh  -%n");
			  }

			// Acquire running cdf for the conditional p(k|M) which we will nedd in the next bit
			// running_mode_cdf[k] =  _MS.mode.lastAcceptedNegLogAllLikelihoods - Math.log(_MS.fluid_state.prob);
			running_mode_cdf[k] =  _MS.mode.lastAcceptedNegLogAllLikelihoods - _MS.fluid_state.log_prob;
			// Get moments for later normalisation constant estimation
			if (getMoments && realisation_number>0)
			  {
				// remember sample
				exper_negLogLiks.setQuick(k,realisation_number-1, _MS.mode.lastAcceptedNegLogLikPriorOnly
										  + _MS.mode.lastAcceptedNegLogAllLikelihoods);
				// record the full sample.......
				DoubleMatrix1D row = exper_samples[k].viewRow(realisation_number-1);
				row.assign(_MS.mode.last_accepted_state);

				// covariance. ptr
				_exper_SS =	exper_SS[k];
				_exper_S = exper_S[k];
				for (i=0;i<_MS.mode.last_accepted_state.length;i++)
				  {
					val = _MS.mode.last_accepted_state[i];
					_exper_S.setQuick(i,_exper_S.getQuick(i) + val);
					// covariance. Need upper triangle only
					for (j=i; j<_MS.mode.last_accepted_state.length; j++)
					  {
						_exper_SS.setQuick(i,j, _exper_SS.getQuick(i,j) + val*_MS.mode.last_accepted_state[j]);
					  }
				  }
			  }
		  }

		// (2) Sampling over the conditional p(k|M)

		// Now sample the model-label from the conditional.
		// (A Metropolis-Hastings step with unit acceptance probability)
		// The cdf computes the  negloglikelihood from seismic, isopach, plus prior fluid combo terms.
		// nedd to construct this carefully, watching out for overflow etc
		negLogLikMin = Double.MAX_VALUE;
		for (k=0;k<n_modes;k++) negLogLikMin = Math.min(negLogLikMin, running_mode_cdf[k]);
		// get normalisation constant
		double norm=0;
		for (k=0;k<n_modes;k++) 
		  {
			running_mode_cdf[k] = Math.exp(-running_mode_cdf[k] + negLogLikMin); 
			norm += running_mode_cdf[k];
		  }
        // normalise and turn into a pdf
		running_mode_cdf[0] /= norm;
		for (k=1;k<n_modes;k++) 
		  {
			running_mode_cdf[k] = running_mode_cdf[k-1] + running_mode_cdf[k]/norm;
		  }
		// make last cdf entry exact
		running_mode_cdf[n_modes-1] = 1.0;
		// At this point running_mode_cdf is an exact cdf for the possible states k=0...n_modes-1

		// and sample!
		double cdf = _MS.Ugen.nextDouble();
		// find associated local mode
		for (k=0;k<n_modes;k++)
		  {
			if (cdf<=running_mode_cdf[k]) break;
		  }
		// The sampled model now has index k. Reset mode
		_MS.setMode(_MS.localModes.get(k));
		// restore Layer properties
		restoreStateFromLastAcceptedState();

		// indicate chosen mode & give some feedback
		if (_m._SI.verbosity>3) 
		  {
			_EH.printf("%4d %3d %2d         | ",mc_cntr, k, _MS.mode.fs);
			for (i=0;i<_m.n_layers;i++) 
			  {
				_EH.printf("%1d ",(_MS.layer_present[i] ? 1 : 0));
			  }
			_EH.printf("%5d | %10.6f %10.6f %10.6f | %14.8f %s %3d%n",_MS.layer_code,
							_MS.mode.lastAcceptedNegLogAllLikelihoods+
							_MS.mode.lastAcceptedNegLogLikPriorOnly,
							_MS.mode.lastAcceptedNegLogLikPriorOnly,
							_MS.mode.lastAcceptedNegLogAllLikelihoods,
							rel_log_lik[k], ((accept_flag[k])? "Ah" : "Rh"), k);
		  }

		// Produce output only for second half of chain
		if (realisation_number>0)
		  {
			// Apply decimation of output to reduce correlations if necessary
			if (!_m._SI.decimateRandomWalkSamplerOutput || ((realisation_number-1) % decimation_constant == 0))
			  {
				write_state=true;
				realisation_label_number++;
			  }
			else
			  {
				write_state=false;
			  }
			// increment fluid state prob cntr
			N_fs[_MS.mode.fs]++;
			// ensure first member of second half of chain is set and dumped out
			boolean isNewState;
			// force new state so properties in output trace are remapped
			// isNewState = (realisation_number==1)? true : newState;
			isNewState=true;
			// dump model
			_m._SI.dumpRealisation(isNewState, realisation_label_number, write_state && allowTraceOutput,false);
			// and seismics if required. Need to recompute synthetic seismic if dump needed
			if (_m._SI.dumpSyntheticSeismic) 
			  {
				_MS.negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
				_MS.negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
				_m.SM.dumpSynthetics(isNewState, realisation_label_number, write_state && allowTraceOutput);
			  }
			// if doing diagnostics, get scalar estimands
			if ((_m._SI.performGelmanConvergenceDiagnostics || _m._SI.performRafteryConvergenceDiagnostics)
				&& write_state) 
			  {
				_MS.aquireScalarEstimands(realisation_label_number-1);
			  }
		  }
	  }

	// after finishing multiple parallel chains, estimate normalisation constants via normal approximation
	if (getMoments)
	  {
		for (k=0;k<n_modes;k++)
		  {
			_MS.setMode(_MS.localModes.get(k));

			// (1) form approximatiing normal distribution parameters
			
			// get means
			_exper_S =	exper_S[k];
			for (i=0;i<_MS.mode.last_accepted_state.length;i++)
			  {
				_exper_S.setQuick(i , _exper_S.getQuick(i)/((double)realisation_number));
			  }
		    // covariance. ptr (formula 14.1.7. in NR). Has to be single pass so swallow the roundoff error
			_exper_SS =	exper_SS[k];
			// and covariance
			for (i=0;i<_MS.mode.last_accepted_state.length;i++)
			  {
				val = _exper_S.getQuick(i);
				for (j=i;j<_MS.mode.last_accepted_state.length;j++)
				  {
					_exper_SS.setQuick(i,j,  (_exper_SS.getQuick(i,j) - ((double)realisation_number) * val * _exper_S.getQuick(j))
									   /((double)realisation_number-1.0));
					// symmetric entry 
					_exper_SS.setQuick(j,i, _exper_SS.getQuick(i,j));
				  }
			  }
			// get covariance inverse of _exper_SS, to facilitate approximations evaluation
			Algebra alg = new Algebra();
			DoubleMatrix2D Cinv =  alg.inverse(_exper_SS);

			// mean of the approximation posterior......
			DoubleMatrix1D mu_approx = _exper_S;
			DoubleMatrix1D tmp = new DenseDoubleMatrix1D(mu_approx.size());
			DoubleMatrix1D ans = tmp.like();

			int Nsamples_sofar = realisation_number;

			// (2) attempt to get normalisation constant. Sum over importance ratios etc
			double sum=0.0;
			for (mc_cntr=0;mc_cntr<Nsamples;mc_cntr++)
			  {
				DoubleMatrix1D row = exper_samples[k].viewRow(mc_cntr);
				for (j=0;j<mu_approx.size();j++) 
				  {
					tmp.setQuick(j,row.getQuick(j) - mu_approx.getQuick(j));
				  }
				// form C^1.tmp -> ends up in-place in "ans"
				Cinv.zMult(tmp,ans,1.0,0.0,false);
	
				// then tmp^T.C^-1.tmp/2.0
				double negLogLikApproximatePosterior = tmp.zDotProduct(ans)/2.0;

				val = exper_negLogLiks.getQuick(k,mc_cntr) - negLogLikApproximatePosterior;
				
				sum+=Math.exp(val);
			  }
			sum/=(double)Nsamples_sofar;
			_EH.println("Norm. const = "+ sum);
		  }
	  }
	if (_m._SI.verbosity>2) 
	  {
		for (k=0;k<n_modes;k++)
		  {
			_EH.println("Using "+(USE_POSTERIOR_WEIGHTS_IN_MODE_JUMP_PROPOSALS?"posterior":"hybrid")+
							  " proposal weights, acceptance rate of mode "+k+" = " + 
							  (double)(N_jumps_accepted[k])/(double)(realisation_number));
		  }
		for (k=0;k<_m.FS.possibleFluidStates.size();k++)
		  {
			_EH.printf("Prob of state fs=%d is %5g%n",k,N_fs[k]/(double)realisation_number+1);
		  }
	  }
  }

  /** 
	  Once the _MS.mode reference is set, this method restores
	  all the internal variables in Model.LayerProperties[] to the
	  state contained in mode.last_accepted_state[]
  */
  private void restoreStateFromLastAcceptedState()
  {
	for (int i=0; i<_MS.mode.last_accepted_state.length; i++)
	  {
		// restore the property for this index
		_MS.setValuesFromIndex(i,_MS.mode.last_accepted_state[i]);
	  }
	// Ensure the times are ordered
	_MS.orderTimes();
	// must now restore internal vp_effs etc as these will be written out
	// exception will never be thrown, as we've been here before
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	
	// and restore reflection coefficients
	ArrayUtils.copy(_MS.mode.last_accepted_R_pp,_m.SM.R_pp);   
	if (_m._SI.doSimpleFluidUnitsAndContacts) 
	  {
		ArrayUtils.copy(_MS.mode.last_accepted_R_pp_fc,_m.SM.R_pp_fc);
		//ArrayUtils.copy(_MS.mode.last_accepted_R_far_fc,_m.SM.R_far_fc);
		ArrayUtils.copy(_MS.mode.last_accepted_t_fc,_m.SM.t_fc);
	  }
  }



  /**
	 Alternating random walk sampler.
  */
  void sampleFromLocalModesUsingModelJumpingRandomWalkSampler()
  {
	int k,realisation_number=0;
	int realisation_label_number=0, decimation_constant=1;
	boolean accept=false, newState = false, write_state=false;
	boolean modeJumpsMustLeaveOldMode;
	boolean useStrictPosteriorWeights = false;
	double log_r=0.0,log_r_old=1.0;
	long timeStart=0;
	int Nmodes,Nsamples;
	int[] N_fs;
	int N_mode_jumps_accepted=0, N_jumps_accepted=0;
	int N_mode_jumps_attempted=0, N_jumps_attempted=0;
	boolean mappingSucceeded,allowOutput,allowTraceOutput=true;
	double test_deviate_for_mode_jump;
	double mode_jump_proposal_rate;
	Algebra alg=null;
	double debug_HM_sum=0;   // FIXME HARMONIC MEANS
	boolean dumpState;

	// It's impossible unless we have some modes to start with
	if (_MS.localModes.size() < 1)
	  {
		_EH.println("No local modes to base sampler on. Optimiser must have failed. No RWS sampling possible!");
		return;
	  }

	mode_jump_proposal_rate = MODEL_JUMPING_PROPOSAL_RATE;      // (0.5) = fraction of time to spend trying to jump modes
	modeJumpsMustLeaveOldMode = MODE_JUMPS_MUST_LEAVE_OLD_MODE; // (true) force mode jump if we try a mode jump
	Nmodes = _MS.localModes.size();
	if (Nmodes==1) mode_jump_proposal_rate = 0.0;
    // set up vector of last states for each mode and other useful stuff
	if (useMALAsampler) alg = new Algebra();
	for (k=0;k<Nmodes;k++)
	  {
		_MS.mode = _MS.localModes.get(k);
		// set up initial state of MCMC chain
		_MS.mode.last_accepted_state =  _MS.mode.QP_optimum.toArray();
		_MS.mode.last_accepted_R_pp = new double[_m.SM.Nstacks][_m.n_layers];
		//_MS.mode.last_accepted_R_far = new double[_m.n_layers];
		if (_m._SI.doSimpleFluidUnitsAndContacts) 
		  {
			_MS.mode.last_accepted_R_pp_fc = new double[_m.SM.Nstacks][_m.n_layers];
			_MS.mode.last_accepted_t_fc = new double[_m.Nvintages][_m.n_layers];
		  }
		if (useMALAsampler)
		  {
			_MS.mode.RtranspInv = alg.inverse(_MS.mode.Rtransp);
		  }
		// Useful test of cross correlations, possible sparsity gains etc
		// testSparsityOfRmatrices();
		//debugGuessPosteriorErrorBars();
	  }
	alg=null;


	if (allowRandomBigSteps) setUpBigStepProposals();

	// place to stuff temporary samples as candidates
	_MS.property_sample = new double[_MS.max_mode_size];

	// output will be decimated according to scaling for approppriate sampler
	if (_m._SI.decimateRandomWalkSamplerOutput) 
	  {
		if (useMALAsampler)
		  { // Roberts' convergence time estimate
			decimation_constant = (int)(Math.pow((double)_MS.max_mode_size,0.33333));
		  }
		else
		  { // Gelman's rule-of-thumb (assumes optimal Hopping sampling)
			decimation_constant = (int)(((double)_MS.max_mode_size)/GELMAN_EFFICIENCY_CONSTANT);
		  }
	  }

	// If we expressly suppress realisations, forbid output
	if (!_m._SI.dumpRealisations) allowTraceOutput = false;

	// writing may be suppressed if this is a calibration run
	if (_MS.CM != null && _MS.CM.RafteryStyleDiagnostics)
	  {
		if (_MS.CM.isRafteryCalibrationRun())
		  {
			allowTraceOutput = false;
		  }
		else
		  {
			// on second pass, will need to thin out, or on first pass, if we aren't
			// doing Raftery-auto-decimation, this will be a factor of unity
			decimation_constant *= _MS.CM.getRafteryWorstCaseThinningRecommendation(); 
		  }
	  }

	// Number of samples in MCMC chain
	Nsamples=_m._SI.n_samples_per_trace_undecimated * decimation_constant;
	// but this may perhaps be overridden
	if (_m._SI.haveSpecifiedPreDecimationNsamples) 
	  {
		Nsamples=_m._SI.n_samples_per_trace_undecimated;
	  }

	if (_m._SI.verbosity>3) 
	  {
		_EH.println("Sample fluid-state mode | layers present, code | full-logLik prior-logLik seis+iso+depthSC-logLik | "+
						  "| log-lik-ratio success%n");
	  }
	// fluid state prob cntrs
	N_fs = new int[_m.FS.possibleFluidStates.size()];
	
	last_accepted_mode = null;
	timeStart = System.currentTimeMillis();
	// Now generate samples
	for (mc_cntr=1;mc_cntr<=2*Nsamples;mc_cntr++) // number of samples
	  {
		// flag to signify new transition in MCMC chain
		newState = false;

		// Compute sample number
		realisation_number = mc_cntr-Nsamples;

		if (_m._SI.verbosity>3) _EH.printf("%4d ",mc_cntr);
		// do a hop or a jump according to this deviate....
		test_deviate_for_mode_jump = _MS.Ugen.nextDouble();

		if ((test_deviate_for_mode_jump < mode_jump_proposal_rate) || (mc_cntr == 1))
		  {
			// Mode jumps (only required if there is more than 1 mode)

			// At initialisation (mc_cntr==1), use the posterior weights to pick the first mode
			useStrictPosteriorWeights = (USE_POSTERIOR_WEIGHTS_IN_MODE_JUMP_PROPOSALS || (mc_cntr==1));
			_MS.proposeNewMode(useStrictPosteriorWeights, modeJumpsMustLeaveOldMode,last_accepted_mode);
			
			// perform a state-jump based on mappings as per Andrieu et al
			
			// initialisation chunk........
			allowOutput=true;
			if (mc_cntr==1)
			  {
				// Nothing to form mode-jump mapping with. Draw a state using formHoppingProposalDraw()
				// force the old state to be the proposed state. This will in turn force acceptance.
				last_accepted_mode = _MS.mode;
				
				// make a draw from the mode
				formHoppingProposalDraw();
				computeHoppingLikelihood(allowOutput);
				takeCopyOfModeVariables();
				allowOutput=false;
			  }
			// ...initialisation chunk
			
			mappingSucceeded = formMappedModeJumpProposalDraw();

			if (mappingSucceeded)
			  {
				// compute synthetics etc
				computeHoppingLikelihood(allowOutput);
				
				log_r = computeNewStateImportanceRatio(modeJumpsMustLeaveOldMode,useStrictPosteriorWeights);

				// (1c) Apply Metropolis-Hastings rule
				if (_m._SI.verbosity>3) 
				  {
					_EH.printf(" Q=%10.6f | %14.8g ",negLogQtheta, log_r); // the log(likelihood ratio)
				  }
			  }
			// Metropolis criterion
			// make comparison with care to avoid overflow in wierd cases
			if (mappingSucceeded && ((log_r <= 0) || _MS.Ugen.nextDouble() <= Math.exp(-log_r)))
			  {
				// proposal accepted.
				if (_m._SI.verbosity>3) _EH.print("Am%n");
				fs_old = _MS.fs;
				newState = true;
				// remember state
				last_accepted_mode = _MS.mode;
				takeCopyOfModeVariables();
				if (realisation_number>0) N_mode_jumps_accepted++;
			  }
			else
			  {
				if (_m._SI.verbosity>3) _EH.print("Rm%n");
			  }
			// count stats only for second half of chain
			if (realisation_number>0) N_mode_jumps_attempted++;
		  }
		else
		  {
			// Within-mode hops

			// Step 2. Keep fluid-state/mode constant. Do Langevin diffusion on current exposed state
			// 2a. keep mode fixed. Reset this if last mode-change proposal failed
			_MS.setMode(last_accepted_mode);
			allowOutput=true;


			// 2b. make Hop		
			if (allowRandomBigSteps && _MS.mode.bigStepInfo!=null && 
				(_MS.Ugen.nextDouble() <= randomBigStepsProbability))
			  {
				formBigStepProposalDraw();
			  }
			else
			  {
				formHoppingProposalDraw();					
			  }

			// 2c. get likelihood
			log_r = computeHoppingLikelihood(allowOutput);
			
			if ( _m._SI.verbosity>3) 
			  {// for absence of Q(theta)
				_EH.print("    -         | ");
			  }

			// (7) Apply Metropolis rule
			log_r_old = last_accepted_mode.lastAcceptedNegLogAllLikelihoods + 
			  last_accepted_mode.lastAcceptedNegLogLikPriorOnly;

			// Test to see if we accept the new state using the Metropolis rule
			if (_m._SI.verbosity>3) 
			  {
				_EH.printf("%14.8f ",log_r - log_r_old - (useMALAsampler ? q_MALA:0));
				// the log(likelihood ratio)
			  }
			// make comparison with care to avoid overflow in wierd cases
			accept=false;
			if (useMALAsampler)
			  {
				// _EH.print("test_MALA = "+(log_r <= (log_r_old+q_MALA))+ " test = " +(log_r <= log_r_old));
				if ((log_r <= (log_r_old+q_MALA)) || 
					(_MS.Ugen.nextDouble() <= Math.exp(-(log_r - log_r_old)+q_MALA)))
				  {
					accept=true;
				  }
			  }
			else
			  {
				if ((log_r <= log_r_old) || _MS.Ugen.nextDouble() <= Math.exp(-(log_r - log_r_old)))
				  {
					accept=true;
				  }
			  }
			if (accept)
			  {
				if (_m._SI.verbosity>3) _EH.print("Ah%n");
				fs_old = _MS.fs;
				newState = true;
				// remember state
				takeCopyOfModeVariables();
				if (realisation_number>0) N_jumps_accepted++;
			  }
			else
			  {
				if (_m._SI.verbosity>3) _EH.printf("Rh%n");
			  }
			// count stats only for second half of chain
			if (realisation_number>0) N_jumps_attempted++;
		  } 
		// Timing. Average 10 samples.
		if ((_m._SI.decimateRandomWalkSamplerOutput && mc_cntr==Math.min(2*Nsamples,10*decimation_constant)) || 
			(!_m._SI.decimateRandomWalkSamplerOutput && mc_cntr==Math.min(2*Nsamples,10))) 
		  {
			if (_m._SI.verbosity>2)
			  {
				_EH.println("time per (written) sample approx. "+decimation_constant*
								  (System.currentTimeMillis() - timeStart)/(double)mc_cntr+" ms");
				_EH.println("time per (predecimation) sample approx. %6.4f ms",
								  (System.currentTimeMillis() - timeStart)/(double)mc_cntr);
			  }
		  }

		// Produce output only for second half of chain
		if (realisation_number>0) 
		  {
			// Apply decimation of output to reduce correlations if necessary
			if (!_m._SI.decimateRandomWalkSamplerOutput || ((realisation_number-1) % decimation_constant == 0))
			  {
				write_state=true;
				realisation_label_number++;
			  }
			else
			  {
				write_state=false;
			  }
			// fluid state cntr
			N_fs[fs_old]++;
			// ensure first member of second half of chain is set and dumped out
			if (realisation_number==1)
			  {
				dumpState = true;
				// Set all internal properties into a valid (reasonably probable) 
				// state before dumping into output trace
				// N.B. This causes a chunk of the preceding chain to get "erased", but this won't
				// matter statistically
				restoreStateFromLastAcceptedState();
			  }
			else
			  {
				dumpState = newState;
			  }

			// dump model
			_m._SI.dumpRealisation(dumpState, realisation_label_number, write_state && allowTraceOutput,false);

			// and seismics if required
			_m.SM.dumpSynthetics(dumpState, realisation_label_number, write_state && allowTraceOutput);
			// get convergence stats if required
			if ((_m._SI.performGelmanConvergenceDiagnostics || _m._SI.performRafteryConvergenceDiagnostics) 
				&& write_state)
			  {
 				_MS.aquireScalarEstimands(realisation_label_number-1);
			  }
			// FIXME test of harmonic mean summing for Bayes factors....
			if (false && fs_old==0) 
			  {
				debug_HM_sum+=1.0/Math.exp(-last_accepted_mode.lastAcceptedNegLogAllLikelihoods);
				if (_m._SI.verbosity>3) _EH.println("Harmonic_mean= "+1.0/(debug_HM_sum/N_fs[0]));
			  }
		  }
	  }
	if (_m._SI.verbosity>2 && Nsamples>0) 
	  {
		_EH.printf("Using "+(USE_POSTERIOR_WEIGHTS_IN_MODE_JUMP_PROPOSALS?"posterior":"hybrid")+
						 " proposal weights, average acceptance rate of mode-jump proposals is %5g%n",
						 (double)(N_mode_jumps_accepted)/Math.max(N_mode_jumps_attempted,1));
		_EH.printf("Within-mode jumping acceptance rate = %5g%n",
						 ((N_jumps_attempted>0) ? (double)N_jumps_accepted/(double)N_jumps_attempted : 1));
		double norm =(double)(N_jumps_attempted+ N_mode_jumps_attempted);
		for (k=0;k<_m.FS.possibleFluidStates.size();k++)
		  {
			_EH.printf("Prob of state fs=%d is %5g%n",k,((norm > 0) ? N_fs[k]/norm : 0));
		  }
	  }	 
  }

  // Store useful things about current state
  private void takeCopyOfModeVariables()
  {
	System.arraycopy(_MS.property_sample,0,_MS.mode.last_accepted_state,
					 0,_MS.mode.last_accepted_state.length);
	ArrayUtils.copy(_m.SM.R_pp,_MS.mode.last_accepted_R_pp);
	if (_m._SI.doSimpleFluidUnitsAndContacts) 
	  {
		ArrayUtils.copy(_m.SM.R_pp_fc,_MS.mode.last_accepted_R_pp_fc);
		ArrayUtils.copy(_m.SM.t_fc,_MS.mode.last_accepted_t_fc);
	  }
	_MS.mode.lastAcceptedNegLogAllLikelihoods = negLogLikSeisIsopachDepthVS;
	_MS.mode.lastAcceptedNegLogLikPriorOnly = negLogLikPriorOnly;
  }













  // Test to see if more efficient hopping might be possible if R is sparse. Do
  // this if it looks like a good thing to do

  // Observations: we don't seem to get much more than 50% sparsity for average
  // problems with a 10% cutoff wrt the diagonal,
  // therefore sparsifying L doesn't seems worthwhile.....
  private void testSparsityOfRmatrices()
  {
	int i,j,N,sparse_count=0;
	double thresh;
	DoubleMatrix2D L = _MS.mode.Rtransp;

	N = L.rows();
	// direct test on L
	for (i=0;i<N;i++)
	  {
		thresh = 0.1*L.getQuick(i,i);
		for (j=0;j<=i;j++)
		  {
			if (Math.abs(L.getQuick(i,j))<thresh) sparse_count++;
		  }
	  }
	_EH.println("Estimated sparsity = "+
					  (((double)sparse_count)/(N*N)));
	// eigenvalues....
	Algebra alg = new Algebra();
	DoubleMatrix2D Ltr =  L.viewDice();
	DoubleMatrix2D Capp = alg.mult(L,Ltr);
	EigenvalueDecomposition eD = new EigenvalueDecomposition(Capp); 
	DoubleMatrix1D eValues = eD.getRealEigenvalues();
	// eValues are in ascending order, and all +ve
	sparse_count=0;
	thresh=0.1*eValues.getQuick(N-1);
	for (i=0;i<N;i++)
	  {
		if (eValues.getQuick(i)<thresh) sparse_count++;
	  }
	_EH.println("Estimated eigenvalue sparsity = "+
					  (((double)sparse_count)/N));
	alg=null;
	Capp=null;
	eD=null;
	Ltr=null;
	eValues=null;
  }

  // Rough compute of expected uncertainties
  private void debugGuessPosteriorErrorBars()
  {
	int i,N;
	DoubleMatrix2D L = _MS.mode.Rtransp;

	N = L.rows();
	DoubleMatrix2D Capp = (new Algebra()).mult(L,L.viewDice());
	for (i=0;i<N;i++)
	  {
		LayerProperties.Property prop = (LayerProperties.Property)((ArrayList) _MS.map_from.get(i)).get(0); 		
		_EH.printf("%-15s %3d %-10g %-6g%n",_m.completePropertyNameMap.get(prop),i,prop.value,Math.sqrt(Capp.getQuick(i,i)));
	  }
  }



  /**
	 Set up some apparatus for model parameters that will need big
	 steps occasionally, to save them doing long random walks in regions
	 that are truncated. So far we have done this for NG. Possibly
	 something analogous could happen for pinched-out times.
  */
 private void setUpBigStepProposals()
  {
	int i,i_m,j,k,sz=0;
	Integer I;
	LayerProperties.Property prop;
	double mean,stddev,tmp,scale;


	// set up vector of last states for each mode and other useful stuff
	for (k=0;k<_MS.localModes.size();k++)
	  {
		_MS.mode = _MS.localModes.get(k);
		_MS.setMode(_MS.mode);
		DoubleMatrix2D Cinv_prior = _MS.Cinv_prior_only[_MS.mode.fs];

		// E.g. for NG>1, The NG is propelled back to the (0,1) interval only by the
		// prior, and if the posterior fixes NG tightly, such a random walk can be slow, leading to
		// long return excursion of the pre-truncated NG above 1.

		_MS.mode.bigStepInfo = new ArrayList();
		HashMap used_indices=new HashMap();
		for (i_m=0;i_m<_m.n_active_layers;i_m++)
		  {
			i = _MS.master_layer_number_sequence[i_m];
			// (1) NG
			prop = _m.LP[i].NG;	
			if (prop!=null && (I = (Integer)_MS.map_to.get(prop))!=null)
			  {
				// if the property in question is a mirror, skip it (we
				// don't want any duplicates
				if (used_indices.get(I) != null) continue;
				else used_indices.put(I,I);
				int indx = I.intValue();
				// suspected posterior mean of NG
				mean = _MS.mode.QP_optimum.getQuick(indx);
				// approximate measure of posterior scale
				// NG posterior variance - obtain from Cholesky factor...
				DoubleMatrix1D L = _MS.mode.Rtransp.viewRow(indx);
				for (j=0,stddev=0.0; j<=indx; j++) stddev += ((tmp=L.getQuick(j))*tmp);
				stddev = Math.sqrt(stddev);
				// does this posterior look like it will allow NG to leak outside [0,1]?
				if (mean+BIG_JUMP_ACTIVATION_MARGIN*stddev>1 || 
					mean-BIG_JUMP_ACTIVATION_MARGIN*stddev<0)
				  {
					// looks like it could cross into the abyss. Mark it.
					scale = 1.0/Math.sqrt(Cinv_prior.getQuick(indx,indx));
					_MS.mode.bigStepInfo.add(new SortableNumbers(scale,indx));
				  }
			  }
			// (2..) Anything else????
		  }

		// after all the additions to the likely-candidates list, check...
		if ((sz=_MS.mode.bigStepInfo.size())==0) 
		  {
			_MS.mode.bigStepInfo=null;
		  }
		else
		  {
			// rescale using the magic Gelman rule
			scale = GELMAN_NUMERATOR_CONSTANT/Math.sqrt((double)sz);
			for (i=0;i<sz;i++)
			  {
				SortableNumbers bigStepInfo = (SortableNumbers)_MS.mode.bigStepInfo.get(i);
				bigStepInfo.value *= scale;
			  }
		  }
		used_indices = null;
		if (_m._SI.verbosity>2 && sz>0)
		  {
			_EH.println("Mode "+k+": Mixing in random big steps for "+sz+" parameter(s)");
		  }
	  }
  }

  void formHoppingProposalDraw()
  {
	int i,j,jj=0,n,N;
	int fail_count, n_trial=0, n_t_trial=0;
	DoubleMatrix1D w;
	double s=0,scale,tscale,scale2=1.0;
	// DoubleMatrix2D Ltransp;
	MatrixInterface2D Ltransp;
	boolean useBlocks=!true;
	boolean sample_succeeded;
	double[] _property_sample;
	double[] _last_accepted_state;

	_property_sample = _MS.property_sample;
	_last_accepted_state = _MS.mode.last_accepted_state;

	// point covariance cholesky factor R and mean at those of current mode
	_MS.mu = _MS.mode.QP_optimum;


	// Ltransp = _MS.mode.Rtransp;
	
	Ltransp = (MatrixInterface2D) _MS.mode.Rtransp;

	N = _MS.mu.size();       // dimension of full property suite
	n = _m.n_active_layers;         // dimension of times subvector

	// sample from doubly-updated Covariance (isopachs plus local times)
	w = new DenseDoubleMatrix1D(N);
	fail_count=0;

	// set proposal scale using MALA or Gelman rules
	if (useMALAsampler)
	  {
		scale = MALA_NUMERATOR_CONSTANT/Math.pow((double)N,1.0/6.0);
		scale2 = scale*scale/2.0;
	  }
	else
	  {
		scale = GELMAN_NUMERATOR_CONSTANT/Math.sqrt((double)N);
	  }
	

	// If we're doing Gelman convergence tests, fling the first scale to an over-dispersed value
	// so the ensemble starting points are overdispersed
	if (mc_cntr==1 && _m._SI.performGelmanConvergenceDiagnostics)
	  {
		scale=2;
	  }
	if (_m._SI.dumpGassmanDiagnostics) _MS.showGassmanDiagnostics(true,true);

	tscale=scale;

	// some proposals may be unphysical, so we have to attempt this possibly several times.
	for (n_trial=0;;n_trial++)
	  {
		// make a sample
		// Variable 0 is t_top, variable n is t_base
		for (n_t_trial=0;;n_t_trial++)
		  {
			// _MS.Ngen.nextDouble();
			for (i=0;i<=n;i++) w.setQuick(i,_MS.Ngen.nextDouble());
			for (i=0;i<=n;i++)
			  {
				s = 0;
				// From inverse of right cholesky decomposition,
				// Ltransp has true left cholesky factor
				for (j=0;j<=i;j++)
				  {
					// notice transpose!
					s += Ltransp.getQuick(i,j) *  w.getQuick(j);
				  }
				// make proposal by small hop from last state
				_property_sample[i] = _last_accepted_state[i] + s * tscale;
				if (useMALAsampler)
				  {
					_property_sample[i] -= 
					  scale2*(_last_accepted_state[i]-_MS.mode.QP_optimum.getQuick(i));
				  }
				// set the property for this index
				_MS.setValuesFromIndex(i,_property_sample[i]);
			  }
			// Ensure the times are ordered
			_MS.orderTimes();

			// FIXME! 
			// DECISION: Shall we constrain the time samples to mutually exclusive "pinchout" 
			// or "mode-domains" states here???
			// Like, e.g.
			if (forbidLayerCodeChangingJumps)
			  {
				// FIXME PONDER: Voronoi constraints alone are enough I think......
				//if ((_MS.layer_code == _MS.mode.layer_code) && isSampleInVoronoiCellOfMode(_MS.mode)) break;
				// NB this breaks reversibility, so effects are uncertain
				// (correct method is to break here, but hit candidate with huge likelihood penalty....
				if (isSampleInVoronoiCellOfMode(_MS.mode)) break;
			  }
			else
			  {
				// No constraints, so proceed
				break;
			  }
		  }

		// now set remaining properties
		for (i=n+1;i<N;i++) 
		  {
			w.setQuick(i,_MS.Ngen.nextDouble());
		  }
		if (useBlocks)
		  {
			// get submatrix of lower rows
			DoubleMatrix2D Ltransp_part = Ltransp.viewPart(n+1,0,N-n-1 /* (N-1)-(n+1)+1 */,N);
			// form L.w product. scale is buried here
			DoubleMatrix1D scale_Lw = Ltransp_part.zMult(w,null,scale,0,false);	
			// map to properties
			for (i=n+1,j=0;i<N;i++,j++)
			  {
				// make proposal from small deviation:
				_property_sample[i] = _last_accepted_state[i] + scale_Lw.getQuick(j);
				if (useMALAsampler)
				  {
					_property_sample[i] -= 
					  scale2*(_MS.mode.last_accepted_state[i]-_MS.mode.QP_optimum.getQuick(i));
				  }				// set the property for this index
				_MS.setValuesFromIndex(i,_property_sample[i]);
			  }
			Ltransp_part = null; // help the GC
			scale_Lw = null;     // help the GC
		  }
		else
		  {
			// form L.w using fast triangular mult
			Ltransp.lowerTriangularMult(_property_sample,w,n+1,N-1);
			
			// scale and add
			for (i=n+1;i<N;i++)
			  {
				_property_sample[i] = _last_accepted_state[i] + scale *_property_sample[i];
				if (useMALAsampler)
				  {
					_property_sample[i] -= 
					  scale2*(_last_accepted_state[i]-_MS.mode.QP_optimum.getQuick(i));
				  }
				// set the property for this index
				_MS.setValuesFromIndex(i,_property_sample[i]);
			  }
		  }
		// now update full layer properties
		sample_succeeded = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
		if (_m._SI.dumpGassmanDiagnostics) 
		  {
			_MS.showGassmanDiagnostics(sample_succeeded,false);
		  }
		// was it unphysical??
		if (!sample_succeeded) 
		  {
			fail_count++;
		  }
		else
		  {
			fail_count=0; // All is well...
		  }

		// a failure shouldn't happen too many times or things are badly suspect
		if (fail_count>_MS.MAX_MC_FAIL_COUNT)
		  {
			_EH.println("Impossible fluid or rock properties in fluid state " +
							  _MS.fs + ". Regression curves, isopachs, seismic, or "
							  + "fluid properties are probably suspect.");
			_MS.showDiagnosticsOfEffectivePropertiesForAllLayers();
			throw new RuntimeException();
		  }
		if (sample_succeeded) 
		  {
			break;
		  }
	  }
	// MALA samplers need to keep track of this factor
	if (useMALAsampler)
	  {
		DoubleMatrix1D m_prop = new DenseDoubleMatrix1D(N);
		// form m-mbar
		for(i=0;i<N;i++) m_prop.setQuick(i,_last_accepted_state[i] - _MS.mode.QP_optimum.getQuick(i));
		// then
		double fac = (1.0-scale2);
		double coeff = (1.0-fac*fac)/scale;
		DoubleMatrix1D wprime = _MS.mode.RtranspInv.zMult(m_prop,null,coeff,0,false);
		for(i=0;i<N;i++) wprime.setQuick(i,wprime.getQuick(i)-fac*w.getQuick(i));
		q_MALA = (w.zDotProduct(w) - wprime.zDotProduct(wprime))/2.0;
		wprime = null;
		m_prop = null;
	  }
	w=null; // help the GC
  }

  private void formBigStepProposalDraw()
  {
	double RANDOM_PRIOR_BASED_BIG_PROB = 0.5;
	
	if ((_MS.mode.bigJumpDirections == null || 
		 (_MS.Ugen.nextDouble() <= RANDOM_PRIOR_BASED_BIG_PROB)))
	  {
		// prior based big-steps if there are no suspicious eigenmodes
		formBigStepPriorBasedDraw();
	  }
	else
	  {
		formBigStepEigenmodeBasedDraw();
	  }
  }

  private void formBigStepPriorBasedDraw()
  {
	int i,j,indx,n,Nbig,N;
	int fail_count, n_trial=0, n_t_trial=0;
	DoubleMatrix1D w;
	double s;
	DoubleMatrix2D Cinv_prior;
	boolean sample_succeeded;
	double[] _property_sample;

	_property_sample = _MS.property_sample;

	// Get inverse prior Covariance for the fluid state of this mode
	Cinv_prior = _MS.Cinv_prior_only[_MS.mode.fs];

	N = _MS.mode.QP_optimum.size(); // size of current parameter vector
	n = _m.n_active_layers;                // dimension of times subvector
	Nbig = _MS.mode.bigStepInfo.size();

	// sample from doubly-updated Covariance (isopachs plus local times)
	// w = new DenseDoubleMatrix1D(N);  // OLD
	w = new DenseDoubleMatrix1D(Nbig);
	fail_count=0;
	// some proposals may be unphysical, so we have to attempt this possibly several times.
	for (n_trial=0;;n_trial++)
	  {
		// make a sample
		// Variable 0 is t_top, variable n is t_base
		if (false)   // may need to do this branch if time big steps allowed
		  {
			for (n_t_trial=0;;n_t_trial++)
			  {
				for (i=0;i<=n;i++) w.setQuick(i,0.0);
				for (i=0;i<=n;i++)
				  {
					s = 0;
					// From inverse of right cholesky decomposition,
					// Cinv_prior has true left cholesky factor
					for (j=0;j<=i;j++)
					  {
						// notice transpose!
						s += Cinv_prior.getQuick(i,j) *  w.getQuick(j);
					  }
					// make proposal by small hop from last state
					_property_sample[i] = _MS.mode.last_accepted_state[i] + s;
					// set the property for this index
					_MS.setValuesFromIndex(i,_property_sample[i]);
				  }
				// Ensure the times are ordered
				_MS.orderTimes();
				if (true) break;
			  }
		  }
		else
		  {
			// map original properties exactly
			for (i=0;i<N;i++) 
			  {
				_property_sample[i] = _MS.mode.last_accepted_state[i];
				// set the property for this index
				_MS.setValuesFromIndex(i,_property_sample[i]);
			  }
			// Ensure the times are ordered
			_MS.orderTimes();
		  }
		// now overwrite properties that need big jumps
		for (i=0;i<Nbig;i++)
		  {
			w.setQuick(i,_MS.Ngen.nextDouble());
			SortableNumbers bigStepInfo = (SortableNumbers)_MS.mode.bigStepInfo.get(i);
			indx = bigStepInfo.index;
			// get hop scale from prior inverse covariance entry....(see setUpBigStepProposals())
			s = bigStepInfo.value * w.getQuick(i);
			// make proposal from deviation:
			_property_sample[indx] = _MS.mode.last_accepted_state[indx] + s;
			// set the property for this index
			_MS.setValuesFromIndex(indx,_property_sample[indx]);
		  }
		// now update full layer properties
		sample_succeeded=_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
		if (_m._SI.dumpGassmanDiagnostics) 
		  {
			_MS.showGassmanDiagnostics(sample_succeeded,false);
		  }
		// was it unphysical??
		if (!sample_succeeded) 
		  {
			fail_count++;
		  }
		else
		  {
			fail_count=0; // All is well...
		  }

		// a failure shouldn't happen too many times or things are badly suspect
		if (fail_count>_MS.MAX_MC_FAIL_COUNT)
		  {
			_EH.println("Impossible fluid or rock properties in fluid state " +
							  _MS.fs + ". Regression curves, isopachs, seismic, or "
							  + "fluid properties are probably suspect.");
			_MS.showDiagnosticsOfEffectivePropertiesForAllLayers(); // FC_OK	
			throw new RuntimeException();
		  }
		if (sample_succeeded) 
		  {
			break;
		  }
	  }
	w=null; // help the GC
  }

  private void formBigStepEigenmodeBasedDraw()
  {
	int i,N,nDirections,iDir;
	int fail_count, n_trial=0;
	double w;
	DoubleMatrix2D Cinv_prior;
	boolean sample_succeeded;
	DoubleMatrix1D dir;

	N = _MS.mode.QP_optimum.size(); // size of current parameter vector
	
	// choose direction uniformly
	nDirections = _MS.mode.bigJumpDirections.size();
	iDir = Math.min(nDirections-1, (int)(_MS.Ugen.nextDouble()*nDirections));
	dir = (DoubleMatrix1D)_MS.mode.bigJumpDirections.get(iDir);

	fail_count=0;
	// some proposals may be unphysical, so we have to attempt this possibly several times.
	for (n_trial=0;;n_trial++)
	  {
		// make a sample
		w = 4*_MS.Ngen.nextDouble();
		for (i=0;i<N;i++)
		  {
			// make sample from hop in line direction of random length
			_MS.property_sample[i] = _MS.mode.last_accepted_state[i] + dir.getQuick(i)*w;
			// set the property for this index
			_MS.setValuesFromIndex(i,_MS.property_sample[i]);
		  }
		// Ensure the times are ordered
		_MS.orderTimes();
		// now update full layer properties
		sample_succeeded=_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
		if (_m._SI.dumpGassmanDiagnostics) 
		  {
			_MS.showGassmanDiagnostics(sample_succeeded,false);
		  }
		// was it unphysical??
		if (!sample_succeeded) 
		  {
			fail_count++;
		  }
		else
		  {
			fail_count=0; // All is well...
		  }
		// a failure shouldn't happen too many times or things are badly suspect
		if (fail_count>_MS.MAX_MC_FAIL_COUNT)
		  {
			_EH.println("Impossible fluid or rock properties in fluid state " +
							  _MS.fs + ". Regression curves, isopachs, seismic, or "
							  + "fluid properties are probably suspect.");
			_MS.showDiagnosticsOfEffectivePropertiesForAllLayers();
			throw new RuntimeException();
		  }
		if (sample_succeeded) 
		  {
			break;
		  }
	  }
  }

  /**
	 Deal with proposed transition from last_accepted_mode to _MS.mode
	 via invertable mappings
  */
  private boolean formMappedModeJumpProposalDraw()
  {
	int i;
	boolean sample_succeeded;

	negLogQtheta=0.0;
	negLogDetCphi=0.0;

	modeHopIsPossible=true;

	if (_MS.mode == last_accepted_mode)
	  {
		// Restore current 'snapshot variables' according to last state for this mode
		ArrayUtils.copy(_MS.mode.last_accepted_state,_MS.property_sample);
		ArrayUtils.copy(_MS.mode.last_accepted_R_pp,_m.SM.R_pp);
		//ArrayUtils.copy(_MS.mode.last_accepted_R_far,_m.SM.R_far);
		if (_m._SI.doSimpleFluidUnitsAndContacts) 
		  {
			ArrayUtils.copy(_MS.mode.last_accepted_R_pp_fc,_m.SM.R_pp_fc);
			//ArrayUtils.copy(_MS.mode.last_accepted_R_far_fc,_m.SM.R_far_fc);
			ArrayUtils.copy(_MS.mode.last_accepted_t_fc,_m.SM.t_fc);
		  }
		for (i=0;i<_MS.mode.last_accepted_state.length ; i++)
		  {
			_MS.setValuesFromIndex(i,_MS.property_sample[i]);
		  }
		// Ensure the times are ordered
		_MS.orderTimes();
	  }
	else
	  {
		// Preserve nested parameters insofar as is possible
		mapModelsViaCommonElements(last_accepted_mode, _MS.mode); 
	  }

	// now update full layer properties
	sample_succeeded=_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
	// was it unphysical??
	if (!sample_succeeded) 	return false;

	// was the hop allowed (are mode integrities preserved)?
	if (!modeHopIsPossible) return false;
	return true;
  }






  private double computeNewStateImportanceRatio(boolean modeJumpsMustLeaveOldMode, boolean useStrictPosteriorWeights)
  {
	int N_old,N_new;
	double negLogLik_ratio;
	double negLogLik_full;
	double u,jacobian;
	double logDetRatio=0;

	// recompute the new likelihood using the variables for the new state.
	// necessary only on initialisation or if the state changes


	// mapping via implied-normal-deviates matching.....
	//jacobian = (_MS.mode.logDetC - last_accepted_mode.logDetCtilde)/2.0;

	// or, via direct mapping
	jacobian = 0.0;


	negLogLik_full = (negLogLikSeisIsopachDepthVS + negLogLikPriorOnly);

  
	negLogLik_ratio = negLogLik_full
	                - (last_accepted_mode.lastAcceptedNegLogAllLikelihoods
					   + last_accepted_mode.lastAcceptedNegLogLikPriorOnly);

	// now compensate for mode-jump factors
	N_old = last_accepted_mode.last_accepted_state.length;  // dimension of current model
	N_new = _MS.mode.last_accepted_state.length;            // dimension of proposed model

	// (1) jacobian
	// jacobian is zero, but put this factor in as a genuflection to the formula.
	negLogLik_ratio -= jacobian;

	// (2) dimension-matching theta parameters
	// (2a) look after determinant issues
	// Ratio of determinants in MH factor (P(new)/P(old)) . 1/q(phi) as they would appear
	// so this will give a { (|C_old|/|C_new| ). |C_q| }^(1/2) type term.
	logDetRatio = last_accepted_mode.logDetCprior - _MS.mode.logDetCprior - negLogDetCphi;
	logDetRatio/=2.0;

	// (2b) then the q(phi) issues for generating missing parameters etc
	// _EH.println("DET(+)="+logDetRatio + " detCphi="+negLogDetCphi);
	negLogLik_ratio -= (logDetRatio + negLogQtheta);


	// Compensate for bias possibly induced by biased mode proposal
	// We use the either the posterior weight or the hybrid weight, depending on what 
	// has been used in ModelSampler.proposeNewMode()
	if (useStrictPosteriorWeights)
	  {
		u = last_accepted_mode.weight/_MS.mode.weight;
		if (modeJumpsMustLeaveOldMode && (_MS.mode.weight<1.0))
		  {
			u *= (1.0-last_accepted_mode.weight)/(1.0-_MS.mode.weight);
		  }
	  }
	else
	  {
		u = last_accepted_mode.hybrid_weight/_MS.mode.hybrid_weight;
		if (modeJumpsMustLeaveOldMode && (_MS.mode.hybrid_weight<1.0))
		  {
			u *= (1.0-last_accepted_mode.hybrid_weight)/(1.0-_MS.mode.hybrid_weight);
		  }
	  }

	// Ratio of terms from prior weights
	// Old 
	// u *= (_MS.mode.prior_weight/last_accepted_mode.prior_weight);
	// new: this term comes from strict prior, not proposal mechanism
	// usually, mode.prior_weight = mode.FS.prob, but not if multimodes per fluid state
	// means the FS.prob has been farmed out over various modes.
	u *= (_MS.mode.FS.prob/last_accepted_mode.FS.prob);



	negLogLik_ratio -= Math.log(u);


	return negLogLik_ratio;
  }

  private double computeHoppingLikelihood(boolean allowOutput)
  {
	int i;
	double w_sum_sq=0.0, negLogLik;
	double negLogLik_iso, negLogLik_depth, negLogLik_seis, negLogVStackLikelihood;
	double negLogLik_PStime; //JChen

	// (1) SEISMIC mismatch term
	negLogLik_seis = _MS.negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
	_MS.negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
	negLogLik = negLogLik_seis;

	
	
	// (2) Add in logLikelihood contribution from ISOPACHS if they're required
	negLogLik_iso = _MS.computeIsopachNegLogLikelihood(false);
	negLogLik += negLogLik_iso;

	// (3) Add in logLikelihood contribution from DEPTH self-consistency if required
	negLogLik_depth= _MS.computeConsistentDepthNegLogLikelihood(false); 
	negLogLik += negLogLik_depth; 

	// Add in logLikelihood contribution from vertical stacks
	if (_m._SI.doStackedModel) 
	  {
		negLogVStackLikelihood = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
		negLogLik += negLogVStackLikelihood;
	  }

	// 4D ensure variable concoction is 4D-plausible
	negLogLik += _MS.compute4DPenaltyLikelihood();

	// Added by JChen : JG check this what about 4D????// FIXME// FIXME// FIXME// FIXME// FIXME	
	// (4) Add in loglikelihood contribution from PS EVENT REGISTRATION TIME
	//negLogLik_PStime = _m.SM.computePSTimeLikelihood();
	//negLogLik += negLogLik_PStime; 

	// Sum contribution of all above
	negLogLikSeisIsopachDepthVS = negLogLik;

	// (5) Also the contribution from the raw PRIOR
	// (i.e. before isopachs were used in the linearised update)
	negLogLikPriorOnly = _MS.computePriorNegLogLikelihood();

	negLogLik += negLogLikPriorOnly;


	// (6) Output
	// This output appears to be a cpu eater, so only switch it on if absolutely required
	if (allowOutput && _m._SI.verbosity>3) 
	  {
		_EH.printf("%2d %2d      | ",_MS.fs, _MS.localModes.indexOf(_MS.mode));
		for (i=0;i<_m.n_layers;i++) 
		  {
			_EH.printf("%1d ",_MS.layer_present[i] ? 1 : 0);
		  }
		_EH.printf("%5d ",_MS.layer_code);
		_EH.printf("| %10.6f %10.6f %10.6f | ",negLogLik, negLogLikPriorOnly, negLogLikSeisIsopachDepthVS); 
		//_EH.printf("| %10.6f %10.6f %10.6f %10.6f %10.6f| ",negLogLik, negLogLikPriorOnly, negLogLikSeisIsopachDepthVS,negLogLik_seis,negLogLik_iso); //
		
	  }


	return negLogLik;
  }
  

  /** 
	  A method to jump from an oldMode to a newMode, mapping common
	  elements in the respective model vectors directly where
	  possible. Parameter i in the old model is mapped to j in the
	  new. Non-common elements have to be simulated from a proposal;
	  at first pass we use the prior. Non-common elements are often
	  fluid densities etc. Follows treatment in Andrieu's signal
	  processing paper
  */
  private void mapModelsViaCommonElements(LocalMode oldMode,LocalMode newMode)
  {
	int i,j,N_new,N_old;
	ArrayList vNew,vOld;
	LayerProperties.Property propOld,propNew; 
	double mode_offset;
	double w,covar;
	double detPhi=0;
	double negLogLikQtheta=0.0;
	DoubleMatrix2D R;
	boolean doTranslation;
	boolean doScaledJumps=false; // Retain?????-------> Doesn't seem to help when isopachs are tight
	double scale,Jac=0.0;        // Retain?????
	double twice_midpoint;

	// set up vector of last states for each mode
	N_new = newMode.QP_optimum.size();
	N_old = oldMode.QP_optimum.size();

	Integer targetIndex;

	// for mapping, set fluid state to new mode
	// make sure mode is correct.......
	_MS.setMode(newMode);

	// useful pointers
	cern.colt.map.OpenIntObjectHashMap map_from_in_old_state = _MS.map_from_index[oldMode.fs];
	HashMap map_to_in_new_state = _MS.map_to_index[newMode.fs];

	// choice of Cholesky factor from which to draw non-nested variable deviates
	// R = oldMode.R;
	// or, could have, R_prior_only[fs] is an upper triangular Cholesky factor of Cinv = C^{-1}
	R = _MS.R_prior_only[oldMode.fs];
	
	doTranslation = true; 

	// loop over parameters of old mode
	for (i=0;i<N_old;i++)
	  {
		// find out if this parameter is in the new (target) mode
		vOld = (ArrayList)map_from_in_old_state.get(i); 

		propOld = (LayerProperties.Property) vOld.get(0);
		// is propOld in the target model?
		if ((targetIndex = (Integer)map_to_in_new_state.get(propOld)) != null)
		  {
			// yes it is. Get its index
			j = targetIndex.intValue();
			//_EH.println("COMMON "+i+" "+j);
			// Now map the property old.i to new.j relative to mode locations...?
			mode_offset = newMode.QP_optimum.getQuick(j);
			if (doScaledJumps)
			  {
				scale = newMode.R.getQuick(j,j)/oldMode.R.getQuick(i,i); // FIXME retain?????
				// stretching Jacobian
				Jac += Math.log(scale);  // FIXME retain?????
				// scaled..
				_MS.property_sample[j] = scale * (oldMode.last_accepted_state[i] - oldMode.QP_optimum.getQuick(i)) + mode_offset;
			  }
			else
			  {
				if (doTranslation)
				  {
					// simple mode-peak offset mapping. Jacobian is unity.
					_MS.property_sample[j] = (oldMode.last_accepted_state[i] - oldMode.QP_optimum.getQuick(i)) + mode_offset;
				  }
				else
				  {
					// centre-point reflection mapping. Also unit Jacobian
					twice_midpoint = (oldMode.QP_optimum.getQuick(i) + mode_offset);
					_MS.property_sample[j] = -oldMode.last_accepted_state[i] + twice_midpoint;
				  }
			  }
			// or vanilla - flavour ones (no offsets) ?
			//_MS.property_sample[j] = oldMode.last_accepted_state[i];
			_MS.setValuesFromIndex(j,_MS.property_sample[j]);
			// stretching Jacobian
			Jac += (doScaledJumps) ? Math.log(newMode.R.getQuick(j,j)/oldMode.R.getQuick(i,i)) : 0.0;
		  }
		else
		  {
			//_EH.println("SOURCE ONLY "+i);
			// Here we have an unpaired variable in the source model. We will map it as an independent variable
			// whose covariance etc is approximated as the corresponding diagonal entry of R^-2. We could
			// use either the posterior or prior R..
			// This term ends up on the numerator, hence the -= for w^2.
			covar = 1.0/R.getQuick(i,i);
			// what the jump must have been to generate these source variables
			// .....I'm modestly sure this is right.......
			w = (oldMode.last_accepted_state[i] - oldMode.QP_optimum.getQuick(i))/covar;
			negLogLikQtheta -= w*w;
			// acquire dimension-matching variable covariance.
			detPhi -= Math.log(covar);   // FIXME can precalculate
		  }
	  }



	// now find any variables in the target that are not present in the source
	// useful pointers
	HashMap map_to_in_old_state = _MS.map_to_index[oldMode.fs];
	cern.colt.map.OpenIntObjectHashMap map_from_in_new_state = _MS.map_from_index[newMode.fs];

	R = _MS.R_prior_only[newMode.fs];

	// loop over variables in new model
	for (j=0;j<N_new;j++)
	  {
		vNew = (ArrayList)_MS.map_from.get(j);
		propNew = (LayerProperties.Property) vNew.get(0);
		// is propNew not in the source model?
		if ((targetIndex = (Integer)map_to_in_old_state.get(propNew)) == null)
		  {
			//_EH.println("TARGET ONLY "+j);
			// Here we have an unpaired variable in the target model. We will map it as an independent variable
			// whose covariance etc is approximated as the corresponding diagonal entry of R squared. We could
			// use either the posterior or prior R.
			
			covar = 1.0/R.getQuick(j,j);
			mode_offset = newMode.QP_optimum.getQuick(j);
			w = _MS.Ngen.nextDouble();
			_MS.property_sample[j] = covar * w + mode_offset;
			_MS.setValuesFromIndex(j,_MS.property_sample[j]);
			negLogLikQtheta += w*w;
			// acquire dimension-matching variable covariance.
			detPhi += Math.log(covar);   // FIXME can precalculate
		  }
	  }
	//_EH.println("stretch-jac = "+ Math.exp(Jac));
	// do final 1/2 factor of neg-log-likelihood of dimension-matching variables
	negLogLikQtheta/=2.0;
	// and factor of 2 for phi-variables' covariance
	detPhi *=2.0; 
	// convert to a negative log-likelihood 
	negLogDetCphi = - detPhi;

	// make sure mode is correct.......
	_MS.setMode(newMode);

	// Ensure the times are ordered
	_MS.orderTimes();
	negLogQtheta = negLogLikQtheta;

	// Test if mapping messes up layer_code or mode-membership (defined via Voronoi cells)
	// Mode boundaries cannot be crossed in the same fluid state if forbidLayerCodeChangingJumps == true
	// and act on this if layer_code preservation required
	// Note reversibility issues attached to this method! Best switched off!
	if (forbidLayerCodeChangingJumps && newMode.fs==oldMode.fs )
	  {
		//if (_MS.layer_code != newMode.layer_code || !isSampleInVoronoiCellOfMode(newMode))
		if (!isSampleInVoronoiCellOfMode(newMode))
		  {
			// FIXME PONDER: Voronoi constraints alone are enough I think......
			// We either changed the layer configuration, or
			// jumped into the Voronoi cell of another mode. Sampling of this part of
			// space ought to happen within that other mode.....
			modeHopIsPossible=false;
		  }
	  }
  }

  /**
	 Determine whether current point in _MS.property_sample belongs
	 to the Voronoi cell of newMode, where the cell is defined in the
	 time subspace only, and the cell centres are the points mode.QP_optimum
	 for all the available modes in the array _MS.localModes. The
	 set of comparison modes is restricted to those of the same fluid state.
  */
  private boolean isSampleInVoronoiCellOfMode(LocalMode newMode)
  {
	LocalMode nearestMode,test_mode;
	double sum,min_sum;
	DoubleMatrix1D mode_centre;
	
	min_sum = Double.MAX_VALUE;
	nearestMode = newMode;
	for (int k=0;k<_MS.localModes.size();k++)
	  {
		test_mode = _MS.localModes.get(k);
		if (test_mode.fs != newMode.fs) continue;  // skip modes with different fluid states
		mode_centre = test_mode.QP_optimum;	
		sum=0.0;
		for (int i=0;i<=_m.n_active_layers;i++)
		  {
			sum += MathUtil.sqr(_MS.property_sample[i]-mode_centre.getQuick(i));
		  }
		if (sum<min_sum)
		  {
			min_sum = sum;
			nearestMode = test_mode;
		  }
	  }
	if (nearestMode == newMode) 
	  {
		return true;
	  }
	else 
	  {
		return false;
	  }
  }



  // test routine only......

  private double experimental_examineNonMappingVariableCovarianceDeterminant(LocalMode newMode, ArrayList indices)
  {
	int i,j,ii,jj,N;
	DoubleMatrix2D C;
	LUDecompositionQuick LUD;

	N=indices.size();
	C = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(N,N);
	
	for (i=0;i<N;i++)
	  {
		ii = ((Integer)indices.get(i)).intValue();
		_EH.println("Unmappped indices: "+ii); 
		for (j=0;j<N;j++)
		  {
			jj = ((Integer)indices.get(j)).intValue();
			C.setQuick(i,j,newMode.R.getQuick(ii,jj));
		  }
	  }
	LUD = new LUDecompositionQuick();
	LUD.decompose(C);
	double detR = LUD.det();
	// return determinant of implied covariance
	_EH.println("DETRATIO (cnr): "+ detR*detR);
	//au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("Rsub."+mc_cntr+".txt",C);
	return 2.0*Math.log(detR);
  }

}
