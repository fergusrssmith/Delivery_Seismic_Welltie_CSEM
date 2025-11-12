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



import cern.colt.matrix.*;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.util.*;
import cern.jet.random.ChiSquare;

import au.csiro.JamesGunning.lib.util.MathUtil;

/**
	 A class to draw samples of the model using the MCMC technique, with a Metropolis-Hastings 
	 algorithm based on full Multi--Gaussian approximations to the local modes. See
	 Independence samplers as discussed in the Reports.
*/
public final class IndependenceSampler
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
	 Independence samples must match the layer_code configuration:
	 This flag sets this behaviour. The matching is not insisted upon if
	 seismic is ignored, as there will not be a posterior surface
	 discontinuity in this case
  */
  public boolean MATCH_LAYER_CODES = true;

  /**
	 Independence samples must match the layer_code configuration:
	 if we fail to get the layer code right this many times in a row, 
	 the mode's pinchout configuration is too narrow
	 and inaccessible. The algorithm implemented here is then not 
	 appropriate
  */
  public static final int MAX_MATCH_LAYER_CODE_ATTEMPT_COUNT=1000;

  /**
	 overdisperser flag for proposals
  */
  private boolean useStudentOverDispersion=false;
  /**
	 whether to use Approximate Mixed Gaussian Sampling
  */
  private boolean useApproximateMixedGaussianSampling=false;
  private boolean insistOnMatchedLayerCode;

  private ChiSquare chiSqareGenerator;
  private double studentDOF;

  protected IndependenceSampler(Model m, boolean useApproximateMixedGaussianSampling)
  {
	this._m = m;
	this._MS = m.MS;
	this.useApproximateMixedGaussianSampling = useApproximateMixedGaussianSampling;
	this.insistOnMatchedLayerCode = (_m._SI.ignoreSeismic ? false : MATCH_LAYER_CODES);

	// Over dispersion seems problematic in these models. We require
	// DOF >> d for the overdispersion to be `mild', and this doesn't
	// allow much to happen once d>50 or so. I.e. experiments show this
	// doesn't do anything useful, so leave it switched off
	useStudentOverDispersion = false;
	if (useStudentOverDispersion)
	  {
		studentDOF = 300;
		chiSqareGenerator = new ChiSquare(studentDOF, _MS.ranEngine);
	  }
  }



  /**
	 Draw samples of the model using the MCMC technique, with a Metropolis Hastings algorithm
  */
  void sampleFromLocalModesUsingFullLocalModes()
  {
	int i,k,mc_cntr,realisation_number=0,realisation_label_number=0;
	int decimation_constant=1;
	boolean write_state,newState = false;

	int fs_state=0, N_accepted=1;
	int[] N_fs;
	int k_mode_indx=0, kp_mode_indx=0;
	// old state k terms. m = model (old state)
	double P_k_m=0.0, L_k_m=0.0, Papp_k_m=0.0;
	// state k' terms  mp = model m' (new state)
	double P_kp_mp=0, P_kp_mp_sansdet=0, L_kp_mp=0, Papp_kp_mp=0;
	double Pi_new, Pi_old;   // full posterior densities of old and new states
	double JoldToNew;        // proposal density from old to New state
	double JnewToOld;        // the converse
	double[] pk;             // mode prior logProb
	double[] qk;             // mode proposal density
	double numerator,denominator;
	boolean modeJumpPossible,allowTraceOutput=true;

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

	// Check all modes are accessible
	eliminateFaultyModes();


	// loop over number of samples in MCMC chain
	int Nsamples = 2 * _m._SI.n_samples_per_trace_undecimated*decimation_constant;

	// set up initial array of p_k entries (prior mode weights) and qk's (mode proposal probabilities)
	pk = new double[_MS.localModes.size()];
	qk = new double[_MS.localModes.size()];
	for (k=0;k<_MS.localModes.size();k++)
	  {
		_MS.mode = _MS.localModes.get(k);
		pk[k] = Math.log(_MS.mode.prior_weight);
		// Mode proposal density is prior if using rigorous sampler, or
		// approximate posterior weigth if using useApproximateMixedGaussianSampling
		if (useApproximateMixedGaussianSampling)
		  {
			qk[k] = Math.log(_MS.mode.weight);
		  }
		else
		  {
			qk[k]=pk[k];  // prior
		  }
	  }
	if (_m._SI.verbosity>3) 
	  {
		_m._SI.EH.println("Sample mode fluid-state | layers present, code | prior-logLik iso+seismic-logLik full-logLik | "+
						  "proposal_prob success%n");
	  }
	
	// fluid state prob cntrs
	N_fs = new int[_m.FS.possibleFluidStates.size()];

	for (mc_cntr=1; mc_cntr<=Nsamples; mc_cntr++)
	  {
		// Compute sample number
		realisation_number = mc_cntr-Nsamples/2;
		// flag to allow mode-jumps as well as model jumps. (Must be state-independent/deterministic)
		// Currently alternate even odd steps....
		//	modeJumpPossible = ((mc_cntr & 1)==1);
		modeJumpPossible = true;

		// flag to signify new transition in MCMC chain
		newState = false;

		// Form new candidate
		// Propose state k'
		// alternate possible mode-jumpers with steps within a mode
		if (modeJumpPossible || mc_cntr==1)
		  {
			if (useApproximateMixedGaussianSampling)
			  {
				kp_mode_indx = _MS.proposeModeFromApproximatePosterior();
			  }
			else
			  {
				kp_mode_indx = _MS.proposeModeFromPrior();
			  }
		  }
		else kp_mode_indx = k_mode_indx;

		// Form proposal from posterior approximation and get its proposal prob.
		Papp_kp_mp = formProposalDrawAndGetNegLogPapp(true);
		// must augment prior with determinant terms. (2pi)^d/2 terms cancel.
		Papp_kp_mp += _MS.mode.logDetCtilde/2.0;
		// The prior prob. of this proposal (chisq bit)
		P_kp_mp_sansdet = computeProposalPriorNegLogProb();
		// must augment prior with determinant terms. (2pi)^d/2 terms cancel.
		P_kp_mp = P_kp_mp_sansdet + _MS.logDetC_prior_only[_MS.fs]/2.0;
		// likelihood
		L_kp_mp = computeProposalNegLogLikelihood();

		// initialise 'old' state k
		if (mc_cntr==1) 
		  {
			k_mode_indx = kp_mode_indx;
			L_k_m = L_kp_mp;
			P_k_m = P_kp_mp;
			Papp_k_m = Papp_kp_mp;
		  }
		// (a) New state m'_{k'} cf old state m_k

		// Apply Metropolis-Hastings rule 
		Pi_new = pk[kp_mode_indx] - L_kp_mp - P_kp_mp;
		JoldToNew = ( modeJumpPossible? qk[kp_mode_indx] : 0.0) - Papp_kp_mp;  
		// (a) note the qk[k_mode_indx] will usually cancel with JnewToOld, but we put it here
		// to be more explicit about the proposal density
		Pi_old = pk[k_mode_indx] - L_k_m - P_k_m;
		JnewToOld = ( modeJumpPossible? qk[k_mode_indx] : 0.0) - Papp_k_m;  
		// see comment (a) above
		
		// the fraction alpha = exp(numerator-denominator), is that which appears in the MH acceptance rule;
		numerator = Pi_new + JnewToOld;
		denominator = Pi_old + JoldToNew;

		// some feedback
		if (_m._SI.verbosity>3)  
		  {
			_m._SI.EH.printf("%4d %3d %2d  | ",mc_cntr, kp_mode_indx, _MS.fs);
			for (i=0;i<_m.n_layers;i++) 
			  {
				_m._SI.EH.printf("%1d ",(_MS.layer_present[i] ? 1 : 0));
			  }
			_m._SI.EH.printf("%5d ",_MS.layer_code);
			_m._SI.EH.printf("| %10.6f %10.6f %10.6f | %10.6f %10.6f ",
							P_kp_mp, L_kp_mp, P_kp_mp + L_kp_mp,
							Papp_kp_mp, numerator - denominator);
		  }


		// Apply acceptance rule: accept with probability min(1,alpha)
		// make comparison with care to avoid overflow in wierd cases
		if (useApproximateMixedGaussianSampling ||
			(numerator >= denominator) || _MS.Ugen.nextDouble() <= Math.exp(numerator - denominator))
		  {
			// proposal accepted. Update everything
			if (_m._SI.verbosity>3) 
			  {
				_m._SI.EH.print("A%n");
			  }
			fs_state=_MS.fs;
			newState = true;

			L_k_m = L_kp_mp;
			P_k_m = P_kp_mp;
			Papp_k_m = Papp_kp_mp;
			k_mode_indx = kp_mode_indx;
			
			_MS.mode.lastAcceptedNegLogAllLikelihoods = L_kp_mp;
			_MS.mode.lastAcceptedNegLogLikPriorOnly = P_kp_mp_sansdet;

			if (realisation_number>0) N_accepted++;
		  }
		else
		  {
			if (_m._SI.verbosity>3) _m._SI.EH.print("R%n");
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
			N_fs[fs_state]++;
			// ensure first member of second half of chain is set and dumped out
			boolean dumpState = (realisation_number==1)? true : newState;
			// dump model
			_m._SI.dumpRealisation(dumpState,realisation_label_number,write_state && allowTraceOutput,false);
			// FIXME REMOVE NEXT LINE
			// Util.printfileZO("error."+realisation_number+".txt", _m.SM.seismic_error_vector,_m.SM.seismic_error_vector.length);
			// and seismics if required
			_m.SM.dumpSynthetics(dumpState,realisation_label_number,write_state && allowTraceOutput);
			// if doing diagnostics, get scalar estimands
			if (_m._SI.performGelmanConvergenceDiagnostics ||
				_m._SI.performRafteryConvergenceDiagnostics) 
			  {
				_MS.aquireScalarEstimands(realisation_label_number-1);
			  }
		  }
	  }	 
	if (_m._SI.verbosity>2) 
	  {
		double norm =(double)((realisation_number>0)?realisation_number:1);
		_m._SI.EH.println("Acceptance rate = "+Math.min(1.0,(double)N_accepted/norm));
		for (k=0;k<_m.FS.possibleFluidStates.size();k++)
		  {
			_m._SI.EH.printf("Prob of state fs=%d is %6g%n",k,(N_fs[k]/norm));
		  }
	  }
	pk=null; // GC
	qk=null; // GC
  }

  /**
	 Some modes occasionally won't be able to be generated from the independence sampler.
	 the signature of this is that formProposalDrawAndGetNegLogPapp(false) == Double.MAX_VALUE
	 (occurs after many abortive attempts to get the right pinchout configuration)
	 Here we find these and eliminate them
  */

  private void eliminateFaultyModes()
  {
	int k,j=0;
	LocalMode mode;
	boolean mustRenormalise=false;

	// loop over current modes
	for (k=0;k<_MS.localModes.size();k++)
	  {
		mode = _MS.localModes.get(k);
		_MS.setMode(mode);
		// check if mode is faulty
		if (formProposalDrawAndGetNegLogPapp(false) == Double.MAX_VALUE)
		  {
			// if so eliminate it...
			j++;
			_MS.localModes.remove(k--);	  
			mustRenormalise=true;
		  }
	  }
	if (mustRenormalise) 
	  {
		if (_m._SI.verbosity>2) 
		  {
			_m._SI.EH.println(j+ " spurious modes for independence sampler found & eliminated");
		  }
		_MS.localModes.trimToSize();
		_MS.reNormaliseLocalModes(false);
	  }
  }


  private double updateScaledlogLikelihoodArray(int k, double negLogLik, double[] llk, double Lmin)
  {
	if (Lmin <= negLogLik)
	  { 
		// min doesn't change: update of only one number in scaled-array needed
		llk[k] = Math.exp(Lmin - negLogLik);
	  }
	else 
	  {
		// Lmin has to change, rescaling of whole array required
		double tmp = Math.exp(negLogLik - Lmin);
		Lmin = negLogLik;
		for (int i=0; i<_MS.localModes.size(); i++) llk[i] *= tmp;
		// this is the new element setting the scaling
		llk[k]=1.0;
	  }
	return Lmin;
  }

  /**
	 Generate sample from local mode. Sometimes the modes are badly configured
	 and we can't get realisations to satisfy the pinchout criterion. If so, throw
	 and exception after MAX_MATCH_LAYER_CODE_ATTEMPT_COUNT failed attempts, if
	 throwExceptionOnFaultyMode=true. Otherwise, return Double.MAX_VALUE. 

	 If the draw is successful, will return the usual negative Log proposal probability
  */
  protected double formProposalDrawAndGetNegLogPapp(boolean throwExceptionOnFaultyMode)
  {
	int i,j,jj=0,n,N;
	int fail_count;
	double[] w;
	double s;
	DoubleMatrix2D Ltransp;
	boolean sample_succeeded;

	// Point covariance cholesky factor R and mean at those of current mode
	// NOTE! We have two choices for the mean here:
	// (1) as mean, use the extrapolated quadratic approximation mean, which may not generate
	// a valid pinchout configuration at this point. A large fraction of realisations
	// may however be valid (the mode determination should have tested this requirement)
	// Caveat - the quadratic_mean may be a long way from a valid pinchout region (i.e. where the
	// layer code will match the mode.layer_code), so we can spend ages here trying to get a valid 
	// candidate. This corresponds to....
	//_MS.mu = _MS.mode.quadratic_mean;
	
	// (2) simply use the point at the pinchout boundary as the mode mean. The shape will be totally
	// wrong, but will at least `cover' the valid pinchout region with modest probability. We then trust in the
	// MCMC Hastings terms to correct for the poor shape approximation. This corresponds to
	// _MS.mu = _MS.mode.QP_optimum;
	// Note that even these cases result in getting locked in some instances
	_MS.mu = _MS.mode.QP_optimum;

	Ltransp = _MS.mode.R;

	N = _MS.mu.size();         // dimension of full property suite
	//	n = _m.n_layers;     // dimension of times subvector
	n = _m.n_active_layers;     // dimension of times subvector
	_MS.property_sample = new double[N];
	double w_sum_sq;            // Will keep track of proposal log-likelihood
	double chisqDeviate=1.0;    // Will produce Gaussian samples by default

	// sample from doubly-updated Covariance (isopachs plus local times)
	w = new double[N];
	fail_count=0;
	do
	  {
		// make a sample
		// Variable 0 is t_top, variable n is t_base
		int matchLayerCodeAttemptCount=0;
		do
		  {
			if (useStudentOverDispersion)
			  {
				chisqDeviate =  chiSqareGenerator.nextDouble();
				// See Gelman "Bayesian Data Analysis" p 481 to get multivariate t samples;
				chisqDeviate = Math.sqrt(studentDOF/chisqDeviate);
			  }

			w_sum_sq=0.0;
			for (i=0;i<=n;i++) 
			  {
				w[i] = _MS.Ngen.nextDouble();
				w_sum_sq += w[i]*w[i];    // amounts to proposal negLogLik term
			  }
			for (i=0;i<=n;i++)
			  {
				s = 0;
				// From inverse of right cholesky decomposition,
				// Ltransp has true left cholesky factor, but transposed.
				for (j=0;j<=i;j++)
				  {
					// notice transpose!
					s += Ltransp.getQuick(j,i) *  w[j];
				  }
				// stretch to make into t-deviate
				s *= chisqDeviate;
				// add back the mean
				s += _MS.mu.getQuick(i);
				_MS.property_sample[i]=s;
				// s now has the random deviate.
				// set the property for this index
				_MS.setValuesFromIndex(i,s);
			  }
			// Ensure the times are ordered
			_MS.orderTimes();
			// if layer codes don't need to match we can bail out
			if (!insistOnMatchedLayerCode) break;

			if (_m._SI.verbosity>4)
			  {
				if (_MS.layer_code != _MS.mode.layer_code)
				  { 
					_m._SI.EH.println("Attempt "+jj+" : layer_code = "+_MS.layer_code + " != " +  _MS.mode.layer_code);
					jj++;
				  }
				else jj=0;
			  }
			// debugging stuff to see how stringent the layer-code matching is
			/*
			  if (_MS.layer_code != _MS.mode.layer_code) 
			  {
				_m._SI.EH.println("LC: miss "+_MS.mode.layer_code);
			  }
			else
			  {
				_m._SI.EH.println("LC: success "+_MS.mode.layer_code);
			  }
			  */
			matchLayerCodeAttemptCount++;

			// ** Confine to same layer configuration as mode is defined for ?? FIXME? 
			// But give up if this proves too hard (1000 attempts)
			}
		while (_MS.layer_code != _MS.mode.layer_code && 
			   matchLayerCodeAttemptCount < MAX_MATCH_LAYER_CODE_ATTEMPT_COUNT);

		// throw or return error if mode is faulty
		if (matchLayerCodeAttemptCount == MAX_MATCH_LAYER_CODE_ATTEMPT_COUNT)
		  {
			if (throwExceptionOnFaultyMode)
			  {
				throw new RuntimeException("Proving near-impossible to match layer_code="
										   +_MS.mode.layer_code+" with independence sampler."
										   +"Try using the Random walk sampler");
			  }
			else
			  {
				// flag this with magic answer
				return Double.MAX_VALUE;
			  }
		  }

		// if not faulty, now set remaining properties
		for (i=n+1;i<N;i++) 
		  {
			w[i] = _MS.Ngen.nextDouble();
			w_sum_sq += w[i]*w[i];    // amounts to proposal negLogLik term
		  }
		for (i=n+1;i<N;i++)
		  {
			s = 0;
			// From inverse of right cholesky decomposition,
			// Ltransp has true left cholesky factor, but transposed.
			for (j=0;j<=i;j++)
			  {
				// notice transpose!
				s += Ltransp.getQuick(j,i) *  w[j];
			  }
			// stretch to make into t-deviate
			s *= chisqDeviate;
			// add back the mean
			s += _MS.mu.getQuick(i);
			_MS.property_sample[i]=s;
			// s now has the random deviate.
			// set the property for this index
			_MS.setValuesFromIndex(i,s);
		  }
		// now update full layer properties
		sample_succeeded=_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
		// was it unphysical??
		if (!sample_succeeded) 
		  {
			fail_count++;
			_m._SI.EH.println("failure..."+fail_count);
		  }
		else
		  {
			fail_count=0; // All is well...
		  }

		// a failure shouldn't happen too many times or things are badly suspect
		if (fail_count>_MS.MAX_MC_FAIL_COUNT)
		  {
			throw new RuntimeException("Impossible fluid or rock properties in fluid state " +
									   (_MS.fs - 1) + ". Regression curves, isopachs, seismic, or "
									   + "fluid properties are probably suspect");
		  }
	  } while (!sample_succeeded);

	// approximate Prior + isopach (+Seismic) term
	double negLogProbProposed;
	if (!useStudentOverDispersion)
	  {
		negLogProbProposed = w_sum_sq/2.0;
	  }
	else
	  {
		// correct for t-distribution stretching....
		negLogProbProposed = 0.5*(studentDOF+N)*Math.log(1.0+w_sum_sq/studentDOF);
	  }

	// Must include determinant contribution?
	// negLogProbProposed	+= _MS.mode.logDetCtilde/2.0;


   	return negLogProbProposed;
  }

  private double computeProposalNegLogLikelihood()
  {
	double negLogLik;
	double negLogLik_iso,negLogLik_seis,negLogLik_depth;

	// (1a) SEISMIC mismatch term
	negLogLik_seis =  _MS.negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
	_MS.negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
	negLogLik = negLogLik_seis;

	// (1b) Add in logLikelihood contribution from ISOPACHS if they're required
	negLogLik_iso = _MS.computeIsopachNegLogLikelihood(false);
	negLogLik += negLogLik_iso;

    // (1c) Add in logLikelihood contribution from DEPTH self-consistency if required
	negLogLik_depth = _MS.computeConsistentDepthNegLogLikelihood(false);
	negLogLik += negLogLik_depth;

	// (1d) And from vertical stacking
	if (_m._SI.doStackedModel)
	  {
		double negLogLik_Vstack = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false);
		negLogLik += negLogLik_Vstack;
	  }

	return negLogLik;
  }
 
  private double computeProposalPriorNegLogProb()
  {
	return _MS.computePriorNegLogLikelihood();
  }
}








