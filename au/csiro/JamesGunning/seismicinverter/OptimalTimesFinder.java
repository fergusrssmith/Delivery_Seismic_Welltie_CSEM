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

import au.csiro.JamesGunning.lib.Minpack.optimization.*;
import au.csiro.JamesGunning.lib.util.MathUtil;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.Numerical.*;


import java.util.*;

import cern.colt.matrix.*;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.*;

/**
	 Class to find optimal (maximum-likelihood) layer times, with a
	 selection of other layer properties also possibly variable.  The
	 model top and base times are presumed fixed (no reflections, so no
	 contribution to posterior except through the prior)

	 The function to minimise is the sum of the negative log-likelihood error
	 terms from the prior, the seismic mismatch, and the isopach requirement.
  */	 
public final class OptimalTimesFinder extends Object implements Uncmin_methods 
{
  /**
	 optimising function penalty for Illegal rock parameters
  */
  public static final double GASSMAN_FAIL_PENALTY = 1e10;
  
  Model _m;
  int n;                        // n=dimensionality of times subvector of full model vector
  int n_sub;                    // dimension of current times sub-vector (=n if no-pinchout)
  int n_nontime;                // dimension of non-time subvector
  int n_tot;                    // sum of n_sub+n_nontime

  boolean reParameteriseToAdditiveTimes = false; // whether to impose t_i=sum_j^i tau_j type mapping
  SortableNumbers[] additiveMapsArray;           // mapping array to use in additive reparametriser


  double[] fullParameters;      // Contracting full suite of parameters, comprising union of times and non-time parameters
  double[] times;               // working vector for full suite of layer times
  double[] fullParameters_best; // Contracting vector: will contain the local optimal point at exit
  double[] grad_best;           // the local gradient at the local optimum at exit
  
  double t_scale;               // rough time scale for scaling

  private boolean neverApplyTimeOffset=false; // will force any offsetting to zero, no matter what the next option says
  private boolean applyTimeOffset = true;  // whether to apply time offset to get scaling better....
  private double time_offset=0.0;

  private int it_count=0;

  private TimePinchMap[] timeMapping;          // map of subset of times to actual times

  private int[] latch_count;
  boolean usePinchoutDetection = false;
  private byte[] layerIsPinchedOut;              // flag used to detect pinchouts on layer i=1...n
                                                  // 0=none, -1=full pinchout, +1= thinlayer
  private byte[]  latch_sense_last_evaluation;    // solver uses this to detect pinchout solutions. 

  // If usePinchoutDetection==true, any two times found to be 
  // closer together than ModelSampler.TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER on more than
  // LATCH_COUNT_THRESHOLD consecutive evaluations of the function f_to_minimise are taken to be
  // sufficient evidence that the solution is a pinchout, and the solution is aborted and re-solved for with
  // the pinchout enforced by reparametrising f_to_minimise(). A flakey trick, but works.
  // Good values for LATCH_COUNT_THRESHOLD seem to be about n.
  private int LATCH_COUNT_THRESHOLD;           

  protected boolean abortNow = false;

  protected double negLogIsopachLikelihood;    // terms available via getters() to get pieces of overall posterior
  protected double negLogDepthLikelihood;      // terms available via getters() to get pieces of overall posterior
  protected double negLogVStackLikelihood=0;   // these are set by each f_to_minimise call
  protected double negLogSeisLikelihood;       
  protected double negLogPriorLikelihood;     
  protected double negLogTruncationLikelihood; // term to force solutions away from problem time boundaries

  // Handles for non-Time properties if used
  protected ArrayList nonTimeProperties;

  // Control variables for solver. Needed at all levels of the pinchout recursion, so allocated once here
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



  /**
	 Constructor.
	 nonTimeProperties must be an ArrayList of PropertyMapper objects, pointing to
	 suitable layer properties to optimise over
  */

  public OptimalTimesFinder(Model m, ArrayList nonTimeProperties, boolean reParameteriseToAdditiveTimes)
  {
	this._m = m;
	this.reParameteriseToAdditiveTimes = reParameteriseToAdditiveTimes;

	// times = new double[_m.n_layers] is original definition
	// this.n is dimensionality of times vector
	this.n = _m.n_active_layers-1;

	times = new double[n+1];

	// Array of non-time properties
	this.nonTimeProperties = nonTimeProperties;

	grad_best = new double[n+1];

	usePinchoutDetection = _m._SI.usePinchoutDetection;
	latch_sense_last_evaluation = new byte[n+1];
	layerIsPinchedOut = new byte[n+1];
	latch_count = new int[n+1];

	// initialise sane time mapping
	formFullTimeMapping();

	// Reasonable guess for the pinchout detection repetition threshold number
	LATCH_COUNT_THRESHOLD = 2;
  }

  /**
	 Find optimal solution, perhaps starting at a new set of starting times, as given by
	 start_times
  */
  
  public double findOptimalSolution(double[] timesPlusOtherParameters)
  {
	int k;
	double[] sub_times;
	double f;

	setFlags(false);
	
	// size of time sub problem. On entry, this must be _m.n_layers-1
	n_sub = _m.n_active_layers-1;
	setProblemDimensions();

	// set up full parameter vector (incoming fullParameters must be of length _m.n_layers+n_nontime
	if (timesPlusOtherParameters.length!=n_tot+1) throw new RuntimeException("findOptimalSolution() called with wrong dimension");
	// elements [1...(_m.n_layers-1)] in fullParameters are the times.
	fullParameters = (double[])timesPlusOtherParameters.clone();
	fullParameters_best = (double[])timesPlusOtherParameters.clone();
	time_offset = timesPlusOtherParameters[1];

	// Copy time part of full parameters into times subvector
	System.arraycopy(timesPlusOtherParameters,1,times,1,n);	
	sub_times =  (double[])times.clone();
	f = findBestSubTimes(sub_times);
	
	// test if we aborted prematurely in order to find Pinchouts more efficiently
	while (wasAborted())
	  {
		// try again in new mode		
		sub_times = getSubTimesVectorForCurrentPinchState();
		if (_m._SI.verbosity>3)
		  {
			_m._SI.EH.print("               New starting subvector times ");
			for (k=1;k<sub_times.length;k++) _m._SI.EH.printf("%9.6g ",sub_times[k]);
			_m._SI.EH.print("%n");
		  }				
		setFlags(true); // resets abortNow flag
		f = findBestSubTimes(sub_times);  
		// At success, fullParameters has optimal solution, with no time-offsetting
		// and these have also been written back into times[]
	  }

	// Ensure final state is correct.
	// Write ***non-time*** properties back to calling vector.
	// Here fullParameters may be contracted if there are pinchouts, but write
	// the solution to the proper place [n+1.....n+n_nontime] of timesPlusOtherParameters
	System.arraycopy(fullParameters,n_sub+1,timesPlusOtherParameters,n+1,n_nontime);

	// reset flags so default (no-pinchout) mode is used by any future calls
	setFlags(false);
	// restore vanilla flavour time mapping
	formFullTimeMapping();

	// copy optimal ***times*** back to calling vector & current state
	System.arraycopy(times,1,timesPlusOtherParameters,1,n);	
	// and get respectable copy of entire uncontracted optimal parameters
	System.arraycopy(timesPlusOtherParameters,0,fullParameters_best,0,timesPlusOtherParameters.length);
	
	// Remap fullParameters to additive reparameterisation if required
	// from hereon, time-offsetting must be off, as we are outside the optimisation loop
	toggleTimeOffsetting(false);
	constructAdditiveTimesFromTrueTimes(timesPlusOtherParameters);

	// And finally, ensure these values are reflected in current state. Evaluate the function
	f_to_minimize(timesPlusOtherParameters);

	return f;
  }


  /*
	Discretised line derivative test for new starting points. Evaluates posterior density
	along straight line in parameter space (for all known modes with same fluid state) from
	in-situ starting point in _m.MS.property_sample to the mode MAP point, to check
	if such a path is purely "downhill". If so, we conjecture the start point in in the
	basin of attraction of the test mode, and the optimiser needn't be initiated
  */
  boolean testForPureDescentPath(double[] timesPlusOtherParameters,ArrayList modeIdentifiers)
  {
	int N_STEPS_IN_TEST_PATH = 30;
	double MINIMA_THRESHOLD = 0.2;
	int i,j,k,dim;
	double wt,fold=0;
	boolean pathDescends=false,good_state,performMinimaTest=false;
	double[] test_params, start_params=null;
	DoubleMatrix1D cmp_params;
	double deriv,val;
	ModelSampler.ModeIdentifier mi;
	double f=Double.MAX_VALUE,fmin=Double.MAX_VALUE;

	boolean debug = false;                      // Set false for production code

	double[] fvals=null,xvals=null;

	if (debug)
	  {
		fvals = new double[N_STEPS_IN_TEST_PATH]; // FIXME REMOVE
		xvals = new double[N_STEPS_IN_TEST_PATH]; // FIXME REMOVE
	  }

	if (modeIdentifiers.size()>0)
	  {
		dim = _m.MS.mu.size();

		// on entry, _m.MS.property_sample contains starting point for optimiser, correctly set up from
		// setOptimiserInitialValues().
		test_params = new double[dim];

		// loop over known modes with same fluid state
		// If there is a pure descent path to any of these, return true

		for (k=0;k<modeIdentifiers.size();k++)
		  {
			mi = (ModelSampler.ModeIdentifier)(modeIdentifiers.get(k));
			if ((mi.fs == _m.MS.fs) && (mi.mode != null)) 
			  {
				// compare only to other like fluid states
				cmp_params = mi.mode.QP_optimum;
				if (start_params==null) start_params = (double[])_m.MS.property_sample.clone();

				// We assume this is true until shown otherwise
				pathDescends = true;

				for (i=0;i<N_STEPS_IN_TEST_PATH;i++)
				  {
					wt = (double)i/(N_STEPS_IN_TEST_PATH-1.0);
					// set parameter vector along joining line from start configuration to QP optimum
					for (j=0; j<dim; j++)
					  {
						val = test_params[j] = wt * cmp_params.getQuick(j) + (1-wt)*start_params[j];
						_m.MS.setValuesFromIndex(j,val);
						_m.MS.property_sample[j]=val; // needed for prior
					  }
					_m.MS.orderTimes();
					good_state = _m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	
					if (good_state)
					  {
						f = _m.MS.debugComputeNegLogPosterior();
					  }
					else 
					  {
						// Catastrophic. Indicate ascent, and restore initial state, then return, as below
						pathDescends = false;
						break;
					  }
					if (debug)
					  {
						fvals[i] = f;
						xvals[i] = i;
					  }
					if (i>0)
					  {
						fmin = Math.min(fmin,f);
						deriv = (f-fold);
						if (deriv>MINIMA_THRESHOLD) 
						  {
							// Obvious uphill step
							pathDescends = false;
						    if (!debug) break; // needn't complete path if not in debugging mode
						  }
						else if (deriv>0)
						  {
							// A slight rise: we may be comparing to an optimum that isn't perfect...
							// Defer the check of this...
							performMinimaTest = true;
						  }
					  }
					fold=f;
				  }
				if (debug)
				  {
					String name = new String("DescentTestProfile.fs="+mi.fs+"."+k+"."+pathDescends+".txt");
					Util.printfileZO(name,xvals,fvals,fvals.length);
				  }

				// We require a pure descent path to any known mode. If demonstrated, we are done. No
				// need to examine other modes
				// Restore condition, then return true
				if (pathDescends) break;
			  }
		  }
		if (start_params==null)
		  {
			// no comparisons were made. Bail out
			return false;
		  }

		if (performMinimaTest)
		  {
			// value of f is at last point - supposedly the minimum
			if ((f-fmin) > MINIMA_THRESHOLD) pathDescends = false;
		  }

		// restore starting conditions
		for (j=0; j<dim; j++)
		  {
			_m.MS.setValuesFromIndex(j,start_params[j]);
			_m.MS.property_sample[j] = start_params[j]; // needed for prior
		  }
		_m.MS.orderTimes();
		good_state = _m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	


		// GC
		test_params = null;
		start_params=null;

		return pathDescends;
	  }
	else
	  {
		// no modes to compare to. Can't be a descent path
		return false;
	  }
  }


  // The actual minimiser call. times[] here is always the set of actual times, but possibly reduced
  // in dimensionality if pinchouts are imposed.
  private double findBestSubTimes(double sub_times[])
  {
	int i;
	
	setProblemDimensions();
	
	System.arraycopy(sub_times,1,fullParameters,1,n_sub);

	if (toggleTimeOffsetting(true))
	  {
		for (i=1;i<=n_sub;i++) fullParameters[i]-=time_offset;
	  }

	// construct fullParameters from additive reparameterisation if required
	reparametriseTimePropertiesUsingAdditiveScheme(fullParameters);

	// add default (in-situ) non-time properties to parameter vector
	fillParameterVectorWithNonTimeProperties(fullParameters);

	// recompute effective properties for all layers
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false);

	// needed workspace
	double[][] Hessian_wk = new double[n_tot+1][n_tot+1];
	double[] Hessian_diag_wk = new double[n_tot+1];
	// fullParameters_best = new double[n_tot+1];
	double[] grad_sub_best = new double[n_tot+1];
	double[] typical_values = new double[n_tot+1];

	// Typical scaling. take care here as it affects convergence
	for (i=1;i<=n_sub;i++)
	  {	
		typical_values[i] = getCharacteristicLayerTimeForThisIndex(i);
		// For some reason, Gauss-Newton updates don't like too-small time-scales, which
		// may be imposed by requiring tight isopach constraint: I don't understand this yet :-d
		// So if using Gauss-Newton, don't make typical scales any smaller than _m.SM.delta_t_indep_samples
		// if (!_m._SI.useSecantOptimisation) typical_values[i] = _m.SM.delta_t_indep_samples;
	  }
	if (reParameteriseToAdditiveTimes) 
	  {
		// parameters 2...n_sub are additive times parametrised as squares. Scale suitably
		for (i=2;i<=n_sub;i++)
		  {	
			typical_values[i] = Math.sqrt(typical_values[i]);
		  }
	  }
  
	// use default values also for scaling here....
	fillParameterVectorWithNonTimeProperties(typical_values);
	// under certain transformations, fullParameters[i] has default value 0.0, and scale 1,  so...
	for (i=n_sub+1;i<=n_tot;i++) if (typical_values[i]==0.0) typical_values[i]=1.0;

	fscale[1]=(double)n_tot;  // assumes chisq \sim n_tot if properly scaled

	// Tolerances...
	double epsilon  = 1.12e-16; // 1.12e-16;
	// Notes: gradient tolerance weakening does not seem to help much
	gradient_tol[1] = Math.pow(epsilon,1.0/3.0);
	gradient_tol[1] = 0.01;
	max_step_size[1] = 2.0*Math.sqrt(n_tot); // Default is 0.0, this nonzero guess is a tight "factor of 2 per dimension" ansatz
	step_tol[1]=Math.sqrt(epsilon);
	// step_tol[1]=1e-6;
	trust_radius[1]=-1;  // -1 is default

	// Flags for methods/outputs...
	method[1]=1;          // 1 = line-search
	// NB. UNCMIN can be a bit wierd under BFGS if scaling is dodgy. I believe scaling issues are
	// fixed so BFGS (Secant) is safe
	is_exp[1] = 1;        // _m._SI.useSecantOptimisation ? 1 : 0;
	msg[1]=8;             // msg=8 seems to suppress most of the noise... :-). msg=0,1 is very noisy
	n_good_digits[1]=-1;  // -1 gives defaults
	max_its[1]=200;
	is_gradient_supplied[1]=0;
	is_hessian_supplied[1]=0;

	

	// Solver run-time state variables
	// (1) iteration counts. Don't reset if this is a second pass with pinchouts enforced
	if (!areThereAnyPinchouts()) it_count=0;
	
	// (2) are non-time variables at their means? Check this and set flags. This affects efficiency

	// Now do it !
	Uncmin_f77.optif9_f77(n_tot, fullParameters, this,
						  typical_values,fscale,method,is_exp,msg,n_good_digits,
						  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
						  gradient_tol,max_step_size,step_tol,
						  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
	// if a premature abort, get out quickly
	if (abortNow) return 0;

	if (info[1]>=5)
	  {
		_m._SI.EH.println("UNCMIN: Solution did not converge!!!");
	  }
	if (_m._SI.verbosity>3) 
	  {
		_m._SI.EH.print("Minimiser its = "+it_count+" subvector solution-times ");
		for (i=1;i<=n;i++) _m._SI.EH.printf("%9.6g ",times[i]);
		_m._SI.EH.print("%n");
	  }

	// Copy final solution to fullParameters vector
	System.arraycopy(fullParameters_best,1,fullParameters,1,n_tot);
	// and map these back to the original times vector
	mapSubTimesVectorToCompleteTimes(fullParameters,false);

	Hessian_wk = null;
	Hessian_diag_wk = null;
	grad_sub_best = null;
	typical_values =null;

	return (double) f_best[1];
  }

  protected boolean toggleTimeOffsetting(boolean on)
  {
	if (on) 
	  {
		applyTimeOffset=true;
	  }
	else
	  {
		applyTimeOffset=false;
	  }
	if (neverApplyTimeOffset) applyTimeOffset=false;
	return applyTimeOffset;
  }

  private void setProblemDimensions()
  {
	// size of non-time variables vector
	n_nontime = (nonTimeProperties!=null)? nonTimeProperties.size(): 0;
	// full problem size
	n_tot = n_sub + n_nontime;
  }
 
  // Reset n_nontime, n_sub & n_tot fields
  void setStandardProblemDimensions()
  {
	n_nontime = (nonTimeProperties!=null)? nonTimeProperties.size(): 0;
	// update full problem size to keep things coherent
	n_sub = n;
	n_tot = n_sub + n_nontime;
  }

  /*
	Remap times to additive scheme. On input fullParameters[1...n]
	contain the layer times in ordinary ordering. After this call completes
	fullParameters[1] is unchanged,
	additiveMapsArray[0...nsub-1] is an array with..
	additiveMapsArray[i].index   giving the time index of the next layer time
	fullParameters[i]=T_i, with property <br>
	t_{index[i-1]}  = t_{index[0]} + sum_{j=2}^{i} fullParameters[j]^2 <br>
	with also<br>
	 t_{index[0]} = fullParameters[1]  (+ time_offset if required )<br>
  */
  private void reparametriseTimePropertiesUsingAdditiveScheme(double[] fullParameters)
  {
	int i;
	if (reParameteriseToAdditiveTimes)
	  {
		// times are in elements fullParameters[1...n_sub], and there are no pinchouts
		additiveMapsArray = new SortableNumbers[n_sub];
		double[] times_copy = new double[1+n_sub];
		System.arraycopy(fullParameters,1,times_copy,1,n_sub);
		for (i=1;i<=n_sub; i++)
		  {
			additiveMapsArray[i-1] = new SortableNumbers(fullParameters[i],i);
		  }
		// sort into desired order
		Arrays.sort(additiveMapsArray);
		for (i=1;i<=n_sub; i++)
		  {
			if (i==1) 
			  {
				fullParameters[1]=additiveMapsArray[0].value;
				//if (applyTimeOffset) fullParameters[1] -= time_offset;
			  }	
			else
			  {
				fullParameters[i]=Math.sqrt(times_copy[additiveMapsArray[i-1].index] - 
											times_copy[additiveMapsArray[i-2].index]);
			  }
		  }
		times_copy = null;
	  }
  }

  /*
	Remap addtive times in fullParameters[1..n] to true times, in place,  using additive scheme 
  */
  private void constructAdditiveTimesFromTrueTimes(double[] fullParameters)
  {
	if (reParameteriseToAdditiveTimes)
	  {
		// times are in elements fullParameters[1...n_sub], and there are no pinchout remappings
		double[] times_copy = new double[1+n_sub];
		System.arraycopy(fullParameters,1,times_copy,1,n_sub);
		for (int i=1;i<=n_sub; i++)
		  {
			if (i==1) 
			  {
				//fullParameters[i]=additiveMapsArray[0].value;
				fullParameters[i]=times_copy[additiveMapsArray[0].index];
				if (applyTimeOffset) fullParameters[i]-=time_offset;
			  }
			else
			  {
				fullParameters[i]=Math.sqrt(times_copy[additiveMapsArray[i-1].index] - 
											times_copy[additiveMapsArray[i-2].index]);
			  }
		  }
		times_copy = null;
	  }
  }


  /*
	fill fullParameters[n_sub+1,...] with non-time properties
  */
  void fillParameterVectorWithNonTimeProperties(double[] fullParameters)
  {
	if (nonTimeProperties!=null)
	  {  
		int i,j,k;
		PropertyMapper PM;
		for (i=n_sub+1, j=0; j<n_nontime; i++,j++)
		  {
			// get the Mapper for this property in the parameters
			PM = (PropertyMapper) nonTimeProperties.get(j);
			// in the linear algebra setup of modelSampler, this property has index PM.propertyMasterIndx
			// in the HashMaps map_to, map_from etc.
			// The propertyMasterIndx will point to an Array (usually length 1) of
			// properties synchronised if necessary. The first in the list is sufficient to give
			// us the current value
			ArrayList list = (ArrayList) _m.MS.map_from.get(PM.propertyMasterIndx);

			LayerProperties.Property prop = (LayerProperties.Property)list.get(0);
			// The real value is prop.value, but we may want to transform it before stuffing
			// into the fullParameters list...
			fullParameters[i] = PM.getTransformedValue(prop.value);
			// and keep current state snapshot up-to-date
			_m.MS.property_sample[PM.propertyMasterIndx] = prop.value;

			// repeat above on any coupled properties
			if (PM.isACoupledMapper)
			  {
				// compute coupled value in "transformed" space
				double base_value = fullParameters[i];
				for (k=0; k<PM.coupledProperties.length; k++)
				  {
					// The real value is prop.value, but we may want to transform it before stuffing
					// into the fullParameters list...
					PropertyMapper PM_coup = PM.coupledProperties[k];
					ArrayList list_coup = (ArrayList) _m.MS.map_from.get(PM_coup.propertyMasterIndx);
					prop = (LayerProperties.Property)list_coup.get(0);
					// compute coupled value in "transformed" space
					double coupled_value = PM.coupledCoeffs[k] * base_value + PM.coupledConst[k];
					// compute true value
					double true_value = PM_coup.getTrueValue(coupled_value);
					// and set the coupled Layer properties suitably
					_m.MS.setValuesFromIndex(PM_coup.propertyMasterIndx,true_value);
					// and keep current state snapshot up-to-date
					_m.MS.property_sample[PM_coup.propertyMasterIndx] = prop.value;
				  }
			  }
		  }
	  }
  }


  /*
	fill fullParameters[n_sub+1,...] with non-time properties
  */
  private void fillNonTimePropertiesFromParameterVector(double[] fullParameters)
  {
	if (nonTimeProperties!=null)
	  {  
		int i,j,k;
		double true_value;
		PropertyMapper PM;
		for (i=n_sub+1, j=0; j<n_nontime; i++,j++)
		  {
			// get the Mapper for this property in the parameters
			PM = (PropertyMapper) nonTimeProperties.get(j);
			// in the linear algebra setup of modelSampler, this property has index 
			// PM.propertyMasterIndx in the HashMaps map_to, map_from etc.
			true_value = PM.getTrueValue(fullParameters[i]);
			// and set the Layer properties suitably
			_m.MS.setValuesFromIndex(PM.propertyMasterIndx,true_value);
			// and update the current 1D vector snapshot
			_m.MS.property_sample[PM.propertyMasterIndx] = true_value;

			if (PM.isACoupledMapper)
			  {
				// compute coupled value in "transformed" space
				double base_value = fullParameters[i];
				for (k=0; k<PM.coupledProperties.length; k++)
				  {
					// The real value is prop.value, but we may want to transform it before stuffing
					// into the fullParameters list...
					PropertyMapper PM_coup = PM.coupledProperties[k];
					ArrayList list_coup = (ArrayList) _m.MS.map_from.get(PM_coup.propertyMasterIndx);
					LayerProperties.Property prop = (LayerProperties.Property)list_coup.get(0);
					// compute coupled value in "transformed" space
					double coupled_value = PM.coupledCoeffs[k] * base_value + PM.coupledConst[k];
					// compute true value
					true_value = PM_coup.getTrueValue(coupled_value);
					// and set the coupled Layer properties suitably
					_m.MS.setValuesFromIndex(PM_coup.propertyMasterIndx,true_value);
					// and keep current state snapshot up-to-date
					_m.MS.property_sample[PM_coup.propertyMasterIndx] = prop.value;
				  }
			  }
		  }
	  }
  }

  // Return if there are any pinchouts
  private boolean areThereAnyPinchouts()
  {
	int i;
	for (i=1;i<=n;i++)
	  {
		if (layerIsPinchedOut[i]!=0) return true;
	  }
	return false;
  }
	
  /** 
	  A debugging method to examine the posterior surface near the ML
	  points. The surface generally looks smooth and quadratic, but definite steps are
	  evident where pinchouts come into play.
  */
  private void debugPrintXsections(String prefix)
  {
	int i,j,N=100;
	int NSTEPS=50;
	double DT_STEP = 0.002;
	double NONTIME_FRACTIONAL_RANGE = 0.2;
	double[] fvals=new double[N+1];
	double[] d1fvals=new double[N];
	double[] d2fvals=new double[N-1];


	_m._SI.EH.print("Solution: ");
	for (i=1;i<=fullParameters_best.length-1;i++) _m._SI.EH.printf("%11.8g ",fullParameters_best[i]);
	_m._SI.EH.print("%n");

	for (int k=1;k<=n_tot;k++)
	  {
		// restore values
		if (k>n) _m._SI.EH.println("Non-time-property "+(k-n));
		System.arraycopy(fullParameters_best,1,fullParameters,1,fullParameters_best.length-1);
		boolean zero_mean = (fullParameters_best[k]==0.0);
		for (j=0,i=-NSTEPS;i<=NSTEPS;i++,j++)
		  {
			if (k<=n)
			  { // times
				fullParameters[k] = fullParameters_best[k]*(1+DT_STEP*((double)i/(double)NSTEPS));
			  }
			else
			  {
				// non-time properties
				if (!zero_mean)
				  {
					fullParameters[k] = fullParameters_best[k]*(1+NONTIME_FRACTIONAL_RANGE*((double)i/(double)NSTEPS));
				  }
				else
				  {
					fullParameters[k] = NONTIME_FRACTIONAL_RANGE*((double)i/(double)NSTEPS);
				  }
			  }

			fvals[j]=f_to_minimize(fullParameters);

			if (j>0) d1fvals[j-1]=fvals[j]-fvals[j-1];
			if (j>1) d2fvals[j-2]=fvals[j]+fvals[j-2]-2*fvals[j-1];
		  }
		String key = new String("."+_m.MS.fs+"."+_m.MS.localModes.size()+"."+k+".txt");
		Util.printfileZO(prefix+key,fvals,fvals.length);
		_m._SI.EH.println("Wrote "+prefix+key);
		/*
		  Util.printfileZO("d1_"+prefix+key, d1fvals, d1fvals.length);
		  _m._SI.EH.println("Wrote d1_"+prefix+k ey);
		  Util.printfileZO("d2_"+prefix+key, d2fvals, d2fvals.length);
		  _m._SI.EH.println("Wrote d2_"+prefix+ke y);
		*/
		// restore fullParameters....
		System.arraycopy(fullParameters_best,1,fullParameters,1,fullParameters_best.length-1);
		// restore full system configuration
		f_to_minimize(fullParameters);
	  }
  }

  

  /**
	 The solver throws execution to here everytime it gets a major new
	 point. This is where we can monitor  things, like, e.g. convergence towards
	 a pinchout-like solution (a crease in the optimising surface)
  */
  public void updateNewPoint(double fnew, double[] xnew)
  {
	int l,j,k;
	// monitor pinchouts
	if (usePinchoutDetection && !reParameteriseToAdditiveTimes) 
	  {
		countTimeLatchingOccurrences();
	  }

	// Feedback a snapshot of new point if desired
	if (_m._SI.verbosity>3)
	  {
		_m._SI.EH.printf("It: %4d f=%11.8g layer_code %d parameters: ",it_count, fnew, _m.MS.layer_code);
		// time parameters
		for (l=1;l<=n_sub;l++) _m._SI.EH.printf("%11.8g ",xnew[l]);
		// non-time parameters
		if (n_nontime>0)
		  {
			_m._SI.EH.printf(" | ",xnew[l]);
			for (j=0,l=n_sub+1;l<=n_tot;j++,l++) 
			  {
				PropertyMapper PM;
				// get the Mapper for this property in the parameters
				PM = (PropertyMapper) nonTimeProperties.get(j);
				if (!PM.isACoupledMapper)
				  {
					// in the linear algebra setup of modelSampler, this property has index 
					// PM.propertyMasterIndx in the HashMaps map_to, map_from etc.
					_m._SI.EH.printf("%11.8g ",PM.getTrueValue(xnew[l]));
				  }
				else
				  {
					// Just print the master parameter for coupled Mappers, with _ suffix
					// in the linear algebra setup of modelSampler, this property has index 
					// PM.propertyMasterIndx in the HashMaps map_to, map_from etc.
					_m._SI.EH.printf("%11.8g_ ",PM.getTrueValue(xnew[l]));	
				  }
			  }
		  }
		_m._SI.EH.println("");
		
		_m._SI.EH.print("          complete time parameters: ");
		for (l=1;l<=n;l++) _m._SI.EH.printf("%11.8g ",times[l]);
		_m._SI.EH.println("");
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
	Arrays.fill(latch_count,0);
	if (setOrReset)
	  {
		abortNow=false;
	  }
	else
	  {
		it_count=0;
		Arrays.fill(layerIsPinchedOut,(byte)0);
		Arrays.fill(latch_sense_last_evaluation,(byte)0);
	  }
  }


  /**
	 For a given state of layer-pinchouts, the function can be parametrised in fewer
	 dimensions. E.g. if we know min f(x,y) falls on x=y, we need minimise only f(x,x).
	 Under pinchouts, 'clusters' of adjacent time variables get locked together, so we
	 need to produce a new state vector which points only to one variable in each cluster, for
	 reduced dimensionality in the minimisation routine.

	 t_1.... t_n is the full suite of times that will be broken up into clusters.

	 Once we have l = 1...n_sub clusters, the times are set by setting

	 times[TimePinchMap[l].timeIndx[kk]] = t_sub[l] + TimePinchMap[l].timeOffset[kk]
	 
	 for each kk=0..(TimePinchMap[l].length-1)
  */
  private double[] getSubTimesVectorForCurrentPinchState()
  {
	int nFreeInCluster;
	int nTimesAccountedFor=0;
	int i,j,k,kk,l;
	double tmin,tmax, toffset;
	double[] t_sub;

	n_sub=0;

	if (!usePinchoutDetection)
	  {
		return times;
		// FIXME to be clean should make sure timeMapping is coherent with this 
	  }
	else
	  {
		// Pass 1. Get the size of the subvector required by cluster counting
		for (i=1;i<=n;i++) 
		  {
			 n_sub++;
			 if (layerIsPinchedOut[i]!=0)
			   {
				 // Beginning of a isolated sequence of linked layers. Set the time
				 // Search ahead until the sequence of layers ends
				 for (j=i+1;j<n;j++) if (layerIsPinchedOut[j]==0) break;
				 // Ok now the end of the sequence of linked layer is indexed j-1
				 // This means the times t_i, ... t_j are in a "cluster" with j-i+1 members
				 // Kick encompassing loop to start at j+1...
				 i=j;   // for(;;i++) will kick this up to j+1
			   }
		  }
		// create the new timeMapping array associated with the subvector t_sub
		timeMapping = new TimePinchMap[n_sub+1];
		// create the subvector (unit offset, so of dimension n_sub+1)
		t_sub = new double[n_sub+1];
		l=1;  // index to run over t_sub entries as we find them
		// populate it with suitable values from current state
		tmin = tmax = times[1];
		// Now break the time sequence into clusters. layers are Free (F), or Pinched-out (P)
		// layerIsPinchedOut[j]==0 means F,   layerIsPinchedOut[j]==1 means P
		// Clusters have to be formed as per this kind of pattern: [F][F][F][PPPPF][PF][F][PF]
		// A cluster is defined as the largest group we can bracket with only one "free" layer
		// where the right-most member in the cluster is "free"
		for (i=1,nTimesAccountedFor=0; i<=n && (l<=n_sub); i++,l++) 
		  {
			for (nFreeInCluster=0,j=i;j<=n;j++) 
			  {
				if (layerIsPinchedOut[j]==0) nFreeInCluster++;
				if (nFreeInCluster==2) 
				  {
					// if we span a pinchout layer or group, the clustering should group like this
					// F=Free, P=Pinched:    [F][PF], or [F][PPPPF]
					if (j-i>1 && layerIsPinchedOut[i]==0)
					  {
						// This deals with the case FP...F, from which we want to break out the cluster [F]P...F
						j=i+1;  // so this cluster is [F], the next cluster will be the big one [P...F]
					  }
					nFreeInCluster=1;
					break; // end-of-cluster index = j-1
				  }
			  }
			// Ok here our cluster spans the indices i...j-1
			nTimesAccountedFor += (j-i);
			if (nFreeInCluster==1)
			  {
				// Beginning of a isolated sequence of linked layers, indexed i....j-1
				// Latch times sensibly over this subgroup of linked layers i...j-1
				// This means the times t_i, ... t_{j-1} are in a "cluster", with j-i members
				// The most important thing is to keep everything in the cluster 
				// behind the most exposed layer 
				// which has tmax = max{t_i, ... t_{j-1}}
				tmin = tmax = times[i];
				for (k=i+1; k<=j-1; k++)
				  {
					if (times[k] > tmax) tmax = times[k];
					if (times[k] < tmin) tmin = times[k];					
				  }
				// now latch all times wrt tmax. Let the time that has tmax be the new sub parameter
				t_sub[l]=tmax;
				timeMapping[l] = new TimePinchMap(j-i);
				// _m._SI.EH.println("TIME MAPPING("+l+") size = "+(j-i)); // FIXME DELETE
				// We want to preserve the mapping t_sub[l] <-> tmax, so hang everything in
				// this cluster off t_sub[l].
				// Latch all times wrt tmax, preserving the ordering of times[i...j] in the
				// interval (tmax,tmin). The interval (tmax,tmin) is mapped to 
				// (t_sub[l] , t_sub[l] - _m.MS.TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER)
				for (k=i,kk=0; k<=j-1; k++,kk++)
				  {
					timeMapping[l].timeIndx[kk]=k;
					if (tmax!=tmin) 
					  {						
						toffset = -(_m.MS.TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER
						          * (tmax - times[k])/(tmax-tmin));
					  }
					else
					  {
						toffset=0;
					  }
					timeMapping[l].timeOffset[kk] = toffset;
				  }
				// and kick encompassing loop to start at j...
				i=j-1;   // for(;;i++) will kick this up to j
			  }
			else
			  {
				// something's wrong
				_m._SI.EH.println("%nWARNING!!!!!: Clustering algorithm buggy FIXME");
			  }	
		  }
		if (nTimesAccountedFor != n)
		  {
			// something's wrong
			_m._SI.EH.println("%nWARNING!!!!!: Clustering algorithm buggy. Clusters incomplete. FIXME."
							  +" n-times in clusters = "+nTimesAccountedFor+" n-times="+n);
		  }	
	  }
	return t_sub;
  }

  // Make vanilla flavour time Mapping for no-pinchout state
  private void formFullTimeMapping()
  {
	int i;

	timeMapping = new TimePinchMap[n+1];
	for (i=1;i<=n;i++) 
	  {
		timeMapping[i] = new TimePinchMap(1);
		timeMapping[i].timeIndx[0]=i;
		timeMapping[i].timeOffset[0]=0.0;
	  }
	n_nontime = (nonTimeProperties!=null)? nonTimeProperties.size(): 0;
	// update full problem size to keep things coherent
	n_sub = n;
	n_tot = n_sub + n_nontime;
  }

  /* 
	 map the Additive times back to true ones using the current additiveMapsArray
   */

  protected void mapSubTimesVectorToCompleteTimes(double[] parameters, boolean mapParametersInPlaceSansOffsets)
  {
	int i,j;
	double t_i;

	for (i=1; i<=n_sub; i++)
	  {
		if (reParameteriseToAdditiveTimes)
		  {
			// use additive times reparameterisation
			if (i==1)
			  {
				// Fix the smallest time....
				times[additiveMapsArray[0].index] = parameters[1];
			  }
			else
			  {
				// compute all subsequent (later) times.
				// Make time overlaps strictly impossible: This mapping will then hold....
				// t_i = t_{i-1}+parameters[i]*parameters[i] + TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER;
				// so t_i > t_{i-1} for all parameters[i]
				times[additiveMapsArray[i-1].index] = times[additiveMapsArray[i-2].index]
				                                     + parameters[i]*parameters[i] 
				                                     + _m.MS.TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER;
			  }
		  }
		else
		  {
			// direct times are parameters & pinchouts allowed, so have to use
			// timeMapping
			t_i = parameters[i];
			for (j=0; j<timeMapping[i].size(); j++)
			  { 
				times[timeMapping[i].timeIndx[j]] = t_i + timeMapping[i].timeOffset[j];
			  }
		  }
	  }
	// restore times to true values if offset
	if (applyTimeOffset && !mapParametersInPlaceSansOffsets)
	  {
		for (i=1; i<=n; i++) times[i] += time_offset;
	  }

	// and write these values into the in-situ fields
	for (i=1; i<=n; i++)
	  {
		_m.MS.setValuesFromIndex(i,times[i]);
		// update current snapshot too
		_m.MS.property_sample[i] = times[i];
	  }
	// overwrite times element if required also. This won't be necessary if there
	// is no reparameterisation (reParameteriseToAdditiveTimes==false)
	if (mapParametersInPlaceSansOffsets && reParameteriseToAdditiveTimes)
	  {
		for (i=1; i<=n; i++)
		  {
			parameters[i] = times[i];
		  }
	  }
  }

  /* 
	 Get the characteristic time-scale associated with a particular subtimes parameter
	 
   */

  private double getCharacteristicLayerTimeForThisIndex(int subtimes_indx)
  {
	int i,layer_no,i_layer;
	double dt, dt_iso_scale;

	if (reParameteriseToAdditiveTimes)
	  {
		// use additive times reparameterisation
		layer_no = additiveMapsArray[subtimes_indx-1].index;
	  }
	else
	  {
		// direct times are parameters & pinchouts allowed, so have to use timeMapping
		layer_no = timeMapping[subtimes_indx].timeIndx[0];
	  }

	// Now, at this point, layer_no is the true (master) layer index of the time
	// which is mapped by subtimes_indx in the times subvector.
	// If this layer or the one above has an isopach, we should get the
	// characteristic times associated with the isopach constraints for these layers
	// If there is no isopach, or the constraint is very weak, the 
	// characteristic time scale should be the delta_t_indep_samples scale associated
	// with the seismic model (usually about 10 ms or so)

	// Get isopach time scales.
	dt=Double.MAX_VALUE;
	if (!(_m.MPD.n_isopachs==0 || _m._SI.ignoreIsopachs))
	  {
		dt=Double.MAX_VALUE;
		// run over this layer and the one above
		for (i=Math.max(0,layer_no-1);  i<=Math.min(_m.n_active_layers-1,layer_no); i++)
		  {
			i_layer = _m.MS.master_layer_number_sequence[i];
			ModelPriorDescriptions.MetaInfo MI = _m.MPD.MetaInfo[i_layer];
			LayerProperties LPI = _m.LP[i_layer];
			if ((MI.dz!=null) && !(MI.dz.value==0 && MI.sigma_dz.value==0)) // layer has an real isopach
			  {
				// approximate time-thickness scale of isopach (from severity of isopach constraint)
				// factor of 2 from two-way times
				dt_iso_scale =  2.0*MI.sigma_dz.value/LPI.vp_eff;
				dt = Math.min(dt,dt_iso_scale);
			  }
		  }
	  }
	// return a sensible time scale
	if (dt==0.0 || dt==Double.MAX_VALUE)
	  {
		dt = _m.SM.delta_t_indep_samples;
	  }
	
	// compare to natural seismic time scale. Choose the least (fastest scale-of-variation of the objective)
	dt = Math.min(dt,_m.SM.delta_t_indep_samples);

	return dt;
  }

	

  /**
	 Here parameters[] is the union of the times and non-times properties we wish to
	 optimise over. The times are the components parameters[1....n_sub], and the
	 non-time variables are parameters[n_sub+1...n_tot]. The time component is the 
	 set of independent `cluster times' used in smart pinchout detection.
  */

  public double f_to_minimize(double parameters[])
  {
	int i,N;
	double negLogLik;

	// testForNaN(_m.LP[1].t.value);
	  
	// (1) Constraint on times to fit into total reservoir interval
	testAndImposeTimeTruncation(parameters,true);

	// map the sub-vector parameters[1...n_sub] onto times[]
	mapSubTimesVectorToCompleteTimes(parameters,false);

	
	  
	// Copy times into LayerProperties data structure...
	// the first indices in the reduced model maps seen
	// by setValuesFromIndex() are times
	// (Remember Uncmin works with unit-offsets....)
	// These times then have to be ordered to account for pinchouts etc by
	// orderTimes()

	if (_m._SI.verbosity>5)
	  { 
		_m._SI.EH.print("     Layer times: ");
		for (i=1;i<=n;i++) _m._SI.EH.printf("%11.8g ",times[i]); 
		_m._SI.EH.print("%n");
	  }
	// order the times
	_m.MS.orderTimes();
	
	if (_m._SI.verbosity>5)
	  { 
		_m._SI.EH.print("Post-ordering  t: ");
		for (i=1;i<=n;i++) _m._SI.EH.printf("%11.8g ",_m.LP[i].t.value);
		_m._SI.EH.print("%n");		
	  }

	// update all non-time properties
	fillNonTimePropertiesFromParameterVector(parameters);

	if (nonTimeProperties!=null && nonTimeProperties.size()>0)
	  {
		// recomputing effective properties for all layers is necessary
		if (true)
		  {			
			// Gassman -handling case
			if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false))
			  {
				_m._SI.EH.println("Optimiser hit illegal Gassman state...");
				return GASSMAN_FAIL_PENALTY;
			  }
		  }
		else
		  {
			// Gassman catastrophic case
			_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false);
		  }
	  }
	

	// (1) Acquire logLikelihood contribution from seismic mismatch
	negLogLik =  _m.SM.computeSyntheticAndSeismicError();
	if (_m._SI.verbosity>5)
	  {
		_m._SI.EH.println("Likelihood: Seismic mismatch contribution: %-11.8g ",negLogLik);
	  }
	negLogSeisLikelihood = negLogLik;

	// (2) Acquire logLikelihood contribution from prior and isopachs

	// USE EXACT ISOPACH LIKELIHOOD WITH PRIOR
	// Acquire logLikelihood contribution from isopachs if we're not ignoring them
	// and there are some isopachs to attend to
	negLogIsopachLikelihood = _m.MS.computeIsopachNegLogLikelihood(false);
	if (_m._SI.verbosity>5)
	  {
		_m._SI.EH.println("Likelihood: Isopach contribution: %-11.8g ",negLogIsopachLikelihood);
	  }
	negLogLik += negLogIsopachLikelihood;

	// (3) Contribution from depth self-consistency
	negLogDepthLikelihood = _m.MS.computeConsistentDepthNegLogLikelihood(false); 
	if (_m._SI.verbosity>5)
	  {
		_m._SI.EH.println("Likelihood: Depth self-consistency contribution: %-11.8g ",negLogIsopachLikelihood);
	  }
	negLogLik += negLogDepthLikelihood;

	// (3b) Contribution from vertical stacking self-consistency
	if (_m._SI.doStackedModel)
	  {
		negLogVStackLikelihood = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
		if (_m._SI.verbosity>5)
		  {
			_m._SI.EH.println("Likelihood: Vertical stacking contribution: %-11.8g ",negLogVStackLikelihood);
		  }
		negLogLik += negLogVStackLikelihood;
	  }

	//  (4) Also the contribution from the raw prior 
	// (i.e. before isopachs were used in the linearised update)
	// Acquire logLikelihood contribution from prior.
	// The times are the variables indexed 0...(_m.n_active_layers-1)
	// top and base are fixed at their modes. 	

	double sum=0.0;
	DoubleMatrix1D mu_prior_only =  _m.MS.mbar_prior_only[_m.MS.fs];
	DoubleMatrix2D Cinv_prior_only = _m.MS.Cinv_prior_only[_m.MS.fs];
	N = mu_prior_only.size();
	
	// Compute prior contribution. There may be efficiencies to be gained here (only a 
	// subvector varies), but for now, be slow but safe.
	if (true)
	  {
		negLogPriorLikelihood = _m.MS.computePriorNegLogLikelihood();
	  }
	else
	  {
		// old........

		// Have to recompute everything. Form current state vector
		// The leading part will have changes (times), but also a few values
		// lower down in the vector will change 

		DoubleMatrix1D tmp = new DenseDoubleMatrix1D(N);
		DoubleMatrix1D ans = tmp.like();
		for (i=0;i<mu_prior_only.size();i++) 
		  {
			tmp.setQuick(i,_m.MS.property_sample[i] - mu_prior_only.getQuick(i));
		  }
		// form C^1.tmp -> ends up in-place in "ans"
		Cinv_prior_only.zMult(tmp,ans,1.0,0.0,false);
		// then tmp^T.C^-1.tmp/2.0
		sum = tmp.zDotProduct(ans);
		tmp = null;
		ans = null; 

		// Add to this the  determinant term...
		// (note prior fluid combo term cancels in MCMC method so omitted here)
		// FIXME: I think I can omit the determinant here. Yes????
		// negLogPriorLikelihood = sum/2.0 + _m.MS.logDetC_prior_only[_m.MS.fs]/2.0;
		negLogPriorLikelihood = sum/2.0;
	  }



	negLogLik += negLogPriorLikelihood;
	
	if (_m._SI.verbosity>5)
	  {
		_m._SI.EH.println("Likelihood: Prior contribution: %11.8g ",
						  negLogPriorLikelihood);
	  }


	if (_m._SI.verbosity>5)  _m._SI.EH.println("f="+negLogLik);

	it_count++;

	// Add truncation error too
	negLogLik+= negLogTruncationLikelihood;

	//if (_m._SI.verbosity>2)  _m._SI.EH.println("IT_count= "+it_count+" "+negLogLik);

	return negLogLik;
  }


  /** 
	  Tests whether optimiser has times at unacceptable values.
	  All times must lie in the interval <br>
	  ( _m.MPD.MetaInfo[0].t.value < t <  _m.MPD.t_base.value ) <br>
	  I.e. prior mean value of model-region top, and base
  */
  protected boolean testAndImposeTimeTruncation(double[] parameters, boolean offsettingPossible)
  {
	int i;
	double t_i=0.0,t_iminus1=0.0;
	negLogTruncationLikelihood=0.0;

	boolean truncated = false;
	for (i=1;i<=n_sub;i++) 
	  {
		if (reParameteriseToAdditiveTimes)
		  {
			if (i==1)
			  {
				t_i = t_iminus1 = parameters[i];
			  }
			else
			  {
				// additive non-negative mapping. Pinchouts never allowed so
				// t_{i} = times[i] always.
				t_i = t_iminus1+ parameters[i]*parameters[i];
				t_iminus1 = t_i;
			  }
		  }
		else
		  {
			t_i = parameters[i];
		  }

		if (applyTimeOffset && offsettingPossible) t_i+=time_offset;

		if (t_i<=_m.MPD.MetaInfo[0].t.value) 
		  {
			// form penalty
			negLogTruncationLikelihood += MathUtil.sqr(_m.MPD.MetaInfo[0].t.value - t_i);
			// truncate??
			if (!reParameteriseToAdditiveTimes || i==1)
			  {
				parameters[i] = _m.MPD.MetaInfo[0].t.value;	
				if (applyTimeOffset && offsettingPossible) parameters[i]-=time_offset;
			  }
			else
			  {
				parameters[i] =0.0;
			  }
			truncated = true;
			if (_m._SI.verbosity>4) _m._SI.EH.println("Upper T truncation imposed!");
		  }
		if (t_i>=_m.MPD.t_base.value)
		  {
			// form penalty
			negLogTruncationLikelihood += MathUtil.sqr(_m.MPD.t_base.value - t_i);
			// truncate??
			if (reParameteriseToAdditiveTimes)
			  {
				parameters[i] =0.0;
			  }
			else
			  {
				parameters[i] = _m.MPD.t_base.value; 	
				if (applyTimeOffset && offsettingPossible) parameters[i]-=time_offset;
			  }
			truncated = true;
			if (_m._SI.verbosity>4) _m._SI.EH.println("Bottom T truncation imposed!");
		  }
	  }
	return truncated;
  }

  /**
	 method to detect if two layers are looking to have the same layer times (a pinchout), 
	 called for each new point in the minimiser outer loop.
  */
  private void  countTimeLatchingOccurrences()
  {
	int i;
	byte[] latch_sense;
	boolean[] couldBeLatched;
	
	latch_sense = new byte[n+1];
	Arrays.fill(latch_sense,(byte)0);
	couldBeLatched = new boolean[n+1];
	Arrays.fill(couldBeLatched,false);

	for (i=1;i<n;i++) 
	  {
		if (Math.abs(times[i] - times[i+1]) <= _m.MS.TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER)
		  {
			couldBeLatched[i] = true;
			latch_sense[i] = (times[i+1] > times[i]) ? (byte)1 : (byte)(-1);
		  }
		if (couldBeLatched[i])
		  {
			// do these tests only if we're not latched already
			if (layerIsPinchedOut[i]==0)			
			  {
				// These times {t[i+1], t[i]} are latched just as last time...
				if (latch_sense_last_evaluation[i] == latch_sense[i])
				  {
					// ...and with same sense
					// Keep counting...
					latch_count[i]++;
				 
					// Trap the solver by setting the flag layerIsPinchedOut if it looks like we need a pinchout solution
					// Do this if this layer is not already pinched out
					if (latch_count[i] > LATCH_COUNT_THRESHOLD)
					  {
						if (_m._SI.verbosity>3)
						  {
							_m._SI.EH.println("Iteration "+it_count+", pinchout enforced on layer " + i);
						  }
						layerIsPinchedOut[i] = latch_sense[i];
						// Set flag to abort from the solver
						abortNow=true;
					  }
				  }
			  }
			latch_sense_last_evaluation[i] = latch_sense[i];
		  }
		else
		  {
			// Doesn't appear to be latched. Reset the counter if we're not trapped yet
			if (layerIsPinchedOut[i]==0)
			  {
				latch_count[i]=0;
				latch_sense_last_evaluation[i] = 0;
			  }
		  }
	  }
	couldBeLatched=null;
	latch_sense=null;
  }

  /**
	 Returns component of likelihood due to negLogIsopachLikelihood
	 to be used only after a successful f_to_minimize() call
  */
  public double getNegLogIsopachLikelihood()
  {
	return negLogIsopachLikelihood;
  }
  /**
	 Returns component of likelihood due to negLogDepthLikelihood
	 to be used only after a successful f_to_minimize() call
  */
  public double getNegLogDepthLikelihood()
  {
	return negLogDepthLikelihood;
  }
  /**
	 Returns component of likelihood due to negLogVerticalStackingLikelihood
	 to be used only after a successful f_to_minimize() call
  */
  public double getNegLogVStackLikelihood()
  {
	return negLogVStackLikelihood;
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

  // test for debugging
  boolean testForNaN(double x)
  {
	Double T = new Double(x);
	if(T.isNaN())
	  {
		return true;
	  }
	T=null;
	return false;
  }
}
