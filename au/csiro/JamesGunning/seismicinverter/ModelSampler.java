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

import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.Minpack.optimization.*;
import au.csiro.JamesGunning.lib.util.MathUtil;

import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.stat.util.MCMCConvergenceMonitor;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.Numerical.*;


import java.util.*;
import static java.lang.Math.*;
import cern.jet.math.Functions;
import cern.jet.math.*;
import cern.colt.matrix.*;

import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix1D;


import cern.colt.list.*;
import cern.colt.matrix.linalg.*;
import cern.jet.random.*;

import cern.jet.random.engine.MersenneTwister;

/**
 Mapping from layer properties to smallest possible set for the current trace,
 for all fluid combinations, incorporating  synchronication etc
*/

public final class ModelSampler
{
  /**
	  pointer to underlying model
  */
  Model _m;

  /*
	  shortcut to seismic inverter too
  */
  SeismicInverter _SI;

  /**
	array of hashMaps of indices for properties
	one for each fluid state
  */
  HashMap map_to_index[];
  cern.colt.map.OpenIntObjectHashMap map_from_index[];

  /**
	 Array of master layer numbers needed if layer suppresion is invoked
	 time_synch_master_layer_number[i] is the master layer number of layer i.
	 If time_synch_master_layer_number[i] != i the layer is time-synched.
  */
  int[] time_synch_master_layer_number;

  /**
	 Array of master layer numbers.
	 master_layer_number_sequence[0...length-1] is the sequence of 
	 independent master layers
  */
  int[] master_layer_number_sequence;

  /**
	 Arrays of Properties to be mapped from a certain
	 index value when synchronisation applies
  */
  /**
	 time of layer
  */
  ArrayList t_array[];
  /**
	  model depth of layer
  */
  ArrayList d_array[];

  /**
	  rock_curves depth of layer
  */
  ArrayList d_curves_array[];

  /**
	Low freq interval velocity
  */
  ArrayList LFIV_array[];

  /**
	 Net-to-gross (fraction of "sand", or
	 `sand-end' rock member, susceptible of fluid substitution)
  */
  ArrayList NG_array[];
  /**
	 sand porosity
  */
  ArrayList phi_s_array[];
  /**
	 mud density ("mud" is shale/impermeable rock)
  */
  ArrayList rho_m_array[];
  /**
	pure sand p-wave velocity
  */
  ArrayList vp_s_array[];
  /**
	 pure sand s-wave velocity
  */
  ArrayList vs_s_array[];
  /**
	pure mud p-wave velocity
  */
  ArrayList vp_m_array[];
  /**
	 pure mud s-wave velocity
  */
  ArrayList vs_m_array[];


  // Hard rock properties
  ArrayList NG_hard_array[];
  ArrayList vp_hard_array[];
  ArrayList vs_hard_array[];
  ArrayList rho_hard_array[];


  /**
	  Density of reference fluid (brine?)
  */
  ArrayList rho_b_array[];

  /**
	  Velocity of reference fluid (brine?)
  */
  ArrayList vp_b_array[];

  /**
	 Arrays used for synchronised fluid properties, if required
  */
  ArrayList vp_brine_array[];
  ArrayList vp_gas_array[];
  ArrayList vp_lsg_array[];
  ArrayList vp_oil_array[];
  ArrayList rho_brine_array[];
  ArrayList rho_gas_array[];
  ArrayList rho_lsg_array[];
  ArrayList rho_oil_array[];
  ArrayList S_brine_array[];
  ArrayList S_gas_array[];
  ArrayList S_lsg_array[];
  ArrayList S_oil_array[];

  /**
   float fraction in reservoir rock
  */
  ArrayList float_fraction_array[];

  /**
	 Seismic A_AVO array for synchronisation
  */
  ArrayList A_AVO_array[];

  /**
	 Flag to indicate success of optimisation pass
  */
  boolean noModesFound = false;
  /** 
	 Working array of prior covariance matrices for each
	 fluid state. Always has the pure prior
  */
  DoubleMatrix2D[] Cinv_prior_only;
  /**
	 Working array of prior covariance matrices for each
	 fluid state. Used for general work.
  */
  DoubleMatrix2D[] Cinv_pr;
  /**
	 Working array of covariance matrices for each
	 fluid state, combining prior plus isopachs
  */
  DoubleMatrix2D[] Cinv_prior_plus_isopachs;

  /**
	 Working array of prior mean
  */
  DoubleMatrix1D[] mbar_prior_only;
  /**
	 Working array of corresponding prior mean
  */
  DoubleMatrix1D[] mbar_prior_plus_isopachs;
  /**
	 Working array of corresponding Cholesky factor of prior inverse C
  */
  DoubleMatrix2D[] R_prior_only;
  /**
	 Working array of Cholesky right-factor of prior C (no isopachs)
	 or R_prior_plus_isopachs[fs]^T.R_prior_plus_isopachs[fs] = Ctilde[fs] after isopachs are included
  */
  DoubleMatrix2D[] R_prior_plus_isopachs;
  /**
	Working array of log (determinant of prior Covariance C matrix)
  */
  double[] logDetC_prior_only;

  /**
	Working array of determinants of prior Covariance C updated with isopachs
  */
  double[] detC_prior_plus_isopachs;

  /**
	 Working array of relative time-thickness uncertainties, first index
	 is fluid state, second is master layer number. These entries are the layer
	 time-thickness std deviation, divided by the mean time-thickness.
	 High values indicate  high probability of layer pinchout. The
	 values are set after prior and isopach inclusion, but before seismic
  */
  DoubleMatrix2D time_uncertainty;
  /**
	 Working array of master-layer time uncertainties, first index is fluid
	 state, second is layer number. These entries are the master-layer time
	 std deviation. The values are set after prior and isopach
	 inclusion, but before seismic
  */
  DoubleMatrix2D t_layer_uncertainty;


  /**
	 Working array of local modes for sampling from
  */
  ArrayList<LocalMode>localModes;

  /** 
	  Working timeFinding object for locating and sampling modes
  */
  OptimalTimesFinder timeFinder;

  /** 
	  Working global genetic timeFinding object for locating modes
  */
  GeneticOptimalTimesFinder gTimeFinder;

  /**
	 Useful list of approximate std. deviations of layer times, in order,
	 for each fluid state. The are n_layers + 1 entries.
	 Last entry in time variances is for t_base
  */
  TreeMap[] time_std_devs;

  /**
	 Boolean flag to overcome impossible Gassman substitutions by adding
	 a large penalty to overall objective
  */
  private boolean useGassmanPenaltyTrick = true;


  boolean pinchoutFreeSolutionLikely = false;

  /**
	 Value of penalty applied int case of Gassman failure
  */
  private double GassmanPenaltyValue = 1000;
  /**
	 When mode enumerating, after Newton-updating for seismic over all
	 parameters, modes with identical fluid combos and pinchout
	 configurations are compared by means of their QP_optimum
	 vectors. If, for every element in this vector,
	 |mode1.QP_optimum(j) - mode2.QP_optimum(j)| < MODE_COMPARISON_EPSILON * scale, 
	 they are declared to be identical, and the duplication is eliminated. The scale is set by
	 mode1.R(j,j)
  */
  public static double MODE_COMPARISON_EPSILON=0.5;

  /**
	 Range of mode negative log-posteriors to accept as valid modes (any
	 mode with neg-log-posterior greater than the best neg-log-posterior 
	 *plus* this number is eliminated)
	 This number is ln(100), so modes of less than 1% approximate 
	 probability wrt the most likely mode
	 are culled.
  */
  private static double MODE_NLL_ELIMINATION_MARGIN = 4.605;

  /**
	 Epsilon to detect optimiser locking into "groove" at NG=1 or 0 caused
	 by likelihood driving to an endpoint, but concave prior pulling it back outside [0,1]
  */
  private static double NG_EPSILON=0.005;


  /**
	 Mode finding always searches along lines from current point p to,
	 say, a Gauss-Newton point p'. If Gassman-substitution is impossible at p', the
	 method will go as far along the current search direction (joining
	 p and p') as Gassman is legal. Finding the point along this line
	 nearest to p' where Gassman is still legal is done by bisection,
	 with this limit on iterations
  */
  private static int GASSMAN_BISECTION_MAX_ITS = 40;

  /**
	 Mode finding always searches along lines from current point p to,
	 say, a Gauss-Newton point p'. If Gassman is impossible at p', the
	 method will go as far along the current search direction (joining
	 p and p') as Gassman is legal. Finding the point along this line
	 nearest to p' where Gassman is still legal is done by
	 bisection. Bisection stops when we are within this fractional
	 amount of the "Gassman-fails" edge, on the legal side. This
	 number should be comfortably larger than the DESIGN_MATRIX_EPS variable
	 below, since finite difference perturbations may be computed
	 around this point.
  */
  private static double GASSMAN_BISECTION_EPS = 0.01;

  /**
	 epsilon to trap when backtracking in mode hunting has become degenerate
  */
  private static final double BACKTRAP_EPSILON = 1e-6;

  /**
	 Maximum allowable successive failures for fluid substitution
	 Successive attempts that fail more than this number of
	 times almost certainly indicate unphysical rock or fluid data.
  */
  static public int MAX_MC_FAIL_COUNT = 100;

  /**
	 Standardised deviate threshold to determine whether a layer seems 
	 likely to pinchout, and therefore whether to seek a local mode
	 in the likelihood surface assuming time-starting points that impose
	 the pinchout of the layer in question. This is a threshold for the
	 layer time-thickness std deviation over the mean layer time-thickness.
	 If it's Gaussian, a factor of less than 0.33 indicates 99% chance of not
	 being pinched out. 0.5 Corresponds to 95%.
  */
  static public final double LAYER_PINCHOUT_THRESHOLD = 0.5;

  /**
	 This switch will try a full-parameter Gauss-Newton first, for all models, before
	 trying subspace BFGS. If it succeeds, the subspace BFGS is never called. If it fails, the
	 subspace BFGS is tried on subsequent calls.
  */
  static public final boolean ALWAYS_ATTEMPT_GAUSS_NEWTON_FIRST = true;

  /**
	 scaled seismic error mismatch threshold at which to start firing off expensive global optimisers, 
	 presuming the mismatch is due to being trapped in poor local minima when mode hunting
  */
  static public final double UNACCEPTABLE_SCALED_SEISMIC_MISMATCH_GLOBAL_OPTIMISER_THRESHOLD = 1.2;

  /**
	 Dimensionless value at which to start firing off expensive global optimisers:
	 xi is a measure of the strength of the prior relative to the seismic likelihood. 
	 A single-reflection problem with a Ricker-2 wavelet is unimodal for xi<0.5, assuming 
     an S/N ratio of 4:1, and a prior
	 guess close to the true global minima.
  */
  static public final double XI_THRESHOLD_FOR_MULTIMODALITY = 0.6;
  /**
	 Mode searcher will often try this many initial
	 combinations of values in order to get a good match
  */
  static public int N_RANDOM_STARTING_TRIALS = 50;

  /**
	 Mode searcher will retry this many optimised initial
	 combinations before giving up when seeking a mode from
	 fixed time values
  */
  static public int MAX_NUM_TRIAL_GUESSES = 10;

  /**
	 Mode searcher will retry this many Newton refinements
	 of a mode if faulty Newton convergence is detected. If still
	 faulty, we give up.
  */
  static public int MAX_NUM_FAULTY_NEWTON_REFINEMENTS = 10;

  /**
	 Threshold at which to make frequently-hammered sparse matrices dense, 
	 for efficiency reasons.
  */
  static public double SPARSITY_THRESHOLD = 0.05;

  /**
	 Difference in mode max. posterior values less than this (for the same
	 fluid type) is taken as evidence the modes are in fact the same,
	 and they are amalgamated. This kind of energy barrier is easily
	 overcome in diffusion (easy tunneling)
  */
  static public double MODE_SEPARATION_EPSILON = 0.5;


  // standard Gaussian and Uniform generators
  MersenneTwister ranEngine;
  AbstractDistribution Ngen;
  AbstractDistribution Ugen;



  // WORKING variables

  // working reference to current array of layer properties
  LayerProperties LPI=null;

  // working references for prior distributions
  DoubleMatrix2D Cinv = null;
  DoubleMatrix1D mu = null;

  // current sample of the reduced property vector when doing MCMC sampling
  double[] property_sample;

  // current values of scaled isopach errors working vector
  private double[] isopach_errors;

  // flag to signal application of depth consistency
  boolean applyDepthConsistency = true;

  // current values of scaled consistent-depth errors working vector
  private double[] consistent_depth_errors;

  // current values of scaled vertical-stacking errors working vector
  double[] vertical_stacking_errors;

  // constant used for scaling depth-self-consistent likelihood function
  double sigma_consistent_depth = 1;

  // current local mode in multimode MCMC sampling
  LocalMode mode=null;

  // fluid state working variable
  int fs;  
  /**
	 Fluid state working variable
	 Never alter the value of this variable, but it
	 can be read to determine the current configuration of fluids
  */
  FluidStates.FluidState fluid_state=null;

  // working maps only, re-referenced for each fs
  HashMap map_to;
  cern.colt.map.OpenIntObjectHashMap map_from;
  // array versions of map_from, for speed.
  LayerProperties.Property[][] map_from_arrays;

  // variable to recognise Synchronised property ArrayLists
  static private Class arraylistClass=null;

  // working MutableInteger
  MutableInteger mutI = new MutableInteger();

  /**
	 Epsilon used to multiply depth of last layer to obtain sensible 
	 sigma_consistent_depth value used in depth-consistency likelihood
  */
  private static final double  EPSILON_FOR_DEPTH_SELF_CONSISTENCY = 1e-3;

  /**
	 Default sigma_consistent_depth value used in depth-consistency likelihood, used
	 only when scaling seems to fail
  */
  private static final double  DEFAULT_SIGMA_DEPTH_SELF_CONSISTENCY = 1;

  /**
	 Time epsilon used to determine that two times are the same when
	 comparing minima. try 0.1 ms
  */
  private static final double  TIME_EPSILON_FOR_MINIMA_MATCH 
	= 0.1*FundamentalScalingConstants.BASE_SCALING_TIME;    // usually 0.1 ms

  /**
	 Time epsilon used to lock layer times together in pinchout situation.
	 Must be less than TIME_EPSILON_FOR_MINIMA_MATCH, but not too much
  */
  static final double  TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER = 0.5*TIME_EPSILON_FOR_MINIMA_MATCH;  // usually 0.05 ms

  /**
	 Single floating point epsilon small for typical times
  */
  static final double  DEPSILON = 1e-7*FundamentalScalingConstants.BASE_SCALING_TIME;

  /**
	 Epsilon to be used in finite differences. Forward
	 difference suggest this ought to be O(sqrt(10*f_accuracy),
	 here we will imagine 12 sig. figures in the fwd model, which allows 
	 (pessimistically) for  loss of 4 figures
	 from numerical conspiracies, if eps_double=10^-16.
	 (All properties are perturbed by this fractional amount, or this amount
	 absolutely if the property is 0 or 1)
  */
  private static final double DESIGN_MATRIX_EPS = 1e-6;

  // Make sure time perturbations are sane etc
  private static final double SEISMIC_DESIGN_MATRIX_EPS = DESIGN_MATRIX_EPS;

  private static final double SEISMIC_DESIGN_MATRIX_TIME_EPS = DESIGN_MATRIX_EPS*FundamentalScalingConstants.BASE_SCALING_TIME;
  private static final double TEN_MACHINE_EPS = 1e-15;


  /**
	 State variable to indicate dispoaition of pinchouts
  */

  boolean[] layer_present;
  /**
	 Binary encoding of exposed layers as a state variable
  */
  int layer_code;

  int _verbosity;

  // generally useful flag for overhead communications
  private boolean abortNow=false;

  /**
	 RandomWalkSampler to be used if desired
  */
  RandomWalkSampler randomWalkSampler;

  /**
	 IndependenceSampler to be used if desired
  */
  IndependenceSampler independenceSampler;


  // Variable for Newton optimisations
  /**
	 Maximum number of Newton steps used in Newton methods of computing the posterior 
	 covariance and final solution
  */
  public static int MAX_NEWTON_STEPS = 10;
  /**
	 Epsilon for Newton method terminations. Usually f=O(d), and we stop optimising when
	 |delta f | < 0.0001 |f|. This is the 0.0001
  */
  public static double NEWTON_EPS = 0.0001;
  /**
	 There is an epsilon for termination of Newton methods associated with floating point
	 rounding too
  */
  public static double EPS_FLOATINGPOINT=1e-8;

  private double fNewNewton; // value at new point
  private double fNewNewtonOld; // value at old point
  private double delta_fNewNewton; // change from old point to new point

  /**
	 ConvergenceMonitor object for managing multiple runs
	 and monitoring convergence
  */
  MCMCConvergenceMonitor CM;

  /**
	 create instance of a Model, details to be filled in with setUpFully() when known
  */

  /**
	 seismic chisq error variable to be used or accessed by the samplers
  */
  double negLogLikSeisOnly;


  /**
	 hidden seismic chisq error variable to be used or accessed by the samplers
  */
  double negLogLikHiddenSeisOnly;


  /** 
	  useful terms for marginal model likelihood
  */
  double logDetCisopachAnd2PiPieces=0.0;
  double logDetCvsAnd2PiPieces=0.0;


  // Useful flag to communicate status of mode searches and reasons for failure
  // 0 = success
  // 1 = failure due to Gassman illegality
  // 2 = failure due to Newton method not converging
  private int modeSearchSuccessFlag;  

  //private int debugCntr=0;

  ModelSampler(Model m)
  {
	this._m = m;
	this._SI = m._SI;
	_verbosity = m._SI.verbosity;

	t_array = new ArrayList[m.n_layers];
	d_array = new ArrayList[m.n_layers];
	d_curves_array = new ArrayList[m.n_layers];
	LFIV_array = new ArrayList[m.n_layers];
	NG_array = new ArrayList[m.n_layers];
	phi_s_array = new ArrayList[m.n_layers];
	rho_m_array = new ArrayList[m.n_layers];
	vp_s_array = new ArrayList[m.n_layers];
	vs_s_array = new ArrayList[m.n_layers];
	vp_m_array = new ArrayList[m.n_layers];
	vs_m_array = new ArrayList[m.n_layers];
	vp_b_array = new ArrayList[m.n_layers];
	rho_b_array = new ArrayList[m.n_layers];


	if (_SI.MF.ENABLE_HARD_ROCK_MODEL) // FIXME Save mem if possible?
	  {
		NG_hard_array = new ArrayList[m.n_layers];
		rho_hard_array = new ArrayList[m.n_layers];
		vp_hard_array = new ArrayList[m.n_layers];
		vs_hard_array = new ArrayList[m.n_layers];
	  }
	if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		float_fraction_array = new ArrayList[m.n_layers];
	  }

	A_AVO_array = new ArrayList[1];  // for A and B if needed

	// we need a Class identifier to detect synchronised properties
	try
	  {
		arraylistClass = Class.forName("java.util.ArrayList");
	  }
	catch (ClassNotFoundException ex) {}; // never happens

	// deviate generators
	if (_SI.useSeed && (_SI.seed < 0))
	  {
		// use System nano time, randomising seed across threads or processes...
		ranEngine = new MersenneTwister((int)System.nanoTime());
	  }
	else
	  {
		ranEngine = new MersenneTwister(_SI.seed);
		// increment seed to prevent blatant systematics across traces (but retain determinism)
		_SI.seed++;
	  }
	Ngen = new Normal(0,1.0,ranEngine);
	Ugen = new Uniform(0,1.0,ranEngine);

	// pinchout state array
	layer_present = new boolean[m.n_layers];

	// smp blas. Have to know the number of cpus.
	int cpu_s;
	// NB: as of java 1.4, we can get this with Runtime.availableProcessors()
	// cpu_s = Runtime.getRuntime().availableProcessors();
	// for now, leave it at 1.
	cpu_s = 1;
	SmpBlas.allocateBlas(cpu_s, SeqBlas.seqBlas);
  }

  void completeSetUpUsingFluidContactInfo()
  {
	// Allocatons for possibly synchronised fluids...
	if (_SI.doSimpleFluidUnitsAndContacts)
	  {
		vp_brine_array = new ArrayList[_m.n_layers];
		vp_gas_array = new ArrayList[_m.n_layers];
		vp_lsg_array = new ArrayList[_m.n_layers];
		vp_oil_array = new ArrayList[_m.n_layers];
		rho_brine_array = new ArrayList[_m.n_layers];
		rho_gas_array = new ArrayList[_m.n_layers];
		rho_lsg_array = new ArrayList[_m.n_layers];
		rho_oil_array = new ArrayList[_m.n_layers];
		S_brine_array = new ArrayList[_m.n_layers];
		S_gas_array = new ArrayList[_m.n_layers];
		S_lsg_array = new ArrayList[_m.n_layers];
		S_oil_array = new ArrayList[_m.n_layers];
	  }
  }

  void setUpLayerSynching()
  {
	int i,j,k;
	boolean is_time_synched;
	// pre-initialise master layer numbers
	time_synch_master_layer_number = new int[_m.n_layers];
	master_layer_number_sequence = new int[_m.n_active_layers];
	for (i=0,j=0,k=0;i<_m.n_layers;i++)
	  {
		is_time_synched = (i==0) ? false : !_m.MPD.MetaInfo[i-1].isActive;
		if (!is_time_synched) 
		  {
			j=i;
			master_layer_number_sequence[k++] = i; // list of master layers
		  }
		time_synch_master_layer_number[i]=j;       // master layer of layer i
	  }
  }

  /**
	 set up mappings for true free variables to linear algebra indices
  */
  void setUpMappings()
  {
	int time_synch_indx=0;
	short i_h,N_possible_fluids;
	boolean is_time_synched=false;

	int i,j,j_s,indx,t_indx,synch_indx, tmp_indx;
	ModelPriorDescriptions.MetaInfo MI,MI_synch,MI_tmp;
	ModelPriorDescriptions.PriorDetails P_h = null;
	ModelPriorDescriptions.PriorDetails PD_sat = null;
	double sigma_S_h = 0;

	RockandFluidData.RegressionSet CS=null;
	RockandFluidData.MaterialProperties fl_props=null;
	RockandFluidData.MaterialProperties hydro_props=null;
	LayerProperties.Property S_h = null;
	LayerProperties.Property rho_h = null;
	LayerProperties.Property vp_h = null;
    String hydro_name = null;
	boolean is_synched=false, must_synch = false;
	boolean overwriteVpFl,overwriteRhoFl;
	String name=null;
	short fluid_type;

	// Fluid contacts
	FluidUnitsAndContacts.FluidUnit fu;
	int fluid_synch_indx = 0;
	boolean fluid_synched_image = false;
	FluidUnitsAndContacts _FUC;
	// references...
	ArrayList[] vp_h_array=null;
	ArrayList[] rho_h_array=null;
	ArrayList[] S_h_array=null;

	_FUC = _m.fluidUnitsAndContacts;
	if (_m._SI.doSimpleFluidUnitsAndContacts)
	  {
		_m.fluidUnitsAndContacts.fc_props = new ArrayList<LayerProperties.Property>(0);
	  }
	// Seismic AVO array common to all models
	A_AVO_array[0] = new ArrayList(1);

	// Hashmaps for getting/retreiving linear algebra indices from properties etc
	map_to_index = new HashMap[_m.FS.possibleFluidStates.size()];
	map_from_index = new cern.colt.map.OpenIntObjectHashMap[_m.FS.possibleFluidStates.size()];

	setUpLayerSynching();
	updateIsPermeableFlagAndNIsopachs();  // might not be needed here, but fast & harmless I think
	_FUC.updateMasterFluidLayerTable();


	// loop over fluid combinations
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		map_to_index[fs] = new HashMap();
		map_from_index[fs] = new cern.colt.map.OpenIntObjectHashMap();

		map_to = map_to_index[fs];
		map_from = map_from_index[fs];

		if (_verbosity>3) _SI.EH.print("%nVariable ordering for fluid state " + fs + "%n");

		indx=0;   // will serve as indexing for non-time variables
		t_indx=0; // will serve as index for time variables
		// count number of time variables required (almost certainly _m.n_layers+1 unless some layers are 
		// pinchouted by the flag sigma_t==0
		indx = _m.n_active_layers+1;

		for (i=0;i<_m.n_layers;i++)
		  {
			LPI = _m.LP[i];
			MI = _m.MPD.MetaInfo[i];
			
			boolean isMasterLayer = (i==(_m.master_depth_layer_number-1));
			
			// zero thickness imposed by synching bounding times
			// if layer above is of zero thickness (not Active), the times must be synched
			is_time_synched = (i==0)? false : !_m.MPD.MetaInfo[i-1].isActive;

			// only unsynchronised rock properties get to enter the reduced model
			// but we will keep lookups to stuff since it's handy
			if (is_time_synched)
			  {
				time_synch_indx = time_synch_master_layer_number[i];
			  }
			if (_m.MPD.synchMap!=null && (_m.MPD.synchMap.containsKey(i)))
			  {
				// layer number of master layer (this layer is an image, except for fluids)
				synch_indx = _m.MPD.synchMap.get(i);
				is_synched = true;
			  }
			else
			  {
				is_synched = false;
				synch_indx = i;
			  }
			// Set up new set of ArrayLists of synchronised layers if required
			// These are the same for all fluid combos, so only do it once
			if (fs==0) setUpSynchPropertyArrays(MI,i,is_time_synched,is_synched);


			// layer to set properties from
			MI_synch = _m.MPD.MetaInfo[synch_indx];

			// Get indexes for properties etc

			// times if variable
			LPI.t.setDistr(MI.t.value,MI.sigma_t.value);
			if (!is_time_synched)
			  {
				// add new master time 
				// addTerm(t_indx,LPI.t);
				addSynchTerm(false, i, t_indx, LPI.t, t_array);
				t_indx++;
			  } 
			else
			  {
				// add synchronised time
				addSynchTerm(true, time_synch_indx, t_indx, LPI.t, t_array);
			  }


			// Model depth. If synchronised, synchronised to master layer number
			tmp_indx = MI.depth.master_layer_number;
			must_synch = (MI.depth.master_layer_number != i);
			MI_tmp = _m.MPD.MetaInfo[MI.depth.master_layer_number];
			LPI.d.setDistr(MI_tmp.depth.value,MI_tmp.sigma_depth.value);
			if (MI_tmp.sigma_depth.value > 0 && (MI.isActive || MI.hasDepthDependents || isMasterLayer))
			  {
				addSynchTerm(must_synch, tmp_indx, indx, LPI.d, d_array);
				if (!must_synch) indx++;
			  }

			// Possibly distinct rock_curves depth. If distinct, it has sigma=0, so
			// code chunk is as above, with sigma <- 0
			// If not distinct, Ensure LPI.d_curves references LPI.d.
			if (MI.depth != MI.depth_rock_curves)
			  {
				tmp_indx = MI.depth_rock_curves.master_layer_number;
				must_synch = (MI.depth_rock_curves.master_layer_number != i);
				MI_tmp = _m.MPD.MetaInfo[MI.depth_rock_curves.master_layer_number];
				LPI.d_curves.setDistr(MI_tmp.depth_rock_curves.value, 0);
			  }
			else
			  {
				// Ensure LPI.d_curves references LPI.d.
				LPI.d_curves = null; // GC
				LPI.d_curves = LPI.d;
			  }


			// LFIV: if synchronised, synchronised to master layer number
			tmp_indx = MI.LFIV.master_layer_number;
			MI_tmp=_m.MPD.MetaInfo[tmp_indx];
			must_synch = (tmp_indx != i);
			LPI.LFIV.setDistr(MI_tmp.LFIV.value, MI_tmp.sigma_LFIV.value);

			if (MI_tmp.sigma_LFIV.value > 0 && (MI.isActive || MI.hasLFIVDependents))
			  {
				addSynchTerm(must_synch, tmp_indx, indx, LPI.LFIV, LFIV_array);
				if (!must_synch) indx++;
			  }
			
			// MI.is_permeable updated as above...
			
			if (MI.isActive || (_m.MPD.synch_copies!=null && _m.MPD.synch_copies[i]>0))
			  {
				// impermeable rocks...(must not be pure permeable rock mixture)
				if (!MI.is_permeable || !(MI.sigma_NG.value==0 && MI.NG.value == 1.0))
				  {
					// impermeable rock properties. May be synchronised
					CS = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(MI.mud_rock_member);

					if (CS.Vp_versus_depth.error > 0)  // Vp needed
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.vp_m, vp_m_array);
						if (!is_synched) indx++;
					  }
					if (CS.Vs_versus_Vp.error > 0)  // Vs needed
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.vs_m, vs_m_array);
						if (!is_synched) indx++;
					  }
					if (CS.Rho_or_Phi_versus_Vp.error > 0)  // rho needed
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.rho_m, rho_m_array);
						if (!is_synched) indx++;
					  }
				  }

				// contributions from hard rock
				if (_SI.MF.ENABLE_HARD_ROCK_MODEL)
				  {
					if (_m.needHardRockTerms) // prelim check to speed things up
					  {
						boolean isPureSand = (MI_synch.NG.value==1.0) && (MI_synch.sigma_NG.value==0.0);
						if (MI.has_hard_rock_mixture && !isPureSand && !(MI_synch.sigma_NG_hard.value==0 && MI_synch.NG_hard.value == 0.0))
						  {
							// Hard-rock NG
							LPI.NG_hard.setDistr(MI_synch.NG_hard.value,MI_synch.sigma_NG_hard.value);
							if (MI.sigma_NG_hard.value > 0)
							  {
								addSynchTerm(is_synched, synch_indx, indx, LPI.NG_hard, NG_hard_array);
								if (!is_synched) indx++;
							  }

							// impermeable rock properties. May be synchronised
							if (MI.hard_shale_rock_member != null)
							  {
								// fixed rock type
								CS = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(MI.hard_shale_rock_member);
							  }
							else // presume  fluid_dep_hard_shale_rock_members!=null by construction
							  {
								fluid_type = ((FluidStates.FluidState)_m.FS.possibleFluidStates.get(fs)).types[i];
								String s = MI.fluid_dep_hard_shale_rock_members[fluid_type];
								if (s==null)
								  {
									throw new RuntimeException("Required hard rock type, layer i="+(i+1)+", "+MI.layer_name+
															   " doesn't exist for fluid state fs="+fs+" layer_fluid_type="+fluid_type);
								  }
								CS = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(s);
								if (CS==null)
								  {
									throw new RuntimeException("Required hard rock type "+s+", layer i="+(i+1)+", "+MI.layer_name+
															   "not furnished in XML: fluid state fs="+fs+" layer_fluid_type="+fluid_type);
								  }
							  }
						
							if (CS.Vp_versus_depth.error > 0)  // Vp needed
							  {
								addSynchTerm(is_synched, synch_indx, indx, LPI.vp_hard, vp_hard_array);
								if (!is_synched) indx++;
							  }
							if (CS.Vs_versus_Vp.error > 0)  // Vs needed
							  {
								addSynchTerm(is_synched, synch_indx, indx, LPI.vs_hard, vs_hard_array);
								if (!is_synched) indx++;
							  }
							if (CS.Rho_or_Phi_versus_Vp.error > 0)  // rho needed
							  {
								addSynchTerm(is_synched, synch_indx, indx, LPI.rho_hard, rho_hard_array);
								if (!is_synched) indx++;
							  }
						  }
					  }
				  }

				// contributions from permeable rock and fluids........
		
				if (MI.is_permeable)
				  {
					// Sand & maybe fluid properties will enter the state vector

					LPI.NG.setDistr(MI_synch.NG.value,MI_synch.sigma_NG.value);
					if (MI.sigma_NG.value > 0)
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.NG, NG_array);
						if (!is_synched) indx++;
					  }
					  

					// floating grain fraction if required/enabled
					if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
					  {
						// floating grain if synchronised, synchronised to master layer number
						tmp_indx = MI.float_fraction.master_layer_number;
						MI_tmp=_m.MPD.MetaInfo[tmp_indx];
						must_synch = (tmp_indx != i);
						LPI.float_fraction.setDistr(MI_tmp.float_fraction.value, MI_tmp.sigma_float_fraction.value);
						if (MI_tmp.sigma_float_fraction.value > 0 && (MI.isActive || MI.hasFloat_fractionDependents))
						  {
							addSynchTerm(must_synch, tmp_indx, indx, LPI.float_fraction, float_fraction_array);
							if (!must_synch) indx++;
						  }
					  }


					// permeable rock properties
					CS = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(MI.sand_rock_member);

					if (CS.Vp_versus_depth.error > 0)  // error on Vp curve
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.vp_s, vp_s_array);
						if (!is_synched) indx++;
					  }
					if (CS.Vs_versus_Vp.error > 0)  // error on Vs curve
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.vs_s, vs_s_array);
						if (!is_synched) indx++;
					  }
					if (CS.Rho_or_Phi_versus_Vp.error > 0)  // error on rho curve
					  {
						addSynchTerm(is_synched, synch_indx, indx, LPI.phi_s, phi_s_array);
						if (!is_synched) indx++;
					  }

					// fluid terms. Compute these, but add to the model only if layer is active
					// (They can never be a synch master)

					
					// Brine fluid density will enter even if there are no hydros

					// Brine Density.
					// (we model phi_s, and overall density has a phi_s * rho_b term)
					// if there is no areal variation in mean or std deviation, 
					// set distribution from fluid properties entries
					fl_props = (RockandFluidData.MaterialProperties)
					  _m.RFD.FluidTypes.get(MI.brine_name);
				
					LPI.rho_b.setDistr((MI.rho_b==null) ? fl_props.rho : MI.rho_b.value, 
									   (MI.sigma_rho_b==null) ? fl_props.sigma_rho : MI.sigma_rho_b.value); 
					if (LPI.rho_b.sigma>0 && MI.isActive)
					  {
						addTerm(indx,LPI.rho_b);
						indx++;
					  }
					// Brine Velocity
					// this is needed b/c The K_fluid = rho_b vp_b^2 enters Gassman
					// and *implies* a certain Kdry. Some Kdry values are illegal, so a certain
					// subspace of the brine {rho_b,vp_b} space is illegal and has to be caught
					LPI.vp_b.setDistr((MI.vp_b==null) ? fl_props.vp : MI.vp_b.value, 
									  (MI.sigma_vp_b==null) ? fl_props.sigma_vp : MI.sigma_vp_b.value); 
					if (LPI.vp_b.sigma>0 && MI.isActive)
					  {
						addTerm(indx,LPI.vp_b);
						indx++;
					  }

					// Here either 1 or 3 possible hydrocarbon cases, 
					// for enumeration-mode, or contacts-mode respectively
					// Here 3 = FluidStates.FluidTypes.Nfluids-1
					N_possible_fluids = _SI.doSimpleFluidUnitsAndContacts ? (short)3: (short)1;
					for (i_h=1; i_h<=N_possible_fluids; i_h++) // either {1}, or {1,2,3} = { oil, lsg, gas}
					  {
						overwriteRhoFl = overwriteVpFl = false;
						if (_SI.doSimpleFluidUnitsAndContacts)
						  {
							fluid_type = i_h;
							if (!_FUC.fluidRequired[i][i_h]) continue; // skip rest of loop
						  }
						else
						  {
							// fluid required by current fluid state
							fluid_type = ((FluidStates.FluidState)_m.FS.possibleFluidStates.get(fs)).types[i];
						  }
						switch (fluid_type)
						  {
						  case FluidStates.FluidTypes.brine:
							break;

						  case FluidStates.FluidTypes.oil:
							P_h = MI.P_oil;
							if (P_h.value>0 || _SI.doSimpleFluidUnitsAndContacts)
							  {
								sigma_S_h = (MI.sigma_S_oil!=null) ? MI.sigma_S_oil.value : 0;
								LPI.S_oil.setDistr(MI.S_oil.value, sigma_S_h);
								if (MI.vp_oil_VA) 
								  {
									overwriteVpFl  = true; 
									LPI.vp_oil.setDistr(MI.vp_oil.value, (MI.sigma_vp_oil!=null) ? MI.sigma_vp_oil.value : 0 );    // JGCO2 
								  }
								if (MI.rho_oil_VA) 
								  {
									overwriteRhoFl  = true; 
									LPI.rho_oil.setDistr(MI.rho_oil.value, (MI.sigma_rho_oil!=null) ? MI.sigma_rho_oil.value : 0 );    // JGCO2 
								  }
								S_h = LPI.S_oil;
								rho_h = LPI.rho_oil;
								vp_h = LPI.vp_oil;
								name = MI.oil_name;
								vp_h_array = vp_oil_array;
								rho_h_array = rho_oil_array;
								S_h_array = S_oil_array;
							  }
							break;

						  case FluidStates.FluidTypes.lsg:
							P_h = MI.P_lsg;
							if (P_h.value>0 || _SI.doSimpleFluidUnitsAndContacts)
							  {
								sigma_S_h = (MI.sigma_S_lsg!=null) ? MI.sigma_S_lsg.value : 0;
								LPI.S_lsg.setDistr(MI.S_lsg.value, sigma_S_h);
								if (MI.vp_lsg_VA) 
								  {
									overwriteVpFl  = true; 
									LPI.vp_lsg.setDistr(MI.vp_lsg.value, (MI.sigma_vp_lsg!=null) ? MI.sigma_vp_lsg.value : 0 );    // JGCO2 
								  }
								if (MI.rho_lsg_VA) 
								  {
									overwriteRhoFl  = true; 
									LPI.rho_lsg.setDistr(MI.rho_lsg.value, (MI.sigma_rho_lsg!=null) ? MI.sigma_rho_lsg.value : 0 );    // JGCO2 
								  }
								S_h = LPI.S_lsg;
								rho_h = LPI.rho_lsg;
								vp_h = LPI.vp_lsg;
								name = MI.lsg_name;
								vp_h_array = vp_lsg_array;
								rho_h_array = rho_lsg_array;
								S_h_array = S_lsg_array;
							  }
							break;

						  case FluidStates.FluidTypes.gas:
							P_h = MI.P_gas;
							if (P_h.value>0 || _SI.doSimpleFluidUnitsAndContacts)
							  {
								sigma_S_h = (MI.sigma_S_gas!=null) ? MI.sigma_S_gas.value : 0;
								LPI.S_gas.setDistr(MI.S_gas.value, sigma_S_h);
								if (MI.vp_gas_VA) 
								  {
									overwriteVpFl  = true; 
									LPI.vp_gas.setDistr(MI.vp_gas.value, (MI.sigma_vp_gas!=null) ? MI.sigma_vp_gas.value : 0 );    // JGCO2 
								  }
								if (MI.rho_gas_VA) 
								  {
									overwriteRhoFl  = true; 
									LPI.rho_gas.setDistr(MI.rho_gas.value, (MI.sigma_rho_gas!=null) ? MI.sigma_rho_gas.value : 0 );    // JGCO2 
								  }
								S_h = LPI.S_gas;
								rho_h = LPI.rho_gas;
								vp_h = LPI.vp_gas;
								name = MI.gas_name;	
								vp_h_array = vp_gas_array;
								rho_h_array = rho_gas_array;
								S_h_array = S_gas_array;
							  }
							break;
						  }

						if (fluid_type != FluidStates.FluidTypes.brine)
						  {
							// fluid substitution can occur if the P(hydro) > 0 ?
							// ref fluid uncertainties can also only enter if hydro's are possible
							if (P_h.value>0 || _SI.doSimpleFluidUnitsAndContacts)
							  {
								// Ditto for hydrocarbon
								hydro_props = (RockandFluidData.MaterialProperties)
								  _m.RFD.FluidTypes.get(name);
								boolean fluidsSynched = false;
								fluid_synched_image = false;
								if (_SI.doSimpleFluidUnitsAndContacts && _FUC.fluidsSynchedForFUowningLayer(i))
								  {
									fluidsSynched = true;
									fu = _FUC.getFUofThisLayer(i);
									fluid_synch_indx = fu.masterFluidLayer[i_h];
									if (fluid_synch_indx < 0)
									  {
										// should never happen but be safe...
										_SI.EH.printf("WARNING: No intelligible master layer for "+
													  "fluid %d in layer %d! Removing synchronisation for safety!",i_h,i+1);
										fluid_synch_indx = i;
									  }
									fluid_synched_image = (fluid_synch_indx!=i); // true on a fluid "image" layer
								  }
								// vp for hydrocarb
								// use fluid defaults if not overriden before  JGCO2
								if (!overwriteVpFl) vp_h.setDistr(hydro_props.vp, hydro_props.sigma_vp);

								if (MI.isActive && (hydro_props.sigma_vp>0 || fluid_synched_image))
								  {
									if (fluidsSynched)
									  { 
										handleFluidSynchArrayAllocs(fluid_type,fluid_synch_indx,'v');
										addSynchTerm(fluid_synched_image, fluid_synch_indx, indx, vp_h, vp_h_array);
										if (!fluid_synched_image) indx++; 
									  }
									else
									  {
										addTerm(indx,vp_h);
										indx++;
									  }
								  }
								// density for hydrocarb
								if (!overwriteRhoFl) rho_h.setDistr(hydro_props.rho, hydro_props.sigma_rho);
								if (MI.isActive && (hydro_props.sigma_rho>0 || fluid_synched_image))
								  {
									if (fluidsSynched)
									  { 
										handleFluidSynchArrayAllocs(fluid_type,fluid_synch_indx,'r');
										addSynchTerm(fluid_synched_image, fluid_synch_indx, indx, rho_h, rho_h_array);
										if (!fluid_synched_image) indx++; 
									  }
									else
									  {
										addTerm(indx,rho_h);
										indx++;
									  }
								  }

								// Saturation of hydrocarb

								if (MI.isActive && (sigma_S_h>0 || fluid_synched_image))
								  {
									// if a hydrocarbon can possibly be there
									// with uncertain saturation, add S
									if (fluidsSynched)
									  { 
										handleFluidSynchArrayAllocs(fluid_type,fluid_synch_indx,'S');
										addSynchTerm(fluid_synched_image, fluid_synch_indx, indx, S_h, S_h_array);
										if (!fluid_synched_image) indx++; 
									  }
									else
									  {
										addTerm(indx,S_h);
										indx++;
									  }
								  }
							  }
						  }
					  } // possible fluids in layer, either 1 (enumeration mode), up to 3 (simple contacts mode)
				  } // if() contributions from permeable rocks and fluids
			  } // if () contribution from active -layer rocks or a rock-synch master
		  } // loop over layers

		// Finally, add terms for the base time.
		// base times if variable
		LPI.t_base.setDistr(_m.MPD.t_base.value,_m.MPD.sigma_t_base.value);
		if (_m.MPD.sigma_t_base.value > 0)
		  {
			addTerm(t_indx,LPI.t_base);
			t_indx++;
		  }


		// Add to the map any additional variables

		// (1) internal variables of the Seismic model that need to be stochastic
		_m.SM.A.setDistr(_m.SM.A_avo.value, _m.SM.sigma_A_avo.value);
		if (_m.SM.sigma_A_avo.value > 0)
		  {
			addSynchTerm(false,0,indx,_m.SM.A,A_AVO_array);
			indx++;
		  }
		if (_m.SM.A_B_avo_synchronised)
		  {
			_m.SM.B.setDistr(_m.SM.A_avo.value, _m.SM.sigma_A_avo.value);
		  }
		else
		  {
			_m.SM.B.setDistr(_m.SM.B_avo.value, _m.SM.sigma_B_avo.value);
		  }
		if (_m.SM.B.sigma > 0)
		  {
			if (!_m.SM.A_B_avo_synchronised) 
			  {
				addTerm(indx,_m.SM.B);
				indx++;				
			  }
			else
			  {
				addSynchTerm(true,0,indx,_m.SM.B, A_AVO_array);
				// if (!_m.SM.A_B_avo_synchronised) indx++;   // no need to increment
			  }
		  }

		// (2) Simple contacts parameters
		if (_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			for (i=0;i<	_m.fluidUnitsAndContacts.fluidUnits.length; i++)
			  {
				fu = _m.fluidUnitsAndContacts.fluidUnits[i];
				for (j=0; j< fu.contactInfo.length; j++)
				  {
					FluidUnitsAndContacts.ContactInfo CI = 	fu.contactInfo[j];
					CI.depth.setDistr(CI.depthPD.value, CI.sigma_depthPD.value);
					if (CI.sigma_depthPD.value > 0)
					  {
						addTerm(indx,CI.depth);
						indx++;
						// Set up useful ArrayList for the simple contact model for
						// sensitivity calculation (there will be only 1 fluid state)
						_m.fluidUnitsAndContacts.fc_props.add(CI.depth);
					  }
				  }
			  }
		  }
	  } // loop over fluid states
  }

  void updateIsPermeableFlagAndNIsopachs()
  {
	int i;
	ModelPriorDescriptions.MetaInfo MI;
	ModelPriorDescriptions _MPD = _m.MPD;

	_MPD.n_isopachs=0;
	// loop on layers
	for (i=0;i<_m.n_layers;i++)
	  {
		MI = _MPD.MetaInfo[i];
		
		// Determine the rock mixture: update permeable_flag
		MI.is_permeable = ((MI.sand_rock_member != null)
						   && (MI.sigma_NG.value>0 || MI.NG.value>0));
		if (MI.dz != null && !((MI.dz.value==0) && (MI.sigma_dz.value==0))) _MPD.n_isopachs++;
	  }
  }

  /**
	 Add a term to the index mappings that can never be synchronised
  */
  private void addTerm(int indx, LayerProperties.Property prop)
  {
	ArrayList v = new ArrayList(1);
	v.add(prop);
	v.trimToSize();
	map_from.put(indx,v);

	map_to.put(prop,new Integer(indx));
	if (_verbosity>3) printVariableOrdering(indx,prop,false);
  }

  private void printVariableOrdering(int indx, LayerProperties.Property prop,boolean is_a_synch_image)
  {
	String name = (String)_m.completePropertyNameMap.get(prop);
	if (name!=null)
	  {
		_SI.EH.printf("%2d %s",indx, name);
		_SI.EH.println(((is_a_synch_image)? " (synch-image)" : ""));
	  }
  }

  private String getPropertyName(LayerProperties.Property prop)
  {
	return (String)_m.completePropertyNameMap.get(prop);
  }


  /**
	 Add a term to the index mappings that may potentially be synchronised
	 (e.g. rock properties)
  */
  private void addSynchTerm(boolean is_synched, int synch_indx, int indx,
							LayerProperties.Property prop, ArrayList[] arr)
  {
	if (!is_synched)
	  {
		// only need to grow list when on first fluid combo:
		// Things that are synched are the same for all fluid combos
		if (fs==0) arr[synch_indx].add(prop);
		// But since the indx may vary with fluid combo,
		// we need fluid-combo-specific index mappings
		map_from.put(indx,arr[synch_indx]);
		map_to.put(prop,new Integer(indx));
		if (_verbosity>3) printVariableOrdering(indx,prop,false);
	  }
	else
	  {
		// This layer is synched. It's non-fluid properties mirror those of the master layer.
		// What's the master layer index? It's the index corresponding to the masterproperty,
		// which is the first in the list arr[]. Get the masterproperty, and
		// then we can find what it's original index in the linear algebra was via the
		// reverse lookup HashMap.
		ArrayList AL = arr[synch_indx];
		// It's just conceivable we have synchronised to a layer where this property is
		// fixed, but in this layer it isn't, so therefore we should bail out...
		if (AL.size()==0) return;
		// but otherwise, carry on
		LayerProperties.Property masterP = (LayerProperties.Property)(AL.get(0));
		Integer masterIndex = (Integer) map_to.get(masterP);
		// add this new property to the list of synchronised guys, but only for the first fluid state
		// as we explained above
		if (fs==0) arr[synch_indx].add(prop);
		// and add lookup index of image-property
		map_to.put(prop,masterIndex);
		if (_verbosity>3) printVariableOrdering(masterIndex.intValue(),prop,true);
	  }
  }


  /**
	 Suite of rock properties that may be synchronised between layers
  */
  private void setUpSynchPropertyArrays(ModelPriorDescriptions.MetaInfo MI,int i,boolean is_time_synched,boolean is_synched)
  {
	if (!is_time_synched) 
	  {
		t_array[i] = new ArrayList(1);
	  }
	if (!is_synched)
	  {
		// model depth
		d_array[i] = new ArrayList(1);
		// trend depth
		if (MI.depth != MI.depth_rock_curves)
		  {
			d_curves_array[i] = new ArrayList(1);
		  }
		else
		  {
			d_curves_array[i] = d_array[i];
		  }
		LFIV_array[i] = new ArrayList(1);
		NG_array[i] = new ArrayList(1);
		phi_s_array[i] = new ArrayList(1);
		rho_m_array[i] = new ArrayList(1);
		vp_s_array[i] = new ArrayList(1);
		vs_s_array[i] = new ArrayList(1);
		vp_m_array[i] = new ArrayList(1);
		vs_m_array[i] = new ArrayList(1);
		vp_b_array[i] = new ArrayList(1);
		rho_b_array[i] = new ArrayList(1);
		if (_SI.MF.ENABLE_HARD_ROCK_MODEL)
		  {
			if (_m.needHardRockTerms) // FIXME Save mem if possible?
			  {
				NG_hard_array[i]  = new ArrayList(1);
				rho_hard_array[i] = new ArrayList(1);
				vp_hard_array[i]  = new ArrayList(1);
				vs_hard_array[i]  = new ArrayList(1);
			  }
		  }
		if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			float_fraction_array[i] = new ArrayList(1);
		  }	
	  }
  }
  /**
	 Similar allocation handling for Synchronised fluid properties. For (zero_offset) layer i.
	 Call using p = 'v','r' or 'S' to allocate vp, rho, S accordingly.
  */
  private void handleFluidSynchArrayAllocs(short fluid_type,int i,char p)
  {
	switch (fluid_type)
	  {	
	  case FluidStates.FluidTypes.brine:
		break;

	  case FluidStates.FluidTypes.oil:
		if ((p=='v') && vp_oil_array[i]==null) vp_oil_array[i] = new ArrayList(1);
		if ((p=='r') && rho_oil_array[i]==null)  rho_oil_array[i] = new ArrayList(1);
		if ((p=='S') && S_oil_array[i]==null)  S_oil_array[i] = new ArrayList(1);
		break;
	  case FluidStates.FluidTypes.lsg:
		if ((p=='v') && vp_lsg_array[i]==null) vp_lsg_array[i] = new ArrayList(1);
		if ((p=='r') && rho_lsg_array[i]==null)  rho_lsg_array[i] = new ArrayList(1);
		if ((p=='S') && S_lsg_array[i]==null)  S_lsg_array[i] = new ArrayList(1);
		break;
	  case FluidStates.FluidTypes.gas:
		if ((p=='v') && vp_gas_array[i]==null) vp_gas_array[i] = new ArrayList(1);
		if ((p=='r') && rho_gas_array[i]==null)  rho_gas_array[i] = new ArrayList(1);
		if ((p=='S') && S_gas_array[i]==null)  S_gas_array[i] = new ArrayList(1);
		break;
	  }
  }

  /**
	 Function to set up prior Gaussian matrix for each fluid combination
  */
  void setUpPriorMatrices()
  {
	int i,j,indx;
	ModelPriorDescriptions.MetaInfo[] Marr;
	ModelPriorDescriptions.MetaInfo MI;
	LayerProperties LPI=null;
	short fluid_type;

	RockandFluidData.RegressionSet CS_m=null, CS_s=null, CS_h=null;

	LayerProperties.Property S = null;
	LayerProperties.Property rho_h = null;
	LayerProperties.Property vp_h = null;

	boolean is_time_synched=false,is_synched=false;
	boolean depth_is_free=true;

	double A=0,B=0,C=0,D=0,sigma=0,vp_bar=0,vs_bar=0,phirho_bar=0,float_fraction_bar=0.0;
	HashMap RockTrends;
	
	short N_possible_fluids, i_h;
	FluidUnitsAndContacts _FUC;
	FluidUnitsAndContacts.FluidUnit fu;
	int fluid_synch_indx = 0;
	boolean fluid_synched_image = false;

	_FUC = _m.fluidUnitsAndContacts;

	RockTrends = _m.RFD.RockTrends;


	Cinv_pr = new DoubleMatrix2D[_m.FS.possibleFluidStates.size()];
	Cinv_prior_only = new DoubleMatrix2D[_m.FS.possibleFluidStates.size()];
	mbar_prior_only = new DoubleMatrix1D[_m.FS.possibleFluidStates.size()];
	logDetC_prior_only = new double[_m.FS.possibleFluidStates.size()];
	mbar_prior_plus_isopachs = new DoubleMatrix1D[_m.FS.possibleFluidStates.size()];		

	// working array for isopach errors
	isopach_errors = new double[_m.MPD.n_isopachs];
	// working array for self-consistent depths
	consistent_depth_errors = new double[_m.n_layers-1];

	// and for layer time-thickness standardised uncertainties
	time_uncertainty = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(_m.FS.possibleFluidStates.size(), _m.n_active_layers);
	// and raw layer time uncertainties
	t_layer_uncertainty = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(_m.FS.possibleFluidStates.size(), _m.n_active_layers+1);

	// loop over fluid combinations
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		map_to = map_to_index[fs];
		map_from = map_from_index[fs];

		int dim = map_from.size();

		Cinv_pr[fs] = new au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D(dim,dim);
		mbar_prior_only[fs] = new cern.colt.matrix.impl.DenseDoubleMatrix1D(dim);
		Cinv = Cinv_pr[fs];
		mu = mbar_prior_only[fs];

		indx=0;
		for (i=0;i<_m.n_layers;i++)
		  {
			LPI = _m.LP[i];
			Marr = _m.MPD.MetaInfo;
			MI = _m.MPD.MetaInfo[i];

			if (_m.MPD.synchMap!=null && (_m.MPD.synchMap.containsKey(i)))
			  {
				is_synched = true;
			  }
			else
			  {
				is_synched = false;
			  }

			is_time_synched = (i==0)? false : !Marr[i-1].isActive;
			// Get indices for properties etc

			// times if variable (may be synched)
			// if (MI.sigma_t.value > 0)
			if (!is_time_synched)
			  {
				addPriorMGDiagEntry(LPI.t, 1.0, MI.sigma_t.value);
				addPriorMGMuEntry(LPI.t,MI.t.value);
			  }

			// Rock trend curves we need for this layer
			CS_m = (RockandFluidData.RegressionSet)RockTrends.get(MI.mud_rock_member);
			if (MI.is_permeable)
			  {
				CS_s = (RockandFluidData.RegressionSet)RockTrends.get(MI.sand_rock_member);
			  }
			if (_SI.MF.ENABLE_HARD_ROCK_MODEL)
			  {
				if (_m.needHardRockTerms && MI.has_hard_rock_mixture)
				  {
					if (MI.hard_shale_rock_member != null)
					  {
						// fixed rock type
						CS_h = (RockandFluidData.RegressionSet)RockTrends.get(MI.hard_shale_rock_member);
					  }
					else // presume  fluid_dep_hard_shale_rock_members!=null by construction
					  {
						// fluid-dependent rock type
						fluid_type = ((FluidStates.FluidState)_m.FS.possibleFluidStates.get(fs)).types[i];
						String s = MI.fluid_dep_hard_shale_rock_members[fluid_type];
						if (s==null)
						  {
							throw new RuntimeException("Required hard rock type, layer i="+(i+1)+", "+MI.layer_name+
													   " doesn't exist for fluid state fs="+fs+" layer_fluid_type="+fluid_type);
						  }
						CS_h = (RockandFluidData.RegressionSet)RockTrends.get(s);
					  }
				  }
				else CS_h = null;				
			  }

			// Model depth. If distinct to trend depth, add a separate model-depth term.
			if (MI.depth != MI.depth_rock_curves)
			  {
				// May be synchronised 
				// with the depth from another layer
				depth_is_free = (MI.depth.master_layer_number == i);
				if (MI.sigma_depth.value > 0)
				  {
					if (depth_is_free)
					  {
						addPriorMGDiagEntry(LPI.d, 1.0, MI.sigma_depth.value);
						addPriorMGMuEntry(LPI.d, MI.depth.value);
					  }
				  }
			  }

			// Trend depths, and coupled terms. 
			// If trend-depth uncertain, add this. May be synchronised with the depth from another layer
			// Remember here, MI.depth_rock_curves references MI.depth if there is no distinct trend depth entry
			depth_is_free = (MI.depth_rock_curves.master_layer_number == i);
			if (MI.sigma_depth_rock_curves.value > 0)
			  {
				if (depth_is_free)
				  {
					addPriorMGDiagEntry(LPI.d_curves, 1.0, MI.sigma_depth_rock_curves.value);
					addPriorMGMuEntry(LPI.d_curves, MI.depth_rock_curves.value);
				  }

				// Possible terms from sand
				if (MI.is_permeable && (CS_s.Vp_versus_depth.error > 0))
				  {
					if (!is_synched) 
					  {
						// extra diag terms from sand
						A = CS_s.Vp_versus_depth.intercept;
						B = CS_s.Vp_versus_depth.depth_coeff;
						sigma = CS_s.Vp_versus_depth.error;
						addPriorMGDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves, B*B,sigma);
						// extra off-diag terms from sand
						addPriorMGSystemOffDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves, LPI.vp_s,-B, sigma);
					  }
				  }

				// Possible terms from impermeable rocks
				// (1) shale
				sigma = CS_m.Vp_versus_depth.error;
				if (sigma>0)
				  {
					// extra diag terms from mud
					A = CS_m.Vp_versus_depth.intercept;
					B = CS_m.Vp_versus_depth.depth_coeff;
					if (!is_synched) 
					  {
						addPriorMGDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves, B*B, sigma);
						// extra off-diag terms from mud
						addPriorMGSystemOffDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves,
													 LPI.vp_m,-B, sigma);
					  }
				  }
				// (2) hard-rock impermeable rock contributions: Here, just coupling terms with depth
				if (_m.needHardRockTerms && MI.has_hard_rock_mixture)
				  {
					sigma = CS_h.Vp_versus_depth.error;
					if (sigma>0)
					  {
						// extra diag terms from hard-rock
						A = CS_h.Vp_versus_depth.intercept;
						B = CS_h.Vp_versus_depth.depth_coeff;
						if (!is_synched) 
						  {
							addPriorMGDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves, B*B, sigma);
							// extra off-diag terms from hard-rock
							addPriorMGSystemOffDiagEntry(_m.LP[MI.depth_rock_curves.master_layer_number].d_curves,
														 LPI.vp_hard,-B, sigma);
						  }
					  }
				  }			
			  }

			// also LFIV. May be synchronised
			// if (!is_synched && MI.sigma_LFIV.value > 0)
			// is not synchronised if the master_layer_number is the same as this layer
			boolean LFIV_is_free = (MI.LFIV.master_layer_number == i);
			if (MI.sigma_LFIV.value > 0)
			  {
				if (LFIV_is_free) 
				  {
					addPriorMGDiagEntry(LPI.LFIV, 1.0, MI.sigma_LFIV.value);
					addPriorMGMuEntry(LPI.LFIV, MI.LFIV.value);
				  }

				// sand
				if (MI.is_permeable && (CS_s.Vp_versus_depth.error > 0))
				  {
					// extra diag terms from sand
					A = CS_s.Vp_versus_depth.intercept;
					B = CS_s.Vp_versus_depth.depth_coeff;
					C = CS_s.Vp_versus_depth.LFIV_coeff;
					sigma = CS_s.Vp_versus_depth.error;
					if (!is_synched) 
					  {
						addPriorMGDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV, C*C, sigma);
						// extra off-diag terms from sand
						// LFIV x vp_s
						addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV,LPI.vp_s,-C, sigma);
						if (MI.sigma_depth_rock_curves.value > 0)
						  {
							// LFIV x depth
							addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV,
														 _m.LP[MI.depth_rock_curves.master_layer_number].d_curves,B*C, sigma);
						  }
					  }
				  }
				// mud
				sigma = CS_m.Vp_versus_depth.error;
				if (sigma>0)
				  {
					// extra diag terms from mud
					A = CS_m.Vp_versus_depth.intercept;
					B = CS_m.Vp_versus_depth.depth_coeff;
					C = CS_m.Vp_versus_depth.LFIV_coeff;
					if (!is_synched) 
					  {
						addPriorMGDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV, C*C, sigma);
						// LFIV x vp_m. 
						// Off-diag entries may still apply if at least one of the pair of properties referenced is real
						// & LPI.vp_m must be real unless the whole layer is synched
						addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV, LPI.vp_m,-C, sigma);
						if (MI.sigma_depth_rock_curves.value > 0)
						  {
							// LFIV x depth
							addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV,
														 _m.LP[MI.depth_rock_curves.master_layer_number].d_curves, B*C, sigma);
						  }
					  }
				  }
				// hard-rock impermeable rock contributions: Here, just coupling terms with LFIV
				if (_m.needHardRockTerms && MI.has_hard_rock_mixture)
				  {
					sigma = CS_h.Vp_versus_depth.error;
					if (sigma>0)
					  {
						// extra diag terms from hard-rock
						A = CS_h.Vp_versus_depth.intercept;
						B = CS_h.Vp_versus_depth.depth_coeff;
						C = CS_h.Vp_versus_depth.LFIV_coeff;
						if (!is_synched) 
						  {
							addPriorMGDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV, C*C, sigma);
							// LFIV x vp_m. 
							// Off-diag entries may still apply if at least one of the pair of properties referenced is real
							// & LPI.vp_hard must be real unless the whole layer is synched
							addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV, LPI.vp_hard,-C, sigma);
							if (MI.sigma_depth_rock_curves.value > 0)
							  {
								// LFIV x depth
								addPriorMGSystemOffDiagEntry(_m.LP[MI.LFIV.master_layer_number].LFIV,
															 _m.LP[MI.depth_rock_curves.master_layer_number].d_curves, B*C, sigma);
							  }
						  }
					  }
				  }		
			  }


			// other standard *** impermeable rock contributions

			sigma = CS_m.Vp_versus_depth.error;
			A = CS_m.Vp_versus_depth.intercept;
			B = CS_m.Vp_versus_depth.depth_coeff;
			C =  CS_m.Vp_versus_depth.LFIV_coeff;
			vp_bar = A + B * Marr[MI.depth_rock_curves.master_layer_number].depth_rock_curves.value
			  + C * Marr[MI.LFIV.master_layer_number].LFIV.value;

			if (!is_synched && sigma > 0)  // Vp entry has 3 terms
			  {
				addPriorMGDiagEntry(LPI.vp_m, 1.0, sigma);
				addPriorMGMuEntry(LPI.vp_m, vp_bar);

				sigma = CS_m.Vs_versus_Vp.error;
				if (!is_synched && sigma>0)
				  {
					B = CS_m.Vs_versus_Vp.slope;
					addPriorMGDiagEntry(LPI.vp_m, B*B, sigma);
					addPriorMGSystemOffDiagEntry(LPI.vp_m, LPI.vs_m, -B, sigma);
				  }

				sigma = CS_m.Rho_or_Phi_versus_Vp.error;
				if (!is_synched && sigma>0)
				  {
					B = CS_m.Rho_or_Phi_versus_Vp.slope;
					addPriorMGDiagEntry(LPI.vp_m, B*B, sigma);
					addPriorMGSystemOffDiagEntry(LPI.vp_m, LPI.rho_m, -B, sigma);
				  }
			  }

			A = CS_m.Vs_versus_Vp.intercept;
			B = CS_m.Vs_versus_Vp.slope;
			sigma = CS_m.Vs_versus_Vp.error;
			vs_bar = A + B*vp_bar;
			if (!is_synched && sigma > 0)  // Vs entry has only 1 term
			  {
				addPriorMGDiagEntry(LPI.vs_m, 1.0, sigma);
				addPriorMGMuEntry(LPI.vs_m, vs_bar);
			  }

			sigma = CS_m.Rho_or_Phi_versus_Vp.error;
			A = CS_m.Rho_or_Phi_versus_Vp.intercept;
			B = CS_m.Rho_or_Phi_versus_Vp.slope;
			phirho_bar = A + B*vp_bar;
			if (!is_synched && sigma > 0)  // rho_m entry has 1 term
			  {
				sigma = CS_m.Rho_or_Phi_versus_Vp.error;
				A = CS_m.Rho_or_Phi_versus_Vp.intercept;
				B = CS_m.Rho_or_Phi_versus_Vp.slope;
				addPriorMGDiagEntry(LPI.rho_m, 1.0, sigma);
				addPriorMGMuEntry(LPI.rho_m, phirho_bar);
			  }
			  
			// Possible hack to make things go if rock props are desired to be fixed. Needed even if layer synched
			if (CS_m.Vp_versus_depth.error == 0)
			  {
				if ((MI.sigma_depth_rock_curves.value == 0) && (MI.sigma_LFIV.value == 0))
				  {
					// Can handle this case OK
					LPI.vp_m.setDistr(vp_bar, 0);
					if (CS_m.Vs_versus_Vp.error == 0) LPI.vs_m.setDistr(vs_bar, 0); // Can handle this case OK					
					if (CS_m.Rho_or_Phi_versus_Vp.error == 0) LPI.rho_m.setDistr(phirho_bar, 0); // Can handle this case OK					
				  }
				else
				  {
					throw new RuntimeException("Check regressions for rock "+MI.mud_rock_member+
											   ": can handle case sigma_vp=0 only if sigma_d(maybe rock_curves) and sigma_LFIV=0 also.");
				  }
			  }
			else if (CS_m.Vs_versus_Vp.error == 0)
			  {
				throw new RuntimeException("Check regressions for rock "+MI.mud_rock_member+
										   ": can handle case sigma_vs=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
			  }
			else if (CS_m.Rho_or_Phi_versus_Vp.error == 0)
			  {
				throw new RuntimeException("Check regressions for rock "+MI.mud_rock_member+
										   ": can handle case sigma_rho_m=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
			  }	
		  


			// hard-rock impermeable rock contributions: FIXME! Below is OK. Possible coupling terms with depth, LFIV are above.
			if (_m.needHardRockTerms && MI.has_hard_rock_mixture)
			  {
				// Net to Gross.
				if (!is_synched && MI.sigma_NG_hard.value > 0) // NG_hard has 1 term
				  {
					addPriorMGDiagEntry(LPI.NG_hard, 1.0, MI.sigma_NG_hard.value);
					addPriorMGMuEntry(LPI.NG_hard, MI.NG_hard.value);					
				  }

				// standard terms corresponding to *** above
				sigma = CS_h.Vp_versus_depth.error;
				A = CS_h.Vp_versus_depth.intercept;
				B = CS_h.Vp_versus_depth.depth_coeff;
				C =  CS_h.Vp_versus_depth.LFIV_coeff;
				vp_bar = A + B * Marr[MI.depth_rock_curves.master_layer_number].depth_rock_curves.value
				  + C * Marr[MI.LFIV.master_layer_number].LFIV.value;
				if (!is_synched && sigma > 0)  // Vp entry has 3 terms
				  {
					addPriorMGDiagEntry(LPI.vp_hard, 1.0, sigma);
					addPriorMGMuEntry(LPI.vp_hard, vp_bar);

					sigma = CS_h.Vs_versus_Vp.error;
					if (!is_synched && sigma>0)
					  {
						B = CS_h.Vs_versus_Vp.slope;
						addPriorMGDiagEntry(LPI.vp_hard, B*B, sigma);
						addPriorMGSystemOffDiagEntry(LPI.vp_hard, LPI.vs_hard, -B, sigma);
					  }

					sigma = CS_h.Rho_or_Phi_versus_Vp.error;
					if (!is_synched && sigma>0)
					  {
						B = CS_h.Rho_or_Phi_versus_Vp.slope;
						addPriorMGDiagEntry(LPI.vp_hard, B*B, sigma);
						addPriorMGSystemOffDiagEntry(LPI.vp_hard, LPI.rho_hard, -B, sigma);
					  }
				  }

				A = CS_h.Vs_versus_Vp.intercept;
				B = CS_h.Vs_versus_Vp.slope;
				sigma = CS_h.Vs_versus_Vp.error;
				vs_bar = A + B*vp_bar;
				if (!is_synched && sigma > 0)  // Vs entry has only 1 term
				  {
					addPriorMGDiagEntry(LPI.vs_hard, 1.0, sigma);
					addPriorMGMuEntry(LPI.vs_hard, vs_bar);
				  }

				sigma = CS_h.Rho_or_Phi_versus_Vp.error;
				A = CS_h.Rho_or_Phi_versus_Vp.intercept;
				B = CS_h.Rho_or_Phi_versus_Vp.slope;
				phirho_bar = A + B*vp_bar;
				if (!is_synched && sigma > 0)  // rho_hard entry has 1 term
				  {
					addPriorMGDiagEntry(LPI.rho_hard, 1.0, sigma);
					addPriorMGMuEntry(LPI.rho_hard, phirho_bar);
				  }

				// Possible hack to make things go if rock props are desired to be fixed. Needed also if synched
				if (CS_h.Vp_versus_depth.error == 0)
				  {
					if ((MI.sigma_depth_rock_curves.value == 0) && (MI.sigma_LFIV.value == 0))
					  {
						// Can handle this case OK
						LPI.vp_hard.setDistr(vp_bar, 0);
						if (CS_h.Vs_versus_Vp.error == 0) LPI.vs_hard.setDistr(vs_bar, 0); // Can handle this case OK					
						if (CS_h.Rho_or_Phi_versus_Vp.error == 0) LPI.rho_hard.setDistr(phirho_bar, 0); // Can handle this case OK					
					  }
					else
					  {
						throw new RuntimeException("Check regressions for rock "+MI.hard_shale_rock_member+
												   ": can handle case sigma_vp=0 only if sigma_d(maybe rock_curves) and sigma_LFIV=0 also.");
					  }
				  }
				else if (CS_h.Vs_versus_Vp.error == 0)
				  {
					throw new RuntimeException("Check regressions for rock "+MI.hard_shale_rock_member+
											   ": can handle case sigma_vs=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
				  }
				else if (CS_h.Rho_or_Phi_versus_Vp.error == 0)
				  {
					throw new RuntimeException("Check regressions for rock "+MI.hard_shale_rock_member+
											   ": can handle case sigma_rho_m=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
				  }	
			  } // _m.needHardRockTerms && MI.has_hard_rock_mixture
			  

			// contributions from permeable rock and fluids........

			if (MI.is_permeable)
			  {
				// Sand & maybe fluid properties will enter the state vector

				// Floating grain effects
				float_fraction_bar = 0.0;
				if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
				  {
					// Full capture model with cross terms in regressions
					
					// is not synchronised if the master_layer_number is the same as this layer
					boolean float_fraction_is_free = (MI.float_fraction.master_layer_number == i);
					float_fraction_bar = MI.float_fraction.value;
					if (MI.sigma_float_fraction.value > 0)
					  {
						if (float_fraction_is_free) 
						  {
							addPriorMGDiagEntry(LPI.float_fraction, 1.0, MI.sigma_float_fraction.value);
							addPriorMGMuEntry(LPI.float_fraction, float_fraction_bar);
						  }
						if (CS_s.Vp_versus_depth.error > 0)
						  {
							// extra diag terms from sand
							A = CS_s.Vp_versus_depth.intercept;
							B = CS_s.Vp_versus_depth.depth_coeff;
							D = CS_s.Vp_versus_depth.float_fraction_coeff;
							sigma = CS_s.Vp_versus_depth.error;
							if (!is_synched) 
							  {
								addPriorMGDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction, D*D, sigma);
								// extra off-diag terms from vp trend
								// float_fraction x vp_s
								addPriorMGSystemOffDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction,LPI.vp_s,-D, sigma);
								if (MI.sigma_depth_rock_curves.value > 0)
								  {
									// float_fraction x depth
									addPriorMGSystemOffDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction,
																 _m.LP[MI.depth_rock_curves.master_layer_number].d_curves,B*D, sigma);
								  }
								if (MI.sigma_LFIV.value > 0)
								  {
									// float_fraction x LFIV
									addPriorMGSystemOffDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction,
																 _m.LP[MI.LFIV.master_layer_number].LFIV,C*D, sigma);
								  }
							  }
						  }
						if (CS_s.Rho_or_Phi_versus_Vp.error > 0)
						  {
							// porosity relation (since Permeable)
							A = CS_s.Rho_or_Phi_versus_Vp.intercept;
							B = CS_s.Rho_or_Phi_versus_Vp.slope;
							D = CS_s.Rho_or_Phi_versus_Vp.float_fraction_coeff;
							sigma = CS_s.Rho_or_Phi_versus_Vp.error;
							if (!is_synched) 
							  {
								addPriorMGDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction, D*D, sigma);
								// extra off-diag terms from vp trend
								// float_fraction x phi
								addPriorMGSystemOffDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction,LPI.phi_s,-D, sigma);
								// float_fraction x vp (vp guaranteed to be stochastic)
								addPriorMGSystemOffDiagEntry(_m.LP[MI.float_fraction.master_layer_number].float_fraction,
															 _m.LP[MI.depth_rock_curves.master_layer_number].vp_s,B*D, sigma);
							  }
						  }
					  }
				  }

				// Terms from Permeable rock regressions (vp_s, vs_s, phi_s etc )
				sigma = CS_s.Vp_versus_depth.error;
				A = CS_s.Vp_versus_depth.intercept;
				B =  CS_s.Vp_versus_depth.depth_coeff;
				C =  CS_s.Vp_versus_depth.LFIV_coeff;
				vp_bar = A + B * Marr[MI.depth_rock_curves.master_layer_number].depth_rock_curves.value 
				  + C * Marr[MI.LFIV.master_layer_number].LFIV.value; 
					
				if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
				  {
					D = CS_s.Vp_versus_depth.float_fraction_coeff;
					vp_bar += D * Marr[MI.float_fraction.master_layer_number].float_fraction.value;
				  }
				if (!is_synched && sigma > 0)  // Vp entry has 3 terms
				  {
					addPriorMGDiagEntry(LPI.vp_s, 1.0, sigma);
					addPriorMGMuEntry(LPI.vp_s, vp_bar);

					sigma = CS_s.Vs_versus_Vp.error;
					if (!is_synched && sigma>0)
					  {
						B = CS_s.Vs_versus_Vp.slope;
						addPriorMGDiagEntry(LPI.vp_s, B*B, sigma);
						addPriorMGSystemOffDiagEntry(LPI.vp_s, LPI.vs_s, -B, sigma);
					  }

					sigma = CS_s.Rho_or_Phi_versus_Vp.error;
					if (!is_synched && sigma>0)
					  {
						B = CS_s.Rho_or_Phi_versus_Vp.slope;
						addPriorMGDiagEntry(LPI.vp_s,  B*B, sigma);
						addPriorMGSystemOffDiagEntry(LPI.vp_s, LPI.phi_s, -B, sigma);
					  }
				  }
				  
				sigma = CS_s.Vs_versus_Vp.error;
				A = CS_s.Vs_versus_Vp.intercept;
				B = CS_s.Vs_versus_Vp.slope;
				vs_bar = A+B*vp_bar;
				if (!is_synched && sigma > 0)  // Vs entry has only 1 term
				  {
					addPriorMGDiagEntry(LPI.vs_s, 1.0, sigma);
					addPriorMGMuEntry(LPI.vs_s, vs_bar);
				  }

				sigma = CS_s.Rho_or_Phi_versus_Vp.error;
				A = CS_s.Rho_or_Phi_versus_Vp.intercept;
				B = CS_s.Rho_or_Phi_versus_Vp.slope;
				phirho_bar = A + B*vp_bar;
				if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL) 
				  {
					C = CS_s.Rho_or_Phi_versus_Vp.float_fraction_coeff;
					phirho_bar += C* float_fraction_bar;
				  }
				if (!is_synched && sigma > 0)  // phi_s entry has 1 term
				  {
					addPriorMGDiagEntry(LPI.phi_s, 1.0, sigma);
					addPriorMGMuEntry(LPI.phi_s, phirho_bar);
				  }
				// Possible hack to make things go if rock props are desired to be fixed. Needed even if synched
				if (CS_s.Vp_versus_depth.error == 0)
				  {
					if ((MI.sigma_depth_rock_curves.value == 0) && (MI.sigma_LFIV.value == 0))
					  {
						// Can handle this case OK
						LPI.vp_s.setDistr(vp_bar, 0);
						if (CS_s.Vs_versus_Vp.error == 0) LPI.vs_s.setDistr(vs_bar, 0); // Can handle this case OK						
						if (CS_s.Rho_or_Phi_versus_Vp.error == 0) LPI.phi_s.setDistr(phirho_bar, 0); // Can handle this case OK						
					  }
					else
					  {
						throw new RuntimeException("Check regressions for rock "+MI.sand_rock_member+
												   ": can handle case sigma_vp=0 only if sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
					  }
				  }
				else if (CS_s.Vs_versus_Vp.error == 0)
				  {
					throw new RuntimeException("Check regressions for rock "+MI.sand_rock_member+
											   ": can handle case sigma_vs=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
				  }
				else if (CS_s.Rho_or_Phi_versus_Vp.error == 0)
				  {
					throw new RuntimeException("Check regressions for rock "+MI.sand_rock_member+
											   ": can handle case sigma_phi_s=0 only if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
				  }					
				  

				// Net to Gross.
				if (!is_synched && MI.sigma_NG.value > 0) // NG has 1 term
				  {
					addPriorMGDiagEntry(LPI.NG, 1.0, MI.sigma_NG.value);
					addPriorMGMuEntry(LPI.NG, MI.NG.value);
				  }



				// fluid terms. Independent of synchronisation
				// Here either 1 or 3 possible hydrocarbon cases, 
				// for enumeration-mode, or contacts-mode respectively
				// Here 3 = FluidStates.FluidTypes.Nfluids-1
				N_possible_fluids = _SI.doSimpleFluidUnitsAndContacts ? (short)3: (short)1;
				for (i_h=1; i_h<=N_possible_fluids; i_h++) // either {1}, or {1,2,3}
				  {
					if (_SI.doSimpleFluidUnitsAndContacts)
					  {
						fluid_type = i_h;
						if (!_FUC.fluidRequired[i][i_h]) continue; // skip rest of loop

						fluid_synched_image = false;
						if (_FUC.fluidsSynchedForFUowningLayer(i))
						  {
							fu = _FUC.getFUofThisLayer(i);
							fluid_synch_indx = fu.masterFluidLayer[i_h];
							fluid_synched_image = (fluid_synch_indx!=i); // true on a fluid "image" layer
						  }
					  }
					else
					  {
						// fluid required by current fluid state
						fluid_type = ((FluidStates.FluidState)_m.FS.possibleFluidStates.get(fs)).types[i];
					  }

					switch (fluid_type)
					  {
					  case FluidStates.FluidTypes.brine:
						S = null;
						vp_h = null;
						rho_h = null;
						break;

					  case FluidStates.FluidTypes.oil:
						S = LPI.S_oil;
						vp_h = LPI.vp_oil;
						rho_h = LPI.rho_oil;
						break;

					  case FluidStates.FluidTypes.lsg:
						S = LPI.S_lsg;
						vp_h = LPI.vp_lsg;
						rho_h = LPI.rho_lsg;
						break;

					  case FluidStates.FluidTypes.gas:
						S = LPI.S_gas;
						vp_h = LPI.vp_gas;
						rho_h = LPI.rho_gas;
						break;
					  }
					if (i_h==1)
					  {
						// add these only once
						if (map_to.get(LPI.vp_b)!=null && !fluid_synched_image)  // the ref fluid vp is uncertain
						  {
							addPriorMGDiagEntry(LPI.vp_b, 1.0, LPI.vp_b.sigma);
							addPriorMGMuEntry(LPI.vp_b, LPI.vp_b.mean);
						  }
						if (map_to.get(LPI.rho_b)!=null && !fluid_synched_image)  // the new brine rho is uncertain
						  {
							addPriorMGDiagEntry(LPI.rho_b, 1.0, LPI.rho_b.sigma);
							addPriorMGMuEntry(LPI.rho_b, LPI.rho_b.mean);
						  }
					  }
					if (map_to.get(S)!=null && !fluid_synched_image)  // there must have been a hydrocarbon of uncertain Sat.
					  {
						// if a hydrocarbon can possibly be there
						// with uncertain saturation, add S
						addPriorMGDiagEntry(S, 1.0, S.sigma);
						addPriorMGMuEntry(S, S.mean);
					  }
					if (map_to.get(vp_h)!=null && !fluid_synched_image)  // there must have been a hydrocarbon of uncertain vp
					  {
						// if a hydrocarbon can possibly be there
						// with uncertain vp add vp
						addPriorMGDiagEntry(vp_h, 1.0, vp_h.sigma);
						addPriorMGMuEntry(vp_h, vp_h.mean);
					  }
					if (map_to.get(rho_h)!=null && !fluid_synched_image)  // there must have been a hydrocarbon of uncertain rho
					  {
						// if a hydrocarbon can possibly be there
						// with uncertain saturation, add rho_h
						addPriorMGDiagEntry(rho_h, 1.0, rho_h.sigma);
						addPriorMGMuEntry(rho_h, rho_h.mean);
					  }
				  }

			  } // if() contributions from permeable rocks and fluids
		  } // loop over layers

		// Finally, add terms for the base time.
		// base times if variable. Never synchronised
		if (_m.MPD.sigma_t_base.value > 0)
		  {
			addPriorMGDiagEntry(LPI.t_base, 1.0, _m.MPD.sigma_t_base.value);
			addPriorMGMuEntry(LPI.t_base, _m.MPD.t_base.value);
		  }

		// Internal stochastic variables
		// (1)  of the Seismic model
		if (_m.SM.sigma_A_avo.value > 0)
		  {
			addPriorMGDiagEntry(_m.SM.A, 1.0, _m.SM.sigma_A_avo.value);
			addPriorMGMuEntry(_m.SM.A, _m.SM.A_avo.value);
		  }
		// Check if B_avo is synched to A
		if (!_m.SM.A_B_avo_synchronised &&  (_m.SM.sigma_B_avo.value > 0))
		  {
			addPriorMGDiagEntry(_m.SM.B, 1.0, _m.SM.sigma_B_avo.value);
			addPriorMGMuEntry(_m.SM.B, _m.SM.B_avo.value);
		  }

		// (2) Simple contacts parameters
		if (_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			for (i=0;i<	_m.fluidUnitsAndContacts.fluidUnits.length; i++)
			  {
				fu = _m.fluidUnitsAndContacts.fluidUnits[i];
				for (j=0; j< fu.contactInfo.length; j++)
				  {
					FluidUnitsAndContacts.ContactInfo CI = 	fu.contactInfo[j];
					if (CI.sigma_depthPD.value > 0)
					  {
						addPriorMGDiagEntry(CI.depth, 1.0,CI.sigma_depthPD.value);
						addPriorMGMuEntry(CI.depth, CI.depthPD.value);
					  }
				  }
			  }
		  }
	  } // loop over fluid states


	if (_SI.doStackedModel) 
	  {
		// working array for vertical stacking
		_m.stackedModelInfo.setUpConstraints();
		modifyPriorMatricesForVerticalStacking();
		vertical_stacking_errors = new double[_m.stackedModelInfo.n_constraints];
	  }
  }
  
  /** 
	 Modifies prior covariances to reflect terms from VerticalStacking
	 Not hard-rock complaint yet....
  */
  void modifyPriorMatricesForVerticalStacking()
  {
	int i,j,i_b,jp,k;
	StackedModelInfo _smi;
	StackedModelInfo.StackedPropertyMappingDetails[] _l_maps;	
	Integer I,J;
	double sigma=0.0;

	_smi = _m.stackedModelInfo;
	_l_maps = _m.stackedModelInfo.l_maps;

	// loop over fluid combinations
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		map_to = map_to_index[fs];
		Cinv = Cinv_pr[fs];
		mu = mbar_prior_only[fs];

		// make copy 
		DoubleMatrix1D mu_tmp = mu.like();
		// form mu_tmp = C_p^{-1}.mu
		Cinv.zMult(mu,mu_tmp);

		// sum over target layer range
		for (i_b=0; i_b<_l_maps.length; i_b++)
		  {
			StackedModelInfo.StackedPropertyMappingDetails map = _l_maps[i_b];
			j = map.base_layer-1; // indx of base layer

			LayerProperties LPI = _m.LP[j];

			for (k=0; k< map.target_layers.length; k++)
			  {
				jp = map.target_layers[k] - 1;  // index of comparison layer

				LayerProperties LPJ  = _m.LP[jp];
			
				// vp_m
				if ((I = (Integer)(map_to.get(LPI.vp_m)))!=null && ((J = (Integer)(map_to.get(LPJ.vp_m)))!=null) 
					&&  _smi.vp_m.doConstraint && _smi.vp_m.isApplicable[j])  // the vp_m is uncertain
				  {
					sigma = _smi.vp_m.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }   
				// vp_s
				if ((I = (Integer)(map_to.get(LPI.vp_s)))!=null && ((J = (Integer)(map_to.get(LPJ.vp_s)))!=null) 
					&&  _smi.vp_s.doConstraint && _smi.vp_s.isApplicable[j])  // the vp_s is uncertain
				  {
					sigma = _smi.vp_s.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// vs_m
				if ((I = (Integer)(map_to.get(LPI.vs_m)))!=null && ((J = (Integer)(map_to.get(LPJ.vs_m)))!=null) 
					&&  _smi.vs_m.doConstraint && _smi.vs_m.isApplicable[j])  // the vs_m is uncertain
				  {
					sigma = _smi.vs_m.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// vs_s
				if ((I = (Integer)(map_to.get(LPI.vs_s)))!=null && ((J = (Integer)(map_to.get(LPJ.vs_s)))!=null) 
					&&  _smi.vs_s.doConstraint && _smi.vs_s.isApplicable[j])  // the vs_s is uncertain
				  {
					sigma = _smi.vs_s.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// rho_m
				if ((I = (Integer)(map_to.get(LPI.rho_m)))!=null && ((J = (Integer)(map_to.get(LPJ.rho_m)))!=null) 
					&&  _smi.rho_m.doConstraint && _smi.rho_m.isApplicable[j])  // the rho_m is uncertain
				  {
					sigma = _smi.rho_m.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// phi_s
				if ((I = (Integer)(map_to.get(LPI.phi_s)))!=null && ((J = (Integer)(map_to.get(LPJ.phi_s)))!=null) 
					&&  _smi.phi_s.doConstraint && _smi.phi_s.isApplicable[j])  // the phi_s is uncertain
				  {
					sigma = _smi.phi_s.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// NG
				if ((I = (Integer)(map_to.get(LPI.NG)))!=null && ((J = (Integer)(map_to.get(LPJ.NG)))!=null) 
					&&  _smi.NG.doConstraint && _smi.NG.isApplicable[j])  // the NG is uncertain
				  {
					sigma = _smi.NG.sigma;
					addPriorMGSystemSymmetricDifferenceEntry(I,J,sigma);
					// logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
				// thickness
				// this is nonlinear, so cannot be incorporated directly into prior
				if (false)
				  {
					if (_smi.thickness.doConstraint && _smi.thickness.isApplicable[j])  // the thickness is uncertain
					  {
						i++;
						sigma = _smi.thickness.sigma;
						// logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
				  }
			  }
		  }

		// now update mean
		// At this point Cp^{-1} has been updated to  Cp^{-1} + C_VS^{-1}
		// Need to update mu to mu <- (Cp^{-1} + C_VS^{-1})^{-1} .  Cp^{-1} . mu
		// Do this by the linear problem for x: (Cp^{-1} + C_VS^{-1}) .  x = Cp^{-1} . mu = mu_tmp

		// copy Cinv
		DoubleMatrix2D Cinv_cp = Cinv.like().assign(Cinv);
		LUDecompositionQuick LUD = new LUDecompositionQuick();
		// decompose destructively
		LUD.decompose(Cinv_cp);
		if (LUD.isNonsingular())
		  {
			//Util.debugPrintMatrixToFile("mu.txt",mu);
			LUD.solve(mu_tmp);
			//Util.debugPrintMatrixToFile("mu_p.txt",mu_tmp);
			mbar_prior_only[fs] = mu_tmp;
			//System.exit(0);
		  }
		else
		  {
			_SI.EH.die_err("Warning!: Vertical stacking prior-mean update impossible. Specification is impossible!");
		  }
		Cinv_cp = null;
		LUD = null;
	  }

	// restore maps to safe state.
	map_to = map_to_index[0];
	map_from = map_from_index[0];
  }




  private void addPriorMGDiagEntry(LayerProperties.Property prop,double C_inv_numerator,double sigma)
  {
	Integer I;
	int indx;
	double var,entry;

	// get the correct index for this property
	I = (Integer)map_to.get(prop);
	if (I==null) return;
	indx = I.intValue();

	// add the diag entry in C^{-1}
	var = MathUtil.sqr(sigma);
	if (var==0)
	  {
		System.out.println("PROBLEM1!");
	  }
	entry = Cinv.getQuick(indx,indx);
	entry += C_inv_numerator/var;
	Cinv.setQuick(indx,indx,entry);

  }

  private void modifyPriorMGMuEntry(LayerProperties.Property prop, double mu_entry)
  {
	Integer I;
	int indx;
	double var,entry;

	// get the correct index for this property
	I = (Integer)map_to.get(prop);
	if (I==null) return;
	indx = I.intValue();
	mu.setQuick(indx,mu_entry);
  }

 private void addPriorMGMuEntry(LayerProperties.Property prop, double mu_entry)
  {
	Integer I;
	int indx;
	double var,entry;

	// get the correct index for this property
	I = (Integer)map_to.get(prop);
	if (I==null) return;
	indx = I.intValue();

	// add the entry in the RHS mu
	entry = mu.getQuick(indx);
	// add mu entry by inspection
	entry += mu_entry;
	mu.setQuick(indx,entry);
  }

  private void addPriorMGSystemOffDiagEntry(LayerProperties.Property prop1,LayerProperties.Property prop2,
											double numerator, double sigma)
  {
	Integer I;
	int indx1,indx2;
	double var,entry;

	// get the correct indices for this property
	I = (Integer)map_to.get(prop1);
	if (I==null) return;
	indx1 = I.intValue();

	I = (Integer)map_to.get(prop2);
	if (I==null) return;
	indx2 = I.intValue();

	// add the off-diag entry in C^{-1}
	var = MathUtil.sqr(sigma);
	if (var==0) _SI.EH.println("PROBLEM2 in addPriorMGSystemOffDiagEntry()!");

	// add to existing entry
	entry = Cinv.getQuick(indx1,indx2);
	entry += numerator/var;
	Cinv.setQuick(indx1,indx2,entry);

	// and the symmetric entry....
	Cinv.setQuick(indx2,indx1,entry);
  }


  /**
	 For difference terms (vertical stacking) in prior.
	 Adds a {{1,-1},{-1,1}} / sigma^2 block entry to elements I,J of Cinv
  */
  private void addPriorMGSystemSymmetricDifferenceEntry(Integer I, Integer J,double sigma)
  {
	int i,j;
	double dx;

	i = I.intValue();
	j = J.intValue();
	
	dx = 1.0/(sigma*sigma);
					
	// add in block entries
	Cinv.setQuick(i,i,dx + Cinv.getQuick(i,i));
	Cinv.setQuick(j,j,dx + Cinv.getQuick(j,j));
	Cinv.setQuick(i,j,-dx + Cinv.getQuick(i,j));
	Cinv.setQuick(j,i,-dx + Cinv.getQuick(j,i));
  }

  /**
	 This method does a quick check to confirm that Gassman
	 substutions derived from the prior for all fluid modes of
	 interest are possible.  If running with checkPriorSanity set in
	 the SeismicInverter, only error spiel is produced (this will then
	 continue to other traces).  If running with
	 checkPriorSanity=false, an exception will be thrown if the basic
	 prior model is unacceptable, as this precludes sensible
	 inversion.
  */
  void testGassmanSanity()
  {
	int i;
	DoubleMatrix1D mu;
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		setMapsFromFluidState(fs);
		mu = mbar_prior_only[fs];
		for (i=0;i<mu.size();i++) 
		  {
			setValuesFromIndex(i,mu.getQuick(i));
		  }
		_m.computeEffectivePropertiesForAllLayers(false,true,false,false,false);
	  }
  }

  /**
	 This method does a quick check to confirm that Gassman
	 substutions at the state vector m_test for the current fluid state
	 interest are possible. Returns true or false accordingly
  */
  boolean isGassmanSubstitutionPossible(DoubleMatrix1D m_test)
  {
	int i;
	boolean possible;
	
	for (i=0;i<m_test.size();i++) 
	  {
		setValuesFromIndex(i,m_test.getQuick(i));
	  }

	// order the times
	orderTimes();
	possible = _m.computeEffectivePropertiesForAllLayers(false,false,true,false,false);

	return possible;
  }
  
  // method that assumes x is a legal point
  private void setVariablesAtThisLegalGassmanPoint(DoubleMatrix1D x)
  {
	isGassmanSubstitutionPossible(x);
  }

  /**
	 Finds optimal point in parameter space along the line joing xold to xnew that
	 yields a valid Gassman substitution. xnew is tried first. In case of failure, the
	 nearest valid point to xnew along the line is found by bisection
  */
  DoubleMatrix1D getBestLegalGassmanPointAlongCurrentLine(DoubleMatrix1D xold, DoubleMatrix1D xnew)
  {
	int j,k;
	DoubleMatrix1D xmid;
	DoubleMatrix1D xright=xnew;
	double norm=1.0,fM=1,fL=0,fR=1,u;

	boolean debug=false;

	// test xnew to see if it's legal
	if (isGassmanSubstitutionPossible(xnew))  
	  {
		return xnew;
	  }
	else
	  {
		// find the best point along the bisection line
		xmid = xold.like();
		if (debug) _SI.EH.println("BIS: %10.8g%n",fM);
		// find using bisection
        for (j=0;j<GASSMAN_BISECTION_MAX_ITS;j++) 
		  {
			// construct interpolater
			fM = (fL+fR)/2.0;
			u = (1-fM);
			for(k=0;k<xold.size();k++)
			  {				
				xmid.setQuick(k, u*xold.getQuick(k) + fM*xright.getQuick(k));
			  }
			if (isGassmanSubstitutionPossible(xmid)) // FC_FIXME?
			  {			
				fL=fM;
			  }
			else
			  {
				fR=fM;
			  } 
			if (debug) _SI.EH.println("BIS: %10.8g ",fM);
			norm/=2.0;
			if (norm < GASSMAN_BISECTION_EPS) 
			  {
				// ensure we are at least GASSMAN_BISECTION_EPS safely away from the
				// failure point...
				fM = Math.max(0,fM - GASSMAN_BISECTION_EPS);
				u = (1-fM);
				for (k=0;k<xold.size();k++)
				  {
					xmid.setQuick(k, u*xold.getQuick(k) + fM*xright.getQuick(k));
				  }
				orderTimes();
				_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
				if (debug) _SI.EH.printf("BIS: %10.8g%nBIS:%n",fM);
				return xmid;
			  }
		  }
	  }
	// never gets here
	return xold;
  }

  /**
	 Apply depth self consistency only if all depths have intrinsic freedom (non-zero std dev's)
  */
  private void determineWhetherToApplyDepthConsistency()
  {
	int i;
	ModelPriorDescriptions.MetaInfo MI;

	// assume true unless we prove it shouldn't be done
	applyDepthConsistency = true;
	for (i=0;i<_m.n_layers;i++)
	  {
		MI = _m.MPD.MetaInfo[i];
		// Depth may be synchronised with the depth from another layer
		boolean depth_is_free = (MI.depth.master_layer_number == i);
		if (MI.sigma_depth.value ==0 && depth_is_free)
		  {
			// it's a free layer but has no depth freedom. Can't work
			applyDepthConsistency = false;
		  }
	  }
	// we may have insisted on neglecting it also
	if (_SI.ignoreDepthConsistency) applyDepthConsistency = false;

	// Choose also a sane value for the depth self-consistent error
	// make it the mean depth of the last layer / 1000
	if (applyDepthConsistency)
	  {
		sigma_consistent_depth = _m.MPD.MetaInfo[_m.n_layers-1].depth.value * 
		  EPSILON_FOR_DEPTH_SELF_CONSISTENCY;
		// protection
		if (sigma_consistent_depth<=0) sigma_consistent_depth = DEFAULT_SIGMA_DEPTH_SELF_CONSISTENCY;
	  }
	if (_verbosity>3) _SI.EH.println("Depth self-consistency (sigma_d = %-6.4g) applied = %s",
									 sigma_consistent_depth, applyDepthConsistency);

  }

  void incorporateLinearisedIsopachAndOrDepthConstraints()
  {
	int i,j,jj,k,n;
	int maxNewtonSteps;
	int nt,n_constraints,nd_start,nvs_start;
	ModelPriorDescriptions.MetaInfo MI;
	au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D X;
	DenseDoubleMatrix1D d;
	double DZ_mean, DZ_new, t_base, t_base_new, val, newval, rhs,tmp;
	Integer I;
	boolean debug=false;
	double[] base_depth_errors;
	double[] base_vertical_stacking_errors=null;
	R_prior_only = new DoubleMatrix2D[_m.FS.possibleFluidStates.size()];
	R_prior_plus_isopachs = new DoubleMatrix2D[_m.FS.possibleFluidStates.size()];
	Cinv_prior_plus_isopachs = new DoubleMatrix2D[_m.FS.possibleFluidStates.size()];
	double[] base_errors = (double [])isopach_errors.clone();



	determineWhetherToApplyDepthConsistency();

	base_depth_errors = (double [])consistent_depth_errors.clone();
	if (_SI.doStackedModel)
	  {
		base_vertical_stacking_errors = (double [])vertical_stacking_errors.clone();
	  }
	
	if (!(_m.MPD.n_isopachs==0 || _SI.ignoreIsopachs) || applyDepthConsistency || _SI.doStackedModel)
	  {
		detC_prior_plus_isopachs = new double[_m.FS.possibleFluidStates.size()];
	  }
		
	time_std_devs = new TreeMap[_m.FS.possibleFluidStates.size()];
	// loop over fluid combinations, obtaining X matrix for isopachs if any
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		time_std_devs[fs] = new TreeMap();

		// (A) If there are no isopachs, or we're ignoring them, do this
		// (later on we must have the Cholesky factors, determinants etc)
		// Also, we can only get away with this is we are ignoring depth self-consistency and vertical stacking
		if ((_m.MPD.n_isopachs==0 || _SI.ignoreIsopachs) && !applyDepthConsistency && !_SI.doStackedModel)
		  {
			// There will be no isopach update, so these three matrices are the same
			Cinv_prior_plus_isopachs[fs] =  Cinv_pr[fs];
			Cinv_prior_only[fs] =  Cinv_pr[fs];
			// as will the prior mean
			mbar_prior_plus_isopachs[fs] = mbar_prior_only[fs];

			// And we need the Cholesky factor
			ModifiedCholeskyRightDecomposition CDR 
			  = new ModifiedCholeskyRightDecomposition(Cinv_pr[fs],true); 

			R_prior_plus_isopachs[fs] =  CDR.getR();

			R_prior_only[fs] = R_prior_plus_isopachs[fs].copy();
			// R_prior_only is now an upper triangular Cholesky factor of Cinv_pr

			CDR.getCholeskyRightInverse(R_prior_plus_isopachs[fs]);  
			// R_prior_plus_isopachs is now an upper triangular Cholesky factor of C_pr = Cinv_pr^-1
			// viz C = {\tilde C} = R_prior_plus_isopachs[fs]^T . R_prior_plus_isopachs[fs]


			// get and save the log of its determinant
			logDetC_prior_only[fs]=0.0;
			for (i=0;i<R_prior_plus_isopachs[fs].rows();i++) 
			  {
				logDetC_prior_only[fs] += Math.log((tmp=R_prior_plus_isopachs[fs].getQuick(i,i)) * tmp);
			  }

			map_to =  map_to_index[fs];
			aquireLayerTimeThicknessStandardisedUncertainties();
			// everything in (B) below is redundant since no isopach-updating need be done
			CDR=null;
			if (debug && fs==0)   // FIXME COMMENT OUT
			  {
				Util.debugPrintMatrixToFile("Cinv.txt",Cinv_pr[fs]);
				Util.debugPrintMatrixToFile("mbar.txt",mbar_prior_only[fs]);
			  }
			continue;  // will force skip of rest of loop!!!!!!!
		  }


		// (B) But if there are isopachs, we must do the standard linearised update
		//     as follows....

		if (_verbosity>3) _SI.EH.println("Isopachs for fluid state "+fs);


		fluid_state = (FluidStates.FluidState) _m.FS.possibleFluidStates.get(fs);
		n = mbar_prior_only[fs].size();
		map_to =  map_to_index[fs];
		map_from =  map_from_index[fs];
		
		// take copy of original prior matrix before it gets corrupted
		Cinv_prior_only[fs] =  Cinv_pr[fs].copy();

		// let the initial expansion point for the isopach likelihood be the prior mean
		// DoubleMatrix1D expansion_point = (DoubleMatrix1D) mbar_prior_only[fs].copy();
		DoubleMatrix1D expansion_point = new DenseDoubleMatrix1D(mbar_prior_only[fs].size());
		expansion_point.assign(mbar_prior_only[fs]);
		// set up X, d of linearised joint constraint (union of isopachs + depth-consistency)
		n_constraints = (_SI.ignoreIsopachs ? 0 : _m.MPD.n_isopachs);

		if (applyDepthConsistency) n_constraints += (_m.n_layers-1);

		if (_SI.doStackedModel) n_constraints += _m.stackedModelInfo.n_constraints;

		// NB logical entry to this piece of code will ensure n_constraints will always be non-zero
		X = new au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D(n_constraints,n);
		d = new DenseDoubleMatrix1D(n_constraints);

		// starting index of depth self-consistency block in linear algebra
		nd_start = (_SI.ignoreIsopachs ? 0 : _m.MPD.n_isopachs);
		// starting index of vertical-stack block in linear algebra
		nvs_start = nd_start + ((applyDepthConsistency)? (_m.n_layers-1) : 0);

		// Working array.
		ArrayList<LayerProperties.Property> prop_arrayl = new ArrayList<LayerProperties.Property>();
		
		// This update can be iterated: 3 times seems to produce ok results.
		// I've noticed single updates are occasionally flaky, especially with depth-self-consistency switched on
		maxNewtonSteps=MAX_NEWTON_STEPS;
		for (int NewtonIts=0; NewtonIts<maxNewtonSteps; NewtonIts++)
		  {
			// set all layer-property values to the expansion point &
			// compute all effective properties
			for (k=0;k<n;k++)
			  {
				val = expansion_point.getQuick(k);
				setValuesFromIndex(k,val);
			  }
			// order the times
			orderTimes();  // FOR_FC
			_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK

			// (1) compute base isopach errors
			computeIsopachNegLogLikelihood(false);

			// info about isopach errors before we do the linear update
			if (!_SI.ignoreIsopachs && _verbosity>3)
			  {
				_SI.EH.print("   Pre-update mean isopach deviates:  ");
				for (k=0;k<_m.MPD.n_isopachs;k++)
				  {
					_SI.EH.printf("%10.4g ",isopach_errors[k]);
				  }
				_SI.EH.print("%n");
			  }

			// Take copy of base isopach-errors errors for finite differences
			System.arraycopy(isopach_errors,0,base_errors,0,isopach_errors.length);

			// (2) Base self-consistent depths
			if (applyDepthConsistency)
			  {
				if (_verbosity>3)
				  {
					_SI.EH.print("   Pre-update depth-error deviates: ");					
					computeConsistentDepthNegLogLikelihood(true);
					_SI.EH.print("%n");
				  }
				else 
				  {
					computeConsistentDepthNegLogLikelihood(false);
				  }
				// Take copy of base isopach-errors errors for finite differences
				System.arraycopy(consistent_depth_errors,0,base_depth_errors,0,consistent_depth_errors.length);
			  }
			// (3) Base vertical stacking data
			if (_SI.doStackedModel)
			  {
				if (_verbosity>3)
				  {
					_SI.EH.print("   Pre-update vertical-stacking-error deviates: ");					
					_m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(true);
					_SI.EH.print("%n");
				  }
				else 
				  {
					_m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false);
				  }
				// Take copy of base isopach-errors errors for finite differences
				System.arraycopy(vertical_stacking_errors,0,
								 base_vertical_stacking_errors,0,vertical_stacking_errors.length);
			  }

			// Start assembling the "data" vector d. First block is isopachs
			for (j=0;j<nd_start;j++) d.setQuick(j,-base_errors[j]);
			// next block is depth-consistency
			if (applyDepthConsistency) 
			  {
				for (j=0,jj=nd_start;j<_m.n_layers-1; j++,jj++) 
				  {
					d.setQuick(jj,-base_depth_errors[j]);
				  }
			  }
			// next block for vertical stacking
			if (_SI.doStackedModel)
			  {
				for (j=0,jj=nvs_start; j<vertical_stacking_errors.length; j++,jj++) 
				  {
					d.setQuick(jj,-vertical_stacking_errors[j]);
				  }
			  }

			// Loop over layers seeking all properties that will affect the isopachs or depth self-consistency.
			// Loop over columns of X, which are properties that we want an update over
			// We could do this using the map_from_index() HashMap, but we want more explicit
			// control over *which* properties we want to use in the design matrix X.
			// Se we'll loop over all layers, and all properties, being careful not to
			// create any columns in X from "mirror" properties (i.e. only their masters)
			// We'll use the following array to detect any "mirrors"
			int[] used_indices = new int[map_to.size()];
			Arrays.fill(used_indices,-1);

			for (i=0;i<_m.n_layers;i++)
			  {
				MI = _m.MPD.MetaInfo[i];
				LayerProperties LPI = _m.LP[i];

				// Here's a list of all the things that *might* affect the general error term
				// of this layer or any mirrored layers. 
				LayerProperties.Property[] prop_array =
			    	{ LPI.d, LPI.d_curves, LPI.LFIV, LPI.t, (i==_m.n_layers-1) ? LPI.t_base : null,
					  LPI.rho_b, LPI.vp_b,
					  LPI.rho_oil, LPI.rho_lsg, LPI.rho_gas,
					  LPI.vp_oil, LPI.vp_lsg, LPI.vp_gas,
					  LPI.S_oil, LPI.S_lsg, LPI.S_gas,
					  LPI.NG, LPI.vp_s, LPI.vs_s, LPI.phi_s,
					  LPI.vp_m, LPI.vs_m, LPI.rho_m, 
					  ((_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL) ? 
					   LPI.float_fraction : null),
					  ((i==0) ?_m.SM.A : null), 
					  (((i==0) && !_m.SM.A_B_avo_synchronised) ?_m.SM.B : null)
					};

				// form union of properties master list via ArrayList methods...
				prop_arrayl.clear();
				prop_arrayl.addAll((List<LayerProperties.Property>)Arrays.asList(prop_array));
				if (_SI.MF.ENABLE_HARD_ROCK_MODEL && _m.needHardRockTerms)
				  {
					prop_arrayl.addAll((List<LayerProperties.Property>)
									   Arrays.asList(LPI.NG_hard, LPI.vp_hard, LPI.vs_hard, LPI.rho_hard));
				  }
				// Get effect of systematic properties on first pass (i==0) only
				if (_m._SI.doSimpleFluidUnitsAndContacts && i==0) prop_arrayl.addAll(_m.fluidUnitsAndContacts.fc_props);

				// so find the derivative of the logLikelihood isopach & depth term wrt each of these parameters
				// by finite differences
				for (k=0; k<prop_arrayl.size(); k++)
				  {
					LayerProperties.Property prop = prop_arrayl.get(k);
					if (prop!=null && (I = (Integer)map_to.get(prop))!=null)
					  {
						// if the property in question is a mirror, skip it (we
						// don't want any duplicate columns in X)

						int iVal = I.intValue();
						if (used_indices[iVal] != -1) continue;
						else used_indices[iVal] = iVal;

						// This property is indexed, so is stochastic, and might contribute to X
						// get the column index
						int indx = I.intValue();
						// Skip any column indices out of range
						if (indx>=n) continue;

						// Add a perturbation to the property for obtaining derivatives via finite differences
						val = expansion_point.getQuick(indx);
						if (prop.normalised)
						  {
							if ((newval = val*(1.0+DESIGN_MATRIX_EPS)) > 1.0)
							  {
								// too close to one, so prop -> (1+eps)prop doesn't work.
								// try  prop -> prop/(1+eps)
								newval = val/(1.0+DESIGN_MATRIX_EPS);
							  }
						  }
						else
						  {
							newval = val*(1.0+DESIGN_MATRIX_EPS);
						  }
						if (newval == 0.0)
						  {
						    // If prop was zero, (1+eps)prop=0 so we'll change this to prop->eps
							// Will very likely only ever happen with normalised properties, so scaling is ok
							newval=DESIGN_MATRIX_EPS;
						  }
						setValuesFromIndex(indx,newval);
						if (_SI.doSimpleFluidUnitsAndContacts)
						  {
							// If there are contacts, effective properties can in principle always change
							orderTimes();  // FOR_FC
							// Recompute the effective properties for all layers, since 
							// others might be synchronised to this one. A change in times 
							// won't affect properties, so we can skip this if (prop==LPI.t)
							_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
						  }
						else if (prop!=LPI.t && prop!=LPI.t_base)
						  {
							// Recompute the effective properties for all layers, since 
							// others might be synchronised to this one. A change in times 
							// won't affect properties, so we can skip this if (prop==LPI.t)
							_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
						  }
						// compute new isopach errors (ends up in isopach_errors[])
						computeIsopachNegLogLikelihood(false);
						for (j=0;j<nd_start;j++)
						  {
							double deriv = (isopach_errors[j] - base_errors[j])/(newval-val);
							X.setQuick(j,indx,deriv);
							// and update RHS of likelihood constraint (sum of terms in X.m_0)
							d.setQuick(j, d.getQuick(j) + deriv * val);
						  }
						// depth self-consistency entries
						if (applyDepthConsistency) 
						  {
							computeConsistentDepthNegLogLikelihood(false);
							for (j=0,jj=nd_start; j<_m.n_layers-1; j++,jj++)
							  {
								double deriv = (consistent_depth_errors[j] - base_depth_errors[j])/(newval-val);
								X.setQuick(jj,indx,deriv);
								// and update RHS of likelihood constraint (sum of terms in X.m_0)
								d.setQuick(jj, d.getQuick(jj) + deriv * val);
							  }
						  }
						// vertical stacking entries
						if (_SI.doStackedModel)
						  {
							_m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false);
							for (j=0,jj=nvs_start; j<vertical_stacking_errors.length; j++,jj++)
							  {
								double deriv = (vertical_stacking_errors[j] - base_vertical_stacking_errors[j])/(newval-val);
								X.setQuick(jj,indx,deriv);
								// and update RHS of likelihood constraint (sum of terms in X.m_0)
								d.setQuick(jj, d.getQuick(jj) + deriv * val);
							  }
						  }
						// Here FC_FIXME -> get suspended_depth_top Depth sensitivity matrix?????

						// restore value
						setValuesFromIndex(indx,val);
						if (_SI.doSimpleFluidUnitsAndContacts)
						  {
							// Restore all times, restore ordering, recompute all props
							for (j=0;j<=_m.n_active_layers;j++) setValuesFromIndex(j,expansion_point.getQuick(j));
							orderTimes(); // FOR_FC
							_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
						  }
						else
						  {
							if (prop!=LPI.t && prop!=LPI.t_base)
							  {
								// Restore all layers' effective properties, since others might be synchronised to this one
								// A change in times won't affect properties, so we can skip this if (prop==LPI.t)
								_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
							  }
							else
							  {
								// if we've perturbed a time, restore all of them & make sure ordering is ok
								// number of times in model
								// nt = _m.n_active_layers-1;						  
								for (j=0;j<=_m.n_active_layers;j++) setValuesFromIndex(j,expansion_point.getQuick(j));
								orderTimes();
							  }
						  }
					  } // property contributes (non-trivial column of X)
				  } // loop over properties
			  } // loop over layers
			// now we can update the prior matrices etc to get the new prior which honours isopachs
			updateMGSystemForIsopachsAndDepthConsistency(X,d,expansion_point,NewtonIts);
			// and we may want to iterate this with a better expansion point
			expansion_point = mbar_prior_plus_isopachs[fs];

			// Rough termination criterion
			if ((NewtonIts>0) && (2.0*Math.abs(delta_fNewNewton) < NEWTON_EPS * (Math.abs(fNewNewton) + EPS_FLOATINGPOINT)))
			  {
				if (_verbosity>3) 
				  {
					_SI.EH.println("Isopachs & depth criterion: Newton terminated after "+NewtonIts+" iterations");
				  }
				break;
			  }
			else
			  {
				// not converged. If we haven't done this by maxNewtonSteps something's awry
				if (NewtonIts==(maxNewtonSteps-1) && maxNewtonSteps>1)
				  {
					if (_verbosity>3) 
					  {
						_SI.EH.println("WARNING: Isopachs & Depth criterion: unsatisfactory Newton termination criterion");  
					  }
					// Probably the mode's not worth trusting, but continue anyway
				  }
			  }
		  }
		aquireLayerTimeThicknessStandardisedUncertainties();

	  } // loop over fluid combos

	// finally, update the fluid-probabilities
	_m.FS.normaliseFluidProbabilities();
	_m.FS.trimPossibleCombinationsFromThresholds(false); 
	if (_verbosity>2)
	  {
		_SI.EH.print("Updated fluid combination probabilities from isopachs"+
						(applyDepthConsistency ? " and depth-self-consistency..." : "..."));
		_m.FS.printFluidCombos();
	  }
  }

  
  private boolean aquireLayerTimeThicknessStandardisedUncertainties()
  {
	int i,j,k;
	double val,s;
	double R_SPARSITY_THRESHOLD = 0.2;
	DoubleMatrix2D R;
	DoubleMatrix1D m_tilde = mbar_prior_plus_isopachs[fs];
	Integer I;

	// get time-covariance submatrix
	// compute this from Cholesky factor *after* isopach, depth consistency etc is incorporated
	R = R_prior_plus_isopachs[fs];
	// Don't compute full covariance matrix; only the diagonal as a vector (Ctimes_diag)
	// and the upper diagonal (covariances for adjacent master layers)
	double[] Ctimes_diag = new double[_m.n_active_layers+1];
	double[] Ctimes_upper_diag = new double[_m.n_active_layers];

	// compute covariance matrix diagonal and upper diagonal for times
	if (R.cardinality() < R_SPARSITY_THRESHOLD * R.size())
	  {
		// Efficient code for sparse R
		int kMin,colIMin;
		MatrixColumnInfo info = new MatrixColumnInfo(R);
		int[] colMins = info.colmins;
		for (i=0;i<=_m.n_active_layers;i++)
		  {
			colIMin = colMins[i];
			for (j=Math.max(0,(i-1)); j<=i; j++)
			  {
				val=0.0;
				kMin = Math.max(colIMin, colMins[j]);
				for (k=kMin;k<=j;k++) val += R.getQuick(k,i) * R.getQuick(k,j);
				if (j==i) Ctimes_diag[i] = val;
				if (j==i-1) Ctimes_upper_diag[j] = val;
			  }
			// set  master-layer time standard deviations
			t_layer_uncertainty.setQuick(fs, i, Math.sqrt(Ctimes_diag[i]));
		  }
		info.colmins = null;
		info = null;
	  }
	else
	  {
		// Compute as if R were dense
		for (i=0;i<=_m.n_active_layers;i++)
		  {
			for (j=Math.max(0,(i-1)); j<=i; j++)
			  {
				val=0.0;
				for (k=0;k<=j;k++) val += R.getQuick(k,i) * R.getQuick(k,j);
				if (j==i) Ctimes_diag[i] = val;
				if (j==i-1) Ctimes_upper_diag[j] = val;
			  }
			// set master-layer time standard deviations
			t_layer_uncertainty.setQuick(fs, i, Math.sqrt(Ctimes_diag[i]));
		  }
	  }

	// accumulate layer-time std deviations, since we have updated the Covariance
	// and Cholesky factors in updateMGSystemForIsopachsAndDepthConsistency()
	if (_verbosity>3) _SI.EH.printf("Fluid state %d, master-layer time uncertainties (post isopach+depth): ",fs);
	for (i=0;i<=_m.n_layers;i++)
	  {
		// last entry in time variances is for t_base
		if ((I = (Integer)map_to.get((i==_m.n_layers) ? _m.LP[i-1].t_base : _m.LP[i].t))!=null)
		  {
			j = I.intValue();  // index in linear algebra
			s = Math.sqrt(Ctimes_diag[j]);  // time std deviation
			boolean cont=true;
			SortableNumbers SN = new SortableNumbers(s,i);
			do
			  {
				// There's a small chance of duplication, so test to be sure...
				if (time_std_devs[fs].get(SN)==null)
				  {
					time_std_devs[fs].put(SN, new Integer(i));
					cont=false;
					break;
				  }
				// if it's already there, wind up the variance a bit to break the duplication
				// this ensures , when there are duplicates, upper layers are set first
				SN.changeValue((s*=1.00001));
			  }	while(cont);	
			if (_verbosity>3) _SI.EH.printf("%10.4g ",s);
		  }
	  }
	if (_verbosity>3) _SI.EH.println("");

	// Assess relative time-horizon separations

	// and the relative layer time-thickness uncertainty

	if (_verbosity>3) _SI.EH.println("Normalised master-layer-time-thickness uncertainties of mode:");

	for (i=0;i<_m.n_active_layers;i++)
	  {
		s = Math.sqrt(Ctimes_diag[i] + Ctimes_diag[i+1] 
					  - 2.0* Ctimes_upper_diag[i]);
		// normalise
		val = Math.abs( m_tilde.getQuick(i+1) - m_tilde.getQuick(i) );
		time_uncertainty.setQuick(fs, i, s/((val>0.0) ? val : 1.0));
		if (_verbosity>3) _SI.EH.printf("%10.4g ",time_uncertainty.getQuick(fs,i));
	  }
	if (_verbosity>3) _SI.EH.println("");
	
	Ctimes_diag = null;
	Ctimes_upper_diag = null;
	return true;
  }
  /**
	 Same as aquireLayerTimeThicknessStandardisedUncertainties() 
	 but for layer upper model depths
  */
  boolean aquireLayerTopDepthStandardisedUncertainties()
  {
	int i,j,k;
	double val,s;
	double R_SPARSITY_THRESHOLD = 0.2;
	DoubleMatrix2D R;
	DoubleMatrix1D m_tilde = mbar_prior_plus_isopachs[fs];
	Integer I;

	// get time-covariance submatrix
	// compute this from Cholesky factor *after* isopach, depth consistency etc is incorporated
	R = R_prior_plus_isopachs[fs];
	// Don't compute full covariance matrix; only the diagonal as a vector (Ctimes_diag)
	// and the upper diagonal (covariances for adjacent master layers)
	double[] Ctimes_diag = new double[_m.n_active_layers+1];
	double[] Ctimes_upper_diag = new double[_m.n_active_layers];

	// compute covariance matrix diagonal and upper diagonal for times
	if (R.cardinality() < R_SPARSITY_THRESHOLD * R.size())
	  {
		// Efficient code for sparse R
		int kMin,colIMin;
		MatrixColumnInfo info = new MatrixColumnInfo(R);
		int[] colMins = info.colmins;
		for (i=0;i<=_m.n_active_layers;i++)
		  {
			colIMin = colMins[i];
			for (j=Math.max(0,(i-1)); j<=i; j++)
			  {
				val=0.0;
				kMin = Math.max(colIMin, colMins[j]);
				for (k=kMin;k<=j;k++) val += R.getQuick(k,i) * R.getQuick(k,j);
				if (j==i) Ctimes_diag[i] = val;
				if (j==i-1) Ctimes_upper_diag[j] = val;
			  }
			// set  master-layer time standard deviations
			t_layer_uncertainty.setQuick(fs, i, Math.sqrt(Ctimes_diag[i]));
		  }
		info.colmins = null;
		info = null;
	  }
	else
	  {
		// Compute as if R were dense
		for (i=0;i<=_m.n_active_layers;i++)
		  {
			for (j=Math.max(0,(i-1)); j<=i; j++)
			  {
				val=0.0;
				for (k=0;k<=j;k++) val += R.getQuick(k,i) * R.getQuick(k,j);
				if (j==i) Ctimes_diag[i] = val;
				if (j==i-1) Ctimes_upper_diag[j] = val;
			  }
			// set master-layer time standard deviations
			t_layer_uncertainty.setQuick(fs, i, Math.sqrt(Ctimes_diag[i]));
		  }
	  }

	// accumulate layer-time std deviations, since we have updated the Covariance
	// and Cholesky factors in updateMGSystemForIsopachsAndDepthConsistency()
	if (_verbosity>3) _SI.EH.printf("Fluid state %d, master-layer time uncertainties (post isopach+depth): ",fs);
	for (i=0;i<=_m.n_layers;i++)
	  {
		// last entry in time variances is for t_base
		if ((I = (Integer)map_to.get((i==_m.n_layers) ? _m.LP[i-1].t_base : _m.LP[i].t))!=null)
		  {
			j = I.intValue();  // index in linear algebra
			s = Math.sqrt(Ctimes_diag[j]);  // time std deviation
			boolean cont=true;
			SortableNumbers SN = new SortableNumbers(s,i);
			do
			  {
				// There's a small chance of duplication, so test to be sure...
				if (time_std_devs[fs].get(SN)==null)
				  {
					time_std_devs[fs].put(SN, new Integer(i));
					cont=false;
					break;
				  }
				// if its already there, wind up the variance a bit to break the duplication
				// this ensures , when there are duplicates, upper layers are set first
				SN.changeValue((s*=1.00001));
			  }	while(cont);	
			if (_verbosity>3) _SI.EH.printf("%10.4g ",s);
		  }
	  }
	if (_verbosity>3) _SI.EH.println("");

	// Assess relative time-horizon separations

	// and the relative layer time-thickness uncertainty

	if (_verbosity>3) _SI.EH.println("Normalised master-layer-time-thickness uncertainties of mode:");

	for (i=0;i<_m.n_active_layers;i++)
	  {
		s = Math.sqrt(Ctimes_diag[i] + Ctimes_diag[i+1] 
					  - 2.0* Ctimes_upper_diag[i]);
		// normalise
		val = Math.abs( m_tilde.getQuick(i+1) - m_tilde.getQuick(i) );
		time_uncertainty.setQuick(fs, i, s/((val>0.0) ? val : 1.0));
		if (_verbosity>3) _SI.EH.printf("%10.4g ",time_uncertainty.getQuick(fs,i));
	  }
	if (_verbosity>3) _SI.EH.println("");
	
	Ctimes_diag = null;
	Ctimes_upper_diag = null;
	return true;
  }

  private void updateMGSystemForIsopachsAndDepthConsistency(DoubleMatrix2D X, DoubleMatrix1D d, 
															DoubleMatrix1D expansion_point, int iteration)
  {
	// Use the standard linearised posterior bayesian Updating formulas to get the posterior
	// after accounting for isopachs and/or depth-self-consistency
	// in the writeup (standard inverse theory) to get new MG distribution
	int i,j,k,n,n_params;
	double val,rhs,s,A_term,det;
    ModifiedCholeskyRightDecomposition CDR;
	DoubleMatrix2D Xtilde,Rinv;
	DoubleMatrix1D m_tilde;

	DoubleMatrix2D Cinv;
	DoubleMatrix1D mbar = mbar_prior_only[fs];

	// must start with a fresh copy of the unpolluted prior
	Cinv = Cinv_pr[fs];
	Cinv.assign(Cinv_prior_only[fs]);


	m_tilde = mbar_prior_only[fs].copy();
	boolean sample_succeded;

	boolean debug=false; // FIXME turn off & ultimately remove


	// This routine will overwrite Cinv with Ctilde.
	// Before we do that, we'll need the determinant of Cinv
	// Get this via the diagonal of the Cholesky decomposition
	CDR = new ModifiedCholeskyRightDecomposition(Cinv);
	Rinv = CDR.getR();
	R_prior_only[fs] = Rinv.copy(); // we will need duplicate
	// At this point, R_prior_only[fs] is an upper triangular Cholesky factor of Cinv = C^{-1}

	// the diagonal contains all we need for the determinant
	n_params = Rinv.rows();
	double[] diagRinv_pr = new double[n_params];
	for (j=0;j<n_params;j++)
	  {diagRinv_pr[j] = Rinv.getQuick(j,j);}
	Rinv = null;

	n = X.rows();


	if (debug && fs==0)   // FIXME COMMENT OUT
	  {
		Util.debugPrintMatrixToFile("Cinv.txt",Cinv);
		Util.debugPrintMatrixToFile("mbar.txt",mbar_prior_only[fs]);
		Util.debugPrintMatrixToFile("d.txt",d);
	  }
	// form Xtilde = X^T and A_term = d^T.d
	// (A_term is acquired here for updating overall fluid probability update factor)
	A_term = d.zDotProduct(d,0,d.size());

	// form rhs = C^-1.mbar and  update A_term += mbar^T.C^-1.mbar
	SmpBlas.smpBlas.dsymv(false,1.0,Cinv,mbar,0.0,m_tilde);
	A_term += mbar.zDotProduct(m_tilde,0,n_params);
	// At this point m_tilde = C^-1.m_bar
	
	if (debug && fs==0)  // FIXME COMMENT OUT
	  {
		Util.debugPrintMatrixToFile("X.txt",X);
	  }
	// Form C^-1   <-   X^T.X + C^-1
	DoubleMatrix1D col_i,col_j;
	for (i=0 ; i < Cinv.rows() ; i++)
	  {
		for (j=0 ; j <= i; j++)
		  {
			val = Cinv.getQuick(i,j);
			col_i = X.viewColumn(i);
			col_j = X.viewColumn(j);
			val += col_i.zDotProduct(col_j,0,n);

			Cinv.setQuick(i,j,val);
			if (i!=j) Cinv.setQuick(j,i,val); // use symmetry to set this one
		  }
	  }
	// form rhs = (X^T.d + C^-1.m_bar) = X^T.d + m_tilde *in place*
	// (at this point m_tilde = C^-1.m_bar. Add  X^T.d bit
	SmpBlas.smpBlas.dgemv(true,1.0,X,d,1.0,m_tilde);

	if (false)
	  {
		// test of gradient
		// at this point m_tilde contains xi = X^T.d + C^-1.mbar 
		DoubleMatrix1D gradF = m_tilde.copy();
		// form -gradF = xi - {\tilde C}_{s}^{-1}.m_0
		SmpBlas.smpBlas.dgemv(false,-1.0,Cinv,expansion_point,1.0,gradF);
		Util.debugPrintMatrixToFile("gradF_iso.txt",gradF);
	  }

	// Now at this point Cinv contains {\tilde C}_m^{-1}
	// we may very likely need a copy of {\tilde C}_m^{-1}
	Cinv_prior_plus_isopachs[fs] = Cinv;
			
	// And m_tilde contains (X^T.d + C^-1.mbar) "the RHS"
	DoubleMatrix1D RHS = m_tilde.copy();

	// get right-Cholesky factor of Cinv
	CDR = new ModifiedCholeskyRightDecomposition(Cinv); 
	R_prior_plus_isopachs[fs] = CDR.getR();
	
	

	CDR.getCholeskyRightInverse(R_prior_plus_isopachs[fs]);
	// At this point R_prior_plus_isopachs is the transpose of the upper triangular R Cholesky
	// factor of {\tilde C}, viz {\tilde C} = R_prior_plus_isopachs[fs]^T . R_prior_plus_isopachs[fs]
	
	Rinv=R_prior_plus_isopachs[fs];
	int nC = Cinv.rows();
	// form {\tilde C}.m_tilde =  (Rinv^T.Rinv).m_tilde *in place*
	// (1) form m_tilde <- Rinv.m_tilde *in place*
	for (i=0;i<nC;i++)
	  {
		s = 0;
		for (j=i;j<nC;j++)
		  {
			s += Rinv.getQuick(i,j) * m_tilde.getQuick(j);
		  }
		// write over m_tilde in place in last element no longer needed)
		m_tilde.setQuick(i,s);
	  }
	

	// (2) now Rinv^T.m_tilde and the determinant ratio |Ctilde|/|C_m|
	// and also the (log)determinant of the new prior Ctilde
	// det=1.0;
	logDetC_prior_only[fs]=0.0;
	detC_prior_plus_isopachs[fs]=0.0;
	for (i=nC-1;i>=0;i--)
	  {
		s = 0;
		for (j=0;j<=i;j++)
		  {
			// notice transpose!
			s += Rinv.getQuick(j,i) * m_tilde.getQuick(j);
		  }
		// write over m_tilde in place in ith element as no longer needed)
		m_tilde.setQuick(i,s);
		// update the term bit too, to subtract off (mbar^T.C_m^-1+d^T.Ciso^-1.X).m_tilde = RHS^T.m_tilde
		A_term -= s * RHS.getQuick(i);
		// and also the determinant ratio. No need for logs as these should numerically stabilise...
		// det *=  Rinv.getQuick(i,i) * diagRinv_pr[i];
		// we'll need the determinant of the prior covariance (subtract since R is an inverse)
		logDetC_prior_only[fs] -= Math.log(MathUtil.sqr(diagRinv_pr[i]));
		// and also the determinant of the updated prior covariance {\tilde C}
		detC_prior_plus_isopachs[fs] += Math.log(MathUtil.sqr(R_prior_plus_isopachs[fs].getQuick(i,i)));
	  }

	// Note: update overall fluid probability update factor. This is left till later..
	// fluid_state.prob *= Math.exp(-A_term/2.0) * det;

	diagRinv_pr=null;
	RHS=null;



	// Test Gassman validity of new point and get best point along line from
	// expansion_point to m_tilde (hopefully will be the end, m_tilde)
	m_tilde = getBestLegalGassmanPointAlongCurrentLine(expansion_point, m_tilde);

	// update now the new multiGaussian means
	mbar_prior_plus_isopachs[fs] = m_tilde;


	// Util.debugPrintMatrixToFile("m_tilde.txt",m_tilde);

	//  Now compute mean isopach & depth deviations to verify that linearisation is valid
	
	// set all times to updated (prior+isopach) means
	// _SI.EH.print("Updated mean times to compute isopachs");
	for (i=0;i<n_params;i++) setValuesFromIndex(i,m_tilde.getQuick(i));

	orderTimes();  // FOR_FC
	// now update full layer properties
	sample_succeded = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
		
	// was it unphysical??
	fNewNewton=0.0;
	if (sample_succeded)
	  {
		// All is well...
		if (!_SI.ignoreIsopachs && _m.MPD.n_isopachs>0)
		  {
			// compute the scaled isopach error vector
			fNewNewton += computeIsopachNegLogLikelihood(false);
			if (_verbosity>3)
			  {
				_SI.EH.print("   Post-update mean isopach deviates: ");
				for (k=0;k<_m.MPD.n_isopachs;k++)
				  {
					_SI.EH.printf("%10.4g ",isopach_errors[k]);
				  }
				_SI.EH.print("%n");
			  }
		  }
		// Check depth self-consistencies too if required
		if (applyDepthConsistency) 
		  {
			fNewNewton += computeConsistentDepthNegLogLikelihood(false);
			if (_verbosity>3)
			  {
				_SI.EH.print("   Post-update depth-error deviates: ");
				for (k=0; k<_m.n_layers-1; k++)
				  {
					_SI.EH.printf("%10.4g ",consistent_depth_errors[k]);
				  }
				_SI.EH.print("%n");
			  }
		  }
		if (_SI.doStackedModel)
		  {
			fNewNewton += _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false);
			if (_verbosity>3)
			  {
				_SI.EH.print("   Post-update vertical-stacking-error deviates: ");
				for (k=0; k<vertical_stacking_errors.length; k++)
				  {
					_SI.EH.printf("%10.4g ",vertical_stacking_errors[k]);
				  }
				_SI.EH.print("%n");
			  }
		  }
	  }
	else
	  {
		// a failure shouldn't happen or things are badly suspect
		// getBestLegalGassmanPointAlongCurrentLine() call will almost certainly prevent this
		throw new RuntimeException("Mean isopach deviations force impossible fluid properties."
								   + "in fluid state " +
								   (fs - 1) +
								   ". Regression curves, isopachs, or fluid properties " +
								   "are probably suspect");
	  }

	// finally, incorporate prior term
	fNewNewton += computePriorNegLogLikelihood(m_tilde);

	// set Newton variables
	delta_fNewNewton = fNewNewton - fNewNewtonOld;
	fNewNewtonOld = fNewNewton;

	if (debug && fs==0)  // COMMENT OUT
	  {
		Util.debugPrintMatrixToFile("L.txt",Rinv);		
		Util.debugPrintMatrixToFile("mu_tilde.txt",m_tilde);
	  }
  }

  /**
	 Draw samples out of the prior model which incorporates prior + isopachs
  */
  void drawSetOfPriorSamples(int nSamplesPerFluidCombo, boolean showIsopachErrors)
  {
	int i_s,k=1;

	// loop over fluid combinations, obtaining fully defined prior MG system
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		// loop for nSamplesPerFluidCombo samples
		for (i_s=0;i_s<nSamplesPerFluidCombo;i_s++)
		  {
			drawPriorSample(showIsopachErrors);
			if (_SI.dumpRealisations) _SI.dumpRealisation(true,k,true,false);
		    negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
			negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
			_m.SM.dumpSynthetics(true,k++,true);
		  }
	  }
  }

  /**
	 Method to draw samples from MG distribution incorporating priors,
	 depth self-consistency, & isopach requirements, for given
	 fluid-state fs
	 @param showIsopachErrors Prints isopach standardised errors if true
  */
  private void  drawPriorSample(boolean showIsopachErrors)
  {
	int i,j,k,i_s,fail_count;
	boolean sample_succeded;
	DoubleMatrix2D Ltransp;
	Ltransp = R_prior_plus_isopachs[fs];
	mu = mbar_prior_plus_isopachs[fs];
	int n = Ltransp.rows();

	setMapsFromFluidState(fs);

	// make vector of normal deviates
	DoubleMatrix1D w = new DenseDoubleMatrix1D(n);
	// multiply by Cholesky "L" to get correlated true properties
	// The R factor we have is the inverse of the R from R.R^T = C^{-1}, so is
	// actually of form (R^{-1})^T . R^{-1} = C. So if C = L. L^T,
	// we have in place R = L^T. So x = L.w = R^T.w is the random vector we want

	fail_count=0;
	do
	  {
		// make a sample
		for (i=0;i<n;i++) w.setQuick(i, Ngen.nextDouble());
   		for (i=0;i<n;i++)
		  {
			double s = 0;
			// From inverse of right cholesky decomposition,
			// Ltransp has true left cholesky factor, but transposed.
			for (j=0;j<=i;j++)
			  {
				// notice transpose!
				s += Ltransp.getQuick(j,i) *  w.getQuick(j);
			  }
			// add back the mean
			s += mu.getQuick(i);
			// s now has the random deviate.
			// get the property for this index
			setValuesFromIndex(i,s);
		  }
		orderTimes();
		// now update full layer properties
		sample_succeded = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
		// was it physical??
		if (sample_succeded)
		  {
			// All is well...
			fail_count=0;
			if (showIsopachErrors && _m.MPD.n_isopachs>0)
			  {
				computeIsopachNegLogLikelihood(true);
				_SI.EH.printf("Fluid state %d, layers: ",fs);
				for (i=0;i<_m.n_layers;i++) 
				  {
					_SI.EH.printf("%1d ",(layer_present[i] ? 1 : 0));
				  }
				_SI.EH.printf("code %2d isopach deviates:",layer_code); 
				for (k=0;k<_m.MPD.n_isopachs;k++)
				  {
					_SI.EH.printf("%10.4g ",isopach_errors[k]);
				  }
				_SI.EH.print("%n");
			  }
		  }
		// a failure shouldn't happen too many times or things are badly suspect
		if (fail_count>MAX_MC_FAIL_COUNT)
		  {
			throw new RuntimeException("Impossible fluid or rock properties in fluid state " +
									   (fs - 1) +
									   ". Regression curves, isopachs, or fluid properties " +
									   "are probably suspect");
		  }
	  } while (sample_succeded==false);

	// if we got here, the sample is ok. Order the times
	// Impose time truncations
	orderTimes();
  }














  /**
	 Run over sampled times to ensure that they are monotonic.
	 private orderTimes(). Returns sum of time-crossovers
	 as an error measure (not used, but may be useful)
  */
  double orderTimes()
  {
	int i;
	boolean[] isSet=new boolean[_m.n_layers+1]; // defaults to false
	double t_below;
	LayerProperties.Property t_this;
	double crossover_time_sum=0.0;
	Iterator I = time_std_devs[fs].values().iterator();
	int J = time_std_devs[fs].size();

	// fragment to show in-situ layer times...
	if (_verbosity>5)
	  {
		_SI.EH.print("Pre-ordering  layer times: ");
		for (i=0;i<=_m.n_layers;i++) 
		  {
			_SI.EH.printf("%9.7g ",
							((i==_m.n_layers) ?  _m.LP[i-1].t_base.value :  _m.LP[i].t.value));
		  }
		_SI.EH.println("");
	  }

	while (I.hasNext())
	  {
		int indx = ((Integer)I.next()).intValue();
		// the layer time in question
		t_this = (indx==_m.n_layers) ?  _m.LP[indx-1].t_base :  _m.LP[indx].t;
		// look over existing fixed times for crossovers...
		for (i=indx-1;i>=0;i--)
		  {
			if (isSet[i])
			  {
				if (_m.LP[i].t.value > t_this.value)
				  {
					// we hit a cross over...Fix this at the lowest layer
					crossover_time_sum += _m.LP[i].t.value - t_this.value;
					t_this.value = _m.LP[i].t.value;
					isSet[indx]=true;
					break;
				  }
				else
				  {// layer above is ok, so we can bail out
					break;
				  }
			  }
		  }
		// look down for crossovers...
		for (i=indx+1;i<=_m.n_layers;i++)
		  {
			t_below = (i==_m.n_layers) ?  _m.LP[i-1].t_base.value :  _m.LP[i].t.value;
			if (isSet[i])
			  {
				if (t_below < t_this.value)
				  {
					// we hit a cross over...Fix this at the lowest layer
					crossover_time_sum += t_this.value - t_below;
					t_this.value = t_below;
					isSet[indx]=true;
					break;
				  }
				else
				  {// layer below is ok, so we can bail out
					break;
				  }
			  }
		  }
		// if no cross overs, the current value is legal. Indicate so with the flag
		isSet[indx]=true;
	  }


	// set the state of the layer_present state variable
	for (i=0;i<_m.n_layers-1;i++)
	  {
		layer_present[i] = (_m.LP[i+1].t.value > _m.LP[i].t.value); 
	  }
	layer_present[_m.n_layers-1]= (_m.LP[_m.n_layers-1].t_base.value > _m.LP[_m.n_layers-1].t.value);
	
	// compute a pseudo-binary representation of the layers present
	layer_code=0;
	for (i=0;i<_m.n_layers;i++) if (layer_present[i]) layer_code += (1<<i);


	if (_verbosity>5)
	  {
		_SI.EH.print("Post-ordering layer times: ");
		for (i=0;i<=_m.n_layers;i++) 
		  {
			_SI.EH.printf("%9.7g ",
							((i==_m.n_layers) ?  _m.LP[i-1].t_base.value :  _m.LP[i].t.value));
		  }
		_SI.EH.println("");
	  }

	isSet=null;
	return crossover_time_sum;
  }

  void computeMostLikelyPriorFluidState()
  {
	int fs_max=0;
	double prob_max = -1;
	for (fs = 0 ; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		fluid_state = (FluidStates.FluidState) _m.FS.possibleFluidStates.get(fs);
		if (fluid_state.prob > prob_max)
		  {
			fs_max = fs;
			prob_max = fluid_state.prob;
		  }
	  }
	_m.FS.fs_max = fs_max;
  }

  void drawMostLikelyPrior(boolean showIsopachErrors)
  {
	int k;
	boolean sample_succeeded=false;
	ModelPriorDescriptions.MetaInfo MI;
	LayerProperties LPI;
 
	// loop over fluid combinations, obtaining most likely fluid combo
	// sets _m.FS.fs_max
	computeMostLikelyPriorFluidState();

	// set the current fluid state to the most likely
	fs = _m.FS.fs_max;
	fluid_state = (FluidStates.FluidState) _m.FS.possibleFluidStates.get(_m.FS.fs_max);
	if (_verbosity>2)
	  {
		_SI.EH.println("Creating model for most likely fluid state: (fluid state %d),",_m.FS.fs_max);
		_SI.EH.println("incorporating linearised isopach and depth-consistency terms if applicable");
	  }

	// set properties at mean from isopach and depth-consistency update
	mu = mbar_prior_plus_isopachs[_m.FS.fs_max];
	map_from = map_from_index[_m.FS.fs_max];
	map_to = map_to_index[_m.FS.fs_max];
	int n = mu.size();
	// set relevant variable properties at their modes...
	for (k=0;k<n;k++)
	  {
		double s = mu.getQuick(k);
		setValuesFromIndex(k,s);
	  }
	orderTimes();
	// now update full layer properties for effective vp's etc
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK
  
	if (_SI.dumpGassmanDiagnostics) 
	  {
		showGassmanDiagnostics(false,true);
		showGassmanDiagnostics(false,false);
	  }

	// All is well...
	if (showIsopachErrors && _verbosity>2)
	  {
		if (_m.MPD.n_isopachs>0)
		  {
			_SI.EH.printf("Isopach errors: ",_m.FS.fs_max);
			for (k=0;k<_m.n_layers;k++)
			  {
				MI = _m.MPD.MetaInfo[k];
				LPI = _m.LP[k];
				if ((MI.dz != null) && !(MI.dz.value==0 && MI.sigma_dz.value==0)) // layer has an non trivial isopach
				  {
					// layer time thickness...
					double t_base = (k==_m.n_layers-1) ? LPI.t_base.value : _m.LP[k+1].t.value;
					// and => thickness. Factor of 1/2 since times are two-way times
					double DZ = LPI.vp_eff * 0.5 * (t_base - LPI.t.value);
					double dev = (DZ-MI.dz.value)/MI.sigma_dz.value;
					_SI.EH.printf("%10.4g ",dev);
				  }
			  }
			_SI.EH.print("%n");	
		  }

		if (applyDepthConsistency)
		  {
			_SI.EH.printf("Depth-self-consistency errors: ",_m.FS.fs_max);
			computeConsistentDepthNegLogLikelihood(true);
			_SI.EH.print("%n");	
		  }

		if (_SI.doStackedModel)
		  {
			_SI.EH.printf("Vertical-stacking errors: ",_m.FS.fs_max);
			_m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(true);
			_SI.EH.print("%n");	
		  }

		_SI.EH.print("Layer times: ");
		for (k=0;k<=_m.n_layers;k++)
		  {
			_SI.EH.printf("%10.8g ",((k<_m.n_layers)? _m.LP[k].t.value :_m.LP[k-1].t_base.value ));
		  }
		_SI.EH.print("%n");
	  }
  }

  void resetModelPriorTimesFromInSituState()
  {
	int i,i_inc,k,fs_l;
	double dz,dt;
	boolean descending,is_time_synched;
	double[] adjusted_times;
	ModelPriorDescriptions.MetaInfo MI;

	// The following wraps up these two loops: 
	// From master layer downwards, then master layer upwards...

	adjusted_times = new double[_m.n_layers+1];
	adjusted_times[_m.master_depth_layer_number-1] = _m.LP[_m.master_depth_layer_number-1].t.value;
	for (i_inc=1; i_inc>=-1; i_inc -=2)
	  {
		for (i = ((i_inc==1) ? _m.master_depth_layer_number-1 : _m.master_depth_layer_number-2); 
			 (i_inc==1)? (i<_m.n_layers) : (i>=0) ;        i += i_inc)			 
		  { 
			descending = (i_inc == 1);
			LPI = _m.LP[i];
			MI = _m.MPD.MetaInfo[i];
			// default thickness from existing times is: 
			dz = LPI.suspended_depth_bottom - LPI.suspended_depth_top; 
			// but use target thickness! if supplied 
			if (MI.dz != null) dz = MI.dz.value;  
			dt = (LPI.vp_eff > 0) ? (2.0*dz/LPI.vp_eff) : 0.0;
			if (descending) 
			  {
				adjusted_times[i+1] = adjusted_times[i] +dt;
			  }
			else
			  {
				adjusted_times[i] = adjusted_times[i+1] -dt;
			  }
		  }
	  }
	if (_verbosity>2) 
	  {
		_SI.EH.println("UPDATED layers times from depths:");
		_SI.EH.println("                            t_old      t_new");
		for (i=0;i<_m.n_layers;i++)
		  {
			_SI.EH.printf("Layer %5d:        %15.6g %15.6g%n",(i+1),
						  _m.MPD.MetaInfo[i].t.value,adjusted_times[i]);
		  }
		    _SI.EH.printf("Layer %5d, t_base %15.6g %15.6g%n",_m.n_layers,
						  _m.MPD.t_base.value,adjusted_times[i]);
	  }


	// rejig priors, prior vectors etc. Same fragment as from setUpPriorMatrices()
	// loop over fluid combinations & reset mean for each.
	for (fs_l = 0 ; fs_l < _m.FS.possibleFluidStates.size();fs_l++)
	  {
		map_to = map_to_index[fs_l];
		mu = mbar_prior_only[fs_l];
		for (i=0;i<_m.n_layers;i++)
		  {
			is_time_synched = (i==0)? false : !_m.MPD.MetaInfo[i-1].isActive;
			if (!is_time_synched)
			  {
				// modify prior mean only
				modifyPriorMGMuEntry(_m.LP[i].t, adjusted_times[i]);
				_m.MPD.MetaInfo[i].t.value = adjusted_times[i];
			  }
		  }
		// tbase
		if (_m.MPD.sigma_t_base.value > 0)
		  {
			modifyPriorMGMuEntry(_m.LP[_m.n_layers-1].t_base, adjusted_times[_m.n_layers]);
			_m.MPD.t_base.value = adjusted_times[i];
		  }
	  }
	adjusted_times = null;
	// point again at best fluid state
	fs = _m.FS.fs_max;
	mu = mbar_prior_plus_isopachs[_m.FS.fs_max];
	map_to = map_to_index[_m.FS.fs_max];
  }


  /**
	 Given a local mode, sets all HashMap references 
	 {fluid_state, map_from, map_to} to point to the correct place. These
	 references are heavily used in subsequent operations on the mode
  */
  void setMode(LocalMode new_mode)
  {
	mode = new_mode;

	// set correct fluid state and working variables accordingly
	fs = mode.fs;
	fluid_state = (FluidStates.FluidState) _m.FS.possibleFluidStates.get(fs);
	map_from = map_from_index[fs];
	map_to = map_to_index[fs];
  }

  /**
	 Given a fluid state, sets all HashMap references 
	 {fluid_state, map_from, map_to} to point to the correct place. These
	 references are heavily used in subsequent operations on the mode
  */
  void setMapsFromFluidState(int fs)
  {
	// set correct fluid state and working variables accordingly
	fluid_state = (FluidStates.FluidState) _m.FS.possibleFluidStates.get(fs);
	map_from = map_from_index[fs];
	map_to = map_to_index[fs];
  }



  void enumerateAllLocalModes()
  {
	int i,j,k,kk, warnCount=0;
	double[] timesPlusOtherParameters = null;
	double norm,prior_norm,weight_norm,prior_weight_norm,tmp;
	DoubleMatrix2D R_times=null;
	DoubleMatrix1D t_uncert,dt_uncert;
	ArrayList modeIdentifiers;
	ArrayList nonTimeProperties;

	boolean useOldParameters = true;       // Old style only ever used shale properties, and for all NG values
	boolean searchAlongRockCurves = false;  // Experimental idea - doesn't work as well as one might hope

	// ArrayList containing full suite of local modes, and all their
	// requisite properties: fluid state, marginal probability, local mean, local covariance etc


	localModes = new ArrayList<LocalMode>();

	// Map of local minima found, used to prevent adding duplicate modes
	modeIdentifiers = new ArrayList();

	// Define a new subvector of internal-reflector times only, to optimise over.
	// IE exclude top time and base time. This gives N = _m.n_active_layers-1 unknowns
	// but times has to be unit-offset, i.e. times[1...N], 
	// (so times.length = N+1 =_m.n_active_layers if there are no non-time parameters)
	// Hence the first 1...N values in timesPlusOtherParameters are times.
	// times = new double[_m.n_layers];
	nonTimeProperties = new ArrayList();
	timeFinder = new OptimalTimesFinder(_m,nonTimeProperties,
										_SI.reParameteriseToAdditiveTimes);

	for (fs = 0; fs < _m.FS.possibleFluidStates.size();fs++)
	  {
		if (_verbosity>2) _SI.EH.println("Finding modes for fluid state "+fs+"...");
		setMapsFromFluidState(fs);

		mu = mbar_prior_plus_isopachs[fs];
		dt_uncert = time_uncertainty.viewRow(fs);
		t_uncert = t_layer_uncertainty.viewRow(fs);

		//.........	ADDING OTHER PROPERTIES
		/*
		  The mode-finding method seeks for optimal times, with (perhaps)
		  everything else fixed.  It can also optimise over selected
		  non-time properties. Here we add some good choices to the
		  parameter vector. Sand/Shale velocities (vp), and NG seem good choices.
		*/
		warnCount=0;
		if (!_SI.timesOnlyInOptimiser)
		  {
			Integer I,Is,Iden;
			HashMap used_indices = new HashMap();
			// must cleanly assemble extra properties for each fluid state
			nonTimeProperties.clear();
			LayerProperties LPI=null;
			for (i=0;i<_m.n_layers;i++)
			  {
				LPI = _m.LP[i];
				// rock properties.
				// (1) dominant rock p-wave velocity and density
				// [ old style always used shale props even in sand-rich layers ]
				if ((_m.MPD.MetaInfo[i].NG.value <= 0.5) || useOldParameters) 
				  {
					// shale dominated	case	
					if (useOldParameters && (_verbosity>3) && (warnCount==0)) 
					  {
						_SI.EH.println("Using old optimiser: subspace sweep on shale props only");
						warnCount++;
					  }
					if (searchAlongRockCurves) // TURNED OFF. See above. Doesn't work so well. Also not upgraded for hard rocks
					  {
						// Add linked property mapping vp/vs/rho along maximum-prior prob. line
						I =  (Integer) map_to.get(LPI.vp_m);                   // Vp
						if (I != null && used_indices.get(I) == null)
						  {
							Is =  (Integer) map_to.get(LPI.vs_m);              // Vs
							if (Is != null && used_indices.get(Is) == null)
							  {
								Iden =  (Integer) map_to.get(LPI.rho_m);       // rho
								if (Iden != null && used_indices.get(Iden) == null)
								  {
									// Ok, vp and vs and rho are available free properties
									// Regression data for this shale....
									RockandFluidData.RegressionSet CS;
									CS  = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(_m.MPD.MetaInfo[i].mud_rock_member);
										
									// Add ordinary ABS mappers for vp, vs and rho , but couple vs and rho to vp
										
									// parent
									PropertyMapper pMapper_vp = new PropertyMapper(I.intValue(),PropertyMapper.ABS);
									//children
									PropertyMapper pMapper_vs = new PropertyMapper(Is.intValue(),PropertyMapper.ABS);
									PropertyMapper pMapper_rho = new PropertyMapper(Iden.intValue(),PropertyMapper.ABS);
									// entries
									pMapper_vp.makeCoupled(2);
									pMapper_vp.coupledProperties[0] = pMapper_vs;
									pMapper_vp.coupledCoeffs[0] = CS.Vs_versus_Vp.slope;
									pMapper_vp.coupledConst[0] = CS.Vs_versus_Vp.intercept;
									pMapper_vp.coupledProperties[1] = pMapper_rho;
									pMapper_vp.coupledCoeffs[1] = CS.Rho_or_Phi_versus_Vp.slope;
									pMapper_vp.coupledConst[1] = CS.Rho_or_Phi_versus_Vp.intercept;										
									used_indices.put(I,I);
									used_indices.put(Is,Is);
									used_indices.put(Iden,Iden);
									// add property
									nonTimeProperties.add(pMapper_vp);
								  }
							  }
						  }
					  }
					else
					  {		
						// Individual + independent rock properties in optimiser

						// (1) Impermeable rock p-wave velocities
						// get actual in-situ property prop. It's value is prop.value
						I =  (Integer) map_to.get(LPI.vp_m);
						if (I!=null && used_indices.get(I) == null)
						  {
							used_indices.put(I,I);
							// ordinary ABS mapper
							nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.ABS));
							// logarithmic type mapper prototype, if you wanted to do this kind of thing
							// nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.LOGARITHM,mu.getQuick(I.intValue()),0));
						  }
						// (2) Impermeable rock density
						I =  (Integer) map_to.get(LPI.rho_m);
						if (I!=null && used_indices.get(I)==null)
						  {
							used_indices.put(I,I);
							nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.ABS));
						  }
					  }
				  }
				else
				  {
					// sand dominated case:
					if (!useOldParameters)
					  {
						// There are sometimes problems here when the optimiser hits illegal Gassman states...
						// If serious, set true above to false. This omits sand properties from the optimising space.
						if (searchAlongRockCurves)  // LINKING TURN OFF
						  {
							// Coupled vp and vs Mapper properties in optimiser helps stay away from Gassman boundaries
							// This couples the vp,vs,phi along the maximum prior line in {vp,vs,phi} space

							// (1) Vp and Vs
							I =  (Integer) map_to.get(LPI.vp_s);           // Vp
							if (I != null && used_indices.get(I) == null)
							  {
								Is =  (Integer) map_to.get(LPI.vs_s);      // Vs
								if (Is != null && used_indices.get(Is) == null)
								  {
									Iden =  (Integer) map_to.get(LPI.phi_s);      // phi_s
									if (Iden != null && used_indices.get(Iden) == null)
									  {
										// Ok, vp and vs and phi are available free properties
										// Regression data for this sand...
										RockandFluidData.RegressionSet CS;
										CS  = (RockandFluidData.RegressionSet)_m.RFD.RockTrends.get(_m.MPD.MetaInfo[i].sand_rock_member);
										
										// Add ordinary ABS mappers for vp and vs, but couple vs to vp
										
										// parent
										PropertyMapper pMapper_vp = new PropertyMapper(I.intValue(),PropertyMapper.ABS);
										//children
										PropertyMapper pMapper_vs = new PropertyMapper(Is.intValue(),PropertyMapper.ABS);
										PropertyMapper pMapper_phi = new PropertyMapper(Iden.intValue(),PropertyMapper.ABS);
										// entries
										pMapper_vp.makeCoupled(2);
										pMapper_vp.coupledProperties[0] = pMapper_vs;
										pMapper_vp.coupledCoeffs[0] = CS.Vs_versus_Vp.slope;
										pMapper_vp.coupledConst[0] = CS.Vs_versus_Vp.intercept;
										pMapper_vp.coupledProperties[1] = pMapper_phi;
										pMapper_vp.coupledCoeffs[1] = CS.Rho_or_Phi_versus_Vp.slope;
										pMapper_vp.coupledConst[1] = CS.Rho_or_Phi_versus_Vp.intercept;										
										used_indices.put(I,I);
										used_indices.put(Is,Is);
										used_indices.put(Iden,Iden);
										// add property
										nonTimeProperties.add(pMapper_vp);
									  }
								  }
							  }
						  }
						else
						  {
							// Individual + independent rock properties in optimiser

							// (1) Impermeable rock p-wave velocities
							// get actual in-situ property prop. It's value is prop.value
							I =  (Integer) map_to.get(LPI.vp_s);
							if (I != null && used_indices.get(I) == null)
							  {
								used_indices.put(I,I);
								// ordinary ABS mapper
								nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.ABS));
								// logarithmic type mapper prototype, if you wanted to do this kind of thing
								// nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.LOGARITHM,mu.getQuick(I.intValue()),0));
							  }
							// (2) Permeable rock porosity
							I =  (Integer) map_to.get(LPI.phi_s);
							if (I!=null && used_indices.get(I)==null)
							  {
								used_indices.put(I,I);
								nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.ABS));
							  }
						  }
					  }
				  }
				
				// (3) Net-to-Gross ??
				// Handle this as a logit map acting on [0,1]
				I =  (Integer) map_to.get(LPI.NG);
				if (I!=null && used_indices.get(I)==null)
				  {
					used_indices.put(I,I);
					nonTimeProperties.add(new PropertyMapper(I.intValue(),PropertyMapper.LOGIT,0.0,1.0));
				  }
			  }
			used_indices = null; // help GC
		  }

		// set dimensions right // JG 2013
		timeFinder.setStandardProblemDimensions();


		// create full parameter vector for optimiser
		timesPlusOtherParameters = new double[_m.n_active_layers + nonTimeProperties.size()];

		
		// and full snapshot vector for the current realisation (includes constant properties)
		property_sample = new double[mu.size()];

		// test if pinchout-free solution likely, used for branching directly to Gauss-Newton
		pinchoutFreeSolutionLikely = true;
		if (!ALWAYS_ATTEMPT_GAUSS_NEWTON_FIRST)
		  {
			// Assume pinchoutFreeSolutionLikely is true, so Gauss-Newton will be tried first.
			// Don't test conditions that may set it false.
			for (i=1;i<_m.n_active_layers-1;i++)
			  {
				if (dt_uncert.getQuick(i) > LAYER_PINCHOUT_THRESHOLD)
				  {
					pinchoutFreeSolutionLikely = false;
					break;
				  }
			  }
		  }

		// Here we need to loop over a good set of starting configurations
		// These fall into two groups:

		// (1) Starting points are the times + other properties at the
		// prior+isopach posterior mean obtained by the linearisation update
		// INITIAL layer times guess
		for (i=0;i<_m.n_active_layers;i++) timesPlusOtherParameters[i] = mu.getQuick(i);

		// JG 2013: ensure timesPlusOtherParameters fully populated & consistent
		pickNonTimeStartingValues(0,timesPlusOtherParameters);
		timeFinder.fillParameterVectorWithNonTimeProperties(timesPlusOtherParameters);

		// make first search for each fluid mode free of descent test
		searchForMode(modeIdentifiers,timesPlusOtherParameters,false,null);



		// (2) Set of `doublets' suspected from relative time_thickness errors after prior+isopach
		// inclusion. A 'doublet' is a starting configuration with one particular pinchout enforced:
		// there are clearly a suite of 'doublets' for any model with more than 3 layers. From the
		// prior+isopach uncertainty of the layer times we can compute how probable a pinchout seems
		// in the absence of seismic data, and loop over the probable pinchout configurations to look
		// for solutions.
		// We loop over possible pinchouts in internal layers only
		for (j=1;j<=_m.n_active_layers-2;j++)
		  {
			if (dt_uncert.getQuick(j) > LAYER_PINCHOUT_THRESHOLD)
			  {
				// Here we have a doublet formed by master layers j and j+1
				// reduced-model indexes (as in mu, timesPlusOtherParameters) are j and j+1

				// This layer seems `likely' to pinch out. Seek for a new
				// mode assuming starting conditions that impose the said pinchout
				// Starting times for previous layers
				for (k=1;k<j;k++) timesPlusOtherParameters[k] = mu.getQuick(k);
				// Set the starting times for the bounding layers in the middle of the layer,
				// weighted by the relative uncertainties in each bounding time
				timesPlusOtherParameters[j] = (t_uncert.getQuick(j+1) * mu.getQuick(j) + t_uncert.getQuick(j) *mu.getQuick(j+1))/
				  (t_uncert.getQuick(j) + t_uncert.getQuick(j+1));
				// ensure pinchout will not be crossed by finite difference minimiser method calls
				// (on first attempt anyway)
				timesPlusOtherParameters[j+1] =  timesPlusOtherParameters[j] - 2.0 * SEISMIC_DESIGN_MATRIX_TIME_EPS;
				// following layers
				for (k=j+2;k<=_m.n_active_layers-1;k++) timesPlusOtherParameters[k] = mu.getQuick(k);

				// now get mode
				searchForMode(modeIdentifiers,timesPlusOtherParameters,_SI.useDescentTest,null);
			  }
		  }

		// (3) Event-detected starting points.
		// Also, skip if we ignore seismic (there's no point detected events that aren't in the likelihood).
		if (_SI.useEventDetection && !_SI.ignoreSeismic)
		  {
			// detects events on near data
			if (_m.SM.nearstack_data!=null)
			  {
				EventDetector ED = new EventDetector(_m);

				double[] eventTimes = ED.detectLayerTimes();
				if (eventTimes != null)
				  {
					for (k=0;k<=_m.n_active_layers;k++)
					  {
						timesPlusOtherParameters[k] = eventTimes[k];
					  }
					searchForMode(modeIdentifiers,timesPlusOtherParameters,false,null);
				  }
			  }
		  }


		// (4) Experimental. Random time-starts drawn from prior
		if (false)
		  {
			if (_verbosity>3) _SI.EH.println("Random starts for fluid state "+fs+"...");
			for (j=1;j<=50;j++)
			  {
				// don't thrash top or bottom (no reflections)
				for (k=1;k<=_m.n_active_layers-1;k++) 
				  {
					timesPlusOtherParameters[k] = mu.getQuick(k) + Ngen.nextDouble()*t_uncert.getQuick(k);
				  }
				// now get mode
				searchForMode(modeIdentifiers,timesPlusOtherParameters,_SI.useDescentTest,null);
			  }
		  }

		// (5) Finally, if these local approaches haven't yielded anything plausible, try a Genetic algorithm
		testForAndInitiateGeneticModeFindingIfRequired(localModes, modeIdentifiers, timesPlusOtherParameters);
	  } // loop over fluid states

 	// GC
	t_uncert = null;
	dt_uncert = null;
	timesPlusOtherParameters =null;
	property_sample = null;

	// Final crisis-saving placeholder if all the optimisations fail for stupid models
	// set a warning flag noModesFound: Pollute output files seismic_chisquare with marker if noModesFound=true
	_SI.EH.println("%nFound "+localModes.size()+" mode(s)");
	noModesFound = false;
	if (localModes.size() == 0)
	  {
		noModesFound = true;
		if (_verbosity>1) _SI.EH.println("WARNING! No modes found! Model or data very likely dodgy!. "+
					   "Using prior of most-likely fluid mode as a placeholder!!!");
		drawMostLikelyPrior(false); // sets fs, layer_code, mu = mbar_prior_plus_isopachs[_m.FS.fs_max];
		// All entries below are really token placeholders. No chance of any meaning..just allows code to press on
		ModifiedCholeskyRightDecomposition CDR = new ModifiedCholeskyRightDecomposition(Cinv_prior_only[fs]);
		// overwrite it on R
		DoubleMatrix2D Rinv = CDR.getR();
		// now take inverse (in place)
		CDR.getCholeskyRightInverse(Rinv);
		LocalMode mode = new LocalMode(1.0,1e10,fs,fluid_state,layer_code,mu.copy(), Rinv);
		mode.Rtransp = mode.R.viewDice();
		mode.seismic_chisq = 1e10;
		mode.logDetCprior = logDetC_prior_only[fs];
		mode.prior_weight = 1.0;
		mode.logMargModelLikelihood = 1.0;
		int _n = Cinv_prior_only[fs].size();
		mode.CtildeInvDiag = new double[_n];
        for(k=0;k<_n;k++) mode.CtildeInvDiag[k] = Cinv_prior_only[fs].getQuick(k,k);
		localModes.add(mode);
	  }


	if (_verbosity>2) _SI.EH.println("%nFound "+localModes.size()+" mode(s)");

	// (A) Now sort array of local modes into a useful order
	// First, normalise quasi-posterior weights
	localModes.trimToSize();
	// loop over modes to establish normalisation
	weight_norm = Double.MAX_VALUE;
	for (k=0; k< localModes.size(); k++) 
	  {
		// weight_norm = min((localModes.get(k)).weight, weight_norm);
		weight_norm = min((localModes.get(k)).getPriorFSAdjustedWeight(), weight_norm);
	  }
	// normalise relative 'posterior' probabilities to be less than 1
	// Accumulate also un-normalised marginal model likelihood
	double least_nLL = Double.MAX_VALUE;
	for (k=0; k< localModes.size(); k++)
	  {
		mode = localModes.get(k);
		// normalise weight so most probable state has probability of 1
		// tmp = mode.weight - weight_norm;
		tmp = mode.getPriorFSAdjustedWeight() - weight_norm;
		// retain record of most likely mode in posterior
		least_nLL = min(least_nLL, tmp);
		mode.weight = tmp;
		if (_verbosity>3) _SI.EH.println("Initial mode list entry: "+k+" weight = "+tmp);
	  }
	// drop modes that seem impossible wrt best mode, where
	// the best mode has the negloglikelihood least_nLL
	for (k=kk=0; k<localModes.size(); k++,kk++) 
	  {
		if ((localModes.get(k)).weight > least_nLL + MODE_NLL_ELIMINATION_MARGIN)
		  {
			if (_verbosity>3) _SI.EH.println("Dropping mode "+kk);
			localModes.remove(k--);
		  }
	  }

	// New code (Aug2010) to rebalance prior weights to compensate for
	// unbalanced numbers of local modes per fluid state
	setModeWeights();

	// equalise the number of modes per fluid combo if required
	if (_SI.equaliseModeCountAcrossFluidCombinations)
	  {
		// ensure equal number of modes per fluid combo....
		int[] n_modes = new int[_m.FS.possibleFluidStates.size()];
		Arrays.fill(n_modes,0);
		// count no. of modes per fluid state
		for (k=0; k< localModes.size(); k++) 
		  {
			mode = localModes.get(k);
			n_modes[mode.fs]++;
		  }
		int least_number_of_modes_per_fs = localModes.size();
		for (k=0;k<_m.FS.possibleFluidStates.size();k++)
		  {
			if (least_number_of_modes_per_fs > n_modes[k])
			  {
				least_number_of_modes_per_fs = n_modes[k];
			  }
		  }

		// Now eliminate modes so each fluid state retains only the best least_number_of_modes_per_fs modes
		Collections.sort(localModes, new LocalMode.LocalModeSorter(true,false));

		if (least_number_of_modes_per_fs<1) least_number_of_modes_per_fs=1;
		int fs_old=0,n_this_mode=0;
		for (k=0; k< localModes.size(); k++) 
		  {
			mode = localModes.get(k);
			if (k==0) fs_old = mode.fs;
			if (mode.fs==fs_old)
			  {
				n_this_mode++;
				// strip out excess one if too many
				if (n_this_mode>least_number_of_modes_per_fs)
				  {
					localModes.remove(k);
					k--;
				  }
			  }	
			else
			  { // just jumped to a new fs
				fs_old = mode.fs;
				n_this_mode=1;
			  }
		  }
		localModes.trimToSize();
	  }
	// redo normalisation to counter dropped modes, fluid-state equalisation etc
	reNormaliseLocalModes(true);
	// shan't need this anymore
	modeIdentifiers=null;

	// Show mode profile fence sections if required
	if (checkModeJoiningProfiles(_SI.printModeProfiles))
	  {
		// recompute mode weights from scratch again
		setModeWeights();
		// have to do this
		reNormaliseLocalModes(true);
	  }

	// print cross sections if required
	if (_SI.printXsections) 
	  {
		property_sample = new double[getMaxModeSize()];
		for (k=0; k< localModes.size(); k++)
		  {
			setMode(localModes.get(k));
			debugPrintXsections("PosteriorXsection.fs"+fs+".mode"+k,mode);
			// computeBigJumpDirectionVectorsFromPrincipalComponents(mode);
		  }
		property_sample = null;
	  }
  }

  /**
	 Set mode weights based on current set of local optima in LocalModes
	 Rebalances prior_weights and mode weights depending on number of
	 local modes per fluid state
  */
  void setModeWeights()
  {
	int k,i,j,size;
	double least_nLL,tmp,weight_norm;
	double[] normFactors;
	int[] nModesPerFS;

	size = localModes.size();

	// recompute estimated weights using MLL
	normFactors = new double[_m.FS.possibleFluidStates.size()];
	// first, safe rescaling
	for (k=0, least_nLL = Double.MAX_VALUE; k< size; k++)  
	  { 
		mode = localModes.get(k);
		least_nLL = min(least_nLL, mode.logMargModelLikelihood);
	  }
	// sum of MML probabilities (rescaled for safety) over modes shared by each fluid state
	for (k=0; k< size; k++) 
	  {
		mode = localModes.get(k);
		normFactors[mode.fs] += exp(-(mode.logMargModelLikelihood - least_nLL));
	  }
	// Now reweight mode weights : hopefully underflow/overflow safe here
	for (k=0; k< size; k++) 
	  {
		mode = localModes.get(k);
		// relative contribution (factor 0<tmp<1) of this mode among the group with the same fs
		if (normFactors[mode.fs]>0) 
		  tmp = exp(-(mode.logMargModelLikelihood - least_nLL))/normFactors[mode.fs];
		else tmp=1.0;
		// Decent looking formula for mode probability
		// m_prob = tmp * exp(-mode.logMargModelLikelihood) * mode.FS.prob;
		// -log() of m_prob, above
		mode.weight = -log(tmp) + mode.logMargModelLikelihood - log(mode.FS.prob);
	  }
	// Again, loop over modes to establish normalisation
	weight_norm = Double.MAX_VALUE;
	for (k=0; k< size; k++) 
	  {
		weight_norm = min(localModes.get(k).weight, weight_norm);
	  }
	for (k=0; k< size; k++)
	  {
		mode = localModes.get(k);
		// normalise weight so most probable state has probability of 1
		mode.weight -= weight_norm;
	  }

	// Now reallocate prior weights
	nModesPerFS = new int[_m.FS.possibleFluidStates.size()];
	for (k=0; k< size; k++) 
	  {
		mode = localModes.get(k);
		nModesPerFS[mode.fs]++;
	  }
	if (false)
	  {
		// reallocate prior weight from fluid-state probability 
		// equally across modes found in this fluid state. Guard against /0
		for (k=0; k< size; k++) 
		  {
			mode = localModes.get(k);
			mode.prior_weight = mode.FS.prob/max(1,nModesPerFS[mode.fs]);
		  }
	  }
	else
	  {
		// partition-up mode.FS.prob into portions proportional to the
		// MML weight, for efficiency
		for (i=0; i< size; i++) 
		  {
			mode = localModes.get(i);
			for (j=i+1;j<size;j++)
			  {
				if ((localModes.get(j)).fs != mode.fs) break;
			  }
			// same-fluid-state  modes span the range [i,j-1] now
			// total prob of members
			for(k=i,tmp=0.0; k<=(j-1); k++) tmp += exp(-localModes.get(k).weight);
			// farmed out
			for(k=i; k<=(j-1); k++)
			  {
				mode = localModes.get(k);
				mode.prior_weight = mode.FS.prob*exp(-mode.weight)/tmp;
			  }
			i=j-1;// to start at new chunk above
		  }
	  }

	normFactors=null; // GC
	nModesPerFS=null;
  }

  void dumpMAPmodels()
  {
	int i,k,l,size,kMax=0,n;
	double wt,optimum=0; 
	LocalMode lM;

	// Set to false to get an su file with a dump of all modes, keyed by mark=0....size-1
	boolean dumpBestModeOnly = true;

	property_sample = new double[getMaxModeSize()];

	// find most likely mode
	if (_SI.drawMinChisqModels) optimum=Double.MAX_VALUE;
	size = localModes.size();
	for (k=0; k< size; k++)
	  {
		lM = localModes.get(k);
		if (_SI.drawMinChisqModels)
		  {
			// lowest chisq
			wt = lM.chisq;
			if (wt < optimum) 
			  {
				kMax = k;
				optimum = wt;
			  }
		  }
		else
		  {
			// highest weight (MML)
			wt = lM.weight;
			if (wt >  optimum) 
			  {
				kMax = k;
				optimum = wt;
			  }
		  }
	  }
	if (_verbosity>3) _SI.EH.println("Dumping "+(_SI.drawMinChisqModels ? "min-chisq" : "MAP")+" point for mode "+kMax);

	for (k=0,l=0; l < (dumpBestModeOnly ? 1 : (size+1)); k++,l++)
	  {
		if (dumpBestModeOnly || (l==size)) k = kMax;
		else k=l;

		// set mode to most likely mode
		setMode(localModes.get(k));
		// central point 
		DoubleMatrix1D centre_point = mode.QP_optimum;
	
		n = centre_point.size();
		for (i=0;i<n;i++)
		  {
			double val = centre_point.getQuick(i);
			setValuesFromIndex(i,val);
			property_sample[i]=val; // needed for prior
		  }
		orderTimes();
		_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false);

		// Must compute seismic to enforce Refl coefficients etc at correct places
		debugComputeNegLogPosterior();

		// dump model. Last pass is to restore state
		if (l<size) _SI.dumpRealisation(true,l,true,true);
	  }
	// GC
	property_sample = null;
  }


  // Return length of longest mode vector in list of modes localModes
  int getMaxModeSize()
  {
	int max_mode_size=0;
	LocalMode mode;
	for (int k=0;k<localModes.size();k++)
	  {
		mode = localModes.get(k);
		max_mode_size = (k==0) ? mode.QP_optimum.size() : 
		  max(mode.QP_optimum.size(), max_mode_size);
	  }
	return max_mode_size;
  }

  /**
	 Method to perform normalisation
	 N.B. if convertFromLogWeights==false, 
	ensure mode.weights are *actual* weights (not -log(weights))
  */
  void  reNormaliseLocalModes(boolean convertFromLogWeights)
  {
	int k,j,size;
	double norm,prior_norm;
	String fstr = ""; 

	size = localModes.size();


	for (k=0,norm=prior_norm=0.0; k< size; k++) 
	  {
		mode = localModes.get(k);
		norm += ((convertFromLogWeights) ? exp(-mode.weight) : mode.weight);
		prior_norm += mode.prior_weight;
	  }

	// convert weight probabilities to real ones now
	double sum=0.0,prior_sum=0.0;
	if (_verbosity>2) _SI.EH.println("%nNormalised local modes summary:");
	for (k=0; k< size; k++) 
	  {
		mode = localModes.get(k);
		// normalise weight so sum of states is unity
		// mode.weight's are henceforth absolute +ve probabilities
		mode.weight = ((convertFromLogWeights) ? exp(-mode.weight)/norm : mode.weight/norm);
		mode.prior_weight /= prior_norm;
		sum += mode.weight;
 		prior_sum += mode.prior_weight;
 		// compute both CDF's
 		mode.CDF = sum;
 		mode.prior_CDF = prior_sum;
		// and the half-weighted hybrid weight
		mode.hybrid_weight = (mode.prior_weight + mode.weight)/2.0;
		if (_verbosity>2) 
		  {
			if (_SI.drawMAPmodels || _SI.drawMinChisqModels)
			  {
				FluidStates.FluidState f = _m.FS.possibleFluidStates.get(mode.fs);
				fstr = f.fState;
			  }
			_SI.EH.println("Normalised mode entry: mode number %2d ,fs= %2d (%s) ,dimension= %3d ,layer_code= %3d "
							  +" ,weight= %-10.6g"
							  +" ,CDF= %-10.6g ,-log(MML)= %-14.6g ,-chisq(raw)= %-14.6g, prior_weight= %-10.6g"
						   +" ,prior_CDF= %-10.6g",k, mode.fs, fstr,
						   mode.QP_optimum.size(), mode.layer_code,
						   mode.weight, mode.CDF,
						   mode.logMargModelLikelihood, mode.chisq,
						   mode.prior_weight, mode.prior_CDF);
		  }		
	  }

	_SI.EH.println("");
	// set last mode CDF to be 1.0 exactly (precaution against roundoff)
	mode.CDF = mode.prior_CDF = 1.0;


	if (_verbosity>2)
	  {
		// Add also preliminary guess of fluid probabilities, on a layer basis....
		int i;
		double prob;
		boolean[] fsIsInteresting = {false,false,false,false};
		for (i=0;i<_m.n_layers;i++)
		  {
			for (k=0; k< size; k++) 
			  {
				mode = localModes.get(k);
				FluidStates.FluidState fs = (FluidStates.FluidState)_m.FS.possibleFluidStates.get(mode.fs);
				if (FluidStates.FluidTypes.getFluidChar(fs.types[i]) != FluidStates.FluidTypes.brine)
				  {
					fsIsInteresting[fs.types[i]] = true;
				  }
			  }
			for (j=1;j<FluidStates.FluidTypes.Nfluids;j++)  // examine oil,lsg,gas
			  {
				if (fsIsInteresting[j])
				  {
					// sum probabilities of modes with this matching fluid type in this layer
					prob=0;
					for (k=0; k< size; k++) 
					  {
						mode = localModes.get(k);
						FluidStates.FluidState fs = (FluidStates.FluidState)_m.FS.possibleFluidStates.get(mode.fs);
						if (fs.types[i] == j) prob += mode.weight;
					  }
					_SI.EH.println("MAP-estimated fluid probability: layer %3d type %s ,pr= %-10.6g",
								   // (i+1),"HI",prob);//FluidStates.FluidTypes.getFluidChar(fs.types[i]),prob);
								   (i+1),FluidStates.FluidTypes.getFluidChar((short)j),prob);
					// clear type
					fsIsInteresting[j]=false;
				  }
			  }
			// clear types
			for(j=0;j<4;j++) fsIsInteresting[j]=false;
		  }
		fsIsInteresting=null;
	  }
  }	


  private void testForAndInitiateGeneticModeFindingIfRequired(ArrayList<LocalMode> localModes, 
															  ArrayList modeIdentifiers, 
															  double[] timesPlusOtherParameters)
  {
	int i,k;
	double best_seismic_chisq=Double.MAX_VALUE;
	double xi,xi_max=0.0,lambda;
	DoubleMatrix1D t_uncert;
	boolean multimodalityLikelyFromWeakPrior=false;

	// loop over current modes with this fluid state to see if genetic mode-finding needs to be
	// fired off
	for (k=0;k<localModes.size();k++)
	  {
		mode = localModes.get(k); 
		// only interested in matches for current fluid state
		if (mode.fs!=fs) continue;
		best_seismic_chisq = min(mode.seismic_chisq,best_seismic_chisq);
		if (_verbosity>3) _SI.EH.println("Best_seismic_chisq = %g%n",best_seismic_chisq);
	  }


	// Examine also the standard deviation of the prior+isopach times relative to the scale of the
	// wavelet correlation length. This is a measure of how far the mode search is allowed to stray
	// from the prior-mean-times in order to find a mode, and how many "loops" it might have traversed
	// in doing so, and therefore the likelihood of there being considerable multimodality in the problem.
	// I.e. if the prior is quite tight, we won't have gone very far, and therefor won't be likely
	// to be loop-skipping. We seek an indication of the converse.
	// For a single reflection, the secondary minima is at dt = 0.6 lambda
	t_uncert = t_layer_uncertainty.viewRow(fs);  // array of standard deviations of times from prior +isopach
	double fBE = _m.SM.getWaveletBandEdgeFrequency();
	for (i=1;i<_m.n_active_layers;i++)
	  {
		// scaled measure assuming 4:1 S/N ratio
		xi = 14.0*MathUtil.sqr(t_uncert.getQuick(i)*fBE);
		xi_max = max(xi_max,xi);
		if (_verbosity>3) _SI.EH.printf("GA: layer "+i+" xi-factor = %g%n",xi);
	  }
	// here's the magic criterion for guessing multimodality from the prior's alone.
	// If we ignore the seismic, the surface will be unimodal, so take this into account
	if (xi_max > XI_THRESHOLD_FOR_MULTIMODALITY && !_SI.ignoreSeismic)
	  {
		multimodalityLikelyFromWeakPrior = true;
	  }
	t_uncert=null;


	/*
	  Ok, if the best match we've achieved so far is still hopeless,
	 maybe we got caught in a dodgy local minima. Or maybe the prior
	 is so weak there could yet be a better minima elsewhere.  Fire of
	 the global optimiser, as the prior wasn't good enough to point us
	 at a good starting place using the efficient local minimiser
	*/
	if ((best_seismic_chisq> UNACCEPTABLE_SCALED_SEISMIC_MISMATCH_GLOBAL_OPTIMISER_THRESHOLD
		 || multimodalityLikelyFromWeakPrior) && _SI.useGeneticAlgorithmsIfNeeded)
	  {
		if (_verbosity>2) 
		  {
			_SI.EH.println("Suspect strong multimodality: Starting Genetic solver...");
		  }
		gTimeFinder = new GeneticOptimalTimesFinder(_m);

		// solve for the best global times using other properties at their modes
		gTimeFinder.findOptimalSolution();

		// Now start local mode hunting from this solution. The best times[] solution will have been left in
		// the first entries of property_sample[] by the Evolutionary optimiser
		for (k=1;k<_m.n_active_layers;k++) timesPlusOtherParameters[k] = property_sample[k];
		searchForMode(modeIdentifiers,timesPlusOtherParameters,true,null);
	  }
  }


  private void searchForMode(ArrayList modeIdentifiers, 
							 double[] timesPlusOtherParameters, 
							 boolean performDescentTest,
							 LocalMode goodStartingMode)
  {
	int i,k,n,i_trial;
	double f=0.0;
	LocalMode new_mode=null;
	double[] timesPlusOtherParameters_cpy;
	boolean hasFaultyNewtonConvergence=false;
	ModeIdentifier mi=null, best_mi=null;
	long runTime,updateTime=0,tmpTime;
	boolean branchToGaussNewton = false;

	runTime = System.currentTimeMillis();

	timesPlusOtherParameters_cpy = (double[]) timesPlusOtherParameters.clone();
	n = mu.size();	

	// loop over minimisation attempts, with different initial conditions
	// attempt (i_trial=0):        everything set at the modal values as starting point, except passed-in times
	// attempt (i_trial=1,2,3...): random guesses from the prior conditioned upon the input times
	// We bail out at the first attempt that doesn't fail.
	modeSearchSuccessFlag=0;   // Flag to indicate success until overwritten with failure indication

	for (i_trial=0;i_trial<MAX_NUM_TRIAL_GUESSES;i_trial++)
	  {
		if (i_trial>0 && modeSearchSuccessFlag==1) 
		  {
			// As of 27/5/02, shouldn't need random guesses, so itrial==0 nearly always succeeds.
			// If not, there's probably a bug, so raise a warning
			_SI.EH.println("WARNING: mode-means are not Gassman-safe!!! (Probably a bug)");
		  }
		// later on in the loop, timesPlusOtherParameters[] gets polluted, so we have to restore it here
		System.arraycopy(timesPlusOtherParameters_cpy,0,timesPlusOtherParameters,0,timesPlusOtherParameters.length);
		
		// Set initial guesses for the optimiser.
		// First guess is at mode mu, and is guaranteed 'Gassman safe' 
		// for fluid state 0, for i_trial==0
		if (!setOptimiserInitialValues(i_trial,timesPlusOtherParameters,goodStartingMode)) 
		  {
			if (goodStartingMode == null) 
			  {
				continue;  // usual mode of operation
			  }
			else
			  { 
				return; // couldn't set starting conditions sensibly from suggested mode. Give up
			  }
		  }


		// For subsequent attempts, if starting point is obviously purely "uphill" from best mode
		// found so far, we save time by skipping the optimisation: A descent method seems 
		// very likely to end up where we were before,
		// so there's no point in starting the optimiser...

		if (performDescentTest && timeFinder.testForPureDescentPath(timesPlusOtherParameters,modeIdentifiers))
		  {
			if (_verbosity>2) _SI.EH.println("Skipping mode search by descent test");
			// It's a pure descent to the known mode: don't bother pressing on
			return;
		  }

		// ---------- Branch of BFGS vs Gauss-Newton methods here! -----------------
		branchToGaussNewton = (pinchoutFreeSolutionLikely && (i_trial==0));
		// branchToGaussNewton=false;
		if (_SI.doSimpleFluidUnitsAndContacts)
		  {
			// Avoid getting into first branch below if at all possible
			branchToGaussNewton = true;
		  }
		if (!branchToGaussNewton)
		  {
			// optimise layer-boundary times using BFGS; expect possible pinchouts
			if (_SI.doSimpleFluidUnitsAndContacts)
			  {
				// FIXME
				_SI.EH.printf("FIXME: findOptimalSolution() methods not yet working for contacts problems! Try switching contacts off!%n");
				System.exit(1);
			  }
			f = timeFinder.findOptimalSolution(timesPlusOtherParameters); 
			// Call above will return with timesPlusOtherParameters[] at optimum best values
			if (_verbosity>3)
			  {
				_SI.EH.printf("Best match at overall f = %10.8g, prior = %10.8g, isopach = %10.8g, "+
								"depth-self-cons = %10.8g, vert-stack = %10.8g, seismic NLL =  %10.8g%n",
								f, timeFinder.getNegLogPriorLikelihood(),
								timeFinder.getNegLogIsopachLikelihood(),
								timeFinder.getNegLogDepthLikelihood(),
							    timeFinder.getNegLogVStackLikelihood(),
								timeFinder.getNegLogSeisLikelihood());
				_SI.EH.print("Best match times: ");
				for (k=1;k<=_m.n_active_layers-1;k++) _SI.EH.printf("%9.6g ",timesPlusOtherParameters[k]);
				_SI.EH.print("%n");
			  } 
		  }
		else
		  {
			// Pinchout-free solution: try pure Gauss-Newton on first attempt
			if (_verbosity>3)
			  {
				_SI.EH.print("Branching to Gauss-Newton; best-match times are start of Gauss-Newton updates...%n");
			  }
			System.arraycopy(timesPlusOtherParameters,0,timesPlusOtherParameters_cpy,0,timesPlusOtherParameters.length);
			// remap times to true times if needed
			timeFinder.mapSubTimesVectorToCompleteTimes(timesPlusOtherParameters,true);
			f = timeFinder.f_to_minimize(timesPlusOtherParameters);
			// Call above will return with timesPlusOtherParameters[] at starting values
			if (_verbosity>3)
			  {
				_SI.EH.printf("Starting Gauss-Newton at overall f = %10.8g, prior = %10.8g, isopach = %10.8g, "+
								"depth-self-cons = %10.8g, vert-stack = %10.8g, seismic NLL =  %10.8g%n",
								f, timeFinder.getNegLogPriorLikelihood(),
								timeFinder.getNegLogIsopachLikelihood(),
								timeFinder.getNegLogDepthLikelihood(),
							    timeFinder.getNegLogVStackLikelihood(),
								timeFinder.getNegLogSeisLikelihood());
				_SI.EH.print("Active layer times: ");
				for (k=1;k<=_m.n_active_layers-1;k++) _SI.EH.printf("%9.6g ",timesPlusOtherParameters[k]);
				_SI.EH.print("%n");
			  }
		  }


		
		// if the optimiser ended up with some times at the region boundary, reject this mode as unsuitable
		if (timeFinder.testAndImposeTimeTruncation(timesPlusOtherParameters,true))
		  {
			if (false)
			  {
				if (_verbosity>3) // FIXME do we really want this? what do we really want??????
				  // i.e. we may enforce a layer absence in the priors that looks illegal here, but in fact is ok.....
				  {
					_SI.EH.println("Mode *rejected* due to time truncations");
					return;
				  }
			  }
		  }

		// remap times to true times if needed
		timeFinder.mapSubTimesVectorToCompleteTimes(timesPlusOtherParameters,true);

		// Scan over existing modes to test if we have found this mode already. 
		// Duplicates are frequently found by starting at different points.
		// We scan over the array of existing mode - identifiers
		// ( do this before wasting time on the getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods() call)
		for (i=0;i<modeIdentifiers.size();i++)
		  {
			// assume it matches until we refute this
			boolean matches=true;
			mi = (ModeIdentifier) modeIdentifiers.get(i);
			if (mi.fs != fs) 
			  {
				matches=false; continue;
			  }
			if (mi.layer_code != layer_code)
			  {
				matches=false; continue;
			  }
			if (_verbosity>3) 
			  {
				_SI.EH.printf("Scaled mode-time comparisons for history list: mode "+i+", ");
			  }
			for (k=1;k<=_m.n_active_layers-1;k++) 
			  {
				if (_verbosity>3)
				  {
					_SI.EH.printf("%9.6g ",Math.abs(mi.timesPlusOtherParameters[k] 
													  - timesPlusOtherParameters[k])
									/TIME_EPSILON_FOR_MINIMA_MATCH);
				  }
				if (Math.abs(mi.timesPlusOtherParameters[k] - timesPlusOtherParameters[k]) 
					> TIME_EPSILON_FOR_MINIMA_MATCH) 
				  {
					matches = false;
					break;
				  }
			  }
			if (_verbosity>3) _SI.EH.println("");

			// if we survived all these tests, I'm sure the mode is the same
			if (matches)
			  {
				if (_verbosity>3)
				  {
					_SI.EH.println("Mode *rejected* since an identical one already found");
				  }
				return;
			  }
		  }
		// The mode is  unique. Better record it
		if (_verbosity>3)
		  {
			_SI.EH.println("Mode *accepted*");
		  }


		mi = new ModeIdentifier(fs,layer_code,timesPlusOtherParameters,f);	
		modeIdentifiers.add(mi);

		

		// GAUSS NEWTON UPDATES

		// Do an inverse theory mode update based on a linearisation of some
		// subset of the full suite of properties
		// about the point with layer times initially set to what the findOptimalSolution() optimisation
		// yielded, and other properties at their prior+isopach posterior mean		
		tmpTime = System.currentTimeMillis();
		new_mode = getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods(mi);
		updateTime += System.currentTimeMillis()-tmpTime;

		if (new_mode != null)
		  {
			// search seemed to work and not produce an illegal state. 
			// Save link to full mode spec in mode identifier
			mi.mode = new_mode;

			// Clean up and bail out
			break;
		  }
		else
		  {
			// otherwise, the loop will try again, from different starting points...
			if (_verbosity>3)
			  {
				switch (modeSearchSuccessFlag)
				  {
				  case 1:
					_SI.EH.println("Extrapolated mode from seismic update forces illegal Gassman state. "
										+"Trying again...");
					break;
				  case 2:
					_SI.EH.println("Newton method did not converge. Trying again...");
					break;
				  }
			  }
		  }
	  }
	// if, after multiple attempts, we still dip out, give it away. Later on we'll detect if 
	// no decent modes were found at all
	if (new_mode==null) return;

	tmpTime = System.currentTimeMillis();
	for (i_trial=0; i_trial<MAX_NUM_FAULTY_NEWTON_REFINEMENTS; i_trial++)
	  {
		hasFaultyNewtonConvergence = modeHasFaultyNewtonConvergence(new_mode);
		if (hasFaultyNewtonConvergence)
		  {
			if (_verbosity>3) 
			  {
				_SI.EH.println("Mode doesn't look concave-up: retrying Newton method..");
			  }
			// mode doesn't look suitably concave up. Try again from point detected
			new_mode = getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods(mi);
			// if new attempt at Newton gets lost, give up
			if (new_mode==null) 
			  {
				modeSearchSuccessFlag=2;
				return;
			  }
		  }
		else break; // get out first time things look ok
	  }
	if (hasFaultyNewtonConvergence) // might as well give up
	  {
		if (_verbosity>3) _SI.EH.println("Unable to escape from faulty Newton convergence. Giving up.");
		modeSearchSuccessFlag=2;
		new_mode=null;
		return;
	  }

	// if, after multiple attempts, we still dip out, give it away. Later on we'll detect if 
	// no decent modes were found at all
	if (new_mode==null) return;
	
	// Now add this local approximation to the list of local modes
	// Save the mode's negative log(likelihood) at this point, and we 
	// do proper normalisation later
	if (isUpdatedModeUnique(new_mode))
	  {
		localModes.add(new_mode);
	  }

	// finally, dump synthetic seismic if required. Label it with the mode number 
	// which is currently modeIdentifiers.size()
	if (_SI.dumpOptimalSeismic) 
	  {
		if (_verbosity>4)
		  {
			_SI.EH.printf("Dumping optimal seismic, header-labelled %d%n",modeIdentifiers.size());
		  }
		_m.SM.dumpSynthetics(true,modeIdentifiers.size(),true);
		// dump model too 
		if (_SI.dumpRealisations) _SI.dumpRealisation(true,1,true,false);
	  }



	updateTime += System.currentTimeMillis()-tmpTime;

	// Check mode Xsections if desired
	// not here: this is now done after the mode list is massaged + filtered etc
	// if (_SI.printXsections) debugPrintXsections("posteriorXsection.hunt.fs"+fs+".mode"+localModes.size(),new_mode);
	  

	if (_verbosity>3)
	  {
		_SI.EH.println("Added mode of relative weight "+Math.exp(-new_mode.weight));
		_SI.EH.print("%n");
		_SI.EH.print("Mode (QP) times: ");
		for (k=0;k<=_m.n_active_layers;k++)
		  {
			_SI.EH.printf("%10.7g ",new_mode.QP_optimum.getQuick(k));
		  }
		_SI.EH.print("%n");
		_SI.EH.print("Mode (QP) non-time parameters: ");
		for (;k<new_mode.QP_optimum.size();k++)
		  {
			_SI.EH.printf("%10.7g ",new_mode.QP_optimum.getQuick(k));
		  }
		_SI.EH.print("%n");
		_SI.EH.println("Update-time/total mode-hunt time = %6.3g ms",
						(double)updateTime/(System.currentTimeMillis()-runTime));
		// show also isopach mismatch
		if (_m.MPD.n_isopachs>0 && !_SI.ignoreIsopachs)
		  {
			_SI.EH.print("Isopach mismatches at mode with seismic: ");
			for (k=0;k<_m.MPD.n_isopachs;k++)
			  {
				_SI.EH.printf("%10.4g ",isopach_errors[k]);
			  }
			_SI.EH.print("%n");
		  }
	  }
  }

  /**
	 check if mode seems to be the same as one we already have in the collection
   */
  private boolean isUpdatedModeUnique(LocalMode new_mode)
  {
	int i,j,k;
	double scale;
	boolean matches=false;
	DoubleMatrix2D L;

	// Scan over existing modes to test if we can find a "similar" one
	for (i=0;i<localModes.size();i++)
	  {
		// assume it matches until we refute this
		matches=true;
		mode = localModes.get(i);
		if (mode.fs != new_mode.fs) 
		  {
			matches=false;
			continue; // won't match this one
		  }
		if (mode.layer_code != new_mode.layer_code)
		  {
			matches=false;
			continue; // won't match this one
		  }
		// Ok if we got here, we're comparing modes with the same layer config and fluids
		// Now check if optimum seems to be within a standard deviation or so of comparable mode.
		// Since they have the same fluid configuration, the parameters in the model have the
		// same ordering
		L = new_mode.R;
		for (j=0; j<mode.QP_optimum.size(); j++)
		  {
			// measure of scale is the diagonal entry in the covariance
			// We get this by a sum over the row
			for(k=0,scale=0.0;k<=j;k++) scale+=L.getQuick(j,k);
			scale = MODE_COMPARISON_EPSILON*Math.sqrt(scale);
			//_SI.EH.println("SCALE TEST: "+j+" "+(Math.abs(mode.QP_optimum.getQuick(j) - new_mode.QP_optimum.getQuick(j))/scale)); 
			if (Math.abs(mode.QP_optimum.getQuick(j) - new_mode.QP_optimum.getQuick(j)) > scale) 
			  {
				matches=false;
				break;
			  }
		  }
		if (matches) break; // one match is good enough
	  }
	if (matches)
	  {
		if (_verbosity>3)
		  {
			_SI.EH.println("Updated mode *rejected* since an identical one already found");
		  }
	  }
	return !matches;
  }

  /*
	 Set initial guesses of parameters in optimiser that are to be
	 varied. Returns false if initial parameters render Gassman
	 substitution unphysical. For i_trial==0 this should be impossible
	 if the prior is 'Gassman-safe' (all subsequent linearisation updates
	 of the prior have been guaranteed to be Gassman-safe).
  */
  private boolean setOptimiserInitialValues(int i_trial, double[] timesPlusOtherParameters, LocalMode goodStartingMode)
  {
	int i,k;
	
	// Set initial layer times guess as per passed in parameters in timesPlusOtherParameters
	for (k=1;k<=_m.n_active_layers-1;k++) 
	  {
		setValuesFromIndex(k,timesPlusOtherParameters[k]);
		property_sample[k]=timesPlusOtherParameters[k];
	  }
	setValuesFromIndex(0,mu.getQuick(0));                     // top layer time always fixed
	property_sample[0]=mu.getQuick(0);
	setValuesFromIndex(_m.n_active_layers,mu.getQuick(_m.n_active_layers)); // base layer time always fixed
	property_sample[_m.n_active_layers] = mu.getQuick(_m.n_active_layers);

	// Ensure the times are ordered
	orderTimes();
	
	// First attempt, on first fluid combo, search from the prior mean
	// (should always succeed for i_trial==0)

	if (!pickNonTimeStartingValues(i_trial,timesPlusOtherParameters)) return false;

	// First attempt; possibly override with smart initial guess pulled from previous modes on fluid state 0
	if (i_trial==0)
	  {
		if ((fs>0 && _SI.useSmartModeFindingStartingPoints && !pinchoutFreeSolutionLikely) // First fluid-cases sweep
			|| ((goodStartingMode != null) && pinchoutFreeSolutionLikely))                    // 2nd sweep if used
		  {
			// Operate only:
			// (1) for fluid modes >0
			// (2) if useSmartModeFindingStartingPoints mode applies
			// (3) using BFGS (Gauss-Newton applies for tighter isopachs - previous solutions may not be useful
			// OR:
			// Second sweep through fluid states: used best mode found so far (will be Gauss-Newton).
			overrideStartingValuesFromPreviousModes(timesPlusOtherParameters, goodStartingMode); 
		  }
		if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false)) // FC_OK
		  {
			// if override doesn't work, go back to 'Gassman-safe' prior
			pickNonTimeStartingValues(0,timesPlusOtherParameters);
		  }
	  }

	// Give some feedback
	if (_verbosity>3)
	  {
		_SI.EH.printf("%nFluid state fs=%d, starting values ",fs);
			
		for (k=1;k<=_m.n_active_layers-1;k++) _SI.EH.printf("%11.8g ",timesPlusOtherParameters[k]);
		// add feedback for non-time properties
		if (timeFinder.nonTimeProperties!=null && timeFinder.nonTimeProperties.size()>0)
		  {  
			_SI.EH.print(" | ");
			PropertyMapper PM;
			for (k=0; k< timeFinder.nonTimeProperties.size(); k++)
			  {
				// get the Mapper for this property in the parameters
				PM = (PropertyMapper) timeFinder.nonTimeProperties.get(k);
				// print the master parameter even for coupled Mappers 
				// if (PM.isACoupledMapper)
				{
				  // and the arraylist this points to in the maps
				  ArrayList list = (ArrayList) map_from.get(PM.propertyMasterIndx); 
				  // and the first property in the arraylist will give us the value
				  LayerProperties.Property prop = (LayerProperties.Property)list.get(0);
				  _SI.EH.printf("%11.8g ",prop.value);
				}
			  }
		  }
		_SI.EH.print("%n");
	  }

	return true;
  }

  /** 
	  Function to pick starting properties for minimiser, conditional to the given times
	  for the i_trial-th iteration. On first attempt, everything will be set to the mean, subsequently
	  we will try better guesses based on proximity to the seismic.

	  (As of rewrites 27/5/2002, code should never get  into (i_trial>0) branch)
  */
  private boolean pickNonTimeStartingValues(int i_trial, double[] timesPlusOtherParameters)
  {
	int i,j,k,n,N,n_trial;
	DoubleMatrix2D R;
	double s,dt,negLogLik_iso, negLogLik_seis, negLogLik_prior;
	double negLogLik_depth = 0.0, negLogLik=0.0, negLogLik_best=0.0;
	boolean sample_succeeded;
	

	if (i_trial==0)
	  {
		// On first attempt, set remaining variable properties at their 
		// modes...(i.e. the post-isopach update mode, 
		// mu = mbar_prior_plus_isopachs if isopachs are required, otherwise the prior mu)
		// The various getBestLegalGassmanPointAlongCurrentLine() calls involved in the
		// establishment of the mu vector should guarantee this configuration
		// is 'Gassman-safe'
		for (k=_m.n_active_layers+1;k<mu.size();k++)
		  {
			property_sample[k]=mu.getQuick(k);
			setValuesFromIndex(k,property_sample[k]);
		  }
	  }
	else
	  {
		// draw a large number of (say 50ish) samples and keep the best state

		// basic MG system after isopach plus depth self-consistency updating (if any)
		mu = mbar_prior_plus_isopachs[fs];
		R = R_prior_plus_isopachs[fs];

		N = mu.size();         // dimension of full property suite
		n = _m.n_active_layers;       // dimension of times subvector
		double[] best_property_sample = new double[N]; // will keep track of best state

		// sample from doubly-updated Covariance (isopachs plus local times)
		double[] w = new double[N];
		int fail_count=0;
		// Variable 0 is t_top, variable n is t_base
		// solve for random deviates w[i] that would lead to the current starting times
		for (i=0; i<n; i++) 
		  {
			property_sample[i] = timesPlusOtherParameters[i];
			for (s = timesPlusOtherParameters[i]-mu.getQuick(i),j=0;j<=i-1;j++) s -= R.getQuick(j,i) *  w[j];
			w[i] = s/R.getQuick(i,i);
		  }
		// base time has been set at mean, but the times vector does not have this
		property_sample[n] = mu.getQuick(n);
		for (s=0.0,j=0;j<=n-1;j++) s -= R.getQuick(j,n) *  w[j];
		w[n] = s/R.getQuick(n,n);		

		if (_SI.dumpGassmanDiagnostics) showGassmanDiagnostics(true,true);

		// try a range of other values
		for (n_trial=0; n_trial<N_RANDOM_STARTING_TRIALS; n_trial++)
		  {
			do
			  {
				// make a sample
				// Set remaining (non-time) properties
				for (i=n+1;i<N;i++) 
				  {
					w[i] = Ngen.nextDouble();
				  }
				for (i=n+1;i<N;i++)
				  {
					s = 0;
					// From inverse of right cholesky decomposition,
					// Ltransp has true left cholesky factor, but transposed.
					for (j=n+1;j<=i;j++)
					  {
						// notice transpose!
						s += R.getQuick(j,i) *  w[j];
					  }
					// add back the mean
					s += mu.getQuick(i);
					property_sample[i]=s; // record it
					// s now has the random deviate.
					// set the property for this index
					setValuesFromIndex(i,s);
				  }
				// now update full layer properties
				sample_succeeded = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
				if (_SI.dumpGassmanDiagnostics) showGassmanDiagnostics(sample_succeeded,false);
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
				if (fail_count>MAX_MC_FAIL_COUNT)
				  {
					if(useGassmanPenaltyTrick)
					  {
						// FIXME replace the following with something clever  .... :-)
						throw new RuntimeException("Impossible fluid or rock properties in fluid state " +
											   (fs - 1) + ". Regression curves, isopachs, seismic, or "
													   + "fluid properties are probably suspect");
					  }
					else
					  {
						throw new RuntimeException("Impossible fluid or rock properties in fluid state " +
											   (fs - 1) + ". Regression curves, isopachs, seismic, or "
													   + "fluid properties are probably suspect");
					  }
				  }
			  } while (!sample_succeeded);

			// OK, at this point, the sample is valid. See how good it is

		    // (1a) SEISMIC mismatch term
			negLogLik_seis = negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
			negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
			negLogLik = negLogLik_seis;

		    // (1b) Add in logLikelihood contribution from ISOPACHS if they're required
			negLogLik_iso = computeIsopachNegLogLikelihood(false);
			negLogLik += negLogLik_iso;

			// (1c) Add in logLikelihood contribution from DEPTH self-consistency if required
			negLogLik_depth = computeConsistentDepthNegLogLikelihood(false);
			negLogLik += negLogLik_depth;

			// (1d) Add in logLikelihood contribution from Vert stacking
			if (_SI.doStackedModel)
			  {
				double negLogLik_stack= _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
				negLogLik += negLogLik_stack;
			  }

		    // (1e) Also the contribution from the raw PRIOR
		    // (i.e. before isopachs were used in the linearised update)
		    // Acquire logLikelihood contribution from prior
		    // The logLikelihood here must span all the variables
			DoubleMatrix2D Ltrans = R_prior_only[fs]; 
			DoubleMatrix1D mu_prior_only =  mbar_prior_only[fs];

			negLogLik_prior = 0.0;
			for (i=0;i<N;i++) 
			  {
				s=0.0;
				for (j=0;j<=i;j++)
				  {
					dt = (property_sample[j] - mu_prior_only.getQuick(j));
					// remember L is stored in transposed form
					s += Ltrans.getQuick(j,i) * dt;
				  }
				negLogLik_prior += s * s;
			  }
			negLogLik_prior/=2.0;
			negLogLik += negLogLik_prior;

			if (_verbosity>4) _SI.EH.println("Initial guess : nLL="+negLogLik);
			// Keep a record of the best state, and keep this up to date
			if (n_trial==0 || negLogLik < negLogLik_best)
			  {
				// record this.
				negLogLik_best = negLogLik;
				if (_verbosity>3) _SI.EH.println("  Ntrial="+n_trial+
														  ", storing improved guess for : nLL="+negLogLik);
				System.arraycopy(property_sample,0,best_property_sample,0,N);
			  }
		  }
		// at the end of this loop, best_property_sample has what we want.
		// Restore it
		for (i=n+1;i<N;i++)
		  {
			setValuesFromIndex(i,best_property_sample[i]);
			property_sample[i] = best_property_sample[i];
		  }
		best_property_sample=null;
	  }
	

	// now update full layer properties
	// The various getBestLegalGassmanPointAlongCurrentLine() calls involved in the
	// establishment of the mu vector shouyld guarantee this configuration
	// is 'Gassman-safe' for i_trial==0
	return _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
  }

  /*
	If we are optimising on fluid states fs>0, use the solutions found
	for fs=0 (or suggestedMode if non-null) as a guide to good
	starting values. Overwrites existing starting values only for
	first attempt though (i_trial==0).
  */
  private boolean overrideStartingValuesFromPreviousModes(double[] timesPlusOtherParameters, LocalMode suggestedMode)
  {
	int k,indx;
	boolean found_useful_previous_mode = false;
	LayerProperties.Property prop;
	ArrayList list;
	Integer I;
	boolean useSuggestedMode = false;
	

	if (suggestedMode != null) useSuggestedMode = true;
		
	if (useSuggestedMode)
	  {
		found_useful_previous_mode = true;
		mode = suggestedMode;
	  }
	else
	  {
		// search for a mode that has the same layer configuration as
		// already suggested
		for (k=0; k< localModes.size(); k++) 
		  {
			mode = localModes.get(k);
			if (mode.fs==0 && mode.layer_code == layer_code)
			  {
				found_useful_previous_mode = true;
				break;
			  }
		  }
	  }
	if (found_useful_previous_mode)
	  {
		// ok we have a mode from fluid state 0 with the same pinchout sense. 
		// Set initial layer times guess
		// (1) Set times from the mode 'mean'.....
		for (k=0;k<_m.n_active_layers;k++) 
		  {
			timesPlusOtherParameters[k] = mode.QP_optimum.getQuick(k);
			setValuesFromIndex(k,timesPlusOtherParameters[k]);
			property_sample[k]=timesPlusOtherParameters[k];
		  }
		setValuesFromIndex(0,mu.getQuick(0));                     // top layer time always fixed
		property_sample[0]=mu.getQuick(0);
		setValuesFromIndex(_m.n_active_layers,mu.getQuick(_m.n_active_layers)); // base layer time always fixed
		property_sample[_m.n_active_layers]=mu.getQuick(_m.n_active_layers);
		// Ensure the times are ordered
		orderTimes();

		// (2) set other properties. Different indexing makes this a bit messy
		// Loop over non-time properties in source mode
		// and set those of the current fluid state to those found in source-mode
		// at the optimum. Properties in the new fluid state model vector not
		// present in the source-mode model vector are left at their values as set
		// up in the preceding pickNonTimeStartingValues() call.
		cern.colt.map.OpenIntObjectHashMap map_from_suggested_mode = map_from_index[mode.fs];

		for (k=_m.n_active_layers+1; k<map_from_suggested_mode.size(); k++) 
		  {
			list = (ArrayList)map_from_suggested_mode.get(k);
			// and the first property in the arraylist will give us the value
			prop = (LayerProperties.Property)list.get(0);
			I = (Integer)map_to.get(prop);
			if (I!=null)
			  {
				// the properties are also in the current vector.
				// set the new property values accordingly
				indx = I.intValue();
				property_sample[indx] = mode.QP_optimum.getQuick(k);
				setValuesFromIndex(indx, property_sample[indx]);
			  }
		  }
	  }

	return found_useful_previous_mode;
  }
		





  LocalMode getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods(ModeIdentifier mi)
  {
	int i,j,jj,k,n,nt,n_rows,n_rows_seismic,iVal;
	int jStartPOgroup,jEndPOgroup,jSnapPOgroup=0,indx_t;
	ModelPriorDescriptions.MetaInfo MI;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D X;
	DenseDoubleMatrix1D d=null, e=null;
	double val, newval, mean_error, delta=0, t=0, deriv;
	double negLogLik_iso=0, negLogLik_vert_stack=0, negLogPrior=0;
	double[] base_errors=null, base_isopach_errors = null, base_vertical_stacking_errors=null, new_errors;
	Integer I;	
	boolean perturbationIsNonZero = true;         // switch to save us if things are really wierd (e.g. multi pinchouts)
	boolean debug=false; // FIXME turn off & ultimately remove
	int base_layer_code;
	nt = _m.n_active_layers-1;
	n = mbar_prior_plus_isopachs[fs].size();
	LayerProperties.Property prop_t;

	// Finite difference gradient proves useful to compute here, leading to better Newton directions than
	// previous methods
	boolean COMPUTE_FD_GRADIENT = true;

	double negLogPost=0,base_negLogPost=0,base_negLogPrior=0;


	// let the initial expansion point for the isopach likelihood be the current optimal state.
	// The optimiser will have filled property_sample[] with all the values at the
	// optimum it found
	DoubleMatrix1D expansion_point = new DenseDoubleMatrix1D(property_sample);

	DoubleMatrix1D gradF=null;

	if (COMPUTE_FD_GRADIENT) gradF = expansion_point.like();

	// define a new mode with dummy entries at this point
	LocalMode new_mode = new LocalMode(0.0,mi.chisq,fs,fluid_state,mi.layer_code,null,null);

	
	// Working array
	ArrayList<LayerProperties.Property> prop_arrayl = new ArrayList<LayerProperties.Property>();

	// This update can be iterated, but a few times seems enough??
	abortNow=false;  // early escape device for next loop
	int maxNewtonSteps = 100; //10


	for (int NewtonIts=0;NewtonIts<maxNewtonSteps;NewtonIts++)
	  {
		if (_verbosity>3) _SI.EH.println("Gauss-Newton updates, iteration "+NewtonIts+":");

		// set all layer-property values to the expansion point
		// (which, on first iteration, is the posterior mean after isopach inclusion)
		for (k=0;k<n;k++)
		  {
			val = expansion_point.getQuick(k);
			setValuesFromIndex(k,val);
		  }
		if (_verbosity>3) 
		  {
			_SI.EH.printf("Gauss-Newton expansion_point: ");
			for (k=0;k<n;k++)
			   {
				_SI.EH.printf("%10.6g ",expansion_point.getQuick(k));
			  }
			_SI.EH.println("");
		  }
		orderTimes();		
		new_mode.layer_code = base_layer_code = layer_code;

		// now update full layer properties for effective vp's etc
		// Abort and reject mode if the state looks illegal
		if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false))  // FC_OK
		  {
			modeSearchSuccessFlag=1;
			return null;
		  }

		// Acquire logLikelihood contribution from seismic (+maybe isopach+ maybe vertical stack) mismatch
		n_rows_seismic = _m.SM.getErrorVector().length;
		if (true)
		  {
			// base errors are concatenation of seismic and isopachs and vertical stacking
			negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
			negLogLik_iso = computeIsopachNegLogLikelihood(false);
			negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
			n_rows =  n_rows_seismic + isopach_errors.length;
			base_errors =  _m.SM.getErrorVector();
			// get baseline isopach errors
			base_isopach_errors = (double[]) isopach_errors.clone();

			if (_SI.doStackedModel)
			  {
				negLogLik_vert_stack = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
				n_rows += vertical_stacking_errors.length;
				// get baseline vertical stacking errors
				base_vertical_stacking_errors = (double[]) vertical_stacking_errors.clone();
			  }
		  }

		if (COMPUTE_FD_GRADIENT) 
		  {
			base_negLogPrior = computePriorNegLogLikelihood(expansion_point);
			base_negLogPost = negLogLikSeisOnly + negLogLik_iso + negLogLik_vert_stack + base_negLogPrior;
		  }

		if (debug && fs==0) Util.printfileZO("errors.txt",base_errors,_m.MPD.n_isopachs);
		

		// now construct design matrix X and data d from finite differences
		// We also need the normalised errors e for the line search later on
		X = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(n_rows,n);
		e = new DenseDoubleMatrix1D(n_rows);
		d = new DenseDoubleMatrix1D(n_rows);
		for (j=0;j<n_rows_seismic;j++) d.setQuick(j,-base_errors[j]);
		for (j=0;j<n_rows_seismic;j++) e.setQuick(j,base_errors[j]);

		// initialise for isopachs too
		// for (jj=0;j<n_rows;j++,jj++)
		for (jj=0; jj<base_isopach_errors.length; j++,jj++)
		  {
			d.setQuick(j,-base_isopach_errors[jj]);
			e.setQuick(j,base_isopach_errors[jj]);
		  }
		if (_SI.doStackedModel)
		  {
			// initialise for vertical stacking too
			for (jj=0; jj<base_vertical_stacking_errors.length; j++,jj++)
			  {
				d.setQuick(j,-base_vertical_stacking_errors[jj]);
				e.setQuick(j,base_vertical_stacking_errors[jj]);
			  }
		  }



		
		// Loop over columns of X, which are properties that we want an update over
		// We could do this using the map_from_index() HashMap, but we want more explicit
		// control over *which* properties we want to use in the design matrix X.
		// Se we'll loop over all layers, and all properties, being careful not to
		// create any columns in X from "mirror" properties (i.e. only their masters)
		// We'll use the following array to detect any "mirrors"
		int[] used_indices = new int[map_to.size()];
		Arrays.fill(used_indices,-1);

		for (i=0;i<_m.n_layers;i++)
		  {
			MI = _m.MPD.MetaInfo[i];
			LPI = _m.LP[i];
		
			// This is a full list of properties that affect seismic & other constraints
			// Most are layer properties, but there are internal Seismic-model variables too
			// (see below)
			LayerProperties.Property[] prop_array =
			    {LPI.t,(i==_m.n_layers-1)?LPI.t_base:null,
				 LPI.d, LPI.LFIV,
				 LPI.rho_b, LPI.vp_b,
				 LPI.rho_oil, LPI.rho_lsg, LPI.rho_gas,
				 LPI.vp_oil, LPI.vp_lsg, LPI.vp_gas,
				 LPI.S_oil, LPI.S_lsg, LPI.S_gas,
				 LPI.NG, 
				 LPI.vp_s, LPI.vs_s, LPI.phi_s,
				 LPI.vp_m, LPI.vs_m, LPI.rho_m, 
				 ((_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL) ? 
				  LPI.float_fraction : null)
				};

			// form union of properties master list via ArrayList methods...
			prop_arrayl.clear();
			prop_arrayl.addAll((List<LayerProperties.Property>)Arrays.asList(prop_array));
			if (_SI.MF.ENABLE_HARD_ROCK_MODEL && _m.needHardRockTerms)
			  {
				prop_arrayl.addAll((List<LayerProperties.Property>)
								   Arrays.asList(LPI.NG_hard, LPI.vp_hard, LPI.vs_hard, LPI.rho_hard));
			  }
			// Get effect of systematic properties on first pass (i==0) only
			if (i == 0)
			  {
				prop_arrayl.add(_m.SM.A);  // seismic AVO systematics
				if (!_m.SM.A_B_avo_synchronised) prop_arrayl.add(_m.SM.B);
				if (_m._SI.doSimpleFluidUnitsAndContacts) prop_arrayl.addAll(_m.fluidUnitsAndContacts.fc_props);
			  }

			// Find the derivative of the logLikelihood Seismic term wrt each 
			// of these parameters by finite differences
			
			for (k=0; k<prop_arrayl.size(); k++)
			  {
				// LayerProperties.Property prop = prop_array[k];
				LayerProperties.Property prop = prop_arrayl.get(k);
				if (prop!=null && (I = (Integer)map_to.get(prop))!=null)
				  {
					// if the property in question is a mirror, skip it (we
					// don't want any duplicate rows in X)

					// get the column index
					int indx = I.intValue();
					if (used_indices[indx] != -1) 
					  {
						continue;
					  }
					else
					  {
						used_indices[indx] = indx;
					  }
                    // this property is indexed, so is part of the model, and might contribute to X

					// Skip any synched mirror properties (but not the masters) 
					// (column indexes out of range)
					if (indx>=n) continue;   

					// indx is column index in X						
					val = expansion_point.getQuick(indx);
					// Add a perturbation to the property for obtaining derivatives via finite differences
					perturbationIsNonZero = true;

					// time perturbations
					// t_base never gives reflections....
					// if ((prop == LPI.t) || (prop == LPI.t_base))   // -> see also (A)**** below
					if (prop == LPI.t)
					  {
						// Time perturbation has to be very subtle and careful
						// Also, no backwards diffs on base_time
						// (A)*****  if statement for timeperts that have t_base (buggy?)
						// if ((i>0) && (i<(_m.n_layers)) && (prop!=LPI.t_base) && (_m.LP[i].t.value <= _m.LP[i-1].t.value))
						if ((i>0) && (i<(_m.n_layers)) && (_m.LP[i].t.value <= _m.LP[i-1].t.value))
						  {
							// layer indexed [i-1] must be pinched out. 
							// => Do backward difference
							// to make perturbation maintain the pinchout state (layer code)
							// _SI.EH.warn("DOING BACKWARD DIFF!");
							// find extent of pinchout (i.e. what group of layers are involved), and safe delta if needed
							delta = -SEISMIC_DESIGN_MATRIX_TIME_EPS;
							for (j=i-1;j>=0;j--) 
							  {
								if (_m.LP[j].t.value < _m.LP[i].t.value) 
								  {
									delta = Math.max(delta, -0.99*(_m.LP[i].t.value - _m.LP[j].t.value));
									break;
								  }
							  }
							if (abs(delta) < TEN_MACHINE_EPS*_m.LP[i].t.value) // usually about 1e-15
							  {
								// configuration reconstruction probably unreliable. Set delta=0 and proceed
								delta = 0;
								perturbationIsNonZero = false;
							  }
							jStartPOgroup = j+1;
							for (j=i+1;j<=_m.n_layers;j++)
							  {
								t = (j<_m.n_layers) ? _m.LP[j].t.value : _m.LP[_m.n_layers-1].t_base.value;
								if (t > _m.LP[i].t.value) break;
							  }
							jEndPOgroup = j-1;
							// pinchout group now spans [jStartPOgroup,....,jEndPOgroup]
							// loop over time-ordering priority list (in increasing order) to find the first-snapped time
							// in the group (which we call jSnapPOgroup)
							Iterator timeStdIterator = time_std_devs[fs].values().iterator();
							while (timeStdIterator.hasNext())
							  {
								// layer number
								jSnapPOgroup = ((Integer)timeStdIterator.next()).intValue();
								if (jSnapPOgroup >= jStartPOgroup && jSnapPOgroup <= jEndPOgroup) break;
							  }
							if (i != jSnapPOgroup)
							  {
								delta = 0;          // derivatives should be zero through action of layer-ordering
								newval = val;       // and our desire to move in a direction that maintains the layer_config
								perturbationIsNonZero = false;
							  }
							else
							  {
								// OK, now *all* times from [jSnapPOgroup,...,jEndPOgroup] must be backwards-differenced
								// to preserve the layercode/pinchout configuration
								// set all lower values to have same backward difference
								newval = val+delta;
								for (j=jSnapPOgroup; j<=jEndPOgroup; j++)
								  {
									prop_t = (j<_m.n_layers) ? _m.LP[j].t : _m.LP[_m.n_layers-1].t_base;
									indx_t = ((Integer)map_to.get(prop_t)).intValue();
									// Now perturb the time value. Time properties require special consideration...
									setValuesFromIndex(indx_t,newval);
								  }
							  }
						  }
						else
						  {
							// Do a forward difference, but be careful not to cross the next layer
							// or any subsequent pinchouts
							// find next non-pinched out time
							delta = SEISMIC_DESIGN_MATRIX_TIME_EPS;
							if (prop == LPI.t)
							  {					  
								double t_i = _m.LP[i].t.value;
								for (j=i+1;j<=_m.n_layers;j++)
								  {
									t = (j<_m.n_layers) ? _m.LP[j].t.value : _m.LP[_m.n_layers-1].t_base.value;
									if (t > t_i) 
									  {
										if (j==i+1)
										  {
											// most common; no in-between layers
											delta = Math.min(delta,0.99*(t - t_i));
										  }
										else
										  {
											// choose delta between last layers t_i....t_{j-1}...t_half....t
											double t_half = 0.5*(_m.LP[j-1].t.value + t);
											delta = Math.min(delta,(t_half - t_i));
										  }
										break;
									  }		
								  }	
								if (delta<0) 
								  {
									if (_verbosity>1) 
									  {
										_SI.EH.println("Bad fwd difference delta in "
													   +"getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods(). FIXME!!!");				
									  }
								  }
								else if (delta < TEN_MACHINE_EPS*t_i) // usually about 1e-15
								  {
									// configuration reconstruction probably unreliable. Set delta=0 and proceed
									delta = 0;
								  }
								if (delta==0) perturbationIsNonZero = false;
							  }
							// else delta = SEISMIC_DESIGN_MATRIX_TIME_EPS for t_base, always valid
							newval = val+delta;
							// Now perturb the time value. Time properties require special consideration...
							setValuesFromIndex(indx,newval);
						  }

						// orderTimes() can alter more than the current time,
						// so we have to set all of them to be sure
						orderTimes();
						if (base_layer_code!=layer_code)
						  {
							// time perturbations were ill-constructed. A bug. Fix it
							_SI.EH.println("Time perturbations in getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods() "+
										   "alter pinchout configuration on layer "+(i+1)+"!! Sensitivity may be wrong. Fix this up!?");
							_SI.EH.print(" time config to perturb: ");
   							for (k=0;k<=_m.n_layers;k++)
							  {
								_SI.EH.printf("%10.8g ",((k<_m.n_layers) ? _m.LP[k].t.value : _m.LP[_m.n_layers-1].t_base.value));
							  }
							_SI.EH.println("");
							// Desperate. Set zero difference, restore state
							delta=0;
							newval = val;
							// Set value..
							setValuesFromIndex(indx,newval);
							// Reorder
							orderTimes();
							if (base_layer_code!=layer_code)
							  {
								// time perturbations hopeless. A bug. Fix it
								_SI.EH.die_err("Fatal time perturbations bug. Fix this up!?");
							  }
							perturbationIsNonZero = false;
						  }
						if (_SI.doSimpleFluidUnitsAndContacts)
						  {
							// time perturbation will affect depths and => internal layer FC's
							_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
						  }
					  }
					else
					  {
						// Non-time properties
						if (prop.normalised) // (0<x<1)
						  {
							// this is probably better for normalised (0<x<1) properties
							if ((newval = val+SEISMIC_DESIGN_MATRIX_EPS) > 1.0)
							  {
								// too close to one, so prop -> prop+eps  doesn't work.
								// try  prop -> prop -eps
								newval = val-SEISMIC_DESIGN_MATRIX_EPS;
							  }
						  }
						else
						  {
							newval = val*(1.0+SEISMIC_DESIGN_MATRIX_EPS);
						  }
						if (newval == 0.0)
						  {
							// If prop was zero, (1+eps)*prop=0 so we'll change this to prop->eps
							// Will very likely only ever happen with normalised properties, so scaling is ok
							newval=SEISMIC_DESIGN_MATRIX_EPS;
						  }
						// Now perturb the property's value. Time properties require special consideration...
						setValuesFromIndex(indx,newval);
						// We'll have to update all layers, since others might be synchronised to this one
						// A change in times won't affect properties, so we can skip this if (prop==LPI.t)
						// Abort and reject mode if the state looks illegal
						if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false)) // FC_OK	
						  {
							// If this failed, we are sitting at the edge of the physically real
							// region, and the forward difference took us over the edge.
							// If we go in opposite directions, the backward direction points should 
							// be sensible. Let's try the backwards difference.
							// set newval to a backwards difference point...
							newval = 2.0*val - newval;
							setValuesFromIndex(indx,newval);
							if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false)) // FC_OK
							  {
								// This should never happen. We may be sitting at the edge of the physically real
								// region, but at least one of the forward
								// or backward direction points should be sensible
								_SI.EH.warn("Gassman failed both forwards and backward differences. A disaster.");
								modeSearchSuccessFlag=1;
								return null;
							  }
						  }
					  }



					// compute new vector of seismic errors under this perturbation
					negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
					new_errors =  _m.SM.getErrorVector();
					
					if (newval==val)
					  {
						// perturbations were ill-constructed. A bug. Fix it
						if (prop == LPI.t)
						  {
							// Zero perturbations are now legit for times.
							// _SI.EH.warn("Warning: Time perturbations in getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods() "
							//										   +"badly constructed. Sensitivity has been set to 0. Fix this bug!?");
						  }
						else
						  {
							_SI.EH.warn("Warning: Property (not time) perturbations in "
										   +"getLocalModeByUpdatingMGSystemUsingSeismicAndOtherLikelihoods() "
										   +"badly constructed. Sensitivity has been set to 0. Fix this bug!?");
						  }
					  }
					// Protection when finite-differences don't seem to work (e.g. multi-pinchouts)
					double d_prop = (perturbationIsNonZero ? (newval-val) : Double.MAX_VALUE);

					// Fill in sensitivity matrix

					// Seismic
					for (j=0;j<n_rows_seismic;j++)
					  {
						deriv = (new_errors[j] - base_errors[j])/d_prop;
						X.setQuick(j,indx,deriv);
						// and update RHS of likelihood constriant
						d.setQuick(j, d.getQuick(j) + deriv * val);
					  }
					// Isopachs: continue on through system, using isopach errors
					negLogLik_iso = computeIsopachNegLogLikelihood(false);
					new_errors = isopach_errors;
					// for (jj=0;j<n_rows;j++,jj++)  // jj is index of isopach error, j continues on through rows of X
					for (jj=0; jj<new_errors.length; j++,jj++)  // jj is index of isopach error, j continues on through rows of X
					  {
						deriv = (new_errors[jj] - base_isopach_errors[jj])/d_prop;
						X.setQuick(j,indx,deriv);
						// and update RHS of likelihood constraint: d = {\tilde X}.m_0 - e(m_0)
						d.setQuick(j, d.getQuick(j) + deriv * val);
					  }
					// Vertical stacking
					if (_SI.doStackedModel)
					  {
						negLogLik_vert_stack = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
						new_errors = vertical_stacking_errors;
						for (jj=0; jj<new_errors.length; j++,jj++)  // jj is index of vertical_stacking_error, j continues on through rows of X
						  {
							deriv = (new_errors[jj] - base_vertical_stacking_errors[jj])/d_prop;
							X.setQuick(j,indx,deriv);
							// and update RHS of likelihood constraint: d = {\tilde X}.m_0 - e(m_0)
							d.setQuick(j, d.getQuick(j) + deriv * val);
						  }
					  }

					if (COMPUTE_FD_GRADIENT)
					  {
						expansion_point.setQuick(indx,newval);
						negLogPrior = computePriorNegLogLikelihood(expansion_point);						
						negLogPost = negLogLikSeisOnly + negLogLik_iso + negLogLik_vert_stack +  negLogPrior;
						gradF.setQuick(indx,(negLogPost - base_negLogPost)/d_prop);
						// restore
						expansion_point.setQuick(indx,val);	
					  }
					// restore value
					if (prop==LPI.t || prop==LPI.t_base)
					  {
						// if we've perturbed a time(s), restore all of them & make sure ordering is ok
						for (j=0;j<=nt+1;j++) setValuesFromIndex(j,expansion_point.getQuick(j));
						orderTimes();
						if (_SI.doSimpleFluidUnitsAndContacts)
						  {
							// restore internal states
							if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false)) // FC_OK	
							  {
								modeSearchSuccessFlag=1;
								return null;
							  }
						  }
					  }
					else
					  {
						setValuesFromIndex(indx,val);
						// Restore all layers' effective properties, since others might be synchronised to this one
						// A change in times won't affect properties, so we can skip this if (prop==LPI.t)
						// Abort and reject mode if the state looks illegal
						if (!_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false)) // FC_OK	
						  {
							modeSearchSuccessFlag=1;
							return null;
						  }
					  }
				  } // property contributes (if prop!=null etc)
			  }// loop over possible properties
		  } // loop over layers

		used_indices=null;
		base_isopach_errors = null;
		base_vertical_stacking_errors = null;

		// now we can perform the linearised Bayesian update
		if (!updateMGSystemForSeismic(X,d,e,new_mode,expansion_point,gradF))
		  {
			modeSearchSuccessFlag=1;
			return null;
		  }
		// update the expansion point to the current mode
		expansion_point.assign(new_mode.QP_optimum);

		// check to see if we detected a degenerate (albeit successful) backtrack
		if (abortNow) 
		  {
			modeSearchSuccessFlag=0; // I think this is the right thing to do
			if (_verbosity>3)
			  {
				_SI.EH.println("Newton terminated due to degenerate line-search");
			  }
			break;
		  }

		// termination criterion
		/* FIXME we have to get Newton to behave better....*/
		if ((2.0*Math.abs(delta_fNewNewton) < NEWTON_EPS * (Math.abs(fNewNewton) + EPS_FLOATINGPOINT)))
		  {
			modeSearchSuccessFlag=0;
			if (_verbosity>3)
			  {
				_SI.EH.println("Newton terminated after "+NewtonIts+" iterations");
			  }
			// we should...FIXME we have to get Newton to behave better...
			break;
		  }
		else
		  {
			// not converged. If we haven't done this by maxNewtonSteps something's awry
			if (NewtonIts==(maxNewtonSteps-1))
			  {
				modeSearchSuccessFlag=2;
				if (_verbosity>2) 
				  {
					_SI.EH.println("Unsatisfactory Newton termination criterion");  
				  }
				// and probably the mode's not worth trusting....
				return null;
			  }
		  }
	  }

	// for efficiency reasons, change mode.R to a dense matrix if it isn't sufficiently sparse
	double sparsity;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D Rnew;
	
	if (new_mode.R == null)
	  {
		_SI.EH.println("No decomposition found");
	  }
	sparsity = new_mode.R.cardinality()/(double)new_mode.R.size();
	if (sparsity > SPARSITY_THRESHOLD)
	  {
		Rnew = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(new_mode.R.rows(),new_mode.R.columns());
		Rnew.assign(new_mode.R);
		new_mode.R = null; // discard old sparse version
		new_mode.R = Rnew;
	  }
   	if (false)
	  {
		new_mode.Rtransp = new_mode.R.viewDice();
	  }
	else
	  {
		// make new_mode.Rtransp efficient for lower tri mult in RWS...
		new_mode.Rtransp = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(new_mode.R.rows(),
																					   new_mode.R.columns());
		new_mode.Rtransp.assign(new_mode.R.viewDice());
	  }

	// save also scaled seismic error
	// m.MS.negLogLikSeisOnly = 0.5*sum {errors^2}, so undo factor of 2
	new_mode.seismic_chisq = 2.0*negLogLikSeisOnly/
	  (double)(_m.SM.seismic_error_vector.length * _m.SM.delta_t_ratio_correction_factor);

	// these are useful for --MAP runs
	// computed using relation 
	// chisq = negLogLikSeisOnly + negLogIsopachLikelihood + negLogDepthLikelihood + negLogVStackLikelihood + negLogPrior
	// in updateMGSystemForSeismic()
	new_mode.lastAcceptedNegLogAllLikelihoods = new_mode.chisq - negLogPrior;
	new_mode.lastAcceptedNegLogLikPriorOnly = negLogPrior;
	// these pieces are constant for a particular mode, so we initialise here
	new_mode.normalisationConstants = _m.SM.logDetCnoiseAnd2PiPieces        // normalisation from seismic
	                                  + logDetCisopachAnd2PiPieces          // normalisation from isopachs
	                                  + logDetCvsAnd2PiPieces               // normalisation from vertical stacking *likelihood* (thicknesses)
	                                   - 0.5 * (new_mode.logDetCtilde  - logDetC_prior_only[fs]);  



	// and logDetCprior
	new_mode.logDetCprior = logDetC_prior_only[fs];

	// prior weight depends only on fluid state prior probability. Use actual probability for this
	// (its floating point range is not large)
	new_mode.prior_weight = fluid_state.prob;

	// store also unnormalised marginal model likelihood
	// before reweightings etc occur
	new_mode.logMargModelLikelihood = new_mode.weight;

	// Now discount mode weight by fluid state prior probability
	// Overall weight must be in proportion to fluid_state combination probability
	// adding this to weight = log(mode_weight) effects this change
	// (we are working with the -ve log of the mode weight...)

	new_mode.weight -= Math.log(fluid_state.prob);

	// and layer code
	new_mode.layer_code = layer_code;

	return new_mode;
  }

  private boolean updateMGSystemForSeismic(DoubleMatrix2D X, DoubleMatrix1D d, DoubleMatrix1D e, LocalMode new_mode,
										   DoubleMatrix1D current_expansion_point, DoubleMatrix1D grad)
  {
	// Use the standard inverse theory updating formulas ala 
	// Tarantola, but with X the
	// design matrix formed from the seismic likelihood,
	// to get new MG distribution
	// Supplied gradient grad may be used for a line search if Newton direction doesn't work

	int i,ii,j,k,n_params;
	double val,rhs,s,chisq=0;
	double weight,logDetCtilde;
	boolean sample_succeeded, closelyMatchedLayerTimesAreLikely,useQP=false;
	DoubleMatrix1D m_tilde;
	ModifiedCholeskyRightDecomposition CDR;
	DoubleMatrix2D Rinv;
	DoubleMatrix2D Cinv;
	DoubleMatrix1D mbar;
	double[] CtildeInvDiag;
	boolean debug = false;     // FIXME turn off & ultimately remove
	IntArrayList pinchedOrSqueezedLayers = new IntArrayList();

	char[] NGstate = new char[_m.n_layers];

	// allow line search to use steepest descent direction if Newton is bad
	boolean useSteepestDescentIfRequired = true;

	// if no supplied gradient direction, recompute it
	boolean gradientSupplied = (grad!=null);

	DoubleMatrix1D gradF=null;

	n_params = mbar_prior_only[fs].size();

	// Force use of dense matrices, since following operations are n^3 anyway...
	Cinv = new DenseDoubleMatrix2D(n_params,n_params);
	
	// Perform updates on pure Gaussian prior. Likelihood incorporates isopachs + seismic
	mbar = mbar_prior_only[fs];
	// new matrix for the local mode to hold the new mode Cholesky R ultimately		
	Cinv.assign(Cinv_prior_only[fs]);  // copy Cinv_pr[fs] into Cinv
  
	m_tilde = new DenseDoubleMatrix1D(mbar.size()); 

	if (debug && fs==0)   // FIXME COMMENT OUT
	  {
		m_tilde.assign(mbar);                         // copy mbar_prior_plus_isopachs into mtilde
		Util.debugPrintMatrixToFile("Cinv.txt",Cinv);
		Util.debugPrintMatrixToFile("mbar.txt",mbar);
		Util.debugPrintMatrixToFile("d.txt",d);
		Util.debugPrintMatrixToFile("X.txt",X);
	  }


	if (!gradientSupplied)	  
	  {
		// Numerical finite difference gradient isn't available, use the theory to get m_tilde:
		// (1) form m_tilde = C^-1.mbar. Can use symmetric matrix-vector multiply dsymv
		SeqBlas.seqBlas.dsymv(false,1.0,Cinv,mbar,0.0,m_tilde);
		
		// at this point m_tilde = C^-1.m_bar
		// (2) Form m_tilde = xi = X^T.d + C^-1.mbar   *in place*
		SeqBlas.seqBlas.dgemv(true,1.0,X,d,1.0,m_tilde);
	  }

	// form {\tilde Cinv} = (X^T.X+C^-1) *in place* (Cinv is created with C^-1 already in place)
	SeqBlas.seqBlas.dgemm(true,false,1.0,X,X,1.0,Cinv);

	// Now at this point Cinv contains {\tilde C}_{s}^{-1}

	if (gradientSupplied)
	  {
		// form m_tilde from supplied gradient here
		// Roundoff errors studies show m_tilde is much better formed from the numerical
		// finite difference gradient (yields better Newton directions).
		// m_tilde = xi = -gradF + {\tilde C}_{s}^{-1}.m_0
		gradF = grad;
		m_tilde.assign(grad);
		SeqBlas.seqBlas.dgemv(false,1.0,Cinv,current_expansion_point,-1.0,m_tilde);
	  }

	if (debug && fs==0)  // FIXME COMMENT OUT
	  {		
		Util.debugPrintMatrixToFile("rhs."+fs+".txt",m_tilde);
		Util.debugPrintMatrixToFile("Ctilde_inv."+fs+".txt",Cinv);
	  }

	// Now at this point Cinv contains {\tilde C}_{s}^{-1}
	// And m_tilde contains (X^T.d + C^-1.mbar) "the RHS"
	DoubleMatrix1D RHS = m_tilde.copy();

	if (!gradientSupplied)
	  {
		// Gradient calculation. Some appreciable roundoff error here 
		// due to formula for gradF containing cancelling terms.
		// At this point m_tilde = xi = X^T.d + C^-1.mbar 
		gradF = m_tilde.copy();
		// form gradF = - xi + {\tilde C}_{s}^{-1}.m_0
		SeqBlas.seqBlas.dgemv(false,1.0,Cinv,current_expansion_point,-1.0,gradF);
		// Util.debugPrintMatrixToFile("gradF.txt",gradF);
		// Util.debugPrintMatrix("gradF.txt",gradF,10);		
	  }

	// store diagonal entries before we destroy them...
	CtildeInvDiag = new double[n_params];
	for (i=0; i<n_params; i++) CtildeInvDiag[i] = Cinv.getQuick(i,i);



	// get right-Cholesky factor of Cinv
	CDR = new ModifiedCholeskyRightDecomposition(Cinv);
	// overwrite it on R
	Rinv = CDR.getR();





	if (debug && fs==0)  // FIXME COMMENT OUT
	  {		
		Util.debugPrintMatrixToFile("Rtilde_inv.txt",Rinv);
	  }

	// now take inverse (in place)
	CDR.getCholeskyRightInverse(Rinv);



	// At this point Rinv is the transpose of the upper triangular R Cholesky
	// factor of {\tilde C}, viz {\tilde C} = Rinv^T . Rinv
	// reconstructLayerBasedCovariance(Rinv);// FIXME COMMENT OUT

	// form {\tilde C}.m_tilde =  (Rinv^T.Rinv).m_tilde *in place*. 
	// Later on new_mode.R will point to Rinv if the step improves things
	// (1) form m_tilde -> Rinv.m_tilde *in place*
	for (i=0;i<n_params;i++)
	  {
		s = 0;
		for (j=i;j<n_params;j++)
		  {
			s += Rinv.getQuick(i,j) * m_tilde.getQuick(j);
		  }
		// write over m_tilde in place, as last element is no longer needed
		m_tilde.setQuick(i,s);
	  }
	
   
	// (2) now Rinv^T.m_tilde
	// and also the (log)determinant of the new prior Ctilde
	// NO LONGER NEED det_ratio=1.0;
	logDetCtilde = 0.0;
	for (i=n_params-1;i>=0;i--)
	  {
		s = 0;
		for (j=0;j<=i;j++)
		  {
			// notice transpose!
			s += Rinv.getQuick(j,i) * m_tilde.getQuick(j);
		  }
		// write over m_tilde in place in ith element as no longer needed)
		m_tilde.setQuick(i,s);
		// we'll need the determinant of the prior covariance (subtract since R is an inverse)
		// and also the determinant of the updated prior covariance {\tilde C}
		logDetCtilde += Math.log((s=Rinv.getQuick(i,i)) * s);
	  }

	if (debug && fs==0)   // FIXME COMMENT OUT
	  {
		Util.debugPrintMatrixToFile("m_tilde.pre-backtrack.txt",m_tilde);
	  }

	// Now compute new mean seismic error to verify that linearisation is valid
	// set all layer-property values to the Newton point m_tilde
	for (k=0;k<n_params;k++) setValuesFromIndex(k,m_tilde.getQuick(k));
	orderTimes();

	// Now the new mode may have forced a pinchout or changed the pinchout configuration.	
	
	// little fragment of code to detect minima living near pinchout...
	// Time tests...
	closelyMatchedLayerTimesAreLikely = false;
	for (j=1;j<_m.n_active_layers-1;j++) 
	  {
		i = master_layer_number_sequence[j];
		ii = master_layer_number_sequence[j+1];
		if (Math.abs(_m.LP[ii].t.value - _m.LP[i].t.value)<= TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER)
		  {
			pinchedOrSqueezedLayers.add(j+1);
			closelyMatchedLayerTimesAreLikely = true;
		  }
	  }
	// NG tests.
	if (false && !closelyMatchedLayerTimesAreLikely)
	  {
		// test also for NG near endpoints...
		Integer I;
		for (i=0,k=0;i<_m.n_active_layers;i++) 
		  {
			LayerProperties.Property prop = _m.LP[i].NG;	
			if (prop!=null)
			  {
				if ((I = (Integer)map_to.get(prop)) != null)
				  {
					j = I.intValue();
					double NG = current_expansion_point.getQuick(j);
					double NGnew = m_tilde.getQuick(j);
					if (((NG<=NG_EPSILON) && (NGnew<=NG_EPSILON))
						|| ((NG >= (1.0-NG_EPSILON)) && (NGnew >= (1.0-NG_EPSILON))))
					  {
						closelyMatchedLayerTimesAreLikely= true;
						break;
					  }
				  }
			  }
		  }
	  }

    if (((layer_code != new_mode.layer_code || closelyMatchedLayerTimesAreLikely)))
	  {
		// Suspect that updated curvature is about right, but new_mode mean is wrong. 
		// (The minima is a pinchout & not parabolic. There may be a big jump across to the
		// new mode, so the extrapolated mode.mean could be very wrong. Use the
		// updated curvature but get the mode location from Quadratic programming.
		if (_verbosity>3) 
		  {
			_SI.EH.print("Seismic update forces closely matched layer times or pinchouts: layer(s)");
			for (i=0;i<pinchedOrSqueezedLayers.size();i++) _SI.EH.print(" "+pinchedOrSqueezedLayers.get(i));
			_SI.EH.println(". Obtaining quadratic programming optimum...");
		  }
		// we will need the projected new mean for the Quadratic programming stuff
		new_mode.QP_optimum = m_tilde;
		// quadratic programming (QP) will set new_mode.QP_optimum to optimal point in current layer config.
		updateModeUsingConstrainedQuadraticProgramming(new_mode,current_expansion_point,m_tilde,Rinv);
		// set new point to optimal point along line from current_expansion_point 
		// to Quadratic programming update QP_optimum (hopefully QP_optimum if Gassman doesn't fail)
		m_tilde = getBestLegalGassmanPointAlongCurrentLine(current_expansion_point,m_tilde);
		new_mode.QP_optimum = m_tilde;
		useQP=true;
	  }
	else
	  {
		// No pinchout flip forced going from current_expansion_point to Newton update m_tilde.
		// Set mode at best Gassman point along current direction
		m_tilde = getBestLegalGassmanPointAlongCurrentLine(current_expansion_point,m_tilde);
		// Actual quadratic mode can point to the mode mean
		// and point mode members at the bits we derived
		new_mode.QP_optimum = m_tilde;
		new_mode.quadratic_mean = new_mode.QP_optimum;
	  }
	pinchedOrSqueezedLayers=null; // GC

	// BACKTRACKING

	if (false)
	  {
		// We need an equivalent of the Wolfe conditions for this problem, for
		// the unconstrained as well as the QP case. Here's some of the pieces
		// (1) We compute  grad(f_k)^T.p_k described in Nocedal & Wright, "Numerical optimization", 
		// Procedure 3.1. The "c" and "alpha" in N & W will probably be 0.001 & 1
		// This is as yet incomplete....

		// we need a gradient for the Wolfe conditions in the backtracking method
		DoubleMatrix1D m_grad = m_tilde.copy();
		// subtract: grad = m_tilde-mbar
		m_grad.assign(mbar_prior_plus_isopachs[fs], Functions.minus);
		// Form grad = -Cp^{-1}.(m-mbar) + X^T.e
		SeqBlas.seqBlas.dsymv(false,-1.0,Cinv_prior_plus_isopachs[fs],m_grad,0,m_grad);
		// SeqBlas.seqBlas.dsymv(false,0,Cinv_prior_plus_isopachs[fs],m_grad,0,m_grad);
		SeqBlas.seqBlas.dgemv(true,1.0,X,e,0,m_grad);
		// Form step vector
		DoubleMatrix1D delta_m = m_tilde.copy();
		delta_m.assign(current_expansion_point, Functions.minus);
		// and the dot product (grad(f).p_k)
		_SI.EH.println("QP Wolfe = "+SeqBlas.seqBlas.ddot(m_grad,delta_m));

	  }


	double frac=0,frac_start,frac_end,frac_step;
	DoubleMatrix1D xmid;
	DoubleMatrix1D projDir=null;
	double negLogIsopachLikelihood, negLogDepthLikelihood, negLogVStackLikelihood=0.0, negLogPrior, f_current;
	// value at current point
	setVariablesAtThisLegalGassmanPoint(current_expansion_point);
	orderTimes();
	// acquire total error
	negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
	negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
	negLogIsopachLikelihood = computeIsopachNegLogLikelihood(false);
	// Add in logLikelihood contribution from DEPTH self-consistency if required
	negLogDepthLikelihood = computeConsistentDepthNegLogLikelihood(false); 
	// Add in logLikelihood contribution from vertical stacks
	if (_SI.doStackedModel) negLogVStackLikelihood = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 

		
	negLogPrior = computePriorNegLogLikelihood(current_expansion_point);
	f_current = negLogLikSeisOnly + negLogIsopachLikelihood + negLogDepthLikelihood + negLogVStackLikelihood + negLogPrior;
	if (_verbosity>3) 
	  {
		_SI.EH.printf((useQP?"QP_":"")+"Backtrack:               Pre-update overall f = %-10.8g, prior = %-10.8g, isopach = %-10.8g, "+
					  "depth self-consistency = %-10.8g, vert-stack = %-10.8g, seismic NLL =  %-10.8g, layer code=%d%n",
					  f_current, negLogPrior, negLogIsopachLikelihood,
					  negLogDepthLikelihood, negLogVStackLikelihood, negLogLikSeisOnly, layer_code);
	  }
	// Now backtrack from m_tilde back to current_expansion_point.
	// This is a terrible approximate method: it ought to be the Wolfe conditions, but this is a first pass
	// the loop *must* contain the null step frac=0 also (must satisfy frac < BACKTRAP_EPSILON)

	boolean usingSteepestDescent = false;
	double steepestDescentScale;

	// Steps to move toward projected point (either Newton step or QP step)
	boolean successfulBacktrack = false;	
	double fStepMin=Double.MAX_VALUE,fracMin=1.0;
	xmid = current_expansion_point.copy();

	frac_start=1.0;   // should be the Newton point on a smooth surface
	frac_end=0.0;    // If -ve, we may even go backward if e.g. the QP functions 
	                 // have screwed up (i.e. line is uphill...). Sigh.. Desperate.
	frac_step=-0.1;

	if (_verbosity>3) 
	  {
		// _SI.EH.printf("Backtracking: fraction, Post-update overall f, prior,      isopach, depth self-consistency, vert-stack, seismic NLL, layer code%n");
		_SI.EH.printf("Backtracking: %-10s %-25s %-10s %-10s %-10s %-10s %-10s %-11s %-20s%n",
					  "fraction","Post-update-overall-f","prior","isopach","depth-SC","vert-stack","seismic-NLL","layer-code","NG-state");
	  }

	boolean exitSearch=false;

	PASS:
	for (byte pass=0;pass<=3;pass++)
	  {
		// pass = 0,1,2 = {Newton, steepest descent (SD), subspace SD}	
		// pass = 3 = possible last exitSearch to populate model after pass=2
		for (frac=frac_start; frac>=frac_end; frac+=frac_step)
		  {	
			// backtrack
			for (k=0;k<xmid.size();k++)
			  {
				xmid.setQuick(k,frac*m_tilde.getQuick(k)+(1.0-frac)*current_expansion_point.getQuick(k));
			  }
			xmid = getBestLegalGassmanPointAlongCurrentLine(current_expansion_point,xmid);

			negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
			negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
			negLogIsopachLikelihood = computeIsopachNegLogLikelihood(false);
			// Add in logLikelihood contribution from DEPTH self-consistency if required
			negLogDepthLikelihood = computeConsistentDepthNegLogLikelihood(false);
			// Add in logLikelihood contribution from vertical stacks
			if (_SI.doStackedModel) 
			  {
				negLogVStackLikelihood = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false);  
			  }
			negLogPrior = computePriorNegLogLikelihood(xmid);
			chisq = negLogLikSeisOnly + negLogIsopachLikelihood + negLogDepthLikelihood + negLogVStackLikelihood + negLogPrior;
			if (_verbosity>3) 
			  {
				/*
				  _SI.EH.printf("Backtracking: %-10s %-25s %-10s %-10s %-10s %-10s %-10s %-11s %-20s%n",
				  "fraction","Post-update-overall-f","prior","isopach","depth-SC","vert-stack","seismic-NLL","layer-code","NG-state");
				*/
				for (i=0;i<_m.n_layers;i++) NGstate[i] = _m.LP[i].NGstate();
				_SI.EH.printf("              %-10.4g %-25g %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g %11d %s%n",
							  frac, chisq, negLogPrior, negLogIsopachLikelihood,
							  negLogDepthLikelihood, negLogVStackLikelihood, negLogLikSeisOnly, layer_code, new String(NGstate));
			  }

			// QUASI EXHAUSTIVE LINE SEARCH

			// Final single passthrough to populate model with best solution
			if (exitSearch) 
			  {
				if (fStepMin >= f_current || frac<BACKTRAP_EPSILON) successfulBacktrack = false;
				m_tilde.assign(xmid);
				break PASS;
			  }

			// store min
			if (chisq<=fStepMin)
			  {
				fStepMin = chisq;
				fracMin = frac;
			  }

			// Last point checks
			// check if fwd direction has been any good in Newton direction. 
			// Stop when frac=0 and a degenerate line search (best=(frac=0))
			if ( Math.abs(frac) < BACKTRAP_EPSILON)   // nearest-point in line search
			  {
				// Testing steps
				if (pass==0)
				  {
					// End of Newton pass: check for suspiciously small decrease in Newton direction. If so, flip to SD
					if (2.0*Math.abs(fStepMin - f_current) < NEWTON_EPS * (Math.abs(f_current) + EPS_FLOATINGPOINT))
					  {
						// flip to pure steepest descent direction
						usingSteepestDescent = true;
						if (_verbosity>3) _SI.EH.print("Flipping line search to steepest descent...%n");
					  }
					else
					  {
						exitSearch = true; // prepare to jump out
						successfulBacktrack = true;
					  }
				  }
				if (pass==1)
				  {
					if (fracMin == frac)  // degenerate line search
					  {
						successfulBacktrack = false;
						exitSearch = false;
					  }
					else exitSearch = true;
				  }
				if (pass==2)
				  {
					exitSearch = true;    // will guarantee finishing
					if (fracMin == frac)  // degenerate line search
					  {
						successfulBacktrack = false;
					  }
					else successfulBacktrack = true;
					  
				  }
				if (!exitSearch && usingSteepestDescent && useSteepestDescentIfRequired)
				  {
					// Setting up steps for SD passes
					if (_verbosity>3) _SI.EH.printf("steepest descent pass %d...%n",(pass+1));
					if (pass==0)
					  {
						// Set up steepest descent  
						// form new mtilde...
						projDir = gradF.like();
					  }
					// point downhill:
					for (k=0;k<n_params;k++) projDir.setQuick(k,-gradF.getQuick(k));
					if (pass==0)
					  {
						// freeze times on first pass - work in subspace
						for (k=0;k<=_m.n_active_layers;k++) projDir.setQuick(k,0.0);
					  }
				
					steepestDescentScale = computeSteepestDescentScale(projDir,gradF,Cinv, (pass==0));
					for (k=0;k<n_params;k++)
					  {
						m_tilde.setQuick(k,steepestDescentScale*projDir.getQuick(k) + current_expansion_point.getQuick(k));
					  }
					// project along pinchout and NG disconuities if required
					// experimentalProjectSDdirectionOntoQuadraticProgrammingConstraint(new_mode,current_expansion_point,m_tilde);

					// debugShowDeltaMComponents(new_mode,current_expansion_point,m_tilde);

					// force restart to loop, operating from 1.0...0.0 in Steepest descent direction
					fStepMin = Double.MAX_VALUE;
					fracMin = 1.0;
				  }
			  }
			
			// force re-pass on last time
			if (exitSearch) // last entry in loop
			  {
				frac = fracMin - frac_step; // so increment in loop brings it back to fracMin....
			  }
		  }
	  }

	// flag a degenerate (albeit successful) backtrack
	if (Math.abs(frac)<BACKTRAP_EPSILON) abortNow= true;
	  
	if (debug && fs==0)   // FIXME COMMENT OUT
	  {
		Util.printfileZO("new_errors.txt",_m.SM.getErrorVector(),_m.SM.getErrorVector().length);
		Util.debugPrintMatrixToFile("m_tilde.txt",m_tilde);
		Util.debugPrintMatrixToFile("R_iso+seis.txt", R_prior_plus_isopachs[fs]);
		Util.debugPrintMatrixToFile("R_tilde.txt",Rinv);
	  }

	// Mode weight must be proportional to likelihood fits, but
	// also should contain a covariance (from integrating to get the marginal) factor
	// Let's assert that everything is taken into account by the following Ansatz:
	// (start with the -ve log of the mode chisq...)

	// First cut at -log(marginal model likelihood)
	weight = chisq  
	  + _m.SM.logDetCnoiseAnd2PiPieces        // normalisation from seismic
	  + logDetCisopachAnd2PiPieces            // normalisation from isopachs
	  + logDetCvsAnd2PiPieces                 // normalisation from vertical stacking *likelihood* (thicknesses)
	  - 0.5 * (logDetCtilde  - logDetC_prior_only[fs]);  // note (1/2)log(2Pi) factors cancel here
	//_SI.EH.printf("WEIGHT: %8g %8g %8g %8g %8g %8g %8g %8g%n",weight,chisq,_m.SM.logDetCnoiseAnd2PiPieces,logDetCisopachAnd2PiPieces,logDetCvsAnd2PiPieces,- 0.5 * (logDetCtilde  - logDetC_prior_only[fs]),logDetCtilde,logDetC_prior_only[fs]);

	// set the status variables
	fNewNewton = chisq;
	delta_fNewNewton = (chisq - new_mode.chisq);

	// Has the step been an improvement? If so, update the mode variables
	if (chisq <= new_mode.chisq || successfulBacktrack)
	  {
		new_mode.logDetCtilde = logDetCtilde;
		new_mode.R = Rinv;
		new_mode.CtildeInvDiag = CtildeInvDiag;
		new_mode.QP_optimum = m_tilde; 
		new_mode.weight = weight;
		new_mode.chisq = chisq;
		new_mode.layer_code = layer_code;
	  }
	// but if it was a degenerate backtrack, try this
	if (!successfulBacktrack)
	  {
		new_mode.logDetCtilde = logDetCtilde;
		new_mode.R = Rinv;
		new_mode.CtildeInvDiag = CtildeInvDiag;
		new_mode.QP_optimum = m_tilde;
		new_mode.weight = weight;
		new_mode.chisq = chisq;
	  }

	RHS=null;
	CDR = null;



	if (debug)
	  {
		// form {\tilde Cinv} = (Rinv^T.Rinv) *in place* 
		DoubleMatrix2D Cpost = Rinv.like();
		SeqBlas.seqBlas.dgemm(true,false,1.0,Rinv,Rinv,0,Cpost);
		Util.debugPrintMatrixToFile("Cpost.txt",Cpost);
	  }
	
	// if we got here, things look ok
	return true;
  }

  private double computeSteepestDescentScale(DoubleMatrix1D delta_m, DoubleMatrix1D gradF, DoubleMatrix2D Cinv, boolean hasZeroTimeComponent)
  {
	int i;
	double alpha_safe, gradMax=0;
	double kappa,alpha,beta,disc;

	// Use parabolic approx to surface:
	// delta_chisq = a |delta_m . grad F| +0.5 a^2 delta_m^T.Cinv.delta_m = -beta a +kappa a^2
	// solve for minima in a, or a reqired to generate delta_chisq = DELTA_CHISQ_TARGET
	// use double the minima range, for wise bracketing (2 a_min)

	double DELTA_CHISQ_TARGET = 1.0;

	// form kappa = 0.5 gradF^T.Cinv.gradF	
	DoubleMatrix1D y = gradF.like();
	SeqBlas.seqBlas.dgemv(false,1.0,Cinv,delta_m,0,y);
	kappa = 0.5 * SeqBlas.seqBlas.ddot(delta_m,y);
	// beta = -|delta_m.gradF|
	beta = -SeqBlas.seqBlas.ddot(delta_m,gradF);

	disc = Math.sqrt(Math.max(0, beta*beta - 4 * kappa * DELTA_CHISQ_TARGET));
	alpha = (beta - disc)/(2.0*((kappa>0) ? kappa : 1));

		
	if (!hasZeroTimeComponent)
	  {
		// For safety, compute also this crude bound: 
		// seek changes in times no bigger than FundamentalScalingConstants.BASE_SCALING_TIME
		// along the line 0<a<1, timeparams go like dt = a gradF
		// alpha_safe is then, roughly DELTA_CHISQ_SD/max(gradient)^2
		for (i=0; i<=_m.n_active_layers; i++)
		  {
			gradMax = Math.max(gradMax, Math.abs(gradF.getQuick(i)));
		  }
		alpha_safe = 0.1 * FundamentalScalingConstants.BASE_SCALING_TIME/((gradMax>0)? gradMax : 1.0);

		alpha = Math.min(2.0*alpha, alpha_safe);
	  }
	y=null;
	return alpha;

  }

  // Prior contribution for current fluid mode based on state in current_point
  double computePriorNegLogLikelihood(DoubleMatrix1D current_point)
  {
	int j;
	double chisq;

	DoubleMatrix1D mu_prior_only =  mbar_prior_only[fs];

	if (mu_prior_only.size() != current_point.size() )
	  {
		throw new RuntimeException("Bug when invoking computePriorNegLogLikelihood(). FIXME!!!!"); // 
	  }

 	DoubleMatrix1D delta_m = new DenseDoubleMatrix1D(mu_prior_only.size());

	for (j=0;j<mu_prior_only.size();j++) 
	  {
		delta_m.setQuick(j,current_point.getQuick(j)  - mu_prior_only.getQuick(j));
	  }	

	chisq = computePriorNegLogLikelihood_kernel(delta_m);
	delta_m = null; // GC

	return chisq;
  }

  /**
	 Prior contribution for current fluid mode based on state in Class variable property_sample
	 which must have been set correctly for this call to be meaningful
  */
  double computePriorNegLogLikelihood()
  {
	int j;
	double chisq;

	DoubleMatrix1D mu_prior_only =  mbar_prior_only[fs];

 	DoubleMatrix1D delta_m = new DenseDoubleMatrix1D(mu_prior_only.size());

	for (j=0;j<mu_prior_only.size();j++) 
	  {
		delta_m.setQuick(j,property_sample[j] - mu_prior_only.getQuick(j));
	  }
	chisq = computePriorNegLogLikelihood_kernel(delta_m);
	delta_m = null; // GC

	if (_SI.verbosity > 5) _SI.EH.println("Prior negLog= %10.4g ",chisq);

	return chisq;
  }

  // The prior evaluation kernel, returning (delta_m^T.C^-1.delta_m)/2
  private double computePriorNegLogLikelihood_kernel(DoubleMatrix1D delta_m)
  {
	double chisq;

 	DoubleMatrix1D ans = delta_m.like();

	// form C^1.delta_m -> ends up in-place in "ans"
	Cinv_prior_only[fs].zMult(delta_m,ans);

	// then delta_m^T.C^-1.delta_m/2.0
	chisq = delta_m.zDotProduct(ans)/2.0;

	ans = null; // GC

	return chisq;
  }


  /** 
	  If the inverse theory update using the seismic changes the time
	  ordering such that the updated mode has a different pinchout
	  configuration to the original expansion point, then the
	  quadratic approximation used to compute the update point no
	  longer has validity (we crossed a major step in the surface due
	  to the pinchout). The best point is that obtained from the new
	  quadratic form, but constrained to have the same pinchout
	  configuration as the original expansion point. This will lie
	  along the lines { t_{i} = t_{i+1} }, for indices i corresponding
	  to pinchout layers. This is a quadratic programming problem, and
	  this function returns the solution we need.

	  On entry, predicted distribution has approximate shape
	  ((m - m_new)^T.C^{-1}.(m - m_new)
	  where C = R_new^T. Rnew 

	  This update may also be operating in a reduced space of the t's if the prior
	  has synchronised some of the layers. Allow for this.
  */
  void updateModeUsingConstrainedQuadraticProgramming(LocalMode new_mode,DoubleMatrix1D expansion_point,
													  DoubleMatrix1D m_new, DoubleMatrix2D R_new)
  {
	boolean expansionStateIsPinchFree=false;
	int i,j,jj,k,N,N_k,all_layers_present_layer_code;
	int[] pinchout_indices;
	double[] pinchout_epsilon;

	/*
	  (1) Compute the desired pinchout configuration and the set {s_k}
	  of indices corresponding to pinched-out layers (i.e. solution
	  is within epsilon of t_{s_k} = t_{s_k+1}, for the set {s_k},
	  k=0...(N_k-1) 
	*/
	N_k=0;
	pinchout_indices = new int[_m.n_active_layers];
	pinchout_epsilon = new double[_m.n_active_layers];
	all_layers_present_layer_code =  (1<<_m.n_layers) - 1;
	// On entry new_mode.layer_code has the layer code for the state when the optimiser finished, i.e.
	// corresponding to the expansion point at first updating iteration
	if (_verbosity>3) _SI.EH.print("QP_constraints: t parameters (layer-indx): ");
	if (new_mode.layer_code == all_layers_present_layer_code)
	  {
		// Expansion point is in a non-pinched out state
		expansionStateIsPinchFree=true;
		// Current system configuration (_m.LP.t) is in the state
		// that produced some pinchouts. Let's find out what layers look to be pinching out
		for (i=1;i<_m.n_active_layers-1;i++) 
		  {
			// OLD
			// if (_m.LP[i+1].t.value == _m.LP[i].t.value)
			j = master_layer_number_sequence[i];
			jj = master_layer_number_sequence[i+1];
			if (_m.LP[jj].t.value <= _m.LP[j].t.value+TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER)
			  {
				pinchout_indices[N_k]=i;
				pinchout_epsilon[N_k]=DEPSILON;
				if (_verbosity>3) _SI.EH.printf("%d ",j);
				N_k++;  
			  }
		  }
	  }
	else
	  {
		// Expansion point is already in a pinched out state
		expansionStateIsPinchFree=false;
		// and update tries to force a different pinchout state
		// We want the pinchout_indices for the old state, so restore
		// times back to them & get pinchout indices
		for (i=0;i<=_m.n_active_layers-1;i++) setValuesFromIndex(i,expansion_point.getQuick(i));
		orderTimes();
		for (i=1;i<_m.n_active_layers-1;i++) 
		  {
			j = master_layer_number_sequence[i];
			jj = master_layer_number_sequence[i+1];
			if (_m.LP[jj].t.value == _m.LP[j].t.value)
			  {
				pinchout_indices[N_k]=i;
				pinchout_epsilon[N_k]=0.0;
				if (_verbosity>3) _SI.EH.printf("%d ",j);
				N_k++;  
			  }
		  }
	  }
	if (_verbosity>3) _SI.EH.printf("%n");

	/*
	  (2) Ok now we know the constraints are t[pinchout_indices[k]] = t[pinchout_indices[k]+1]
	  for k=0...N_k-1.

	  (2a) Form Linear Quadratic programming equations for Lagrange multiplier coefficients
	*/
	// But first, the original extrapolated mean will be useful, so copy it
	new_mode.quadratic_mean = m_new.like(); 
	new_mode.quadratic_mean.assign(m_new);
	if (!true)
	  {
		// Legacy method
		double sum=0.0;
		DoubleMatrix2D A = new DenseDoubleMatrix2D(N_k,N_k);
		DoubleMatrix1D rhs = new DenseDoubleMatrix1D(N_k);
		DoubleMatrix1D lambda;

		boolean RHSisNonZero=false;
		for (i=0;i<N_k;i++)
		  {
			for (j=0;j<=i;j++)
			  {
				sum = covarianceMatrixEntry(pinchout_indices[i],pinchout_indices[j],R_new);
				sum -= covarianceMatrixEntry(pinchout_indices[i],pinchout_indices[j]+1,R_new);
				sum -= covarianceMatrixEntry(pinchout_indices[i]+1,pinchout_indices[j],R_new);
				sum += covarianceMatrixEntry(pinchout_indices[i]+1,pinchout_indices[j]+1,R_new);
				A.setQuick(i,j,sum);
				if (i!=j) A.setQuick(j,i,sum);
			  }
			rhs.setQuick(i,m_new.getQuick(pinchout_indices[i]+1)-m_new.getQuick(pinchout_indices[i])+
					     pinchout_epsilon[i]);
			if (m_new.getQuick(pinchout_indices[i]+1) != m_new.getQuick(pinchout_indices[i])) RHSisNonZero=true;
		  }
		// (2b) Now solve them. 
		// now solve. This will set m_new to be on the quadratic programming boundary.
		if (RHSisNonZero)
		  {
			/* new */
			LUDecompositionQuick LUD = new LUDecompositionQuick();
			LUD.decompose(A);
			if (LUD.isNonsingular())
			  {
				LUD.solve(rhs);
				lambda=rhs;
			  }
			else
			  {
				_SI.EH.warn("Warning!: QP linear algebra system is singular! Fix this bug!" +
							"I'll set the lambda solution to zero and continue...");
				lambda = rhs;
				lambda.assign(0);
			  }
			// And recompute QP solution from formula in writeup
			for (i=0,k=0;i<m_new.size();i++)
			  {
				sum = m_new.getQuick(i);
				for (j=0;j<N_k;j++)
				  {
					sum += lambda.getQuick(j) * (covarianceMatrixEntry(i,pinchout_indices[j],R_new)
												 - covarianceMatrixEntry(i,pinchout_indices[j]+1,R_new));
				  }
				m_new.setQuick(i,sum);
			  }	
			// alg = null;
			LUD = null;
		  }
	  }
	else
	  {
		// Generic formulation for time pinchouts, grooves at NG=0,1 etc
		// Lagrange problem has form 
		//      m_QP = argmin { (m_QP -m_new)^T.C^{-1}.(m_QP-m_new) +(F.m_QP -m_c)^T.lambda }
		// lambda is a vector of Lagrange multipliers
		// F is a sparse matrix that embeds the pinchout bits, plus other linear constraints
		// (1) for pinchouts, e.g. F has rows with  {1,-1} adjacent pairs, m_c=0
		// (2) for NG grooves, F has 1's and m_c=0,1
		int Nc = N_k;          // number of constraints
		N = m_new.size();  // dimension of system

		// NG detection stuff here
		int N_NG=0;
		int[] NG_sat_indices = new int[_m.n_active_layers];
		Integer I;
		if (_verbosity>3) _SI.EH.print("QP_constraints: NG parameters (layer,indx): ");
		for (i=0,k=0;i<_m.n_active_layers;i++) 
		  {
			LayerProperties.Property prop = _m.LP[i].NG;	
			if (prop!=null)
			  {
				if ((I = (Integer)map_to.get(prop)) != null)
				  {
					j = I.intValue();
					double NG = expansion_point.getQuick(j);
					double NGnew = m_new.getQuick(j);
					if (((NG<=NG_EPSILON) && (NGnew<=NG_EPSILON))
						|| ((NG >= (1.0-NG_EPSILON)) && (NGnew >= (1.0-NG_EPSILON))))
					  {
						// code NG=0 as -j, NG=1 as +j (We know j is *never* 0)
						NG_sat_indices[N_NG] = (NG<NG_EPSILON) ? -j: j; 
						if (_verbosity>3) _SI.EH.printf("(%d,%d) ", i,j);
						N_NG++;
					  }
				  }
			  }
		  }
		Nc+=N_NG;  // total num constraints
		if (_verbosity>3) _SI.EH.printf("%n");

		// other constraints here....

		DoubleMatrix1D m_c = new DenseDoubleMatrix1D(Nc);  // constraint RHS
		DoubleMatrix2D F = new SparseDoubleMatrix2D(Nc,N);

		DoubleMatrix2D RFt = new DenseDoubleMatrix2D(N,Nc);
		DoubleMatrix2D RtRFt = new DenseDoubleMatrix2D(N,Nc);
		DoubleMatrix2D C_c = new DenseDoubleMatrix2D(Nc,Nc);

		// Fill F and mu_c matrices for set of linear constraints:
		// (1) Pinchout entities
		for (i=0;i<N_k;i++)
		  {
			F.setQuick(i,pinchout_indices[i],1.0);
			F.setQuick(i,pinchout_indices[i]+1,-1.0);
			m_c.setQuick(i,0.0);
		  }
		// (2) NG grooves: attempt to deal with NG=0 or 1 saturated optima (V-shapes)
		for (i=N_k,j=0;i<Nc;i++,j++)
		  {
			F.setQuick(i,Math.abs(NG_sat_indices[j]),1.0);
			// m_c.setQuick(i,(NG_sat_indices[j] > 0)? 1.0 : 0.0);
			m_c.setQuick(i,(NG_sat_indices[j] > 0)? 1.0-NG_EPSILON : NG_EPSILON);
		  }


		// Get QP solution
		// form m_c <- F.m_new -m_c
		SeqBlas.seqBlas.dgemv(false,1.0,F,m_new,-1.0,m_c);
		// form RFt = R_new.F^T 
		SeqBlas.seqBlas.dgemm(false,true,1.0,R_new,F,0,RFt);
		// then  RtRFt = C.F^T = R_new^T.R_new.F^T = R_new^T . RFt 
		SeqBlas.seqBlas.dgemm(true,false,1.0,R_new,RFt,0,RtRFt);
		// then C_c = F.R_new^T.R_new.F^T = F.RtRFt
		SeqBlas.seqBlas.dgemm(false,false,1.0,F,RtRFt,0,C_c);

		// Now Solve C_c.rhs = m_c, 
		// i.e. rhs = C_c^{-1}.m_c = {F.C.F^T}^{-1}. (F.m_new -m_c) in original quantities
		LUDecompositionQuick LUD = new LUDecompositionQuick();
		LUD.decompose(C_c);
		if (LUD.isNonsingular())
		  {
			LUD.solve(m_c);
		  }
		else
		  {
			_SI.EH.warn("Warning!: Generalised QP linear algebra system is singular! Fix this bug!" +
						"I'll set the lambda solution to zero and continue...");
			m_c.assign(0);
		  }
		// Final adjustement to new QP solution
		// form m_c <- F.m_new -m_c
		SeqBlas.seqBlas.dgemv(false,-1.0,RtRFt,m_c,1.0,m_new);
		// QP solution is now in m_new

	  }

	// else: (RHS==0) degenerate solution. all lambda's=0 => m_new unchanged
	if (expansionStateIsPinchFree)
	  {
		// enforce time ordering in case roundoff error didn't quite produce it
		for (i=1;i<_m.n_active_layers-1;i++) 
		  {
			if (m_new.getQuick(i) >= m_new.getQuick(i+1))
			  {
				m_new.setQuick(i,m_new.getQuick(i+1)-DEPSILON);
			  }
		  }
	  }
	else
	  {
		// make sure times sit exactly on constraints prescribed
		for (i=0;i<N_k;i++)
		  {
			if (m_new.getQuick(pinchout_indices[i]) != m_new.getQuick(pinchout_indices[i]+1))
			  {
				m_new.setQuick(pinchout_indices[i]+1,
							   m_new.getQuick(pinchout_indices[i]));
			  }
		  }
	  }

	/*
	  FIXME....should apply a test here to see if quadratic_mean is too far from QP QP_optimum to be able to generate viable
	  candidates with the correct pinchout configuration or not......
	  // Example: Test whether QP mode is not too far from quadratic point
	  // Get quantity ((m_Q - m_{QP})^T.C^{-1}.(m_Q - m_{QP}) and see if it's outside N(n,sqrt(n))
	  double chisq=0.0;
	  double y[];
	  N = m_new.size();
	  y = new double[N];
	  for (i=0;i < N; i++)
	  {
	  for (sum=0.0, j=0; j<=i-1 ; j++) sum += R_new.getQuick(j,i) * y[j];
	  y[i] = ((m_new.getQuick(i) - new_mode.quadratic_mean.getQuick(i)) - sum)/R_new.getQuick(i,i);
	  chisq += y[i]*y[i];
	  }
	  _SI.EH.println("chisq-test = "+chisq);
	  y = null;
	*/
	// Simple feedback to reflect above worries
	if (_verbosity>3) 
	  {
		// _SI.EH.println("QP DEBUG CNT: "+debugCntr);
		_SI.EH.print("QP current-point times: ");
		for (i=1;i<_m.n_active_layers-1;i++) _SI.EH.printf("%10.7g ",expansion_point.getQuick(i));
		_SI.EH.println("");
		_SI.EH.print("QP     new-point times: ");
		for (i=1;i<_m.n_active_layers-1;i++) _SI.EH.printf("%10.7g ",m_new.getQuick(i));
		_SI.EH.println("");
		_SI.EH.print("QP time-separations from Quadratic optimum (std-dev-normalised):  ");
		for (i=1;i<_m.n_active_layers-1;i++) 
		  {
			_SI.EH.printf("%10.7g ",((m_new.getQuick(i) - new_mode.quadratic_mean.getQuick(i))/
											Math.sqrt(covarianceMatrixEntry(i,i,R_new))));
		  }
		_SI.EH.println("");
	  }
	// clean up
	pinchout_indices=null;
	//	A = null;
	//	rhs=null;
	// lambda=null;
  }

  /**
	 Returns true if adjustments made
  */
  void experimentalProjectSDdirectionOntoQuadraticProgrammingConstraint(LocalMode new_mode,DoubleMatrix1D expansion_point,
																		   DoubleMatrix1D m_new)
  {
	int i,j,jj,Ng,N,N_k;
	byte[] pinchoutIndicator;
	double delta;

	/*
	  (1) Compute the desired pinchout configuration and the set {s_k}
	  of indices corresponding to pinched-out layers (i.e. solution
	  is within epsilon of t_{s_k} = t_{s_k+1}, for the set {s_k},
	  k=0...(N_k-1) 
	*/
	N_k=0;

	pinchoutIndicator = new byte[_m.n_active_layers];
	// On entry new_mode.layer_code has the layer code for the state when the optimiser finished, i.e.
	// corresponding to expansion_point at first updating iteration
	for (i=1;i<_m.n_active_layers-1;i++) 
	  {
		j = master_layer_number_sequence[i];
		jj = master_layer_number_sequence[i+1];
		if (_m.LP[jj].t.value <= _m.LP[j].t.value+TIME_EPSILON_DETECT_PINCHOUT_IN_MINIMISER)
		  {
			pinchoutIndicator[j]=1;         // within "epsilon"
			N_k++;  
			if (_m.LP[jj].t.value <= _m.LP[j].t.value) pinchoutIndicator[j] = 2; // actually strictly pinched 
		  }
	  }
	  
	/*
	  (2) Ok now we know the constraints are t[pinchout_indices[k]] = t[pinchout_indices[k]+1]
	  for k=0...N_k-1.

	  Here's the recipe for the nearest point to m_tilde on the constraint subspace:
	  Components that are pinch-free are preserved.
	  Components in a pinch group get assigned the average of the change for that group
	*/	
	for (i=1;i<_m.n_active_layers-1;i++) 
	  {
		if (true)
		  {
			Ng=1;
			for (j=i;j<_m.n_active_layers-1;j++) 
			  {
				if (pinchoutIndicator[j] > 0 ) Ng++; 
				else break;
			  }
			// group to average now spans i....(j-1) inclusive
			if (Ng>1)
			  {
				// get average
				delta=0;
				for (jj=i;jj<=j;jj++) delta += m_new.getQuick(jj) - expansion_point.getQuick(jj);
				delta/=Ng;
				// and assign
				for (jj=i;jj<=j;jj++) 
				  {
					m_new.setQuick(jj, expansion_point.getQuick(jj) + delta);
				  }
				// now start next one after group
				i=j-1;
			  }
		  }
		else
		  {
			_m._SI.EH.println("Shutting off time changes");
			// Test! Completely switch off all time-changes
			m_new.setQuick(i, expansion_point.getQuick(i));	
		  }
	  }

	if (false)
	  {
		// Experimental attempt to reduce scale-prejudice to times....
		double SCALING=1.0;
		LayerProperties.Property prop;
		Integer I;
		for (i=0;i<_m.n_active_layers;i++) 
		  {
			prop = _m.LP[i].NG;	
			if (prop!=null)
			  {
				if ((I = (Integer)map_to.get(prop)) != null)
				  {
					j = I.intValue();
					delta = m_new.getQuick(j) - expansion_point.getQuick(j);
					m_new.setQuick(j,expansion_point.getQuick(j) + SCALING * delta);
				  }
			  }
		  }
		if ((prop = _m.LP[_m.n_layers-1].t_base) != null)
		  {
			if ((I = (Integer)map_to.get(prop)) != null)
			  {
				j = I.intValue();
				delta = m_new.getQuick(j) - expansion_point.getQuick(j);
				m_new.setQuick(j,expansion_point.getQuick(j) + SCALING * delta);
			  }
		  }
	  }



	// Experimental attempt to deal with NG=0 or 1 saturated optima (V-shapes)
	if (false)
	  {
		LayerProperties.Property prop;
		Integer I;
		for (i=0;i<_m.n_active_layers;i++) 
		  {
			prop = _m.LP[i].NG;	
			if (prop!=null)
			  {
				double NG = prop.value;
				if ((NG<=NG_EPSILON) || (NG >=(1.0-NG_EPSILON)))
				  {
					if ((I = (Integer)map_to.get(prop)) != null)
					  {
						_m._SI.EH.println("LOCKING NG...layer indx "+(i));
						j = I.intValue();
						// lock to descent direction onto this value already
						m_new.setQuick(j,(NG<=NG_EPSILON) ? 0 : 1.0);
					  }
				  }
			  }
		  }
	  }

	return;
  }

  void debugShowDeltaMComponents(LocalMode new_mode,DoubleMatrix1D expansion_point,
								 DoubleMatrix1D m_new)
  {
	int i,j;

	
	for (i=0;i<_m.n_active_layers;i++) 
	  {
		LayerProperties.Property prop = _m.LP[i].t;	
		if (prop!=null)
		  {
			double t = prop.value;
			Integer I;
			if ((I = (Integer)map_to.get(prop)) != null)
			  {
				j = I.intValue();
				_m._SI.EH.println("t component "+i+" delta = "+(m_new.getQuick(j) - expansion_point.getQuick(j)));
			  }
		  }
	  }
	for (i=0;i<_m.n_active_layers;i++) 
	  {
		LayerProperties.Property prop = _m.LP[i].NG;	
		if (prop!=null)
		  {
			double NG = prop.value;
			Integer I;
			if ((I = (Integer)map_to.get(prop)) != null)
			  {
				j = I.intValue();
				_m._SI.EH.println("NG component "+i+" delta = "+(m_new.getQuick(j) - expansion_point.getQuick(j)));
			  }
		  }
	  }
  }

  /*
	Get the entry C_{p,q} in the covariance matrix of the mode new_mode, where
	C has the decomposition C = R^T . R
  */
  private double covarianceMatrixEntry(int p, int q, DoubleMatrix2D R)
  {
	double sum=0.0;
	int i,pp,qq;
	// since C is symmetric, there's an efficient way to get this
	if (p<=q)
	  {
		pp=p;
		qq=q;
	  }
	else
	  {
		pp=q;
		qq=p;
	  }
	for (i=0;i<=pp;i++)
	  {
		sum += R.getQuick(i,pp) * R.getQuick(i,qq);
	  }
	return sum;
  }

  /**
	 Draw samples of the model using the MCMC technique, choosing algorithm suitably
  */
  void sampleFromLocalModes() 
  {

	// convergence diagnostics 
	if (_SI.performRafteryAutoDecimation || _SI.performGelmanConvergenceDiagnostics ||
		_SI.performRafteryConvergenceDiagnostics) 
	  {
		// Make a new convergence monitor
		// if its not Gelman, it must be Raftery....
		// (performGelmanConvergenceDiagnostics==false) => Raftery
		CM = new MCMCConvergenceMonitor(_SI.EH,_verbosity,
										_SI.performGelmanConvergenceDiagnostics,
										_SI.performRafteryConvergenceDiagnostics);
		if (_SI.performRafteryAutoDecimation) CM.setRafteryCalibrationRun(true);
	  }

	// single or double loop over samplers, depending on whether auto-decimation is demanded
	// first run is calibration, second the actual output
	for (int run=0; run<2; run++)
	  {
		if (_SI.UseRandomWalkSampler)
		  {
			// Metropolis random walker
			randomWalkSampler = new RandomWalkSampler(_m);
			if (_SI.performGelmanConvergenceDiagnostics ||
				_SI.performRafteryConvergenceDiagnostics) 
			  {
				// this call wraps randomWalkSampler.sampleFromLocalModesUsingRandomWalkSampler
				// with extra functionality to monitor convergence
				performConvergenceDiagnosticsWhileSamplingFromLocalModes(CM);
			  }
			else
			  {
				randomWalkSampler.sampleFromLocalModesUsingRandomWalkSampler();
			  }
		  }
		else
		  {
			// Independence sampler based on full mode approximations
			independenceSampler = new IndependenceSampler(_m,_SI.useApproximateMixedGaussianIndependenceSampler);
			if (_SI.performGelmanConvergenceDiagnostics ||
				_SI.performRafteryConvergenceDiagnostics) 
			  {
				// this call is a loop around independenceSampler.sampleFromLocalModesUsingRandomWalkSampler
				// with extra functionality to monitor convergence
				performConvergenceDiagnosticsWhileSamplingFromLocalModes(CM);
			  }
			else
			  {
				independenceSampler.sampleFromLocalModesUsingFullLocalModes();
			  }
		  }
		// second pass, decimation statistics available. Switch on writing (i.e. calibration off)
		if (_SI.performRafteryAutoDecimation) 
		  {
			CM.setRafteryCalibrationRun(false);
		  }
		else break;  // only single pass needed if we are not doing performRafteryAutoDecimation
	  }
  } 

  /**
	 Perform Convergence diagnostics in desired style, either Gelman or Raftery, as
	 per the state variables set in the constructor.
  */
  void performConvergenceDiagnosticsWhileSamplingFromLocalModes(MCMCConvergenceMonitor CM)
  {
	if (CM.GelmanStyleDiagnostics)
	  {
		// NB Raftery style may be done as well if CM.RafteryStyleDiagnostics==true
		CM.initialiseGelmanStyleConvergenceDiagnostics(_m.n_active_layers-1,_SI.n_samples_per_trace_undecimated);
		int ensemble_cntr;
		// run over ensemble of walks
		for (ensemble_cntr=0; ensemble_cntr<CM.N_INDEPENDENT_CHAINS; ensemble_cntr++)
		  {
			if (_SI.UseRandomWalkSampler)
			  {
				randomWalkSampler.sampleFromLocalModesUsingRandomWalkSampler();
			  }
			else
			  {
				independenceSampler.sampleFromLocalModesUsingFullLocalModes();
			  }
		  }
		CM.finishGelmanStyleConvergenceDiagnostics(_SI.n_samples_per_trace_undecimated);
	  }
	else if (CM.RafteryStyleDiagnostics)
	  {
		CM.initialiseRafteryStyleConvergenceDiagnostics(_m.n_active_layers-1,_SI.n_samples_per_trace_undecimated);
		// run over ensemble of walks
		if (_SI.UseRandomWalkSampler)
		  {
			randomWalkSampler.sampleFromLocalModesUsingRandomWalkSampler();
		  }
		else
		  {
			independenceSampler.sampleFromLocalModesUsingFullLocalModes();
		  }
		CM.finishRafteryStyleConvergenceDiagnostics();
	  }
  }

  void aquireScalarEstimands(int realisation_number)
  {
	int k,j,time_indx;
	// variables in a rejected realisation. The trace has the true state variables.
	// Index of INTERNAL layer times in output trace...
	time_indx = _SI.MF.reducedPropertyList.indexOf("t") * _m.n_layers;
	for (k=1; k<_m.n_active_layers; k++)
	  {
		j = master_layer_number_sequence[k];
		CM.setEstimand(k-1,(double)_SI.tr_model.get_data_value(time_indx+j));
	  }
	CM.completeScalarEstimandsUpdate(realisation_number);
  }


  /**
	 From a full set of local modes, picks one with probability proportional
	 to the mode.pdf value, which is set from the *posterior* weights. 
	 Returns the index of the mode selected.
	 (using uniform deviate and find where it slots in the CDF array... )
  */
  int proposeModeFromApproximatePosterior()
  {
	int k;
	double cdf;

	cdf = Ugen.nextDouble();
	// find associated local mode
	for (k=0;k<localModes.size();k++)
	  {
		mode = localModes.get(k);
		if (cdf<=mode.CDF) break;
	  }
	// set correct fluid state and working variables accordingly
	fs = mode.fs;
	setMapsFromFluidState(fs);
	return k;
  }



  /**
	 From a full set of local modes, picks one with probability proportional
	 to the mode.pdf value. Returns the index of the mode selected.
	 (using uniform deviate and find where it slots in the CDF array... )
	 If modeJumpsMustLeaveOldMode is set, will force random pick of new mode if possible

	 If useEstimatedPosteriorWeight==true, mode selection is based on cdf constructed from estimated 
	 posterior weights *only*. If false, mode selection based on equal probablity of
	 using posterior weight and prior weight, i.e. weight = 0.5(posterior_weight + prior_weight)
  */
  int proposeNewMode(boolean useEstimatedPosteriorWeight, 
					 boolean proposedModeMustLeaveOldMode, LocalMode oldMode)
  {
	int k,size;
	double sum_prior=0, sum_post=0, sum=0, cdf,weight;

	/* form cdf totals maybe omitting current mode */

	size = localModes.size();
	for (k=0,sum=0.0; k<size; k++)
	  {
		mode = localModes.get(k);
		if (mode==oldMode && proposedModeMustLeaveOldMode) continue;
		sum += (useEstimatedPosteriorWeight ? mode.weight : mode.hybrid_weight);
	  }
	// draw random deviate
	cdf = sum * Ugen.nextDouble();
	for (k=0,sum=0.0; k<size; k++)
	  {
		mode = localModes.get(k);
		if (mode==oldMode && proposedModeMustLeaveOldMode) continue;
		sum += (useEstimatedPosteriorWeight ? mode.weight : mode.hybrid_weight);
		if (cdf<=sum) break;
	  }
	
	// set correct fluid state and working variables accordingly
	fs = mode.fs;
	setMapsFromFluidState(fs);
	return k;
  }

  /**
	 From a full set of local modes, picks one with probability proportional
	 to the mode.prior_pdf value. Returns the index of the mode selected.
	 (using uniform deviate and find where it slots in the CDF array... )
  */
  int proposeModeFromPrior()
  {
	int k;
	double cdf;

	cdf = Ugen.nextDouble();
	// find associated local mode
	for (k=0;k<localModes.size();k++)
	  {
		mode = localModes.get(k);
		if (cdf<=mode.prior_CDF) break;
	  }
	// set correct fluid state and working variables accordingly
	fs = mode.fs;
	setMapsFromFluidState(fs);
	return k;
  }

		

  /**
	 Diagnostic tool for when sampler gets into a state where most samples seem unphysical
  */
  void showDiagnosticsOfEffectivePropertiesForAllLayers()
  {
	_SI.EH.println("Persistent fluid substitution errors. Some diagnostics...");

	_m.computeEffectivePropertiesForAllLayers(false,false,false,true,false);

// 	for (int k=0;k<_m.n_layers;k++)
// 	  {
// 		// Some MC combinations might just be unphysical
// 		_SI.EH.println("Layer ... "+ k);
// 		if (!_m.LP[k].computeEffectiveProperties(fluid_state.types[k],false))
// 		  {
// 			// Failed. Show some diagnostics.
// 			_SI.EH.println("Fluid substitution error: diagnostics...");
// 			_m.LP[k].computeEffectiveProperties(fluid_state.types[k],true);
// 		  }
// 	  }
  }

  /**
	 Useful method for gleaning information about Gassman substiution problems
  */
  void showGassmanDiagnostics(boolean sample_succeeded, boolean setUp)
  {
	boolean success;
	// here we print out a table of all the relevant Gassman substituions
	// and can use this to diagnose problems
	if (setUp)
	  {
		_SI.EH.println("%nGassman_test: layer-indx phi_s vp_s vs_s K_satb rho_satb%n");
	  }
	else
	  {
		_m.computeEffectivePropertiesForAllLayers(false,false,false,false,true);
//      for (int k=0;k<_m.n_layers;k++)
// 		  {
// 			// Some MC combinations might just be unphysical
// 			success = _m.LP[k].computeEffectiveProperties(fluid_state.types[k],false);
// 			_SI.EH.println("Gassman_test: %d %10.4f %10.4f %10.4f %s",k, _m.LP[k].phi_s.value, 
// 							  _m.LP[k].vp_s.value,
// 							  _m.LP[k].vs_s.value,
// 							  (success?"ok":"bad"));
// 		  }
	  } 
  }



  double computeIsopachNegLogLikelihood(boolean force_calculation)
  {	
	int i,k;
	double chisq=0.0, negLogLik=0.0;

	i=0;  // runs along isopach error vector isopach_errors[]
	logDetCisopachAnd2PiPieces = 0.0;

	if (force_calculation || !(_m.MPD.n_isopachs==0 || _SI.ignoreIsopachs))
	  {
		for (k=0;k<_m.n_layers;k++)
		  {
			ModelPriorDescriptions.MetaInfo MI = _m.MPD.MetaInfo[k];
			LayerProperties LPI = _m.LP[k];
			if (MI.dz!=null && !(MI.dz.value==0 && MI.sigma_dz.value==0)) // layer has an non-trivial isopach
			  {
				// layer time thickness...
				double t_base = (k==_m.n_layers-1) ? LPI.t_base.value : _m.LP[k+1].t.value;
                // and => thickness. Factor of 1/2 since times are two-way
				double DZ = LPI.vp_eff * 0.5 * (t_base - LPI.t.value);
				// double sdz = (MI.sigma_dz.value > 0) ? MI.sigma_dz.value : 0.5; // FIXME remove
				isopach_errors[i] = (DZ-MI.dz.value)/MI.sigma_dz.value;
				chisq += isopach_errors[i] * isopach_errors[i];
				
				// make efficient somehow....
				logDetCisopachAnd2PiPieces += MI.log_sigma_dz; // = Math.log(MI.sigma_dz.value);
				logDetCisopachAnd2PiPieces += MathUtil.HALF_LOG_TWOPI;

				i++;
			  }
		  }
		negLogLik += chisq/2.0;
	  }
	if (_SI.verbosity>5) _SI.EH.println("Isopach negLog= %10.4g ",negLogLik);

	return negLogLik;
  }

  /** 
	  Likelihood function for imposing depth self consistency. Computed
	  after times are truncated and ordered, so depths are non-negative.
	  LayerProperty effective properties must be current and valid.
  */

  double computeConsistentDepthNegLogLikelihood(boolean showErrors)
  {	
	int i;
	double chisq=0.0, negLogLik=0.0;

	if (applyDepthConsistency)
	  {
		for (i=0;i<_m.n_layers-1;i++)
		  {
			LayerProperties LPI = _m.LP[i];
			LayerProperties LPIplus1 = _m.LP[i+1];
			// layer time thickness...
			// double t_base = (i==_m.n_layers-1) ? LPI.t_base.value : LPIplus1.t.value;
			double t_base = LPIplus1.t.value;
			// and => thickness. Factor of 1/2 since times are two-way
			double DZ = 0.5 * LPI.vp_eff * (t_base - LPI.t.value);
			consistent_depth_errors[i] = (LPIplus1.d.value - (LPI.d.value+DZ))/sigma_consistent_depth;
			chisq += consistent_depth_errors[i] * consistent_depth_errors[i];
			if (showErrors) _SI.EH.printf("%10.4g ",consistent_depth_errors[i]);
		  }
		negLogLik += chisq/2.0;

		// if (_SI.verbosity>4) _SI.EH.println("Depth self-consistency negLog= %10.4g ",negLogLik);
	  }
	return negLogLik;
  }



  void setValuesFromIndex(int i,double value)
  {
	int size;
	LayerProperties.Property prop;

	ArrayList list =(ArrayList) map_from.get(i);

	size = list.size();
	for (int j=0;j<size;j++)
	  {
		prop = (LayerProperties.Property)list.get(j);
		if (prop==null)
		  {
			throw new RuntimeException("Error finding property for index "+i+", fluid combo "+fs);
		  }
		// and set its value suitably
		prop.setPosValue(value);
	  }
  }

  /** 
	 Simple class used to keep key identifiers of local minima found in the
	 likelihood surface, used to prevent incorporating duplicate
	 modes into the list of local modes
  */
  static class ModeIdentifier
  {
	/**
	   fluid state
	*/
	int fs;
	/**
	   layers-present code-identifier
	*/
	int layer_code;
	/**
	   suite of parameters at minima
	*/
	double[] timesPlusOtherParameters;

	/**
	   maximum of mode at minima
	*/
	double chisq;
	
	/**
	   Ref. to full mode associated with the successful update of this mode
	*/
	LocalMode mode=null;

	public  ModeIdentifier(int fs, int layer_code, double[] timesPlusOtherParameters, double chisq)
	{
	  this.fs = fs;
	  this.layer_code = layer_code;
	  this.timesPlusOtherParameters = (double [])timesPlusOtherParameters.clone();
	  this.chisq = chisq;
	}
  }







  /*************************************/   

  // Useful debugging functions

  /************************************/

  /** 
	  A method to examine the posterior surface near the ML (maximum
	  likelihood) points, by dumping files with Xsections for selected
	  properties, with all other properties fixed at their ML
	  values. The surface generally looks smooth and quadratic, but
	  definite steps are evident where pinchouts come into play, plus
	  kinks at time pinchouts.  Xsections for time properties are
	  shown for 
	  (a) if scaled==false:<br> 
	  +/- N*DT_STEP/2 time range about the ML point, and for
	  non-time properties, for +/- 20% about the ML point,<br>
	  (b) if scaled==true:<br>
	  over a range (-2<chi<2) where logPost ~ chi^2/2, so if mode
	  is truly quadratic, parabolas with variation of +2 above ML point 
	  should be evident.<br> Xsections
	  are writtenas as ascii files, with the supplied prefix, viz
	  prefix.property_name.layer_number.txt. Must be called with a
	  completely defined mode, where newmode.QP_optimum is expected to
	  be the ML point.
  */
  private void debugPrintXsections(String prefix, LocalMode newmode)
  {
	int i,ii,j,jj,k,kk,N;
	int NSTEPS=50;
	int n;
	Integer I;
	double DT_STEP = 0.002;
	double NONTIME_FRACTIONAL_RANGE = 0.2;
	boolean scaled,good_state=true;
	double val, newval, ans, scale=1;
	double[] fvals;
	double[] xvals;
	ModelPriorDescriptions.MetaInfo MI;
	HashMap used_indices=new HashMap();
	DoubleMatrix1D centre_point=null;
	String name;

	N=2*NSTEPS;
	fvals=new double[N+1];
	xvals=new double[N+1];

	// Working array
	ArrayList<LayerProperties.Property> prop_arrayl = new ArrayList<LayerProperties.Property>();

	// point everything at mode we set.
	setMode(newmode);
	// central point through which to draw Xsections
	centre_point = mode.QP_optimum;

	// set scaled=true if you want scaled Xsections for all variables (-2<chisq<2)
	scaled = true;

	if (scaled)
	  {
		NONTIME_FRACTIONAL_RANGE = _SI.XsectionsDeviation; // user option....typically 2
	  }
	// _SI.verbosity=6;

	n = centre_point.size();
	// property_sample = new double[n];
	for (k=0;k<n;k++)
	  {
		val = centre_point.getQuick(k);
		setValuesFromIndex(k,val);
		property_sample[k]=val; // needed for prior
	  }
	orderTimes();
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	


	// loop over properties that are interesting to draw graphs of
	for (i=0;i<_m.n_layers;i++)
	  {
		MI = _m.MPD.MetaInfo[i];
		LPI = _m.LP[i];

		// Here's a list of all the things that you *might* 
		// want to see Xsections of
		/*
		LayerProperties.Property[] prop_array =
		  { LPI.d, LPI.LFIV, LPI.t, (i==_m.n_layers-1) ? LPI.t_base : null,
			LPI.rho_b, LPI.vp_b,
			LPI.rho_oil, LPI.rho_lsg, LPI.rho_gas,
			LPI.vp_oil, LPI.vp_lsg, LPI.vp_gas,
			LPI.S_oil, LPI.S_lsg, LPI.S_gas,
			LPI.NG, LPI.vp_s, LPI.vs_s, LPI.phi_s,
			LPI.vp_m, LPI.vs_m, LPI.rho_m}; */
		// But we'll do this subset, plus systematics below......
		LayerProperties.Property[] prop_array =
		  { LPI.t, (i==_m.n_layers-1) ? LPI.t_base : null,
			LPI.NG,LPI.vp_s,LPI.vs_s,LPI.vp_m,LPI.vs_m,LPI.phi_s,LPI.rho_m,
			((_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL) ? 
			 LPI.float_fraction : null),
			((i==0) ?_m.SM.A : null), 
			(((i==0) && !_m.SM.A_B_avo_synchronised) ?_m.SM.B : null)};

		// form union of properties master list via ArrayList methods...
		prop_arrayl.clear();
		prop_arrayl.addAll((List<LayerProperties.Property>)Arrays.asList(prop_array));
		if (_SI.MF.ENABLE_HARD_ROCK_MODEL && _m.needHardRockTerms)
		  {
			prop_arrayl.addAll((List<LayerProperties.Property>)
							   Arrays.asList(LPI.NG_hard, LPI.vp_hard, LPI.vs_hard, LPI.rho_hard));
		  }
		if (_m._SI.doSimpleFluidUnitsAndContacts) prop_arrayl.addAll(_m.fluidUnitsAndContacts.fc_props);

		for (k=0; k<prop_arrayl.size(); k++)
		  {
			LayerProperties.Property prop = prop_arrayl.get(k);
			if (prop!=null && (I = (Integer)map_to.get(prop))!=null)
			  {
				// if the property in question is a mirror, skip it (we've already done it)
				if (used_indices.get(I) != null) continue;
				else used_indices.put(I,I);
				
				// This property is indexed, so is stochastic.
				// get the column index
				int indx = I.intValue();
				// Skip any column indices out of range (shouldn't happen)
				if (indx>=n) continue;

				// we need the name of this property
				name = getPropertyName(prop);
				// get the base value of the property
				val = centre_point.getQuick(indx);
				if (scaled) scale = 1.0/Math.sqrt(mode.CtildeInvDiag[indx]);				
				boolean zero_mean = (val==0.0);
				// run over range for val. NB Step (ii == NSTEPS+1) is to restore value
				for (j=0,ii=-NSTEPS;ii<=NSTEPS+1;ii++,j++)
				  {
					if (ii==NSTEPS+1) 
					  {
						newval = val; // restore on last loop pass
					  }
					else
					  {
						if (prop==LPI.t || prop==LPI.t_base)
						  { 
							// times
							if (scaled) 
							  newval = val+scale*NONTIME_FRACTIONAL_RANGE*((double)ii/(double)NSTEPS);
							else
							  newval = val*(1+DT_STEP*((double)ii/(double)NSTEPS));
						  }
						else
						  {
							// non-time properties
							if (!zero_mean)
							  {
								if (scaled) 
								  newval = val+ scale*NONTIME_FRACTIONAL_RANGE*((double)ii/(double)NSTEPS);
								else
								  newval = val*(1+NONTIME_FRACTIONAL_RANGE*((double)ii/(double)NSTEPS));
							  }
							else
							  {
								if (scaled) 
								  newval = scale*NONTIME_FRACTIONAL_RANGE*((double)ii/(double)NSTEPS);
								else
								  newval = NONTIME_FRACTIONAL_RANGE*((double)ii/(double)NSTEPS);
							  }
						  }
					  }
					property_sample[indx]=newval;     // must set this for prior evaluation
					setValuesFromIndex(indx,newval);

					// set times
					if (prop==LPI.t || prop==LPI.t_base)
					  {
						// (A) set all times cleanly from ML point
						for (jj=0;jj<=_m.n_layers;jj++) 
						  {
							if (jj==indx) continue; // skip this time
							setValuesFromIndex(jj,centre_point.getQuick(jj));
							property_sample[jj]=centre_point.getQuick(jj);
						  }
						// set new time point
						setValuesFromIndex(indx,newval);
						property_sample[indx]=newval;     // must set this for prior evaluation

						// must get times in right order
						orderTimes();
					  }

					if (_SI.doSimpleFluidUnitsAndContacts || (prop!=LPI.t && prop!=LPI.t_base))
					  {
						// Recompute the effective properties for all layers
						// Required if doing contacts, but if not, is not needed for times
						good_state = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false);
					  }
					else good_state=true;
					
					if (good_state)
					  {
						ans = debugComputeNegLogPosterior();
					  }
					else ans=0.0;
					// store result
					if (ii != NSTEPS+1) 
					  {
						fvals[j] = ans;
						xvals[j] = newval;
						// Amusing hack to show profiles for depth...
						// if (prop==LPI.t || prop==LPI.t_base)
						//   {
						//    xvals[j] = (prop==LPI.t) ? LPI.suspended_depth_top : LPI.suspended_depth_bottom;
						//   }
					  }
					// NB Restoring happens at step (ii == NSTEPS+1)
				  }
				String key = new String(prefix+"."+name+"."+i+".txt");
				Util.printfileZO(key,xvals,fvals,fvals.length);
			  }
		  }
	  }
	// GC
	// property_sample=null;
  }


  /** 
	  A method to examine the posterior surface near the ML (maximum
	  likelihood) points and produce a set of directions which look suspicious, and
	  likely candidates for "big jumps" in the RWS sampler.

	  We do this by exploring points along the pricipal axis of the reqquested mode's
	  quadratic approximation, at positions 
	  -N_STD_DEVIATIONS,-N_STD_DEVIATIONS+1....N_STD_DEVIATIONS
	  along each pricipal axis.


	  Must be called with a completely defined mode, where
	  newmode.QP_optimum is expected to be the ML point.
  */
  private void computeBigJumpDirectionVectorsFromPrincipalComponents(LocalMode newmode)
  {
	int N_STD_DEVIATIONS = 4;
	double THRESHOLD = 0.7;
	int i,ii,j,jj,k,N,nXsectionSteps;
	int n,iVal;
	boolean good_state=true;
	double val, newval, ans, mode_ans, scale,lambda;
	// double[] fvals;
	// double[] xvals;
	DoubleMatrix1D centre_point=null;
	String name=null;
	ArrayList directions = new ArrayList();

	nXsectionSteps = N_STD_DEVIATIONS;
	N = 2*nXsectionSteps;

	// point everything at mode we set.
	setMode(newmode);
	// central point through which to draw Xsections
	centre_point = mode.QP_optimum;
	n = centre_point.size();

	// set centre point and evaluate
	for (k=0;k<n;k++)
	  {
		val = centre_point.getQuick(k);
		setValuesFromIndex(k,val);
		property_sample[k]=val; // needed for prior
	  }
	orderTimes();
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	
	mode_ans = debugComputeNegLogPosterior();



	// Require eigensystem of mode covariance
	DoubleMatrix2D R = mode.R;
	// Util.debugPrintMatrixToFile("R.txt",R);
	// form symmetric permuted covariance C = R^T.R
	DoubleMatrix2D C = R.like();
	SeqBlas.seqBlas.dgemm(true,false,1,R,R,0,C);

	// Util.debugPrintMatrixToFile("C.txt",C);

	// now decompose
	EigenvalueDecomposition eDecomp = new EigenvalueDecomposition(C);
	DoubleMatrix2D Lambda = eDecomp.getD();
	DoubleMatrix2D U = eDecomp.getV();
	if (false)
	  {
		DoubleMatrix2D Utmp = U.like().assign(U);
		DoubleMatrix1D vec;
		for(i=0;i<n;i++)
		  {
			DoubleMatrix1D v = Utmp.viewColumn(i);
			v.assign(Mult.mult(Math.sqrt(Math.max(Lambda.getQuick(i,i),0.0)))); // scale eigenvector by sqrt(eigenvalue) 
		  }
		Util.debugPrintMatrixToFile("Uscaled.txt",Utmp);
	  }

	//Util.debugPrintMatrixToFile("D.txt",eDecomp.getRealEigenvalues());
	DoubleMatrix1D eM = centre_point.like();
	DoubleMatrix1D eWork = centre_point.like();
	DoubleMatrix1D mu_post = centre_point;

	for (i=0; i<n; i++)
	  {
		// eM = perturbations in principal directions space
		eM.assign(0);
		scale = Math.sqrt(Math.max(Lambda.getQuick(i,i),0.0));
		for (j=0,ii=-nXsectionSteps; ii<=nXsectionSteps; ii++,j++)
		  {
			lambda = N_STD_DEVIATIONS*((double)ii/(double)nXsectionSteps);
			eM.setQuick(i,lambda*scale);
			U.zMult(eM,eWork,1,0,false);
			// now eWork  = U.eM

			// now form P.eWork and write into currentParams[]
			for (jj=0;jj<n;jj++) 
			  {
				val = eWork.getQuick(jj) + centre_point.getQuick(jj);
				property_sample[jj] = val;
				setValuesFromIndex(jj,val);    
			  }
			// must get times in right order
			orderTimes();
			good_state = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
			if (good_state)
			  {
				ans = debugComputeNegLogPosterior();
			  }
			else ans=0.0;// call function with perturbation
			// store result
			// fvals[j] = ans;
			// xvals[j] = lambda;
			// compute to quadratic approximation from mode centre:
			// perfect quadtratic (prior) is: ans = mode_ans + lambda*lambda/2.0
			if (good_state && (ans < mode_ans + THRESHOLD * lambda*lambda/2.0))
			  {
				_SI.EH.println("BIG STEP MODE "+i);
				DoubleMatrix1D vec = new DenseDoubleMatrix1D(n);
				vec.assign(U.viewColumn(i));
				vec.assign(Mult.mult(scale)); // scale eigenvector by sqrt(eigenvalue) 
				directions.add(vec);          // add to list
				// can break out of loop: one outlier is enough
				break;
			  }
		  }
	  }
	// assign directions to mode
	if (directions.size()>0)
	  {
		mode.bigJumpDirections = directions;
	  }
	else directions=null;
  }



  /** 
	  A method to examine the posterior surface near the ML (maximum
	  likelihood) points, by dumping files with Xsections *along principal 
	  components direction*. The surface generally looks smooth and quadratic, but
	  definite steps are evident where pinchouts come into play, plus
	  kinks at time pinchouts. Xsections
	  are writtenas as ascii files, with the supplied prefix, viz
	  prefix.property_name.layer_number.txt. Must be called with a
	  completely defined mode, where newmode.QP_optimum is expected to
	  be the ML point.
  */
  private void debugPrintPrincipalXsections(String prefix, LocalMode newmode)
  {
	int i,ii,j,jj,k,N;
	int N_XSECTION_STEPS = 20;
	double N_STD_DEVIATIONS = 4;
	int n,iVal;
	boolean good_state=true;
	double val, newval, ans, scale,lambda;
	double[] fvals;
	double[] xvals;
	DoubleMatrix1D centre_point=null;
	String name=null;

	N=2*N_XSECTION_STEPS;
	fvals=new double[N+1];
	xvals=new double[N+1];

	// point everything at mode we set.
	setMode(newmode);
	// central point through which to draw Xsections
	centre_point = mode.QP_optimum;

	n = centre_point.size();

	// Require eigensystem of mode covariance
	DoubleMatrix2D R = mode.R;
	// Util.debugPrintMatrixToFile("R.txt",R);
	// form symmetric permuted covariance C = R^T.R
	DoubleMatrix2D C = R.like();
	SeqBlas.seqBlas.dgemm(true,false,1,R,R,0,C);

	// Util.debugPrintMatrixToFile("C.txt",C);

	// now decompose
	EigenvalueDecomposition eDecomp = new EigenvalueDecomposition(C);
	DoubleMatrix2D Lambda = eDecomp.getD();
	DoubleMatrix2D U = eDecomp.getV();
	//Util.debugPrintMatrixToFile("U.txt",U);
	//Util.debugPrintMatrixToFile("D.txt",eDecomp.getRealEigenvalues());
	DoubleMatrix1D eM = centre_point.like();
	DoubleMatrix1D eWork = centre_point.like();
	DoubleMatrix1D mu_post = centre_point;

	// X section traces
	DoubleMatrix2D traces =new DenseDoubleMatrix2D(N+1,n);

	for (i=0; i<n; i++)
	  {
		// eM = perturbations in principal directions space
		eM.assign(0);
		scale = Math.sqrt(Math.max(Lambda.getQuick(i,i),0.0));
		for (j=0,ii=-N_XSECTION_STEPS; ii<=N_XSECTION_STEPS; ii++,j++)
		  {
			lambda = N_STD_DEVIATIONS*((double)ii/(double)N_XSECTION_STEPS);
			eM.setQuick(i,lambda*scale);
			U.zMult(eM,eWork,1,0,false);
			// now eWork  = U.eM

			// now form P.eWork and write into currentParams[]
			for (jj=0;jj<n;jj++) 
			  {
				val = eWork.getQuick(jj) + centre_point.getQuick(jj);
				property_sample[jj] = val;
				setValuesFromIndex(jj,val);    
				traces.setQuick(j,jj,val);
			  }
			// must get times in right order
			orderTimes();
			good_state = _m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK	
			if (good_state)
			  {
				ans = debugComputeNegLogPosterior();
			  }
			else ans=0.0;// call function with perturbation
			// store result
			fvals[j] = ans;
			xvals[j] = lambda;
		  }
		String key = new String(prefix+"."+i+".txt");
		Util.printfileZO(key,xvals,fvals,fvals.length);
		Util.debugPrintMatrixToFile("trace."+i+".txt",traces);
	  }
  }


  /** 
	  A method to examine the posterior surface near the MAP (maximum
	  aposteriori likelihood) points. The Newton solver is sometimes
	  fooled by wierd pinch out points, and gets stuck at little
	  "ledges", close to a steep downhill slope.  To detect this we
	  will check the posterior surface at points where the time values
	  are sampled at -1,-0.5,+0.5,+1 std deviations away from the MAP
	  point, taken one at a time (like a very simple Box Benkhen
	  scheme). Points of posterior probability lower than some Newton
	  threshold less than the MAP point are taken as evidence of the
	  above pathology, and flagged.  If so, we return true, to
	  indicate wierdness of the surface (probable false convergence)
	  Must be called with a completely defined mode, where
	  newmode.QP_optimum is expected to be the MAP point.
  */
  private boolean modeHasFaultyNewtonConvergence(LocalMode newmode)
  {
	int i=0,ii,i_m,j,jj,j_m,k,kk,n;
	int NSTEPS=2;  // left and right of MAP point
	int N_MAP_PROFILE_STEPS=20;
	double FRACTION_OF_STDDEV_STEP=0.5;

	double val=0, newval, ans, scale=1;
	DoubleMatrix1D centre_point=null;
	double MAPnegLogPost;
	boolean looksFaulty = false;

	// point everything at mode we set.
	setMode(newmode);
	// central point through which to draw Xsections
	centre_point = mode.QP_optimum;

	n = centre_point.size();
	// set MAP point
	for (k=0;k<n;k++)
	  {
		val = centre_point.getQuick(k);
		setValuesFromIndex(k,val);
		property_sample[k]=val; // needed for prior
	  }
	orderTimes();
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	

	MAPnegLogPost = debugComputeNegLogPosterior();
	if (_SI.verbosity>4) _SI.EH.println("MAP centre-point  %10.4g ",MAPnegLogPost);

	// loop over times, including base time
	for (i_m=0; i_m<=_m.n_active_layers; i_m++)
	  {
		// i = master_layer_number_sequence[i_m];
		// get the base value of the property
		val = centre_point.getQuick(i_m);
		// (a) formally correct with proper invCovariance
		scale = FRACTION_OF_STDDEV_STEP/Math.sqrt(mode.CtildeInvDiag[i_m]);  // try +/- 0.5, +/-1 std dev
		boolean zero_mean = (val==0.0);
		// run over -2,-1,0,1,2 sample points. NB Step (ii == NSTEP+1) is to restore value
		for (j=0,ii=-NSTEPS;ii<=NSTEPS+1;ii++,j++)
		  {
			if (ii==0) continue; // skip MAP point
			if (ii==NSTEPS+1) 
			  {
				newval = val; // restore on last loop pass
			  }
			else
			  {
				// times
				newval = val+scale*(double)ii;
			  }

			// Always reset all times, because of polluting effect of orderTimes()
			for (k=0;k<=_m.n_active_layers;k++)
			  {
				if (k==i_m) continue;  // index we are testing
				setValuesFromIndex(k,centre_point.getQuick(k));
			  }
			setValuesFromIndex(i_m,newval);
			property_sample[i_m]=newval;     // must set this for prior evaluation

			// must get times in right order
			orderTimes();

			// For fluid contacts, use this...
			if (_SI.doSimpleFluidUnitsAndContacts)
			  {
				// time perturbation will affect depths and => internal layer FC's
				_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
			  }

			ans = debugComputeNegLogPosterior();
			// examine result
			if (ii != NSTEPS+1) 
			  {
				if (_SI.verbosity>3) 
				  {
					_SI.EH.println("MAP check (master)t.%d, %10.4g std-devs: %10.4g ",
									  i_m, ii*FRACTION_OF_STDDEV_STEP, ans-MAPnegLogPost);
				  }
				if (ans-MAPnegLogPost < -NEWTON_EPS)
				  {
					looksFaulty=true;
					break;
				  }
			  }
			// restore times 
			// ordering times can mess things up, so restore all of them.
			for (jj=0;jj<=_m.n_active_layers;jj++) 
			  {
				setValuesFromIndex(jj,(newval=centre_point.getQuick(jj)));
				property_sample[jj]=newval;
			  }
			orderTimes();
		  }
		if (looksFaulty) break;
	  }
	if (!looksFaulty && _SI.doSimpleFluidUnitsAndContacts)
	  {
		// clean up after last orderTimes()
		// time perturbation will affect depths and => internal layer FC's
		_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
	  }
	if (looksFaulty)
	  {
		if (_SI.useMAPcheck)
		  {
			// N.B. we have broken out on master layer index i_m
			double min=Double.MAX_VALUE, minval=0;
			if (_SI.verbosity>3) _SI.EH.print("MAP Profile: ");
			for (j=0,ii=-N_MAP_PROFILE_STEPS;ii<=N_MAP_PROFILE_STEPS;ii++,j++)
			  {
				// restore times at centre point first
				// (ordering times can mess things up)
				for (j_m=0; j_m<=_m.n_active_layers; j_m++) 
				  {
					// jj = master_layer_number_sequence[j_m];
					if (j_m == i_m) continue;  // can skip faulty layer
					setValuesFromIndex(j_m,(newval=centre_point.getQuick(j_m)));
					property_sample[j_m]=newval;
				  }
				// then offset faulty-layer time
				newval = val+scale*(double)ii;
				property_sample[i_m]=newval;     // must set this for prior evaluation
				setValuesFromIndex(i_m,newval);
				orderTimes();
				// For fluid contacts, use this...
				if (_SI.doSimpleFluidUnitsAndContacts)
				  {
					// time perturbation will affect depths and => internal layer FC's
					_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
				  }
				ans = debugComputeNegLogPosterior();
				if (ans < min)
				  {
					min = ans;
					minval = newval;
				  }
				if (_SI.verbosity>3) _SI.EH.printf("%10.4g ",ans-MAPnegLogPost);
			  }
			if (_SI.verbosity>3) _SI.EH.println("");
			// restore things at the culprit point
			property_sample[i_m]=minval;     // must set this for prior evaluation
			for (j_m=0; j_m<=_m.n_active_layers; j_m++) 
			  {
				// jj = master_layer_number_sequence[j_m];
				// if (jj==i) continue;  // can skip faulty layer
				if (j_m == i_m) continue;  // can skip faulty layer
				setValuesFromIndex(j_m,centre_point.getQuick(j_m));
			  }
			setValuesFromIndex(i_m,minval);
			orderTimes();
			// For fluid contacts, use this...
			if (_SI.doSimpleFluidUnitsAndContacts)
			  {
				// time perturbation will affect depths and => internal layer FC's
				_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false); // FC_OK
			  }
			ans = debugComputeNegLogPosterior();
			if (_SI.verbosity>3) _SI.EH.println("MAP Improved point found with f = %10.8g",ans);
		  }
		else
		  {
			if (_SI.verbosity>3) _SI.EH.println("MAP Check fails, but continuing on blindly (MAPcheck disabled)");
			looksFaulty = false;
		  }	
	  }
	return looksFaulty;
  }

  /** 
	  A method to examine the posterior surface along the lines
	  joining the ML (maximum likelihood) points of each mode, for
	  those subsets of the set of localModes which have the same fluid
	  state. Modes are examines pairwise within these subsets, and if
	  the posterior profile does not show a "barrier" (i.e. a "hill to
	  get over") along the line, one of the modes is very likely
	  spurious. The mode with the maximum -ve log posterior is then
	  eliminated, and the process iterated until no more eliminations
	  can be performed in all fluid states.

	  Returns true if any modes were removed and some renormalisation required
  */
  private boolean checkModeJoiningProfiles(boolean printAsciiProfiles)
  {
	int i,ii,j,jj,k,kk,k_mark,n;
	int NSTEPS = 25;                          // number of points along line joining mode centers
	double MODE_ELIMINATION_WEIGHT_EPS = 0.5;
	boolean good_state=true;
	double max,val,frac,step=1.0/(NSTEPS-1);
	double[] xvals = new double[NSTEPS];
	double[] fvals = new double[NSTEPS];
	DoubleArrayList fvalsList = new DoubleArrayList(fvals);
	DoubleMatrix1D start_point=null;
	DoubleMatrix1D end_point=null;

	boolean mustRenormalise=false;
	boolean[] markedForDeletion = null;
	LocalMode m_k,m_kk;

	// loop over modes in localModes suite to find groups of modes with same fluid state
	markedForDeletion = new boolean[localModes.size()];
	for (i=0;i<localModes.size();)
	  {
		// ok we can only draw lines between modes in the same dimension space, so
		// the lines have to be for modes with the same fluid combination.
		// Get this range
		// point everything at mode we set.
		setMode(localModes.get(i));
		for (j=i+1;j<localModes.size();j++)
		  {
			if ((localModes.get(j)).fs != fs) break;
		  }
		// same-fluid-state  modes span the range [i,j-1] now
		if (j-1>i)
		  {
			property_sample = new double[mode.QP_optimum.size()];
			// more than two modes with this fluid state. We can draw a Profile
			// loop over all possible pairs k,kk in the range [i,j-1]
			for (k=i;k<j-1;k++)
			  {
				if (markedForDeletion[k]) continue;
				m_k = localModes.get(k);
				setMode(m_k);				
				// points through which to draw Xsections
				start_point = mode.QP_optimum;
				n = start_point.size();
				for (kk=k+1; kk<=j-1; kk++)
				  {
					if (markedForDeletion[kk]) continue;
					m_kk = localModes.get(kk);
					end_point = m_kk.QP_optimum;

					// draw Xsection
					for (ii=0,frac=0.0;ii<NSTEPS;ii++,frac+=step)
					  {
						// interpolate model vector
						for (jj=0;jj<n;jj++)
						  {
							val = (1.0-frac)*start_point.getQuick(jj) + frac*end_point.getQuick(jj);
							setValuesFromIndex(jj,val);
							property_sample[jj]=val; // needed for prior
						  }
						orderTimes();
						good_state = _m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FC_OK	
						if (good_state)
						  {
							val = debugComputeNegLogPosterior();
						  }
						else val=0.0;
						// store result
						xvals[ii] = k+frac;
						fvals[ii] = val;
					  }
					// print profile to file if required
					if (printAsciiProfiles) 
					  {
						String name = new String("Profile.fs="+fs+"."+k+"-"+(kk)+".txt");
						Util.printfileZO(name,xvals,fvals,fvals.length);
					  }

					// Possibilities for failure

					// (1) We can have instances of two modes that ought to be the same, differing by only
					// floating point errors.
					max = cern.jet.stat.Descriptive.max(fvalsList);
					if ((max-fvals[0] < MODE_SEPARATION_EPSILON) && (max-fvals[NSTEPS-1] < MODE_SEPARATION_EPSILON))
						//	&& !(max==fvals[0] || max==fvals[NSTEPS-1]))
					  {
						// One of the modes looks spurious. 
						// check if mode weights differ significantly. If so, eliminate on this basis
						if (m_k.weight < MODE_ELIMINATION_WEIGHT_EPS * m_kk.weight)
						  {
							k_mark = k;  // mode number to remove
						  }
						else if (m_kk.weight < MODE_ELIMINATION_WEIGHT_EPS * m_k.weight)
						  {
							k_mark=kk;
						  }
						else
						  {
							// Eliminate the on the basis of the
							// most pinched-out one (often happens for priors with enforced pinchouts)
							if (m_k.layer_code < m_kk.layer_code)
							  { 
								k_mark = k; // mode number to remove
							  }
							else
							  {
								k_mark = kk; // mode number to remove
							  }
						  }
						// localModes.remove(k_mark);
						markedForDeletion[k_mark] = true;
						if (_SI.verbosity>3) 
						  _SI.EH.println("Mode-fence profile test: Eliminated, on weight basis, spurious (near-match) mode "+k_mark);
						if (k_mark==k) break; // can skip rest of kk loop of we lose the start index
					  }
					else
					  {
						// (2) test if either endpoint mode looks spurious
						if (max==fvals[0])
						  {
							// starting mode looks spurious. Eliminate it, and break from kk loop.
							// localModes.remove(k);
							markedForDeletion[k] = true;
							if (_SI.verbosity>3) 
							  _SI.EH.println("Mode-fence profile test: Eliminated spurious (start) mode "+k);
							break;
						  }
						else if (max==fvals[NSTEPS-1])
						  {
							// end mode looks spurious. Eliminate it
							// localModes.remove(kk);
							markedForDeletion[kk] = true;
							if (_SI.verbosity>3) 
							  _SI.EH.println("Mode-fence profile test: Eliminated spurious (end) mode "+(kk));
						  }
					  }
					
				  }
			  }
			property_sample=null; // GC
		  }
		// point i to start of new fluid state group
		i=j;
	  }

	// Now loop in reverse sequence through modes, deleting marked ones
	for (i=localModes.size()-1;i>=0;i--)
	  {
		if (markedForDeletion[i]) 
		  {
			localModes.remove(i);
			// _SI.EH.println("Mode-fence profile test: removing mode "+i);
			mustRenormalise=true;
		  }
	  }
	localModes.trimToSize();

	markedForDeletion = null; // GC
	return mustRenormalise;
  }

  double debugComputeNegLogPosterior()
  {
	int i,j;
	double dt, w_sum_sq=0.0, negLogLik, term;
	double negLogLik_iso, negLogLik_depth, negLogLik_seis,negLogLik_stack;
	double negLogLikPriorOnly, negLogLikSeisIsopachDepthVS;
	

	// (1) SEISMIC mismatch term
	negLogLik_seis = negLogLikSeisOnly = _m.SM.computeSyntheticAndSeismicError();
	negLogLikHiddenSeisOnly = _m.SM.hiddenSeismicError;
	negLogLik = negLogLik_seis;

	// (2) Add in logLikelihood contribution from ISOPACHS if they're required
	negLogLik_iso = computeIsopachNegLogLikelihood(false);
	negLogLik += negLogLik_iso;

	// (3) Add in logLikelihood contribution from DEPTH self-consistency if required
	negLogLik_depth= computeConsistentDepthNegLogLikelihood(false); 
	negLogLik += negLogLik_depth;

	// (4) Add in logLikelihood contribution from stacking self-consistency if required
	if (_SI.doStackedModel) 
	  {
		negLogLik_stack= _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
		negLogLik += negLogLik_stack;
	  }



	negLogLikSeisIsopachDepthVS = negLogLik;


	// (4) Also the contribution from the raw PRIOR
	negLogLikPriorOnly = computePriorNegLogLikelihood();
	
	negLogLik += negLogLikPriorOnly;

	return negLogLik;
  }







































  // testing routine
  private void debugShowMatrices()
  {
	int i,j,n;

	n = Cinv.rows();
	_SI.EH.print("%nMatrix:%n");
	_SI.EH.print("   ");
	// col numbers
	//for (i = 0;i<n;i++)	_SI.EH.printf("%9d",p.add(i));
	_SI.EH.print("%n");
	for (i = 0;i<n;i++)
	  {
		// row number
		for (j = 0;j<n;j++)
		  {
			_SI.EH.printf("%9.3g ",Cinv.getQuick(i,j));
		  }
		_SI.EH.printf("  %8g%n",mu.getQuick(i));
	  }
	_SI.EH.print("%n");
  }



  // testing routine
  private void debugShowMatrix(String name, DoubleMatrix2D A)
  {
	int i,j,n;

	n = A.rows();
	_SI.EH.printf("%n%s%n",name);
	for (i = 0;i<n;i++)
	  {
		for (j = 0;j < A.columns();j++)
		  {
			_SI.EH.printf("%10.4g ",A.getQuick(i,j));
		  }
		_SI.EH.print("%n");
	  }
  }

  // testing routine
  private void debugShowMatrix(String name, DoubleMatrix1D A)
  {
	int i,j,n;

	n = A.size();
	_SI.EH.printf("%n%s%n",name);
	for (i = 0;i<n;i++)
	  {
		_SI.EH.printf("%10.4g ",A.getQuick(i));
		_SI.EH.print("%n");
	  }
  }


  // testing routine
  private void debugTestMatrixSizes(int fs)
  {
	int i;

	LayerProperties.Property prop;
	ArrayList AL = null;

	_SI.EH.printf("%nFluid combo %d:%n",fs);

	for (i=0;i<map_from.size();i++)
	  {
		AL = (ArrayList) map_from.get(i);
		prop = (LayerProperties.Property)AL.get(0);
		_SI.EH.printf("%2d, list %2d, mean=%6.2f, val=%6.2f%n",
						i, AL.size(), prop.mean, prop.value);
	  }
  }

  /**
	 Special class for setting up useful info about sparse matrices
  */

  public class MatrixColumnInfo implements cern.colt.function.IntIntDoubleFunction
  {
	/**
	   colmins[j] = smallest row index that is non-zero for column j
	*/
	int[] colmins;

	/**
	   Setup the information colmins for the (usually sparse) matrix A.
	   colmins[j] = smallest row index that is non-zero for column j
	   It can then be accessed directly via this.colmins
	*/
	public MatrixColumnInfo(DoubleMatrix2D A)
	{
	  int n = A.rows();
	  colmins = new int[n];
	  Arrays.fill(colmins,n-1);
	  A.forEachNonZero(this);
	}

	/**
	   This method is called only within the instantiator
	   forEachNonZero() call on the sparse matrix.  It 
	   computes the colmins values.
	*/
	public double apply(int i, int j, double value) 
	{
	  if (i<colmins[j]) colmins[j]=i;
	  return value;
	}
  }
}
