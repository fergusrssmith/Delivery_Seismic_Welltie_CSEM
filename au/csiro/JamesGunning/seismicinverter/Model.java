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
import au.csiro.JamesGunning.lib.ErrorHandler;
import java.util.HashMap;

  /**
	 Full suite of parameters for Bayesian seismic inverse calculations,
	 grouped into convenient classes
  */

public final class Model
{
  int n_layers;
  int n_active_layers;

  ModelPriorDescriptions MPD;
  RockandFluidData RFD;
  LayerProperties[] LP;
  HashMap completePropertyNameMap=null;

  OneDSeismicConvolutionalModel SM;
  FluidStates FS;
  ModelSampler MS;  
  StackedModelInfo stackedModelInfo;
  FluidUnitsAndContacts  fluidUnitsAndContacts;
  boolean needHardRockTerms=false;
  boolean needFloatingGrainTerms=false;
  int master_depth_layer_number=1;

  /**
	 Each instance of Model has a parent SeismicInverter, from which
	 we control all the options, error-handling etc
  */
  SeismicInverter _SI;

  static int debugCounter=0; // FIXME remove!

  /**
	 create instance of a Model, details to be filled in with setUpFully() when known
  */
  protected Model(SeismicInverter SI)
  {
	this._SI = SI;
  }

  /**
	 set up model details fully once number of layers are known
  */
  protected void setUpFully(int n_layers,ErrorHandler EH,int verbosity)
  {
	// Set model dimension here
	// suitable defaults
	this.n_layers = this.n_active_layers = n_layers;	
	MPD = new ModelPriorDescriptions(this);	
	RFD = new RockandFluidData();

	FS = new FluidStates(this);
	MS = new ModelSampler(this);

	// we need a set of LayerProperties for each layer
	// set up properties values for each layer, plus a map of property names
	completePropertyNameMap = new HashMap();
	LP = new LayerProperties[n_layers];
	for (int i=0;i<n_layers;i++)
	  {
		boolean needsTBase = (i==n_layers-1);
		this.LP[i] = new LayerProperties(this,i,needsTBase);
	  }
  }

  protected void setUpFluidContactInfoInLayers()
  {
	int i;

	if (_SI.doSimpleFluidUnitsAndContacts)
	  {
		for (i=0;i<n_layers;i++)
		  {
			LP[i].setUpFluidContactInfo(MPD.MetaInfo[i].is_permeable);			 
		  }
	  }
  }

 /**
	 Given a model trace tr,
	 generate stochastic forward models from it
  */
  protected void resetPriorTimesFromDepthsIfNeeded()
  {	
	if (_SI.overwritePriorTimesFromDepths)
	  {
		MS.drawMostLikelyPrior(true);
		MS.resetModelPriorTimesFromInSituState();
	  }
  }


  /**
	 Given a model trace tr,
	 generate stochastic forward models from it
  */
  protected void generateForwardModelsFromPrior()
  {	
	if (_SI.makeSyntheticAndWavelet || _SI.drawPriorPlusIsopachModel)
	  {
		// makes single realisation of most likely model
		MS.drawMostLikelyPrior(true);
		SM.generateSyntheticTrace();
		SM.dumpSynthetics(true,1,true);
		// dump model
		if (_SI.dumpRealisations) 
		  {
			_SI.dumpRealisation(true,1,true,false);
		  }
		if (_SI.verbosity>2) 
		  {
			if (_SI.ignoreIsopachs)
			  {
				_SI.EH.println("Wrote MAP model from prior only");
			  }
			else
			  {
				_SI.EH.println("Wrote MAP model from prior + isopachs only");
			  }
		  }
	  }
	else
	  {
		MS.drawSetOfPriorSamples(_SI.n_samples_per_trace_undecimated,true);
	  }
  }
  
  /**
	 Sets up the prior model. Enumerates all local modes using various heuristics.
	 Samples from posterior using information from these local modes.
  */
  protected void drawSamples()
  {
	MS.enumerateAllLocalModes();
	if (_SI.drawMAPmodels || _SI.drawMinChisqModels)
	  {
		MS.dumpMAPmodels();
	  }
	// Sample only if other outputs not requested
	// 	if (!_SI.dumpOptimalSeismic && !_SI.drawMAPmodels)
	if (!_SI.dumpOptimalSeismic)
	  {
		MS.sampleFromLocalModes();
	  }

  }


  void determineIfHardRockModelNeeded()
  {
	// hard rocks
	needHardRockTerms = false;

	if (_SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		ModelPriorDescriptions.MetaInfo MI;
		for (int i=0;i<n_layers;i++)
		  {
			MI = MPD.MetaInfo[i];
			if ((MI.hard_shale_rock_member != null) || (MI.fluid_dep_hard_shale_rock_members !=null)) 
			  {
				if (MI.has_hard_rock_mixture)
				  {
					// it's a non-trivial mixture that may have hard-rock terms
					needHardRockTerms = true; 
					break;
				  }
			  }
		  }	
	  }
  }


  void determineIfFloatingGrainModelNeeded()
  {
	// hard rocks
	needFloatingGrainTerms  =false;

	if (_SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		ModelPriorDescriptions.MetaInfo MI;
		for (int i=0;i<n_layers;i++)
		  {
			if (MPD.MetaInfo[i].needs_float_capability)
			  {
				needFloatingGrainTerms = true; 
				break;
			  }
		  }	
	  }
  }



  void setUpReducedModel()
  {
	performBasicSanityCheck();
	// Possible fluid states might vary if NG varies etc. Update isPermeable flags
	MS.updateIsPermeableFlagAndNIsopachs();	
	FS.setUpFluidStates();

	// if trends vary, we need this. Do not force remapping, since we have already done this in parsing,
	// and trends will remap internally if they need to (areal variation)
	RFD.resetTrendData(this,false);
	RFD.checkRockSpeedRegressions(this,false,_SI.rockandfluids_xml_filename);

	// Put together indexing for linear algebra
	MS.setUpMappings();
	MS.setUpPriorMatrices();
	MS.testGassmanSanity();
	// Set up Stacked Velocity model detail if used
	if (_SI.doStackedModel)
	  {
		// Set up exclusion zone if required		
		SM.computeExclusionZones();
	  }
	// get design Matrix X and d for isopachs, depth self-consistency
	MS.incorporateLinearisedIsopachAndOrDepthConstraints();
  }

  /**
	 Suite of basic sanity checks.
  */
  private void performBasicSanityCheck()
  {
	double t_start = MPD.MetaInfo[0].t.value;
	double t_end = MPD.t_base.value;

	// (1) Span of model is positive
	if (t_end-t_start < 0) 
	  {
		throw new RuntimeException("Trace " + _SI.n_traces + 
								   ": mean time-span of modelling region ("+
								   t_start+","+t_end+") is negative!!!");
	  }
	// (2) add whatever you like here......

  }
  

  /**
	 Computes full suite of effective layer properties: to
	 be used after making a change in the linear algebra
	 representation.
	 If must_throw_exception==true, will throw an exception if something goes
	 wrong, otherwise will return false.
	 Usual call has all arguments false: these are set only for special purpose feedback.
  */
  boolean computeEffectivePropertiesForAllLayers(boolean must_throw_exception, boolean testGassmanSanity,
												 boolean showIfGassmanSubstitutionPossible, boolean showDiagnostics,
												 boolean showGassmanDiagnostics)
  {
	int i,i_inc;
	boolean success,descending;
	double t_layer_base,suspended_depth_hook;

	// all depths hung or supported from master_depth_layer_number (unit offset)
	// I.e. if master_depth_layer_number==1, this is the top layer.
	LP[master_depth_layer_number-1].suspended_depth_top  = LP[master_depth_layer_number-1].d.value;
	suspended_depth_hook = LP[master_depth_layer_number-1].suspended_depth_top;

	// The following wraps up these two loops: 
	// From master layer downwards, then master layer upwards...

	// for (i = master_depth_layer_number-1; i<(n_layers-1); i++) ;
	// for (i = master_depth_layer_number-2; i>=0 ; i--);
	for (i_inc=1; i_inc>=-1; i_inc -=2)
	  {
		for (i = ((i_inc==1) ? master_depth_layer_number-1 : master_depth_layer_number-2); 
			 (i_inc==1)? (i<n_layers) : (i>=0) ; 
			 i += i_inc)			 
		  { 
			// Main call: check if combinations are unphysical
			// Update suspended depths
			t_layer_base = (i < (n_layers-1)) ? LP[i+1].t.value : LP[n_layers-1].t_base.value;
			descending = (i_inc == 1);
			
			if (descending)
			  {			
				// Going down: 
				// 0 <= master_depth_layer_number-1 <= i <= n_layers-1
				if (i > master_depth_layer_number-1) suspended_depth_hook = LP[i-1].suspended_depth_bottom; // bottom of layer above
			  }
			else
			  {
				// Going up (master_depth_layer_number-2 >= i >= 0)
				suspended_depth_hook = LP[i+1].suspended_depth_top; // top of layer below
			  }
			success = LP[i].computeEffectiveProperties(MS.fluid_state.types[i],false,t_layer_base,suspended_depth_hook,descending);
			
			if (!showGassmanDiagnostics)
			  {
				if (!success)
				  {
					if (showIfGassmanSubstitutionPossible) return false;  // early bailout if a failure
					if (!must_throw_exception) 
					  {
						return false;
					  }
					else
					  {
						// Appropriate messages for diverse contexts which call this function
						if (testGassmanSanity && _SI.checkPriorSanity) // only occurs when called from testGassmanSanity()
						  {
							// a poor prior. Be noisy about it
							_SI.EH.println("Fluid substitution error, layer "+(i+1)+": Requested layer state is unphysical. "+
										   "Your prior model is shonky for layer " + i +", fluid state " + MS.fs+". Some diagnostics...");
						  }
						else
						  {
							//  probably just unlucky. but something might be spastic. Show some diagnostics.
							_SI.EH.println("Fluid substitution error: layer "+(i+1)+": some diagnostics...");
						  }
						LP[i].computeEffectiveProperties(MS.fluid_state.types[i],true,t_layer_base,suspended_depth_hook,descending);  // call for gory detail

						if (!(showDiagnostics || (testGassmanSanity && _SI.checkPriorSanity)))
						  {
							// Throw an error
							throw new RuntimeException("Requested layer state is unphysical. Your "+
													   (testGassmanSanity ? "prior ":"")+
													   "model is shonky for layer " + (i+1) +
													   ", fluid state " + MS.fs);	
						  }
					  }
				  }
			  }
			else
			  {
				_SI.EH.println("Gassman_test: %d %10.4f %10.4f %10.4f %s",i, LP[i].phi_s.value, 
							   LP[i].vp_s.value,LP[i].vs_s.value,(success?"ok":"bad"));
			  }
		  }
	  }

	return true;
  }



  void initialiseLayerSynching()
  {
	int i;
	ModelPriorDescriptions.MetaInfo MI;

	// time synching
	
	// top layer must always be active
	MPD.MetaInfo[0].isActive = true;
	n_active_layers = 1;

	if (false)
	  {
		// Emergency code to clean up models with pinchouts that should have been forced etc
		_SI.EH.println("FIXME EMERGENCY CODE!!");
		boolean valid=false;
		for (i=1;i<n_layers;i++)
		  {
			if (MPD.MetaInfo[i].t.value > MPD.MetaInfo[0].t.value) valid=true;
		  }
		if (!valid)
		  {
			for (i=1;i<n_layers;i++)
			  {
				MI = MPD.MetaInfo[i];
				if (MI.sigma_dz!= null) MI.sigma_dz.value=0;
				if (MI.dz!=null) MI.dz.value=0;
			  }
		  }
	  }

	for (i=1;i<n_layers;i++)
	  {
		// Time sunching criterion
		MI = MPD.MetaInfo[i];

		if (_SI.DETECT_LAYER_SYNCHING && (MI.dz!=null) && (MI.dz.value==0) 
			&& (MI.sigma_dz!= null) && (MI.sigma_dz.value == 0))
		  {
			MI.isActive = false;
		  }
		else
		  {
			MI.isActive = true;
			n_active_layers++;
		  }
	  }

	if (n_active_layers<1) throw new RuntimeException("Must have at least one active layer in model (not all sigma_dz=0)");
  }

  /**
	 Simple precaluations to save log() calls, set up safe synthetic times etc
  */
  
  void computeEfficiencyInitialisations()
  {
	int i;
	for (i=0;i<n_layers;i++)
	  {
		ModelPriorDescriptions.MetaInfo MI = MPD.MetaInfo[i];
		if ((MI.sigma_dz != null) && (MI.sigma_dz.value > 0)) 
		  {
			MI.log_sigma_dz = Math.log(MI.sigma_dz.value);
		  }
	  }
	SM.precomputeLogNoise();	
  }
  
  void computeSyntheticWindowInitialisations()
  {
	int i;
	double t;

	// compute safe ranges for synthetics from means
	t = MPD.MetaInfo[1].t.value;
	for (i=2;i<n_layers;i++) t = Math.min(t,MPD.MetaInfo[i].t.value);
	MPD.min_interesting_reflector_time_mean = t; 

	t = MPD.MetaInfo[n_layers-1].t.value;
	for (i=n_layers-2; i>=0 ; i--) t = Math.max(t,MPD.MetaInfo[i].t.value);
	MPD.max_interesting_reflector_time_mean = t;
  }
}
