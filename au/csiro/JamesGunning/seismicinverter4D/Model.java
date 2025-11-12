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
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;

import java.util.HashMap;
import java.util.List;
import java.util.Arrays;
import static java.lang.Math.*;

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
  LayerProperties[][] L4P; // 4D replacement of above

  // 4D stuff
  TimeLapseInfo timeLapseInfo = null;
  HashMap<LayerProperties.Property,String>completePropertyNameMap = null;
  OneDSeismicConvolutionalModel SM;
  FluidStates FS;
  ModelSampler MS;  
  StackedModelInfo stackedModelInfo;
  FluidUnitsAndContacts  fluidUnitsAndContacts;
  boolean needHardRockTerms=false;
  boolean needFloatingGrainTerms=false;
  boolean is4D = false;
  boolean need4DTerms = false;
  int Nvintages=1;
  int MaxOffsetsPerVintage=1;
  int[][] stackVintageLookup=null;
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
	completePropertyNameMap = new HashMap<LayerProperties.Property,String>();
	LP = new LayerProperties[n_layers];
	for (int i=0;i<n_layers;i++)
	  {
		boolean needsTBase = (i==n_layers-1);
		this.LP[i] = new LayerProperties(this,i,needsTBase,(byte)0); // presume starting vintage
	  } 
  }

  // call this once Nvintages is known
  protected void complete4DsetUp()
  {
	MS.complete4DMSSetUp();
  }

  void setUpLP4DInfo()
  {	
	int i;
	for (i=0;i<n_layers;i++) LP[i].init4DInfo();	
  }

  void setUpStackVintageXref()
  {
	int i,j,k;

	// stack <-> vintage cross indexing
	stackVintageLookup= new int[Nvintages][MaxOffsetsPerVintage];
	for(i=0;i<Nvintages;i++) Arrays.fill(stackVintageLookup[i],-1);  // use -1 as flag for not-present
	for(k=0;k<SM.Nstacks;k++)
	  {
		StackData stack = SM.stacks[k];
		i = stack.vintage - 1;    // vintages are unit-offset
		j = stack.offset_index;
		if(stack.DSI==null) stackVintageLookup[i][j] = k;
	  }
  }

  void complete4DLPinfo()
  {
	int i,k;
	for (k=1;k<Nvintages;k++)
	  {
		for (i=0;i<n_layers;i++)
		  {
			L4P[k][i].RFD_PermeableRegSet = LP[i].RFD_PermeableRegSet;
		  }
	  }
  }

  protected boolean reallocateModel()
  {
	int i,k;

	// snapshot of existing properties
	LayerProperties[] _LP = LP;
	L4P = new LayerProperties[Nvintages][n_layers];
	// copy existing base model
	for (i=0;i<n_layers;i++) L4P[0][i] = _LP[i];  
	// later vintages
	for(k=1;k<Nvintages;k++)
	  {
		for (i=0;i<n_layers;i++)
		  {
			boolean needsTBase = (i==n_layers-1);
			L4P[k][i] = new LayerProperties(this,i,needsTBase,(byte)k);
		  }
	  }
	
	// Also complete Map for 4D stochastic properties
	if (is4D)
	  {
		if (Nvintages > 1)
		  {
			if (timeLapseInfo.DeltaP4D!=null && (timeLapseInfo.DeltaP4D.length != (Nvintages-1)))
			  {
				return false;
			  }
		  
			// already alloced & set up at parsing ...
		
			// and reverse map for other uses e.g. model-vector feedback
			for(k=1;k<Nvintages;k++)
			  {
				completePropertyNameMap.put(timeLapseInfo.DeltaP4D[k-1],"delta_p,"+(k-1));
			  }
			for(k=1;k<Nvintages;k++)
			  {
				completePropertyNameMap.put(timeLapseInfo.DeltaT[k-1],"time_shift,"+(k-1));
			  }
		  }
	  }
	// set up stack-vintage cross-refs
	setUpStackVintageXref();

	// complete model allocations for the sampler
	MS.complete4DMSSetUp();

	// fluid property allocs
	// for case no contacts..
	MS.completeSetUpSansFluidContactInfo();
	// with contacts...
	MS.completeSetUpUsingFluidContactInfo();
	
	return true;
  }


  protected void setUpFluidContactInfoInLayers()
  {
	if (_SI.doSimpleFluidUnitsAndContacts)
	  {
		int i,i_v;

		for (i_v=0; i_v < Nvintages; i_v++)
		  {
			LayerProperties[] LP  = L4P[i_v];
			for (i=0;i<n_layers;i++)
			  {
				LP[i].setUpFluidContactInfo(MPD.MetaInfo[i].is_permeable);			 
			  }
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
		MS.drawSetOfPriorSamplesPerFluidState(_SI.n_samples_per_trace_undecimated,true);
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

  void determineIf4DModelTermsNeeded()
  {
	need4DTerms =  false;
	if (is4D)
	  {
		need4DTerms = true; 
	  }
  }


  void setUpReducedModel()
  {
	performBasicSanityCheck();
	// Possible fluid states might vary if NG varies etc. Update isPermeable flags
	// also N isopachs
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
  


  /*
	copy,set later-vintage properties from vintage=0 (base vintage)
  */
  boolean setLaterVintagePropertiesFromBaseVintage()
  {
	int i,i_v;
	double delta_p,zref;
	
	LayerProperties[] LPB = L4P[0];
	LayerProperties[] LP = null;
	zref = LPB[master_depth_layer_number-1].d.value;
	for (i_v=1; i_v < Nvintages; i_v++)
	  {
		LP = L4P[i_v];
		delta_p = timeLapseInfo.DeltaP4D[i_v-1].value;
		for (i=0; i<n_layers;i++)
		  {
			LP[i].setPropertiesFromPreviousVintage(i_v,LPB[i],delta_p);
		  }
		// propagate ref depth
		LP[master_depth_layer_number-1].d.value = zref;		
	  }
	// order contact depths if required
	orderFluidContactsAcrossVintages();

	return true;
  }



  /*
	Order fluid contacts to sensible order, e.g. depth of
	GO then OB must be increasing, for each fluid unit. With priors that are
	sensible, we might still get cross-overs in the tail area, e.g.
	z_GO < z_OB, so we have to truncate somehow.
  */
  final void orderFluidContacts()
  {
	if (!_SI.doSimpleFluidUnitsAndContacts) return;
	else if (!fluidUnitsAndContacts.contactsVariable) return; // no variable contacts
	else
	  {
		int i,j,k,Nc;
		double z,zlast=0;
		for (i=0;i<	fluidUnitsAndContacts.fluidUnits.length; i++)
		  {
			FluidUnitsAndContacts.FluidUnit fu = fluidUnitsAndContacts.fluidUnits[i];
			for(j=0;j<Nvintages;j++)
			  {
				FluidUnitsAndContacts.ContactInfo[] CIJ = fu.contactInfo[j];
				Nc = CIJ.length;
				if (Nc<2) continue;
				for (k=0; k< Nc; k++)
				  {
					FluidUnitsAndContacts.ContactInfo CI = 	CIJ[k];
					z = CI.depth.value; 
					//_SI.EH.printf("Curr contacts: %d %d %7g%n",j,k,z);
					if ((k > 0) && (z < zlast))
					  {
						CI.depth.setPosValue(zlast);
						if (_SI.verbosity>5) _SI.EH.printf("Snapped contacts: %d %d %7g %7g%n",j,k,z,zlast);
					  }
					zlast = z;
				  }			  
			  }
		  }
	  }
  }
  
  /*
	Order fluid contacts to desired ordering (Tobit-wise) in
	vintage-sequence.  Eg model may have specified that a GB interface
	rises with depth as vintages are acquired...
  */
  final void orderFluidContactsAcrossVintages()
  {
	if (!_SI.doSimpleFluidUnitsAndContacts) return;
	else if (!fluidUnitsAndContacts.contactsVariable) return; // no variable contacts
	else
	  {
		int i,j,k;
		for (i=0;i<	fluidUnitsAndContacts.fluidUnits.length; i++)
		  {
			FluidUnitsAndContacts.FluidUnit fu = fluidUnitsAndContacts.fluidUnits[i];
			FluidUnitsAndContacts.ContactInfo[] CIJm1 = fu.contactInfo[0];
			for(j=1;j<Nvintages;j++)
			  {
				FluidUnitsAndContacts.ContactInfo[] CIJ = fu.contactInfo[j];
				for (k=0; k< CIJ.length; k++)
				  {
					FluidUnitsAndContacts.ContactInfo CIJk = 	CIJ[k]; // fu.contactInfo[j][k] 
					byte code = CIJk.depthPD.previousVintageOrderingCode;
					if (code != 0)
					  {
						LayerProperties.Property P = CIJk.depth;
						LayerProperties.Property Pprev = CIJm1[k].depth;
						switch (code)
						  {
						  case 1:							
							if (P.value > Pprev.value) P.setPosValue(Pprev.value);
							break;
						  case 2:
							if (P.value != Pprev.value) P.setPosValue(Pprev.value);
							break;
						  case 3:
							if (P.value < Pprev.value) P.setPosValue(Pprev.value);
							break;
						  }
					  }
				  }
				CIJm1 = CIJ;			  
			  }
		  }
	  }
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
	int i,i_inc,i_v,iMD;
	boolean success,descending;
	double t_layer_base,suspended_depth_hook,dz=0,t,dt;
	LayerProperties[] LP,LP0;
	LayerProperties.Property P;

	if (is4D)
	  {
		setLaterVintagePropertiesFromBaseVintage();
	  }
 
	// Order within-fluid_unit fluid contacts as well!
	orderFluidContacts();
	LP0 = L4P[0];      // base vintage
	iMD = master_depth_layer_number-1;
	for (i_v=0; i_v < Nvintages; i_v++)
	  {
		LP = L4P[i_v];

		// 4D time offsets
		if (i_v > 0)
		  {		
			// reference layer time shifted by model time shift
			dt = timeLapseInfo.DeltaT[i_v-1].value;  // the shift
			t = LP0[iMD].t.value + dt;                 // shifted time
			LP[iMD].t.setPosValue(t);                  // set reference layer time, later vintage
		  }

		// all depths hung or supported from master_depth_layer_number (unit offset)
		// I.e. if master_depth_layer_number==1, this is the top layer.
		LP[iMD].suspended_depth_top  = LP[iMD].d.value;
		suspended_depth_hook = LP[iMD].suspended_depth_top;

		// The following wraps up these two loops: 
		// From master layer downwards, then master layer upwards...
		for (i_inc=1; i_inc>=-1; i_inc -=2)
		  {
			descending = (i_inc == 1);
			for (i = (descending ? iMD : master_depth_layer_number-2); 
				 descending? (i<n_layers) : (i>=0) ; 
				 i += i_inc)			 
			  { 
				// Main call: check if combinations are unphysical
				// Update suspended depths
				t_layer_base = (i < (n_layers-1)) ? LP[i+1].t.value : LP[n_layers-1].t_base.value;
				if (i_v>0) dz = LP0[i].suspended_depth_bottom - LP0[i].suspended_depth_top; 
			
				if (descending)
				  {			
					// Going down: 
					// 0 <= master_depth_layer_number-1 <= i <= n_layers-1
					if (i > iMD) suspended_depth_hook = LP[i-1].suspended_depth_bottom; // bottom of layer above
				  }
				else
				  {
					// Going up (master_depth_layer_number-2 >= i >= 0)
					suspended_depth_hook = LP[i+1].suspended_depth_top; // top of layer below
				  }
				
				// Use base time, or layer-thickness (i_v>0) to computeEffectiveProperties()
				success = LP[i].computeEffectiveProperties(MS.fluid_state.types[i],false,
														   t_layer_base , dz,  
														   suspended_depth_hook,descending);


				// set layer top times LP[i].t.value etc immediately for use in contacts model etc
				if (i_v>0)
				  {
					dt = 2.0*(LP[i].suspended_depth_bottom -LP[i].suspended_depth_top)/LP[i].vp_eff;
					if (descending)
					  {
						t = LP[i].t.value + dt;
						P = (i < (n_layers-1)) ? LP[i+1].t : LP[n_layers-1].t_base;
						P.setPosValue(t);
					  }
					else
					  {
						t = t_layer_base - dt;
						LP[i].t.setPosValue(t);
					  }
				  }
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
							LP[i].computeEffectiveProperties(MS.fluid_state.types[i],true,t_layer_base,dz,
															 suspended_depth_hook,descending);  // call for gory detail

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

	for (i=1;i<n_layers;i++)
	  {
		// Time synching criterion
		MI = MPD.MetaInfo[i];

		
		/*
		forced pinchout test hack
		  if(true)
		  {
			if((i<n_layers-1) && (MI.dz!=null) && (MI.sigma_dz!= null))
			  {
			  if(MPD.MetaInfo[i+1].t.value <= MI.t.value)
				{
				MI.dz.value=0;
				MI.sigma_dz.value =0;
				}
              }
		  }
        */
		
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

  /*
	 Simple precalculations to save log() calls, set up safe synthetic times etc
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
