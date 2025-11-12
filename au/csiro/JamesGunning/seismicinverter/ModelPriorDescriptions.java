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

import cern.colt.map.*;

import au.csiro.JamesGunning.lib.ErrorHandler;

/**
   Characteristics of layers to be determined from a model-description xml file
   
   @author james.gunning@csiro.au
   @version 1.0, 9/11/2001 
*/

public final class ModelPriorDescriptions
{
  /**
	 Ref to calling model object
  */
  Model _m;

  /**
   Error Handler
  */
  ErrorHandler EH;

  /**
	 Descriptions of properties in each layer. See MetaInfo class
  */
  MetaInfo MetaInfo[];

  /**
	 Descriptions of prior near-stack noise level for seismic
  */
  public PriorDetails near_stack_noise_rms;
  /**
	 Descriptions of prior far-stack noise level for seismic
  */
  public PriorDetails far_stack_noise_rms;

  /**
	 Descriptions of prior properties for AVO model coefficient A
  */
  public PriorDetails A_avo;
  /**
	 Descriptions of prior properties for AVO model coefficient uncertainty sigma_A
  */
  public PriorDetails sigma_A_avo;

  /**
	 Descriptions of prior properties for AVO model coefficient B
  */
  public PriorDetails B_avo;
  /**
	 Descriptions of prior properties for AVO model coefficient uncertainty sigma_B
  */
  public PriorDetails sigma_B_avo;



  /**
  flags to detect if input traces have times & sigma_times
  */
  boolean prior_times_vary = false;
  boolean prior_sigma_times_vary = false;

  /* synthetic safe bounds */
  double min_interesting_reflector_time_mean;
  double max_interesting_reflector_time_mean;

  /**
	 This will hold the number of isopachs
  */
  int n_isopachs=0;
  /* 
	 stack of layers has a base time and uncertainty in last layer
  */
  PriorDetails t_base;
  PriorDetails sigma_t_base;

  /**
	 HashMap of layer numbers sychronised with an earlier layer
  */
  OpenIntIntHashMap synchMap;

  /**
	 synch_copies[i] contains the number of copies (via synching) of layer i
  */
  int synch_copies[]=null;

  /**
	 List of prior properties, index, layer, and descriptions, which vary with each trace
  */
  protected ArrayList PriorIndexInfo;

  /**
	 A null value threshold to indicate data coming in from prior traces is dummy
  */
  public static final double PRIOR_TRACES_NULL_THRESHOLD = -900;

 
  /**
	 Allocate room for n_layers MetaInfo objects to 
	 summarise information about layer i in MetaInfo[i-1].
  */

  public ModelPriorDescriptions(Model m)
  {
	this._m = m;
	this.EH = m._SI.EH;	
	this.MetaInfo = new ModelPriorDescriptions.MetaInfo[m.n_layers];
  }
  
  /**
	 Class describing generic information about a particular layer.
	 Layers are a mixture of 2 rock types, blended with a Net to gross (fraction of
	 sand_rock_member). Fluid substitution only occurs in the sand_rock_member
  */

  public final static class MetaInfo 
  {
	/**
	   Name of layer (optional)
	*/
	String layer_name = null;               // eg "top seal"
	/**
	   Name of sand_rock_member (optional)
	*/
	String sand_rock_member = null;          // eg "sand A"
	/**
	   Name of mud_rock_member (mandatory)
	*/
	String mud_rock_member = null;         // eg "shale"

	/**
	   Name of hard_shale_rock_member (optional)
	*/
	String hard_shale_rock_member = null;         // eg "shale2"
	/**
	   Fluid case-dependent hard-rock members
	*/
	String[] fluid_dep_hard_shale_rock_members = null;

	/**
	   Whether the rock properties in this layer should 
	   always be the same as those in another
	   layer with the same end-member rock names 
	   (in adjacent layers, reveals fluid effects...)
	*/
	boolean synchronised = false;
	/**
	   flag to indicate whether layer is forcibly pinched out
	*/
	boolean isActive = true;
	/**
	   flag to indicate whether layer has depth-property dependents
	*/
	boolean hasDepthDependents = false;
	/**
	   flag to indicate whether layer has LFIV-property dependents
	*/
	boolean hasLFIVDependents = false;
	/**
	   flag to indicate whether layer has float_fraction-property dependents
	*/
	boolean hasFloat_fractionDependents = false;

	/**
	   Name of in-situ brine to be replaced by hydrocarbon (optional)
	*/
	String brine_name = null;             // eg "brine A"
	/**
	   Name of oil hydrocarbon (optional)
	*/
	String oil_name = null;             // eg "40APIoil"
	/**
	   Name of high-sat gas (optional)
	*/
	String gas_name = null;             // eg "methane"
	/**
	   Name of low-sat-gas (optional)
	*/
	String lsg_name = null;             // eg "low-sat-gas"

	
	/**
	   flag to indicate whether layer is a mixture or pure shale/mud
	*/
	boolean is_permeable=false;

	/**
	   flag to indicate whether layer has possibly a hard-rock mixture 
	*/
	boolean has_hard_rock_mixture = false;
	/**
	   flag to indicate whether layer needs floating grain apparatus 
	*/
	boolean needs_float_capability = false;

	/**
	   flag to indicate whether oil substituion will really be needed
	*/
	boolean oil_substitution_possible=false;
	/**
	   flag to indicate whether gas substituion will really be needed
	*/
	boolean gas_substitution_possible=false;
	/**
	   flag to indicate whether lsg substituion will really be needed
	*/
	boolean lsg_substitution_possible=false;

	/**
	  time of layer top (must be public for field access via Class.getField)
	*/
	public PriorDetails t;
	/**
	  std-dev of time of layer top (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_t;

	/**
	  thickness of layer (isopach) (must be public for field access via Class.getField)
	*/
	public PriorDetails dz;
	/**
	  std-dev of thickness (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_dz;
	/**
	   Logarithm of sigma_dz
	*/
	double log_sigma_dz=0.0;

	/**
	   Default Net-to-gross info {value,varies,index} (must be public for field access via Class.getField)
	*/
	public PriorDetails NG;

	/**
	   Default Net-to-gross sd dev info {value,varies,index} (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_NG;

	/**
	   Default hard-Net-to-gross info {value,varies,index} (must be public for field access via Class.getField)
	*/
	public PriorDetails NG_hard;

	/**
	   Default hard-Net-to-gross std dev info {value,varies,index} (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_NG_hard;

	/**
	   Default float_fraction (if used) 
	*/
	public PriorDetails float_fraction;

	/**
	   Default sigma_float_fraction (if used) 
	*/
	public PriorDetails sigma_float_fraction;

	/**
	   Default P_oil (must be public for field access via Class.getField)
	*/
	public PriorDetails P_oil;

	/**
	   Default P_gas (must be public for field access via Class.getField)
	*/
	public PriorDetails P_gas;

	/**
	   Default P_lsg (must be public for field access via Class.getField)
	*/
	public PriorDetails P_lsg;

	/**
	   Default S_oil (must be public for field access via Class.getField)
	*/
	public PriorDetails S_oil;

	/**
	   Default S_gas (must be public for field access via Class.getField)
	*/
	public PriorDetails S_gas;

	/**
	   Default S_lsg (must be public for field access via Class.getField)
	*/
	public PriorDetails S_lsg;

	/**
	   Default sigma_S_oil (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_S_oil;

	/**
	   Default sigma_S_gas (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_S_gas;
	
	/**
	   Default sigma_S_lsg (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_S_lsg;

	/**
	  mean depth of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails depth;


	/**
	   std dev of depth of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_depth;
	/**
	   std dev of trend depth of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_depth_rock_curves = null;

	/**
	  mean depth_rock_curves of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails depth_rock_curves = null;

	/**
	  low frequency interval velocity of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails LFIV;

	/**
	  std dev of low frequency interval velocity of layer (must be public for field access via Class.getField)
	*/
	public PriorDetails sigma_LFIV;

	/**
	   in-situ reference fluid (probably brine) velocity
	*/
	public PriorDetails vp_b;
	/**
	   in-situ reference fluid (probably brine) std dev of velocity
	*/
	public PriorDetails sigma_vp_b;

	/**
	   in-situ reference fluid (probably brine) density
	*/
	public PriorDetails rho_b;
	/**
	   in-situ reference fluid (probably brine) std dev of density
	*/
	public PriorDetails sigma_rho_b;



	// variable fluids
	/**
	   Oil
	*/
	public PriorDetails vp_oil;         // JGCO2
	public PriorDetails sigma_vp_oil;
	public PriorDetails rho_oil;
	public PriorDetails sigma_rho_oil;
	public boolean vp_oil_VA,rho_oil_VA;
	/**
	   Gas
	*/
	public PriorDetails vp_gas;
	public PriorDetails sigma_vp_gas;
	public PriorDetails rho_gas;
	public PriorDetails sigma_rho_gas;
	public boolean vp_gas_VA,rho_gas_VA;
	/**
	   LSG
	*/
	public PriorDetails vp_lsg;
	public PriorDetails sigma_vp_lsg;
	public PriorDetails rho_lsg;
	public PriorDetails sigma_rho_lsg;
	public boolean vp_lsg_VA,rho_lsg_VA;

	public MetaInfo()
	{
	  /* useful (in)complete list of properties in layers ....
	  t;
	  sigma_t;
	  dz;
	  sigma_dz;
	  NG;
	  sigma_NG;
	  P_oil;
	  P_gas;
	  P_lsg;
	  S_oil;
	  S_gas;
	  S_lsg;
	  sigma_S_oil;
	  vp_oil;  // JGCO2
	  sigma_S_gas;
	  sigma_S_lsg;
	  depth;
	  depth_rock_curves;
	  LFIV;
	  sigma_depth;
	  sigma_LFIV;

	  if (*_FLOATING_GRAIN_MODEL)
	  {
	  float_fraction
	  sigma_float_fraction
	  }

	  vp_b;
	  sigma_vp_b;
	  rho_b;
	  sigma_rho_b;
	  
	  */
	}

  }
  
  void printPriorInfo(boolean isUpdated)
  {
	int i,k;

	String VS="[V]";

	if (isUpdated)
	  {
		EH.println("%nTRACE-UPDATED MODEL AND PRIORS:%n");
	  }
	else
	  {
		EH.println("%nDEFAULT MODEL AND PRIORS:%n");
	  }
	for (i=0;i<MetaInfo.length;i++)
	  {
		ModelPriorDescriptions.MetaInfo MI = MetaInfo[i];
		EH.println("Layer " + (i+1) + ":");
		if (MI.layer_name != null) 
		  EH.println("                   name  = " + MI.layer_name);
		if (MI.sand_rock_member != null) 
		  EH.println("                   reservoir member  = " + MI.sand_rock_member);
		if (MI.mud_rock_member != null) 
		  EH.println("                   non-reservoir member = " + MI.mud_rock_member);
		if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
		  {
			if (MI.hard_shale_rock_member != null) 
			  EH.println("                   non-reservoir hard-rock member = " + MI.hard_shale_rock_member);
		  }
		if (MI.brine_name != null) 
		  EH.println("                   brine        = " + MI.brine_name);
		if (MI.oil_name != null) 
		  EH.println("                   oil          = " + MI.oil_name);
		if (MI.gas_name != null) 
		  EH.println("                   gas          = " + MI.gas_name);
		if (MI.lsg_name != null) 
		  EH.println("                   low sat gas  = " + MI.lsg_name);
		if (MI.t != null) 
		  EH.println("                   t (top) = %-10.5g %s",MI.t.value, MI.t.varies ? VS : "");
		if (MI.sigma_t != null) 
		  {
			EH.println("                   sigma_t = %-10.5g %s" , MI.sigma_t.value,  MI.sigma_t.varies ? VS : "");
		  }
		EH.println("                   forced pinchout = "+!MI.isActive);
		if (MI.dz != null) 
		  EH.println("                   dz = %-10.5g %s" , MI.dz.value,  MI.dz.varies ? VS : "");
		if (MI.sigma_dz != null) 
		  EH.println("                   sigma_dz = %-10.5g %s" , MI.sigma_dz.value,  MI.sigma_dz.varies ? VS : "");
		if (MI.NG != null) 
		  EH.println("                   NG = %-10.5g %s" , MI.NG.value,  MI.NG.varies ? VS : "");
		if (MI.sigma_NG != null) 
		  EH.println("                   sigma_NG = %-10.5g %s" , MI.sigma_NG.value,  MI.sigma_NG.varies ? VS : "");

		if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
		  {
			if (MI.NG_hard != null) 
			  EH.println("                   NG_hard = %-10.5g %s" , MI.NG_hard.value,  MI.NG_hard.varies ? VS : "");
			if (MI.sigma_NG_hard != null) 
			  EH.println("                   sigma_NG_hard = %-10.5g %s" , MI.sigma_NG_hard.value,  MI.sigma_NG_hard.varies ? VS : "");
		  }
		if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			if (MI.float_fraction != null)
			  {
				EH.println("                   float_fraction = %-10.5g %s" , MI.float_fraction.value,  MI.float_fraction.varies ? VS : "");
				if (MI.float_fraction.master_layer_number != i)
				  {
					EH.println("                   --- true float_fraction master = " + 
							   MetaInfo[MI.float_fraction.master_layer_number].float_fraction.value);
				  }
			  }
			if (MI.sigma_float_fraction != null) 
			  {
				EH.println("                   sigma_float_fraction = %-10.5g %s" , MI.sigma_float_fraction.value,  MI.sigma_float_fraction.varies ? VS : "");
				if (MI.float_fraction.master_layer_number != i)
				  {
					EH.println("                   --- true sigma_float_fraction master = " + 
							   MetaInfo[MI.float_fraction.master_layer_number].sigma_float_fraction.value);
				  }
			  }
		  }






		if (MI.rho_b != null) // 200104
		  EH.println("                   rho_brine = %-10.5g %s" , MI.rho_b.value,  MI.rho_b.varies ? VS : "");
		if (MI.sigma_rho_b != null) // 200104
		  EH.println("                   sigma_rho_brine = %-10.5g %s" , MI.sigma_rho_b.value,  MI.sigma_rho_b.varies ? VS : "");

		// Fluid probabilities only meaningful in enumerations-mode...
		if (!_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			if (MI.P_oil != null) 
			  EH.println("                   P_oil = %-10.5g %s" , MI.P_oil.value,  MI.P_oil.varies ? VS : "");
			if (MI.P_gas != null) 
			  EH.println("                   P_gas = %-10.5g %s" , MI.P_gas.value,  MI.P_gas.varies ? VS : "");
			if (MI.P_lsg != null) 
			  EH.println("                   P_lsg = %-10.5g %s" , MI.P_lsg.value,  MI.P_lsg.varies ? VS : "");
		  }
		if (MI.S_oil != null) 
		  EH.println("                   S_oil = %-10.5g %s" , MI.S_oil.value,  MI.S_oil.varies ? VS : "");
		if (MI.S_gas != null) 
		  EH.println("                   S_gas = %-10.5g %s" , MI.S_gas.value,  MI.S_gas.varies ? VS : "");
		if (MI.S_lsg != null) 
		  EH.println("                   S_lsg = %-10.5g %s" , MI.S_lsg.value,  MI.S_lsg.varies ? VS : "");
		if (MI.sigma_S_oil != null) 
		  EH.println("                   sigma_S_oil = %-10.5g %s" , MI.sigma_S_oil.value,  MI.sigma_S_oil.varies ? VS : "");
		// variable fluids (JGCO2)
		if (MI.vp_oil != null) 
		  EH.println("                   vp_oil = %-10.5g %s" , MI.vp_oil.value,  MI.vp_oil.varies ? VS : "");
		if (MI.rho_oil != null) 
		  EH.println("                   rho_oil = %-10.5g %s" , MI.rho_oil.value,  MI.rho_oil.varies ? VS : "");
		if (MI.vp_gas != null) 
		  EH.println("                   vp_gas = %-10.5g %s" , MI.vp_gas.value,  MI.vp_gas.varies ? VS : "");
		if (MI.rho_gas != null) 
		  EH.println("                   rho_gas = %-10.5g %s" , MI.rho_gas.value,  MI.rho_gas.varies ? VS : "");
		if (MI.vp_lsg != null) 
		  EH.println("                   vp_lsg = %-10.5g %s" , MI.vp_lsg.value,  MI.vp_lsg.varies ? VS : "");
		if (MI.rho_lsg != null) 
		  EH.println("                   rho_lsg = %-10.5g %s" , MI.rho_lsg.value,  MI.rho_lsg.varies ? VS : "");


		if (MI.sigma_S_gas != null) 
		  EH.println("                   sigma_S_gas = %-10.5g %s" , MI.sigma_S_gas.value,  MI.sigma_S_gas.varies ? VS : "");
		if (MI.sigma_S_lsg != null) 
		  EH.println("                   sigma_S_lsg = %-10.5g %s" , MI.sigma_S_lsg.value,  MI.sigma_S_lsg.varies ? VS : "");
		if (MI.depth != null) 
		  {
			EH.println("                   depth = %-10.5g %s" , MI.depth.value,  MI.depth.varies ? VS : "");
			if (MI.depth.master_layer_number != i)
			  EH.println("                   --- true depth master = " + 
						 MetaInfo[MI.depth.master_layer_number].depth.value);
		  }
		if ((MI.depth_rock_curves != null) && (MI.depth_rock_curves != MI.depth))
		  {
			EH.println("                   depth_rock_curves = %-10.5g %s" , MI.depth_rock_curves.value,  MI.depth_rock_curves.varies ? VS : "");
			if (MI.depth_rock_curves.master_layer_number != i)
		      EH.println("                   --- true depth_rock_curves master = " + 
						 MetaInfo[MI.depth_rock_curves.master_layer_number].depth_rock_curves.value);
		  }
		if (MI.LFIV != null)
		  {
			EH.println("                   LFIV = %-10.5g %s" , MI.LFIV.value,  MI.LFIV.varies ? VS : "");
			if (MI.LFIV.master_layer_number != i)
			  {
				 EH.println("                   --- true LFIV master = " + 
						 MetaInfo[MI.LFIV.master_layer_number].LFIV.value);
			  }
		  }
		if (MI.sigma_depth != null) 
		  {
			EH.println("                   sigma_depth = %-10.5g %s" , MI.sigma_depth.value,  MI.sigma_depth.varies ? VS : "");
			if (MI.depth.master_layer_number != i)
			  EH.println("                   --- true sigma_depth master = " + 
						 MetaInfo[MI.depth.master_layer_number].sigma_depth.value);
		  }
		if (MI.sigma_LFIV != null) 
		  {
			EH.println("                   sigma_LFIV = %-10.5g %s" , MI.sigma_LFIV.value,  MI.sigma_LFIV.varies ? VS : "");
			if (MI.LFIV.master_layer_number != i)
			  {
				 EH.println("                   --- true sigma_LFIV master = " + 
						 MetaInfo[MI.LFIV.master_layer_number].sigma_LFIV.value);
			  }
		  }


		if (MI.vp_b != null)
		  EH.println("                   vp_b = %-10.5g %s" , MI.vp_b.value,  MI.vp_b.varies ? VS : "");
		if (MI.sigma_vp_b != null)
		  EH.println("                   sigma_vp_b = %-10.5g %s" , MI.sigma_vp_b.value,  MI.sigma_vp_b.varies ? VS : "");
		if (MI.rho_b != null)
		  EH.println("                   rho_b = %-10.5g %s" , MI.rho_b.value,  MI.rho_b.varies ? VS : "");
		if (MI.sigma_rho_b != null)
		  EH.println("                   sigma_rho_b = %-10.5g %s" , MI.sigma_rho_b.value,  MI.sigma_rho_b.varies ? VS : "");
	  }
	EH.println("                   t_base  = %-10.5g %s" , t_base.value,  t_base.varies ? VS : "");
	EH.println("                   sigma_t_base  = %-10.5g %s" , sigma_t_base.value,  sigma_t_base.varies ? VS : "");
  
	// Contact info
	if (_m._SI.doSimpleFluidUnitsAndContacts)
	  {
		FluidUnitsAndContacts FUC;
		FluidUnitsAndContacts.FluidUnit fu;
		FluidUnitsAndContacts.ContactInfo CI;

		EH.printf("%nPrior fluid contacts info:%n");
		FUC = _m.fluidUnitsAndContacts;
		for (i=0;i<FUC.fluidUnits.length; i++)
		  {
			fu = FUC.fluidUnits[i];
			EH.println("   Fluid Unit %3d (fluids synchronised=%s): Layer range span: %d %d",
					   (i+1), (fu.fluids_synched ? "true" : "false"),
					   fu.firstLayerIndx, fu.finalLayerIndx);
			for (k=0;k<fu.contactInfo.length; k++)
			  { 
				CI = fu.contactInfo[k];
				EH.printf("                   Contact %3d type %s%s,       depth = %-10.5g %s%n",(k+1),
						  FluidStates.FluidTypes.getFluidChar(CI.UF),FluidStates.FluidTypes.getFluidChar(CI.LF),
						  CI.depthPD.value,CI.depthPD.varies ? VS : "");

				EH.printf("                                        sigma_depth = %-10.5g %s\n",
						  CI.sigma_depthPD.value,CI.sigma_depthPD.varies ? VS : "");
			  }
		  }
	  }

	if (!isUpdated)
	  {
		EH.println("%nFluid probabilities post-ordering = %s%n",_m._SI.useFluidProbabilitiesAsPostOrdering);
	  }
  }


  /**
   Simple class defining a property's prior default value, and
   whether it is variable with each trace. Useful when parsing the xml files
   After parsing all layers, the //  VariablePropertyDetails will be
   set up to indicate whether the property varies and what index 
   to find it at
  */
  static class PriorDetails 
  {
	/**
	   Value of property (e.g. S_oil=0.2). Must keep this default at zero
	   as this default is used in subsequent tests on new PriorDetails() objects.
	*/
	double value=0.0;

	/**
	   whether it varies and will be expected in the input prior traces
	*/
	boolean varies=false;

	/**
	   whether this property is [0,1] bounded
	*/
	boolean normalised=false;

	/**
	   whether this property is synched with the layer above (LFIV, depth)
	*/
	boolean synched_above=false;
	/**
	   If this property is synched with the layer above (LFIV, depth), the
	   index of the lowest-numbered layer which it mirrors (i.e. the master layer)
	*/
	int master_layer_number;
	/** 
		what index this property has if it varies
	*/
	int index=0;
	
	public PriorDetails(boolean normalised)
	{
	  this.normalised = normalised;
	}

	public PriorDetails(double value)
	{
	  this.value = value;
	}


	void setValue(double value) throws RuntimeException
	{
	  if (value < PRIOR_TRACES_NULL_THRESHOLD)
		{
		  // do nothing. input is dummy
		  return;
		}
	  if (normalised)
		{
		  if (value<0.0) throw new RuntimeException("value < 0 for normalised property");
		  if (value>1.0) throw new RuntimeException("value >1 for normalised property");
		}
	  this.value = value;
	}

  }


  /**
   Simple class defining a property's variablity with each trace,
   and the index to expect it at 
   Useful when parsing the xml files
  */
  protected static class VariablePropertyDetails 
  {
	/**
	   whether it varies and will be expected in the input prior traces
	*/
	boolean varies=false;
	/**
	   what index it ought to have in the prior traces if it varies
	*/
	int index=0;
  }

  /** 
	A useful object for the ArrayList we use to get variable names
  */
  public static final class propertyInfo
  {
	int layer_no;
	PriorDetails PD;

	public propertyInfo(int layer_no,PriorDetails PD)
	{
	  this.layer_no = layer_no;
	  this.PD = PD;
	}

	public propertyInfo() {}
	
  }
}

