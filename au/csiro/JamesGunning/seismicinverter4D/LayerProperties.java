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
import static au.csiro.JamesGunning.lib.util.MathUtil.*;
import static java.lang.Math.*;
import java.util.*;

/** 
	Set of fundamental properties required at each 
	layer for complete description of rock physics
*/

public final class LayerProperties
{
  /**
	 Link to underlying model
  */
  Model _m;

  byte vintage=0;

  /** 
	  time of layer top
  */
  Property t;

  /** 
	  time of layer base. Not used except for last layer
  */
  Property t_base = null;

  /** 
	  depth of layer
  */
  Property d;
  /** 
	  rock-trends (curves) depth of layer
  */
  Property d_curves;

  /**
	 Low freq interval velocity
  */
  Property LFIV;

  /**
	 Net-to-gross (fraction of "sand", or 
	 `sand-end' rock member, susceptible of fluid substitution)
  */
  Property NG;
  /**
	 sand porosity
  */
  Property phi_s;
  /**
	 mud density ("mud" is shale/impermeable rock)
  */
  Property rho_m;
  /**
	 pure sand p-wave velocity
  */
  Property vp_s;
  /**
	 pure sand s-wave velocity
  */
  Property vs_s;

  /**
	 pure mud p-wave velocity
  */
  Property vp_m;
  /**
	 pure mud s-wave velocity
  */
  Property vs_m;

  /**
	 Hard net-to-gross (fraction of hard rock in overall impermeable volume fraction)
  */
  Property NG_hard;
  /**
	 hard rock p-wave velocity
  */
  Property vp_hard;
  /**
	 hard rock  s-wave velocity
  */
  Property vs_hard;
  /**
	 hardrock density 
  */
  Property rho_hard;


  /**
	 Floating grain volume fraction of reservoir rock
  */
  Property float_fraction;

  /** 
	  Current hydrocarbon-fluid (o) saturation
  */
  Property S_oil;
  /** 
	  Current hydrocarbon-fluid (l) saturation
  */
  Property S_lsg;
  /** 
	  Current hydrocarbon-fluid (g) saturation
  */
  Property S_gas;


  /** 
	  P-wave speed of current oil
  */
  Property vp_oil;
  /** 
	  P-wave speed of current lsg
  */
  Property vp_lsg;
  /** 
	  P-wave speed of current gas
  */
  Property vp_gas;

  /** 
	  Density of oil
  */
  Property rho_oil;
  /** 
	  Density of lsg
  */
  Property rho_lsg;
  /** 
	  Density of gas
  */
  Property rho_gas;


  /** 
	  P-wave speed of new "brine"
  */
  Property vp_b;
 
  /** 
	  Density of new "brine"
  */
  Property rho_b;

  /**
	 Reference to regression set, which has grain and reference fluid 
	 info as member variables
	 K_g, rho_g, vp_b_ref, rho_b_ref
  */
  RockandFluidData.RegressionSet RFD_PermeableRegSet=null;

  /**
	 Hydrocarbon saturation (computed by and accessible after a 
	 computeEffectiveProperties(fluid_type) call). Sensible if the 
	 fluid_type is not brine (0)
  */
  double S_h=0;

  /**
	 Hydrocarbon density (computed by and accessible after a 
	 computeEffectiveProperties(fluid_type) call. Sensible if the 
	 fluid_type is not brine (0)
  */
  double rho_h=0;

  /**
	 Effective Oil saturation. Thickness weighted average if there are 
	 sublayers when working with contacts, so overall volume is preserved.
  */
  double S_oil_eff=0;

  /**
	 Effective Gas saturation. Thickness weighted average if there are 
	 sublayers when working with contacts, so overall volume is preserved.
  */
  double S_gas_eff=0;

  /**
	 Effective lsg saturation. Thickness weighted average if there are 
	 sublayers when working with contacts, so overall volume is preserved.
  */
  double S_lsg_eff=0;

  /**
	 Things to be computed
  */
  double rho_eff;
  double vp_eff;
  double vs_eff;

  /**
	 Useful depths
  */
  double suspended_depth_top;
  double suspended_depth_bottom;

  /**
	 indx of this LayerProperty in the array
  */
  int thisLayerIndx;

  /** 
	  Useful map of all properties we need, which will have string keys
  */
  HashMap<String,Property> completePropertyNameMap=null;

  /**
	 Possibly unused class with internal fluid-contact info
  */
  InternalFluidContactInfo fcInfo = null;

 /**
	Possibly unused class with 4D seismic info
  */
  Info4D info4D = null;

  public LayerProperties(Model m, int thisLayerIndx, boolean needsTBase,byte vintage)
  {
	this._m = m;
	this.thisLayerIndx = thisLayerIndx;
	this.vintage = vintage;                // zero-offset index
	t = new Property(0,false);
	if (needsTBase) t_base = new Property(0,false);
	d = new Property(0,false);
	d_curves = new Property(0,false);
	LFIV = new Property(0,false);
	NG = new Property(0,true);
	phi_s = new Property(0,true);
	rho_m = new Property(0,false);
	vp_s = new Property(0,false);
	vs_s = new Property(0,false);
	vp_m = new Property(0,false);
	vs_m = new Property(0,false);
	vp_b = new Property(0,false);
	rho_b = new Property(0,false);

	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		NG_hard = new Property(0,true);
		vp_hard = new Property(0,false);
		vs_hard = new Property(0,false);
		rho_hard = new Property(0,false);
	  }

	if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		float_fraction = new Property(0,true);
	  }

	// the o/l/g properties are fake Properties 
	// which wrap the hydrocarbon 'h' value
	S_oil = new Property(0,true);
	S_lsg = new Property(0,true);
	S_gas = new Property(0,true);

	vp_oil = new Property(0,false);
	vp_lsg = new Property(0,false);
	vp_gas = new Property(0,false);

	rho_oil = new Property(0,false);
	rho_lsg = new Property(0,false);
	rho_gas = new Property(0,false);

	

	// get internal map for ModelFormatting
	assembleStringKeyedHashMap();
	// and reverse map for other uses e.g. model-vector feedback
	addPropertiesToHashMap(_m.completePropertyNameMap);
  }

  // form local HashMap, based on string keys, to enable model formatting as per ModelFormat
  private void assembleStringKeyedHashMap()
  {
	completePropertyNameMap = new HashMap<String,Property>();

	completePropertyNameMap.put("t",t);
	if (t_base != null) completePropertyNameMap.put("t_base",t_base);
	completePropertyNameMap.put("d",d);
	completePropertyNameMap.put("LFIV",LFIV);
	completePropertyNameMap.put("NG",NG);
	completePropertyNameMap.put("phi_s",phi_s);
	completePropertyNameMap.put("rho_m",rho_m);
	completePropertyNameMap.put("vp_s",vp_s);
	completePropertyNameMap.put("vs_s",vs_s);
	completePropertyNameMap.put("vp_m",vp_m);
	completePropertyNameMap.put("vs_m",vs_m);
	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		completePropertyNameMap.put("NG_hard",NG_hard);
		completePropertyNameMap.put("vp_hard",vp_hard);
		completePropertyNameMap.put("vs_hard",vs_hard);
		completePropertyNameMap.put("rho_hard",rho_hard);
	  }
	if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		completePropertyNameMap.put("float_fraction",float_fraction);
	  }
  }

  /**
	 External routines will want to construct a map of the Properties objects in this
	 class. Supply a HashMap to do this, after the object is instantiated
  */
  private void addPropertiesToHashMap(HashMap map)
  {
	map.put(t,"t");
	if (t_base != null) map.put(t_base,"t_base");
	map.put(d,"d");
	map.put(LFIV,"LFIV");
	map.put(NG,"NG");
	map.put(phi_s,"phi_s");
	map.put(rho_m,"rho_m");
	map.put(vp_s,"vp_s");
	map.put(vs_s,"vs_s");
	map.put(vp_m,"vp_m");
	map.put(vs_m,"vs_m");
	map.put(vp_b,"vp_b");
	map.put(rho_b,"rho_b");

	// Add suffixes to later vintage saturations etc to make them distinguishable
	String V = ((!_m.is4D && vintage==0)) ? "" : (","+vintage);
	map.put(S_oil,"S_oil"+V);	
	map.put(S_lsg,"S_lsg"+V);
	map.put(S_gas,"S_gas"+V);

	map.put(vp_oil,"vp_oil"+V);
	map.put(vp_lsg,"vp_lsg"+V);
	map.put(vp_gas,"vp_gas"+V);

	map.put(rho_oil,"rho_oil"+V);
	map.put(rho_lsg,"rho_lsg"+V);
	map.put(rho_gas,"rho_gas"+V);

	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		map.put(NG_hard,"NG_hard");
		map.put(vp_hard,"vp_hard");
		map.put(vs_hard,"vs_hard");
		map.put(rho_hard,"rho_hard");
	  }
	if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		map.put(float_fraction,"float_fraction");
	  }

	// these needed for dumpRealisations() in SeismicInverter, as per ModelFormat
  }

  void setUpFluidContactInfo(boolean isPermeable)
  {	
    if (isPermeable)
	  {
		// use fuller class for permeable layers
		fcInfo = new InternalPermeableFluidContactInfo();
	  }
	else
	  {
		// simple storage only needed for shales
		fcInfo = new InternalFluidContactInfo();
	  }
  }

  void init4DInfo()
  {	
	info4D = new Info4D();	
  }

  /**
	 Gassman substitution master. Forks contact model if needed 
  */
  boolean computeEffectiveProperties(short fluid_type, boolean showDiagnostics,double t_below, 
									 double thickness, double suspended_depth_hook, boolean descending)
  {
	if (!_m.MPD.MetaInfo[thisLayerIndx].isActive)
	  {
		setDefaultInactiveProperties(suspended_depth_hook);
		return true;
	  }
	if (!_m._SI.doSimpleFluidUnitsAndContacts)
	  {
		return computeEffectivePropertiesSansContacts(fluid_type,showDiagnostics,t_below,
													  thickness,suspended_depth_hook,descending);
	  }
	else
	  {
		return computeEffectivePropertiesUsingContacts(showDiagnostics,t_below,
													   thickness,suspended_depth_hook,descending);
	  }
  }

  /** Inactive layer safe defaults 
   */
  void setDefaultInactiveProperties(double suspended_depth_hook)
  {
	suspended_depth_top = suspended_depth_bottom = suspended_depth_hook; // passed in
  }

  /**
	 Gassman substitution for layer, contacts-free model 
  */
  boolean computeEffectivePropertiesSansContacts(short fluid_type, boolean showDiagnostics,
												 double t_below, double thickness, 
												 double suspended_depth_hook, boolean descending)
  {
	double K_b=0,K_fl=0,K_h,vp_h=0;
	double rho_satb,M_satb,mu_satb,K_satb;
	double K_satfl=0,rho_satfl=0,M_satfl=0,mu_satfl;
	double M_m,mu_m;
	double M_eff=0,mu_eff;
	double K_g=0,rho_g;
	double K_b_new;
	double S_L=0;
	double[] vp_vs_rho = new double[3];

	if (fluid_type<0 || fluid_type >FluidStates.FluidTypes.Nfluids-1)
	  {
		throw new RuntimeException(" Fluid type must be one of 0,1,2..." + 
								   (FluidStates.FluidTypes.Nfluids - 1));
	  }

	if (NG.value == 0.0)
	  {
		// it's pure impermeable rock....
		computeImpermeableProps(vp_vs_rho);
		vp_eff = vp_vs_rho[0];
		vs_eff = vp_vs_rho[1];
		rho_eff = vp_vs_rho[2];
		S_h = S_oil_eff = S_gas_eff = S_lsg_eff = 0;
	  }
	else
	  {
		if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			// remap regression properties if float fraction requested was negative
			// this prevent wierd smearing of the distribution onto the plane float_fraction=0
			// caused by the truncation. vp and vs are projected back to the plane parallel
			// to the rock properties {phi_s,vp_s,vs_s} principal component axes
			computeNegativeFloatRemappedLayerProperties();			
		  }

		// set grain props from appropriate rock
		K_g = RFD_PermeableRegSet.K_g;	
		rho_g = RFD_PermeableRegSet.rho_g;	

		// ref brine + matrix density
		// rho_satb = phi_s.value * rho_b.value + (1.0-phi_s.value) * rho_g;

		rho_satb = phi_s.value * RFD_PermeableRegSet.rho_b_ref + (1.0-phi_s.value) * rho_g;
		// moduli of ref-fluid saturated sand-member rock
		M_satb = rho_satb *  sqr(vp_s.value);
		mu_satb =  rho_satb *  sqr(vs_s.value);
		K_satb = M_satb - 1.3333333333333 * mu_satb;

		if (showDiagnostics)
		  {
			_m._SI.EH.println("phi = " + phi_s.value);
			_m._SI.EH.println("vp_s(satb) = " + vp_s.value);
			_m._SI.EH.println("vs_s(satb) = " + vs_s.value);
			_m._SI.EH.println("K_satb = " + K_satb);
			_m._SI.EH.println("rho_satb = " + rho_satb);
			_m._SI.EH.println("rho_b = " + rho_b.value);
			_m._SI.EH.println("vp_b = " + vp_b.value);
		  }
		if (K_satb<0 && phi_s.value!= 0.0) 
		  {			
			return false;
		  }

		mu_satfl = mu_satb;

		// ref fluid bulk mod
		// K_b = rho_b.value * sqr(vp_b.value);		
		K_b = RFD_PermeableRegSet.K_b_ref; 

		if (fluid_type != FluidStates.FluidTypes.brine)
		  {
			if (showDiagnostics)
			  {
				_m._SI.EH.println("Implied Kdry = " + 
								  ((phi_s.value==0.0)? K_g : K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 
															1.0/(((K_g/K_b)-1.0)*phi_s.value)))));
			  }

			if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0) && phi_s.value!=0.0)
			  {
				// Kdry has no positive solution. Reject this				
				if (showDiagnostics)
				  {
					_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
									  " under these conditions (N.B. can be large at small porosities)." +
									  " Rock speeds may be too low in trend curves");
					_m._SI.EH.printf("*Possible* required minimum vp,vs boost as function of porosity...%nphi       boost%n");
					for (double phi_t=0.02; phi_t<0.4; phi_t+=.02)
					  {
						double alpha = Math.sqrt(K_g/((phi_t*((K_g/K_b)-1.0)+1.0) * (sqr(vp_s.value) - (4.0/3.0)* sqr(vs_s.value))
													  * (phi_t * RFD_PermeableRegSet.rho_b_ref + (1.0-phi_t) * rho_g)));
						_m._SI.EH.printf("%g %g%n",phi_t,alpha);
					  }
				  }
				return false;
			  }

			switch (fluid_type)
			  {
			  case FluidStates.FluidTypes.oil:
				S_h = S_oil.value;
				S_oil_eff = S_h; S_gas_eff = S_lsg_eff = 0;
				rho_h = rho_oil.value;
				vp_h = vp_oil.value;
				break;

			  case FluidStates.FluidTypes.lsg:
				S_h = S_lsg.value;
				S_lsg_eff = S_h; S_oil_eff = S_gas_eff = 0;
				rho_h = rho_lsg.value;
				vp_h = vp_lsg.value;
				break;

			  case FluidStates.FluidTypes.gas:
				S_h = S_gas.value;
				S_gas_eff = S_h; S_oil_eff = S_lsg_eff = 0;
				rho_h = rho_gas.value;
				vp_h = vp_gas.value;
				break;			
			  case FluidStates.FluidTypes.brine:
				break; // never happens
			  }

			// hydrocarbon bulk mod.
			K_h = rho_h * sqr(vp_h);
			// new brine bulk mod
			K_b_new = rho_b.value * sqr(vp_b.value);

			// New effective fluid bulk modulus & density...
			// A mixture of hydrocarbon and new brine
			K_fl = 1.0/(S_h/K_h + (1.0-S_h)/K_b_new);
			
			// simple sanity test
			/*
			  if ((1.0/((K_g/K_satb)-1.0) - 1.0/((K_g/K_b)-1.0)*phi_s.value) < 0)
			  {
			  System.out.println("Crazy Gassman brine substitution");
			  return false;	
			  }
			*/
			// do fluid substitution
			if (phi_s.value !=0)
			  {
				K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s.value // FD
										+ 1.0/((K_g/K_satb)-1.0)));
			  }
			else
			  {
				K_satfl = K_g;
			  }
			if (showDiagnostics)
			  {
				_m._SI.EH.println("Hydrocarbon substitution details:"); // FD
				_m._SI.EH.println("K_b = " + K_b);
				_m._SI.EH.println("K_h = " + K_h);
				_m._SI.EH.println("S_h = " + S_h);
				_m._SI.EH.println("K_fl = " + K_fl);
				_m._SI.EH.println("K_satfl = " + K_satfl);
			  }
			if (K_satfl<0) // FD
			  {
				return false;			
			  }
			// new effective density
			rho_satfl = (1.0-phi_s.value) * rho_g + phi_s.value * (S_h*rho_h+(1.0-S_h)*rho_b.value); // FD
			// new effective p-wave modulus
			M_satfl = K_satfl + 1.3333333333333 * mu_satfl; // FD
		  }
		else
		  {
			S_h = S_oil_eff = S_gas_eff = S_lsg_eff = 0;
			// there are no hydrocarbons, Substitute pure new brine
			if ((rho_b.value != RFD_PermeableRegSet.rho_b_ref) || (vp_b.value != RFD_PermeableRegSet.vp_b_ref))
			  {
				// use pure new brine ONLY if different
				if (showDiagnostics)
				  {
					_m._SI.EH.println("Implied Kdry = " + 
									  ((phi_s.value==0.0)? K_g : K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 
																			   1.0/(((K_g/K_b)-1.0)*phi_s.value)))));
				  }
				
				if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0) && phi_s.value!=0.0)
				  {
					// Kdry has no positive solution. Reject this				
					if (showDiagnostics)
					  {
						_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
										  " under these conditions (N.B. can be large at small porosities)." +
										  " Rock speeds may be too low in trend curves");
					  }
					return false;
				  }
				// New fluid h is pure new brine
				// New effective fluid bulk modulus & density...
				K_b_new = rho_b.value * sqr(vp_b.value);

				// New effective fluid bulk modulus & density...
				K_fl = K_b_new;

				// do fluid substitution
				if (phi_s.value !=0)
				  {
					K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s.value 
											+ 1.0/((K_g/K_satb)-1.0)));
				  }
				else
				  {
					K_satfl = K_g;
				  }
				if (showDiagnostics)
				  {
					_m._SI.EH.println("New S=1 brine substitution details:");
					_m._SI.EH.println("K_b_ref = " + K_b);
					_m._SI.EH.println("K_b_new = " + K_fl);
					_m._SI.EH.println("K_satfl = " + K_satfl);
				  }
				if (K_satfl<0) 
				  {
					return false;			
				  }
				// new effective density
				rho_satfl = (1.0-phi_s.value) * rho_g + phi_s.value * rho_b.value;
				// new effective p-wave modulus
				M_satfl = K_satfl + 1.3333333333333 * mu_satfl;
			  }
			else
			  {
				// New brine same as ref brine, no hydrocarbs. No need to compute Gassman
				K_satfl = K_satb;
				rho_satfl = rho_satb;
				M_satfl = M_satb;

				// sanity test
				if (showDiagnostics)
				  {
					_m._SI.EH.println("Implied Kdry = " + K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 1.0/(((K_g/K_b)-1.0)*phi_s.value))));
				  }
				if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0))
				  {
					// We're not fluid-substituting, but Kdry has no positive solution. Reject this
					if (showDiagnostics)
					  {
						_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
										  " under these conditions (N.B. can be large at small porosities)." +
										  " Rock speeds may be too low in trend curves");
					  }
					return false;
				  }
			  }
		  }



		// shale moduli
		computeImpermeableProps(vp_vs_rho);  // stuffs equivalent impermeable rocks props into vp_vs_rho[]={vp,vs,rho}
		mu_m = vp_vs_rho[2] * sqr(vp_vs_rho[1]);
		M_m = vp_vs_rho[2] * sqr(vp_vs_rho[0]);
		// NG weighted average for densities
		rho_eff = NG.value * rho_satfl + (1.0-NG.value)*vp_vs_rho[2];
		// And moduli are averaged using Backus
		if (NG.value<1)
		  {
			// Backus average of p-wave modulus.
			M_eff = 1.0/(NG.value/M_satfl + (1.0-NG.value)/M_m);
			// Backus average of shear modulus
			mu_eff = 1.0/(NG.value/mu_satfl + (1.0-NG.value)/mu_m);
		  }
		else 
		  {
			M_eff = M_satfl;
			mu_eff = mu_satfl;
		  }
		// and resulting equivalent velocities
		vp_eff = Math.sqrt(M_eff/rho_eff);
		vs_eff = Math.sqrt(mu_eff/rho_eff);
	  }

	// Update suspended depths
	updateSuspendedDepths(descending,suspended_depth_hook,t_below, thickness);

	return true;
  }


 /**
	 Gassman substitution for layer,
	 using on-the-fly contact-determined fluid content
  */
  boolean computeEffectivePropertiesUsingContacts(boolean showDiagnostics,
												  double t_below, double thickness, 
												  double suspended_depth_hook, boolean descending)
  {
	int i;
	double K_b=0,K_fl=0,K_h,vp_h=0;
	double rho_satb,M_satb,mu_satb,K_satb;
	double K_satfl=0,rho_satfl=0,M_satfl=0,mu_satfl;
	double M_m,mu_m;
	double M_eff=0,mu_eff;
	double K_g=0,rho_g;
	double K_b_new;
	double S_L=0;
	short fluid_type = FluidStates.FluidTypes.brine;  // default

	boolean exhaustedContacts=true;
	double z_find,t_find,t_contact,z_top=0;
	double[] zBracket = new double[2];
	InternalPermeableFluidContactInfo fcpInfo = null;
	boolean isPermeable=false;

	// _m.debugCounter++;

	// First, set hooked depth.
	if (descending)
	  {
		// going down
		suspended_depth_top = suspended_depth_hook; // passed in
	  }
	else
	  {
		// Going up: 
		suspended_depth_bottom = suspended_depth_hook;					
	  }	

	// set to no contacts as default;
	if (fcInfo!=null) 
	  {
		fcInfo.n_internal_sublayers = 1;
	  }
	ModelPriorDescriptions.MetaInfo _MI = _m.MPD.MetaInfo[thisLayerIndx];
	if (_MI.is_permeable)
	  {
		fcpInfo = (InternalPermeableFluidContactInfo)fcInfo;
		isPermeable=true;
	  }
	S_h = S_oil_eff = S_gas_eff = S_lsg_eff = 0;
	if (NG.value == 0.0)
	  {
		// it's pure impermeable rock....
		vp_eff = vp_m.value;
		vs_eff = vs_m.value;
		rho_eff = rho_m.value;
		if (vp_eff<=0 || vs_eff<=0 || rho_eff<=0) return false;
		// code to track contact position within shales (no effect on likelihood)
		if (descending)
		  {
			z_find = suspended_depth_top;
			t_find = t.value;  // expected to be set at correct value, all vintages
		  }
		else
		  {
			z_find = suspended_depth_bottom;
			t_find = t_below;			
		  }

		// loop over internal contact-separated intervals (maximum of 4!)
		for (i=0;i<4;i++)
		  {
			fcInfo.t_contact[i] = t_find;
			fluid_type = getFluidTypeAt(z_find,zBracket,descending);
			if (fluid_type == -1) return false; // Error. wierdness

			// test for contact within-layer
			// compute thicknesses based on vp_eff
			if (descending)
			  {
				double z_bottom;
				// test depth, assuming fluid holds down to t_below_or_thickness, 
				// (or bottom of layer from previous vintage)
				if (vintage==0) 
				  {
					z_bottom = z_find + 0.5*(t_below - t_find) * vp_eff;
				  }
				else
				  {
					// later vintages from preserved thickness
					z_bottom = suspended_depth_top + thickness; 
				  }
				if (z_bottom < zBracket[1])
				  {
					// layer ends before contact. Simple pure one-fluid layer
					exhaustedContacts = true;
					if (vintage==0) 
					  {
						fcInfo.t_contact[i+1] = t_below;
					  }
					else
					  {
						fcInfo.t_contact[i+1] = t_find + 2.0*(z_bottom - z_find)/vp_eff;
					  }
				  }
				else
				  {
					exhaustedContacts = false;
					// Solve for contact time accordingly:
					// zBracket[1] = z_find + 0.5*(t_contact - t_find) * vp_eff
					t_contact = t_find + 2.0*(zBracket[1] - z_find)/vp_eff;
					// set new time
					t_find = t_contact;
					// set new depth at encountered contact
					z_find = zBracket[1];
				  } 
			  }
			else
			  {
				// test depth, assuming fluid holds up to t.value
				if (vintage==0) 
				  {
					z_top = z_find - 0.5*(t_find - t.value) * vp_eff;
				  }
				else 
				  {
					// later vintages from preserved thickness
					z_top = suspended_depth_bottom - thickness; 
				  }
				if (z_top > zBracket[0])
				  {
					// layer ends before encountering contact. Simple pure one-fluid layer
					exhaustedContacts = true;
					if (vintage==0) 
					  {
						fcInfo.t_contact[i+1] = t.value;
					  }
					else
					  {
						fcInfo.t_contact[i+1] = t_find - 2.0*(z_find -z_top)/vp_eff;					  
					  }
				  }
				else
				  {
					exhaustedContacts = false;
					// Solve for contact time accordingly:
					// zBracket[1] = z_find + 0.5*(t_contact - t_find) * vp_eff
					t_contact = t_find - 2.0*(z_find - zBracket[0])/vp_eff;
					// set new time
					t_find = t_contact;
					// set new depth at encountered contact
					z_find = zBracket[0]; 
				  }
			  }
			if (isPermeable)
			  {
				// JG FIX 6Aug 2016
				int l = fcpInfo.n_internal_sublayers-1;
				fcpInfo.fluid_type[l] = fluid_type;
				fcpInfo.rho_eff[l] = rho_eff;
				fcpInfo.vp_eff[l] = vp_eff;
				fcpInfo.vs_eff[l] = vs_eff;
				fcpInfo.S_h[l] = 0;
				fcpInfo.rho_h[l] = 1.0;// placeholder - no meaning if S_h=0
			  }
			// escape when layer has accounted for all internal fluid contacts
			if (!exhaustedContacts) 
			  {
				fcInfo.n_internal_sublayers++;
			  }
			else
			  {
				break;
			  }		
		  } // possible internal fluids
 
			// flip times array into useful top-down order...
		if (!descending) fcInfo.reverseTimes();
	  } // NG==0
	else
	  {
		//fcpInfo = (InternalPermeableFluidContactInfo)fcInfo;		
	
		if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			// remap regression properties if float fraction requested was negative
			// this prevent wierd smearing of the distribution onto the plane float_fraction=0
			// caused by the truncation. vp and vs are projected back to the plane parallel
			// to the rock properties {phi_s,vp_s,vs_s} principal component axes
			computeNegativeFloatRemappedLayerProperties();			
		  }

		// set grain props from appropriate rock
		K_g = RFD_PermeableRegSet.K_g;	
		rho_g = RFD_PermeableRegSet.rho_g;	

		// ref brine + matrix density
		// rho_satb = phi_s.value * rho_b.value + (1.0-phi_s.value) * rho_g;


		rho_satb = phi_s.value * RFD_PermeableRegSet.rho_b_ref + (1.0-phi_s.value) * rho_g;
		// moduli of ref-fluid saturated sand-member rock
		M_satb = rho_satb *  sqr(vp_s.value);
		mu_satb =  rho_satb *  sqr(vs_s.value);
		K_satb = M_satb - 1.3333333333333 * mu_satb;

		if (showDiagnostics)
		  {
			_m._SI.EH.println("phi = " + phi_s.value);
			_m._SI.EH.println("vp_s(satb) = " + vp_s.value);
			_m._SI.EH.println("vs_s(satb) = " + vs_s.value);
			_m._SI.EH.println("K_satb = " + K_satb);
			_m._SI.EH.println("rho_satb = " + rho_satb);
			_m._SI.EH.println("rho_b = " + rho_b.value);
			_m._SI.EH.println("vp_b = " + vp_b.value);
		  }
		if (K_satb<0 && phi_s.value!= 0.0) 
		  {			
			return false;
		  }

		mu_satfl = mu_satb;

		// ref fluid bulk mod
		// K_b = rho_b.value * sqr(vp_b.value);		
		K_b = RFD_PermeableRegSet.K_b_ref;


		// shale moduli - constant for all internal fluid contacts...
		mu_m = rho_m.value * sqr(vs_m.value);
		M_m = rho_m.value * sqr(vp_m.value);

		// Determine fluid type from contacts

		if (descending)
		  {
			z_find = suspended_depth_top;
			t_find = t.value;
		  }
		else
		  {
			z_find = suspended_depth_bottom;
			t_find = t_below;
		  }
		// loop over internal contact-separated intervals (maximum of 4!)
		for (i=0;i<4;i++)
		  {
			fcpInfo.t_contact[i] = t_find;
			fluid_type = getFluidTypeAt(z_find,zBracket,descending);
			if (fluid_type == -1) return false; // Error. wierdness

			if (fluid_type != FluidStates.FluidTypes.brine)
			  {
				if (showDiagnostics)
				  {
					_m._SI.EH.println("Implied Kdry = " + 
									  ((phi_s.value==0.0)? K_g : K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 
																			   1.0/(((K_g/K_b)-1.0)*phi_s.value)))));
				  }

				if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0) && phi_s.value!=0.0)
				  {
					// Kdry has no positive solution. Reject this				
					if (showDiagnostics)
					  {
						_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
										  " under these conditions (N.B. can be large at small porosities)." +
										  " Rock speeds may be too low in trend curves");
					  }
					return false;
				  }

				switch (fluid_type)
				  {
				  case FluidStates.FluidTypes.oil:
					S_h = S_oil.value;
					S_oil_eff = S_h;
					rho_h = rho_oil.value;
					vp_h = vp_oil.value;
					break;

				  case FluidStates.FluidTypes.lsg:
					S_h = S_lsg.value;
					S_lsg_eff = S_h;
					rho_h = rho_lsg.value;
					vp_h = vp_lsg.value;
					break;

				  case FluidStates.FluidTypes.gas:
					S_h = S_gas.value;
					S_gas_eff = S_h;
					rho_h = rho_gas.value;
					vp_h = vp_gas.value;
					break;			
				  case FluidStates.FluidTypes.brine:
					// never happens
					break; 
				  }

				// hydrocarbon bulk mod.
				K_h = rho_h * sqr(vp_h);
				// new brine bulk mod
				K_b_new = rho_b.value * sqr(vp_b.value);

				// New effective fluid bulk modulus & density...
				// A mixture of hydrocarbon and new brine
				K_fl = 1.0/(S_h/K_h + (1.0-S_h)/K_b_new);
			
				// do fluid substitution
				if (phi_s.value !=0)
				  {
					K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s.value // FD
											+ 1.0/((K_g/K_satb)-1.0)));
				  }
				else
				  {
					K_satfl = K_g;
				  }
				if (showDiagnostics)
				  {
					_m._SI.EH.println("Hydrocarbon substitution details:"); // FD
					_m._SI.EH.println("K_b = " + K_b);
					_m._SI.EH.println("K_h = " + K_h);
					_m._SI.EH.println("S_h = " + S_h);
					_m._SI.EH.println("K_fl = " + K_fl);
					_m._SI.EH.println("K_satfl = " + K_satfl);
				  }
				if (K_satfl<0) // FD
				  {
					return false;			
				  }
				// new effective density
				rho_satfl = (1.0-phi_s.value) * rho_g + phi_s.value * (S_h*rho_h+(1.0-S_h)*rho_b.value); // FD

				// new effective p-wave modulus
				M_satfl = K_satfl + 1.3333333333333 * mu_satfl; // FD
			  }
			else
			  {
				// there are no hydrocarbons, Substitute pure new brine
				// Set S_h in case we have sublayers and they are used for averaging
				// use "new brine" if so..
				S_h=0;

				if ((rho_b.value != RFD_PermeableRegSet.rho_b_ref) || (vp_b.value != RFD_PermeableRegSet.vp_b_ref))
				  {
					// use pure new brine ONLY if different
					if (showDiagnostics)
					  {
						_m._SI.EH.println("Implied Kdry = " + 
										  ((phi_s.value==0.0)? K_g : K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 
																				   1.0/(((K_g/K_b)-1.0)*phi_s.value)))));
					  }
				
					if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0) && phi_s.value!=0.0)
					  {
						// Kdry has no positive solution. Reject this				
						if (showDiagnostics)
						  {
							_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
											  " under these conditions (N.B. can be large at small porosities)." +
											  " Rock speeds may be too low in trend curves");
						  }
						return false;
					  }
					// New fluid h is pure new brine
					// New effective fluid bulk modulus & density...
					K_b_new = rho_b.value * sqr(vp_b.value);

					// New effective fluid bulk modulus & density...
					K_fl = K_b_new;

					// do fluid substitution
					if (phi_s.value !=0)
					  {
						K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s.value 
												+ 1.0/((K_g/K_satb)-1.0)));
					  }
					else
					  {
						K_satfl = K_g;
					  }
					if (showDiagnostics)
					  {
						_m._SI.EH.println("New S=1 brine substitution details:");
						_m._SI.EH.println("K_b_ref = " + K_b);
						_m._SI.EH.println("K_b_new = " + K_fl);
						_m._SI.EH.println("K_satfl = " + K_satfl);
					  }
					if (K_satfl<0) 
					  {
						return false;			
					  }
					// new effective density
					rho_satfl = (1.0-phi_s.value) * rho_g + phi_s.value * rho_b.value;
					// new effective p-wave modulus
					M_satfl = K_satfl + 1.3333333333333 * mu_satfl;
				  }
				else
				  {
					// New brine same as ref brine, no hydrocarbs. No need to compute Gassman
					K_satfl = K_satb;
					rho_satfl = rho_satb;
					M_satfl = M_satb;

					// sanity test
					if (showDiagnostics)
					  {
						_m._SI.EH.println("Implied Kdry = " + K_g/(1.0+1.0/(1.0/((K_g/K_satb)-1.0) - 1.0/(((K_g/K_b)-1.0)*phi_s.value))));
					  }
					if (K_satb < K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0))
					  {
						// We're not fluid-substituting, but Kdry has no positive solution. Reject this
						if (showDiagnostics)
						  {
							_m._SI.EH.println("K_satb must be at least= " + (K_g/(phi_s.value * (K_g/K_b - 1.0) +1.0)) + 
											  " under these conditions (N.B. can be large at small porosities)." +
											  " Rock speeds may be too low in trend curves");
						  }
						return false;
					  }
				  }
			  }

			// NG weighted average for densities
			rho_eff = NG.value * rho_satfl + (1.0-NG.value)*rho_m.value;
			// And moduli are averaged using Backus
			if (NG.value<1)
			  {
				// Backus average of p-wave modulus.
				M_eff = 1.0/(NG.value/M_satfl + (1.0-NG.value)/M_m);
				// Backus average of shear modulus
				mu_eff = 1.0/(NG.value/mu_satfl + (1.0-NG.value)/mu_m);
			  }
			else 
			  {
				M_eff = M_satfl;
				mu_eff = mu_satfl;
			  }
			// and resulting equivalent velocities
			vp_eff = Math.sqrt(M_eff/rho_eff);
			vs_eff = Math.sqrt(mu_eff/rho_eff);

			// test for contact within-layer
			// compute thicknesses based on vp_eff
			if (descending)
			  {
				double z_bottom;
				// test depth, assuming fluid holds down to t_below_or_thickness
				if (vintage==0) 
				  {
					z_bottom = z_find + 0.5*(t_below - t_find) * vp_eff;
				  }
				else
				  {
					// later vintages from preserved thickness
					z_bottom = suspended_depth_top + thickness; 
				  }
				if (z_bottom < zBracket[1])
				  {
					// layer ends before contact. Simple pure one-fluid layer
					exhaustedContacts = true;
					if (vintage==0) 
					  {
						fcpInfo.t_contact[i+1] = t_below;
					  }
					else
					  {
						fcInfo.t_contact[i+1] = t_find + 2.0*(z_bottom - z_find)/vp_eff;
					  }
				  }
				else
				  {
					exhaustedContacts = false;
					// Solve for contact time accordingly:
					// zBracket[1] = z_find + 0.5*(t_contact - t_find) * vp_eff
					t_contact = t_find + 2.0*(zBracket[1] - z_find)/vp_eff;
					// set new time
					t_find = t_contact;
					// set new depth at encountered contact
					z_find = zBracket[1];
				  }
			  }
			else
			  {
				// test depth, assuming fluid holds up to t.value
				if (vintage==0) 
				  {
					z_top = z_find - 0.5*(t_find - t.value) * vp_eff;
				  }
				else
				  {
					// later vintages from preserved thickness
					z_top = suspended_depth_bottom - thickness; 				  
				  }
				if (z_top > zBracket[0])
				  {
					// layer ends before encountering contact. Simple pure one-fluid layer
					exhaustedContacts = true;
					if (vintage==0) 
					  {
						fcpInfo.t_contact[i+1] = t.value;
					  }
					else
					  {
						fcInfo.t_contact[i+1] = t_find - 2.0*(z_find -z_top)/vp_eff;					  
					  }
				  }
				else
				  {
					exhaustedContacts = false;
					// Solve for contact time accordingly:
					// zBracket[1] = z_find + 0.5*(t_contact - t_find) * vp_eff
					t_contact = t_find - 2.0*(z_find - zBracket[0])/vp_eff;
					// set new time
					t_find = t_contact;
					// set new depth at encountered contact
					z_find = zBracket[0]; 
				  }
                //_m._SI.EH.println("IC %g %g\n",z_top , zBracket[0]);
			  }

			// set sub-layer properties
			fcpInfo.fluid_type[i] = fluid_type;
			fcpInfo.rho_eff[i] = rho_eff;
			fcpInfo.vp_eff[i] = vp_eff;
			fcpInfo.vs_eff[i] = vs_eff;
			fcpInfo.S_h[i] = S_h;
			fcpInfo.rho_h[i] = rho_h;

			// escape when layer has accounted for all internal fluid contacts
			if (!exhaustedContacts) 
			  {
				fcpInfo.n_internal_sublayers++;
			  }
			else
			  {
				break;
			  }
		  }
        if (fcpInfo.n_internal_sublayers > 4)
		  {
			_m._SI.EH.println("n_internal_sublayers > 4 in fluid-contact mode. Illegal state. FIXME!");
            _m._SI.EH.println("IC %g %g\n",z_top , zBracket[0]);
			return false;
		  }
		// flip arrays into useful top-down order...
		if (!descending) fcpInfo.reverseArrays();
	
		if (fcpInfo.n_internal_sublayers > 1)
		  {
			// compute "upscaled" layer effective properties as time-interval weighted averages
			double t_thick,t_thick_tot=0;                // time thickness
			double z_thick,z_thick_tot=0;                // sublayer thickness
			double z_hydro_thick, z_hydro_thick_tot=0;   // hydrocarbon sublayer thickness
			double z_thick_max=0;
			short shifted_fluid_type;

			rho_eff = vp_eff = vs_eff = S_h = rho_h = 0.0;
			S_oil_eff = S_gas_eff = S_lsg_eff = 0;
			// stuff number of internal sublayers into bits 12,11
			// will be in [0,1,2,3], i.e. needs 2 bits
			shifted_fluid_type = (short)(fcpInfo.n_internal_sublayers-1);
			shifted_fluid_type <<= 10;  
			fluid_type = shifted_fluid_type;

			for (i=0; i< fcpInfo.n_internal_sublayers; i++)
			  {
				t_thick = fcpInfo.t_contact[i+1] - fcpInfo.t_contact[i]; // time thickness
				z_thick = t_thick * fcpInfo.vp_eff[i];                  // sublayer 2way thickness
				z_hydro_thick = z_thick * fcpInfo.S_h[i];               // hydrocarbon column thickness

				t_thick_tot += t_thick;                                // Norms...
				z_thick_tot += z_thick;
				z_hydro_thick_tot += z_hydro_thick;

				rho_eff += z_thick * fcpInfo.rho_eff[i];
				vp_eff += z_thick;                       // this average ensures overall layer 2way thickness is right
				vs_eff += t_thick * fcpInfo.vs_eff[i];
				S_h += z_hydro_thick;
				rho_h += z_hydro_thick * fcpInfo.rho_h[i];

				shifted_fluid_type = fcpInfo.fluid_type[i];         // will be in [0,1,2,3], i.e. needs 2 bits
				shifted_fluid_type <<= (2*(i+1));                  // bitpairs [10,9] [8,7] [6,5] [4,3]
				fluid_type += shifted_fluid_type;                  // accumulate total internal fluid types

				switch (fcpInfo.fluid_type[i])
				  {
				  case FluidStates.FluidTypes.oil:
					S_oil_eff += z_hydro_thick;
					break;
				  case FluidStates.FluidTypes.lsg:
					S_lsg_eff += z_hydro_thick;
					break;
				  case FluidStates.FluidTypes.gas:
					S_gas_eff += z_hydro_thick;
					break;
				  }

				// use majority fluid_type of non-brine fluids (on thickness criterion)
				// downweighted saturation will cope with majority brine cases
				if ((z_thick > z_thick_max) && (fcpInfo.fluid_type[i] != FluidStates.FluidTypes.brine)) 
				  {
					z_thick_max = z_thick;
					shifted_fluid_type = fcpInfo.fluid_type[i];         // will be in [0,1,2,3], i.e. needs 2 bits
					fluid_type &= (0xFFFF ^ 3);       // AND with 1111..00 (from 0xFFFF ExOr 3 ) to clear bits [2,1]
					fluid_type |= shifted_fluid_type; // stuff majority fluid type into bits [2,1]
				  }
			  }
			if (t_thick_tot > 0)
			  {
				vp_eff /= t_thick_tot;
				vs_eff /= t_thick_tot;
			  }
			if (z_thick_tot > 0)
			  {
				rho_eff /= z_thick_tot;
				S_h /= z_thick_tot;
				S_oil_eff /= z_thick_tot;
				S_gas_eff /= z_thick_tot;
				S_lsg_eff /= z_thick_tot;
			  }
			if (z_hydro_thick_tot>0) 
			  {
				rho_h /= z_hydro_thick_tot;
			  }
		  }
	  }
	if (vp_eff<=0.0 || vs_eff<=0 || rho_eff<=0) return false;
	
	// Update suspended depths
	updateSuspendedDepths(descending,suspended_depth_hook,t_below,thickness);

	// Set majority fluid type
	if (fcpInfo!=null) fcpInfo.fluid_type_eff = fluid_type;

	return true;
  }

  final void updateSuspendedDepths(boolean descending, double suspended_depth_hook,double t_below, double thickness)
  {
	// Update suspended depths
	if (vintage==0)
	  {
		if (descending)
		  {
			// going down
			suspended_depth_top = suspended_depth_hook; // passed in
			suspended_depth_bottom = suspended_depth_top 
			  + 0.5*(t_below - t.value) * vp_eff;
		  }
		else
		  {
			// Going up: 
			suspended_depth_bottom = suspended_depth_hook;					
			suspended_depth_top = suspended_depth_bottom 
			  - 0.5*(t_below- t.value) * vp_eff;
		  }	
	  }
	else
	  {
		// just copy thicknesses from vintage 1
		if (descending)
		  {
			// going down
			suspended_depth_top = suspended_depth_hook; // passed in
			suspended_depth_bottom = suspended_depth_top + thickness;
		  }
		else
		  {
			// Going up: 
			suspended_depth_bottom = suspended_depth_hook;		
			suspended_depth_top = suspended_depth_bottom - thickness;
		  }
	  }
  }


  /** 
	  For impermeable rock mixture, in current state, compute Backus upscaled velocities etc into state vector
	  vp_vs_rho[3]=[vp,vs,rho]
  */
  void computeImpermeableProps(double vp_vs_rho[])
  {
	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		if (NG_hard.value==0.0)
		  {
			vp_vs_rho[0] = vp_m.value;
			vp_vs_rho[1] = vs_m.value;
			vp_vs_rho[2] = rho_m.value;
		  }
		else if (NG_hard.value==1.0)
		  {
			vp_vs_rho[0] = vp_hard.value;
			vp_vs_rho[1] = vs_hard.value;
			vp_vs_rho[2] = rho_hard.value;
		  }
		else
		  {
			// shale moduli
			double M_m = rho_m.value * sqr(vp_m.value);
			double mu_m = rho_m.value * sqr(vs_m.value);
			// hard-rock moduli
			double M_hard = rho_hard.value * sqr(vp_hard.value);
			double mu_hard = rho_hard.value * sqr(vs_hard.value);

			// NG_hard weighted average for density
			vp_vs_rho[2] = NG_hard.value * rho_hard.value + (1.0-NG_hard.value)*rho_m.value;
			// And moduli are averaged using Backus
			double M_eff = 1.0/(NG_hard.value/M_hard + (1.0-NG_hard.value)/M_m);
			// Backus average of shear modulus
			double mu_eff = 1.0/(NG_hard.value/mu_hard + (1.0-NG_hard.value)/mu_m);
			  
			// and resulting equivalent velocities
			vp_vs_rho[0] = Math.sqrt(M_eff/vp_vs_rho[2]);
			vp_vs_rho[1] = Math.sqrt(mu_eff/vp_vs_rho[2]);
		  }
	  }
	else
	  {
		vp_vs_rho[0] = vp_m.value;
		vp_vs_rho[1] = vs_m.value;
		vp_vs_rho[2] = rho_m.value;
	  }
  }


  // on entry here, i_v>=1 (i.e. indexes a later vintage)
  void setPropertiesFromPreviousVintage(int i_v,LayerProperties LP,double delta_p)
  {
	int i_SP,code;
	double c,b,dVp;
	RockandFluidData.RegressionSet _RS=null;
	ModelPriorDescriptions.MetaInfo[] MI;
	double[] stress_path_coeffs_m;
	double[] stress_path_coeffs_h;
	

	_RS = LP.info4D.RFD_ImpermeableRegSet;
	stress_path_coeffs_m = (_RS!=null) ? _RS.stress_path_coeffs : null;
	_RS = LP.info4D.RFD_HardRegSet;
	stress_path_coeffs_h = (_RS!=null) ? _RS.stress_path_coeffs : null;

 
	// preserve NG
	NG.setPosValue(LP.NG.value); 

				
	// stress-path index is vintage-index - 1
	i_SP = i_v-1;	
	// map velocities according to stress correlation coefficient c
	// and rock-specific coefficient R
	MI = _m.MPD.MetaInfo;
	c = MI[LP.thisLayerIndx].U_coeff.value;


	// (1) Vp response
	dVp = stress_path_coeffs_m[i_SP] * c * delta_p;
	vp_m.setPosValue(LP.vp_m.value + dVp);

	// (2) Shear response
	// vs_m.setPosValue(LP.vs_m.value);     // default no vs shift	
	// Modest suggestion - use loading trend slope:
	_RS = LP.info4D.RFD_ImpermeableRegSet;
	b = (_RS!=null) ? _RS.Vs_versus_Vp.slope : 0.0;
	vs_m.setPosValue(LP.vs_m.value + b * dVp);
	

	// (3) Density response
	rho_m.setPosValue(LP.rho_m.value);   // Preserve density under stress

	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		// preserve vol.frac hard
		NG_hard.setPosValue(LP.NG_hard.value);
		if (stress_path_coeffs_h != null)
		  {
			// (1) Vp
			dVp = stress_path_coeffs_h[i_SP] * c * delta_p;
			vp_hard.setPosValue(LP.vp_hard.value + dVp);
			// (2) Shear response
			// vs_hard.setPosValue(LP.vs_hard.value);     // default no vs shift	
			// Modest suggestion - use loading trend slope:
			_RS = LP.info4D.RFD_HardRegSet;
			b = (_RS!=null) ? _RS.Vs_versus_Vp.slope : 0.0;
			vs_hard.setPosValue(LP.vs_hard.value + b * dVp);
			// (3) density
			rho_hard.setPosValue(LP.rho_hard.value);  // Preserve density under stress
		  }
	  }

	// permeable
	_RS = LP.RFD_PermeableRegSet;
	// (1) Vp
	if (_RS!=null)
	  {
		dVp = _RS.stress_path_coeffs[i_SP] * c * delta_p;
	  }
	else 
	  {
		dVp=0.0;
	  }
	vp_s.setPosValue(LP.vp_s.value + dVp);

	// (2) Vs 	
	// vs_s.setPosValue(LP.vs_s.value);    // Preserved Vs
	// Modest suggestion - use loading trend slope:   
	b = (_RS!=null) ? _RS.Vs_versus_Vp.slope : 0.0;
	vs_s.setPosValue(LP.vs_s.value + b * dVp);

	// (3) 
	phi_s.set(LP.phi_s);                // 4D: Neglect loss of porosity

	// 4D FIXME check:
	// copy properties that *do not* appear as free parameters in current map
	// if they will do anything
	// Here, we do only oil & gas.
	ModelPriorDescriptions.MetaInfo _MI = MI[thisLayerIndx];
	if (MI[thisLayerIndx].is_permeable)
	  {
		HashMap _map_to = _m.MS.map_to; 

		if (_MI.S_oil!=null) enforceOrdering(S_oil, LP.S_oil,_MI.S_oil,i_v,_map_to);
		if (_MI.S_gas!=null) enforceOrdering(S_gas, LP.S_gas,_MI.S_gas,i_v,_map_to);

		// copy other properties
		//OIL
		if (_map_to.get(rho_oil) == null)
		  {
			rho_oil.value = LP.rho_oil.value;
		  }
		if (_map_to.get(vp_oil) == null)
		  {
			vp_oil.value = LP.vp_oil.value;
		  }
		// GAS
		if (_map_to.get(rho_gas) == null)
		  {
			rho_gas.value = LP.rho_gas.value;
		  }
		if (_map_to.get(vp_gas) == null)
		  {
			vp_gas.value = LP.vp_gas.value;
		  }
		// LSG
		if (_map_to.get(S_lsg) == null)
		  {
			S_lsg.value = LP.S_lsg.value;
		  }
		if (_map_to.get(rho_lsg) == null)
		  {
			rho_lsg.value = LP.rho_lsg.value;
		  }
		if (_map_to.get(vp_lsg) == null)
		  {
			vp_lsg.value = LP.vp_lsg.value;
		  }
		// BRINE
		if (_map_to.get(rho_b) == null)
		  {
			rho_b.value = LP.rho_b.value;
		  }
		if (_map_to.get(vp_b) == null)
		  {
			vp_b.value = LP.vp_b.value;
		  }
	  }
  }

  private final void enforceOrdering(Property P, Property Pprev,
									 ModelPriorDescriptions.PriorDetails[] _PD,
									 int i_v, 
									 HashMap _map_to)
  {
	if (_map_to.get(P) == null) 
	  {
		P.value = Pprev.value;
	  }
	else 
	  {
		if (i_v < _PD.length)
		  {
			ModelPriorDescriptions.PriorDetails PD = _PD[i_v];
			byte code;
			if ((code = PD.previousVintageOrderingCode) != 0)
			  { 
				switch (code)
				  {
				  case 1:
					// same function as monotonicDeltaProp4Dlist
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
	  }
  }





  /** 
	  Get fluid type from contacts info at requested depth
  */
  short getFluidTypeAt(double z_find,double[] zBracket,boolean descending)
  {
	FluidUnitsAndContacts FUC;
	FluidUnitsAndContacts.FluidUnit fu = null;
	boolean bracketLeftFound;
	boolean bracketRightFound;
	int v,j,k,jl=0,kl=0,jr=0,kr=0;
	double zc;
	short fluid_type;

	FUC = _m.fluidUnitsAndContacts;
	bracketLeftFound=false;
	bracketRightFound=false;

	// get left bracket(search down)
	
	
	// similar chunk to fu = getFUofThisLayer(thisLayerIndx);
	for (j=0;j<FUC.fluidUnits.length; j++)
	  {
		fu = FUC.fluidUnits[j];
		if ((thisLayerIndx >= (fu.firstLayerIndx-1)) && (thisLayerIndx <= (fu.finalLayerIndx-1)))
		  {
			jl = jr = j;
			break;
		  }
	  }
	if (fu == null)
	  {
		_m._SI.EH.println("Contact fluid bug. Layer unidentifiable with fluid unit");
		return -1;
	  }
	// Use class-instance vintage (unit offset)
	FluidUnitsAndContacts.ContactInfo[] _fui = fu.contactInfo[vintage];
	for (k=0;k<_fui.length; k++)
	  {
		zc =  _fui[k].depth.value; 
		if ((z_find > zc) || (descending && z_find==zc))
		  {
			// layer depth z_find below this contact
			bracketLeftFound = true;
			kl=k;
		  }
	  }
  
	// get right bracket (search up)

	// fu = getFUofThisLayer(thisLayerIndx); same as above
	for (k = _fui.length-1; k>=0; k--)
	  {
		zc =  _fui[k].depth.value;
		if ((z_find < zc) || (!descending && z_find==zc))
		  {
			// layer depth z_find above this contact
			bracketRightFound = true;
			jr=j; kr=k;
		  }
	  }
	// Pass back the depth of the right bracket
	if (bracketRightFound)
	  {
		zBracket[1] = FUC.fluidUnits[jr].contactInfo[vintage][kr].depth.value;
	  }
	else
	  {
		zBracket[1] = Double.MAX_VALUE;
	  }
	_fui = FUC.fluidUnits[jl].contactInfo[vintage];
	// OK, so at depth z_find, we have 
	if (bracketLeftFound)
	  {
		// Pass back the depth of the left bracket
		zBracket[0] = _fui[kl].depth.value;
		if (bracketRightFound)
		  {
			// FUC.fluidUnits[jl].contactInfo[kl].depth.value <= z_find < FUC.fluidUnits[jr].contactInfo[kr].depth.value
			fluid_type = _fui[kl].LF;
			if (_fui[kl].LF != FUC.fluidUnits[jr].contactInfo[vintage][kr].UF)
			  {
				if (jl != jr)
				  {
					// comparisons across different fluid units. Use the contact straddling fluid
					// of the Fluid Unit which owns the current layer
					fluid_type = _fui[kl].LF;
				  }
				else
				  {
					_m._SI.EH.println("Contact fluid determination conflict. "
									  +"Are fluid units & contacts sensible? "
									  +"Are priors too wide in model? ");
					return -1;
				  }
			  }
		  }
		else
		  {
			// FUC.fluidUnits[jl].contactInfo[kl].depth.value <= z_find
			fluid_type = _fui[kl].LF;
		  }
	  }
	else
	  {
		// left bracket
		zBracket[0] = -Double.MAX_VALUE;
		// fluid must be upper of right-bracket contact
		// z_find < FUC.fluidUnits[jr].contactInfo[kr].depth.value
		fluid_type = FUC.fluidUnits[jr].contactInfo[vintage][kr].UF;
	  }	

	return fluid_type;
  }


  /**
	 For cases of negative float values, remap internal properties
	 to truncated versions, translating along lines of principal axes of
	 correlation matrix. Applied to all layer properties. Synchronisations
	 go top-down layerwise, so these should (I believe) come along for the ride...:-d.
	 See: setValuesFromIndex() for example of where synchronisations appear
	 as mapping lists etc.

	 Notice also state-machine like use of requested_value members. This ensures full
	 recomputation of the four mapped fields {vp_s,vs_s,phi_s, float_fraction} from the
	 last occasion a legal quartet of values was set with a set of setPosValue() calls, typically
	 from e.g. setValuesFromIndex(). This is so derivatives computed using the finite difference
	 (e.g. Prop=vp_s with index i)

	 setValuesFromIndex(i, v+dv) - setValuesFromIndex(i, v))/dv

	 make sense, without having to call setValuesFromIndex() for the other 3 properties
	 that get remapped.
  */ 
  private void computeNegativeFloatRemappedLayerProperties()
  {
	if (float_fraction.requested_value >= 0)
	  {
		// This chunk is necessary in case the previous call fell into the
		// else branch and screwed up the  (vp_s|phi_s|vs_s).value 
		// entries. This will clean things up.
		vp_s.setPosValueWithoutMemory(vp_s.requested_value); 
		phi_s.setPosValueWithoutMemory(phi_s.requested_value);
		vs_s.setPosValueWithoutMemory(vs_s.requested_value);
	  }
	else
	  {
		double Ap,Bphi,Cphi,Bs,Dp,vp,vp_prime;

		// Extract required regression coefficients from sand regression set for layer
	
		// vp regression
		Dp = RFD_PermeableRegSet.Vp_versus_depth.float_fraction_coeff;
		// phi regression. Double snap
		Cphi = RFD_PermeableRegSet.Rho_or_Phi_versus_Vp.float_fraction_coeff;
		Bphi = RFD_PermeableRegSet.Rho_or_Phi_versus_Vp.slope;
		// shear
		Bs = RFD_PermeableRegSet.Vs_versus_Vp.slope;

		// current vp
		vp = vp_s.requested_value;

		// adjustments along principal axes
		vp_prime = vp - Dp * float_fraction.requested_value;
		vp_s.setPosValueWithoutMemory(vp_prime); 
		phi_s.setPosValueWithoutMemory(phi_s.requested_value + Bphi*(vp_prime - vp) - Cphi*float_fraction.requested_value);
		vs_s.setPosValueWithoutMemory(vs_s.requested_value + Bs*(vp_prime - vp));

		// and finally snap back float fraction
		float_fraction.setPosValueWithoutMemory(0.0);
	  }
  }

  char NGstate()
  {
	return ((NG.value == 0)? '0' : ((NG.value==1) ? '1' : 'x'));
  }


  /**
	 Class to hold info about fluid contacts in this layer. Outer class
	 only needed for impermeable rocks.
  */
  public static class InternalFluidContactInfo 
  {	
	// here, a copy of m.MPD.MetaInfo[i].is_permeable, for this layer
	boolean is_permeable = false;
	// For all kinds of layers
	byte n_internal_sublayers=0;
	// Include in t_contacts layer top and base
	// i.e. spans i=0...n_internal_sublayers
	double[] t_contact = new double[5];
 
	void reverseTimes()
	{
	  byte i,j,n_half;
	  double tmp;

	  // do swaps of contact times
	  n_half = (byte)((n_internal_sublayers+1)/2);
	  for (i=0,j=n_internal_sublayers; i<n_half;i++,j--)
		{
		  // do swaps...
		  tmp = t_contact[i]; t_contact[i] = t_contact[j]; t_contact[j] = tmp; 
		}
	}
 }

  /**
	 Implementation for permeable rocks: uses InternalFluidContactInfo members
	 plus some extra info....
  */
  public static class InternalPermeableFluidContactInfo extends InternalFluidContactInfo
  {		  
	// This holds label of thickest non-brine phase in layer after a
	// computeEffectivePropertiesUsingContacts() call
	short fluid_type_eff = FluidStates.FluidTypes.brine;

	// These span i=0...n_internal_sublayers-1
	short[] fluid_type = new short[4];
	double[] rho_eff = new double[4];
	double[] vp_eff = new double[4];
	double[] vs_eff = new double[4];
	double[] S_h = new double[4];
	double[] rho_h = new double[4];

	public InternalPermeableFluidContactInfo()
	{
	  is_permeable = true;
	}

	void reverseArrays()
	{
      byte i,j,n_half;
	  double tmp;
	  short tmp_s;

	  n_half = (byte)(n_internal_sublayers/2);
	  // do swaps of properties
	  j = n_internal_sublayers;
    
	  j--;
	  for (i=0; i<n_half;i++,j--)
		{
		  // do swaps...
          tmp_s = fluid_type[i]; fluid_type[i] = fluid_type[j]; fluid_type[j] = tmp_s; 
		  tmp = rho_eff[i]; rho_eff[i] = rho_eff[j]; rho_eff[j] = tmp; 
		  tmp = vp_eff[i]; vp_eff[i] = vp_eff[j]; vp_eff[j] = tmp; 
		  tmp = vs_eff[i]; vs_eff[i] = vs_eff[j]; vs_eff[j] = tmp;
		  tmp = S_h[i]; S_h[i] = S_h[j]; S_h[j] = tmp; 
		  tmp = rho_h[i]; rho_h[i] = rho_h[j]; rho_h[j] = tmp;  
		}
	  // do swaps of contact times
	  n_half = (byte)((n_internal_sublayers+1)/2);
	  for (i=0,j=n_internal_sublayers; i<n_half;i++,j--)
		{
		  // do swaps...
		  tmp = t_contact[i]; t_contact[i] = t_contact[j]; t_contact[j] = tmp; 
		}
	}
  }

  /**
	 4D info
  */
  public static class Info4D
  {
	RockandFluidData.RegressionSet RFD_ImpermeableRegSet=null;
	RockandFluidData.RegressionSet RFD_HardRegSet=null;
  }

  /**
	 Property class. set values through setPosValue if at all possible. Internal
	 fields are then guaranteed consistent
  */
  public static final class Property
  {
	double value;
	double requested_value;    // stores last argument to setPosValue()
	boolean normalised=false;
	boolean pos=true;
	double mean;
	double sigma;

	public Property(double value, boolean normalised)
	{
	  this.value = value;
	  this.normalised = normalised;
	}

	public Property(double value, boolean normalised ,boolean pos)
	{
	  this.value = value;
	  this.normalised = normalised;
	  this.pos = pos;
	}

	/**
	   simple setter from Parent property: copy transient quantities only
	*/

	final void set(Property P)
	{
	  value = P.value;
	  requested_value = P.requested_value;
	}


	/**
	   simple setter for positive physical quantitites
	*/

	final void setPosValue(double value)
	{
	  requested_value = value;
	  if (pos)
		{
		  if (value < 0.0)
			{
			  this.value=0.0;
			  return;
			}
		  if (normalised && value>1.0)
			{
			  this.value=1.0;
			}
		  else
			{
			  this.value = value;
			}
		}
	  else
		{
		  if (normalised) this.value = max(min(value,1.0),0);
		  else this.value = value;
		}
	}

	/**
	   simple setter for positive physical quantitites, omitting memory
	   of requested_value
	*/

	final void setPosValueWithoutMemory(double value)
	{
	  if (value < 0.0)
		{
		  this.value=0.0;
		  return;
		}
	  if (normalised && value>1.0)
		{
		  this.value=1.0;
		}
	  else
		{
		  this.value = value;
		}
	}
  
	/**
	   simple setter for arbitrary physical quantitites
	*/
	final void setValue(double value)
	{
	  requested_value = value;
	  if (normalised)
		{
		  this.value = max(min(value,1.0),0);
		}
	  else
		{
		  this.value = value;
		}
	}


	/**
	   simple setter for mean and variance, insisting on +ve mean
	*/

	final void setDistr(double mean,double sigma) throws RuntimeException
	{
	  if (sigma < 0.0)
		{
		  throw new RuntimeException("tried to set sigma < 0 for property");
		}
	  if (mean < 0.0)
		{
		  throw new RuntimeException("tried to set sigma < 0 for property");
		}
	  if (sigma == 0)
		{
		  // better set value & requested_value as well..
		  this.value = mean;
		  this.requested_value = mean;
		}
	  this.mean = mean;
	  this.sigma = sigma;
	}

	final void setArbDistr(double mean,double sigma) throws RuntimeException
	{
	  if (sigma < 0.0)
		{
		  throw new RuntimeException("tried to set sigma < 0 for property");
		}
	  if (sigma == 0)
		{
		  // better set value & requested_value as well..
		  this.value = mean;
		  this.requested_value = mean;
		}
	  this.mean = mean;
	  this.sigma = sigma;
	}
  }
}

