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

import au.csiro.JamesGunning.lib.util.MathUtil;

/**
   Basic class for data and fields in Rock Physics for inversion
   
   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/8/2001 
*/

public final class RockandFluidData
{
  // generic curve set with which to form hashmap of rock types
  RockandFluidData.RegressionSet curveSet;

  // generic material_property object with which to form hashmap of fluid and grain types
  RockandFluidData.MaterialProperties material_props;  

  // HashMap of fluid types to look up Fluid props
  //  HashMap<RockandFluidData.MaterialProperties> FluidTypes = null;
  HashMap FluidTypes = null;

  HashMap GrainTypes=null;                 // HashMap of grain types to look up Grain props

  HashMap RockTrends=null;                  // HashMap of rock names to look up RegressionSets

  /**
	 Constructor ???
  */
  public RockandFluidData()
  {
  }

 
  /**
	 Run over all rocks in the HashMap and set up their regression-internal variables.
	 If force==true, all trend curve internal variables (for all three trends) will be remapped.
	 If force==false, they will remap only if they need to (the RockandFluidData.Regression.regressionsVary
	 flag is true)
  */
  void resetTrendData(Model m,boolean force)
  {
	RockandFluidData.RegressionSet reg_set;

	Collection rocks = (Collection)(m.RFD.RockTrends.values());
	// loop over rocks
	for (Iterator i = rocks.iterator(); i.hasNext();) 
	  {
		reg_set = (RockandFluidData.RegressionSet)i.next();
		reg_set.Vp_versus_depth.mapInternalProperties(force);
		reg_set.Vs_versus_Vp.mapInternalProperties(force);
		reg_set.Rho_or_Phi_versus_Vp.mapInternalProperties(force);
	  }
  }

  /**
	 The vs vs vp regression cannot be any old thing. This checks whether
	 unphysical values are likely to be drawn from the vp and vs regressions.
	 If runCompleteCheck==true, we do all rocks.

	 If runCompleteCheck==false, we do only the rocks whose trends vary in space

	 Some properties are specified in the XML filename, but dodgy regressions can result
	 if the trends vary spatially and the parameters read in for the trends are unphysical.
	 Errors thrown by this method are then due to either the XML filename or the
	 read-in parameters from the prior traces.
  */
  //JChen fixed BUGS in this subroutine  
  void checkRockSpeedRegressions(Model m, boolean runCompleteCheck, String filename) throws RuntimeException
  {
     RockandFluidData.RegressionSet reg_set;
     RockandFluidData.Regression vs,vp,phiReg,rhoReg;
     ModelPriorDescriptions.MetaInfo MI;
     ModelPriorDescriptions.PriorDetails loading_depth,sigma_loading_depth;
     HashMap visited_rock_types;
     String rock_name=null;
     double oldDepth=0.0,oldLFIV=0.0, vp_bar, vs_bar, chisq, Vp;
     double MAX_NORMAL_DEVIATE=2.0;
     int j;
     boolean isRockPermeable,checksVaryWithLayer;

     // keep a list of all the rocks we have seen so far...
     // loop over all rock types we seem to need from the meta-Layer data
     visited_rock_types = new HashMap();
     for (int i=0;i<m.MPD.MetaInfo.length;i++) {
        MI = m.MPD.MetaInfo[i];

        // Use correct depth in regressions
        if (MI.depth_rock_curves != MI.depth) {
           loading_depth = MI.depth_rock_curves;
           sigma_loading_depth = MI.sigma_depth_rock_curves;
        }
        else {
           loading_depth = MI.depth;
           sigma_loading_depth = MI.sigma_depth;
        }

        checksVaryWithLayer = (i==0 || ((loading_depth.value != oldDepth ) || (MI.LFIV.value != oldLFIV )));

        /* OLD
        for (isRockPermeable=true, j=0 ; j<=1 ; isRockPermeable=false,j++)
        {
           // do sand-end-members then mud members
           rock_name = (isRockPermeable) ? MI.sand_rock_member : MI.mud_rock_member;
        */

        ROCK_TYPES:
        for (j=0 ; j<=6 ; j++) {
           // do sand-end-members, then mud, then hard_rock members
           switch (j) {
           case 0:
              rock_name = MI.sand_rock_member ; break;
           case 1:
              rock_name = MI.mud_rock_member ; break;
           case 2: 
              rock_name = MI.hard_shale_rock_member ; break;
           case 3:
           case 4:
           case 5:
           case 6:
              if (MI.fluid_dep_hard_shale_rock_members == null) {
                 break ROCK_TYPES;
              }
              else {
                 rock_name = MI.fluid_dep_hard_shale_rock_members[j-3]; // maps to indexes 0,1,2,3
              }
           default:
        }
        isRockPermeable = (j==0);

        //All rocks with non--null names need rock vs vp trends these curves must be sane
        if (rock_name != null) {
           // get its regression set
           reg_set = (RockandFluidData.RegressionSet) m.RFD.RockTrends.get(rock_name);
           vp = reg_set.Vp_versus_depth;
           vs = reg_set.Vs_versus_Vp;

           // Early check on trivial cases where regressions amount to simple constants
           if (vp.error==0) {
              if (!((sigma_loading_depth.value == 0) && (MI.sigma_LFIV.value == 0))) {
                 throw new RuntimeException("Rock " + rock_name + " is only allowed to have " +
                    "sigma_vp=0 if sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
              }
           }
           if ((vs.error==0) && !((vp.error==0) && (sigma_loading_depth.value == 0) 
           && (MI.sigma_LFIV.value == 0))) {
              throw new RuntimeException("Rock " + rock_name + " is only allowed to have " +
                 "sigma_vs=0 if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
           }
           if (reg_set.Rho_or_Phi_versus_Vp.error==0 && !((vp.error==0) 
           && (sigma_loading_depth.value == 0) && (MI.sigma_LFIV.value == 0))) {
              throw new RuntimeException("Rock " + rock_name + " is only allowed to have " +
	      "sigma_phi=0 or sigma_rho_m=0 if sigma_vp=0, sigma_d(maybe rock_curves)=0 and sigma_LFIV=0 also.");
           }

           if (runCompleteCheck || (vp.regressionsVary || vs.regressionsVary)) { 
              // mean p-wave velocity at this layer (BUGS found by JChen)
              //vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value;
              vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value
                     + vp.float_fraction_coeff*MI.float_fraction.value;
              // least likely case of velocity ratio vp/vs
              vs_bar = vs.intercept + vs.slope * vp_bar;
              // Check the likelihood of violating vp > sqrt(2) * vs, with vs fixed at the prior mean, and
              // vp varying. To do this, evaluate chisq for vs fixed at vs_bar, with vp varying, at 
              // physical sanity limit
              // vp=sqrt(2) * vs. If this is highly priorly improbable, chisq will be larger than say, 3.
              // Mode hunting usually fixes Vs, and requires rock physics to be always sane when varying Vp, 
              // with Vs fixed at the mean vs_bar.
              // If the prior keeps strongly away from this unphysical region, we are safe. Test this point.
              Vp = 1.414 * vs_bar;
              chisq = MathUtil.sqr((Vp -vp_bar)/vp.error)
                    + MathUtil.sqr((vs_bar - (vs.intercept + vs.slope * Vp))/vs.error);
              if (checksVaryWithLayer && (m._SI.verbosity > 1)) {
                 m._SI.EH.println("      Rock '%s', layer "+(i+1)
                 +", loading depth = %6f, <v_p> = %6f, <v_s> = %6f, chisq vp/vs sanity test: chisq = %6f",
                 rock_name, loading_depth.value, vp_bar, vs_bar, chisq);
              }
              if (chisq < MAX_NORMAL_DEVIATE) {
                 m._SI.EH.print("%nRock physics errors for rock "+rock_name+"...%n");
                 m._SI.EH.print("%nRegression allows vp/vs ratio to get too low.");
                 m._SI.EH.print("%nAt vs_bar = " + vs_bar+" (computed from vp_bar = " + vp_bar
                 + "), chisq = "+chisq+ " at physical limit point vp=sqrt(2)*vs_bar");
                 throw new RuntimeException("Rock " + rock_name + " has unphysical vs, vp regressions "
                 +"at current trace. See relevant entries in "+filename);
              }
           }
           // (1) check the the porosity regression check
           if (rock_name == MI.sand_rock_member) {
              phiReg = reg_set.Rho_or_Phi_versus_Vp;
              if (runCompleteCheck || phiReg.regressionsVary) {
                 // need vp_bar... (BUGS - fixed by JChen)
                 //vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value;
                 vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value
                        + vp.float_fraction_coeff*MI.float_fraction.value;

                 double phi_bar = phiReg.intercept + phiReg.slope * vp_bar;
                 if (checksVaryWithLayer && (m._SI.verbosity > 1)) {
                    m._SI.EH.println("      Rock '"+rock_name+"', <phi>=%6f",phi_bar);
                 }
                 if (phi_bar <= 0.0) {
                    m._SI.EH.println("<Vp> = "+vp_bar+", mean porosity = "+phi_bar+" for rock "+rock_name);
                    throw new RuntimeException("Rock '" + rock_name + "' has unphysical phi vs vp regressions "
                    +"at current trace. See relevant entries in "+filename);
                 }
                 if (phi_bar > 0.999) {
                    m._SI.EH.println("<Vp> = "+vp_bar+", mean porosity = "+phi_bar+" for rock "+rock_name);
                    throw new RuntimeException("Rock '" + rock_name + "' has unphysical phi vs vp regressions "
                    +"at current trace. See relevant entries in "+filename);
                 }
              }
           }
           // (2) convert GGG power law regressions to linear forms
           // Then record that we have converted this rock
           if (visited_rock_types.get(rock_name)==null) {
              if (rock_name == MI.mud_rock_member) {
                 rhoReg = reg_set.Rho_or_Phi_versus_Vp;
                 if (RockandFluidDataParser.ASSUME_DENSITY_REGRESSIONS_ARE_POWER_LAWS) {
                    // convert GGG regressions to linear forms on mudrocks the first time the mud-rock is 
                    // examined, or if the trends change here rho(vp) = rho_0 vp^exponent is the GGG power
                    // law form

                    // GGG regression
                    if (runCompleteCheck || rhoReg.intercept_PD.varies || rhoReg.slope_PD.varies) {
                       // need vp_bar.... (BUGS - fixed by JChen)
                       //vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value;
                       vp_bar = vp.intercept+ vp.depth_coeff* loading_depth.value + vp.LFIV_coeff* MI.LFIV.value
                              + vp.float_fraction_coeff*MI.float_fraction.value;
                       double exponent,rho_0,rho_bar;

                       // we must have read in a fresh value, or it contains the original unconverted values,
                       // so...
                       rho_0 = rhoReg.intercept_PD.value; // assume first entry is the rho_0
                       exponent = rhoReg.slope_PD.value;  // and second is the power exponent
                       // over write with nearest linear approximation centred on mean p-wave velocity
                       rhoReg.intercept = rho_0 * Math.pow(vp_bar,exponent)*(1.0-exponent);
                       rhoReg.slope = exponent * rho_0 * Math.pow(vp_bar,exponent-1.0);
                       rho_bar = rho_0 * Math.pow(vp_bar,exponent);
                       if (checksVaryWithLayer && (m._SI.verbosity > 1)) {
                          m._SI.EH.println("      Rock '"+ rock_name+
                          "', GGG linearised factors, (using <Vp> = %6f), "
                          +"intercept = %6f, slope = %6g, yields <rho> = %6g",
                          vp_bar, rhoReg.intercept, rhoReg.slope, rho_bar);
                       }
                    }
                 }
              }
              visited_rock_types.put(rock_name,rock_name);
              }
           }
        }
        oldDepth = loading_depth.value;
        oldLFIV = MI.LFIV.value;
     }
     visited_rock_types = null;
  }

  /**
	 Simple trend plotter. Type must be {vp,vs,rho,phi}
  */
  boolean plotTrend(Model m, String type,String rock_name)
  {
	int i,j,n=1;
	double x,xD,v,y;
	
	ModelPriorDescriptions.MetaInfo MI,MIJ;
	ModelPriorDescriptions.PriorDetails loading_depth;
	RockandFluidData.Regression vs,vp,phiReg,rhoReg;

	if (RockTrends.get(rock_name) == null)
	  {
		m._SI.EH.println("Rock " + rock_name+" is not in XML!");
		return false;
	  }
	else
	  {
		RockandFluidData.RegressionSet CS;
		CS = (RockandFluidData.RegressionSet)RockTrends.get(rock_name);
		RockandFluidData.Regression regVp,reg;
		n = m.n_layers;
		vp = CS.Vp_versus_depth;
		RockandFluidData.MaterialProperties grain_props = null;
		if (GrainTypes.get(rock_name) != null) 
		  {
			grain_props = (RockandFluidData.MaterialProperties)GrainTypes.get(rock_name);
		  }

		if (type.equals("vp"))
		  {
			for (i=0;i<n;i++)
			  {
				// j=-1;last layer,   j=0 this layer,  j=1;extend bottom layer
				for (j=((i==0)?0:-1) ; j<=((i==(n-1))?1:0) ; j++)
				  {
					// vertical lines bridging transitions
					MIJ = m.MPD.MetaInfo[Math.min(i+j,n-1)];
					MI = m.MPD.MetaInfo[i];
					xD = MI.depth.value + ((j==1 && MI.dz!=null) ? MI.dz.value : 0);
					y = vp.intercept + vp.depth_coeff* MIJ.depth_rock_curves.value
					  + vp.LFIV_coeff* MIJ.LFIV.value
					  + vp.float_fraction_coeff * MIJ.float_fraction.value;				
					m._SI.EH.println("%g %g",xD, y);
				  }
			  }
		  }
		else if (type.equals("vs"))
		  {
			vs = CS.Vs_versus_Vp;
			for (i=0;i<n;i++)
			  {
				for (j=((i==0)?0:-1) ; j<=((i==(n-1))?1:0) ; j++)
				  {
					MIJ = m.MPD.MetaInfo[Math.min(i+j,n-1)];
					MI = m.MPD.MetaInfo[i];			
					xD = MI.depth.value + ((j==1 && MI.dz!=null) ? MI.dz.value : 0);
					v = vp.intercept+ vp.depth_coeff* MIJ.depth_rock_curves.value
					  + vp.LFIV_coeff* MIJ.LFIV.value	
					  + vp.float_fraction_coeff * MIJ.float_fraction.value;				
					y = vs.intercept+ vs.slope * v;
					m._SI.EH.println("%g %g",xD, y);
				  }
			  }
		  }
		else if (type.equals("rho"))
		  {
			rhoReg = CS.Rho_or_Phi_versus_Vp;
			for (i=0;i<n;i++)
			  {
				for (j=((i==0)?0:-1) ; j<=((i==(n-1))?1:0) ; j++)
				  {
					MIJ = m.MPD.MetaInfo[Math.min(i+j,n-1)];
					MI = m.MPD.MetaInfo[i];			
					xD = MI.depth.value + ((j==1 && MI.dz!=null) ? MI.dz.value : 0);
					v = vp.intercept+ vp.depth_coeff* MIJ.depth_rock_curves.value
					  + vp.LFIV_coeff* MIJ.LFIV.value
					  + vp.float_fraction_coeff * MIJ.float_fraction.value;				
					y = rhoReg.intercept+ rhoReg.slope * v;  // valid shale or sand
					if (grain_props != null)
					  {
						// Reservoir: y = poro. Add FG component
						y += rhoReg.float_fraction_coeff * MIJ.float_fraction.value;
						// y = true regression (fluid) poro now.
						// m._SI.EH.println("%g %g",m._SI.p.add(xD).add(y));
						y = (1-y)*grain_props.rho // rock matrix + floats
						  +  y*CS.rho_b_ref; // ref fluid
					  }

					m._SI.EH.println("%g %g",xD, y);
				  }
			  }
		  }
		else if (type.equals("phi"))
		  {
			if (grain_props == null)
			  {
				m._SI.EH.println("Illegal phi curve type for non reservoir rock: "+rock_name);
				return false;
			  }
			phiReg = CS.Rho_or_Phi_versus_Vp;
			for (i=0;i<n;i++)
			  {
				for (j=((i==0)?0:-1) ; j<=((i==(n-1))?1:0) ; j++)
				  {
					MIJ = m.MPD.MetaInfo[Math.min(i+j,n-1)];
					MI = m.MPD.MetaInfo[i];			
					xD = MI.depth.value + ((j==1 && MI.dz!=null) ? MI.dz.value : 0);
					v = vp.intercept+ vp.depth_coeff* MIJ.depth_rock_curves.value
					  + vp.LFIV_coeff* MIJ.LFIV.value	
					  + vp.float_fraction_coeff * MIJ.float_fraction.value;				
					y = phiReg.intercept+ phiReg.slope * v 
					  + phiReg.float_fraction_coeff * MIJ.float_fraction.value;
					m._SI.EH.println("%g %g",xD, y);
				  }
			  }
		  }
		else if (type.equals("imp"))
		  {
			rhoReg = CS.Rho_or_Phi_versus_Vp;
			for (i=0;i<n;i++)
			  {
				for (j=((i==0)?0:-1) ; j<=((i==(n-1))?1:0) ; j++)
				  {
					MIJ = m.MPD.MetaInfo[Math.min(i+j,n-1)];
					MI = m.MPD.MetaInfo[i];				
					xD = MI.depth.value + ((j==1 && MI.dz!=null) ? MI.dz.value : 0);
					v = vp.intercept+ vp.depth_coeff* MIJ.depth_rock_curves.value
					  + vp.LFIV_coeff* MIJ.LFIV.value	
					  + vp.float_fraction_coeff * MIJ.float_fraction.value;				
					// y = density (shale) or poro (sand)
					y = rhoReg.intercept+ rhoReg.slope * v;
					if (grain_props != null)
					  {
						// y = porosity => convert to density
						// Add on float component to porosity
						y += rhoReg.float_fraction_coeff * MIJ.float_fraction.value;
						y = (1-y)*grain_props.rho // rock matrix
						  + MIJ.float_fraction.value * grain_props.rho // float grains
						  + y*CS.rho_b_ref;       // ref fluid
					  }				
					// impedance = Vp * density
					m._SI.EH.println("%g %g",xD, y*v);
				  }
			  }
		  }
		else
		  {
			m._SI.EH.println("Illegal curve type "+type+": must be one of {vp,vs,rho,phi}");
			return false;
		  }
	  }
	return true;
  }


  /** 
	  class to represent a set of rock curves for each facies
  */
  static class RegressionSet
  {
	String ref_brine_name = null;
	/**
	   reference fluid density
	*/
	double rho_b_ref;
	/**
	   reference fluid vp
	*/
	double vp_b_ref;
	/**
	   reference fluid K_b
	*/
	double K_b_ref;
	/**
	   Grain bulk K modulus
	*/
	double K_g;
	
	/**
	   Grain density
	*/
	double rho_g;

	RockandFluidData.Regression Vp_versus_depth;      // linear regression  of vp vs t to be read from xml file
	RockandFluidData.Regression Vs_versus_Vp;     // linear regression  of vs vs vp to be read from xml file
	RockandFluidData.Regression Rho_or_Phi_versus_Vp;    // linear regression  of density vs vp to be read from xml file

	/* 
	   array of Nvintages-1 stress path coeffs
	*/
	double[] stress_path_coeffs = null;

	public RegressionSet()
	{
	  Vp_versus_depth = new RockandFluidData.Regression();
	  Vs_versus_Vp = new RockandFluidData.Regression();
	  Rho_or_Phi_versus_Vp = new RockandFluidData.Regression();
	}

	void setBrineConstants(RockandFluidData.MaterialProperties material_props)
	{
	  rho_b_ref = material_props.rho;
	  vp_b_ref = material_props.vp;
	  K_b_ref = rho_b_ref * vp_b_ref * vp_b_ref;	
	}
  }


  /**
	 general set of coefficients for trend regressions	 
  */
  static class Regression
  {
	/*
	  used only for terms with phi_float:
	  Vp = A + B d +C LFIV +D float_fraction
	  phi = A+B vp + C float_fraction
	*/
	double depth_coeff=0.0;
	double LFIV_coeff=0.0;
	double float_fraction_coeff=0.0; 

	/* generics */
	double intercept;
	double slope;
	double error;

	boolean regressionsVary = false;    // switched on if regressions vary in space

	// And if any of the relevent coefficients vary in space, we'll need these..
	ModelPriorDescriptions.PriorDetails intercept_PD = new ModelPriorDescriptions.PriorDetails(0.0); 
	ModelPriorDescriptions.PriorDetails slope_PD = new ModelPriorDescriptions.PriorDetails(0.0); 
	ModelPriorDescriptions.PriorDetails error_PD = new ModelPriorDescriptions.PriorDetails(1.0); 

	// used only for Vp = A + B d +C LFIV + D float_fraction regression
	ModelPriorDescriptions.PriorDetails depth_coeff_PD = new ModelPriorDescriptions.PriorDetails(0.0);  
	// used only for Vp = A + B d +C LFIV + D float_fraction regression
	ModelPriorDescriptions.PriorDetails LFIV_coeff_PD = new ModelPriorDescriptions.PriorDetails(0.0);   
	// used only for Vp = A + B d +C LFIV + D float_fraction regression
	ModelPriorDescriptions.PriorDetails float_fraction_coeff_PD = new ModelPriorDescriptions.PriorDetails(0.0); 

	/**
	   This maps the PriorDetails values to the internal fields
	*/
	void mapInternalProperties(boolean forceRemapping)
	{
	  // perform remapping only if regressions vary or we force them
	  if (forceRemapping || regressionsVary)
		{
		  intercept = intercept_PD.value;
		  depth_coeff = depth_coeff_PD.value;
		  LFIV_coeff = LFIV_coeff_PD.value;
		  float_fraction_coeff = float_fraction_coeff_PD.value;
		  slope = slope_PD.value;
		  error = error_PD.value;
		}
	}
  }
  /**
	 Class to represent a fluid or rock's useful sonic properties
  */
  static class MaterialProperties
  {
	/** 
		Name of fluid or rock
	*/
	String name;

	/**
	   P-wave velocity at in-situ conditions
	*/
	double vp=0; 
	boolean vp_varies=false;
	/**
	   P-wave velocity std dev at in-situ conditions
	*/
	double sigma_vp=0;    
	boolean sigma_vp_varies=false;
	/**
	   S-wave velocity at in-situ conditions: likely not used for fluids
	*/
	double vs=0;     // likely not be used for fluids
	/**
	   S-wave velocity std dev at in-situ conditions: likely not used for fluids
	*/
	double sigma_vs=0;     // likely not be used for fluids
	/**
	   Density at in-situ conditions
	*/
	double rho=0;
	boolean rho_varies=false;
	/**
	   Density std dev at in-situ conditions
	*/
	double sigma_rho=0;
	boolean sigma_rho_varies=false;
  }
}




