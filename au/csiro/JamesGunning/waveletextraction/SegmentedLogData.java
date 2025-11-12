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
import static java.lang.Math.*;


import cern.colt.list.*;
import cern.colt.matrix.*;
import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.stat.segmentation.*;
import au.csiro.JamesGunning.lib.util.*;

import static au.csiro.JamesGunning.lib.util.MathUtil.sqr;

/**
   Data class for storing upscaled 
   log information corresponding to segmentation. 
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class SegmentedLogData
{ 
  private ErrorHandler EH;
  private int verbosity;
  /**
	 Useful ref to parent wellinfo
  */
  GeneralWellInfo _WI = null;
  /**
	 MD: measured depth of top of data-block with
	 properties as below
  */
  double[] Zm;
  /**
	 True depth of top of data-block with
	 properties as below
  */
  double[] Z;
  /**
	 Segmented layer top times: Will contain times of layer tops:
	 filled in by the time to depth conversion routine anew for every
	 realisation of the t2d mapping. The indexing has exactly the same
	 meaning as for Z[].
  */
  double[] T;
  /**
	 Backus-averaged Vp (p-wave velocity) data.
  */
  double[] Vp;
  /**
	 Harmonic averaged Vp (p-wave velocity) data.
  */
  double[] VpHarm;
  /**
	 Vs (s-wave velocity) data
  */
  double[] Vs;
  /**
	 -Possible- gamma data
  */
  double[] gamma;
  /**
	 Density data
  */
  double[] rho;
  /**
	 -Possible- petrophysics indicator data
  */
  byte[] facies_indicator;
  /**
	 Segmented zero--offset reflection coefficient
  */
  double[] R_near;
  /*
	Segmented O(sin^2(theta)) reflection coefficient
  */
  double[] R1_AVO;
  double[] Vsq_AVO;
  /*
	Useful arrays for PS data
  */
  double[] Rps_AVO;
  double[] Rps_AVO2;

  /**
	 Fields used if Shuey equation for reflectivity used
  */
  double[] R2_Shuey_AVO;
  /**
	 Fields used if isotropic Zoeppritz  equation for reflectivity used
  */
  double[] vp1,vp2,vs1,vs2,rho1,rho2;

  /**
	 Anisotropy fields
  */
  double[] DeltaR_anis=null;

  /**
	 Segmented delta (v_s*rho) term for output traces only
	 (not used in computation)
  */
  double[] deltaZs;
  /**
	 Segmented delta (rho) term for output traces only
	 (not used in computation)
  */
  double[] deltaRho;
  /**
	 reference to log Suite of interest
  */
  LogData lData;
  /**
	 segmentation points that will be generated
  */
  IntArrayList segIndexValues;

  /**
	 Units conversion factor when comparing segmented log velocities
	 to stack velocities.
	 Usually non-unity if there is a ft<->m issue
  */
  double lengthUnitsConversionFactor=1.0;

  /**
	 Units conversion factor when comparing segmented log velocities 
	 to time-to-depth checkshot velocities.
	 Usually non-unity if there is a ft<->m issue
  */
  double lengthUnitsConversionFactorCheckshot=1.0;

  /**
	 Density conversion from logs to SI
  */
  double densityUnitsConversionFactor=1.0;

  SaturationModelInfo SMI=null;

 
  /**
	 Constructor.
  */
  public SegmentedLogData(GeneralWellInfo WI, ErrorHandler EH, int verbosity)
  {
	this._WI = WI;
	this.EH = EH;
	this.verbosity = verbosity;
  }
  /**
	 Segment this suite of data and upscale properties using suitable
	 averages. Use fancy rigorous ML blocking if desired.
	 The Vdenom is only for debugging purposes, and can be set to 0 without harm
  */
  public void segmentLogSuite(WaveletExtractor _WE, LogData logData,
							  boolean useMLblocking, boolean dumpPandSreflectivity, 
							  double Vdenom)
  {
	int i,j,n,N,cnt,indx;
	double[] data;
	double vp,vs;
	double MpEff,rhoEff,muEff,vpHarmEff;
	double gammaEff=0;
	double rhoBar,vpBar,vsBar;
	double delta_rho,delta_vp,delta_vs,ratio,beta,tmp;
	double Dr,Dp,Ds;
	boolean badData=false;
	double DT_scaling,DTS_scaling;

	// reference for other methods
	lData = logData;

	// Ok first of all get a decent effective p-impedance data suite Imp = rho*vp = rho/DT
	// This is what the segmentation will be based on

	N = logData.DT.length;
	data = new double[N];
	for (i=0; i<N; i++)
	  {
		data[i] = logData.rho[i]/logData.DT[i];
	  }
	// Now segment this data
	DataSegmentor segmentor = new DataSegmentor(data);

	if (useMLblocking)
	  {
		IntArrayList[] segIndexValuesArray;
		segIndexValuesArray = segmentor.dynamicProgrammingMLSegmentWithDecimationIfRequired(false);
		segIndexValues = segIndexValuesArray[0];
		if (segIndexValues==null) 
		  {
			// data size probably too big: use simpler approx method
			EH.println("Data too large (n="+N+") for ML blocking: using approximate method...");
			segIndexValues = segmentor.aggregativeSegment(true,1.0);
		  }		
	  }
	else
	  {
		segIndexValues = segmentor.aggregativeSegment(false,_WE.aggregativeBlockingCutoff);
	  }

	// Now we do some classical Backus style upscaling to get the new blocked logs
	// These all assume the isotropic model of each layer
	// n = size of upscaled data
	n = segIndexValues.size();
	if (logData.haveGamma) gamma = new double[n-1];
	rho = new double[n-1];
	Vp = new double[n-1];
	VpHarm = new double[n-1]; // for slowness average
	Vs = new double[n-1];
	Zm = new double[n];
	T = new double[n];
	R_near = new double[n];
	R1_AVO = new double[n];
	Vsq_AVO = new double[n];

	if (_WE.doShuey) 
	  {
		R2_Shuey_AVO = new double[n];
	  }
	if (_WI.needPSreflectivities)
	  {
		Rps_AVO = new double[n];
		Rps_AVO2 = new double[n];
	  }
	LinearInterpolation MDToTVDInterp = null;

	if (_WE.doAnisotropicModel) 
	  {
		facies_indicator = new byte[n];
		DeltaR_anis = new double[n]; // FIXME - ANISO
		// MD to TVD mapping from survey. MUST BE MONOTONIC
		MDToTVDInterp = 
		  new LinearInterpolation(logData.survey.Zm, logData.survey.Z, logData.survey.Z.length, _WE.EH);
	  }

	if (dumpPandSreflectivity)
	  {
		deltaZs = new double[n];
		deltaRho = new double[n];
	  }

	DT_scaling = logData.getSonicLogScaling();

	for (i=0;i<n-1;i++)
	  {
		indx = segIndexValues.getQuick(i);
		// measured depths of layer tops. Set these to midpoints of changepoints intervals
		Zm[i] = (logData.Zm[indx]+logData.Zm[indx+1])/2.0 ;

		// now upscale properties......

		rhoEff = MpEff = muEff = vpHarmEff = 0.0;
		if (logData.haveGamma) gammaEff=0.0;

		// Form averages within each segment
		//System.out.println(i+" "+indx);
		for (cnt=0,j = indx+1; j <= segIndexValues.getQuick(i+1) && i<n; j++,cnt++)
		  {
			if (logData.rho[j] <= FundamentalConstants.BAD_NUMBER
				|| logData.DT[j] <= FundamentalConstants.BAD_NUMBER
				|| logData.DTshear[j] <= FundamentalConstants.BAD_NUMBER)
			  {
				badData = true;
				break;
			  }
			// we could be fussy here about block weights, but logs are 99.9% equispaced...
			// gamma. Cosmetic, so be less fussy than DT etc
			if (logData.haveGamma && (logData.gamma[j] > FundamentalConstants.BAD_NUMBER))
			  {
				gammaEff += logData.gamma[j];
			  }

			// density
			rhoEff += logData.rho[j];

			// velocities
			vp = DT_scaling/logData.DT[j];
			vs = DT_scaling/logData.DTshear[j];
			
			// p-modulus
			MpEff += 1.0/(logData.rho[j]*vp*vp);

			// s-modulus
			muEff += 1.0/(logData.rho[j]*vs*vs);

			// slowness
			vpHarmEff += 1.0/vp;
		  }
		if (badData)
		  {
			rho[i] = Vp[i] = Vs[i] = VpHarm[i] = FundamentalConstants.BAD_NUMBER;
			if (logData.haveGamma) gamma[i] = FundamentalConstants.BAD_NUMBER;
		  }
		else
		  {
			// arithmetic av. of gamma
			if (logData.haveGamma) gamma[i] = gammaEff/(double)cnt;

			// density is arithmetic average
			rhoEff /= (double)cnt;
			rho[i] = rhoEff;

			// Slowness is Harm average
			VpHarm[i] = (double)cnt/vpHarmEff;

			// Backus average of p-wave modulus.
			MpEff = (double)cnt/MpEff;
			Vp[i] = Math.sqrt(MpEff/rhoEff); 

			// And the same for s-waves...
			muEff = (double)cnt/muEff;
			Vs[i] = Math.sqrt(muEff/rhoEff);

			// reflection coefficients
			if (i>0)
			  {
				rhoBar = 0.5*(rho[i]+rho[i-1]);
				vpBar = 0.5*(Vp[i]+Vp[i-1]);
				vsBar = 0.5*(Vs[i]+Vs[i-1]);
				delta_rho = rho[i]-rho[i-1];
				delta_vp = Vp[i]-Vp[i-1];
				delta_vs = Vs[i]-Vs[i-1];
				Dr = delta_rho/rhoBar;
				Dp = delta_vp/vpBar;
				Ds = delta_vs/vsBar;

				R_near[i] = 0.5*(Dr + Dp);

				// Linearized reflectivity requires:
				tmp = Dr + 2.0*Ds;
				beta = vsBar/vpBar;
				// ratio = delta_vp/(2.0*vpBar) - (2.0*vsBar*vsBar*tmp/(vpBar*vpBar));
				ratio = Dp/2.0 - (2.0*sqr(beta)*tmp);
				// NB, in Mavko, angle is average of theta_1 & theta_2, in Castagna, it is theta_1. Sigh.

				R1_AVO[i] = ratio;
				Vsq_AVO[i] = sqr(Vp[i-1]);   // will be dimension-corrected to match stack params later!
				if (_WE.doShuey)
				  {
					R2_Shuey_AVO[i] = Dp/2.0; // delta_vp/(2.0*vpBar);					 
				  }
				if (_WI.needPSreflectivities)
				  {
					// Dimensionless linearized expression
					Rps_AVO[i] = -0.5*(Dr + 2.0*beta*tmp);  
					Rps_AVO2[i] = ((1.0 + beta*(8 + 9.0*beta))*Dr  +  beta*(16.0 + 24.0*beta)*Ds)/12.0;  
				  }
				if (dumpPandSreflectivity)
				  {
					deltaRho[i] = 0.5*Dr;
					deltaZs[i] = 0.5*(Dr + Ds);
				  }
				if (_WE.doAnisotropicModel) 
				  {
					// compute this once and for all - later models may ignore it if m.hasAnisotropy == false
					double zTVD, deltaThomsen1, deltaThomsen2, gammaShaleThreshold, densityShaleThreshold;
					zTVD = MDToTVDInterp.interpolate(Zm[i], i!=0);

					boolean isShale_i_minus_1, isShale_i;
					
					if (logData.haveGamma)
					  {
						isShale_i = _WE.m.petro.isShale(logData.haveGamma,zTVD,gamma[i],rho[i]);
						isShale_i_minus_1 = _WE.m.petro.isShale(logData.haveGamma,zTVD,gamma[i-1],rho[i-1]);
					  }
					else
					  {
						// gamma won't be used...
						isShale_i = _WE.m.petro.isShale(logData.haveGamma,zTVD,0,rho[i]);
						isShale_i_minus_1 = _WE.m.petro.isShale(logData.haveGamma,zTVD,0,rho[i-1]);
					  }
					// version 1: 0,1, indicator
					if (_WE.m.petro.version==1) facies_indicator[i] = isShale_i ? (byte)1 : 0;  // 0,1, indicator

					deltaThomsen1  = isShale_i_minus_1 ? 1 : 0;
					deltaThomsen2  = isShale_i ? 1 : 0;

					/*
					if (logData.haveGamma && _WE.m.petro.GRThreshProvided)
					  {
						gammaShaleThreshold = _WE.m.petro.GR_intercept + _WE.m.petro.GR_slope * zTVD;
						if (_WE.m.petro.densityThreshProvided)
						  {
							// joint criterion
							densityShaleThreshold = _WE.m.petro.density_intercept + _WE.m.petro.density_slope * zTVD;
							deltaThomsen1  = ((gamma[i-1] > gammaShaleThreshold) && (rho[i-1] > densityShaleThreshold)) ? 1 : 0;
							deltaThomsen2  = ((gamma[i] > gammaShaleThreshold) && (rho[i] > densityShaleThreshold)) ? 1 : 0;
						  }
						else
						  {
							deltaThomsen1  = ((gamma[i-1] > gammaShaleThreshold)) ? 1 : 0;
							deltaThomsen2  = ((gamma[i] > gammaShaleThreshold)) ? 1 : 0;
						  }
					  }
					else
					  {
						densityShaleThreshold = _WE.m.petro.density_intercept + _WE.m.petro.density_slope * zTVD;
						deltaThomsen1  = (rho[i-1] > densityShaleThreshold) ? 1 : 0;
						deltaThomsen2  = (rho[i] > densityShaleThreshold) ? 1 : 0;
					  }
					*/
					// p 68 Mavko, Rock physics HB, to O(theta^2)
					DeltaR_anis[i] = vpBar * vpBar * 0.5*(deltaThomsen2 - deltaThomsen1); 
				  }
				// DEBUG Angle indications *if in same units*. But units conversion not available at this point
				// so don't trust this. Uncomment to see
				/*				
				{ 
				  double theta = Math.toDegrees(vpBar/Math.sqrt(Vdenom));
				  System.out.println("ANGLE "+theta);
				}
				*/
				
			  }
		  }
	  }
	// last layer bottom...
	Zm[n-1] = logData.Zm[segIndexValues.getQuick(n-1)];

	// debugPrintSegmentationResults(); // for parsing only
  }
  
  /**
	 Recompute log properties. If after==false, at saturation 
	 _WI.saturationInfo.Sb, else
	 at model saturation from _WI.saturationInfo.Sa or 
	 current model variable usually in _WI.saturationInfo.satPropMapper
  */
  public void recomputeSegProperties(Model m,boolean dumpPandSreflectivity, boolean doShuey,boolean after)
  {
	int N;
	if ( SMI==null || SMI.sensLookups==null) return;

	N = SMI.sensLookups.size();
	if (N==0) return;

	int i,j,i_last=-2;
	double delta_sat=0;
	double rho_i=0,Vp_i=0,Vs_i=0,VpHarm_i;
	double rho_im1,Vp_im1,Vs_im1,VpHarm_im1;
	double rhoBar,vpBar,vsBar;
	double delta_rho,delta_vp,delta_vs,ratio,beta,tmp;
	double Dr,Dp,Ds;
	boolean isLastChangingSegment, isLowerSegmentOfPair=false;
	SaturationInfo _sInfo = _WI.saturationInfo;
	double Smon=0,Sbase=0;
	TimeLapseInfo _TLI = m.timelapseinfo;				
	double[] vp_vs_rho = new double[3];
	double[] vp_vs_rhoM = new double[3];

	// n = size of upscaled data
	// n = rho.length+1;
	ArrayList<SegSensLookup> _sensLookups = SMI.sensLookups;

	// recompute reflectivities for segments that "change" only
	for (j=0;j<N;j++)
	  {
		SegSensLookup SSL = _sensLookups.get(j);
		// variable segment is i = SMI.i_s
		// affects refl coefficients top & bottom of this segment. Update both of these
		isLastChangingSegment = (j==(N-1)); 
		for(i=SSL.i_seg; i<= SSL.i_seg+1; i++)
		  {
			// skip lower of adjacent affected segments - we'll get it next time round
			if (!isLastChangingSegment)
			  {
				isLowerSegmentOfPair = (_sensLookups.get(j+1).i_seg == i);
				if (isLowerSegmentOfPair) continue; 
			  }
			else 
			  {
				isLowerSegmentOfPair = (i == (SSL.i_seg+1));
			  }
			if (i==i_last+1)
			  {
				// segment above altered as well. Use its altered props
				rho_im1 = rho_i;
				Vp_im1 = Vp_i;
				Vs_im1 = Vs_i;
			  }
			else
			  {
				rho_im1 = rho[i-1];
				Vp_im1 = Vp[i-1];
				Vs_im1 = Vs[i-1];
			  }
			// saturation change. Might just be fixed by input
			if (after && !(isLastChangingSegment && isLowerSegmentOfPair))
			  {
				Sbase = _sInfo.Sb;
				if (_sInfo.satPropMapper!=null)
				  {
					Smon = _sInfo.satPropMapper.getTrueValue(_sInfo.satProp.value);
				  }
				else
				  {
					Smon = _sInfo.satProp.value;
				  }
				delta_sat = Smon - Sbase;
			  }
			else
			  {
				delta_sat = 0;  // assume Sb corresponds to condition of supplied logs
			  }
			if(false)
			  { // OLD
				// Altered segment properties under saturation shift
				Vp_i = Vp[i] + SSL.dvp * delta_sat;
				Vs_i = Vs[i] + SSL.dvs * delta_sat;
				rho_i = rho[i] + SSL.drho * delta_sat;
			  }
			else
			  {
				// Approximate Gassmann to capture main nonlinearity
				Vp_i = Vp[i];
				Vs_i = Vs[i];
				rho_i = rho[i];
				rho_i*=densityUnitsConversionFactor; // to SI
				if(delta_sat != 0.0)
				  {
					double vp,vs,rho_s,phi_s,M_b,mu_b,NG = SSL.NG;
					//undo rock mixing, assume brine sand, to get
					if(NG > 0.1) // safe threshold
					  {
						rho_s = (rho_i - (1.0-NG)*_TLI.rho_m)/NG; // sand density
						M_b = NG/(1.0/(rho_i*sqr(Vp_i)) - (1.0-NG)/(_TLI.rho_m*sqr(_TLI.vp_m)));  // brinesand M
						mu_b = NG/(1.0/(rho_i*sqr(Vs_i)) - (1.0-NG)/(_TLI.rho_m*sqr(_TLI.vs_m))); // brinesand mu
						vp = sqrt(M_b/rho_s);
						vs = sqrt(mu_b/rho_s);
						// porosity assuming normal brine-adjusted conditions
						phi_s = (_TLI.rho_g-rho_s)/(_TLI.rho_g -_TLI.rho_b);
						if((rho_s<0) || (M_b<0) || (mu_b<0) || (phi_s<0))
						  {
							EH.die_err("Warning: strange log data/time-lapse info predicts "
									   +"rho_sand<0, M_b<0, mu_b<0, or phi_sand<0"
									   +"segment {Vp=%f Vs=%f rho=%f}\n",Vp_i,Vs_i,rho_i);
						  }
						// base
						GeophysicsUtils.GassmannBackus(vp,vs,phi_s,         // new rock, inferred porosity
													   _TLI.vp_m,_TLI.vs_m,_TLI.rho_m, // stock shale
													   NG,                         // segment NG
													   Sbase,                          // base sat
													   _TLI.vp_h,_TLI.rho_h,_TLI.vp_b, _TLI.rho_b, _TLI.K_g, _TLI.rho_g, vp_vs_rho); // stock fluids

						// mon
						GeophysicsUtils.GassmannBackus(vp,vs,phi_s,         // new rock, inferred porosity
													   _TLI.vp_m,_TLI.vs_m,_TLI.rho_m, // stock shale
													   NG,                         // segment NG
													   Smon,                          // mon sat
													   _TLI.vp_h,_TLI.rho_h,_TLI.vp_b, _TLI.rho_b, _TLI.K_g, _TLI.rho_g, vp_vs_rhoM); // stock fluids
						// shift by Gassman effect, downweighted by in-segment fraction. 
						Vp_i  += SSL.wt * (vp_vs_rhoM[0] - vp_vs_rho[0]);
						Vs_i  += SSL.wt * (vp_vs_rhoM[1] - vp_vs_rho[1]);
						rho_i += SSL.wt * (vp_vs_rhoM[2] - vp_vs_rho[2]);
					  }
				  }
			  }
			
			// reflection coefficients
			if (i>0)
			  {
				rhoBar = 0.5*(rho_i+rho_im1);
				vpBar = 0.5*(Vp_i+Vp_im1);
				vsBar = 0.5*(Vs_i+Vs_im1);
				delta_rho = rho_i-rho_im1;
				delta_vp = Vp_i-Vp_im1;
				delta_vs = Vs_i-Vs_im1;
				Dr = delta_rho/rhoBar;
				Dp = delta_vp/vpBar;
				Ds = delta_vs/vsBar;
				
				R_near[i] = 0.5*(Dr + Dp);
				
				// Linearized reflectivity requires:
				tmp = Dr + 2.0*Ds;
				beta = vsBar/vpBar;
				// ratio = delta_vp/(2.0*vpBar) - (2.0*vsBar*vsBar*tmp/(vpBar*vpBar));
				ratio = Dp/2.0 - (2.0*sqr(beta)*tmp);
				
				R1_AVO[i] = ratio;
				Vsq_AVO[i] = sqr(Vp_im1);   // will be dimension-corrected to match stack params later!
				if (doShuey)
				  {
					R2_Shuey_AVO[i] = Dp/2.0; // delta_vp/(2.0*vpBar);					 
				  }
				if (_WI.needPSreflectivities)
				  {
					// Dimensionless linearized expression
					Rps_AVO[i] = -0.5*(Dr + 2.0*beta*tmp);  
					Rps_AVO2[i] = ((1.0 + beta*(8 + 9.0*beta))*Dr  +  beta*(16.0 + 24.0*beta)*Ds)/12.0;  
				  }
				if (dumpPandSreflectivity)
				  {
					deltaRho[i] = 0.5*Dr;
					deltaZs[i] = 0.5*(Dr + Ds);
				  }
			  }
			i_last = i;
		  } // i_s and i_s+1 reflections, top & bottom of affected segment
	  }
  }

  /**
	 Print to stdout successive streams of the raw data and the
	 upscaled (segmented) data. Suitable for graphing, debugging etc
  */
  public void debugPrintSegmentationResults()
  {
	int i,j;
	double[] data=null, upscaledData=null;
	double val,DT_scaling=1.0;
	String p="";

	
	DT_scaling = lData.getSonicLogScaling();

	// print out original and upscaled versions of the three
	// salient quantities for graphing.....
	for (i=0;i<3;i++)
	  {
		// point property where we need it
		switch (i)
		  {	  
		  case 0:
			data = lData.rho;
			upscaledData = rho;
			p="rho ";
			break;
		  case 1:
			data = lData.DT;
			upscaledData = Vp;
			p="Vp ";
			break;
		  case 2:
			data = lData.DTshear;
			upscaledData = Vs;
			p="Vs ";
			break;
		  }
		// Diagnostic output
		// print orig data
		for (j=0; j<data.length; j++) 
		  {
			val = (i==0)? data[j] : DT_scaling/data[j];
			System.out.println(p+j+" "+val);
		  }
		System.out.println(p);
	
		// now upscaled segments
		for (j=0; j<segIndexValues.size()-1; j++)
		  {
			System.out.println(p+(segIndexValues.getQuick(j)-0.5)+" "+upscaledData[j]);
			System.out.println(p+(segIndexValues.getQuick(j+1)-0.5)+" "+upscaledData[j]);
		  }
	  }
  }

  /**
	 Attempt to autodetect velocity units. We expect sonic logs possibly in uS/ft, but
	 stacking velocities may be in metric. Before entering this function, the Vp
	 file will be in (say) ft/s (probably XX)
  */
  void autoDetectVunits(WaveletExtractor _WE, double stackingVelocity, double meanTimeToDepthVelocity)
  {
	int i;
	int cnt=0;
	double Vpbar=0;
	//	 Factor to get AVO velocity ratios with logs <-> Stack velocities right
	double unitsConversionFactor=1.0;

	// compute average Vp of blocked log segments in own units
	for (i=0;i<Vp.length;i++)
	  {
		if (Vp[i]>0)
		  {
			Vpbar += Vp[i];
			cnt++;
		  }
	  }
	Vpbar /= Math.max(cnt,1);


	if (WaveletExtractor.AUTO_DETECT_VELOCITY_UNITS)
	  {
		// now do units detection. Two bits

		// (1) Logs vs stacking velocity

		// Here are some common possibilities

		// Compare Vpbar to stacking V. Should be within a factor of two if units are commensurable
		// If the units are the same, we need to do nothing.

		// But there are two mixed units cases we ought to consider....
		if (((Vpbar*FundamentalConstants.METRES_PER_FOOT) < 2.0* stackingVelocity) && 
			((Vpbar*FundamentalConstants.METRES_PER_FOOT) > 0.5* stackingVelocity))
		  {
			// Usual case: logs in uS/ft, stackingVelocity in m/s:
			// we will have to hit Vsq_AVO with a factor to get it into m^2/s^2
			unitsConversionFactor = FundamentalConstants.METRES_PER_FOOT*FundamentalConstants.METRES_PER_FOOT;
			lengthUnitsConversionFactor = FundamentalConstants.METRES_PER_FOOT;
			if (verbosity > 2)
			  {
				EH.println("Taking sonic logs in uS/ft, stacking velocity in m/s");
			  }
		  }
		if ((Vpbar < 2.0*(stackingVelocity*FundamentalConstants.METRES_PER_FOOT)) && 
			(Vpbar > 0.5*(stackingVelocity*FundamentalConstants.METRES_PER_FOOT)))
		  {
			// Unusual case: logs in uS/m, stackingVelocity in ft/s:
			// we will have to hit Vsq_AVO with a factor to get it into m^2/s^2
			unitsConversionFactor = 1.0/(FundamentalConstants.METRES_PER_FOOT*FundamentalConstants.METRES_PER_FOOT);
			lengthUnitsConversionFactor = 1.0/FundamentalConstants.METRES_PER_FOOT;
			if (verbosity > 2)
			  {
				EH.println("Taking sonic logs in uS/m, stacking velocity in ft/s");
			  }
		  }


		// Updata the Vsq_AVO term if required (the only one that depend on units)
		if (unitsConversionFactor != 1.0)
		  {
			// Vsq_AVO[i] = Vp[i-1] * Vp[i-1];
			// has to be in units comparable to Vdenom
			for (i=0;i<Vsq_AVO.length;i++)
			  {
				Vsq_AVO[i] *= unitsConversionFactor;
			  }

			 if (DeltaR_anis != null)
			   {
				 for (i=0;i<DeltaR_anis.length;i++)
				   {
					 DeltaR_anis[i] *= unitsConversionFactor;
				   }
			   }
		  }

		// (2) Logs vs checkshot

		// Compare Vpbar to stacking V. Should be within a factor of two if units are commensurable
		// If the units are the same, we need to do nothing.

		// But there are two mixed units cases we ought to consider....
		if (((Vpbar*FundamentalConstants.METRES_PER_FOOT) < 2.0* meanTimeToDepthVelocity) && 
			((Vpbar*FundamentalConstants.METRES_PER_FOOT) > 0.5* meanTimeToDepthVelocity))
		  {
			// Usual case: logs in uS/ft, meanTimeToDepthVelocity in m/s:
			// we will have to hit Vsq_AVO with a factor to get it into m^2/s^2
			lengthUnitsConversionFactorCheckshot = FundamentalConstants.METRES_PER_FOOT;
			if (verbosity > 2)
			  {
				EH.println("Taking sonic logs in uS/ft, checkshot velocities in m/s");
			  }
		  }
	
		if ((Vpbar < 2.0*(meanTimeToDepthVelocity*FundamentalConstants.METRES_PER_FOOT)) && 
			(Vpbar > 0.5*(meanTimeToDepthVelocity*FundamentalConstants.METRES_PER_FOOT)))
		  {
			// Unusual case: logs in uS/m, meanTimeToDepthVelocity in ft/s:
			// we will have to hit Vsq_AVO with a factor to get it into m^2/s^2
			lengthUnitsConversionFactorCheckshot = 1.0/FundamentalConstants.METRES_PER_FOOT;
			if (verbosity > 2)
			  {
				EH.println("Taking sonic logs in uS/ft, checkshot velocities in ft/s");
			  }
		  }
	  }
  }
  /**
	 Attempt to autodetect density units: either kg/m^3 or gm/cc
  */
  void autoDetectRhounits()
  {
	int i;
	int cnt=0;
	double rhobar=0;

	// compute average rho of blocked log segments in own units
	for (i=0;i<rho.length;i++)
	  {
		if (rho[i]>0)
		  {
			rhobar += rho[i];
			cnt++;
		  }
	  }
	rhobar /= Math.max(cnt,1);
	if(rhobar<10.0)
	  {
		//	assume gm/cc
		densityUnitsConversionFactor = 1e3;
	  }
  }


  /**
	 Print to filename a simple ASCII version of the segmented data.
	 To be called after segmentLogSuite(data) has been invoked, so the
	 blocking has been done, and after the t2dMap and survey objects
	 are populated

	 t2dMap should be passed in as _WE.m.timeToDepthMaps[i]
	 survey should be passed in as _WE.logData[i].survey
  */
  public void writeOutBlockedLogFile(String filename, TimeToDepthMap t2dMap, boolean registrationErrorsAreCommon)
  {
	int i,j,js,je,jen,k;
	double[] data=null, upscaledData=null;
	double val, Z,DT_scaling;

	boolean doTimes = true;
	double t = FundamentalConstants.BAD_NUMBER;
	double zTVDmin=0,zTVDmax=0;
	LinearInterpolation tInterp=null;
	boolean writeSeismicMinicubeTime = true;
	double MLoffset_val=0;


	DT_scaling = lData.getSonicLogScaling();

	if (t2dMap.Z == null) doTimes=false;
	if (doTimes)
	  {
		// time to TVD mapping 
		tInterp = new LinearInterpolation(t2dMap.Z, t2dMap.t, t2dMap.Z.length,EH);
		// but don't ask for times outside this range....
		zTVDmin = t2dMap.Z[0];
		zTVDmax = t2dMap.Z[t2dMap.Z.length-1];
		// minicube seismic scale
		if (registrationErrorsAreCommon)
		  {
			MLoffset_val = t2dMap.systematicTProp.value;		
		  }
		else
		  {
			// 2 different offsets. Can't absorb both into an "effective" t2d conversion
			writeSeismicMinicubeTime = false;
		  }
	  }



	// Write as a pseudo GEOEAS file
	PrintWriter out;
	try 
	  {
		out = new PrintWriter(new FileOutputStream(filename));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("Writing segmented log data: couldn't open " + filename);
	  }
	if (verbosity>1) EH.println("Writing segmented log data as "+filename);

	out.printf("Segmented log data in GEOEAS Format%n");
	// out.printf("10%n");
	// Added time...
	out.printf("%d%n", (writeSeismicMinicubeTime)? 12 : 11);
	out.printf("MD%n");
	out.printf("TVD%n");
	out.printf("t_checkshot(ML)%n");
	if (writeSeismicMinicubeTime) out.printf("t_seismic%n");
	out.printf("rho%n");
	out.printf("rho_blocked%n");
	out.printf("Vp%n");
	out.printf("Vp_blocked%n");
	out.printf("Vs%n");
	out.printf("Vs_blocked%n");
	out.printf("AI(computed)%n");
	out.printf("AI_blocked(computed)%n");	
	int kBound = Math.min(Math.min(rho.length,Vp.length),Vs.length);
	int jBound = Math.min(Math.min(lData.rho.length,lData.DT.length),lData.DTshear.length);
 		
	LinearInterpolation zInterp = 
	  new LinearInterpolation(lData.survey.Zm, lData.survey.Z, lData.survey.Zm.length,EH);
	int seg_size_minus1 =  segIndexValues.size()-1; 
	for (j=0,k=0; j<jBound && k<seg_size_minus1; j++) 
	  {
		// 			 Get corresponding index in segmented suite
		// The depth lData.Zm[0,1...] has segments for the range of indices
		// [segIndexValues.getQuick(k)+1,segIndexValues.getQuick(k+1)]
		//if (segIndexValues.getQuick(k) <= j)
		for(;k<seg_size_minus1; k++)
		  {
			js = (k>1) ? segIndexValues.getQuick(k) + 1 : 0;
			je = segIndexValues.getQuick(k+1);
			if (js<=j && j<=je) break;
		  }
		// equivalently, but slower...
		// k = getLowerBracketingIndex(j,segIndexValues); 
		if (k<=0 || k>= kBound) continue;

		// Attempt to get Z = TVD from lData.Zm at sample j, md = lData.Zm[j]
		Z = zInterp.interpolate(lData.Zm[j],j!=0);
		// time from TVD
		if (doTimes) 
		  {
			if ((Z>=zTVDmin) && (Z<=zTVDmax))
			  {
				t = tInterp.interpolate(Z,j!=0);
			  }
			else t = FundamentalConstants.BAD_NUMBER;
		  }
		// data with velocities computed from DT's
		// MD rho rho_blocked Vp Vp_blocked Vs Vs_blocked
		// MD, TVD, chksht time
		out.printf("%10.6g %10.6g %10.6g ",
				   lData.Zm[j], Z, t);
		// seismic time
		if (writeSeismicMinicubeTime) out.printf("%10.6g ", t + MLoffset_val);
		// density etc
		out.printf("%10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g%n",
				   lData.rho[j], rho[k],
				   DT_scaling/lData.DT[j], Vp[k],
				   DT_scaling/lData.DTshear[j], Vs[k],
				   lData.rho[j]*DT_scaling/lData.DT[j],
				   rho[k]*Vp[k]);
	  }
	out.close();
	if (verbosity>1) EH.println("Wrote "+filename);
  }

  private int getLowerBracketingIndex(double x, IntArrayList X)
  {
	int k,klo,khi;
	int size = X.size();
	// bounding index in checkshot array.
	if (x <X.getQuick(0)) return -1;
	if (x > X.getQuick(X.size()-1)) return size;
	// bisection to find bracket
	klo = 0;
	khi = size-1;
	while (khi-klo > 1) 
	  {
		k=(khi+klo) >> 1;
		if (X.getQuick(k) > x) khi=k;
		else klo=k;
	  }
	// now X[klo] <= x < X[khi]  && khi=klo+1
	return klo;
  }

  /**
	 Class for extra information about saturation dependencies of properties for later vintages
  */
  static class SaturationModelInfo
  {
	ArrayList<SegSensLookup> sensLookups;

	public SaturationModelInfo()
	{
	  sensLookups  = new 	ArrayList<SegSensLookup>();
	}
  }

  public final static class SegSensLookup
  {
	int i_seg;
	double dvp,dvs,drho;
	double wt,NG;
	public SegSensLookup(int i, double[] slope,double _wt, double _NG)
	{
	  i_seg = i;
	  dvp=slope[0];
	  dvs=slope[1];
	  drho=slope[2];
	  wt=_wt;
	  NG=_NG;
	}
  }
}
