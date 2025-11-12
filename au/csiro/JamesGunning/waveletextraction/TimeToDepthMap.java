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
/**
   Data class for storing time to depth information. Not a very general
   design, but efficient.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class TimeToDepthMap
{
  WaveletExtractor _WE;

  /**
	 Use this switch to ensure sane ordering as the map flexes. The need for this
	 no longer exists since we usually run with WE.constrainIntervalVelocities = true
  */
  public static final boolean APPLY_VELOCITY_BASED_PENALTIES = false;

  /**
	 Threshold time for ordering of time-depth map
  */
  public static final double TIME_ORDERING_THRESHOLD_TIME = 0.001*FundamentalConstants.BASE_SCALING_TIME;
  /** 
	  fixed penalty for misordered times in velocity-based penalties
  */
  public static final double MISORDERING_PENALTY_COEFF_VELOCITY_BASED = 1.0;
  /** 
	  fixed penalty coefficient (MISORDERING_PENALTY_COEFF x |delta_t|) for misordered times
  */
  public static final double MISORDERING_PENALTY_COEFF=100;
  /**
	 true depth TVD
  */
  double[] Z;
  /**
	 corresponding time
  */
  double[] t;

  /** 
	  approximate velocities corresponding to initialised map
  */
  double[] scalingVelocity;

  /** 
	  approximate velocities from upscaled sonic logs in interval.
	  Set to scalingVelocity if logs unavailable
  */
  double[] upscaledLogVelocity;

  /**
	 std deviations
  */
  double[] sigma_t;

  /**
	Array of properties to handle times
  */
  Property[] tProp;

  // Property controlling systematic registration offsets if
  Property systematicTProp = null;

  private boolean lastTimesWereDecentlyOrdered=true;

  private double tStartCheck, tEndCheck;

  private static boolean debug = !false;
  /**
	 Constructor. Sets up initial time to depth map from sorted & combined
	 Checkshot and Marker data
  */
  public TimeToDepthMap(WaveletExtractor WE,TimeToDepthMap.TimeToDepthTriples[] t2d, int well_indx)
  {
	int N;
	
	this._WE = WE;

	N = t2d.length;

	if (_WE.m.registrationErrorsAreCommon)
	  {
		int synchIndx = _WE.wellInfo[well_indx].synchronizeWellIndex;
		systematicTProp = _WE.seismicInfo[synchIndx][0].systematicTProp;
	  }

	Z = new double[N];
	t = new double[N];
	sigma_t = new double[N];
	scalingVelocity = new double[N-1];
	upscaledLogVelocity = new double[N-1];
	// initialise values to mean from checkshots + markers
	tProp = new Property[N];
	for (int i=0; i<N ; i++)
	  {
		Z[i] = t2d[i].Z_mean;
		t[i] = t2d[i].T_mean;
		sigma_t[i] = t2d[i].sigma_T;
		// rough velocity of initial time to depth map
		if (i>0) 
		  {
			if ((Z[i]<=Z[i-1]) || (t[i]<=t[i-1])) 
			  {
				throw new RuntimeException("Time to depth map instantiated in non-monotonic state");
			  }
			// factor in two-way times
			scalingVelocity[i-1] = 2.0*(Z[i]-Z[i-1])/(t[i]-t[i-1]);
			// put in sane defaults for these too...
		  }
		//And array of properties for times
		tProp[i] = new Property(t[i],false,false);
	  }
	// put in sane defaults for upscaled velocity too...
	System.arraycopy(scalingVelocity,0,upscaledLogVelocity,0,N-1);

	//	for (int i=0; i<t2d.length; i++) System.out.println("TimeToDepthMap cnstr Tmean: "+t2d[i].T_mean);
  
	// Safe endpoint bounds on maps... 
	tStartCheck = (_WE.assumeNoReflectionsOutsideLogs)? (_WE.tStart - _WE.tCodaMax) : _WE.tStart;
	tEndCheck = (_WE.assumeNoReflectionsOutsideLogs) ? (_WE.tEnd - _WE.tPrecursorMin) : _WE.tEnd;
  }

  /**
	 Set the accessible t field values 
	 to match the current values in
	 tProp. Sets lastTimesWereDecentlyOrdered to false
	 if times were badly ordered.
	 Returns penalty function reflecting misordering criterion.
	 The penalty is zero if the times are monotonically ordered.

	 FIXME!! check if this is a sane thing to do 
	 Important ordering constraints here:
	 (1) points monotonic and non-identical
  */
  double updateTimesAndReturnMisorderingPenalty()
  {
	int i,N;
	double layerPenalty=0,penaltyCoeff=1.0,penalty=0.0;
	double delta_t;

	N = t.length;

	for (i=0; i<N ; i++)
	  {
		t[i] = tProp[i].value;
		if (debug)
		  {
			if (Double.isNaN(t[i]))
			  {
				_WE.EH.println("Nan's in t2D map!");
			  }
		  }	
	  }
	lastTimesWereDecentlyOrdered = true;
	// Heavy-handed monotonic ordering for now. FIXME. Penalties??????
	for (i=1; i<N ; i++)
	  {
		if (APPLY_VELOCITY_BASED_PENALTIES)
		  {
			// End points. Have to be in range. Initial checks guarantee that in-rangedness
			// is a reasonable thing to enforce
			// (eg first+last checkshots 4 prior-std-deviations outside prior mean...)
			if (i==1 && (t[0] > tStartCheck))
			  {
				delta_t = t[0] - tStartCheck;
				t[0] = tStartCheck;
			  }
			else if (i==(N-1) && (t[N-1] < tEndCheck))
			  {
				// final time
				delta_t = tEndCheck - t[N-1];
				t[N-1] = tEndCheck;
			  }
			else
			  {
				delta_t = t[i]-t[i-1];
			  }
			penaltyCoeff=1.0;
			if (delta_t <= TIME_ORDERING_THRESHOLD_TIME) 
			  {
				// Misordered
				// Absolutely *must* correct ordering in any case
				t[i] = t[i-1] + TIME_ORDERING_THRESHOLD_TIME;
				lastTimesWereDecentlyOrdered = false;
				// falling ramp additional penalty
				penaltyCoeff = (1.0 - MISORDERING_PENALTY_COEFF_VELOCITY_BASED*
								(delta_t-TIME_ORDERING_THRESHOLD_TIME)/FundamentalConstants.BASE_SCALING_TIME);
				delta_t = TIME_ORDERING_THRESHOLD_TIME;
			  }
			// normalised deviate of implied new "effective' velocity cf initialised velocity
			double delta_v = (2.0*(Z[i]-Z[i-1])/delta_t - scalingVelocity[i-1])/scalingVelocity[i-1];
			layerPenalty = penaltyCoeff*delta_v*delta_v;
		  }
		else
		  {
			// Simple Absolute time based penalties
			// End points. Have to be in range. Initial checks guarantee that in-rangedness
			// is a reasnable this to enforce
			// (eg first+last checkshots 4 prior-std-deviations outside prior mean...)
			if (t[0] > tStartCheck)
			  {
				penalty += MISORDERING_PENALTY_COEFF * (t[0] - tStartCheck);
				t[0] = tStartCheck;
				if (_WE.verbosity >1) _WE.EH.println("T2D map startpoint penalty added");
			  }
			// ordering
			if (t[i]<t[i-1]) 
			  {
				// simple misordering
				layerPenalty = MISORDERING_PENALTY_COEFF * Math.abs(t[i] - t[i-1]);
				t[i] = t[i-1] + FundamentalConstants.BASE_SCALING_TIME;
				lastTimesWereDecentlyOrdered = false;
			  }
			else 
			  {
				layerPenalty=0.0;
			  }
			// final time
			if (t[N-1] < tEndCheck)
			  {
				penalty += MISORDERING_PENALTY_COEFF * (tEndCheck - t[N-1]);
				t[N-1] = tEndCheck;
				if (_WE.verbosity >1) _WE.EH.println("T2D map endpoint penalty added");
			  }
		  }
		penalty += layerPenalty;
	  }
	if (_WE.verbosity > 4)
	  {
		if (penalty>0)
		  {
			_WE.EH.println("Misordered times...(orig: new)");
			for (i=0; i<t.length ; i++)
			  {
				_WE.EH.println(tProp[i].value+ " "+t[i]);
			  }
		  }
	  }

	return penalty;
  }
  /**
	 Returns flag to indicate if the last set of times fed to updateTimesAndReturnMisorderingPenalty()
	 (which must be called before this) were badly ordered.
  */
  boolean areCurrentTimesDecentlyOrdered()
  {
	return lastTimesWereDecentlyOrdered;
  }


  /**
	 A set of triples that will 
	 hold the combined marker and checkshot data
  */
  public static class TimeToDepthTriples implements Comparable 
  {
	/**
	   After setUpCheckShotAndMarkerData(), will contain
	   Z (TVD) value
	*/
	double Z_mean;

	/**
	   After setUpCheckShotAndMarkerData(), will contain
	   T value
	*/
	double T_mean;
	/**
	   After setUpCheckShotAndMarkerData(), will contain
	   sigmaT value
	*/
	double sigma_T;

	/**
	   Sorter to ensure we can get data in order of ascending depths
	*/
	public int compareTo(Object o)
	{
	  TimeToDepthTriples n = (TimeToDepthTriples)o;
	  int cmp = (Z_mean<n.Z_mean)? -1 : ((Z_mean>n.Z_mean)? 1 : 0 );
	  return cmp;
	}
  }
}
