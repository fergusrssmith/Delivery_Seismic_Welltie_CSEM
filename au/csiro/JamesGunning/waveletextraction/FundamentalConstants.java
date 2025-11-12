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
	 The inversion classes require some fundamental constants for 
	 time constants etc.
  */

public final class FundamentalConstants
{
  /**
	 All magic time constants in the inverter are scaled according
	 to this time. This is a useful constant for most surface based 
	 geophysical seismic work. Successive doubles of this number
	 should be canonical sampling times for usual geophysical work, e.g.
	 1mS, 2mS, 4mS...
  */ 
  public static final double BASE_SCALING_TIME = 0.001;
  /**
   Sonic Logs usually have microseconds units. Whatever it is, this is the
   factor we use to put the time part of the units into seconds.
  */ 
  public static final double SONIC_LOG_TIME_SCALING = 1e6;
  /**
   Factor for converting ft/s to m/s
  */ 
  public static final double METRES_PER_FOOT = 0.3048;
  /**
   We may not have Vs logs. Fake these by Vs = Vp/factor, where
   this is the factor. A simple hack in lieu of good petrophysics.
  */ 
  public static final double FAKE_VS_FACTOR = 1.7;
  /**
   Signing "Unknown" or "Bad Value" entries in logs is performed
   by setting them to this number (or less):
  */
  public static final double BAD_NUMBER = -999.0;
  /**
	 log values etc have to be +ve to be useful. Hence this.
  */
  public static final double POS_THRESHOLD = 0;
  /**
	 detection constant to detect time value in ms rather than seconds
  */
  public static final double TIME_UNITS_AUTO_DETECTION_TIME_MAX = 20;
  /**
	 Magic number for checking some scales based on reflection coefficients
  */
  public static final double TYPICAL_REFLECTION_COEFF=0.1;
  /**
	 Magic number for setting scale of prior standard deviation of wavelet
	 coefficients from TYPICAL_REFLECTION_COEFF and seismic RMS
  */
   public static final double TYPICAL_PEAK_TO_RMS_SEISMIC_RATIO = 3.0;

  /**
	N(std dvs) above & below t2d uncertainty to allow for fluctuations in getting depth
	associated with tStart, tEnd	
  */
  public static final double T_UNCERTAINTY_N_STD_DEVS_FACTOR = 4.0;

  /**
	 Converts a time to seconds if required
  */
  public static double autoTimeConvertMsToS(double time)
  {
	if (Math.abs(time)>TIME_UNITS_AUTO_DETECTION_TIME_MAX)
	  {
		return 0.001*time;
	  }
	else
	  {
		return time;
	  }
  }
}
