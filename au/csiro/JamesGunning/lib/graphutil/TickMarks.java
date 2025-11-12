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

package au.csiro.JamesGunning.lib.graphutil;

import au.csiro.JamesGunning.lib.*;

/**
   Useful class for generating tick marks for graphs. 
   Kludged from GNU plotuils ideas
   
   @author James.Gunning@csiro.au
   @version 1.0, 1/2/2004

*/
public class TickMarks
{
  /**
	 Useful array of small tick marks. Not used directky in this class, but elsewhere.
  */
  public static double[] LOG_SMALL_TICKS = {0.30103, 0.477121, 0.60206, 0.69897, 0.778151, 0.845098, 0.90309, 0.954243};

  public static int TRIAL_NUMBER_OF_TICK_INTERVALS = 5;
  /**
	 valid interval lengths 
  */
  public static double[] vint = {1.0, 2.0, 5.0, 10.0};
  private int tick_spacing_type;
  private double range;
  private int min_tick_count, max_tick_count;


  /** 
	  Corresponding breakpoints.  The published algorithm uses geometric
	 means, i.e. sqrt(2), sqrt(10), sqrt(50), but using sqrt(10)=3.16...
	 will (if nticks=5, as we choose it to be) cause intervals of length
	 1.5 to yield an inter-tick distance of 0.2 rather than 0.5.  So we
	 could reduce it to 2.95.  Similarly we could reduce sqrt(50) to 6.95
	 so that intervals of length 3.5 will yield an inter-tick distance of
	 1.0 rather than 0.5. 
  */
  public static double[] sqr = {1.414159, 3.16228, 7.07107};



  // Other floating point stuff & macros
  public static double EPS=0.000001;

  private ErrorHandler _EH;
  /**
	 Constructor: provide a standarad error handler
  */
  public TickMarks(ErrorHandler EH)
  {
	this._EH = EH;
  }


  /* Ideas from GNU plotutils:
   * Algorithm SCALE1, for selecting an inter-tick spacing that will yield a
   * good-looking linear-format axis.  The spacing is always 1.0, 2.0, or 5.0
   * times a power of ten.
   *
   * Reference: Lewart, C. R., "Algorithms SCALE1, SCALE2, and
   *	SCALE3 for Determination of Scales on Computer Generated
   *	Plots", Communications of the ACM, 10 (1973), 639-640.  
   *      Also cited as ACM Algorithm 463.
   *
   * We call this routine even when the axis is logarithmic rather than
   * linear.  In that case the arguments are the logs of the actual
   * arguments, so that it computes an optimal inter-tick factor rather than
   * an optimal inter-tick spacing.
   * Returns: Inter-tick spacing
   * if trial_number_of_ticks <=0, defaults are used
   */

  public double getTickSpacing(double min, double max, int trial_number_of_ticks)
	//  int *getTickSpacing_type;	/* Inter-tick spacing type (0, 1, or 2,  i.e. S_ONE, S_TWO, or S_FIVE) */
  {
	int k;
	double nal;
	double a, b;
	double tick_spacing;

	/* compute trial inter-tick interval length */
	if (max==min)
	  {
		tick_spacing_type=0;
		tick_spacing=1;
	  }
	else
	  {
		if (trial_number_of_ticks<=0)
		  {
			a = (max - min) / TRIAL_NUMBER_OF_TICK_INTERVALS;
		  }
		else
		  {
			a = (max - min) / trial_number_of_ticks;
		  }
		a *= (max > min) ? 1.0 : -1.0; /* paranoia, max>min always */
		if (a <= 0.0)			
		  {
			_EH.die_err("error: bad trial inter-tick spacing "+a);
		  }
		nal = Math.floor(Math.log(a)/2.30259);  // log10(a)
		b = a * Math.pow (10.0, -nal);	/* 1.0 <= b < 10.0 */

		/* round to closest permissible inter-tick interval length */
		k = 0;
		do
		  {
			if (b < sqr[k])
			  break;
			k++;
		  }
		while (k < 3);

		tick_spacing = (max > min ? 1.0 : -1.0) * vint[k] * Math.pow (10.0, nal);
		/* for increment type, k=0,1,2 means 1,2,5 times a power of 10 */
		tick_spacing_type = (k == 3 ? 0 : k);
	  }
	range = max - min;
	min_tick_count = (int)(Math.floor((min + EPS * range)/tick_spacing));
	max_tick_count = (int)(Math.ceil((max - EPS * range)/tick_spacing));
	return tick_spacing;
  }

  /**
	 Ruturns tick type: 0,1,2 means 1,2,5 times a power of 10
  */
  public int getTickType()
  {
	int y;
	return tick_spacing_type;
  }
  /**
	 Return tick range
  */
  public double getRange() 
  {
	return range;
  }
  /**
	 Return min tick count
  */
  public int getMinTickCount() 
  {
	return min_tick_count;
  }
  /**
	 Return max tick count
  */
  public int getMaxTickCount() 
  {
	return max_tick_count;
  }
}
