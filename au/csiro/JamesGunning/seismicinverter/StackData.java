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

import au.csiro.JamesGunning.lib.util.*;

public final class StackData
{
  public String filename;
  public String wavelet_filename;
  public double Xmin;
  public double Xmax;
  private double tAbove=Double.NaN;
  private double tBelow=Double.NaN;
  private double T_event;
  private double V_stack;
  protected ModelPriorDescriptions.PriorDetails noise_rms;
  double log_noise_rms;

  boolean isPSdata = false;                       // if PS data in PP time
  double theta_specified = -Double.MAX_VALUE;     // specified stack angle if constant stack
  double sq_theta_specified = -Double.MAX_VALUE;  // specified theta^2 of stack angle theta if constant stack

  // pointer to current data. can be useful for efficiency
  float[] currentTrace;

  public StackData() 
  {
  }
  
  public void setTime(double T_event)
  {
	if (T_event>0)
	  this.T_event = T_event;
	else
	  throw new RuntimeException("Stack event-time must be positive");
  }

  public double getTime() {return T_event;}

  public void setVstack(double V_stack)
  {
	if (V_stack>0)
	  this.V_stack = V_stack;
	else
	  throw new RuntimeException("Stack velocity must be positive");
  }

  /**
	 Set the optional tAbove field by passing in a real number. Passing in Double.NaN will
	 not override the wavelet-dictated  defaults
  */
  public void setTabove(double tAbove)
  {
	if (tAbove!=Double.NaN) this.tAbove = tAbove;
  }
 /**
	Return the stack's tAbove value in seconds
  */
  public double getTabove()
  {
	return tAbove;
  }

  /**
	 Set the optional tBove field by passing in a real number. Passing in Double.NaN will
	 not override the wavelet-dictated  defaults
  */
  public void setTbelow(double tBelow)
  {
	if (tBelow!=Double.NaN) this.tBelow = tBelow;
  }

  /**
	Return the stack's tBelow value in seconds
  */
  public double getTbelow()
  {
	return tBelow;
  }

  public double getVstack() 
  {
	return V_stack;
  }
  
  public void setStackRange(double Xmax, double Xmin)
  {
	if (Xmin>Xmax)
	  throw new RuntimeException("Stack must have Xmin<=Xmax");
	else
	  {
		this.Xmax=Xmax;
		this.Xmin=Xmin;
	  }
  }


  public double getOffsetMeanSquare()
  {
	if (Xmin==Xmax)
	  return Xmin*Xmin;
	else
	  return (Xmax*Xmax*Xmax - Xmin*Xmin*Xmin)/(3.0*(Xmax - Xmin));
  }

  

  public boolean isStackAnglePrescribed()
  {
	return (theta_specified >= 0);
  }

  void computeStackParams()
  {
	double Xms;

	if (theta_specified >= 0)
	  {
		// compute this useful guy
		sq_theta_specified = MathUtil.sqr(theta_specified);		
	  }
  }


  public void setNoise(double noise_rms)
  {
	if (noise_rms<0)
	  throw new RuntimeException("Stack noise-rms must be positive");
	else
	  this.noise_rms.value = noise_rms;
  }
  
  public void computeLogNoise()
  {
	if (noise_rms.value > 0) log_noise_rms = Math.log(noise_rms.value);
  }

  public double getNoise()
  {
	return noise_rms.value;
  }

}

