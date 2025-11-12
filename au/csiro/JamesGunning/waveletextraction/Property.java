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

public class Property
{
  double value;
  boolean normalised=false;
  boolean mustBePositive=false;
  double mean;
  double sigma=0.0;

  public Property(double value, boolean normalised, boolean mustBePositive)
  {
	this.value = value;
	this.mean = value;
	this.normalised = normalised;
	this.mustBePositive = mustBePositive;
  }


  /**
	 simple setter for positive physical quantitites
  */

  void setPosValue(double value)
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
	 simple setter for general quantitities
  */

  void setValue(double value)
  {
	if (normalised)
	  {
		if (value < 0.0)
		  {
			this.value=0.0;
		  }
		else if (value>1.0)
		  {
			this.value=1.0;
		  }
		else this.value=value;
	  }
	else
	  {
		if (mustBePositive)
		  {
			this.value = Math.abs(value);
		  }
		else
		  {
			this.value =value;
		  }
	  }
	return;
  }

  /**
	 Set the mean field by passing in a real number. Passing in Double.NaN will be passed
	 over without action. Also sets the value, as this is a useful default thing to do.
  */
  public void setMean(double mean)
  {
	if (!Double.isNaN(mean))
	  {
		this.mean = mean;
		this.value = mean;
	  }
  }
  
  /**
	 Set the sigma field by passing in a real number. Passing in Double.NaN will be passed
	 over without action
  */
  public void setSigma(double sigma)
  {
	if (!Double.isNaN(sigma)) this.sigma = sigma;
  }

  /**
	 Set the value to the mean
  */
  public void resetToMean()
  {
    value = mean;
  }
}
