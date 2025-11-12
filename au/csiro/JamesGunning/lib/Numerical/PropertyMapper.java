/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
// USE_NEW_COPYRIGHT


package au.csiro.JamesGunning.lib.Numerical;

/**
 * PropertyMapper.java
 *
 * A simple class to take properties that may exist only on finite ranges (maxValue,maxValue)
 * and remap them to (-infty,infty) via certain standard tranforms.
 * The class provides: no-transform, the logit, and the logarithm.
 * The no-transform mapType makes an object which returns transformedValue = trueValue always.
 * The logarithm mapType makes an object with transformedValue = log((trueValue - minValue)/defaultAndScalingValue)
 * I.e. maxValue = infty
 * The logit  mapType makes an object with transformedValue = logit(u) = log(u/(1-u))
 *  where u = ((trueValue - minValue)/(maxValue-minValue))
 *
 * The absolute map simply returns transformedValue = |trueValue|
 * 
 * Created: Mon Apr 22 09:49:17 2002
 *
 * @author James Gunning
 * @version
 */

public class PropertyMapper
{
  /**
	 A convenience constant. If L(u) = log(u/(1-u)), (logit transform) and
	 L ~ N(0,s^2), 
	 with s = LOGIT_NORM_TO_UNIFORM_APPROX_STDDEV,
	 the inverse mapping
	 u=(1+exp(-L)) makes u nearly uniform: u ~ U(0,1).
	 If s is any wider than this, the distribution gets two-humped.
   */
  public static double LOGIT_NORM_TO_UNIFORM_APPROX_STDDEV = 1.4;

  public static final short NO_TRANSFORM = 0;
  public static final short LOGIT = 1;
  public static final short LOGARITHM = 2;
  public static final short ABS = 3;

  /**
	 default value used in scaling exponential transforms
  */
  private double defaultAndScalingValue=1.0;
  /**
	 outer ranges of logit transformed variable
  */
  private double maxValue=0;
  private double minValue=1;

  /**
	 Map type key.<br>
	 0=no transform<br>
	 1=logit;<br>
	 2=logarithm<br>
  */
  private short mapType;

  // reference to actual property in question
  public int propertyMasterIndx;

  // flag to indicate whether this PropertyMapper is the controller of an array of coupled PropertyMappers
  public boolean isACoupledMapper = false;
  public PropertyMapper[] coupledProperties=null;
  public double[] coupledCoeffs=null;
  public double[] coupledConst=null;

  // instantiator for a notransform map
  public PropertyMapper(int propertyMasterIndx)
  { 
	this.propertyMasterIndx = propertyMasterIndx;
	this.mapType = NO_TRANSFORM;
  }

 // instantiator for an abs map
  public PropertyMapper(int propertyMasterIndx, short mapType)
  { 
	if (mapType != ABS) 
	  throw new RuntimeException("Illegal mapType  when creating PropertyMapper");
	this.propertyMasterIndx = propertyMasterIndx;
	this.mapType = ABS;
  }

  // instantiator for a logit or logarithm map
  public PropertyMapper(int propertyMasterIndx, short mapType, double value1, double value2)
  { 
	this.propertyMasterIndx = propertyMasterIndx;
	if (mapType == LOGIT)
	  {
		this.mapType = LOGIT;
		this.minValue = value1;
		this.maxValue = value2;
	  }
	else if (mapType == LOGARITHM)
	  {
		this.mapType = LOGARITHM;
		this.defaultAndScalingValue = value1;
		this.minValue = value2;
	  }
	else
	  {
		throw new RuntimeException("Illegal mapType  when creating PropertyMapper");
	  }
  }

  // instantiator for a logit or logarithm map without propertyMasterIndx
  public PropertyMapper(short mapType, double value1, double value2)
  { 
	this.propertyMasterIndx = propertyMasterIndx;
	if (mapType == LOGIT)
	  {
		this.mapType = LOGIT;
		this.minValue = value1;
		this.maxValue = value2;
	  }
	else if (mapType == LOGARITHM)
	  {
		this.mapType = LOGARITHM;
		this.defaultAndScalingValue = value1;
		this.minValue = value2;
	  }
	else
	  {
		throw new RuntimeException("Illegal mapType  when creating PropertyMapper");
	  }
  }

  public double getTrueValue(double transformedValue)
  {
	switch (mapType)
	  {
	  case NO_TRANSFORM:
		return transformedValue;

	  case LOGIT:
		// inverse standard logit
		double u;
		// care with overflow
		if (transformedValue<0)
		  {
			u = Math.exp(transformedValue);
			u = u/(1.0+u);
		  }
		else
		  {
			u = 1.0/(1.0+Math.exp(-transformedValue));
		  }
		// rescaled
		u = u*(maxValue-minValue) + minValue;
		return u;

	  case LOGARITHM:
		return Math.exp(transformedValue)*defaultAndScalingValue+minValue;
	  case ABS:
		return Math.abs(transformedValue);
	  }
	// never gets here
	return 0;
  }

  public double getTransformedValue(double trueValue)
  {
	switch (mapType)
	  {
	  case NO_TRANSFORM:
		return trueValue;

	  case LOGIT:
		double u = (trueValue - minValue)/(maxValue-minValue);
		// must protect against overflow here
		if (u<=0) 
		  {
			return  Math.log(Double.MIN_VALUE);
		  }
		else if (u>=1) 
		  {
			return  Math.log(Double.MAX_VALUE);
		  }
		else
		  {
			return  Math.log(u/(1.0-u));
		  }
	  case LOGARITHM:
		return Math.log((trueValue-minValue)/defaultAndScalingValue);
	  case ABS:
		return Math.abs(trueValue);
	  }
	// never gets here
	return 0;
  }

  /**
	 Return slope of transform at true trueValue, computed bt finite diffs.
  */
  public double getFDSlope(double trueValue)
  {
	double dx,slope;
	
	if (mapType == ABS) return ((trueValue>=0) ? 1: -1);
	
	dx =  trueValue*1e-7;  // usual sqrt(eps) rule for fwd differences
	
	if (trueValue+dx < maxValue)
	  {
		slope = (getTransformedValue(trueValue+dx) - getTransformedValue(trueValue))/dx;
	  }
	else
	  {
		slope = (getTransformedValue(trueValue-dx) - getTransformedValue(trueValue))/(-dx);
	  }
	return slope;
  }

  /** 
	  Simple allocations for a coupled property mapper. Population of arrays left to explicit external control
  */
  public void makeCoupled(int size)
  {
	if (size>0)
	  {
		isACoupledMapper = true;
		coupledProperties = new PropertyMapper[size];
		coupledCoeffs = new double[size];
		coupledConst = new double[size];
	  }
  }
}
