/*
  Copyright 1999 CERN - European Organization for Nuclear Research.
  Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose 
  is hereby granted without fee, provided that the above copyright notice appear in all copies and 
  that both that copyright notice and this permission notice appear in supporting documentation. 
  CERN makes no representations about the suitability of this software for any purpose. 
  It is provided "as is" without expressed or implied warranty.
*/
package au.csiro.JamesGunning.lib.util;

/**
 * Constructs functions that returns
 * <ul>
 * <li><tt>a * a</tt>
 * </ul> 
 * <tt>a</tt> is variable.
 * Intended to be passed to <tt>matrix.assign(function)</tt> methods.
 */
public final class PowDoubleFunction implements cern.colt.function.DoubleFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  public double exponent;

 /**
   * Insert the method's description here.
   * Creation date: 21 Aug 2002
   */
  protected PowDoubleFunction(final double exponent) 
  {
	this.exponent = exponent;
  }

  /**
   * Returns the result of the function evaluation.
   */
  public final double apply(double a) 
  {
	if (exponent==2.0)
	  {
		return a * a;
	  }
	else if (exponent==3.0)
	  {
		return a * a * a;
	  }
	else
	  {
		return Math.pow(a,exponent);
	  }
  }

  /**
   * <tt>a ^ exponent</tt>.
   */
  public static PowDoubleFunction pow(final double exponent) {
	return new PowDoubleFunction(exponent);
  }
}
