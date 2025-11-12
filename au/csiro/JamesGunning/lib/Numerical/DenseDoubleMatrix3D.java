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


/*
Original notice.......
Copyright 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose 
is hereby granted without fee, provided that the above copyright notice appear in all copies and 
that both that copyright notice and this permission notice appear in supporting documentation. 
CERN makes no representations about the suitability of this software for any purpose. 
It is provided "as is" without expressed or implied warranty.
*/
package au.csiro.JamesGunning.lib.Numerical;

import cern.colt.matrix.*;
/**
   <p>
  Extension of cern.colt.matrix.impl.DenseDoubleMatrix3D to allow some faster matrix access
  <p>
@author James.Gunning@csiro.au
@version 1.0, 09/24/99
*/
public class DenseDoubleMatrix3D extends cern.colt.matrix.impl.DenseDoubleMatrix3D implements MatrixInterface3D 
{
  // variables useful in fast 3D interpolation below
  private double A_neg1=0,A_0=1,A_1=0,A_2=0,wt_last=0;

  public DenseDoubleMatrix3D(double[][][] values)
  {
	super(values);
  }
  public DenseDoubleMatrix3D(int slices,int rows, int columns) 
  {
	super(slices,rows,columns);
  }

  /**
	 Returns a new deep copy of the 3D matrix
  */
  public DenseDoubleMatrix3D(DenseDoubleMatrix3D src) 
  {
	super(src.slices(),src.rows(),src.columns());
	assign(src);
  }

  public DoubleMatrix2D viewSlice(int slice)
  {
	return super.viewSlice(slice);
  }

  public double[] viewElements()
  {
	return elements;
  }

  public int getElementIndex(int slice, int row, int column)
  {
	return (sliceZero + slice*sliceStride + rowZero + row*rowStride + columnZero + column*columnStride);
  }

  public int getColumnStride()
  {
	return columnStride;
  }
  /**
	 4 pt Lagrange interpolation along a stride, at column i, row j. No bounds checking, so
	 be sure to call with safe i,j,k values and 9<=wt<=1. Assumes no transposes etc have
	 taken place, so column stride=1. Check (externally) that all this is true before
	 using this method!
  */
  public double interpolateQuick(int i,int j,int k,double wt)
  {
	int indx;
	
	// Compute the four coefficients A_{-1},A_0,A_1,A_2
	// as per p 879 of Abramowitz & Stegun....
	// NB 1/6 = 0.16666666666666666666;
	if (wt != wt_last)
	  {
		A_neg1 = -wt*(wt-1.0)*((A_1 = wt-2.0))*0.1666666666666666; // -wt*(wt-1.0)*((tmp1 = wt-2.0))*0.1666666666666666;
		A_0 = ((A_2 = wt*wt-1.0))*A_1*0.5; // A_0 = ((tmp2 = wt*wt-1.0))*tmp1*0.5;
		A_1 *= -wt*(wt+1.0)*0.5; // -wt*(wt+1.0)*tmp1*0.5;
		A_2 *= wt*0.1666666666666666; // A_2 = wt*tmp2*0.1666666666666666;
		wt_last = wt;
	  }
	indx = sliceZero + i*sliceStride + rowZero + j*rowStride + columnZero + k; // k*columnStride last term ; take columnStride==1;

	
	return A_neg1 * elements[indx-1] + A_0* elements[indx] + A_1 * elements[indx+1] + A_2* elements[indx+2];
  }
  /**
	 3D interpolation, formed by bilinear interpolation in the i,j indices, with
	 bilinear weights xWt, yWt, and 4 pt Lagrange interpolation along the k direction, with
	 weight wt. In k direction, interpolation is mixture of e[k-1],e[k],e[k+1],e[k+2], with
	 wt the fractional distance between e[k] and e[k+1].

	 No bounds checking, so
	 be sure to call with safe i,j,k values and 9<=wt<=1. Assumes no transposes etc have
	 taken place, so column stride=1. Check (externally) that all this is true before
	 using this method!
  */
  public double interpolate3DQuick(int i,int j,int k,double wt,double xWt,double yWt)
  {
	int indx;
	double a00,a01,a11,a10;

	// Compute the four coefficients A_{-1},A_0,A_1,A_2
	// as per p 879 of Abramowitz & Stegun....
	// NB 1/6 = 0.166666666666666;
	if (wt != wt_last)
	  {
		A_neg1 = -wt*(wt-1.0)*((A_1 = wt-2.0))*0.166666666666666; // -wt*(wt-1.0)*((tmp1 = wt-2.0))*0.166666666666666;
		A_0 = ((A_2 = wt*wt-1.0))*A_1*0.5; // A_0 = ((tmp2 = wt*wt-1.0))*tmp1*0.5;
		A_1 *= -wt*(wt+1.0)*0.5; // -wt*(wt+1.0)*tmp1*0.5;
		A_2 *= wt*0.166666666666666; // A_2 = wt*tmp2*0.166666666666666;
		wt_last = wt;
	  }
	// correct formula is 
	// sliceZero + i*sliceStride + rowZero + j*rowStride + columnZero + k*columnStride 
	// but take columnStride=1 (usual case)
	indx = sliceZero + i*sliceStride + rowZero + j*rowStride + columnZero + k;
	a00 = A_neg1 * elements[indx-1] + A_0* elements[indx] + A_1 * elements[indx+1] + A_2* elements[indx+2];

	// indx = sliceZero + (i+1)*sliceStride + rowZero + j*rowStride + columnZero + k;
	indx += sliceStride;
	a10 = A_neg1 * elements[indx-1] + A_0* elements[indx] + A_1 * elements[indx+1] + A_2* elements[indx+2];

	// indx = sliceZero + (i+1)*sliceStride + rowZero + (j+1)*rowStride + columnZero + k;
	indx += rowStride;
	a11 = A_neg1 * elements[indx-1] + A_0* elements[indx] + A_1 * elements[indx+1] + A_2* elements[indx+2];

	// indx = sliceZero + i*sliceStride + rowZero + (j+1)*rowStride + columnZero + k;
	indx -= sliceStride;
	a01 = A_neg1 * elements[indx-1] + A_0* elements[indx] + A_1 * elements[indx+1] + A_2* elements[indx+2];

	// Bilinear formula for result
	return (1-xWt)*(1-yWt)*a00 + xWt*(1-yWt)*a10 + xWt*yWt*a11 + (1-xWt)*yWt*a01;
  } 
}
