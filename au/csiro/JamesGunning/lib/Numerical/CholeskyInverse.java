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
Additions to the Cholesky component of the CERN library for inverse from 
Cholesky decompositions.

Original notice.......
Copyright 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose 
is hereby granted without fee, provided that the above copyright notice appear in all copies and 
that both that copyright notice and this permission notice appear in supporting documentation. 
CERN makes no representations about the suitability of this software for any purpose. 
It is provided "as is" without expressed or implied warranty.
*/
package au.csiro.JamesGunning.lib.Numerical;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix1D;
/** 
Assumes: for a symmetric, positive definite matrix <tt>A</tt>, the Cholesky decomposition
is a lower triangular matrix <tt>L</tt> so that <tt>A = L*L'</tt>; this has been
constructed vis the CD = CholeskyDecomposition Class.

The matrix CD = getL() can then be fed to this routine....
*/
public class CholeskyInverse extends cern.colt.matrix.linalg.CholeskyDecomposition
{
 public CholeskyInverse(DoubleMatrix2D A)
  {
	super(A);
  }


  /*
	Constructs and returns the inverse of the triangular factor L
	
	@param  L   Square, lower triangular Cholesky factor.
	@returns    <tt>L^-1</tt> in place in original L factor.
	@throws IllegalArgumentException if <tt>L</tt> is not square.
  */

  public void getCholeskyInverseInPlace(DoubleMatrix2D L) 
  {
	double sum;
	int n = L.rows();
	
	if (L.columns() != n) throw new IllegalArgumentException("matrix is not square");
	
	// Main loop.
	for (int i = 0; i < n; i++) 
	  {
		L.setQuick(i,i,1.0/L.getQuick(i,i));
		for (int j = i+1; j< n; j++) 
		  {
			sum=0.0;
			for (int k = i; k< j ; k++) sum-=L.getQuick(j,k)*L.getQuick(k,i);
			L.setQuick(j,i,sum/L.getQuick(j,j));
		  }
	  }
  }
  /**
	 Return a new matrix which is the inverse of the lower triangular matrix L, leaving L intact
	 @param  L   Square, lower triangular Cholesky factor.
	 @returns    <tt>L^-1</tt>
  */
 public DoubleMatrix2D getCholeskyInverse(DoubleMatrix2D L) 
  {
	double sum;
	int n = L.rows();
	DoubleMatrix2D Linv;

	if (L.columns() != n) throw new IllegalArgumentException("matrix is not square");

	Linv = L.like().assign(L);

	// Main loop.
	for (int i = 0; i < n; i++) 
	  {
		Linv.setQuick(i,i,1.0/Linv.getQuick(i,i));
		for (int j = i+1; j< n; j++) 
		  {
			sum=0.0;
			for (int k = i; k< j ; k++) sum-=Linv.getQuick(j,k)*Linv.getQuick(k,i);
			Linv.setQuick(j,i,sum/Linv.getQuick(j,j));
		  }
	  }
	return Linv;
  }
}

