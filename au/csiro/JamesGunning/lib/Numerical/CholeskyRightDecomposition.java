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
public class CholeskyRightDecomposition 
  //extends cern.colt.matrix.linalg.CholeskyDecomposition
{
  private DoubleMatrix2D R;

  private int n;

  private boolean isSymmetricPositiveDefinite;

  /*
	Constructs and returns the inverse of the triangular factor L
	
	@param  L   Square, lower triangular Cholesky factor.
	@returns    <tt>L^-1</tt> in place in original L factor.
	@throws IllegalArgumentException if <tt>L</tt> is not square.
  */

  public CholeskyRightDecomposition(DoubleMatrix2D A) 
  {
	// Property.DEFAULT.checkSquare(A); FIXME!
        
	n = A.rows();
	R = A.like(n,n);
	isSymmetricPositiveDefinite = (A.columns() == n);
        
	//precompute and cache some views to avoid regenerating them time and again
	DoubleMatrix1D[] Rrows = new DoubleMatrix1D[n];
	for (int j = 0; j < n; j++) Rrows[j] = R.viewRow(j);
        
	// Main loop.
	for (int i = n-1; i >=0 ; i--) 
	  {
	  for (int j = i; j>=0; j--) 
		{
		  double s = Rrows[i].zDotProduct(Rrows[j],i+1,n);
		  s = (A.getQuick(i,j) - s);
		  
		  if (i==j)
			{
			  isSymmetricPositiveDefinite = isSymmetricPositiveDefinite && (s > 0.0);
			  R.setQuick(i,i, Math.sqrt(Math.max(s,0.0)));
			}
		  else
			{
			  R.setQuick(j,i,s/R.getQuick(i,i));
			}
		}                
	  for (int k = 0; k < i; k++) 
		{
		  R.setQuick(i,k, 0.0);
		}
	  }
  }

  /** 
	  Returns the right-triangular factor, <tt>R</tt>.
	  @return     <tt>R</tt>
  */
  public DoubleMatrix2D getR() {
	return R;
  }

  /** 
	  Returns whether the matrix <tt>A</tt> is symmetric and positive definite.
	  @return     true if <tt>A</tt> is symmetric and positive definite; false otherwise
  */
  public boolean isSymmetricPositiveDefinite() {
	return isSymmetricPositiveDefinite;
  }


  public void getCholeskyRightInverse(DoubleMatrix2D R)
  {
	double sum;
	int n = R.rows();
	
	if (R.columns() != n) throw new IllegalArgumentException("matrix is not square");
	
	// Main loop.
	for (int j = n-1; j>=0; j--) 
	  {
		R.setQuick(j,j,1.0/R.getQuick(j,j));
		for (int i = j-1; i >=0; i--) 
		  {
			sum=0.0;
			for (int k = i+1; k <= j ; k++) sum-=R.getQuick(i,k)*R.getQuick(k,j);
			R.setQuick(i,j,sum/R.getQuick(i,i));
		  }
	  }
  }
}
