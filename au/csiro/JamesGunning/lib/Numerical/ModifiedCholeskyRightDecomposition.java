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
import cern.colt.matrix.*;

/** 
From a symmetric, positive definite matrix <tt>A</tt>,
constructs a 'safely' positive definite Cholesky Right decomposition
using an algorithm very similar to that of Algorithm 6.5, CH.6,
"Numerical Optimization", Nocedal and Wright (Springer). Guaranteed to rescue
you if the matrix is ill-conditioned or not quite positive definite.
*/
public class ModifiedCholeskyRightDecomposition 
  //extends cern.colt.matrix.linalg.CholeskyDecomposition
{
  private DoubleMatrix2D R;
  private DoubleMatrix2D _A;  // reference to calling A

  private int n;

  private boolean isSymmetricPositiveDefinite;
  private boolean useBlockDiagMethod = false;
  private boolean isDiag = false;
  private cern.colt.list.IntArrayList blockStarts;

  public static int MAX_BLOCK_SIZE_FOR_BLOCK_DIAG_METHOD = 20;

  // variables for debugging....
  private static final boolean debug = false;
  private static final double DEBUG_MAX_ERROR = 1e-10;

 /** 
	  Example test harness
  */
  public static void main(String argv[])
  {
	int i,j,n;
	double lambda;
	ModifiedCholeskyRightDecomposition chol;
	DoubleMatrix2D A,R,C;

	cern.colt.matrix.linalg.Blas blas = cern.colt.matrix.linalg.SeqBlas.seqBlas;
	// diag
	n=10;
	A = new DenseDoubleMatrix2D(n,n);
	C = A.like();
	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		A.setQuick(i,j,(i==j)?1.0:0);
	  }
	}
	R = (new ModifiedCholeskyRightDecomposition(A)).getR();
	blas.dgemm(false,true,1,R,R,-1,A);
	System.out.println("A-R.R^T: where A=diag(1)"+A.toString());

	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		A.setQuick(i,j,(i==j)?100:0);
	  }
	}
	R = (new ModifiedCholeskyRightDecomposition(A)).getR();
	blas.dgemm(false,true,1,R,R,-1,A);
	System.out.println("A-R.R^T: where A=diag(100)"+A.toString());


	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		A.setQuick(i,j,(i==j)?1e-4:0);
	  }
	}
	R = (new ModifiedCholeskyRightDecomposition(A)).getR();
	blas.dgemm(false,true,1,R,R,-1,A);
	System.out.println("A-R.R^T: where A=diag(1e-4)"+A.toString());



	// Guassian covariance is badly conditioned....
	lambda=10000;
	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		A.setQuick(i,j,1e-6*Math.exp(-(i-j)*(i-j)/lambda*lambda));
	  }
	}
	R = (new ModifiedCholeskyRightDecomposition(A)).getR();
	blas.dgemm(false,true,1,R,R,-1,A);
	System.out.println("A-R.R^T: where A=Gaussian-covariance(1e-6,range=10000)"+A.toString());
	double max=0;
	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		max=Math.max(max,Math.abs(A.getQuick(i,j)));
	  }
	}
	System.out.println("max error = "+max);
	
	}




  /**
	Constructs a secure Cholesky Right decomposition of the square symmetric matrix A at 
	initialisation. I.e. the symmetric positive definite A is decomposed as A = R.R^T, where
	R is right upper triangular. Does not corrupt A.
	
	@param  A   Square matrix: Positive definite or close thereto.
	@throws IllegalArgumentException if <tt>A</tt> is not square.
  */

  public ModifiedCholeskyRightDecomposition(DoubleMatrix2D A) 
  {
	computeModifiedCholeskyRightDecomposition(A);
  }

  
  /**
	 Same as main constructor, but can test for block diagonal matrices in O(n) time 
	 and speed things up accordingly. Natural ordering is presumed.
	 
	@param  A   Square matrix: Positive definite or close thereto.
	@param  testForDiagonality: If true, test if matrix is diagonal, and apply speedups if so
	@throws IllegalArgumentException if <tt>A</tt> is not square.
  */

  public ModifiedCholeskyRightDecomposition(DoubleMatrix2D A, boolean testForDiagonality) 
  {
	if (!testForDiagonality) 
	  {
		computeModifiedCholeskyRightDecomposition(A);
	  }
	else
	  {
		// try to be smart
		n = A.rows();
		R = A.like(n,n);		
		_A = A;
		isSymmetricPositiveDefinite = (A.columns() == n);
        
		BlockDiagTest test = new BlockDiagTest(A);
		isDiag = test.isDiag();
		if (!isDiag)
		  {
			int maxBlockSize;
			maxBlockSize = test.getMaxBlockSize();
			if (maxBlockSize < MAX_BLOCK_SIZE_FOR_BLOCK_DIAG_METHOD)
			  {
				// System.out.println("Max block size = "+maxBlockSize);
				useBlockDiagMethod = true;
				// do blockwise Cholesky over all blocks
				int i,j,iS,iE,width,nB;
				nB = test.blockStarts.size();
				blockStarts = test.blockStarts;
				for (i=0;i<nB;i++)
				  {
					iS = blockStarts.getQuick(i);
					iE = (i<(nB-1)) ? blockStarts.getQuick(i+1)-1 : n-1;
					width = iE-iS+1;
					DoubleMatrix2D Asub = A.viewPart(iS,iS,width,width);
					DoubleMatrix2D Rsub = R.viewPart(iS,iS,width,width);
					ModifiedCholeskyRightDecomposition cholSub = new ModifiedCholeskyRightDecomposition(Asub);
					// write R factor of block into overall block
					Rsub.assign(cholSub.getR()); 
				  }
			  }
			else
			  {
				// isn't simple. Do the whole catastrophe
				computeModifiedCholeskyRightDecomposition(A);
			  }
		  }
		else
		  {
			// simple diagonal case. Trivial.
			for (int i=0;i<n;i++)
			  {
				R.setQuick(i,i,Math.sqrt(Math.max(A.getQuick(i,i),0)));
			  }
		  }
	  }
  }

  private void computeModifiedCholeskyRightDecomposition(DoubleMatrix2D A)
  {
	int i,j,k;
	n = A.rows();
	R = A.like(n,n);

	_A = A;

	double max=0,min = Double.MAX_VALUE;  // DEBUG -> CAN REMOVE
	boolean haveRepivoted = false;        // DEBUG -> CAN REMOVE

	isSymmetricPositiveDefinite = (A.columns() == n);
        
	// precompute and cache some views to avoid regenerating them time and again
	DoubleMatrix1D[] Rrows = new DoubleMatrix1D[n];
	for (j = 0; j < n; j++) Rrows[j] = R.viewRow(j);
        
	//parameter initialisation for modified Cholesky
	double gamma = A.getQuick(0,0);
	for(i = 1; i < n; i++)
	  {
		gamma = Math.max(gamma,A.getQuick(i,i));
	  }
	double xi = 0.0;
	for(i = 0; i < n; i++)
	  {
		for(j = i+1; j < n; j++)  //xi from linear search of the upper right triangular
		  {
			xi = Math.max(xi,A.getQuick(i,j));
		  }
	  }
	// get machine epsilon EPSILON_SMALL from the com.imsl.math.Sfun class in the colt distro. (1.11 x 10^-16)
	double delta = MachineConstants.EPSILON_SMALL * Math.max(gamma+xi,1.0);
	double beta = Math.max(gamma, Math.max(xi/Math.sqrt(n*n - ((n>1)?1:0)), MachineConstants.EPSILON_SMALL));
		

	// Main loop.
	for (i = n-1; i >=0 ; i--) 
	  {
		double theta = Double.NEGATIVE_INFINITY;
		for(k = i; k >= 0; k--)
	      {
			theta = Math.max(theta,A.getQuick(i,k) - Rrows[i].zDotProduct(Rrows[k],i+1,n));
	      }
		// JG HACK...
		// theta stuff is wrong. This will work until debugged version is available....
		theta=0.0;
		// .... JGHACK
		for (j = i; j>=0; j--) 
		  {
			double s = Rrows[i].zDotProduct(Rrows[j],i+1,n);
			s = (A.getQuick(i,j) - s);
		  
			if (i==j)
			  {
				isSymmetricPositiveDefinite = isSymmetricPositiveDefinite && (s > 0.0);
				R.setQuick(i,i, Math.sqrt(Math.max(delta,Math.max(Math.abs(s),(theta*theta)/(beta*beta)))));
				if (debug) 
				  {
					if (delta>s || (theta*theta/beta/beta)>s) haveRepivoted = true;
					min = Math.min(min,R.getQuick(i,i)); 
					max=Math.max(max,R.getQuick(i,i));  // DEBUG -> CAN REMOVE
				  }
			  }
			else
			  {
				R.setQuick(j,i,s/R.getQuick(i,i));
			  }
		  }                
		for (k = 0; k < i; k++) 
		  {
			R.setQuick(i,k, 0.0);
		  }
	  }
	if (debug)
	  {
		if (haveRepivoted) System.out.println("CHOLESKY REPIVOTED");
		debugTestCholeskyMultiplication();
	  }
  }

  // Test whether derived Cholesky factor satisfies A = R.R^T to prescribed accuracy
  private void debugTestCholeskyMultiplication()
  {
	cern.colt.matrix.linalg.Blas blas = cern.colt.matrix.linalg.SeqBlas.seqBlas;
	double max=0.0, maxNorm=0.0;
	
	DoubleMatrix2D Ans = _A.copy();
	// maximum element of diag(A)
	maxNorm = 0.0;
	for (int i=0;i<Ans.rows();i++) maxNorm = Math.max(maxNorm,Math.abs(Ans.getQuick(i,i)));

	// form max = Max(A-R.R^T)
	blas.dgemm(false,true,1.0,R,R,-1,Ans);
	for (int i=0;i<Ans.rows();i++) {
	 for (int j=0;j<Ans.rows();j++) {
	   max = Math.max(max,Math.abs(Ans.getQuick(i,j)));
	 }
	}
	if (max/maxNorm>DEBUG_MAX_ERROR)
	  {
		System.err.println("CHOLESKY ERROR= "+max/maxNorm);
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


  /**
    Returns the inverse of the right Cholesky matrix R *in place*. Typical calling sequence is
	ModifiedCholeskyRightDecomposition D = ModifiedCholeskyRightDecomposition(A); <br>
	DoubleMatrix2D R = D.getR();<br>
	D.getCholeskyRightInverse(R);<br>
	
	@param R Square upper triangular matrix (usually obtained from getR())
	@throws IllegalArgumentException if <tt>R</tt> is not square.
  */

  public void getCholeskyRightInverse(DoubleMatrix2D R)
  {
	double sum;
	int n = R.rows();
	
	if (R.columns() != n) throw new IllegalArgumentException("matrix is not square");
	
	if (isDiag)
	  {
		// trivial diagonal
		for (int i=0;i<n;i++)
		  {
			sum = R.getQuick(i,i);
			R.setQuick(i,i,(sum!=0)? 1.0/sum : 0);
		  }
	  }
	else
	  {
		if (useBlockDiagMethod)
		  {	
			// do blockwise Cholesky inverse over all blocks
			int i,iS,iE,width,nB;
			nB = blockStarts.size();
			DoubleMatrix2D Rsub;
			useBlockDiagMethod = false; // switch off so we go into dense branch for blocks...
			for (i=0;i<nB;i++)
			  {
				iS = blockStarts.getQuick(i);
				iE = (i<(nB-1)) ? blockStarts.getQuick(i+1)-1 : n-1;
				width = iE-iS+1;
				Rsub = R.viewPart(iS,iS,width,width);
				// do in-place blockwise Cholesky
				getCholeskyRightInverse(Rsub);
			  }
			useBlockDiagMethod = true;
		  }	
		else
		  {
			// Standard method, main loop.
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
  }

  /**
	 Special class for handling block-diagonal matrices in their natural order
  */

  public class BlockDiagTest implements cern.colt.function.IntIntDoubleFunction
  {
	private boolean isDiag=true;
	private int n;
	private int[] jMax;
	public cern.colt.list.IntArrayList blockStarts;

	public BlockDiagTest(DoubleMatrix2D A)
	{
	  n = A.rows();
	  jMax = new int[n];
	  A.forEachNonZero(this);
	}

	/**
	   This method is called only within the instantiator
	   forEachNonZero() call on the sparse matrix.  It detects
	   non-diagonality, and computes some block size information
	*/
	public double apply(int i, int j, double value) 
	{
	  if (i!=j) isDiag=false;
	  if (j>jMax[i]) jMax[i]=j;
	  return value;
	}
	/**
	   Returns whether matrix is diagonal or not
	*/
	public boolean isDiag()
	{
	  return  isDiag;
	}
	/**
	   Call this to compute maximum in-order block size, and set up block indices blockStarts
	   Algorithm decisions based on the returned maximum block size can then be made, 
	   exploiting the blockStart indices in blockStarts.
	*/
	public int getMaxBlockSize()
	{
	  int i,j,iS,iE;
	  int maxBlockSize=1;
	  blockStarts = new cern.colt.list.IntArrayList();
	  for (i=0;i<n;)
		{
		  iS=i;
		  iE=jMax[iS];
		  blockStarts.add(iS);
		  for (j=iS;j<=iE;j++)
			{
			  if (jMax[j]>iE) iE=jMax[j];
			}
		  // System.out.println("Block "+iS+" "+(j-1));		  
		  i=j;
		  maxBlockSize = Math.max(maxBlockSize,(j-iS));
		}
	  blockStarts.trimToSize();
	  jMax=null;

	  return maxBlockSize;
	}
  }
}

