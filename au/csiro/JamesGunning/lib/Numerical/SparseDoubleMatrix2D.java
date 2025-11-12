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
  Extension of cern.colt.matrix.impl.SparseDoubleMatrix2D to allow some snappy
  matrix multiplication calls of lower triangles x vectors (a DGEMV) 
*/
public class SparseDoubleMatrix2D extends cern.colt.matrix.impl.SparseDoubleMatrix2D  implements MatrixInterface2D
{
  public SparseDoubleMatrix2D(double[][] values)
  {
	super(values);
  }
  public SparseDoubleMatrix2D(int rows, int columns) 
  {
	super(rows,columns);
  }
  public SparseDoubleMatrix2D(int rows, int columns, int initialCapacity, double minLoadFactor, double maxLoadFactor)
  {
	super(rows,columns,initialCapacity,minLoadFactor,maxLoadFactor);
  }
  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i> as the receiver, having the specified number of rows and columns.
   * For example, if the receiver is an instance of type <tt>DenseDoubleMatrix2D</tt> the new matrix must also be of type <tt>DenseDoubleMatrix2D</tt>,
   * if the receiver is an instance of type <tt>SparseDoubleMatrix2D</tt> the new matrix must also be of type <tt>SparseDoubleMatrix2D</tt>, etc.
   * In general, the new matrix should have internal parametrization as similar as possible.
   *
   * @param rows the number of rows the matrix shall have.
   * @param columns the number of columns the matrix shall have.
   * @return  a new empty matrix of the same dynamic type.
   */
  public cern.colt.matrix.DoubleMatrix2D like(int rows, int columns) {
	return new SparseDoubleMatrix2D(rows, columns);
  }


  /** 
	  Sparse lower triangular multiplication of this with rhs, 
	  with answer placed into ans, and
	  multiplication only over rows given (both inclusive). 
	  No protection offered on ranges, for speed reasons
  */
  public void lowerTriangularMult(double[] ans, cern.colt.matrix.DoubleMatrix1D rhs, int rowStart, int rowFinish)
  {
	int i,j,k,offset,r_offset;
	double s=0;

	for (i=rowStart;i<=rowFinish;i++)
	  {
		s=0;
		for (j=0;j<=i;j++)
		  {
			s += getQuick(i,j) * rhs.getQuick(j);
		  }
		ans[i]=s;
	  }
  }
  public DoubleMatrix1D like1D(int size) 
  {
	return new cern.colt.matrix.impl.SparseDoubleMatrix1D(size);
  }
  /*
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic type</i>, sharing the same cells.
   * For example, if the receiver is an instance of type <tt>DenseDoubleMatrix2D</tt> the new matrix must be of type <tt>DenseDoubleMatrix1D</tt>,
   * if the receiver is an instance of type <tt>SparseDoubleMatrix2D</tt> the new matrix must be of type <tt>SparseDoubleMatrix1D</tt>, etc.
   *
   * @param size the number of cells the matrix shall have.
   * @param offset the index of the first element.
   * @param stride the number of indexes between any two elements, i.e. <tt>index(i+1)-index(i)</tt>.
   * @return  a new matrix of the corresponding dynamic type.
   */
  /*
	For now, superclasses's SparseDoubleMatrix1D is ok here

	protected DoubleMatrix1D like1D(int size, int offset, int stride) 
	{
	return new SparseDoubleMatrix1D(size,this.elements,offset,stride);
	}
  */
}
