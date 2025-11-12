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
import cern.colt.list.DoubleArrayList;

/** 
  Extension of cern.colt.matrix.impl.DenseDoubleMatrix2D to allow some snappy
  matrix multiplication calls of lower triangles x vectors (a DGEMV) 
*/
public class DenseDoubleMatrix2D extends cern.colt.matrix.impl.DenseDoubleMatrix2D implements MatrixInterface2D
{
  public DenseDoubleMatrix2D(double[][] values)
  {
	super(values);
  }
  public DenseDoubleMatrix2D(int rows, int columns) 
  {
	super(rows,columns);
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
	return new DenseDoubleMatrix2D(rows, columns);
  }

  /** 
	  Risky exposure of elements to external methods
  */
  public double[] getElements()
  {
	return elements;
  }

  /**
	 Identity
  */
  public static DenseDoubleMatrix2D identity(int N) 
  {
    DenseDoubleMatrix2D I = new DenseDoubleMatrix2D(N,N);
	int i;
	for(i=0;i<N;i++) I.setQuick(i,i,1.0);

	return I;
  }


  /** 
	  Fast lower triangular multiplication of this with rhs, 
	  with answer placed into ans, and
	  multiplication only over rows given (both inclusive). 
	  No protection offered on ranges, for speed reasons
  */

  public void lowerTriangularMult(double[] ans, cern.colt.matrix.DoubleMatrix1D rhs, int rowStart, int rowFinish)
  {
	int i,j,k,offset,r_offset;
	double s=0;
	
	if (rhs instanceof au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D && 
		(((au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D)rhs).getStride()==1))
	  {
		// fast method for dense lower triangular DGEMV
		// usual case I believe, unless *this* is a transpose view (columnStride!=1)
		au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D RHS = 
		  (au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D)rhs;
		double[] r_elems = RHS.getElements();
		r_offset = RHS.getIndex0();

		int jj;
		if (columnStride==1)
		  {
			for (i=rowStart;i<=rowFinish;i++)
			  {
				s=0;
				offset = rowZero + i*rowStride + columnZero;
				// basic
				// for (j=0,jj=offset;j<=i;j++,jj++) s += elements[jj] * RHS.getQuick(j);
				// unrolled
				k = r_offset;
				for (j=0,jj=offset;j<=i%4;j++) s += elements[jj++] * r_elems[k++];
				for (;j<=i;j+=4) 
				  {
					s += elements[jj++] * r_elems[k++]
					  +elements[jj++] * r_elems[k++]
					  +elements[jj++] * r_elems[k++]
					  +elements[jj++] * r_elems[k++];
				  }
				ans[i]=s;
			  }
		  }
		else
		  {
			// *this* is probably a transpose view
			for (i=rowStart;i<=rowFinish;i++)
			  {
				s=0;
				offset = rowZero + i*rowStride + columnZero;
				k = r_offset;
				jj=offset;
				for (j=0;j<=i;j++)
				  {
					//s += elements[offset + j*columnStride] * r_elems[k++];
					s += elements[jj] * r_elems[k++];
					jj += columnStride;
				  }
				ans[i]=s;
			  }
		  }
	  }
	else
	  {
		for (i=rowStart;i<=rowFinish;i++)
		  {
			s=0;
			offset = rowZero + i*rowStride + columnZero;
			for (j=0;j<=i;j++)
			  {
				s += elements[offset + j*columnStride] * rhs.getQuick(j);	
			  }
			ans[i]=s;
		  }
	  }
  }

    /**
     * Construct and returns a new 1-d matrix <i>of the corresponding dynamic type</i>, entirelly independent of the receiver.
     * For example, if the receiver is an instance of type <tt>DenseDoubleMatrix2D</tt> the new matrix must be of type <tt>DenseDoubleMatrix1D</tt>,
     * if the receiver is an instance of type <tt>SparseDoubleMatrix2D</tt> the new matrix must be of type <tt>SparseDoubleMatrix1D</tt>, etc.
     *
     * @param size the number of cells the matrix shall have.
     * @return  a new matrix of the corresponding dynamic type.
     */
    public DoubleMatrix1D like1D(int size) {
	    return new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D(size);
    }
    /**
     * Construct and returns a new 1-d matrix <i>of the corresponding dynamic type</i>, sharing the same cells.
     * For example, if the receiver is an instance of type <tt>DenseDoubleMatrix2D</tt> the new matrix must be of type <tt>DenseDoubleMatrix1D</tt>,
     * if the receiver is an instance of type <tt>SparseDoubleMatrix2D</tt> the new matrix must be of type <tt>SparseDoubleMatrix1D</tt>, etc.
     *
     * @param size the number of cells the matrix shall have.
     * @param zero the index of the first element.
	 * @param stride the number of indexes between any two elements, i.e. <tt>index(i+1)-index(i)</tt>.
	 * @return  a new matrix of the corresponding dynamic type.
	 */
  protected DoubleMatrix1D like1D(int size, int zero, int stride) 
  {
	return new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D(size,this.elements,zero,stride);
  }






  /**
	 Basic timing tests
  */
  public static void main(String[] args)
  {
	int i,j,k,l,trans,N=0,n=1000;
	int rowStart=0,rowEnd=1;
	long st,dt,opC=0;
	double s;
	cern.jet.random.engine.MersenneTwister mT = new cern.jet.random.engine.MersenneTwister();

	DenseDoubleMatrix2D L;


	for (trans=0;trans<=1;trans++)
	  {
		System.out.println("Naive Mflop rate:    CERN-full-block-unrolled flop rate:   Custom tri-mult flop rate (trans="+trans+"): (custom-max error)");
		for (l=40;l<=600;l+=80)
		  {
			N = l;

			rowEnd=N-1;

			double[] ans = new double[N];
			double[] ans2 = new double[N];
			// these two are from au.csiro.JamesGunning.lib.Numerical.*
			DenseDoubleMatrix2D A = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(N,N);
			DenseDoubleMatrix1D r = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D(N);

			L=A;
			if (trans==1) 
			  {
				L = (au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D)A.viewDice();
			  }
			// Randomly fill rhs, L
			L.assign(mT);
			r.assign(mT);

			// naive
			st = System.currentTimeMillis();
			dt=0;
			for (k=0;dt<1000;k++)
			  {
				for (i=rowStart,opC=0;i<=rowEnd;i++)
				  {
					s=0;
					for (j=0;j<=i;j++,opC++)
					  {
						s += L.getQuick(i,j) * r.getQuick(j);
					  }
					ans[i]=s;
				  }
				dt = System.currentTimeMillis()-st;
			  }

			System.out.print("" + (int)(1.0e-3*opC*k/(double)dt));

			// unrolled but double the work (full matrix mult)
			st = System.currentTimeMillis();
			dt=0;
			for (k=0;dt<1000;k++) 
			  {
				L.zMult(r,null,1.0,0,false);				
				dt = System.currentTimeMillis()-st;
			  }
			System.out.print("          " + (int)(1.0e-3*k*2*opC/(double)dt));

			// custom
			st = System.currentTimeMillis();
			dt=0;
		    for (k=0;dt<1000;k++)  
			  {
				L.lowerTriangularMult(ans2,r,rowStart,rowEnd);
				dt = System.currentTimeMillis()-st;
			  }
			System.out.print("         " + (int)(1.0e-3*k*opC/(double)dt));
			for (i=0;i<N;i++) ans2[i]=Math.abs(ans2[i]-ans[i]);	
			System.out.println("         "+cern.jet.stat.Descriptive.max(new DoubleArrayList(ans2))); 
		  }
	  }
  }
}



