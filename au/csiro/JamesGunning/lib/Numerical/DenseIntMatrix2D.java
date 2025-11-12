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
   Simple class for storing contiguous 2D int arrays (COLT-like) for communications with
   fortran
*/
public class DenseIntMatrix2D
{
  int rows;
  int columns;
  /**
	 make public so can be passed to external routines as
	 an unrolled array or pointer
  */
  public int[] elements;

  public DenseIntMatrix2D(int rows, int columns)
  {
	this.rows = rows;
	this.columns = columns;

	elements = new int[rows*columns];
  }


  public DenseIntMatrix2D(DenseIntMatrix2D m)
  {
	this.rows = m.rows;
	this.columns = m.columns;

	elements = new int[rows*columns];
  }

  public DenseIntMatrix2D assign(DenseIntMatrix2D m)
  {
	if ((rows == m.rows) && (columns == m.columns))
	  {
		System.arraycopy(m.elements,0,elements,0,rows*columns);
	  }
	else
	  {
		throw new IllegalArgumentException("DenseIntMatrix2D.assign(): arrays of incommensurate size.");
	  }
	return this;
  }

  public void setQuick(int row, int column, int value) 
  {
	// manually inlined:
	elements[row*columns + column] = value;
  }

  public int getQuick(int row, int column) 
  {
	// manually inlined:
	return elements[row*columns + column]; 
  }
}

