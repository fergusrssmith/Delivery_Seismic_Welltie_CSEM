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
   Class of utility functions to produce sampled points in the unit hypercube of
   arbitrary dimensions and sampling density 
*/

public final class Hypercube
{
  // for Hypercube
  private int dim,Np;
  private int curDim;
  private int curPoint;
  float[] curOffsets;
  float[][] offsets;


  private Hypercube(int dim,int Np)
  {
	this.dim = dim;
	this.Np = Np;
	curOffsets = new float[dim];
	offsets = new float[(int)Math.rint(Math.pow(Np,dim))][dim];
  }

  /* 
	 Generate list of points raster scanning the hypercube of
	 dimension dim, with Np points per side in the unit [0,1]
	 hypercube. The returned array is of size [Np^dim][dim]
  */

  public static float [][] createHypercubeList(int dim, int Np)
  {

	Hypercube hc = new Hypercube(dim,Np);
	  
	return hc.createHyperCube();
  }

  /* 
	 Generate list of points raster scanning the hypercube of
	 dimension dim, with Np_arr[i] points per side for dimension i in
	 the unit [0,1] hypercube. The returned array is of size
	 [\prod_i Nparr[i]][dim]
  */

  public static float [][] createHypercubeList(int Np_arr[])
  {

	Hypercube hc = new Hypercube(Np_arr);
	  
	return hc.createNonUniformHyperCube(Np_arr);
  }

  /**
	 setup hypercube from list of number of samples per 
	 dimensions e.g. Np = {4,4,5...}
  */
  private Hypercube(int Np[])
  {
	int i,size;

	this.dim = Np.length;
	curOffsets = new float[dim];

	for(i=0,size=1; i<this.dim; i++) size *= Np[i];

	offsets = new float[size][dim];
  }

  private float [][] createNonUniformHyperCube(int Np_arr[])
  {
	curDim=-1;

	populate(curOffsets,offsets,Np_arr);

	curOffsets=null;
	return offsets;
  }


  private float [][] createHyperCube()
  {
	curDim=-1;

	populate(curOffsets,offsets);

	curOffsets=null;
	return offsets;
  }


  // recursive method for setting curOffsets 0,1...(Np-1) in a dim-dimensional hypercube
  // the hypercode points are stored in successive rows of offsets. E.g.
  /*  0.0 0.0 0.0 
	  0.0 0.0 0.333 
	  0.0 0.0 0.666 
	  0.0 0.0 1.0 
	  0.0 0.333 0.0 
	  0.0 0.333 0.333 
	  0.0 0.333 0.666 
	  0.0 0.333 1.0 
	  0.0 0.666 0.0 
	  0.0 0.666 0.333 
	  0.0 0.666 0.666 
	  0.0 0.666 1.0 
	  0.0 1.0 0.0 
	  0.0 1.0 0.333 
	  0.0 1.0 0.666 
	  0.0 1.0 1.0 
	  0.333 0.0 0.0 
	  0.333 0.0 0.333 
	  0.333 0.0 0.666 
	  ...
  */

  private void populate(float[] curOffsets, float[][] offsets)
  {
	curDim++;

	for (int i=0;i<Np;i++)
	  {
		curOffsets[curDim] = (Np>1)? (float)i/(float)(Np-1) : 0;
		if (curDim < (dim-1)) 
		  {
			populate(curOffsets,offsets);
		  }
		else
		  {
			// store in array
			System.arraycopy(curOffsets,0,offsets[curPoint],0,dim);
			curPoint++;
			/*			
						System.out.print((curPoint-1)+"    ");
						for (int j=0;j<dim;j++) System.out.print(curOffsets[j]+ " ");
						System.out.println("");
			*/
		  }
	  }
	curDim--;
  } 

  private void populate(float[] curOffsets, float[][] offsets,int Np_arr[])
  {
	int Np;

	curDim++;

	Np = Np_arr[curDim];
	
	for (int i=0;i<Np;i++)
	  {
		curOffsets[curDim] = (Np>1)? (float)i/(float)(Np-1) : 0;
		if (curDim < (dim-1)) 
		  {
			populate(curOffsets,offsets,Np_arr);
		  }
		else
		  {
			// store in array
			System.arraycopy(curOffsets,0,offsets[curPoint],0,dim);
			curPoint++;
			/*			
						System.out.print((curPoint-1)+"    ");
						for (int j=0;j<dim;j++) System.out.print(curOffsets[j]+ " ");
						System.out.println("");
			*/
		  }
	  }
	curDim--;
  } 


  /** 
	  Test driver: make 2 dim sample with 4 points per dim
  */
  public static void main(String argv[])
  {
	int i,j;
	float[][] arr=null;


	if (true)
	  {
		// uniform style
		arr = createHypercubeList(2,2);
	  }
	else
	  {
		// non uniform case: chop 3 dim cube into 3 x 6 x 2 dissection
		int[] N = {3,6,2};
		arr = createHypercubeList(N);
	  }
	for (i=0;i<arr.length;i++)
	  {
		for (j=0;j<arr[0].length;j++) System.out.print(arr[i][j]+" ");
		System.out.println("");
	  }
  }
	  

}
