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

import java.util.*;

/**
   Class of utility functions to produce sampled points in the unit hypercube of
   arbitrary dimensions and sampling density 
*/

public final class Partitions
{
  // for partitions
  ArrayList combolist;
  int Np,Ne;
  boolean allowDuplicates;
  int[] pos;

  private Partitions(int Ne,int Np,boolean allowDuplicates)
  {
	this.Ne = Ne;
	this.Np = Np;
	this.allowDuplicates = allowDuplicates;
	pos = new int[Np];
  }



 public static ArrayList createPartitionList(int Ne, int Np, boolean allowDuplicates)
  {
	int j_p;

	if (Np>Ne) return null;

	Partitions par = new Partitions(Ne,Np,allowDuplicates);

	par.combolist = new ArrayList();

	par.comboLoop(0);
	par.combolist.trimToSize();

	return par.combolist;
  }
  
  /**
	 Recursive function to assess combinations
  */
  private void comboLoop(int j_p)
  {
	int i;
	
	if (j_p < Np)
	  {
		// continue to form nested loops
		for (i = ((j_p==0) ? 0 : pos[j_p-1]+(allowDuplicates?0:1)); i < Ne; i++)
		  {
			pos[j_p] = i;
			comboLoop(j_p+1);
		  }
	  }
	else
	  {
		// add to list
		combolist.add((int[])pos.clone());
	  }
  }

  /** 
	Test driver: show partitions of Ne pts into Np points
  */
  public static void main(String argv[])
  {
	int i,j,Ne=1,Np=1,allow=0;
	if (argv.length!=3)
	  {
		System.out.println("Usage: java au.csiro.JamesGunning.lib.Numerical.Partitions Ne Np [1|0]");
		System.exit(1);
	  }
	else
	  {
		Ne = Integer.parseInt(argv[0]);
		Np = Integer.parseInt(argv[1]);
		allow = Integer.parseInt(argv[2]);
	  }

	ArrayList list = createPartitionList(Ne,Np,(allow==0)?false:true);
	int[] pos;

	System.out.println("Number of combinations = "+list.size()+ 
					   " (nCp = "+cern.jet.math.Arithmetic.binomial(Ne, Np)+", n^p/p!="
					   +Math.pow(Ne,Np)/cern.jet.math.Arithmetic.factorial(Np)+")");
	if (true)
	  {
		for (i=0;i<list.size();i++)
		  {
			pos = (int[]) list.get(i);
			for (j=0;j<pos.length;j++) System.out.print(pos[j]+" ");
			System.out.println("");
		  }
	  }
  }
	  
}
