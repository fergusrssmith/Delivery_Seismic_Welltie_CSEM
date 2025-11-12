/**********************************************************************
Copyright (C) 2010, CSIRO,  all rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of CSIRO nor the names of its
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

package au.csiro.JamesGunning.lib.cgeom;

import java.io.*;

import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.util.*; 


/**
   Data structures useful for Triangulation

   @author James.Gunning@csiro.au
   @version 1.0, 30/4/2010
*/

public class TriangleUtil
{
  public int Nt; // number of triangles
  public int Na; // number of attributes
  public int[][] nodes=null;   // k = nodes[i][j], numbered k=1,2...Nt. nodes[i=0...Nt-1][0,1,2] = nodes numbers of triangle i's vertices
  public double[][] xys=null;  // xys[k-1][0,1] = x,y of vertex k
  public double values[]=null; // atttribute 1 of triangle, values[0...Nt-1]

  /**
	 Return a list RLL[0...(xys.length-1)] = getReverseLookupList() such that

	 L[i] = list of triangle numbers that use node i (i=zero-offset). Triangle numbers are zero offset - for use directly in
	 the array nodes. Typically:

	 RLL = getReverseLookupList();
	 IntArrayList L = RLL[i];
	 for(j=0;j<L.size();j++)
	     triangle_number = L.get(j);
		 node1 = nodes[triangle_number][0];
		 node2 = nodes[triangle_number][1];
		 node3 = nodes[triangle_number][2];
		 do stuff with triangle. It should have a corner labelled i....
  */
  public cern.colt.list.IntArrayList[] getReverseLookupList()
  {
	int i,j;
	int[] npi;

	// make reverse lookup list
	int Nn = xys.length;
	cern.colt.list.IntArrayList[] RLL = new cern.colt.list.IntArrayList[Nn];
	for (i=0;i<Nn;i++) RLL[i] = new cern.colt.list.IntArrayList();
	// loop over triangles
	for (i=0;i<Nt;i++)
	  {
		npi = nodes[i];
		//unit--offset node numbers are npi[0,1,2]
		for(j=0;j<3;j++) RLL[npi[j]-1].add(i);         // RLL[i] now has triangle numbers belonging to node i
	  }
	for (i=0;i<Nn;i++) RLL[i].trimToSize(); 
	return RLL;
  }

  /**
	 Write coords of triangle i's center into xy[0,1]
  */
  public final void getTriangleCenter(int i,double[] xy)
  {
	int[] npi;
	double[] _p0=null,_p1=null,_p2=null;

	npi = nodes[i];
	_p0 = xys[npi[0]-1];
	_p1 = xys[npi[1]-1];
	_p2 = xys[npi[2]-1];
	xy[0] = (_p0[0]+_p1[0]+_p2[0])/3.0;
	xy[1] = (_p0[1]+_p1[1]+_p2[1])/3.0;
  }


  public boolean readFromEleNodeFilenamePair(String filename, int verbosity, ErrorHandler EH)
  {
	int i,j,k,Nmax=0;
	String s,fname = filename;
	BufferedReader in=null;
	String[] sa;
	boolean gotHeader;
	// get elements from filename.ele file
	try
	  {
		fname = Util.fileSansSuffix(fname)+".ele";
		in = new BufferedReader(new FileReader(fname));
	  }
	catch (Exception e) 
	  {
		if (verbosity>0) EH.warn("Can't open file "+fname);
		return false;
	  }
	try
	  {
		i=j=0;
		gotHeader = false;
		while ((s=in.readLine()) != null)
		  {
			j++;
			if (s.startsWith("#")) continue;
			sa = s.trim().split("\\s+");     // split just at white space
			if (!gotHeader && sa.length>=3)
			  {
				try
				  {
					Nt = Integer.parseInt(sa[0]);
					Na = Integer.parseInt(sa[2]);
					nodes=new int[Nt][3];
					values = new double[Nt];
					gotHeader=true;
				  }
				catch (NumberFormatException e)
				  {
					if (verbosity>0) EH.warn("Line "+j+" of file " + fname+" has non-numbers");
					return false;
				  }			
			  }
			else if (gotHeader && (sa.length>4) && (i<Nt))
			  {
				try
				  {
					// columns 2,3,4 are node labels for this ith triangle. column 5 is the value
					for(k=0;k<3;k++) nodes[i][k] = Integer.parseInt(sa[k+1]);
					values[i] = Double.parseDouble(sa[4]);
				  }
				catch (NumberFormatException e)
				  {
					if (verbosity>0) EH.warn("Line "+j+" of file " + fname+" has non-numbers");
					return false;
				  }	
				i++;
			  }	
		  }
		in.close();
		if (!gotHeader) 
		  {
			if (verbosity>0) EH.warn("Couldn't parse "+fname);
			return false;
		  }
		if (verbosity>0) EH.warn("read "+fname);
	  }
	catch (Exception e) 
	  {
		if (verbosity>0) EH.warn("Can't open, decipher or close file "+fname);
		return false;
	  }

	// get x,y's from filename.node file
	try
	  {
		fname = Util.fileSansSuffix(fname)+".node";
		in = new BufferedReader(new FileReader(fname));	  
	  }
	catch (Exception e) 
	  {
		if (verbosity>0) EH.warn("Can't open file "+fname);
		return false;
	  }
	try
	  {
		i=j=0;
	    gotHeader = false;
		while ((s=in.readLine()) != null)
		  {
			j++;
			if (s.startsWith("#")) continue;
			sa = s.trim().split("\\s+");     // split just at white space
			if (!gotHeader && sa.length>=3)
			  {
				try
				  {
					Nmax = Integer.parseInt(sa[0]);
					gotHeader=true;
					xys = new double[Nmax][2];
				  }
				catch (NumberFormatException e)
				  {
					if (verbosity>0) EH.die_err("Line "+j+" of file " + fname+" has non-numbers");
				  }			
			  }
			else if (gotHeader && (sa.length>=4))
			  {
				if (i>=Nmax)
				  {
					if(verbosity>0) EH.warn("Too many lines in file " + fname);
					return false;
				  }


				try
				  {
					// columns 2,3 have x,y values of node i
					for(k=0;k<2;k++) xys[i][k] = Double.parseDouble(sa[k+1]); 
				  }
				catch (NumberFormatException e)
				  {
					if (verbosity>0) 
					  {
						EH.die_err("Line "+j+" of file " + fname+" has non-numbers");
					  }
					return false;
				  }	
				i++;
			  }	
		  }
		in.close();
		if (verbosity>0) EH.warn("read "+fname);
	  }
	catch (Exception e) 
	  {
		if (verbosity>0) EH.warn("Can't open, decipher or close file "+fname);
		return false;
	  }
	return true;
  }	
} 
