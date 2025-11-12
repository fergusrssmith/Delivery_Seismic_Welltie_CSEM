/**
   Simple class to enumerate binary splittings of a one-dimensional interval,
   with splittings increasing monotonically to the right. Generates
   integer lists of segment size, e.g. list of integers like

  8 
  4   4 
  4   2   2 
  2   2   2   2 
  2   2   2   1   1 
  2   2   1   1   1   1 
  2   1   1   1   1   1   1 
  1   1   1   1   1   1   1   1 

  which are the binary partitions of the interval [1..8], with any left segment equal to or
  larger than any segment to its right. Must be fed intervals of powers of 2.

*/


package au.csiro.JamesGunning.lib.Numerical;
import java.util.*;
import cern.colt.list.*;

public class IntervalBisectionEnumerator 
{
  int n;
  int n_pts;
  boolean showTree = false;
  boolean sequenceAsSmallToBig = false;
  ArrayList<int[]>points;
  int[] current_map;

  /**
	 Demo example
  */
  public static void main(String[] args) 
  {
	if (args.length<1) 
	  {
		System.out.println("Usage: IntervalBisectionEnumerator log_2(size) (integer)");
		System.exit(0);
	  }
	int N=1;
	try 
	  {
		N =Integer.parseInt(args[0]);
	  }
	catch (NumberFormatException e)
	  {
		System.out.println("Argument "+args[0]+ " not an integer ");
		System.exit(0);
	  }
	if (N<1)
	  {
		System.out.println("Argument "+args[0]+ " < 1");
		System.exit(0);
	  }
	IntervalBisectionEnumerator en;
	en = new IntervalBisectionEnumerator(true,true,false,N);
	List<int[]> pts = en.getPartitions();

  }

  /**
	 Make a new binary partions enumerator for the interval of length 2^N
     Set showTree=true to print result to stdout
	 Segmentation sequence typically like

	 (A) if (noOrderingRestriction = false)
	 .|       .
	 .|   |   .
	 .|   | | .
	 .| | | | .
	 .| | | ||.
	 .| | ||||.
	 .| ||||||.
	 .||||||||.

	 Normal ordering looks like
	 8 
	 4   4 
	 4   2   2 
	 2   2   2   2 
	 2   2   2   1   1 
	 2   2   1   1   1   1 
	 2   1   1   1   1   1   1 
	 1   1   1   1   1   1   1   1

	 Set  sequenceAsSmallToBig to reverse this (numbers on right >= numbers on left)

	 8 
	 4   4 
	 2   2   4 
	 2   2   2   2 
	 1   1   2   2   2 
	 1   1   1   1   2   2 
	 1   1   1   1   1   1   2 
	 1   1   1   1   1   1   1   1
	 
	 
(B) if (noOrderingRestriction = true) we get the full tree:

	 8 
	 4 4 
	 4 2 2 
	 4 2 1 1 
	 4 1 1 2 
	 4 1 1 1 1 
	 2 2 4 
	 2 2 2 2 
	 2 2 2 1 1 
	 2 2 1 1 2 
	 2 2 1 1 1 1 
	 2 1 1 4 
	 2 1 1 2 2 
	 2 1 1 2 1 1 
	 2 1 1 1 1 2 
	 2 1 1 1 1 1 1 
	 1 1 2 4 
	 1 1 2 2 2 
	 1 1 2 2 1 1 
	 1 1 2 1 1 2 
	 1 1 2 1 1 1 1 
	 1 1 1 1 4 
	 1 1 1 1 2 2 
	 1 1 1 1 2 1 1 
	 1 1 1 1 1 1 2 
	 1 1 1 1 1 1 1 1 

	 This explodes enormously with n! The size recurrence is:
	 Size_0 = 1
	 Size_n = 1 + (Size_{n-1})^2
	 So first few sizes are [1,2,5,26,677,458330...]
  */
  public IntervalBisectionEnumerator(boolean showTree,boolean noOrderingRestriction,
									 boolean sequenceAsSmallToBig,int N)
  {
	int i,w;
	double x=1.0;
 
	n = N;
	points = new ArrayList();
	this.showTree = showTree;
	this.sequenceAsSmallToBig = sequenceAsSmallToBig;
	n_pts = (1<<n);
	if (noOrderingRestriction)
	  {
		enumerate_DP(n);
	  }
	else
	  {
		merge_enumerate(n,n_pts); // works
	  }
  }

  /* 
	 Sort of dynamic programming approach to full tree enumeration
  */
  void enumerate_DP(int n)
  {
	int i,j,k,jj,kk,q,sz;
	ArrayList[] L;
	ArrayList<IntArrayList>Lq;
	ArrayList<IntArrayList>Lminus1;
	IntArrayList Lm1j,Lm1k;
	IntArrayList Lqn;
	IntArrayList ls;

	L = new ArrayList[n+1]; 
	for (q=0;q<=n;q++)
	  {
		L[q] = new ArrayList();
		Lq = L[q];
		// full member

		Lqn = new IntArrayList(1);
		Lqn.add(1<<q);
		Lqn.trimToSize();
		Lq.add(Lqn);
		if (q>0)
		  {
			// Splits
			Lminus1 = L[q-1];
			sz = Lminus1.size();
			for (j=0;j<sz;j++)
			  {
				Lm1j = Lminus1.get(j); //I
				for (k=0;k<sz;k++)
				  {
					Lqn = new IntArrayList();
					Lm1k = Lminus1.get(k);
					for (jj=0;jj<Lm1j.size();jj++) Lqn.add(Lm1j.get(jj));
					for (kk=0;kk<Lm1k.size();kk++) Lqn.add(Lm1k.get(kk));
					Lqn.trimToSize();
					Lq.add(Lqn);
				  }
			  }
		  }
		Lq.trimToSize();
	  }
	
	// Assemble into points
	points = new ArrayList<int[]>();
	q = n;
	for(j=0;j<L[q].size();j++)
	  {
		ls = (IntArrayList)(L[q].get(j));
		sz = ls.size()-1;
		int[] chunks = new int[sz+1];
		for(k=0;k<=sz;k++) chunks[k] = ls.get(k);
		points.add(chunks);
	  }
	points.trimToSize();
	if (showTree)
	  {
		for(j=0;j<points.size();j++)
		  {
			int[] chunks = points.get(j);
			for(k=0;k<chunks.length;k++)
			  {
				System.out.print(chunks[k]+" ");					
			  }
			System.out.println("");
		  }
	  }
  }



  /** 
	  Returns arraylist of segment lengths, each sequence garanteed to be a non-increasing
	  integer sequence, e.g.
	  8 
	  4   4 
	  4   2   2 
	  2   2   2   2 
	  2   2   2   1   1 
	  2   2   1   1   1   1 
	  2   1   1   1   1   1   1 
	  1   1   1   1   1   1   1   1 
  */
  public ArrayList<int[]> getPartitions()
  {
	return points;
  }


  private void merge_enumerate(int n, int size)
  {
	int i,j,l,i_len,j_len;

	current_map = new int[size];
	for (i=0;i<size;i++) current_map[i] = 1;
	dumpCurrentList();
	for (l=1;l<=n;l++)
	  {
		for (i=0;i<size;i++)
		  {
			NEWI:
			if (current_map[i] != -1)
			  {
				i_len = current_map[i];
				for (j=i+i_len;j<size;j++)
				  {
					if (current_map[j] != -1)
					  {
						j_len = current_map[j];
						// merge i & j
						current_map[j] = -1;
						current_map[i] = i_len+j_len;
						dumpCurrentList();
						// continue on past merge						
						i=j+j_len-1; // continue on past merge
						break NEWI;
					  }
				  }
			  }
		  }
	  }
	// reverse order simplest to most complex
	l = points.size();
	int[] tmp;
	for (i=0,j=l-1;i<l/2;i++,j--)
	  {
		tmp = points.get(i);
		points.set(i,points.get(j));
		points.set(j,tmp);
	  }
	if (sequenceAsSmallToBig)
	  {
		int t;
		for (i=0;i<l;i++)
		  {
			tmp = points.get(i);
			for(j=0;j<(tmp.length/2);j++)
			  {
				t = tmp[j];
				tmp[j] = tmp[tmp.length-1-j];
				tmp[tmp.length-1-j] = t;
			  }
		  }
	  }

	// printing for visuals...
	if (showTree)
	  {
		int k;
		int[] chunks;
		for (i=0;i<points.size();i++)
		  {
			chunks = points.get(i);
			System.out.printf(".");
			for (j=0; j<chunks.length; j++)
			  {
				System.out.printf("|");
				for (k=1;k<chunks[j];k++) {System.out.printf(" ");}
			  }
			System.out.println(".");
		  }
		for (i=0;i<points.size();i++)
		  {
			chunks = points.get(i);
			for(j=0;j<chunks.length;j++) System.out.printf("%3d ",chunks[j]);
			System.out.printf("%n");
		  }
	  }
  }

  void dumpCurrentList()
  {
	int i,j;
	
	int L=0;
	for (i=0; i<n_pts; i++) if (current_map[i] != -1) L++;
	int[] segment_lengths = new int[L];

	for (i=0,j=0; i<n_pts; i++)
	  {
		if (current_map[i] != -1) 
		  {
			segment_lengths[j] = current_map[i];
			j++;
		  }
	  }
	points.add(segment_lengths);
  }
}