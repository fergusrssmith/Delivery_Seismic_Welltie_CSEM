package au.csiro.JamesGunning.lib.util;

import java.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
 


/**
   Useful Class for sorting general n-tuples


   @author james.gunning@csiro.csiro.au
   @version 1.0, 2/2/2011
*/
public class TupleComparator implements Comparator<double[]> 
{
  private int[] s;
  private int n;

  public TupleComparator(int[] indices)
  {
	s=indices;
	n=s.length;
  }

  public int compare(double[] c1, double[] c2) 
  { 
	double o1,o2;
	for(int i=0;i<n;i++)
	  {
		o1 = c1[s[i]];
		o2 = c2[s[i]];
		if (o1 == o2)
		  {
			continue;
		  }
		else
		  {
			return (o1<o2) ? 1 : -1;
		  }
	  }
	return 1;
  }
}


