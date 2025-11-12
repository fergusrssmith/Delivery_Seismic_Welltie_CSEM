package au.csiro.JamesGunning.lib.util;

import java.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;



/**
   Useful Class for indexing and sorting doubles


   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/8/2001
*/
public class SortableNumbers implements Comparable 
{
  public int index;
  public double value;
  
  public SortableNumbers(double value, int index)
  {
	this.value = value;;
	this.index = index;
  }

  public SortableNumbers(float value, int index)
  {
	this.value = (double)value;
	this.index = index;
  }

  public void changeValue(double value)
  {
	this.value = value;
  }

  /**
	   Comparator to sort by value	
  */
  public int compareTo(Object o)
  {
    SortableNumbers n = (SortableNumbers)o;
	int cmp = (value<n.value)? -1 : ((value>n.value)? 1 : 0 );
	return cmp;
  }

  // to use this class as keys in a Hashmap, we need:
  // (see the implementation in the Double class)
  public int hashCode()
  {
	long v = Double.doubleToLongBits(value);
               
	return (int)(v^(v>>>32));
  }

  // to use this class as keys in a Hashmap, we need:
  public boolean equals(Object obj)
  {
	if (obj == null) return false;

	if(getClass() == obj.getClass())  
	  {
		if (value == ((SortableNumbers)obj).value)
		  {
			return true;
		  }
		else
		  return false;
	  }
	else
	  return false;
  }



  public static class IndexSorter implements Comparator
  {
	/**
	   Comparator class to sort by index
	*/
	public int compare(Object o1,Object o2)
	{
	  SortableNumbers n1 = (SortableNumbers)o1;
	  SortableNumbers n2 = (SortableNumbers)o2;
	  int cmp = (n2.index<n1.index)? 1 : ((n2.index>n1.index)? -1 : 0 );
	  return cmp;
	}
  }
}


