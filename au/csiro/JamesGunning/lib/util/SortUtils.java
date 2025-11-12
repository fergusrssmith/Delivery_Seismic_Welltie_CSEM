package au.csiro.JamesGunning.lib.util;

public class SortUtils implements Comparable 
{
  public int indx;
  public double x;
	
  public int compareTo(Object o)
  {
    SortUtils n = (SortUtils)o;
	int cmp = (x<n.x)? -1 : ((x>n.x)? 1 : 0 );
	return cmp;
  }
}
