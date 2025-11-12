package au.csiro.JamesGunning.lib.util;
import java.util.Arrays;

/**	
Basic class of static public functions which are useful

@author james.gunning@dpr.csiro.au
@version 1.0, 9/8/2003 */
 
public class ArrayUtils
{

  /**
	 Deep copy of double[][] primitive arrays
  */
  public static final double[][] deepCopy(double[][] x)
  {
	double[][] xCp = (double[][]) x.clone();
	for (int i=0; i<x.length; i++) xCp[i] = (double[]) x[i].clone();
	return xCp;
  }
  
  /**
	 Deep copy of double[][][] primitive arrays
  */
  public static final double[][][] deepCopy(double[][][] x)
  {
	double[][][] xCp = (double[][][]) x.clone();
	for (int i=0; i<x.length; i++) 
	  {
		xCp[i] = (double[][]) x[i].clone();
		for (int j=0; j<x[i].length; j++) 
		  {
			xCp[i][j] = (double[]) x[i][j].clone();
		  }
	  }	
	return xCp; 
  }
  
  /**
	 Deep zero of double[][] primitive arrays
  */
  public static final void zero(double[][] x)
  {
	int i,N;
	N=x.length;
	for (i=0; i<N; i++) Arrays.fill(x[i],0.0);
  }
  
  /**
	 Deep zero of double[][][] primitive arrays
  */
  public static final void zero(double[][][] x)
  {
	int i,j,N;
	double[][] xi;
	for (i=0; i<x.length; i++) 
	  {
		xi = x[i];
		N = x[i].length;
		for (j=0; j<N; j++) Arrays.fill(xi[j],0.0);
	  }	
  }

  /**
	 Full wrapper to full System.arraycopy for doubles[]
	 Arrays must be of equal length
  */
  public static void copy(double[] src, double[] dest)
  {
	System.arraycopy(src,0,dest,0,src.length);
  }
 /**
	 Full wrapper to full System.arraycopy for floats[]
	 Arrays must be of equal length
  */
  public static void copy(float[] src, float[] dest)
  {
	System.arraycopy(src,0,dest,0,src.length);
  }

  /**
	 Full wrapper to full System.arraycopy for floats[][]
	 Arrays must be of equal length
  */
  public static void copy(double[][] src, double[][] dest)
  {
	for (int i=0; i<src.length; i++)
	  {
		System.arraycopy(src[i],0,dest[i],0,src[i].length);
	  }
  }
  
  /**
	 Full wrapper to full System.arraycopy for floats[][]
	 Arrays must be of equal length
  */
  public static void copy(float[][] src, float[][] dest)
  {
	for (int i=0; i<src.length; i++)
	  {
		System.arraycopy(src[i],0,dest[i],0,src[i].length);
	  }
  }

  public static double[] prepend(double x, double[] vals)
  {
	double[] vals_new = new double[vals.length+1];

	System.arraycopy(vals,0,vals_new,1,vals.length);
	vals_new[0]=x;

	return vals_new;
  }

  public static double[] append(double x, double[] vals)
  {
	double[] vals_new = new double[vals.length+1];

	System.arraycopy(vals,0,vals_new,0,vals.length);
	vals_new[vals.length]=x;

	return vals_new;
  }
  /**
	 return range of data over index interval named (inclusive)
  */
  public static double[] subselect(int iStart, int iEnd, double[] vals)
  {

	if (iStart<0)
	  throw new RuntimeException("Index out of bounds: " +iStart);
	if (iEnd>=vals.length)
	  throw new RuntimeException("Index out of bounds: " +iEnd);
	if (iEnd<iStart)
	  throw new RuntimeException("Must have iStart<=iEnd in subselect()");

	double[] vals_new = new double[iEnd-iStart+1];

	System.arraycopy(vals,iStart,vals_new,0,vals_new.length);

	return vals_new;
  }
}
