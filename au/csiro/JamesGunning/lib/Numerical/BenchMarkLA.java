/**
   Primitive arrays DGEMM benchmark. Recent jits claim they know how to do
   bounds checking & optimisation. This will check. Compare to
   cern.colt.matrix.bench.BenchmarkMatrix dgemm dense 1 2.0 0.999 false true 5 10 25 50 100 250 500
*/
package au.csiro.JamesGunning.lib.Numerical;

public class BenchMarkLA
{
  public static void main(String argv[])
  {
	int[] sizes={64,128,256,512};

	int i,I,j,k,n,N;
	
	double[][] a,b,c;
	double[] elements,r_elems,c_elems;
	long t;


	for (I=0;I<sizes.length;I++)
	  {
		N=sizes[I];
		a = new double[N][N];
		b = new double[N][N];
		c = new double[N][N];

		for (i=0;i<N;i++)
		  {
			for(j=0;j<N;j++)
			  {
				a[i][j]=b[i][j]=Math.random();
			  }
		  }
		t=System.currentTimeMillis();	
		n=0;
		do
		  {
			for (i=0;i<N;i++)
			  {
				c_elems = c[i];
				elements = a[i];
				for(j=0;j<N;j++)
				  {
					r_elems = b[j];
					double s=0;
					for(k=0;k<N;k++) 
					  {
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
						++k;
						s += elements[k] * r_elems[k];
					  }
					c_elems[j]=s;
				  }
			  }
			n++;
		  }
		// do ops for at least a second
		while (System.currentTimeMillis()-t<4000);
		System.out.println("N = "+N+", Mflops = "+ ((1e-3)*2*N*N*N*n)/((System.currentTimeMillis()-t)));
	  }
  }
}
