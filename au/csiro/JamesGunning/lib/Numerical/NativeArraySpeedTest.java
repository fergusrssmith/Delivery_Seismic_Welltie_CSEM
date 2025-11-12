package au.csiro.JamesGunning.lib.Numerical;

import java.io.*;

import org.ejml.data.*;
import org.ejml.factory.*;
import org.ejml.alg.dense.mult.*;
import org.ejml.alg.dense.decomposition.*;
import org.ejml.interfaces.decomposition.*;
import org.ejml.ops.*;
 
/**
   Basic class to test native blas1 type speeds in hand-unrolled dgemm type thing.
   Also ejml fragments including Cholesky.
   Recent compilers claim they can do things to improve matters
   here. This *very simple* class tests this claim. 
   Compare to the colt  BenchmarkMatrix class.

   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/5/2003 
*/

public final class NativeArraySpeedTest
{
  public static void main(String[] args)
  {
	int i,j,k,l,ii;
	double[][] A,B;
	double[] At,Att,Bt;
	int[] sizes={5, 20, 50,100, 300, 500,1000,1500,2000,2500};
	long time;
	double perf,nReps,u;

	System.out.printf("Raw speeds, in hand unrolled dgemm%n");
	for (ii=0;ii<sizes.length;ii++)
	  {
		l=sizes[ii];
		A = new double[l][l];
		B = new double[l][l];
		// fill
		for(j=0;j<l;j++) {
		  for(k=0;k<l;k++) {
			A[j][k]=Math.random();
		  }
		}
		// time
		time=System.currentTimeMillis();
		nReps=0;
		do 
		  {
			for(i=0;i<l;i++) {
			  Bt=B[i];
			  At=A[i];
			  for(j=0;j<l;j++) {
				Att=A[j];
				for(k=0,u=0.0;k<l;k++) {
				  // a clever JIT ought to be able to do good stuff here
				  // B = A . A^T
				  u += At[k]*Att[k];
				}
				Bt[j] +=u;
			  }
			}
			nReps++;
		  }
		// at least 1 second per size 
		while ((System.currentTimeMillis()-time)<1000);
		perf = 1e-3*nReps*l*l*l/((double)(System.currentTimeMillis()-time));
		System.out.printf("DGEMM: Size=%4d, Mflops= %7.2f%n",l,perf);

		// time
		/*
		  time=System.currentTimeMillis();
		  nReps=0;
		  do 
		  {
		  for(i=0;i<l;i++) {
		  double[] Al=A[i];
		  double[] Bl=B[i];
		  for(j=0;j<l;j++) {
		  double[] Ar=A[j];
		  double s=0;
		  for(k=0;k<l;k++) {
		  // a clever JIT ought to be able to do good stuff here
		  // B = A . A^T
		  s += Al[k]*Ar[k];
		  }
		  Bl[j]=s;
		  }
		  }
		  nReps++;
		  }
		  // at least 1 second per size
		  while ((System.currentTimeMillis()-time)<1000);
		  perf = 1e-3*nReps*l*l*l/((double)(System.currentTimeMillis()-time));
		  System.out.println("DGEMM: Hard-coded Size="+l+" Mflops= "+perf);
		*/
	  }
	/* 
	   // EJML tests: about 10-20% faster
	for (ii=0;ii<sizes.length;ii++)
	  {
		l=sizes[ii];
		DenseMatrix64F A1 = new DenseMatrix64F(l,l);
		DenseMatrix64F B1 = new DenseMatrix64F(l,l);
		DenseMatrix64F C1 = new DenseMatrix64F(l,l);
		for(j=0;j<l;j++) {
		  for(k=0;k<l;k++) {
			A1.unsafe_set(j,k,Math.random());
			B1.unsafe_set(j,k,Math.random());
		  }
		} 

		// time
		time=System.currentTimeMillis();
		nReps=0;
		do 
		  {
			CommonOps.mult(A1,B1,C1);
			nReps++;
		  }
		// at least 1 second per size 
		while ((System.currentTimeMillis()-time)<1000);
		perf = 1e-3*nReps*l*l*l/((double)(System.currentTimeMillis()-time));
		System.out.printf("EJML DGEMM: Size=%4d, Mflops= %7.2f%n",l,perf);
	  } 
	*/
 
	/* EJML cholesky tests */
	System.out.printf("Cholesky benchmarks, ejml%n");
	for (ii=0;ii<sizes.length;ii++)
	  {
		l=sizes[ii];

		DenseMatrix64F H = new DenseMatrix64F(l,l);
	
		// time
		long t=0;
		time=0;
		nReps=0;
		do 
		  {
			for(i=0;i<l;i++)
			  for(j=0;j<l;j++) 
				H.unsafe_set(i,j,(i==j)?1.0 : 0.01);
	
			t = System.currentTimeMillis();
			CholeskyDecomposition<DenseMatrix64F> chol = DecompositionFactory.chol(l,true);
			
			if( !chol.decompose(H))
			  throw new RuntimeException("Cholesky failed!");
			chol.getT(null);		
			time += (System.currentTimeMillis()-t);
		
			nReps++;
		  }
		// at least 1 second per size	   
		while (time<1000);
		double T= time/(double)nReps;
		// work = (n^3+2n)/3
		perf = 1e-6*0.33333*(Math.pow(l,3.0)+2.0*l)/T;
		System.out.printf("Cholesky: Size=%4d, <T>=%4f Gflops= %7.2f%n",l,T,perf);	
	  }
  }
}
