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

// This is a translation & adaption of toms 695 by Elizabeth Eskow and 
// Robert B. Schnabel. Later revisions included.
//
// Author: James Gunning
// Date:   Feb 14/2006
//
/*
	From a symmetric, positive definite matrix A,
	constructs a 'safely' positive definite Cholesky Left decomposition
	using the methods described in
	"Algorithm 695 , collected algorithms from ACM.
	Transactions on mathematical software,
	Vol. 17, no. 3, September, 1991, pp. 306-312."

	Includes also updates in "A Revised Modified Cholesky Factorization Algorithm"
    Robert B. Schnabel, Elizabeth Eskow.
    SIAM Journal on Optimization, 9,No 4, Pages 1135 - 1148)
*/
//  Original comments.......
//********************************************************************************//
//
//  ALGORITHM 695 , COLLECTED ALGORITHMS FROM ACM.
//  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
//  VOL. 17, NO. 3, SEPTEMBER, 1991, PP. 306-312.
//  
//  
//  Driver for  modified cholesky factorization algorithm.


// *********************************************************************
//  
//  subroutine name: modchl
//  
//  authors :  Elizabeth Eskow and Robert B. Schnabel
//  
//  date    : December, 1988
//  
//  purpose : perform a modified cholesky factorization
//  of the form (Ptranspose)AP  + E = L(Ltranspose),
//  where L is stored in the lower triangle of the
//  original matrix A.
//  The factorization has 2 phases:
//  phase 1: Pivot on the maximum diagonal element.
//  Check that the normal cholesky update
//  would result in a positive diagonal
//  at the current iteration, and
//  if so, do the normal cholesky update,
//  otherwise switch to phase 2.
//  phase 2: Pivot on the minimum of the negatives
//  of the lower gerschgorin bound
//  estimates.
//  Compute the amount to add to the
//  pivot element and add this
//  to the pivot element.
//  Do the cholesky update.
//  Update the estimates of the
//  gerschgorin bounds.
//  
//  n       - dimension of matrix A
//  
//  A       - n*n symmetric matrix (only lower triangular
//  portion of A, including the main diagonal, is used)
//  
//  g       - n*1 work array
//  
//  mcheps - machine precision
//  
//  tau1    - tolerance used for determining when to switch
//  phase 2
//  
//  tau2    - tolerance used for determining the maximum
//  condition number of the final 2X2 submatrix.
//  
//  
//  output  : L     - stored in the matrix A (in lower triangular
//  portion of A, including the main diagonal)
//  
//  P     - a record of how the rows and columns
//  of the matrix were permuted while
//  performing the decomposition
//  
//  E     - n*1 array, the ith element is the
//  amount added to the diagonal of A
//  at the ith iteration
// 
// 
// ***********************************************************************

package au.csiro.JamesGunning.lib.Numerical;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;

import java.io.*;
import java.util.*;
import java.util.Random;

public class ModifiedCholeskyLeftDecomposition 
{
  /**
	 Machine epsilon. Usually about 1.12 x 10^-16
  */
  public final static double mcheps = MachineConstants.EPSILON_SMALL;
  private final static double tau = 4.80621738393e-6; // Math.pow(mcheps,1.0/3.0);
  private final static double taubar = 2.309972554166e-11; // mcheps^2/3
  private final static double MU = 0.1;  

  private int[] P;
  private DoubleMatrix2D _A;
  private int n;   // size of A

  // NEW mcld_globals
  private boolean _preScale=false;

  private boolean phase1 = true;
  private double delta=0;
  private double ming=0;
  private double taugam=0;
  private double taubargam=0;
  private double mugam=0;
  private double[] diagScale=null;

  /** 
	  Example test harness. Supply -s to scale
	  java au.csiro.JamesGunning.lib.Numerical.ModifiedCholeskyLeftDecomposition [-s]
  */
  public static void main(String argv[])
  {
	int i,j,k,n=0;
	double lambda,tmp;
	ModifiedCholeskyLeftDecomposition chol;
	DoubleMatrix2D A=null,R,C ,B;
	boolean useScaling=false;
	ErrorHandler EH;
	boolean readFromFile=false,doRandomTests=false;
	String filename=null;
	cern.colt.matrix.linalg.Blas blas;
	long seed=1;

	EH = new ErrorHandler();
	OptionsParser options = new OptionsParser("ModifiedCholeskyLeftDecomposition",EH);


	options.addOption("-s",null,0,null,null,"use diagonal pre-scaling");
	options.addOption("-f",null,1,"String","filename",
					  "ascii n x n matrix to decompose");	
	options.addOption("-r",null,1,"int","integer-seed",
					  "Do large set of random matrix tests");
	options.parseOptions(argv,true);

	if (options.isOptionSet("-r")) 
	  {
		doRandomTests = true;
		seed = ((Integer)options.getOptionValues("-r")[0]).intValue();
	  }
	if (options.isOptionSet("-f")) 
	  {
		readFromFile = true;
		filename = ((String)options.getOptionValues("-f")[0]);
	  }
	if (options.isOptionSet("-s")) 
	  {
		useScaling=true;
		System.out.println("Using scaling");
	  }

	blas = cern.colt.matrix.linalg.SeqBlas.seqBlas;

	if (!doRandomTests)
	  {
		if (readFromFile)
		  {
			// Read in square ascii matrix. Can be one column or nxn
			RandomAccessFile f=null;
			try
			  {
				f = new RandomAccessFile(filename, "r");	
			  }
			catch (Exception e) 
			  {
				EH.die_err("Can't open file "+filename);
			  }
			cern.colt.list.DoubleArrayList a = new cern.colt.list.DoubleArrayList();
			double val=0;
			String s;
			StringTokenizer st;
			int lineNo=1,nt=0;
			try
			  {
				while ((s=f.readLine()) != null)
				  {
					st = new StringTokenizer(s);
					while (st.hasMoreTokens()) 
					  {
						try
						  {
							val = Double.parseDouble(st.nextToken());
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
						  }
						a.add(val);
						nt++;
					  }
					lineNo++;
				  }
				f.close();
			  }
			catch (Exception e) 
			  {
				EH.die_err("I/O exception on "+filename);
			  }

			n = (int)Math.sqrt(a.size());
			if ((n*n) != a.size())
			  {
				EH.die_err("Matrix in "+filename+" is not square!\n");
			  }
			A = new DenseDoubleMatrix2D(n,n);
			for (i=0,k=0;i<n;i++) {
			  for (j=0;j<n;j++,k++) {
				A.setQuick(i,j,a.getQuick(k));
			  }
			}
		  }
		else
		  {
			// Make up a default: a funky just-positive definite matrix.

			n = 6;
		
			A = new DenseDoubleMatrix2D(n,n);
			double[] scale = new double[n]; 
			double scale_eps = 1.0e-4;  // 1e-12
			double corr_range=10; // 100
			for (i = 0; i < n; i++) 
			  {
				scale[i] = (i<n/3) ? 1.0-(double)i/(10.0*n) : ((i>(0.6666666*n))? 1.0-(double)i/(10.0*n) : scale_eps);
			  }

			for (i=0;i<n;i++) {
			  for (j=0;j<n;j++) {
				A.setQuick(i,j,scale[i]*scale[j]*Math.exp(-Math.abs(((double)(i-j))/corr_range)));
			  }
			}
		  }

		System.out.println("Original A matrix = "+A.toString());



		// take a copy
		C = A.like().assign(A);
		DoubleMatrix2D Aorig = A.like().assign(A);

		// test inverse (sans corruption)...
		DoubleMatrix2D Cinv2 = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(A);
		

		// perform Left Cholesky
		chol = new ModifiedCholeskyLeftDecomposition(A,useScaling,false);
		A = chol.getL();
		int[] p = chol.getPermutation();
		System.out.println("L="+A.toString());
		System.out.println("Permutations:");
		for (i=0;i<n;i++) {
		  System.out.print(p[i]+" ");
		}
		// Check multiplication
		// B = A.A^T
		B = A.like();
		blas.dgemm(false,true,1.0,A,A,0,B);
		// permute
		double max=0;
		for (i=0;i<n;i++) {
		  for (j=0;j<n;j++) {
			tmp = B.getQuick(i,j);
			A.setQuick(p[i],p[j], tmp);
			max=Math.max(max, Math.abs(tmp - C.getQuick(p[i],p[j])));
		  }
		}
		//	System.out.println("\nMultiplied out after perms = "+A.toString());
		System.out.println("\nMax error after multiplying out L.L^T = "+max);

		// check left inverse
		chol = new ModifiedCholeskyLeftDecomposition(C,useScaling,false);
		A = chol.getL();
		B.assign(A);
		chol.getCholeskyLeftInverse(B);
		System.out.println("L^{-1}="+B.toString());
		C=A.like();
		blas.dgemm(false,false,1.0,A,B,0,C);
		System.out.println("\nCheck of L.L^{-1} = I:");
		System.out.println("L.L^{-1}="+C.toString());
		max=0;
		for (i=0;i<n;i++) {
		  for (j=0;j<n;j++) {
			max = Math.max(max, Math.abs(((i==j)?1:0) - C.getQuick(i,j)));
		  }
		}
		System.out.println("Max. inverse discrepancy="+max);

		System.out.println("\nMatrix inverse from cholesky L^{-1}T.L^{-1}:");
		//System.out.println("L^{-1}="+B.toString());
		//System.exit(0);
		blas.dgemm(true,false,1.0,B,B,0,C);
		for (i=0;i<n;i++) {
		  for (j=0;j<n;j++) {
			tmp = C.getQuick(i,j);
			A.setQuick(p[i],p[j], tmp);
		  }
		}
		System.out.println("C^{-1}="+A.toString());

		B = (new cern.colt.matrix.linalg.Algebra()).inverse(Aorig);
		System.out.println("Conventional inverse C^{-1}="+B);

		blas.dgemm(false,false,1.0,A,Aorig,0,C);
		max=0;
		for (i=0;i<n;i++) {
		  for (j=0;j<n;j++) {
			max = Math.max(max, Math.abs(((i==j)?1:0) - C.getQuick(i,j)));
		  }
		}
		System.out.println("Max error, from Cholesky derived inverse, in (C^{-1}.C - I) is "+max);

		System.out.println("Inverse, using getInverseViaCholesky()");

		System.out.println(Cinv2.toString());

		System.out.printf("Eigenvalues: ");
		EigenvalueDecomposition eig = new EigenvalueDecomposition(Aorig); 
		DoubleMatrix1D eR = eig.getRealEigenvalues();
		for (i=0;i<eR.size();i++) System.out.printf("%g ",eR.getQuick(i));
		System.out.println("");
		double eMax=-Double.MAX_VALUE,eMin=Double.MAX_VALUE;
		for (i=0;i<eR.size();i++)
		  {
			double x = Math.abs(eR.getQuick(i));
			eMax = Math.max(eMax,x);
			eMin = Math.min(eMin,x);
		  }	
		System.out.printf("Condition number = %g%n",eMax/eMin);
	  }
	else
	  {
		// random tests
		int m,iRank;
		double maxE=0,maxEF=0,rowscale;
		Random rg = new Random(seed);
		DoubleMatrix2D X=null, L=null, Aorig=null;
		int[] p=null;

		for (iRank=0;iRank<=1;iRank++)
		  {
			// form A=X.X^T, with X = n x m
			// for a variety of cases
			for (n=10; n<=30; n+=5)
			  {
				for (m=n-5; m<=n; m++)
				  {
					X= new DenseDoubleMatrix2D(n,m);
					A = new DenseDoubleMatrix2D(n,n);
					for (i=0;i<n;i++) {
					  rowscale = Math.exp(10.0*rg.nextDouble());
						for (j=0;j<m;j++) {
						  X.setQuick(i,j,rowscale * rg.nextDouble());
						}
					}
					blas.dgemm(false,true,1.0,X,X,0,A);
					Aorig = A.like().assign(A);
					// perform Left Cholesky
					if (true)
					  {
						chol = new ModifiedCholeskyLeftDecomposition(A,useScaling,false);
						L = chol.getL();
						p = chol.getPermutation();
					  }
					else
					  {
						cern.colt.matrix.linalg.CholeskyDecomposition cernChol = new cern.colt.matrix.linalg.CholeskyDecomposition(A);
						L = cernChol.getL();	
						if (!cernChol.isSymmetricPositiveDefinite())
						  {
							System.out.println("\nFailed Pos. Def test");
							continue;
						  }
						p = new int[n];
						for (i=0;i<n;i++) p[i]=i;
					  }
					// Check multiplication
					// B = L.L^T
					B = A.like();
					blas.dgemm(false,true,1.0,L,L,0,B);
					// permute
					double max=0;
					for (i=0;i<n;i++) {
					  rowscale =B.getQuick(i,i);
					  for (j=0;j<n;j++) {
						tmp = B.getQuick(i,j);
						max = Math.max(max, Math.abs(tmp - Aorig.getQuick(p[i],p[j]))/rowscale);
					  }
					}
					System.out.println("\nMax error after multiplying out L.L^T = "+max);
					maxE = Math.max(maxE,max);
					if (n==m) maxEF = Math.max(maxEF,max);
				  }
			  }
		  }
		System.out.println("\nWorst error over test set = "+maxE);
		System.out.println("\nWorst error over full rank test set = "+maxEF);
	  }
  }

  /** 
      <p>
      From a symmetric, positive definite matrix <tt>A</tt>,
      constructs a 'safely' positive definite Cholesky Left decomposition
      using the methods described in
      "Algorithm 695 , collected algorithms from ACM.
      Transactions on mathematical software,
      Vol. 17, no. 3, September, 1991, pp. 306-312."
      Guaranteed to rescue
      you if the matrix is ill-conditioned or not quite positive definite.
      </p>
      <p>
      The decomposition is of form P<sup>T</sup>.A.P = L.L<sup>T</sup>
      </p>
      The permutation P is available through getPermutation(), and the
      factor L in getL()
      </p>
      <p>
      Example Usage:<br>
      <pre>
      // set up some matrix A....,then
      chol = new ModifiedCholeskyLeftDecomposition(A);
	  // or, for wierdly scaled matrices
      chol = new ModifiedCholeskyLeftDecomposition(A,true);
	  L = chol.getL();
      int[] p = chol.getPermutation();
      // form B = L.L<sup>T</sup> using the COLT Blas routines
      B = L.like(); 
      cern.colt.matrix.linalg.Blas blas = cern.colt.matrix.linalg.SeqBlas.seqBlas;
      // form  L.L<sup>T</sup>
      blas.dgemm(false,true,1.0,L,L,0,B);
      // unpermute to form  P.L.L<sup>T</sup>.P<sup>T</sup> which should look like A
      for (i=0;i&lt n;i++) {
        for (j=0;j &lt n;j++) {
          tmp = B.getQuick(i,j);
          A.setQuick(p[i],p[j],tmp);
        }
      }
      System.out.println(A.toString());
      </pre>
    </p>	

	@param  A   Square matrix: Positive definite or close thereto.
	@throws IllegalArgumentException if <tt>A</tt> is not square.
  */
  public ModifiedCholeskyLeftDecomposition(DoubleMatrix2D A)
  {
	// initialise
	int i,j;

	this._A = A;

	n = _A.rows();
	if (_A.columns() != n) throw new IllegalArgumentException("matrix is not square");

	double[] g = new double[n];
	double[] E = new double[n];
	P = new int[n];

	// The toms original uses the same tau for tau1 and tau2
	modchl(n, _A, g, mcheps, tau, taubar, MU, P, E);

	// A final fix of Eskow's routine: erase the right triangle
	for (i=0;i<n;i++) {
	  for (j=i+1;j<n;j++) {
		_A.setQuick(i,j,0);
	  }
	}
  }
  /**
	 Same as ModifiedCholeskyLeftDecomposition(), except if
	 preScale==true, the matrix is scaled by the diagonal before
	 decomposition, and the scale vector is used to rescale L after
	 the decomposition. If scaleWithinMachEpsRange==true, scaling
	 applies to rows j whose diagonal element satisfies |A_jj| >
	 machine_eps * max(A_ii), i.e. within a factor of say
	 10^-16. Smaller numbers are assumed to best noise and are
	 unscaled.  If scaleWithinMachEpsRange==false, no thresholding on
	 scaling applies, and all diagonal elemnets are used. Zeros are
	 scaled as if 1, though.
  */
  public ModifiedCholeskyLeftDecomposition(DoubleMatrix2D A, boolean preScale,boolean scaleWithinMachEpsRange)
  {
	// initialise
	int i,j;
	double s;

	this._A = A;
	this._preScale = preScale;

	n = _A.rows();
	if (_A.columns() != n) throw new IllegalArgumentException("matrix is not square");

	double[] g = new double[n];
	double[] E = new double[n];
	P = new int[n];

	if (preScale)
	  {
		diagScale = new double[n];
		double[] divScale =  new double[n];
		s = 0.0;
		if (scaleWithinMachEpsRange)
		  {
			// scaling seems to work best on numbers 
			// within *relative* dynamic range [1,mcheps]. 
			// Smaller numbers are probably noise, and will
			// wreck the scaling
			for (i=0;i<n;i++) s = Math.max(s,Math.abs(_A.getQuick(i,i)));
			s *= mcheps;
		  }
		for (i=0;i<n;i++) 
		  {
			if (Math.abs(_A.getQuick(i,i)) > s)
			  {
				diagScale[i] = Math.sqrt(Math.abs(_A.getQuick(i,i)));
				divScale[i] = diagScale[i]; 
			  }
			else
			  {
				diagScale[i] = 1.0;
				divScale[i] = 1.0;
			  }
		  }
		for (i=0;i<n;i++) {
		  for (j=0;j<n;j++) {
			_A.setQuick(i,j,_A.getQuick(i,j)/(divScale[i]*divScale[j]));
		  }
		}
		divScale=null;
	  }

	// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("A.txt",_A);

	// The toms original uses the same tau for tau1 and tau2
	modchl(n, _A, g, mcheps, tau, taubar, MU, P, E);

	// A final fix of Eskow's routine: erase the right triangle
	for (i=0;i<n;i++) {
	  for (j=i+1;j<n;j++) {
		_A.setQuick(i,j,0);
	  }
	}
  
	// and a final unscaling if required
	if (preScale)
	  {
		// L=S.Lhat
		for (i=0;i<n;i++) {
		  s = diagScale[P[i]];
		  for (j=0;j<=i;j++) {
			_A.setQuick(i,j, s * A.getQuick(i,j));
		  }
		}
	  }
  }



  /** 
	  After the instantiator, this method can show if a 
	  phase2 calculation (modified matrix=true) was needed.
   */
  public boolean wasModified()
  {
	return !phase1;
  }

  /**
	 Call this to obtain the permutation P performed in the decomposition.
	 If unpermuted, the vector getPermutation() will be {0,1,2,3,,...}
	 The length of the permutation is A.rows().
  */
  public int[] getPermutation()
  {
	return P;
  }

 /**
	 Call this to obtain the Cholesky left factor L in P<sup>T</sup>.A.P = L.L<sup>T</sup>.
	 The permutation operator P is o<btained from getPermutation()
	 The dimensions of getL() are the same as A.
 */
  public DoubleMatrix2D getL()
  {
	return _A;
  }

  /**
	 Translation of Eskow's modified Cholesky routine; zero offset.
	 Private access: do not call from outside
  */
  void modchl(int n,DoubleMatrix2D A,double[] g,
			  double mcheps,double tau,double taubar,double mu, 
			  int[] P,double[] E)
  {
	//*********************************************************************
	//
	//       subroutine name: modchl
	//
	//       authors :  Elizabeth Eskow and Robert B. Schnabel
	//
	//       date    :  August, 1999 (revised from original code of 12/88)
	//
	//       purpose : perform a modified cholesky factorization
	//                 of the form (Ptranspose)AP  + E = L(Ltranspose),
	//       where L is stored in the lower triangle of the
	//       original matrix A.
	//       The factorization has 2 phases:
	//        phase 1: Pivot on the maximum diagonal element.
	//            Check that the normal cholesky update
	//            would result in a positive diagonal
	//            at the current iteration, and
	//            if so, do the normal cholesky update,
	//            otherwise switch to phase 2.
	//        phase 2: Pivot on the minimum of the negatives
	//            of the lower gerschgorin bound
	//            estimates.
	//            Compute the amount to add to the
	//            pivot element and add this
	//            to the pivot element.
	//            Do the cholesky update.
	//            Update the estimates of the
	//            gerschgorin bounds.
	//
	//       input   : ndim    - largest dimension of matrix that
	//                           will be used
	//
	//                 n       - dimension of matrix A
	//
	//                 A       - n*n symmetric matrix (only lower triangular
	//            portion of A, including the main diagonal, is used)
	//
	//                 g       - n*1 work array
	//
	//                 mcheps - machine precision
	//
	//                 mu     - tolerance used for determining when to switch
	//                          phase 2 via lookahead, or due to the  
	//                          relative magnitude of the current 
	//                          min. diag. element.
	//
	//                tau     - tolerance used for determining the maximum
	//                          condition number of the final 2X2 submatrix.
	//          
	//                taubar  - tolerance used for determining the amount
	//                          added to the smallest eigenvalue of the 
	//                          final 2X2 submatrix.
	//
	//
	//
	//       output  : L     - stored in the matrix A (in lower triangular
	//                           portion of A, including the main diagonal)
	//
	//                 P     - a record of how the rows and columns
	//                         of the matrix were permuted while
	//                         performing the decomposition
	//
	//                 E     - n*1 array, the ith element is the
	//                         amount added to the diagonal of A
	//                         at the ith iteration
	//
	//
	//***********************************************************************

	int j,iming,i,imaxd,itemp,jp1=0,k;
	double normj,maxd,mind;
	double delta1,temp,jdmin,tdmin,tempjj;
  
	init(n, A, P, g, E,taubar,mu);
	//  
	//  check for n=1
	//  
	if (n==1) 
	  {
		delta = (taubar * Math.abs(A.getQuick(0,0))) - A.getQuick(0,0);
		if (delta > 0) E[0]=delta;
		if (A.getQuick(0,0) == 0) E[0]=taubar;
		A.setQuick(0,0,  Math.sqrt(Math.max(0,A.getQuick(0,0) +E[0])) ); 
	  }

	for (j = 0;j < n;j++)
	  {
		//  
		//  PHASE 1
		//  
		if (phase1 ) 
		  {
			//  
			//  Find index of maximum diagonal element A(i,i) where i>=j
			//  
			mind = maxd = A.getQuick(j,j) ;
			imaxd = j;
			for (i = j+1;i<n;i++) 
			  {
				if (maxd < A.getQuick(i,i) ) 
				  {
					maxd = A.getQuick(i,i) ;
					imaxd = i;
				  }
				if (mind > A.getQuick(i,i)) mind = A.getQuick(i,i) ;
			  }
			// if (maxd <taubargam || (mind < -MU*maxd))
			if (mind < -mu*maxd)
			  {
				phase1 = false;
			  }
			//  Pivot to the top the row and column with the max diag
			// 
			// System.out.println("j="+j+", pivot="+imaxd);
			if (phase1 ) 
			  {			  
				if (imaxd != j) 
				  {
					//  Swap row j with row of max diag
					for (i=0; i<= j-1;i++) 
					  {
						temp = A.getQuick(j,i) ;
						A.setQuick(j,i,  A.getQuick(imaxd,i) ); 
						A.setQuick(imaxd,i,  temp ); 
					  }
					//  Swap colj and row maxdiag between j and maxdiag				
					for (i = j+1;i<=imaxd-1;i++)
					  {
						temp = A.getQuick(i,j) ;
						A.setQuick(i,j,  A.getQuick(imaxd,i) );                 
						A.setQuick(imaxd,i,  temp ); 
					  }
					//  Swap column j with column of max diag
					for (i = imaxd+1 ; i < n ; i ++) 
					  {
						temp = A.getQuick(i,j) ;
						A.setQuick(i,j,  A.getQuick(i,imaxd) ); 
						A.setQuick(i,imaxd,  temp ); 
					  }
					//  Swap diag elements
					temp = A.getQuick(j,j) ;
					A.setQuick(j,j,  A.getQuick(imaxd,imaxd) ); 
					A.setQuick(imaxd,imaxd,  temp ); 
					//  Swap elements of the permutation vector
					itemp = P[j];
					P[j]=P[imaxd];
					P[imaxd]=itemp;
				  }
				//  Check to see whether the normal cholesky update for this;
				//  iteration would result in a positive diagonal,
				//  and if not then switch to phase 2.
		  
				jp1 = j+1;
				tempjj=A.getQuick(j,j) ;	

				if (tempjj > taubargam) 
				  {  		
					if (j<(n-1))
					  {	
						jdmin=A.getQuick(jp1,jp1) ;
						for (i = jp1 ; i < n ; i ++) 
						  {
							temp = A.getQuick(i,j);
							temp = temp * temp/ tempjj;
							tdmin = A.getQuick(i,i)  - temp;
							jdmin = Math.min(jdmin, tdmin);
						  }
						if (jdmin < -mugam) phase1 = false;
					  }
				  }
				else
				  {
					phase1 = false;
				  }
			  }

			if (phase1) 
			  {
				//  Do the normal cholesky update if still in phase 1
				A.setQuick(j,j,  Math.sqrt(Math.max(0.0,A.getQuick(j,j))) ); 
				if (j<(n-1))
				  {
					tempjj = A.getQuick(j,j) ;
					for (i = jp1 ; i < n ; i ++) 
					  {
						A.setQuick(i,j, A.getQuick(i,j)/tempjj); 
					  }
					for (i=jp1 ; i<n ; i++) 
					  {
						temp=A.getQuick(i,j) ;
						for (k = jp1 ; k <= i ; k++) 
						  {
							A.setQuick(i,k, A.getQuick(i,k) - (temp * A.getQuick(k,j))); 
						  }
					  }
				  }
			  }
			else 
			  {
				//  Calculate the negatives of the lower gerschgorin bounds
				gersch(n,A,j,g);
			  }
		  }

		//  PHASE 2
		// System.out.println("j=%d, phase1=%d\n",j,phase1);	// FIXME
		if (!phase1 && (j == n-1)) 
		  {
			// See algorithm in paper
			delta1= -A.getQuick(n-1,n-1) + Math.max(-tau*A.getQuick(n-1,n-1)/(1.0-tau), 
											taubargam);
			delta = Math.max(delta, delta1);
			E[n-1] = delta;
			A.setQuick(n-1,n-1, Math.sqrt(Math.max(0,A.getQuick(n-1,n-1) + E[n-1]))); 
			return;
		  }
		if (!phase1)
		  {
			if (j < (n-2))
			  {
				//  Find the minimum negative gershgorin bound
				iming=j;
				ming = g[j];
				for (i = j+1 ; i <n ; i ++) 
				  {
					if (ming > g[i]) 
					  {
						ming = g[i];
						iming = i;
					  }
				  }
				//  Pivot to the top the row and column with the
				//  minimum negative gerschgorin bound
				if (iming != j) 
				  {
					//  Swap row j with row of min gersch bound
					for (i = 0 ; i<= j-1;i++) 
					  {
						temp = A.getQuick(j,i) ;
						A.setQuick(j,i,  A.getQuick(iming,i)   );
						A.setQuick(iming,i,  temp  );
					  }
					//  Swap colj with row iming from j to iming
					for (i = j+1 ; i <=iming-1 ; i ++) 
					  {
						temp = A.getQuick(i,j) ;
						A.setQuick(i,j,  A.getQuick(iming,i)   );
						A.setQuick(iming,i,  temp  );
					  }
					//  Swap column j with column of min gersch bound
					for (i = iming+1 ; i < n ; i ++) 
					  {
						temp = A.getQuick(i,j) ;
						A.setQuick(i,j,  A.getQuick(i,iming)   );
						A.setQuick(i,iming,  temp  );
					  }
					//  Swap diagonal elements
					temp = A.getQuick(j,j) ;
					A.setQuick(j,j,  A.getQuick(iming,iming)   );
					A.setQuick(iming,iming,  temp  );
					//  Swap elements of the permutation vector
					itemp = P[j];
					P[j]=P[iming];
					P[iming]=itemp;
					//  Swap elements of the negative gerschgorin bounds vecto
					temp = g[j];
					g[j] = g[iming];
					g[iming]=temp;
				  }
				//  Calculate delta and add to the diagonal.
				//  delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
				//  where normj=sum of |A(i,j)|,for i=1,n,
				//  delta_previous is the delta computed at the previous iter
				//  and taugam is tau1*gamma.
				normj = 0.0;
				for (i = j+1 ; i < n ; i ++) 
				  {
					normj += Math.abs(A.getQuick(i,j));
				  }
				temp = Math.max(normj,taubargam);
				delta1 = temp - A.getQuick(j,j) ;
				temp = 0.0;
				delta1 = Math.max(temp, delta1);
				delta = Math.max(delta1,delta);
				E[j]=delta;
				A.setQuick(j,j, A.getQuick(j,j) + E[j]) ;
				//  Update the gerschgorin bound estimates
				//  (note: g(i) is the negative of the
				//  Gerschgorin lower bound.)
				//  
				if (A.getQuick(j,j)  != normj) 
				  {
					temp = (A.getQuick(j,j) != 0) ? (normj/A.getQuick(j,j)) - 1.0 : -1.0; // assume norm=0 if A.getQuick(j,j)==0
					for (i = j+1 ; i < n ; i ++) 
					  {
						g[i] += Math.abs(A.getQuick(i,j))* temp;
					  }
				  }
				//  Do the cholesky update
				  
				A.setQuick(j,j,  Math.sqrt(Math.max(0,A.getQuick(j,j)))  );
				tempjj = A.getQuick(j,j) ;
				if (tempjj>0)
				  {
					for (i = j+1 ; i < n ; i ++) 
					  {
						A.setQuick(i,j, A.getQuick(i,j)/tempjj);
					  }
					for (i = j+1 ; i < n ; i ++) 
					  {
						temp = A.getQuick(i,j) ;
						for (k = j+1 ; k <= i ; k ++) 
						  {
							A.setQuick(i,k, A.getQuick(i,k) - (temp * A.getQuick(k,j))) ;
						  }
					  }
				  }
				else
				  {
					for (i = j+1 ; i < n ; i ++) A.setQuick(i,j,  0.0 );
				  } 
			  }
			else  
			  {
				fin2by2(n, A, E, j, tau);
				return;
			  }
		  }
	  }
  }

  //  ***********************************************************************
  //  subroutine name : init
  //  
  //  purpose : set up for start of cholesky factorization
  //  
  //  input : n, A, tau1
  //  
  //  sets: phase1    - short value set to 1 if in phase one,
  //  otherwise 0.
  //  delta     - amount to add to Ajj at iteration j
  //  P,g,E - described above in modchl
  //  ming      - the minimum negative gerschgorin bound
  //  gamma     - the maximum diagonal element of A
  //  taugam  - tau1 * gamma
  //  
  //  ***********************************************************************
  void init(int n,DoubleMatrix2D A,int[] P,double[] g,double[] E,
			double taubar,double mu)
  {
	int i;
	double gamma;

	phase1 = true;
	delta = 0.0;
	ming = 0.0;
	for (i=0 ; i<n ; i++) 
	  {
		P[i]=i;
		g[i]=0.0;
		E[i]=0.0;
	  }
	//  
	//  Find the maximum magnitude of the diagonal elements.
	//  
	gamma = 0.0;
	for (i=0 ; i<n ; i++)
	  {
		gamma = Math.max(gamma,Math.abs(A.getQuick(i,i)));
	  }
	taubargam = taubar * gamma;
	mugam = mu * gamma;
  }




  //  ***********************************************************************
  //  
  //  subroutine name : gersch
  //  
  //  purpose : Calculate the negative of the gerschgorin bounds
  //  called once at the start of phase II.
  //  
  //  input   : n, A, j
  //  
  //  output  : g - an n vector containing the negatives of the
  //  Gerschgorin bounds.
  //  
  //  ***********************************************************************
  void gersch(int n, DoubleMatrix2D A, int j, double[] g)
  {
	int i, k;
	double offrow;

	for (i = j ; i < n ; i ++) 
	  {
		offrow = 0.0;
		for (k = j ; k <= i-1 ; k ++) 
		  {
			offrow += Math.abs(A.getQuick(i,k) );
		  }
		for (k = i+1 ; k < n ; k ++) 
		  {
			offrow += Math.abs(A.getQuick(k,i) );
		  }
		g[i] = offrow - A.getQuick(i,i) ;
	  }
  }

  //  ***********************************************************************
  //  
  //  subroutine name : fin2by2
  //  
  //  purpose : Handles final 2X2 submatrix in Phase II.
  //  Finds eigenvalues of final 2 by 2 submatrix,
  //  calculates the amount to add to the diagonal,
  //  adds to the final 2 diagonal elements,
  //  and does the final update.
  //  
  //  input : n, A, E, j, tau,
  //  delta - amount added to the diagonal in the
  //  previous iteration
  //  
  //  output : A - matrix with complete L factor in the lower triangle,
  //  E - n*1 vector containing the amount added to the diagonal
  //  at each iteration,
  //  delta - amount added to diagonal elements n-1 and n.
  //  
  //  ***********************************************************************
  void fin2by2(int n, DoubleMatrix2D A, double[] E, int j, double tau)
  {
	double t1, t2, t3,lmbd1,lmbd2,lmbdhi,lmbdlo;
	double delta1;
	//  
	//  Find eigenvalues of final 2 by 2 submatrix
	//  
	t1 = A.getQuick(n-2,n-2)  + A.getQuick(n-1,n-1) ;
	t2 = A.getQuick(n-2,n-2)  - A.getQuick(n-1,n-1) ;
	t3 = Math.sqrt(Math.max(0,t2*t2 + 4.0*A.getQuick(n-1,n-2) *A.getQuick(n-1,n-2)));
	lmbd1 = (t1 - t3)/2.;
	lmbd2 = (t1 + t3)/2.;
	lmbdhi = Math.max(lmbd1,lmbd2);
	lmbdlo = Math.min(lmbd1,lmbd2);
	//  
	//  Find delta such that:
	//  1.  the l2 condition number of the final
	//  2X2 submatrix + delta*I <= tau
	//  2. delta >= previous delta,
	//  3. lmbdlo + delta >= tau * gamma,
	//  where lmbdlo is the smallest eigenvalue of the final
	//  2X2 submatrix
	//

	// See algorithm in paper
	delta1= Math.max(tau*(lmbdhi-lmbdlo)/(1.0-tau), taubargam);
	delta1 -= lmbdlo;
	delta1 = Math.max(delta1, 0.0);
	delta = Math.max(delta, delta1);

	if (delta > 0.0) 
	  {
		A.setQuick(n-2,n-2, A.getQuick(n-2,n-2) + delta) ;
		A.setQuick(n-1,n-1, A.getQuick(n-1,n-1) + delta) ;
		E[n-2] = delta;
		E[n-1] = delta;
	  }
	//  
	//  Final update
	//  
	A.setQuick(n-2,n-2, Math.sqrt(Math.max(0,A.getQuick(n-2,n-2)))) ;
	A.setQuick(n-1,n-2, (A.getQuick(n-2,n-2) !=0 ) ? A.getQuick(n-1,n-2) /A.getQuick(n-2,n-2) : 0.0 );
	A.setQuick(n-1,n-1, A.getQuick(n-1,n-1)  - (A.getQuick(n-1,n-2) * A.getQuick(n-1,n-2)));
	A.setQuick(n-1,n-1, Math.sqrt(Math.max(0,A.getQuick(n-1,n-1))));
	return;
  }





  //  **********************************************************************
  //  
  //  subroutine name : solve
  //  
  //  purpose :  solves (LLtranspose)P(x)=P(b), where L is stored in
  //  the lower triangle of A, P is the record of the
  //  permutations  performed in forming L,b is the rhs vecto
  //  and x is the result.
  //  L is the result of the modified cholesky factorization
  //  which computes Ptranspose(A+E)P=LLtranspose.
  //  
  //  n     - dimension of matrix
  //  
  //  A     - n*n array (only lower triangle portion,
  //  including the diagonal, is used).
  //  
  //  P     - n*1 int vector that contains a record
  //  of permutations performed in forming L.
  //  i.e. each Pi was initialized to i, and if
  //  row and columns i & j were interchanged when
  //  computing L, the values of Pi and Pj were swapped.
  //  
  //  b     - n*1 double vector which is the
  //  right hand side of the system to be solved.
  //  
  //  y,z   - n*1 double  work vectors.
  //  
  //  output:  x     - n*1 double result vector
  //  
  //  
  //  **********************************************************************
  //
  // Solution checked to be OK; 19 Aug 2003 James Gunning
  /**
	 After the Cholesky Decomposition is instantiated, we can solve
	 the linear problem A.x=b using this call. If x==null when called
	 a new vector is returned with the solution, else the solution is written in place into the
	 provided x
	 
	 @param  A  x: vector of size = A.rows()
	 @param  A  b: rhs vector of size = b.rows()
	 
	 @throws IllegalArgumentException if b or a has the wrong size.
  */
  public DoubleMatrix1D solve(DoubleMatrix1D x, DoubleMatrix1D  b)
  {
	int i,j;
	double sum;
	double[] y,z;


	if (x==null) 
	  {
		x = new DenseDoubleMatrix1D(n);
	  }
	else if (x.size() != n)
	  {
		throw new IllegalArgumentException("X vector has incorrect number of elements");
	  }
	if (b.size() != n)
	  {
		throw new IllegalArgumentException("RHS vector has incorrect number of elements");
	  }

	y = new double[n];
	z = new double[n];

	//  
	//  Solve Ly=Pb
	//  
	y[0] = b.getQuick(P[0])/_A.getQuick(0,0);
	for (i=0 ; i<n ; i++) 
	  {
		sum=0.0;
		for (j=0 ; j<=i-1 ; j++) 
		  {
			sum += (_A.getQuick(i,j)*y[j]);
		  }
		y[i] = (b.getQuick(P[i]) - sum) / _A.getQuick(i,i);
	  }
	//  
	//  Solve Ltranspose z = y
	//  
	z[n-1] = y[n-1]/_A.getQuick(n-1,n-1);
	for (i=n-2; i>=0 ; i--) 
	  {
		sum=0.0;
		for (j=i+1 ; j<n ; j++) 
		  {
			sum += _A.getQuick(j,i) * z[j];
		  }
		z[i] = (y[i] - sum) / _A.getQuick(i,i);
	  }
	//  
	//  x = Ptranspose z
	//  
	for (i=0 ; i<n ; i++) 
	  {
		x.setQuick(P[i], z[i]);
	  }
	y=null;
	z=null;

	return x;
  }

  /**
	 After the Cholesky Decomposition is instantiated, we can solve
	 the linear problem A.x = b using this call. If x==null when called
	 a new vector is returned with the solution, else the solution is written in place into the
	 provided x. Same as other solve methods, but for doubles[].
	 
	 @param  A  x: native-style vector of x.length = A.rows()
	 @param  A  b: native-style vector of b.length = A.rows()
	 
	 @throws IllegalArgumentException if b or a has the wrong size.
  */
  public double[] solve(double[] x, double[] b)
  {
	int i,j;
	double sum;
	double[] y,z;


	if (x==null) 
	  {
		x = new double[n];
	  }
	else if (x.length != n)
	  {
		throw new IllegalArgumentException("X vector has incorrect number of elements");
	  }
	if (b.length != n)
	  {
		throw new IllegalArgumentException("RHS vector has incorrect number of elements");
	  }

	y = new double[n];
	z = new double[n];

	//  
	//  Solve Ly=Pb
	//  
	y[0] = b[P[0]]/_A.getQuick(0,0);
	for (i=0 ; i<n ; i++) 
	  {
		sum=0.0;
		for (j=0 ; j<=i-1 ; j++) 
		  {
			sum += (_A.getQuick(i,j)*y[j]);
		  }
		y[i] = (b[P[i]] - sum) / _A.getQuick(i,i);
	  }
	//  
	//  Solve Ltranspose z = y
	//  
	z[n-1] = y[n-1]/_A.getQuick(n-1,n-1);
	for (i=n-2; i>=0 ; i--) 
	  {
		sum=0.0;
		for (j=i+1 ; j<n ; j++) 
		  {
			sum += _A.getQuick(j,i) * z[j];
		  }
		z[i] = (y[i] - sum) / _A.getQuick(i,i);
	  }
	//  
	//  x = Ptranspose z
	//  
	for (i=0 ; i<n ; i++) 
	  {
		x[P[i]] = z[i];
	  }
	y=null;
	z=null;

	return x;
  }

  /**
	 After the Cholesky Decomposition is instantiated, we can solve
	 the linear problem A.x = b using this call. If x==null when called
	 a new vector is returned with the solution, else the solution is written in place into the
	 provided x. Same as other solve methods, but for floats[].

	 @param  A  x: native-style vector of x.length = A.rows()
	 @param  A  b: native-style vector of b.length = A.rows()

	 @throws IllegalArgumentException if b or a has the wrong size.
  */
  public float[] solve(float[] x, float[] b)
  {
	int i,j;
	double sum;
	double[] y,z;


	if (x==null) 
	  {
		x = new float[n];
	  }
	else if (x.length != n)
	  {
		throw new IllegalArgumentException("X vector has incorrect number of elements");
	  }
	if (b.length != n)
	  {
		throw new IllegalArgumentException("RHS vector has incorrect number of elements");
	  }

	y = new double[n];
	z = new double[n];

	//  
	//  Solve Ly=Pb
	//  
	y[0] = b[P[0]]/_A.getQuick(0,0);
	for (i=0 ; i<n ; i++) 
	  {
		sum=0.0;
		for (j=0 ; j<=i-1 ; j++) 
		  {
			sum += (_A.getQuick(i,j)*y[j]);
		  }
		y[i] = (b[P[i]] - sum) / _A.getQuick(i,i);
	  }
	//  
	//  Solve Ltranspose z = y
	//  
	z[n-1] = y[n-1]/_A.getQuick(n-1,n-1);
	for (i=n-2; i>=0 ; i--) 
	  {
		sum=0.0;
		for (j=i+1 ; j<n ; j++) 
		  {
			sum += _A.getQuick(j,i) * z[j];
		  }
		z[i] = (y[i] - sum) / _A.getQuick(i,i);
	  }
	//  
	//  x = Ptranspose z
	//  
	for (i=0 ; i<n ; i++) 
	  {
		x[P[i]] = (float)z[i];
	  }
	y=null;
	z=null;

	return x;
  }


 /**
	 Static method to returns the inverse of the 
	 original matrix. using the decomposition <p>
	 C=L.L^T, and then<p>
	 C^{-1} = (L^{-1})^T.L^{-1}<p>
	 internally.

	 Does not corrupt supplied matrix.

	 Typical calling sequence is
	 <pre>
	 DoubleMatrix2D Ainv = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(A);
	 </pre>

	 @param L square matrix
	 @throws IllegalArgumentException if <tt>L</tt> is not square.
  */

  public static DoubleMatrix2D getInverseViaCholesky(DoubleMatrix2D C)
  {
	int i,j,n;
	double tmp;
	DoubleMatrix2D A,L,Ccpy;
	ModifiedCholeskyLeftDecomposition chol;
	cern.colt.matrix.linalg.Blas blas = cern.colt.matrix.linalg.SeqBlas.seqBlas;

	n = C.rows();
	
	if (C.columns() != n) throw new IllegalArgumentException("matrix is not square");
	
	Ccpy = C.like().assign(C);
	// chol = new ModifiedCholeskyLeftDecomposition(Ccpy,true,true);
	chol = new ModifiedCholeskyLeftDecomposition(Ccpy);
	if (chol.wasModified()) 
	  {
		// System.out.println("Modified!");
	  }
	L = chol.getL();
	int[] p = chol.getPermutation();

	chol.getCholeskyLeftInverse(L);
	A = C.like();
	blas.dgemm(true,false,1.0,L,L,0,A);
	for (i=0;i<n;i++) {
	  for (j=0;j<n;j++) {
		tmp = A.getQuick(i,j);
		L.setQuick(p[i],p[j], tmp);
	  }
	}
	return L;
  }


  /**
	 Returns the inverse of the left Cholesky matrix L *in place*. Typical calling sequence is
	 <pre>
	 ModifiedCholeskyLeftDecomposition D = ModifiedCholeskyLeftDecomposition(A); 
	 DoubleMatrix2D R = D.getL();
	 D.getCholeskyLeftInverse(R);
	 </pre>

	 Can be called with any Cholesky factor L

	 @param L square lower triangular matrix (usually obtained from getL())
	 @throws IllegalArgumentException if <tt>L</tt> is not square.
  */

  public static void getCholeskyLeftInverse(DoubleMatrix2D L)
  {
	int i,j,k;
	double sum,tmp;
	int n = L.rows();
	
	if (L.columns() != n) throw new IllegalArgumentException("matrix is not square");
	
	// Main loop.
	if (true)
	  {
		for (i=0;i<n;i++) 
		  {
			tmp = L.getQuick(i,i);
			L.setQuick(i,i,(tmp!=0)? 1.0/tmp : 0.0);
			for (j=i+1; j<n; j++) 
			  {
				sum=0.0;
				for (k=i; k<j;k++) sum -= L.getQuick(j,k)*L.getQuick(k,i);
				tmp = L.getQuick(j,j);
				L.setQuick(j,i,(tmp!=0)? sum/tmp: 0.0);
			  }
		  }
	  }
	else
	  {
		double[] s = new double[n];
		// scale..
		for (i=0;i<n;i++) s[i]=Math.sqrt(L.getQuick(i,i));
		for (i=0;i<n;i++) 
		  {
			for (j=0; j<i; j++) 
			  {
				L.setQuick(i,j, L.getQuick(i,j)/(s[i] * s[j]));
			  }
			L.setQuick(i,i,1.0);
		  }
		// System.out.println("scaled L="+L.toString());
		// inverse in scaled form..
		for (i=0;i<n;i++) 
		  {
			tmp = L.getQuick(i,i);
			L.setQuick(i,i,(tmp!=0)? 1.0/tmp : 0.0);
			for (j=i+1; j<n; j++) 
			  {
				sum=0.0;
				for (k=i; k<j;k++) sum -= L.getQuick(j,k)*L.getQuick(k,i);
				tmp = L.getQuick(j,j);
				L.setQuick(j,i,(tmp!=0)? sum/tmp: 0.0);
			  }
		  }
		// unscale...
		for (i=0;i<n;i++) 
		  {
			for (j=0; j<=i; j++)  
			  {
				L.setQuick(i,j, L.getQuick(i,j)/(s[i]*s[j]));
			  }
		  }
		
	  }
  }
}
