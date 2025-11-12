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

package au.csiro.JamesGunning.lib.stat.clustering;

import java.io.*;
import java.util.*;

import cern.colt.map.*;
import cern.jet.random.*;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.stat.*;
import cern.colt.list.*;

import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;

 /**
	WARNING: This class is not fully debugged yet. It is not yet used in <i>delivery</i>
	<p>

	 Clustering algorithms based on Kmeans method. Kmeans algorithms
	 translated from algorithm AS 136 Appl. Statist. (1979) v.28,
	 Nn.1. Sourced from kmns.f as used in the GNU <i>R</i>
	 distribution
	 <p>
	 Recursive clustering to find optimal numbers of clusters taken from ideas 
	 in Ch. 7 (Bock) of	 "Clustering & Classification", Arabie, Hubert, De Soete.
	 <p>
	 Typical use requires m data points in n dimensions, here represented as
	 the<br> double[][] array   data[0....m-1][0....n-1]
	 <p>
	 Inital scaling of data (<i>in-place</i>) is highly recommended: use<br>

	 <pre>
	 KmeansClusterer.extremaScaleDataArray(data);</tt>
	 </pre>		
	 Extracting of number of clusters and their properties is typically done <br>
	 <pre>
	 KmeansClusterer autoclusterer = new KmeansClusterer();
	 autoclusterer.recursivelyClusterUsingKstar2n(data);
	 </pre>	
	 And then properties can be extracted using<br>
	 <pre>
	 getNClusters()<br>
	 getClusterLabels()<br>
	 getClusterCenters()<br>
	 getNumberOfPointsPerCluster()<br>
	 printOutClassification()<br>
	 </pre>
	 
 */

public class KmeansClusterer
{
 
  
  //  Code in this file based on Applied Statistics algorithms 
  //  (C) Royal Statistical Society 1979
  // 
  // 
  //      ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1
  // 
  //      Divide M points in N-dimensional space into K clusters so that
  //      the within cluster sum of squares is minimized.
  // 
  //   Translation from Fortran input file: kmns.f used in the R distribution
  // 
  
  private int[] ic1;
  private int[] ic2;
  private int[] nc;
  private int[] ncp;
  private int[] itran;
  private int[] live;
  private double[] an1;
  private double[] an2;
  private double[] d;
  private double[] wss;
  private int[] ifault;
  private double[] _centers;
  public static double BIG = Double.MAX_VALUE;

  public static boolean VERBOSELY_PRINT_FILE = false;

  private int branchCode,levelCode;
  private HashMap labelsMap;
  private int nClusters;

  private 	double[][] centers2D;

  private 	double[][] _data;

  private int verbosity=0;

  public KmeansClusterer()
  {
  }

  /** 
	  <p>
	  Sample driver with some fake data sowing some obvious clusters.
	  With no arguments, invents some fake data to cluster.
	  With argument \<filename\>, will cluster the data in the file (assumes
	  ascii, N data points, one per row, x,y,z, etc in columns)
	  <p>
	  
	  try, for example :<br>
	  <p>
	  java KmeansClusterer -i junk.txt | awk 'p==1{print}/Classification/{p=1}' | sort -n -k 1 | \<br>
	  awk 'BEGIN{nc=$1} $1!=nc{print "\n"}{nc=$1;print $2,$3}' | graph -TX -S 2 -m 0 <br>
	  <p>
	  See the example usage in the source code of main()
	  <p>
  */
  public static void main(String argv[])
  {
	int i,j,k=0,kk;
	int m,n=0;
	double val;
  	double[][] dataR, centersR;
	int numClusters;
	String filename=null;
	boolean doMonteCarloForKstar=false;


	// Parse basic requirements and extra options
	
	// handle errors: can change how these are handled if desired
	ErrorHandler EH = new ErrorHandler();
	OptionsParser options = new OptionsParser("KmeansClusterer",EH);
	options.addSpielText("Typical autoclustering usage:\njava au.csiro.JamesGunning.lib.stat.clustering.KmeansClusterer -i data.txt | awk 'p==1{print}/Classification/{p=1}' | sort -n -k 1 | awk 'BEGIN{nc=$1} $1!=nc{print \"\n\"}{nc=$1;print $2,$3}' | graph -TX -S 2 -m 0");
	options.addOption("-i",null,1,"String","file.txt","Ascii data file to auto-cluster");
	options.addOption("-v",null,1,"int","level","Execute with stated level of verbosity (to stderr)");
	options.addOption("-m",null,0,null,null,"Do Monte-Carlo studies for Kstar");

	options.parseOptions(argv,true);

	KmeansClusterer clusterer = new KmeansClusterer();

	if (options.isOptionSet("-v")) 
	  clusterer.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();

	if (options.isOptionSet("-i"))
	  {
		filename = ((String)options.getOptionValues("-i")[0]);
	  }
	if (options.isOptionSet("-m"))
	  {
		doMonteCarloForKstar=true;
	  }

	if (filename==null)
	  {
		// make up dummy clustered data & try to recover clustering
		
		m=400; // n points
		n=2;   // dim of pts
		dataR = new double[m][n];

		// choose the number of clusters
		if (doMonteCarloForKstar)
		  {
			// Monte Carlo database for modelling log(kstar) distribution as a function of dim and ndata
			MersenneTwister ranEngine;
			AbstractDistribution Ngen;
			Ngen = new Normal(0,1,new MersenneTwister(MersenneTwister.DEFAULT_SEED));
			//Ngen = new Uniform(-1,1,0);  // can affect the kstar calibration etc. see Boch's discussion
			double kstar=0;
			boolean test=false;
			int Ntrials=10000;
			System.out.println("dims ndata <log(kstar)> var(log(kstar))  (for "+Ntrials+" trials)");
			for (n=1;n<21;n++)  // loop on dimensions
			  {
				for (m=10;m<=1000;m*=2) // loop on ndata
				  {
					dataR = new double[m][n];
					DoubleArrayList samples = new DoubleArrayList();

					for (k=0;k<Ntrials;k++)
					  {
						do 
						  {
							// make N(0,1) data set
							for (i=0;i<m;i++)
							  {
								for (j=0;j<n;j++)
								  {
									dataR[i][j] = Ngen.nextDouble();
								  }
							  }
							//		for (int it=0;it<1000;it++)
							// {
							// cluster it 
							kstar = clusterer.doKmeansDualClustering(dataR);
							//	System.out.println("TT "+it+" "+kstar);
							// }
							// occasionally things don't work: try again.....
							test = (kstar==0) || Double.isNaN(kstar) || Double.isInfinite(kstar);
						  }
						while (test); 
						// System.out.println("Kstar = "+kstar);
						// collect samples
						samples.add(Math.log(kstar));
					  }
					// write out dimension, ndata, mean and variance
					double kbar = Descriptive.mean(samples);
					System.out.println(n+" "+m+" "+kbar+ " "+Descriptive.sampleVariance(samples,kbar));
				  }
			  }
		  }
		else
		  {
			numClusters=2;
			Random RG = new Random();
			for (i=0;i<m;i++)
			  {
				for (j=0;j<n;j++)
				  {
					dataR[i][j] = 0.7*RG.nextGaussian() + 0.25*i/(m/numClusters);
				  }
			  }

			int success = clusterer.doKmeansGeneralClustering(dataR,numClusters);

			if (clusterer.verbosity>0)
			  {
				switch (success)
				  {
				  case 1:
					System.out.println("Empty cluster, try better set of initial centers");
					break;
				  case 2:
					System.out.println("did not converge in iter.max iterations");
					break;
				  case 3:
					System.out.println("number of cluster centres must lie between 1 and m");
					break;
				  }
			  }
			clusterer.printOutClassification(true);
			clusterer.printOutMeans();
		  }
	  }
	else
	  {
		// read data from file given by argv[1]
		ArrayList data = new ArrayList();
		try 
		  {
			RandomAccessFile f = new RandomAccessFile(filename, "r");		
			String s;
			m=0;
			while ( (s=f.readLine()) != null )  
			  {
				m++;
				// actual data
				StringTokenizer st = new StringTokenizer(s);
				// set dimension
				if (m==1) n=st.countTokens();
				double[] row = new double[n];
				for (i=0;i<n;i++)
				  {
					try
					  {
						val = Double.parseDouble(st.nextToken());
					  }
					catch (NumberFormatException e)
					  {
						throw new RuntimeException("Line "+m+" of file " + filename+" has non-numbers");
					  }
					// add data to lists
					row[i] = val;
				  } 
				data.add(row);
			  }
		  }
		catch (Exception e) 
		  {
			throw new RuntimeException("Can't open file "+filename);
		  }
		// reorganise into a double[][] array
		dataR = new double[m][n];
		for (i=0;i<m;i++)
		  {
			for (j=0;j<n;j++)
			  {
				dataR[i][j] = ((double[])data.get(i))[j];
			  }
		  }
		// scale it
		KmeansClusterer.extremaScaleDataArray(dataR);
		clusterer.recursivelyClusterUsingKstar2n(dataR);
		clusterer.printOutClassification(false);
	  }
  }

  /**
	 Scale data in-place using z_5 scaling of Milligan etal section 4.3 of Arabie et al
     (This is a linear mapping to [0,1] based on extrema in each dimension)
  */

  public static RescalingInfo extremaScaleDataArray(double[][] data)
  {
	RescalingInfo rescalingInfo = new RescalingInfo(data);
	rescalingInfo.rescale();
	return rescalingInfo;
  }


  private static double[] getPointsFromDataArray(double[][] data)
  {
	// assume rectangular data
	int i,j,k;
	int n = data[0].length; 
	int m = data.length;
	double[] points = new double[m*n];

	for (i=0,k=0;i<n;i++)
	  {
		for (j=0;j<m;j++,k++)
		  {
			points[k] = data[j][i];
		  }
	  }
	return points;
  }

  private static double[] guessCentersFromDataArray(int numClusters, double[][] data, boolean addNoise)
  {
	int i,j;
	// make centers equidistant points spanning the data range
	int n = data[0].length; 
	int m = data.length;
	double[][] centersR = new double[numClusters][n];
	for (i=0;i<n;i++)
	  {
		double max=-Double.MAX_VALUE;
		double min=Double.MAX_VALUE;
		for (j=0;j<m;j++)
		  {
			max=Math.max(data[j][i],max);
			min=Math.min(data[j][i],min);
		  }
		// for this dimension i, divide into n segments
		for (j=0;j<numClusters;j++)
		  {
			centersR[j][i] = min+j*(max-min)/((numClusters>1)? (numClusters-1) : 1);
			// silly noise to test robustness wrt starting points
			if (addNoise) centersR[j][i] += 0.05*(max-min)*Math.random(); 
		  }
	  }
	return getCentersFromCentersArray(centersR);
  }

 private static double[] getCentersFromCentersArray(double[][] centers)
  {
	return getPointsFromDataArray(centers);
  }

  /**
	 Cluster m points in n dimensions data[0....m-1][0....n-1] into nClusters clusters<br>
	 @return <tt>3</tt> too few points<br>
	 <tt>2</tt> Something didn't work.<br>
	 <tt> 1</tt> Auto-guessed cluster centers yielded an empty cluster<br>
	 <tt>otherwise</tt> Things went OK.
  */
  
  public int doKmeansGeneralClustering(double[][] data, int numClusters)
  {
	int success;

	_data = data;
	nClusters = numClusters;
	int n = data[0].length; 
	int m = data.length;

	double[] points = getPointsFromDataArray(data);
	double[] centers = guessCentersFromDataArray(numClusters,data,true);

	success = doKmeansClusteringFortranStyle(n,m,points,centers,numClusters);

	compute2Dcenters(data);

	return success;
  }

  // split points into 2 clusters and get Kstar2n measure

  private double doKmeansDualClustering(double[][] data)
  {
	int i,ii,j,k;
	int n = data[0].length; 
	int m = data.length;
	double Kstar2n;
	double[] points;
	double[] centers;

	Kstar2n = 0.0;

	int success = doKmeansGeneralClustering(data,2);

	switch (success)
	  {
	  case 3:
		// too few points. Fake cluster termination by zero Kstar
		return Kstar2n;
	  case 2:
		// Something didn't work, so leave as is. Same trick as above for now
		return Kstar2n;
	  case 1:
		// Empty cluster so leave as is. Same trick as above for now
		return Kstar2n;
	  default:
		// Things are probably ok. keep going
	  }
	if (verbosity>0) System.out.println("Clustering success = "+success);

	// form Kstar2n statistic

	centers2D = get2Darrayfrom1Darray(n,2,_centers);
	// form mean
	double[] lumpedMean = new double[n];
	for (j=0;j<n;j++)
	  {
		lumpedMean[j]=0.0;
		for (k=0;k<m;k++) lumpedMean[j] += data[k][j];
		lumpedMean[j]/=(double)m;
	  }

	double norm,num=0.0,den=0.0;
	for (i=0,ii=1;i<2;i++,ii++)
	  {
		for (j=0,norm=0.0;j<n;j++)
		  {
			norm += (centers2D[i][j]-lumpedMean[j])*(centers2D[i][j]-lumpedMean[j]);
		  }
		num += getNumberOfPointsPerCluster()[i] * norm;
		int[] labels = getClusterLabels();
		norm=0;
		for (k=0;k<labels.length;k++)
		  {
			if (labels[k]==ii)
			  {
				// point belongs to ith cluster. Add to norm (sum norm over dimensions for this point k)
				for (j=0;j<n;j++)
				  {
					norm += (centers2D[i][j]-data[k][j])*(centers2D[i][j]-data[k][j]);
				  }
			  }
		  }
		den += norm;
	  }
	Kstar2n = num/den;

	return Kstar2n;	
  }

  /**
	 Recursive cluster-splitting algorithm to find optimal number of
	 clusters using the k<sup>*</sup><sub>2n</sub> statistic computed
	 on split clusters, where the clusters are generated using the
	 Kmeans algorithm, as discussed in Bock, p 421 (section 7.1) of
	 "Clustering and Classification", Arabie et al.  We start from 1
	 cluster (all the data) and recursively split each cluster into 2
	 clusters, testing for significance using
	 k<sup>*</sup><sub>2n</sub> > k<sub>1-<FONT FACE="Symbol">a</FONT></sub>, where 
    <FONT FACE="Symbol">a</FONT>=0.05 is the level of significance
	 used. The threshold values of k<sub>1-<FONT
	 FACE="Symbol">a</FONT></sub> are from Engelman and Hartigan,
	 referenced in Bock's chapter of Arabie et al, p 422.
  */
  public void recursivelyClusterUsingKstar2n(double data[][])
  {
	int i,j,n,label,newLabel;
	
	_data = data;
	branchCode=0;
	levelCode=-1;

	// labels map helps keep track of absolute labels as we chop up the
	//data using recursive cluster-splitting
	labelsMap = new HashMap(data.length);
	for (i=0;i<data.length;i++) labelsMap.put(data[i],new MutableInteger(0));
	 
	// ok now cluster
	recursivelyCluster(data);


	// Patching
	// (1) FIX LABELS
	// remap labels to be sane sequence 1..2..3..
    OpenIntIntHashMap uniqueLabels = new OpenIntIntHashMap();
	n = data[0].length;
	newLabel=1;
	ic1 = new int[data.length];

	for (i=0;i<data.length;i++) 
	  {
		label = ((MutableInteger)labelsMap.get(data[i])).value;
		if (!uniqueLabels.containsKey(label))
		  {
			uniqueLabels.put(label,newLabel);
			newLabel++;
		  }
		// rewrite labels
		ic1[i] = uniqueLabels.get(label);
	  }

	// (2) and update number of points-per-cluster, centers
	nClusters = newLabel-1;
	nc = new int[nClusters];
	for (i=0;i<data.length;i++)
	  {
		nc[ic1[i]-1]++;
	  }
	  
	compute2Dcenters(data);

	_data = data;
  }
  /**
	 A method that prints out classification and data to stdout <i>after</i> classification has been done.
	 Produces nm rows of ascii columns, with (n+1) columns, with first column being the cluster label.
  */
  public void printOutClassification(boolean sort)
  {
	int i,j;
	int n = _data[0].length;
	int m = _data.length;

	int labels[] = getClusterLabels();

	System.out.println("#Classification");
	if (!sort)
	  {
		for (i=0;i<m;i++)
		  {
			System.out.print(labels[i]);
			for (j=0;j<n;j++) System.out.print(" " + _data[i][j]);
			System.out.println("");
		  }
	  }
	else
	  {
		for(int l=1;l<=nClusters;l++)
		  {
			for (i=0;i<m;i++)
			  {
				if (labels[i]!=l) continue;
				System.out.print(labels[i]);
				for (j=0;j<n;j++) System.out.print(" " + _data[i][j]);
				System.out.println("");
			  }
			System.out.println("");
		  }
	  }
  }
  /**
	 A method that prints out cluster centers to stdout <i>after</i> classification has been done.
	 Produces nClusters rows of ascii columns, with n columns, showing mean vector for clusters.
  */
  public void printOutMeans()
  {
	int i,j;
	int n = _data[0].length;
	double[][] centers2D = compute2Dcenters(_data);
	int nClusters = centers2D.length;
	System.out.println("Centers");
	for (i=0;i<nClusters;i++)
	  {
		for (j=0;j<n;j++) System.out.print(" " + centers2D[i][j]);
		System.out.printf("%n%n");
	  }
	
  }

  private double[][] compute2Dcenters(double[][] data)
  {
	int i,j;
	// 
	//      Update cluster centres to be the average of points contained
	//      within them.
	// 
	int n = data[0].length;

	double[][] centers2D = new double[nClusters][n];
	for (i=0;i<data.length;i++) 
	  {
		for (j=0;j<n;j++)
		  {
			centers2D[ic1[i]-1][j] += data[i][j];
		  }
	  }
	for (i=0;i<nClusters;i++) 
	  {
		for (j=0;j<n;j++)
		  {
			centers2D[i][j]/=(double)nc[i];
		  }
	  }
	/*
	System.out.println("Centers");
	for (i=0;i<nClusters;i++)
	  {
	  for (j=0;j<n;j++)
	  {
	  System.out.print(" " + centers2D[i][j]);
		  }
		  System.out.println("");
		  }
	*/
	return centers2D;
  }
  /**
	 Threshold function for Kstar, develop from Montecarlo samples
	 of Kstar from N(0,1) data sets (see main() with -m option).
	 The distribution of log(Kstar) has been fitted with a normal distribution,
	 and for dimensions from 1 to MAX_DIM, the parameters of the normal
	 distribution come from an empirical power law regression.
	 See Bock, top of p 421 (sec. 7.1) of 
	 "Clustering & Classification", Arabie et al. Valid for arbitrary
	 dimensionality

	 @param ndata number of data
	 @param d dimensionality
	 @param alpha in (0,1) degree of confidence in cluster split
  */
  private static double KstarThreshold(int ndata, int d, double alpha)
  {
	double mu,sigma,quantile;
	int MAX_DIM=20;
	// These regression coefficients have been developed in mathematica using NonlinearFit[]
	double[] A_mu={0.565387, -0.731523, -1.28076, -1.64025, -1.90543, -2.11811, -2.29288,
				   -2.44643, -2.58419, -2.7024, -2.80825, -2.90732, -2.99633, -3.08285,
				   -3.15865, -3.2419, -3.3052, -3.38032, -3.43925, -3.49819};
	double[] B_mu={6.36497, 4.361, 4.27349, 4.28706, 4.46836, 4.55634, 4.6789, 4.77622,
				   4.82824, 4.94902, 5.04309, 5.11362, 5.18458, 5.23831, 5.32077, 5.32406,
				   5.41405, 5.43931, 5.50376, 5.54168};
	double[] p_mu={1.11383, 0.727107, 0.635334, 0.58091, 0.5528, 0.52859, 0.512309, 0.496514,
				   0.480607, 0.471683, 0.463622, 0.454792, 0.447857, 0.439681, 0.435049,
				   0.423779, 0.421399, 0.412311, 0.409174, 0.404484};
	double[] A_sigma={-0.00604128, -0.0212603, -0.0286518, -0.0344219, -0.0212596, -0.0297688,
					  -0.031449, -0.0358585, -0.0247034, -0.0366897, -0.0345655, -0.0383664,
					  -0.0369979, -0.0431549, -0.0471914, -0.0465607, -0.0476248, -0.0369593,
					  -0.0364509, -0.0446685};
	double[] B_sigma={1.3756, 0.972581, 0.787681, 0.679559, 0.638031, 0.564507, 0.53122, 0.494283,
					  0.479962, 0.44797, 0.433035, 0.410528, 0.39223, 0.381831, 0.372967, 0.363993,
					  0.350869, 0.342943, 0.332395, 0.328538};
	double[] p_sigma={0.476829, 0.381977, 0.344788, 0.316932, 0.334823, 0.302389, 0.292581,
					  0.275909, 0.293101, 0.26104, 0.261331, 0.246913, 0.243515, 0.228702, 0.21916,
					  0.217616, 0.210884, 0.226878, 0.222654, 0.207285};
	if (ndata<=2) 
	  {
		return 0.5*Double.MAX_VALUE;
	  }
	else
	  {
		// handle large dimensionality case in clunky way
		if (d>MAX_DIM) d=MAX_DIM;
		// parameters of Normal approximation to distribution of log(k_star)
		// Empirical kstar obserations (from Monte-Carlo) 
		// have been fitted with power laws of this form
		mu = A_mu[d-1]+B_mu[d-1]*Math.pow(ndata,-p_mu[d-1]);
		sigma = A_sigma[d-1]+B_sigma[d-1]*Math.pow(ndata,-p_sigma[d-1]);
		// corresponding quantile for alpha
		quantile = mu + cern.jet.stat.Probability.normalInverse(1.0-alpha)*sigma;
		// convert back to kstar
		return Math.exp(quantile);
	  }
  }
  /**
	 Threshold function from Engelman and Hartigan (1969), 
	 discussed in Bock, top of p 421 (sec. 7.1) of 
	 "Clustering & Classification", Arabie et al. Valid only for 1D
  */
  private static double KstarThreshold1D(int ndata, double alpha)
  {
	double logterm;

	if (ndata<=2) 
	  {
		return 0.5*Double.MAX_VALUE;
	  }
	else
	  {
		logterm = (-Math.log(1.0-2.0/Math.PI)+
		  cern.jet.stat.Probability.normalInverse(1.0-alpha)/Math.sqrt(ndata-2)+2.4/(ndata-2));
		return Math.exp(logterm)-1.0;
	  }
  }

  private void recursivelyCluster(double data[][])
  {
	int i,j,k,kk;
	double[][] cdata;
	double Kstar, KstarThresh;
	int n = data[0].length; 
	int m = data.length;
	int lbl;

	levelCode++;


	RescalingInfo rescalingInfo = new RescalingInfo(data);
	rescalingInfo.rescale();
	// Get Kstar
	Kstar = doKmeansDualClustering(data);
	// undo scaling before it propogates
	rescalingInfo.restoreData();

	// use 95% threshold
	KstarThresh = KstarThreshold(m,n,0.05);
	if (verbosity>0) System.out.println("Splitting: ndata= "+m+" in dimensions="+n+" returned kstar = "+Kstar+ " threshold = "+ KstarThresh);
	if (Kstar>KstarThresh)
	  {
		// It is statistically significant to split this cluster
		// System.out.println("Results for attempted cluster on branch "+branchCode+" level " +levelCode);
		// System.out.println("Kstar = "+Kstar);

		int[] labels = (int[])(getClusterLabels().clone());
		int[] NumberOfPointsPerCluster = (int[])getNumberOfPointsPerCluster().clone();

		// update global labels
		lbl = 2*branchCode + 1;
		for (k=0;k<labels.length;k++)
		  {
			((MutableInteger)labelsMap.get(data[k])).value = labels[k]+lbl;
		  }

		if (VERBOSELY_PRINT_FILE || verbosity>1) printToFileDualClusterData("cluster-split-data."+branchCode+".txt",data);

		branchCode++;

		// try to split cluster 1 then cluster 2, recursively
		for (int cNumber=1;cNumber<=2;cNumber++)
		  {
			// System.out.println("Trying to subdivide cluster "+cNumber+", level "+levelCode);
			if (NumberOfPointsPerCluster[cNumber-1]>1)
			  {
				cdata = new double[NumberOfPointsPerCluster[cNumber-1]][n];
				for (k=0,kk=0;k<labels.length;k++)
				  {
					if (labels[k]==cNumber)	
					  {
						cdata[kk] = data[k];
						kk++;
					  }
				  }
				recursivelyCluster(cdata);
			  }
		  }
	  }
	levelCode--;
  }

  private double[][] get2Darrayfrom1Darray(int n,int m, double[] x)
  {
	// assume rectangular data
	int i,j,k;
	double[][] x2d = new double[m][n];
	for (i=0,k=0;i<n;i++)
	  {
		for (j=0;j<m;j++,k++)
		  {
			x2d[j][i] = x[k];
		  }
	  }
	return x2d;
  }

  /**
	 Write out dual split cluster data
  */
  private void printToFileDualClusterData(String name,double[][] data)
  {
	int i,k,n;
	n = data[0].length;

	try
	  {
		File f = new File(name);
		PrintWriter out = new PrintWriter(new FileWriter(f));

		int[] labels = (int[])(getClusterLabels().clone());
		int[] NumberOfPointsPerCluster = (int[])getNumberOfPointsPerCluster().clone();
		
		// write out cluster 1, then cluster 2
		for (int cNumber=1; cNumber<=2; cNumber++)
		  {
			for (k=0;k<labels.length;k++)
			  {
				if (labels[k]==cNumber)	
				  {
					out.print("L="+((MutableInteger)labelsMap.get(data[k])).value+" ");
					for (i=0;i<n;i++)
					  {  
						out.print(data[k][i]+" ");
					  }
					out.println("");
				  }
			  }
			out.println("");
		  }
		out.close();
		System.out.println("Wrote "+name);
	  }
	catch (IOException e) {System.out.println("Couldn't open " + name);};
  }

  /**
	 Error codes on return:
	 1: Empty cluster;
	 2: If the specified number of iterations has been exceeded.
	 3: If there are too few points to cluster into nClusters clusters.
  */

  private int doKmeansClusteringFortranStyle(int n, int m, double[] points, double[] centers, int numClusters)
  {
	int i,j,k=0;

	_centers = centers;
	ic1 = new int[m];
	ic2 = new int[m];
	nc=new int[numClusters];
	ncp = new int[numClusters];
	itran = new int[numClusters];
	live = new int[numClusters];
	an1 = new double[numClusters];
	an2 = new double[numClusters];
	d = new double[m];
	wss = new double[numClusters];
	ifault = new int[1];

	int iterMax=20;


	// On return;
	// ic1 = cluster labels
	// centers = cluster centers
	// wss = within cluster sum-of-squares
	// nc = # points in each cluster

	// do Kmeans clustering
	KmeansClusterer.kmns(points,m,n,centers,numClusters,ic1,ic2,nc,an1,an2,ncp,d,itran,live,iterMax,wss,ifault);

	return ifault[0];
  }

  /**
	 Get number of clusters. Call <i>after</i> clustering is done
	 @return number of clusters automatically chosen.
  */
  public int getNClusters()
  {
	return nClusters;
  }

  /**
	 Get cluster labels. Call <i>after</i> clustering is done
	 @return array of cluster labels [1...nClusters].
  */
  public int[] getClusterLabels()
  {
	return ic1;
  }

  /**
	 Get cluster centers[0...nClusters-1][0...n]. Call <i>after</i> clustering is done
	 @return array of cluster centers.
  */
  public double[][] getClusterCenters()
  {
	return centers2D;
  }
  
  /**
   Get number of points per cluster [0...nCluster-1]. Call <i>after</i> clustering is done
   @return array of cluster counts [1...nClusters].
  */
  public int[] getNumberOfPointsPerCluster()
  {
	return nc;
  }
  /**
	 Error codes on return:
	 ifault[0]=1: Empty cluster;
	 ifault[0]=2: If the specified number of iterations has been exceeded.
	 ifault[0]=3: If there are too few points to cluster into k clusters.
  */
  private static void kmns (double [] a,
						   int m,
						   int n,
						   double [] c,
						   int k,
						   int [] ic1,
						   int [] ic2,
						   int [] nc,
						   double [] an1,
						   double [] an2,
						   int [] ncp,
						   double [] d,
						   int [] itran,
						   int [] live,
						   int iter,
						   double [] wss,
						   int[] ifault)  
  {

	double [] dt= new double[(2)];
	int i= 0;
	int il= 0;
	int j= 0;
	int l= 0;
	int[] indx = new int[1];
	int ij= 0;
	int ii= 0;
	double da= 0.0;
	double temp= 0.0;
	double db= 0.0;
	double dc= 0.0;
	double aa= 0.0;
	ifault[0] = 3;

	// cannot have number of points less than number of clusters
	if (k <= 1 || k >= m)  return;

	ifault[0] = 0;
	boolean test;
	// 
	//      For each point I, find its two closest centres, IC1(I) and
	//      IC2(I).     Assign it to IC1(I).
	// 
	for (i = 1; i <= m; i++) 
	  {
		ic1[(i-1)] = 1;
		ic2[(i-1)] = 2;
		for (il = 1; il <= 2; il++) 
		  {
			dt[(il-1)] = 0.0;
			for (j = 1; j <= n; j++) 
			  {
				da = a[(i-1)+(j-1) * m]-c[(il-1)+(j-1) * k];
				dt[(il-1)] += da*da;
			  }
		  }
		if (dt[0] > dt[1])  
		  {
			ic1[(i-1)] = 2;
			ic2[(i-1)] = 1;
			temp = dt[0];
			dt[0] = dt[1];
			dt[1] = temp;
		  }
		  
		for (l = 3; l <= k; l++) 
		  {
			db = 0.0;
			test=false;
			
			for (j = 1; j <= n; j++) 
			  {
				dc = a[(i-1)+(j-1) * m]-c[(l-1)+(j-1) * k];
				db += dc*dc;
				if (db >= dt[1])  
				  {
					test=true;
					break;
				  }
			  }
			if (test) break;
			
			if (db < dt[0])  
			  {
				dt[1] = dt[0];
				ic2[(i-1)] = ic1[(i-1)];
				dt[0] = db;
				ic1[(i-1)] = l;
			  }
			else
			  {
				dt[1] = db;
				ic2[(i-1)] = l;
			  }
		  }
	  }
	// 
	//      Update cluster centres to be the average of points contained
	//      within them.
	// 
	for (l = 1; l <= k; l++) 
	  {
		nc[(l-1)] = 0;
		for (j = 1; j <= n; j++) 
		  {
			c[(l-1)+(j-1) * k] = 0.0;
		  }
	  }
	for (i = 1; i <= m; i++) 
	  {
		l = ic1[(i-1)];
		nc[(l-1)]++;
		
		for (j = 1; j <= n; j++) 
		  {
			c[(l-1)+(j-1) * k] += a[(i-1)+(j-1) * m];
		  }
	  }
	// 
	//      Check to see if there is any empty cluster at this stage
	// 
	for (l = 1; l <= k; l++) 
	  {
		if (nc[(l-1)] == 0)  
		  {
			ifault[0] = 1;
			return;
		  }
		aa = (double)(nc[(l-1)]);
		for (j = 1; j <= n; j++) 
		  {
			c[(l-1)+(j-1) * k] /= aa;
		  }
		// 
		//      Initialize AN1, AN2, ITRAN & NCP
		//      AN1(L) = NC(L) / (NC(L) - 1)
		//      AN2(L) = NC(L) / (NC(L) + 1)
		//      ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,

		//               = 0 otherwise
		//      In the optimal-transfer stage, NCP(L) stores the step at which
		//      cluster L is last updated.
		//      In the quick-transfer stage, NCP(L) stores the step at which
		//      cluster L is last updated plus M.
		// 
		an2[(l-1)] = aa/(aa+1.0);
		an1[(l-1)] = BIG;
		if (aa > 1.0)  
		  an1[(l-1)] = aa/(aa-1.0);
		itran[(l-1)] = 1;
		ncp[(l-1)] = -1;
	  }
	
	indx[0] = 0;
	test=false;
	
	for (ij = 1; ij <= iter; ij++) 
	  {
		// 
		//      In this stage, there is only one pass through the data.   Each
		//      point is re-allocated, if necessary, to the cluster that will
		//      induce the maximum reduction in within-cluster sum of squares.
		// 
		optra(a,m,n,c,k,ic1,ic2,
					nc,an1,an2,ncp,d,
					itran,live,indx);
		// 
		//      Stop if no transfer took place in the last M optimal transfer
		//      steps.
		// 
		if (indx[0] == m)  
		  {
			test=true;
			break;
		  }
		// 
		//      Each point is tested in turn to see if it should be re-allocated
		//      to the cluster to which it is most likely to be transferred,
		//      IC2(I), from its present cluster, IC1(I).   Loop through the
		//      data until no further change is to take place.
		// 
		qtran(a,m,n,c,k,ic1,ic2,
					nc,an1,an2,ncp,
					d,itran,indx);
		// 
		//      If there are only two clusters, there is no need to re-enter the
		//      optimal transfer stage.
		// 
		if (k == 2)  
		  {
			test=true;
			break;
		  }
		// 
		//      NCP has to be set to 0 before entering OPTRA.
		// 
		for (l = 1; l <= k; l++) 
		  {
			ncp[(l-1)] = 0;
		  }
	  }
	
	// 
	//      If the specified number of iterations has been exceeded, set
	//      IFAULT = 2.   This may indicate unforeseen looping.
	// 
	if (!test) ifault[0] = 2;
	// 
	//      Compute within-cluster sum of squares for each cluster.
	// 
	for (l = 1; l <= k; l++) 
	  {
		wss[(l-1)] = 0.0;
		for (j = 1; j <= n; j++) 
		  {
			c[(l-1)+(j-1) * k] = 0.0;
		  }
	  }
	for (i = 1; i <= m; i++) 
	  {
		ii = ic1[(i-1)];
		for (j = 1; j <= n; j++) 
		  {
			c[(ii-1)+(j-1) * k] += a[(i-1)+(j-1) * m];
		  }
	  }
	for (j = 1; j <= n; j++) 
	  {
		for (l = 1; l <= k; l++) 
		  {
			c[(l-1)+(j-1) * k] = c[(l-1)+(j-1) * k]/(float)(nc[(l-1)]);
		  }
		for (i = 1; i <= m; i++) 
		  {
			ii = ic1[(i-1)];
			da = a[(i-1)+(j-1) * m] - c[(ii-1)+(j-1) * k];
			wss[(ii-1)] += da*da;
		  }
	  }
	// 
	return;
  }

 // 
  //      If cluster L is updated in the last quick-transfer stage, it
  //      belongs to the live set throughout this stage.   Otherwise, at
  //      each step, it is not in the live set if it has not been updated
  //      in the last M optimal transfer steps.
  // 

  private static void optra (double [] a,
					  int m,
					  int n,
					  double [] c,
					  int k,
					  int [] ic1,
					  int [] ic2,
					  int [] nc,
					  double [] an1,
					  double [] an2,
					  int [] ncp,
					  double [] d,
					  int [] itran,
					  int [] live,
					  int[] indx)  
  {

	int l= 0;
	int i= 0;
	int l1= 0;
	int l2= 0;
	int ll= 0;
	int j= 0;
	double de= 0.0;
	double df= 0.0;
	double da= 0.0;
	double db= 0.0;
	double r2= 0.0;
	double rr= 0.0;
	double dc= 0.0;
	double dd= 0.0;
	double al1= 0.0;
	double alw= 0.0;
	double al2= 0.0;
	double alt= 0.0;

	boolean test;

	for (l = 1; l <= k; l++) 
	  {
		if (itran[(l-1)] == 1)  
		  live[(l-1)] = m+1;
	  }


	for (i = 1; i <= m; i++) 
	  {
		indx[0] += 1;
		l1 = ic1[(i-1)];
		l2 = ic2[(i-1)];
		ll = l2;
		// 
		//      If point I is the only member of cluster L1, no transfer.
		// 
		if (nc[(l1-1)] == 1)  
		  {
			//		  Dummy.go_to("Optra",90); /// HERE
			if (indx[0] == m)  return;
			break;
		  }
		// 
		//      If L1 has not yet been updated in this stage, no need to
		//      re-compute D(I).
		// 
		if (ncp[(l1-1)] != 0)  
		  {
			de = 0.0;
			for (j = 1; j <= n; j++) 
			  {
				df = a[(i-1)+(j-1) * m]-c[(l1-1)+(j-1) * k];
				de += df*df;
			  }
			d[(i-1)] = de*an1[(l1-1)];
		  }
		// 
		//      Find the cluster with minimum R2.
		// 
		da = 0.0;
		for (j = 1; j <= n; j++) 
		  {
			db = a[(i-1)+(j-1) * m]-c[(l2-1)+(j-1) * k];
			da += db*db;
		  } 

		r2 = da*an2[(l2-1)];
		for (l = 1; l <= k; l++) 
		  {
			// 
			//      If I >= LIVE(L1), then L1 is not in the live set.   If this is
			//      true, we only need to consider clusters that are in the live set
			//      for possible transfer of point I.   Otherwise, we need to consider
			//      all possible clusters.
			// 
			if (!(i >= live[(l1-1)] && i >= live[(l-1)] || l == l1 || l == ll))
			  {
				rr = r2/an2[(l-1)];
				dc = 0.0;
				
				test=false;
				for (j = 1; j <= n; j++) 
				  {
					dd = a[(i-1)+(j-1) * m]-c[(l-1)+(j-1) * k];
					dc += dd*dd;
					if (dc >= rr)  
					  {
						test=true;
						break;
					  }
				  }
				if (test) break;
				r2 = dc*an2[(l-1)];
				l2 = l;
			  }
		  }

		test=false;
		if (!(r2 < d[(i-1)]))
		  {
			// 
			//      If no transfer is necessary, L2 is the new IC2(I).
			// 
			ic2[(i-1)] = l2;
		  }
		else
		  {
			// 
			//      Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
			//      L2, and update IC1(I) & IC2(I).
			// 
			indx[0] = 0;
			live[(l1-1)] = m+i;
			live[(l2-1)] = m+i;
			ncp[(l1-1)] = i;
			ncp[(l2-1)] = i;
			al1 = (double)(nc[(l1-1)]);
			alw = al1-1.0;
			al2 = (double)(nc[(l2-1)]);
			alt = al2+1.0;
			for (j = 1; j <= n; j++) 
			  {
				c[(l1-1)+(j-1) * k] = (c[(l1-1)+(j-1) * k]*al1-a[(i-1)+(j-1) * m])/alw;
				c[(l2-1)+(j-1) * k] = (c[(l2-1)+(j-1) * k]*al2+a[(i-1)+(j-1) * m])/alt;
			  }

			nc[(l1-1)] = nc[(l1-1)]-1;
			nc[(l2-1)] = nc[(l2-1)]+1;
			an2[(l1-1)] = alw/al1;
			an1[(l1-1)] = BIG;
			if (alw > 1.0)  
			  an1[(l1-1)] = alw/(alw-1.0);
			an1[(l2-1)] = alt/al2;
			an2[(l2-1)] = alt/(alt+1.0);
			ic1[(i-1)] = l2;
			ic2[(i-1)] = l1;
		  }
		if (indx[0] == m)  return;
	  }



	for (l = 1; l <= k; l++) 
	  {
		// 
		//      ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
		//      decreased by M before re-entering OPTRA.
		// 
		itran[(l-1)] = 0;
		live[(l-1)] -= m;
	  }
	return;
  }

 private static void qtran (double [] a,
							int m,
							int n,
							double [] c,
							int k,
							int [] ic1,
							int [] ic2,
							int [] nc,
							double [] an1,
							double [] an2,
							int [] ncp,
							double [] d,
							int [] itran,
							int[] indx)  
  {

	int icoun= 0;
	int istep= 0;
	int i= 0;
	int l1= 0;
	int l2= 0;
	int j= 0;
	double da= 0.0;
	double db= 0.0;
	double dd= 0.0;
	double al1= 0.0;
	double alw= 0.0;
	double al2= 0.0;
	double alt= 0.0;
	double r2= 0.0;
	double de= 0.0;
	icoun = 0;
	istep = 0;

	boolean test;

	for(;;)
	  {
		for (i = 1; i <= m; i++) 
		  {
			icoun++;
			istep++;
			l1 = ic1[(i-1)];
			l2 = ic2[(i-1)];
			// 
			//      If point I is the only member of cluster L1, no transfer.
			// 
			if (nc[(l1-1)] != 1)
			  {
				// 
				//      If ISTEP > NCP(L1), no need to re-compute distance from point I to
				//      cluster L1.   Note that if cluster L1 is last updated exactly M
				//      steps ago, we still need to compute the distance from point I to
				//      cluster L1.
				// 
				if (istep <= ncp[(l1-1)])
				  {
					da = 0.0;
					for (j = 1; j <= n; j++) 
					  {
						db = a[(i-1)+(j-1) * m]-c[(l1-1)+(j-1) * k];
						da += db*db;
					  }
					d[(i-1)] = da*an1[(l1-1)];
				  }
				// 
				//      If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
				//      point I at this step.
				// 
				if (!(istep >= ncp[(l1-1)] && istep >= ncp[(l2-1)]))
				  {
					r2 = d[(i-1)]/an2[(l2-1)];
					dd = 0.0;
					test=false;
					for (j = 1; j <= n; j++) 
					  {
						de = a[(i-1)+(j-1) * m]-c[(l2-1)+(j-1) * k];
						dd += de*de;
						if (dd >= r2)
						  {
							test=true;
							break;
						  }			  
					  }
					if (!test) 
					  {
						// 
						//      Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
						//      L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
						//      updating occurs in this stage, INDX is set back to 0.
						// 
						icoun = 0;
						indx[0] = 0;
						itran[(l1-1)] = 1;
						itran[(l2-1)] = 1;
						ncp[(l1-1)] = istep+m;
						ncp[(l2-1)] = istep+m;
						al1 = (double)(nc[(l1-1)]);
						alw = al1-1.0;
						al2 = (double)(nc[(l2-1)]);
						alt = al2+1.0;
						for (j = 1; j <= n; j++) 
						  {
							c[(l1-1)+(j-1) * k] = (c[(l1-1)+(j-1) * k]*al1-a[(i-1)+(j-1) * m])/alw;
							c[(l2-1)+(j-1) * k] = (c[(l2-1)+(j-1) * k]*al2+a[(i-1)+(j-1) * m])/alt;
						  }
					
						nc[(l1-1)] = nc[(l1-1)]-1;
						nc[(l2-1)] = nc[(l2-1)]+1;
						an2[(l1-1)] = alw/al1;
						an1[(l1-1)] = BIG;
						if (alw > 1.0)  an1[(l1-1)] = alw/(alw-1.0);
						an1[(l2-1)] = alt/al2;
						an2[(l2-1)] = alt/(alt+1.0);
						ic1[(i-1)] = l2;
						ic2[(i-1)] = l1;
					  }
				  }
			  }
			// 
			//      If no re-allocation took place in the last M steps, return.
			//
			if (icoun == m) return;  // 60
		  }
	  }
  }

  /**
	   Class that rescales a rectangular array data[0...m-1][0...n-1] 
	   of m points in n dimensions
	   into the interval [0,1] in each dimension.
	   The <i>in-place</i> rescale for each data point <tt>j</tt> is
	   <pre>
	   data[j][i] <-  (data[j][i] - minimums[i])/(maximums[i]-minimums[i])
	   </pre>
	   from which you can recover the original data by using the publicly accessible
	   minimums[] and maximums[] arrays.
	*/
  public static class RescalingInfo
  {
	/**
	   Array of data minimums for each dimension
	*/
	public double[] minimums;
	/**
	   Array of data maximums for each dimension
	*/
	public double[] maximums;	

	private int n,m;

	private double[][] _data;

	/**
	   Set up a Rescaler for the rectangular array data[0...m-1][0...n-1].
	   (Set up the minimums[] and maximums[] arrays)
	*/

	public RescalingInfo(double[][] data)
	{
	  int i,j;

	  _data = data;

	  // assume rectangular data
	  n = data[0].length; 
	  m = data.length;

	  minimums = new double[n];
	  maximums = new double[n];

	  // set up max and mins.
	  for (i=0;i<n;i++)
		{
		  maximums[i]=-Double.MAX_VALUE;
		  minimums[i]=Double.MAX_VALUE;
		  for (j=0;j<m;j++)
			{
			  minimums[i] = Math.min(minimums[i],data[j][i]);
			  maximums[i] = Math.max(maximums[i],data[j][i]);
			}
		}
	}

	/**
	   Rescale the data array data[0...npoints-1][0....m-1] into [0,1].
	*/
	public void rescale()
	{
	  int i,j;
	  double denom;

	  // scale each dimesnion in turn into the interval [0,1]
	  for (i=0;i<n;i++)
		{
		  denom = (maximums[i]==minimums[i])? 1.0 : maximums[i]-minimums[i];
		  for (j=0;j<m;j++)
			{
			  _data[j][i] = (_data[j][i] - minimums[i])/denom;
			}
		}
	}
	/**
	   Restore the data array data[0...npoints-1][0....m-1] back to its original values
	*/
	public void restoreData()
	{
	  int i,j;
	  double denom;

	  for (i=0;i<n;i++)
		{
		  denom = maximums[i]-minimums[i];
		  for (j=0;j<m;j++)
			{
			  _data[j][i] = _data[j][i]*denom + minimums[i];
			}
		}
	}
  }

} // End class.
