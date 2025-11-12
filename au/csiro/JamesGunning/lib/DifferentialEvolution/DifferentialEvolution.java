package au.csiro.JamesGunning.lib.DifferentialEvolution;

/**
<p> 
   DifferentialEvolution global genetic optimizer ala Storn and Price
   This class is a java implementation of the DifferentialEvolution
   genetic optimiser found at
   http://www.icsi.berkeley.edu/~storn/code.html I think the original ref
   is Storn, Rainer and Price, Kenneth (1995). Differential Evolution - a
   Simple and Efficient Adaptive Scheme for Global Optimization over
   Continuous Spaces. Technical Report TR-95-012, ICSI, March
   1995. PostScript-file downloadable from
   http://www.icsi.berkeley.edu/techreports/1995.abstracts/tr-95-012.html.
   <\p>
   <p>
   There are lots of other papers. Default parameters are chosen
   from a montage assembled from
   (1) http://www.icos.ethz.ch/research/wseas02.pdf,
   (2) "A Parameter Study for Differential Evolution" Gamperle et al (they recommend DE/best/2/bin),
   (3) http://www.icsi.berkeley.edu/~storn/code.html#prac
   and other papers.
   <\p>

   @author James.Gunning@csiro.au
   @version 1.0, 27/12/2002 
*/

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.jet.random.*;
import cern.jet.random.engine.MersenneTwister;
import au.csiro.JamesGunning.lib.ErrorHandler;


public class  DifferentialEvolution
  /**
	 The minimize() methods minimizes a function of n variables
	 A method that computes the function value at any point
	 must be supplied.  (See DE_methods.java)
	 The minimize() method provides the simplest user access to the DE
	 routines.@param  n         The number of arguments of the function to minimize
	 @param  x         The initial estimate of the minimum point
	 @param minclass   A class that implements the DE_methods
	 interface (see the definition in
	 DE_methods.java). The class must define:
	 1.) a method, f_to_minimize, to minimize.
	 f_to_minimize must have the form<br>
	 
	 public static double f_to_minimize(DoubleMatrix1D x)<br>


	 where x is the vector of arguments to the function
	 and the return value is the value of the function
	 evaluated at x.

	 2.) a method to generate an initial generation of
	 candidate points in the space of interest. It must have form<br>

	 create_initial_generation(int NP, DoubleMatrix2D c)<br>
	 where c is a cern style DoubleMatrix2D whose rows run over the
	 initial population, and whose columns follow the parameter index of the
	 x used in f_to_minimize(x). The user must provide NP independent
	 points from which subsequent generations are bred. c is already allocated.

	 3.) A bunch of other methods that can generally be transparent (do-nothing),
	 but which facilitate monitoring of the convergence, and emergency bailout
	 type procedures.
  */
{
  private ErrorHandler EH;
  private DifferentialEvolution DE;
  private MersenneTwister ranEngine;
  private AbstractDistribution Ugen,NPgen;
  private int nfeval=0;          /* number of function evaluations     */
  private DenseDoubleMatrix2D c, d;
  private DenseDoubleMatrix2D pOld, pNew, pSwap;

  // Major control knobs.
  /**
	 Maximum size of population each generation
  */
  public static final int MAXPOP = 500;
  /**
	 maximum number of generations (termination usually earlier than this)
  */
  public static final  int MAXGENERATIONS = 1000;
  /**
	 Weighting factor
  */
  double F;
  /**
	 Crossover constant
  */
  double CR;
  /**
	 Strategy
  */
  int strategy;

  /**
	 Default termination epsilon for variance of population
  */
  public static final double DEFAULT_TERMINATION_EPSILON=1.0e-4;

  /* number of population members       */
  private int   NP=10;
  /* Dimension of parameter vector      */
  private int   D=1;
  /* Seed */ 
  private int seed;

  /* computes the cost variance         */
  private double cvar;

  /* termination tolerance on cost variance */
  private double epsilon;

  /**
	 Make a new DifferentialEvolution object for solving a D dimensional
	 optimization problem
  */
  public DifferentialEvolution(int D)
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();
	if(D<1)
	  {
		throw new RuntimeException("Dimensionality of problem must be a +ve integer");
	  }
	this.D=D;
  }
  /**
	 Make a new DifferentialEvolution object for solving a D dimensional
	 optimization problem, using the error handler supplied.
  */	
  public DifferentialEvolution(int D,int seed,ErrorHandler EH)
  {
	// handle errors: can change how these are handled if desired
	this.EH=EH;
 	if(D<1)
	  {
		throw new RuntimeException("Dimensionality of problem must be a +ve integer");
	  }
	this.D=D;
	this.seed = seed;
  }

  /** 
	  Example test harness
  */
  public static void main(String argv[])
  {
	ExampleDE_methods exampleMethods = new ExampleDE_methods();

	DifferentialEvolution DE = new DifferentialEvolution(10);

	// Control knobs...FIXME check epsilon!!!
	DE.setControlParameters(-1.0,-1.0,-1,-1,-1);

	// the actual minimiser
	DE.minimize(exampleMethods);
  }


  /**
	 This must be called before a minimisation loop to set the major control parameters. See the spiel
	 at http://www.icsi.berkeley.edu/~storn/code.html#basi and read lots of papers if you want to 
	 fool around with these. All parameters ought to be positive, and 0<F<1, 0<CR<1.
	 If you don't know, pass a -ve number for each parameter  and sane defaults will be selected
	 Defaults for the population size are mainly reasonable, but you can enlarge it if desired.

			There are some simple rules which are worth following:
			1)  F is usually between 0.5 and 1 (in rare cases > 1)
			2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first
			3)  To start off NP = (5-10)*D is a reasonable choice. Increase NP if misconvergence
			           happens.                                                                     
			4)  If you increase NP, F usually has to be decreased
			5)  When the DE/best... schemes fail DE/rand... usually works and vice versa

	 Strategy codes (read the papers for what these mean)<br>
	 DE/best/1/exp (strategy == 1) <br>
	 DE/rand/1/exp (strategy == 2) /* strategy DE1 in the techreport <br>
	 DE/rand-to-best/1/exp (strategy == 3) /* similiar to DE2 but generally better <br>
	 DE/best/2/exp (strategy == 4)<br>
	 DE/rand/2/exp (strategy == 5)<br>
	 DE/best/1/bin (strategy == 6)<br> 
	 DE/rand/1/bin (strategy == 7) <br>
	 DE/rand-to-best/1/bin (strategy == 8) <br>
	 DE/best/2/bin (strategy == 9)<br>
	 DE/rand/2/bin default (strategy>=10 || strategy<1)<br>
			
  */
  public void setControlParameters(double F, double CR, int strategy, double epsilon, int nPopulation)
  {
	// Defaults

	// strategy=2;      // DE/rand/1/exp
	this.strategy=9;    // DE/best/2/bin
	this.F=0.5;         // F values near 0.8 seem to cause misconvergence in some seismic problems
	this.CR=0.5;
	this.epsilon = DEFAULT_TERMINATION_EPSILON;


	// User overrides

	if (F>=0 && F<=1) this.F = F;
	if (CR>=0 && CR<=1) this.CR = CR;
	if (strategy>=1) this.strategy = strategy;
	if (epsilon>0) this.epsilon = epsilon;
	if (nPopulation>0) this.NP = nPopulation;
 }

  /**
	Generate minimum of solution using np population members
	np = suggested number of members in each generation.
	methods = class implementing DE_methods.
  */
  public void minimize(DE_methods methods)
  {
	int   i, j, L, n;      /* counting variables                 */
	int   r1, r2, r3, r4;  /* placeholders for random indexes    */
	int   r5;              /* placeholders for random indexes    */

	int   imin;            /* index to member with lowest energy */
	int   refresh=10;      /* refresh rate of screen output      */
	int   gen;   
	int   NP_min;          /* default suggestions */

	double trial_cost;      /* buffer variable                    */
	double inibound_h=1;    /* upper parameter bound              */
	double inibound_l=-1;   /* lower parameter bound              */
	DoubleMatrix1D tmp, best, bestit; /* members  */
	DoubleMatrix1D cost;    /* obj. funct. values                 */
	double cmean;           /* mean cost                          */
	double cmin;            /* help variables                     */
	double tmpVal;

	ranEngine = new MersenneTwister(seed);
	Ugen = new Uniform(0,1,ranEngine);


	// Allocs based on dimension
	cost = new DenseDoubleMatrix1D(MAXPOP);
	tmp = new DenseDoubleMatrix1D(D);
	best = new DenseDoubleMatrix1D(D);
	bestit = new DenseDoubleMatrix1D(D);
	c = new DenseDoubleMatrix2D(MAXPOP,D);
	d = new DenseDoubleMatrix2D(MAXPOP,D);


	/*------Initialization------------------------------------------------*/

	// Size of population
	NP_min=Math.max(6,10*D);  // p 356 S & P's paper. Mustn't have fewer than 6 members
	
	NP=Math.max(NP,NP_min);   // user may have specified a greater population

	if (NP>MAXPOP) 
	  {
		throw new RuntimeException("Dimensionality of problem too big for DifferentialEvolution"
								   +" defaults. Your problem is probably too hard");
	  }

	// Make initial generation
	methods.create_initial_generation(NP,c);

	for (i=0; i<NP; i++)  // loop over initial population to get costs
	  {
		cost.setQuick(i , methods.f_to_minimize(c.viewRow(i))); /* obj. funct. value */
	  }
	cmin = cost.getQuick(0);
	imin = 0;
	for (i=1; i<NP; i++)
	  {
		if (cost.getQuick(i)<cmin)
		  {
			cmin = cost.getQuick(i);
			imin = i;
		  }
	  }

	assignd(best,c.viewRow(imin));            /* save best member ever          */
	assignd(bestit,c.viewRow(imin));          /* save best member of generation */

	pOld = c; /* old population (generation G)   */
	pNew = d; /* new population (generation G+1) */

	/*=======================================================================*/
	/*=========Iteration loop================================================*/
	/*=======================================================================*/

	NPgen = new Uniform(0,NP,ranEngine);

	gen = 0;                          /* generation counter reset */
	while ((gen < MAXGENERATIONS))
	  {
		gen++;
		imin = 0;

		for (i=0; i<NP; i++)         /* Start of loop through ensemble  */
		  {
			do                        /* Pick a random population member */
			  {                         /* Endless loop for NP < 2 !!!     */
				r1 = (int)(NPgen.nextDouble());
			  } while(r1==i);            

			do                        /* Pick a random population member */
			  {                         /* Endless loop for NP < 3 !!!     */
				r2 = (int)(NPgen.nextDouble());
			  } while((r2==i) || (r2==r1));

			do                        /* Pick a random population member */
			  {                         /* Endless loop for NP < 4 !!!     */
				r3 = (int)(NPgen.nextDouble());
			  } while((r3==i) || (r3==r1) || (r3==r2));

			do                        /* Pick a random population member */
			  {                         /* Endless loop for NP < 5 !!!     */
				r4 = (int)(NPgen.nextDouble());
			  } while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do                        /* Pick a random population member */
			  {                         /* Endless loop for NP < 6 !!!     */
				r5 = (int)(NPgen.nextDouble());
			  } while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


			/*=======Choice of strategy===============================================================*/
			/*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
			/*=======DE :  stands for Differential Evolution==========================================*/
			/*=======x  :  a string which denotes the vector to be perturbed==========================*/
			/*=======y  :  number of difference vectors taken for perturbation of x===================*/
			/*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
			/*                                                                                        */
			/*=======There are some simple rules which are worth following:===========================*/
			/*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
			/*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
			/*=======3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=*/
			/*           happens.                                                                     */
			/*=======4)  If you increase NP, F usually has to be decreased============================*/
			/*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/


			/*=======EXPONENTIAL CROSSOVER============================================================*/

			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (strategy == 1) /* strategy DE0 (not in our paper) */
			  {
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D);
				L = 0;
				do
				  {                       
					tmp.setQuick(n , bestit.getQuick(n) + F*(pOld.getQuick(r2,n)-pOld.getQuick(r3,n)));
					n = (n+1)%D;
					L++;
				  } while((Ugen.nextDouble() < CR) && (L < D));
			  }
			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (strategy == 2) /* strategy DE1 in the techreport */
			  {
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D);
				L = 0;
				do
				  {                       
					tmp.setQuick(n, pOld.getQuick(r1,n) + F*(pOld.getQuick(r2,n)-pOld.getQuick(r3,n)));
					n = (n+1)%D;
					L++;
				  } while((Ugen.nextDouble() < CR) && (L < D));
			  }
			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
			/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (strategy == 3) /* similiar to DE2 but generally better */
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				L = 0;
				do
				  {                       
					tmp.setQuick(n, 
								 tmp.getQuick(n) + F*(bestit.getQuick(n) - tmp.getQuick(n)) 
								 + F*(pOld.getQuick(r1,n)-pOld.getQuick(r2,n)));
					n = (n+1)%D;
					L++;
				  } while((Ugen.nextDouble() < CR) && (L < D));
			  }
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (strategy == 4)
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				L = 0;
				do
				  {                           
					tmp.setQuick(n , bestit.getQuick(n) + 
								 (pOld.getQuick(r1,n)+pOld.getQuick(r2,n)-pOld.getQuick(r3,n)-pOld.getQuick(r4,n))*F);
					n = (n+1)%D;
					L++;
				  } while((Ugen.nextDouble() < CR) && (L < D));
			  }
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (strategy == 5)
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				L = 0;
				do
				  {                           
					tmp.setQuick(n , pOld.getQuick(r5,n) + 
								 (pOld.getQuick(r1,n)+pOld.getQuick(r2,n)-pOld.getQuick(r3,n)-pOld.getQuick(r4,n))*F);
					n = (n+1)%D;
					L++;
				  } while((Ugen.nextDouble() < CR) && (L < D));
			  }

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (strategy == 6) 
			  {
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				for (L=0; L<D; L++) /* perform D binomial trials */
				  {
					if ((Ugen.nextDouble() < CR) || L == (D-1)) /* change at least one parameter */
					  {                       
						tmp.setQuick(n, bestit.getQuick(n) + F*(pOld.getQuick(r2,n)-pOld.getQuick(r3,n)));
					  }
					n = (n+1)%D;
				  }
			  }
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (strategy == 7) 
			  {
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				for (L=0; L<D; L++) /* perform D binomial trials */
				  {
					if ((Ugen.nextDouble() < CR) || L == (D-1)) /* change at least one parameter */
					  {                       
						tmp.setQuick(n, pOld.getQuick(r1,n) + F*(pOld.getQuick(r2,n)-pOld.getQuick(r3,n)));
					  }
					n = (n+1)%D;
				  }
			  }
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (strategy == 8) 
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				for (L=0; L<D; L++) /* perform D binomial trials */
				  {
					if ((Ugen.nextDouble() < CR) || L == (D-1)) /* change at least one parameter */
					  {                       
						tmp.setQuick(n, 
									 tmp.getQuick(n) + F*(bestit.getQuick(n) - tmp.getQuick(n)) 
									 + F*(pOld.getQuick(r1,n)-pOld.getQuick(r2,n)));
					  }
					n = (n+1)%D;
				  }
			  }
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (strategy == 9)
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				for (L=0; L<D; L++) /* perform D binomial trials */
				  {
					if ((Ugen.nextDouble() < CR) || L == (D-1)) /* change at least one parameter */
					  {                       
						tmp.setQuick(n, bestit.getQuick(n) + 
									 (pOld.getQuick(r1,n)+pOld.getQuick(r2,n)-pOld.getQuick(r3,n)-pOld.getQuick(r4,n))*F);
					  }
					n = (n+1)%D;
				  }
			  }
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else
			  { 
				assignd(tmp,pOld.viewRow(i));
				n = (int)(Ugen.nextDouble()*D); 
				for (L=0; L<D; L++) /* perform D binomial trials */
				  {
					if ((Ugen.nextDouble() < CR) || L == (D-1)) /* change at least one parameter */
					  {                       
						tmp.setQuick(n ,pOld.getQuick(r5,n) + 
									 (pOld.getQuick(r1,n)+pOld.getQuick(r2,n)-pOld.getQuick(r3,n)-pOld.getQuick(r4,n))*F);
					  }
					n = (n+1)%D;
				  }
			  }


			/*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/

			trial_cost = methods.f_to_minimize(tmp);  /* Evaluate new vector in tmp[] */
			nfeval++;
			if (trial_cost <= cost.getQuick(i))   /* improved objective function value ? */
			  {                                  
				cost.setQuick(i, trial_cost); 
				assignd(pNew.viewRow(i),tmp);
				if (trial_cost<cmin)          /* Was this a new minimum? */
				  {                               /* if so...*/
					cmin=trial_cost;           /* reset cmin to new low...*/
					imin=i;
					assignd(best,tmp);           
				  }                           
			  }                            
			else
			  {
				assignd(pNew.viewRow(i),pOld.viewRow(i)); /* replace target with old value */
			  }
		  }   /* End mutation loop through pop. */
					   
		assignd(bestit,best);  /* Save best population member of current iteration */

		/* swap population arrays. New generation becomes old one */

		pSwap = pOld;
		pOld  = pNew;
		pNew  = pSwap;

		/*----Compute the energy variance, for monitoring and termination criterion -----------*/

		cmean = 0.;          /* compute the mean value first */
		for (j=0; j<NP; j++) cmean += cost.getQuick(j);
		cmean/=(double)NP;

		cvar = 0.;           /* now the variance of the cost  */
		for (j=0; j<NP; j++)
		  {
			cvar += ((tmpVal=cost.getQuick(j) - cmean))*tmpVal;
		  }
		cvar /=(NP-1.0);


		/*----Output part----------------------------------------------------------*/
		/*
		  if (gen%refresh==1)   // display after every refresh generations
		  {
			EH.printf("%n%n Generation=%d  NFEs=%d ",gen,add(nfeval));
			EH.printf("%n NP=%d    F=%-4.2g    CR=%-4.2g   cost-variance= %-10.5g%n",
					 NP,F,CR,cvar);
		  }
					 
		  EH.printf("%d   %-15.10g%n",nfeval,cmin);
		*/

		// Feedback to controlling class
		methods.updateNewPoint(cmin,best,gen,nfeval);

		// JG's extra termination criteria... 
		// For well-scaled cost functions (e.g. parabolic minima f=x^2), this is ok.
		// Usually choose syomething like epsilon = precision^2.
		if (cvar < epsilon) break;

		// Early termination escape if something goes wrong
		if (methods.wasAborted()) break;
	  }

  }

  /**
	 After solving, the population cost variance is accessible here. Its usually interesting
	 as evidence of convergence
  */
  public double getPopulationCostVariance()
  {
	return cvar;
  }
  /**
	 For use after solving, the number of function evaluations.
  */
  public int getNumberOfFunctionEvaluations()
  {
	return nfeval;
  }

  private void  assignd(DoubleMatrix1D a,  DoubleMatrix1D b)
  // vector b to vector a. 
  {
	a.assign(b);
  }
}
