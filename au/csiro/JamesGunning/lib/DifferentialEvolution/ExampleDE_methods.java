/**
   Interface for function to be optimised (minimised using
   DifferentialEvolution. No Hessians or Gradients required.
   The new methods abortNow(), wasAborted(), updateNewPoint(), and
   setFlags() can all be fairly trivial, and do virtually nothing.
*/

package au.csiro.JamesGunning.lib.DifferentialEvolution;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import au.csiro.JamesGunning.lib.ErrorHandler;


public class ExampleDE_methods implements DE_methods 
{
  private ErrorHandler EH;

  private double f_current;
  private DoubleMatrix1D x_current;

  public boolean abortNow(){return false;}
  public boolean wasAborted(){return false;}

  public ExampleDE_methods()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();
  }

  /**
	 example of how the DE_methods class gets feedback from the solver...
  */
  public void updateNewPoint(double f, DoubleMatrix1D x, int itCount, int nfeval)
  {
	f_current = f;
	x_current = x;
	// may want to print these etc....

	/*----Output part----------------------------------------------------------*/
	if (itCount%10==1)   /* display after every refresh=10 generations */
	  {
		EH.println("\n\n\n Best-so-far cost funct. value= %-15.10g, num. func. evaluations=%d\n",
				   f, nfeval);			
		for (int j=0;j<x.size();j++)
		  {
			EH.printf("best[%d]=%-15.10g\n",j,x.getQuick(j));
		  }
	  }
  }

  public void setFlags(boolean setOrReset){}

  /**
	 Simple sum of squares and related beasts. Function we wish to optimise
	 is a function of x[0], x[1]...x[x.size()-1]
  */
  public double f_to_minimize(DoubleMatrix1D x)
  {
	double f,flast=0,chisq=0;

	for (int i=0;i<x.size();i++)
	  {
		f = x.getQuick(i)+i;
		// humps 
		//chisq+=f*f;
		// set of increasingly narrow valleys...
		//chisq+= f*f*i*i;
		// correlated
		//flast+=f;
		//chisq+=f*flast;
		// multimodal
		chisq+=f*f;
	  }
	// chisq = 2+Math.cos(5.0*chisq)+chisq/2.0;
	// return chisq*chisq;
	return chisq;
  }
  
  /**
	 Make initial population. Must have c.rows()>=NP
  */
  public void create_initial_generation(int NP, DoubleMatrix2D c)
  {
	int i,j;

	// toy problem. Initial generation uniformly distributed in a hypercube
	double inibound_l = -10;
	double inibound_h = 10;

	for (i=0; i<NP; i++)  // loop over initial population
	  {
		for (j=0; j<c.columns(); j++) /* spread initial population members */
		  {
			c.setQuick(i,j, inibound_l + Math.random()*(inibound_h - inibound_l));
		  }
	  }
  }
}
