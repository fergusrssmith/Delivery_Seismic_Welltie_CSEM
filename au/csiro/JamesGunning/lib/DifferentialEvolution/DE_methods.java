/**
   Interface for function to be optimised (minimised using
   DifferentialEvolution. No Hessians or Gradients required.
   The new methods abortNow(), wasAborted(), updateNewPoint(), and
   setFlags() can all be fairly trivial, and do virtually nothing.
*/
package au.csiro.JamesGunning.lib.DifferentialEvolution;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;


public interface DE_methods {

  /**
	 Occasionally we can detect if the minimiser is getting lost, or is
	 on an unrecoverable spiral. Provide this method to escape gracefully
	 if so. This can be a simple do-nothing method.
  */
  public boolean abortNow();

  /**
	 If we forced a fall out of the solver by calling abortNow(),
	 typically we set a boolean flag in abortNow() to indicate so. The
	 calling environment may want to interrogate this flag, hence the
	 need for this method. This can be a simple do-nothing method.
  */
  public boolean wasAborted();

  /**
	 We often want to monitor the progress of the solver. Here we have access to the
	 current function value and variables, for printing etc. Do not alter any xnew in this
	 method. This can be a simple do-nothing method.
  */
  public void updateNewPoint(double fnew, DoubleMatrix1D xnew,int itCount,int nfeval);

  /** 
	  The function we want to minimise may have certain flags which controls its behaviour. Here
	  we can set these. This can be a simple do-nothing method.
  */
  public void setFlags(boolean setOrReset);
  
  /**
	 The main method generating the function of d variables x[0]....x[d-1]
	 which we want to minimise. x.size() must equal d.
  */
  public double f_to_minimize(DoubleMatrix1D x);

  /**
	 Mandatory method to generate initial population (of size NP<=C.rows()) into matrix C:
	 rows of C = member index (0...NP-1), columns of C = parameters index.
	 Obviously C.columns() = x.size() in the f_to_minimize() method.
  */
  public void create_initial_generation(int NP,DoubleMatrix2D C);
}
