/**
   The original Uncmin from Steve Verrill has been augmented with extra
   calls to enable graceful behaviour with piecewise smooth functions.
   Any implementation of Uncmin_methods must contain these methods
   in addition to those described in Steve's documentation<br>
   The new methods abortNow(), wasAborted(), updateNewPoint(), and
   setFlags() can all be fairly trivial, and do virtually nothing.
*/

package au.csiro.JamesGunning.lib.Minpack.optimization;

// JG HACK: made this public
public interface Uncmin_methods {

  // JG hack for rapid exit
  /**
	 Occasionally we can detect if the minimiser is getting lost, or is
	 on an unrecoverable spiral. Provide this method to escape gracefully
	 if so. This can be a simple do-nothing method.
  */
  public boolean abortNow();

  // JG hack. test for rapid exit
  /**
	 If we forced a fall out of the solver by calling abortNow(),
	 typically we set a boolean flag in abortNow() to indicate so. The
	 calling environment may want to interrogate this flag, hence the
	 need for this method. This can be a simple do-nothing method.
  */
  public boolean wasAborted();

  // JG hack. Snapshot of current point for tracking (can contain nothing)
  /**
	 We often want to monitor the progress of the solver. Here we have access to the
	 current function value and variables, for printing etc. Do not alter any xnew[] in this
	 method. This can be a simple do-nothing method.
  */
  public void updateNewPoint(double fnew, double[] xnew);

  // JG hack for altering function
  /** 
	  The function we want to minimise may have certain flags which controls its behaviour. Here
	  we can set these. This can be a simple do-nothing method.
  */
  public void setFlags(boolean setOrReset);
  
  public double f_to_minimize(double x[]);
  
  public void gradient(double x[], double g[]);

  public void hessian(double x[], double h[][]);

}
