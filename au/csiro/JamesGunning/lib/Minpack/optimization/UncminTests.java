package au.csiro.JamesGunning.lib.Minpack.optimization;

import au.csiro.JamesGunning.lib.Minpack.optimization.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*;

/**
*
* This class tests the Uncmin_f77 class for 
* scaling issues and piecewise cts functions
* Hack to your needs.
*
*@author James Gunning)
*@version 1 --- Jan 3, 2004
*
*/


public class UncminTests extends Object implements Uncmin_methods 
{
  int n;
  int it=0;
  boolean  showX = false;
  boolean  scaleProperly = false;
  double[] distort;                // to alter scales on half the variables

  // Control variables for solver.
  // see dfault_f77 for example of how this is done in Uncmin....

  double[] f_best = new double[2];             // f[1] has the function minimum at exit
  int[] info = new int[2];                     // status information about the code
  double[] fscale = new double[2];
  int[] method = new int[2];
  int[] is_exp = new int[2];
  int[] msg = new int[2];
  int[] n_good_digits = new int[2];
  int[] max_its = new int[2];
  int[] is_gradient_supplied = new int[2];
  int[] is_hessian_supplied = new int[2];
  double[] trust_radius = new double[2];
  double[] gradient_tol = new double[2];
  double[] max_step_size = new double[2];
  double[] step_tol = new double[2];

  // workspace variables
  private double[][] Hessian_wk;
  private double[] Hessian_diag_wk;
  private double[] fullParameters_best;
  private double[] grad_sub_best;
  private double[] typical_values;
  private double[] fullParameters;



  public static void main(String argv[]) 
  {
	UncminTests test = new UncminTests();
	ErrorHandler EH = new ErrorHandler();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("UncminTests",EH);

	options.addMandatoryArg(1,"dimension","dimension of problems to do");
	options.addOption("-x",null,0,null,null,"Show X section");
	options.addOption("-s",null,0,null,null,"Scale properly!");
	// parse options
	options.parseOptions(argv,true);

	// get dimensionality
	try
	  {
		test.n = Integer.parseInt(options.getMandatoryArg(1));
	  }
	catch (NumberFormatException nfe)
	  {
		options.Usage();
		EH.die_err("Incorrect arguments " + argv[0]);
	  }
	if (options.isOptionSet("-x"))  test.showX = true;
	if (options.isOptionSet("-s"))  test.scaleProperly = true;
	// now optimise
	test.optimise();
  }

  private void optimise()
  {
	int i;

	// optimiser setup
	// Tolerances...
	double epsilon  = 1.12e-16; // default 1.12e-16;
	// FIXME can the gradient tolerance can be weakened??????
	gradient_tol[1] = Math.pow(epsilon,1.0/3.0);
	max_step_size[1] = 0;     // 0.0 default
	step_tol[1]=Math.sqrt(epsilon);
	trust_radius[1]=-1;   // -1 is default

	// Flags for methods/outputs...
	method[1]=1;          // 1 = line-search
	// NB. UNCMIN can be a bit wierd under BFGS if scaling is dodgy. I believe scaling issues are
	// fixed so BFGS (Secant) is safe
	is_exp[1]=1;
	msg[1]=0;             // msg=8 seems to suppress most of the noise... :-). msg=0,1 is very noisy
	n_good_digits[1]=-1;  // -1 gives defaults
	max_its[1]=1000;
	is_gradient_supplied[1]=0;
	is_hessian_supplied[1]=0;

	// workspace
	fullParameters = new double[n+1];
	fullParameters_best = new double[n+1];
	grad_sub_best = new double[n+1];
	typical_values = new double[n+1];
	Hessian_wk = new double[n+1][n+1];
	Hessian_diag_wk = new double[n+1];

	for (i=1;i<=n;i++) typical_values[i]=1;

	// adhoc scale distortion
	distort = new double[n+1];
	for (i=1;i<=n;i++) distort[i] = (i<=n/2) ? 1 : 0.01;
	if (scaleProperly) 
	  {
		System.out.println("Scaling....");
		for (i=1;i<=n;i++) typical_values[i]=1/distort[i];
	  }

	fscale[1]=1;
	it=0;

	Uncmin_f77.optif9_f77(n, fullParameters, this,
						  typical_values,fscale,method,is_exp,msg,n_good_digits,
						  max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
						  gradient_tol,max_step_size,step_tol,
						  fullParameters_best,f_best,grad_sub_best,info,Hessian_wk,Hessian_diag_wk);
	// show X section through worse parameter
	if (showX)
	  {
		int NSTEPS=200;
		int		j=1;
		double diffMax =0,xWorst=0;
		for (i=1;i<=n;i++)
		  {
			if (Math.abs(fullParameters_best[i] - 1) > diffMax)
			  {
				j=i;
				diffMax = Math.abs(fullParameters_best[i] - 1);
				xWorst = fullParameters_best[i];
			  }
		  }
		for (i=-NSTEPS;i<=NSTEPS;i++)
		  {
			fullParameters_best[j]=xWorst+(double)i/NSTEPS;
			System.out.println(fullParameters_best[j]+" "+f_to_minimize(fullParameters_best));
		  }
	  }
  }

  public double f_to_minimize(double[] x)
  {
	int i;
	double xn,xm;
	double f=1,a=1.0,eps=0.01;
	for (i=1;i<=n;i++)
	  {
		xn =  (x[i]-a)*distort[i];
		f += xn*xn + eps*Math.abs(xn);
	  }
	// rotate things a bit....
	for (i=1;i<n;i++)
	  {
		xn =  (x[i]-a)*distort[i];
		xm =  (x[i+1]-a)*distort[i+1];
		f += 0.4*(xn*xm +  eps*Math.abs(xn)*Math.abs(xm));
	  }
	it++;

	return f;
  }
  
  /**
   The solver throws execution to here everytime it gets a major new
   point. This is where we can monitor  things.
  */
  public void updateNewPoint(double fnew, double[] xnew)
  {
	System.out.println("it = "+it +" f = "+fnew);
  }

  public boolean abortNow() {return false;};

  public boolean wasAborted() {return false;};

  public void setFlags(boolean setOrReset) {};

  // null methods required since we can't figure these out analytically
  public void gradient(double times[], double g[])
  {};
  public void hessian(double times[], double h[][]) 
  {};
}
