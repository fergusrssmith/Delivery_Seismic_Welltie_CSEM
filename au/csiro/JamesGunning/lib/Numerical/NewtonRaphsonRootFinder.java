/**********************************************************************
Copyright (C) 2010, CSIRO,  all rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of CSIRO nor the names of its
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

package au.csiro.JamesGunning.lib.Numerical;

import static java.lang.Math.*;

/**
   @author James.Gunning@csiro.au
   @version 1.0, 30/4/2010
*/

public class NewtonRaphsonRootFinder
{
  FunctionAndDerivative F;
  int success_code=0;
  
  public NewtonRaphsonRootFinder(FunctionAndDerivative F)
  {
	this.F = F;
  }

  public int success()
  {
	return success_code;
  }

  /**
	 Find root of function bracketed by x1 x2, within tolerance eps
	 If success, success()==0, else
	 success()=1   => not bracketed at call
	 success()=2   => too many iterations
  */
  public double NewtonRaphsonFind(double x1, double x2, double xacc)
  {
	int j,MAXIT=100;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	success_code = 0;
	fl = F.f(x1);
	fh = F.f(x2);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	  {
		success_code = 1;
		return -Double.MAX_VALUE;
	  }
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
	  xl=x1;
	  xh=x2;
	} else {
	  xh=x1;
	  xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold = abs(x2-x1);
	dx = dxold;

	f = F.f_and_dfdx(rts); 
	df = F.get_dfdx();
	for (j=1;j<=MAXIT;j++) 
	  {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0) || (abs(2.0*f) > abs(dxold*df))) 
		  {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		  } 
		else 
		  {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		  }
		if (abs(dx) < xacc) return rts;
		f = F.f_and_dfdx(rts); 
		df = F.get_dfdx();
		if (f < 0.0)
		  xl=rts;
		else
		  xh=rts;
	  }
	success_code = 2;
	return -Double.MAX_VALUE;
  }
}
