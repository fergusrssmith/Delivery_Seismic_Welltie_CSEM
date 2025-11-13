/**********************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
package au.csiro.JamesGunning.CSEM.OneDinversion;

import java.io.*;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.util.*;

import static au.csiro.JamesGunning.lib.util.MathUtil.*;

/** 
	Illustration of path integral formulation for finding lowest-lying valleys
	between modes of fixed surface. Sandpit code only
*/
public class ModePairConnector
{
  ErrorHandler EH = new ErrorHandler();

  public static void main(String argv[]) 
  {
	ModePairConnector MPC = new ModePairConnector();
	MPC.findModePairConnection();
  }

  void findModePairConnection()
  {
	int i,j,k,l,iM,d=2,n,N;
	N=30;
	double s,frac,L,Lo,rho;

	double[][] m;
	double[] f,grad,dm2,dm2t,gradL;
	double[][] gradf;

	// endpoints
	double[] mA,mB;	
	mA = new double[2];
	mB = new double[2];
	mA[0]=0; mA[1]=1;
	mB[0]=1; mB[1]=0;

	n = d*N;
	m = new double[N+2][d];
	dm2 = new double[N+1];
	dm2t = new double[N+1];
	gradf = new double[N+2][d];
	f = new double[N+2];
	gradL = new double[n];

	// linearly interpolate starting state
	for (i=0;i<=(N+1);i++)
	  {
		frac = i/(double)(N+1);
		for(j=0;j<d;j++)
		  {
			m[i][j] = frac*mB[j] + (1.0-frac)*mA[j];
		  }		
	  }

	// SD stuff
	double[][] mt = new double[N+2][d];
	double EPS=1e-5;
	System.arraycopy(mA,0,mt[0],0,d);
	System.arraycopy(mB,0,mt[N+1],0,d);

	int IMAX=2000;
	for(iM=0; iM<=IMAX; iM++)
	  {
		for (i=0;i<=N;i++)
		  {
			for (s=0,j=0;j<d;j++) s += MathUtil.sqr(m[i+1][j]-m[i][j]);
			dm2[i] = Math.sqrt(s);
		  }
		// get f's and grad components
		for (i=0; i<=N+1; i++) f[i] = testf(m[i],gradf[i]);
		for (i=0,L=0; i<=N; i++) L += (f[i]+f[i+1])*dm2[i];
		//EH.println("f="+L);
		// Get gradL
		for (k=0,i=1; i<=N; i++)
		  {
			for(j=0;j<d;j++)
			  {
				gradL[k] = gradf[i][j] * (dm2[i-1]+dm2[i]);
				gradL[k] += (f[i]+f[i-1]) * (m[i][j]-m[i-1][j])/dm2[i-1];
				gradL[k] -= (f[i]+f[i+1]) * (m[i+1][j]-m[i][j])/dm2[i];
				k++;
			  }
		  }
		// Approximate Hessian ??

		// Optimise
		// SD method
		Lo=L;

		for (rho=1; rho>=EPS; rho/=2.0)
		  {
			for (k=0,i=1;i<=N;i++)
			  {
				for(j=0;j<d;j++)
				  {
					mt[i][j] = m[i][j] -rho*gradL[k];
					k++;
				  }		
			  }
			for (i=0;i<=N;i++)
			  {
				for (s=0,j=0;j<d;j++) s += MathUtil.sqr(mt[i+1][j]-mt[i][j]);
				dm2t[i] = Math.sqrt(s);
			  }
			// get f's and grad components
			for (i=0; i<=N+1; i++) f[i] = testf(mt[i],gradf[i]);
			for (i=0,L=0; i<=N; i++) L += (f[i]+f[i+1])*dm2t[i];

			if (L<Lo) 
			  {
				// copy solution
				for (i=1;i<=N;i++)
				  {
					for(j=0;j<d;j++) m[i][j] = mt[i][j];
				  }
				System.err.println("L= "+L);
				break;
			  }
		  }
		if (rho<=EPS) break;		
	  }
	if (true || iM<=IMAX)
	  {
		//EH.println("Solution:");
		System.err.println("Nits="+iM);
		for (i=1;i<=N;i++)
		  {
			for(j=0;j<d;j++)
			  {
				EH.printf("%g ",m[i][j]);
			  }
			EH.println("");
		  }
		EH.println("");
		for (i=0;i<=N;i++)
		  {
			frac=i/(double)N;
			EH.printf("%g %g%n",frac,Math.log10(1.0+10.0-Math.pow(10,frac)));
		  }
	  }

  }
  // test debug function for findModePairConnection()
  private static double testf(double[] m,double[] grad)
  {
	int j;
	double s=0;
	double A,xb,yb;

	if (!true)
	  {
		xb=0;
		yb=1;
		A=1.0;
		s = sqr(m[0]-xb)+A*sqr(m[1]-yb);
		grad[0] = 2*(m[0]-xb);
		grad[1] = 2*A*(m[1]-yb);
	  }
	else
	  {
		double x,y;
		double R = 1.0+10.0;
		double lt = Math.log(10.0);

		A=50;
		x=Math.pow(10,m[0]);
		y=Math.pow(10,m[1]);

		s = A*sqr(x+y-R)+x*x+y*y;
		grad[0] = A*(x+y-R) * lt*x + 2*x;
		grad[1] = A*(x+y-R) * lt*y + 2*y;
	  }
	return s;
  }
}