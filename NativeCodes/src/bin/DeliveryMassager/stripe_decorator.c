/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
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
#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "commandline/commandline.h" 
#include "NR/nrutil.h"
#include "util.h"
#include "modified_left_Cholesky.h"
#include "clapack_cblas_wrappers.h"


#define ARGS_USAGE ""

int real_no=0;
int seed=-1;

/*
  Test program for some downscaling ideas.....
  Run, then: Use plotutils: graph -TX  -q 1 -B realisation.%d.txt

*/


CLOption options[] = 
  {
	{'s',NULL,1,  &seed, "seed",{{CL_INT_ARG,&seed,"seed",NULL}}},
	CL_END_OPTIONS  /* defined as {0,NULL,..} to mark the end of the list   */
  }; 


struct BAYES_UPDATE
{
  double **C_tilde;         /* the zero-offset posterior covariance matrix */
  double *mu_tilde;        /* the zero-offset posterior mean vector */
  int info;                /* Success flag: 
							  0 = successful, 
							  1 = failed Cp inverse, 
							  2 = failed Cy inverse
							  3 = failed Cp^{-1}+X^T.Cy^{-1}.X inverse
						   */
};


/* 
   Experimental file for downscaling.... not part of main source
*/
int main(int argc, char *argv[])
{ 
  void choldc_double(double **a, int n);
  void writeOut(int Nx,int Nz,double *th);
  struct BAYES_UPDATE updateBayes(double **Cp, double *mu, int Np, double **X, double **Cy, double *y, int Ny, 
								  int CyIsAlreadyInverted, int CyIsIdentityMatrix);

  long idum=-1;
  CLContext   context;
  struct BAYES_UPDATE bup;

  /*
	Initialise the command line context variable. 
	The third 'string' parameter is printed with the usage message
	to indicate what arguments, excluding the options and their arguments,
	are required on the command line.
  */
  
  clInit( &context, options,ARGS_USAGE, NULL );
        
  /*  Parse the command line and if successful, print out its contents.
	  If not successful the function prints out the usage message.
  */
  
  if (clParse( &context, argc, argv ))
    {
	  if (context.argc != 1) 
        {
		  fprintf(stderr, "\n");
		  fprintf(stderr, "Missing some command line arguments\n");
		  clPrintUsage(&context);		  
		  return EXIT_FAILURE;
        }
	  clDone( &context );
	  idum = -(long)fabs(seed);
    } 
  else
    {
	  /*  Unsuccessful parse of the command line. */
	  return EXIT_FAILURE;
    }


  int i,ii,i_r,j,k,kk,l,Nm;
  int Nx, Nz;

  double sum,C0,lx,lz,l_xsq,l_zsq;
  double *NG,*th,**C,*thBar,*w;

  Nx=50;
  Nz=20;
  C0 = 1*1;
  lx=10;
  lz=0.1;

  NG = dvector(0,Nx-1);
  for (i=0;i<Nx;i++) NG[i]= DMAX(0,DMIN(1.0, 0.2+(1.0*i)/Nx));

  Nm=Nx*Nz;
  thBar = dvector(0,Nm-1);
  th = dvector(0,Nm-1);
  C=dmatrix(0,Nm-1,0,Nm-1);
  w = dvector(0,Nm-1);
  l_xsq = lx*lx;
  l_zsq = lz*lz;
  for (k=0,j=0;k<Nz;k++)
	{
	  for (i=0;i<Nx;i++,j++)
		{
		  thBar[j] = NG[i];
		  for (kk=0,l=0;kk<Nz;kk++)
			{
			  for (ii=0;ii<Nx;ii++,l++)
				{
				  C[j][l] = (((j==l)? 1e-6 :0) + C0)*exp(-((i-ii)*(i-ii)/l_xsq + (k-kk)*(k-kk)/l_zsq));
				}
			}
		}
	}
  // Take copy for MCMC later
  double **Cp = dmatrix(0,Nm-1,0,Nm-1);
  memcpy(&Cp[0][0],&C[0][0],Nm*Nm*sizeof(double));

  // get Cholesky of C
  choldc_double(C,Nm);
  // Example: simple realisation from prior
  if (0)
	{ 
	  for (i_r=0;i_r<=0;i_r++)
		{
		  for (i=0;i<Nm;i++,l++) 
			{
			  w[i] = Gaussian_deviate(&idum);
			  for (j=0,sum=0;j<=i;j++) sum += C[i][j]*w[j];
			  th[i] = sum + thBar[i];
			  th[i] = DMAX(0,th[i]);
			}
		  writeOut(Nx,Nz,th);
		}
	  exit(0);
	}

  // MCMC sampling:
  double scale,chisqp,chisqp_old,chisq_lik,lik_old=0,sigma_NG;
  double **Lpost;

  // get L^{-1} of prior
  modified_Cholesky_left_inverse(Nm,C);
  double **Linv = C;

  // Desired accuracy of NG match 
  sigma_NG = .1;


  double *dth = dvector(0,Nm-1);
  double *thp = dvector(0,Nm-1);
  double *D_thp = dvector(0,Nm-1);



  // form effective "design" matrix
  double **X = dmatrix(0,Nx-1,0,Nm-1);
  memset(&X[0][0],0,Nx * Nm * sizeof(double));

  for (k=0,j=0;k<Nz;k++)
	{
	  for (i=0;i<Nx;i++,j++)
		{
		  X[i][j] = 1.0/Nz; 
		} 
	}
  // print_dmatrix0(X,Nx,Nm,"X.txt");

  // form effective "Cy" matrix
  double **Cy = dmatrix(0,Nx-1,0,Nx-1);
  for (i=0;i<Nx;i++)
	{
	  for (j=0;j<Nx;j++) Cy[i][j] = (i==j)? 1.0/(sigma_NG*sigma_NG) : 0;
	}

  // print_dvector0(thBar,Nm,"thBar.txt");

  double *thBarNew = dvector(0,Nm-1);
  for (i=0;i<Nm;i++) thBarNew[i] = thBar[i];

  // get posterior system
  bup = updateBayes(Cp,thBarNew,Nm,X, Cy,NG,Nx,1,0);
  printf("Bayes update info = %d\n",bup.info);


  print_dvector0(thBar,Nm,"thpBar.txt");

  // get Cholesky factor of updated Cp
  choldc_double(Cp,Nm);
  Lpost = Cp;

  // simple realisations from approximate (linearised) Bayes posterior
  if (0)
	{
	  for (i_r=0;i_r<=5;i_r++)
		{
		  for (i=0;i<Nm;i++,l++) 
			{
			  w[i] = Gaussian_deviate(&idum);
			  for (j=0,sum=0;j<=i;j++) sum += Lpost[i][j]*w[j];
			  th[i] = sum + thBarNew[i];
			  th[i] = (th[i]<0)? 0 : th[i];
			}

		  // likelihood check
		  chisq_lik=0;
		  for (i=0;i<Nx;i++)
			{
			  // total column sum of thickness
			  for (k=0,sum=0;k<Nz;k++)
				{
				  sum += th[i+(k*Nx)];
				}
			  sum/=Nz;  // puts on NG scale
			  chisq_lik += DSQR((sum-NG[i])/sigma_NG);
			}
		  chisq_lik/=2.0;
		  printf("Sample -logLik = %lg\n",chisq_lik);
		  writeOut(Nx,Nz,th);
		}
	  exit(0);
	}


  // start at mean 
  for (i=0;i<Nm;i++) th[i] =  thBarNew[i];
  chisqp_old=0;
  scale=2.0/sqrt(Nx*Nz);

  int Ntry=10000;
  l=0;
  for (i_r=0;i_r<Ntry;i_r++)
	{
	  // proposals from Bayes posterior L
	  // dth = scale*L.w
	  // thp = th + dth
	  // D_thp = thp-thBar
	  for (i=0;i<Nm;i++) 
		{
		  w[i] = Gaussian_deviate(&idum);
		  for (j=0,sum=0;j<=i;j++) sum += Lpost[i][j]*w[j];
		  dth[i] = scale*sum;
		  thp[i] = th[i] + dth[i];
		  D_thp[i] = thp[i] - thBar[i];
		}
	  // chisq_prior = || Linv.(thp-thbar) ||_2
	  for (i=0,chisqp=0;i<Nm;i++) 
		{
		  for (j=0,sum=0;j<=i;j++) sum += Linv[i][j]*D_thp[j];
		  chisqp += sum*sum;
		}
	  chisqp/=2;

	  if (0)
		{
		  // hacky rescaling to approx. satisfy likelihood....
		  // NB: breaks reversibility 
		  for (i=0;i<Nx;i++)
			{
			  // total column sum of thickness
			  for (k=0,sum=0;k<Nz;k++)
				{
				  sum += DMAX(0,thp[i+(k*Nx)]);
				}
			  sum/=Nz;  // puts on NG scale
			  float deviate = sigma_NG*Gaussian_deviate(&idum);
			  for (k=0;k<Nz;k++)
				{
				  thp[i+(k*Nx)] *= (1+deviate)*NG[i]/sum;
				}			  
			}
		  // recompute prior....
		  for (i=0;i<Nm;i++) 
			{
			  D_thp[i] = thp[i] - thBar[i];
			}
		  // chisq_prior = || Linv.(thp-thbar) ||_2
		  for (i=0,chisqp=0;i<Nm;i++) 
			{
			  for (j=0,sum=0;j<=i;j++) sum += Linv[i][j]*D_thp[j];
			  chisqp += sum*sum;
			}
		  chisqp/=2;
		}

	  // likelihood
	  chisq_lik=0;
	  for (i=0;i<Nx;i++)
		{
		  // total column sum of thickness
		  for (k=0,sum=0;k<Nz;k++)
			{
			  sum += DMAX(0,thp[i+(k*Nx)]);
			}
		  sum/=Nz;  // puts on NG scale
		  chisq_lik += DSQR((sum-NG[i])/sigma_NG);
		}
	  chisq_lik/=2.0;
	  chisqp += chisq_lik;

	  // MH step
	  if (random1(&idum) < exp(-(chisqp-chisqp_old)) || (i_r==0))
		{
		  // printf("ACC: %d 1\n",i_r);
		  // accept
		  chisqp_old = chisqp; 
		  lik_old=sum;
		  memcpy(th,thp,Nm*sizeof(double));
		  l++;
		}
	  else 
		{
		  // reject
		}
	  // printf("MTH: %d %f\nMTH: %d %f\n",i_r,th[0],i_r,th[1]);
	  //	  printf("MTH: %d %g\n",i_r,chisqp_old);

	  if (i_r%100 == 0 ) 
		{
#if 1
		  // example total NG for column i=8   (i=0...Nx-1)
		  i=0;
		  {
			// total column sum of thickness
			for (k=0,sum=0;k<Nz;k++)
			  {
				sum += DMAX(0,thp[i+(k*Nx)]);
			  }
			sum/=Nz;  // p
		  }
		  printf("MTH: %d %lf\n",i_r,thp[0+(0*Nx)]);
#endif
		  // example total NG for block (0,0)
		  // printf(": %d %lf\n",i_r,DMAX(0,thp[8+(0*Nx)]));

		  writeOut(Nx,Nz,th);
		}
	}
  printf("Accept rate %f\n",l/(float)Ntry);
  return 0;
}

// Use plotutils: graph -TX  -q 1 -B realisation.%d.txt

void writeOut(int Nx,int Nz,double *th)
{
  int i,j,k;
  FILE* f;
  char *s;
  double* accHeight = dvector(0,Nx);
  memset(&accHeight[0],0,Nx*sizeof(double));

  asprintf(&s,"realisation.%d.txt",real_no++);

  f = fopen(s,"w");
  for (k=0,j=0;k<Nz;k++)
	{
	  for (i=0;i<Nx;i++,j++)
		{
		  // th[j] = thickness of (i,k) slab
		  // draw little trapezoid block
#if 1
		  if (k%2==0)
			{
			  if (i<(Nx-1))
				{
				  //				  if (th[j]>0 || th[j+1]>0)
					{
					  fprintf(f,"%d %lg\n",i,accHeight[i]);
					  fprintf(f,"%d %lg\n",i,accHeight[i]+DMAX(th[j],0));
					  fprintf(f,"%d %lg\n",i+1,accHeight[i+1]+DMAX(th[j+1],0));
					  fprintf(f,"%d %lg\n",i+1,accHeight[i+1]);
					  fprintf(f,"%d %lg\n\n",i,accHeight[i]);
					}
				}
			  else
				{
				  // if (th[j]>0)
					{
					  fprintf(f,"%d %lg\n",i,accHeight[i]);
					  fprintf(f,"%d %lg\n",i,accHeight[i]+DMAX(th[j],0));
					  fprintf(f,"%f %lg\n",i+0.1,accHeight[i]+DMAX(th[j],0));
					  fprintf(f,"%f %lg\n",i+0.1,accHeight[i]);
					  fprintf(f,"%d %lg\n\n",i,accHeight[i]);
					}
				}
			}
		  accHeight[i] += DMAX(0,th[j]);

#else

		  if (i<(Nx-1))
			{
			  if (th[j]>0 || th[j+1]>0)
				{
				  fprintf(f,"%d %lg\n",i,k-th[j]/2.0);
				  fprintf(f,"%d %lg\n",i,k+th[j]/2.0);
				  fprintf(f,"%d %lg\n",i+1,k+th[j+1]/2.0);
				  fprintf(f,"%d %lg\n",i+1,k-th[j+1]/2.0);
				  fprintf(f,"%d %lg\n\n",i,k-th[j]/2.0);
				}
			}
		  else
			{
			  if (th[j]>0)
				{
				  fprintf(f,"%d %lg\n",i,k-th[j]/2.0);
				  fprintf(f,"%d %lg\n",i,k+th[j]/2.0);
				  fprintf(f,"%f %lg\n",i+0.1,k+th[j]/2.0);
				  fprintf(f,"%f %lg\n",i+0.1,k-th[j]/2.0);
				  fprintf(f,"%d %lg\n\n",i,k-th[j]/2.0);
				}
			}
#endif
		}
	}
  fclose(f);
  printf("Wrote %s\n",s);
  free(s);
}

void choldc_double(double **a, int n)
{
	void nrerror(char error_text[]);
	int i,j,k;
	double sum,*p;

	p = dvector(0,n-1);

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					nrerror("choldc failed");
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}

	for (i=0;i<n;i++) 
	  {
		a[i][i] = p[i];
		// clean upper tri
		for (j=i+1;j<n;j++) a[i][j]=0;
	  }

	free_dvector(p,0,n-1);
}


/*
  Compute a standard Bayesian update of a model m with 
  multi-Gaussian prior m ~ N(mu,Cp), of dimension Np.

  mu[0...Np-1] is the prior mean
  Cp[0...Np-1][0...Np-1] is the symmetric, positive definite, prior covariance  

  The prior is updated with a negative log-likelihood of form

  (y-X.m)^T.Cy^{-1}.(y-X.m)/2

  where there is a vector of "observations" 

  y[0...Ny-1], 

  which are mutually correlated by the symmetric positive  definite matrix 

  Cy[0...Ny-1][0....Ny-1]

  Some points:
  (1) We often have the trivial case Cy=I, and this can be flagged to increase efficiency
      (set CyIsIdentityMatrix =1).
  (2) Also, we often have the trivial case Cy= diagonal matrix, in which case the user should invert it
      explicitly, pass it in and set CyIsAlreadyInverted=1.
      Again, this increases efficiency significantly.

  In this linear Bayes model, the observations y are modelled by the linear form X.m
  and

  X[0...Ny-1][0....Np-1]  

  is the "design" matrix. This function computes  (in latex-ese)
  
  C_tilde = (Cp^{-1}+X^T.Cy^{-1}.X)^{-1}

  mu_tilde = C_tilde.(Cp^{-1}.mu + X^{T}.Cy^{-1}.y)

  and these appear in the returned BAYES_UPDATE 
  structure in the same storage space as the original Cp, mu

  To track problems, we return an info flag with these meanings:
  info =
  0 = successful, 
  1 = failed Cp inverse, 
  2 = failed Cy inverse
  3 = failed Cp^{-1}+X^T.Cy^{-1}.X inverse

  N.B. blas usage:
  using CblasRowMajor, any lda type parameter is the number of columns in the matrix
  e.g. for C[p,q], ldc=q.

  N.B. The input matrix Cp and vector mu are overwritten *in-place* for the elements of the returned
  structure struct BAYES_UPDATE update
  
  I.e. if you want to retain Cp and mu, make copies before calling this function.


*/
struct BAYES_UPDATE updateBayes(double **Cp, double *mu, int Np, double **X, double **Cy, double *y, int Ny, 
								int CyIsAlreadyInverted, int CyIsIdentityMatrix)
{
  int i,j,info;
  struct BAYES_UPDATE update;
  double* tmp;
  double **Xtmp;
  double *y_tmp, *mu_tmp;

  /* Assume update OK till contradicted */
  update.info = 0;

  /* Invert Cp using Cholesky. I think Cp only needs lower triangle in place */
  info = clapack_wrapper_dpotrf('L',Np,&Cp[0][0],Np);
  if (info!=0)
	{
	  /* failure of Cholesky */
	  update.info = 1;
	  return update;
	}
  info = clapack_wrapper_dpotri('L',Np,&Cp[0][0],Np);
  /* put top triangle correvtly there */
  for (i=0; i<Np; i++)
	for (j=0;j<i;j++)
	  Cp[j][i]=Cp[i][j];
  /* Cp now has its inverse in place */

  /* form mu_tmp = CP^{-1}.mu */
  /*  void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
  */
  mu_tmp=dvector(0,Np-1);
  cblas_wrapper_dgemv('N', Np, Np, 1.0, &Cp[0][0], Np, &mu[0], 1, 0.0, &mu_tmp[0],1);

  /* Form Cy^{-1} */
  if (CyIsIdentityMatrix != 1)
	{
	  if (CyIsAlreadyInverted != 1)
		{
		  /* Invert Cy using Cholesky. I think Cy only needs lower triangle in place */
		  info = clapack_wrapper_dpotrf('L',Ny,&Cy[0][0],Ny);
		  if (info!=0)
			{
			  /* failure of Cholesky */
			  update.info = 2;
			  return update;
			}
		  info = clapack_wrapper_dpotri('L',Ny,&Cy[0][0],Ny);
		}
	  for (i=0; i<Ny; i++)
		for (j=0;j<i;j++)
		  Cy[j][i]=Cy[i][j];

	  /* Cy now has its inverse in place, and its not a trivial identity matrix */
	  /* form y_tmp = Cy^{-1}.y */
	  y_tmp=dvector(0,Ny-1);
	  cblas_wrapper_dgemv('N', Ny, Ny, 1.0, &Cy[0][0], Ny, &y[0], 1, 0.0, &y_tmp[0],1);
	}
  else
	{
	  y_tmp = y;
	}

  /* Form Xtmp = Cy^{-1}.X = {Ny x Np} matrix */
  if (CyIsIdentityMatrix != 1)
	{	
	  /*
		void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
	  */
	  /* Need new storage for intermediate entity */
	  Xtmp = dmatrix(0,Ny-1,0,Np-1);
	  cblas_wrapper_dgemm('N', 
                 'N',Ny,Np,
                 Ny,1.0,&Cy[0][0],
                 Ny,&X[0][0],Np,
                 0.0,&Xtmp[0][0],Np);
	}
  else
	{
	  Xtmp=X;
	}
	  
  /* now form C_tilde^{-1} = Cp^{-1}+X^T.Xtmp */
  cblas_wrapper_dgemm('T',
			  'N',Np,Np,
			  Ny,1.0,&X[0][0],
			  Np,&Xtmp[0][0],Np,
			  1.0,&Cp[0][0],Np);

    
  /* clean up Xtmp as quickly as we can if we need to */
  if (CyIsIdentityMatrix != 1)
	{	
	  free_dmatrix(Xtmp,0,Ny-1,0,Np-1);
	}
  /* invert C_tilde^{-1} back to C_tilde */
  info = clapack_wrapper_dpotrf('L',Np,&Cp[0][0],Np);
  if (info!=0)
	{
	  /* failure of Cholesky */
	  update.info = 3;
	  return update;
	}
  info = clapack_wrapper_dpotri('L',Np,&Cp[0][0],Np);
  for (i=0; i<Np; i++)
	for (j=0;j<i;j++)
	  Cp[j][i]=Cp[i][j];
  /* OK, now we have C_tilde */
  update.C_tilde = Cp;
  

  /* get updated mean */
  /* form mu_tmp + X^T.y_tmp */
  /*  void cblas_dgemv(const enum CBLAS_ORDER Order,
	             const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
	             const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
  */
  cblas_wrapper_dgemv( 'T', Ny, Np, 1.0, &X[0][0], Np, &y_tmp[0], 1, 1.0, &mu_tmp[0],1);
  /* Answer is now in mu_tmp */

  /* now do mu <- Ctilde.mu_tmp to get updated mean */
  cblas_wrapper_dgemv( 'N', Np, Np, 1.0, &Cp[0][0], Np, &mu_tmp[0], 1, 0.0, &mu[0],1);

  if (CyIsIdentityMatrix != 1)
	{
	  free_dvector(y_tmp,0,Ny-1);
	}
  free_dvector(mu_tmp,0,Np-1);

  update.mu_tilde = mu;

  return update;
}
