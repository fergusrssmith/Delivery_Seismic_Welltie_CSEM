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
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "NR/nrutil.h"

#include "util.h"
#include "su_util.h"
#include "upperTriangleMatrices.h"

#define LINE_LENGTH 256

/*
  Simple program for translating ascii analyse-out files to compressed binary form
  Usage:

  analyseOutputCompacter massage-analysis-ascii-file.txt [ - for stdin] > massage-analysis-ascii-file.mab

*/



extern int errno;



int main(int argc, char *argv[])
{
  char           *filename;
  FILE           *in;
  char    line_buffer[LINE_LENGTH];
  int i,j,k,nMatches=0,number_of_columns;
  int nNonCovarianceColumns;
  int n_cov_entries=0,n_prop_entries=0;
  float f;
  short f_sh;
  int goodData;

  // binary extraction
  int is_binary=0;
  int machine_is_little_endian = 1;
  int short_scale = 32767;  // java standard MAX_SHORT


  if (argc!=2)
	{    
	  printf("\nSimple program for translating \"deliveryAnalyser --massage-analyse-ascii ...\"\n"
			 "ascii output files to compressed binary form for use with deliveryMassager\n\n"
			 "Usage: analyseOutputCompacter massage-analysis-ascii-file.txt [ - for stdin] > "
			 "massage-analysis-ascii-file.mab\n\n");
	  exit(1);
	}  

  filename = argv[1];

  if (!strcmp(filename,"-"))
	{
	  in = stdin;
	  fprintf(stderr,"Reading Massage-analysis data from stdin\n");
	}
  else
	{
	  if ((in = fopen(filename, "r")) == NULL)
		{
		  fprintf(stderr,"cannot open file %s\n",filename);
		  fprintf(stderr,"%s\n",strerror(errno));
                       
		  return EXIT_FAILURE;
		}
	  fprintf(stderr,"Reading Massage-analysis data from %s\n",argv[1]);
	}

 
  // Parse:
  // title
  fgets(line_buffer,LINE_LENGTH,in);

  // binary test
  if (!strncmp(line_buffer,"Binary-compressed",17))
	{
	  is_binary = 1;
	  i = 1;
	  if (*(char *)&i == 1) machine_is_little_endian = 1;
	  else machine_is_little_endian = 0;
	}
  printf("Binary-compressed [corr_coeff_scale= %d n_bytes_covariance= 2] "
		 "Massaging analysis summary file\n",short_scale);

  if (fscanf(in,"%d\n",&number_of_columns) != 1)
	{
	  fprintf(stderr,"Can't read no. of columns in file %s\n",filename);
	  return EXIT_FAILURE;
	}
  printf("%d\n",number_of_columns);

  /*
	Count columns titles and compute number of properties in file, plus number of
	covariance entries
  */

  for (i=0;i<number_of_columns;i++)
	{
	  if (fscanf(in,"%s",line_buffer)!=1)
		{
		  fprintf(stderr,"Can't read column title %s in file %s\n",line_buffer,filename);
		  return EXIT_FAILURE;
		}
	  printf("%s\n",line_buffer);
	  // count entries that look like "C_xx"
	  if (!strncmp(line_buffer,"C_",2)) n_cov_entries++;
	}
  // since the file has the upper triangle of the covariance
  // we must have n_cov_entries = n_prop_entries*(n_prop_entries+1)/2
  // Compute n_prop_entries.
  for (n_prop_entries=0; (n_prop_entries*(n_prop_entries+1))/2 < n_cov_entries; n_prop_entries++);

  nNonCovarianceColumns = number_of_columns-n_cov_entries;

  // ptrs to data for parsing convenience
  float** covar = matrix(0,n_prop_entries-1,0,n_prop_entries-1);
  float* sqrt_covar_diag = vector(0,n_prop_entries-1);
  float* non_covariance_entries = vector(0,nNonCovarianceColumns-1);

  for (;;)
	{
	NEWPOINT:
	  // properties
	  goodData = 1;
	  for (i=0; i<nNonCovarianceColumns; i++)
		{
		  if ((nMatches = fscanf(in," %f ",&f)) == 1)
			{
			  if (isnan((double)f)) goodData=0;
			  // write out in big-endian
			  // if (i==0) fprintf(stderr,"BINOUT %g\n",f);
			  non_covariance_entries[i] = f;
			  // fwrite(&f,4,1,stdout);
			}
		  else
			{
			  fprintf(stderr,"Finished processing\n");
			  goto FINISH;
			}
		}
	  // read and save covariances. Last columns are always covariances
	  for (j=0;j<n_prop_entries;j++)
		{
		  for (k=j; k<n_prop_entries; k++)
			{
			  if (fscanf(in," %f ",&f)==1 )
				{
				  covar[j][k]=f;
				  if (isnan((double)f)) goodData=0;
				}
			  else
				{
				  fprintf(stderr,"Can't read covariance number in file");
				  exit(1);
				}
			}
		}

	  if (!goodData)
		{
		  fprintf(stderr,"Warning! NaN's in data at location %f %f; omitting trace!\n",
				  non_covariance_entries[0],non_covariance_entries[1]);
		}
	  else
		{
		  for (j=0;j<n_prop_entries;j++) sqrt_covar_diag[j] = sqrt(covar[j][j]);

		  // if ((non_covariance_entries[0] > 1790)) goto NEWPOINT;
		  // write out non-covariances in binary 
		  for (i=0; i<nNonCovarianceColumns; i++)
			{
			  if (machine_is_little_endian) swap_float_4(&non_covariance_entries[i]);
			  fwrite(&non_covariance_entries[i],4,1,stdout);
			}

		  // write out covariances in binary
		  for (j=0;j<n_prop_entries;j++)
			{
			  for (k=j; k<n_prop_entries; k++)
				{
				  if (j==k)
					{
					  f = covar[j][k];
					  // write 2 byte float for diagonals
					  if (machine_is_little_endian) swap_float_4(&f);
					  fwrite(&f,2,1,stdout);
					}
				  else
					{
					  // scale etc
					  if ((covar[j][j]>0) && (covar[k][k]>0))
						{
						  f = sqrt_covar_diag[j] * sqrt_covar_diag[k];
						  f_sh =(short)(rint(short_scale*covar[j][k]/f));
						}
					  else f_sh=0;
					  if (machine_is_little_endian) swap_short_2(&f_sh);
					  fwrite(&f_sh,2,1,stdout);
					}
				}
			}
		}
	}

 FINISH:

  if (strcmp(filename,"-")) fclose(in);
  free_matrix(covar,0,n_prop_entries-1,0,n_prop_entries-1); covar=NULL;
  free_vector(sqrt_covar_diag,0,n_prop_entries-1);sqrt_covar_diag=NULL;
  free_vector(non_covariance_entries,0,nNonCovarianceColumns-1);non_covariance_entries=NULL;

  return 0;
}
