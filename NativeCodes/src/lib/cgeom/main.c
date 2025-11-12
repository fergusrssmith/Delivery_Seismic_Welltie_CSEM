/*** MAIN.C ***/
/*
 example compile

 gcc -Wall -Wno-unused  -g -p -DDEBUG -I.. -I../../../include -L./Linux-i686-unknown main.c -lcgeom -lm -o voronoiTest 
 ./voronoiTest | gr -B -x -0.2 1.2 -y -0.2 1.2 -C --bitmap-size=1000x1000


*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h> /* isatty */

#include "vdefs.h" 


#define NPOINTS 300
#define NMAX 10000
#define LINE_LENGTH 200

int main(int argc, char *argv[])
{
  struct timeval tv;
  struct timezone tz;
  double told;
  int i,n;
  float *x,*y,dist;
  struct VORONOI_INFO info;
  int triangulate = 0;
  int showlattice=0;
  int dont_show_points=0;
  int random=0;
  int uniform=0;
  int help=0;
  int readstdin=0;
  int nplot=-1;
  float eps_noise=0;

  n = NPOINTS;

  for (i=1;i<argc;i++)
	{
	  if (!strcmp(argv[i],"-r")) random=1;
	  if (!strcmp(argv[i],"-u")) uniform=1;
	  if (!strcmp(argv[i],"-t")) triangulate=1;
	  if (!strcmp(argv[i],"-l")) showlattice=1;
	  if (!strcmp(argv[i],"-d")) dont_show_points=1;
	  if (!strcmp(argv[i],"-h")) help=1;
	  if (!strcmp(argv[i],"-e")) 
		{
		  if ((++i) < argc)
			{
			  eps_noise = atof(argv[i]);
			  if (eps_noise<0)
				{
				  printf("-e eps must supply noise level eps >= 0\n");
				  help=1;
				}			  
			}
		  else
			{
			  printf("option -e must supply number\n");
			  help=1;
			}
		}
	  if (!strcmp(argv[i],"-n")) 
		{
		  if ((++i) < argc)
			{
			  n = atoi(argv[i]);
			  if ((n<1) || (n>NMAX))
				{
				  printf("-n N must supply number N>0 & N<%d\n",NMAX);
				  help=1;
				}
			  
			}
		  else
			{
			  printf("option -n must supply number\n");
			  help=1;
			}
		}
	  if (!strcmp(argv[i],"-p")) 
		{
		  if ((++i) < argc)
			{
			  nplot = atoi(argv[i]);
			  if ((nplot<1) || (nplot>NMAX))
				{
				  printf("-p N must supply number N>0 & N<%d\n",NMAX);
				  help=1;
				}			  
			}
		  else
			{
			  printf("option -p must supply number\n");
			  help=1;
			}
		}
	  if (!strcmp(argv[i],"-")) readstdin=1;
	}



  if (!readstdin && !(random||uniform)) help=1; // don't leave stdin hanging 
	
  if (help)
	{
	  printf("usage: voronoiTest [-ruh]\n");
	  printf("       -r for internally generated random test\n");
	  printf("       -u for internally generated uniform points test (plus noise wiggle)\n");
	  printf("       - read 2-column ascii x y input points from stdin\n");
	  printf("       -t for Delaunay triangulation\n");
	  printf("       -l to show nearest regular lattice (imposes Delaunay)\n");
	  printf("       -d to suppress symbols at points\n");
	  printf("       -e noise       for noise level on uniform grid\n");
	  printf("       -n number      for number of test points\n");
	  printf("       -p number      for number of Delaunay triangles to plot \n");
	  printf("       -h this help\n");
	  exit(0);
	}

  if (showlattice) triangulate=1;

  gettimeofday(&tv,&tz);
  srand((int)(tv.tv_sec));

  if (random)
	{
	  // random points
	  x = (float*)malloc(n*sizeof(float));
	  y = (float*)malloc(n*sizeof(float));
	  for (i=0;i<n;i++)
		{
		  x[i]=rand()/(float)RAND_MAX;
		  y[i]=rand()/(float)RAND_MAX;
		}
	}
  else if (uniform)
	{
	  int j,k,m;
	  float X,Y,Xr,Yr,stretchX,stretchY,ep,theta,cs,sn;
	  m=(int)sqrt(n);
	  ep = eps_noise/m;

	  x = (float*)malloc(n*sizeof(float));
	  y = (float*)malloc(n*sizeof(float));
 
	  // regular guys 
	  stretchX = 2.3;
	  stretchY = 1;
	  // spin about
	  theta = 3.14*rand()/((float)RAND_MAX);
	  cs = cos(theta);
	  sn = sin(theta);
	  for (i=0,k=0; (k<n); i++) 
		{
		  for (j=0;(j<m) && (k<n);j++,k++)
			{
		   	  X=stretchX * (((float)i)/(float)m + ep*rand()/(float)RAND_MAX);			
			  Y=stretchY * (((float)j)/(float)m + ep*rand()/(float)RAND_MAX);
			  Xr = cs*X - sn*Y;
			  Yr = sn*X + cs*Y;
			  x[k] = Xr;
			  y[k] = Yr;
			}
		} 
	}
  else 
	{
	  fprintf(stderr,"Expect stdin ascii pairs...\n");
	  // stdin ascii
	  char line_buffer[LINE_LENGTH];
	  int nAlloc;
	  n = nAlloc = NMAX;
	  x = (float*)malloc(NMAX*sizeof(float));
	  y = (float*)malloc(NMAX*sizeof(float));

	  for (i=0;;i++)
		{
		  if (!(fgets(line_buffer , LINE_LENGTH , stdin)))
			{
			  break;
			}
		  else
			{
			  if (i >= nAlloc)
				{
				  nAlloc*=2;
				  x = (float*)realloc(x,nAlloc*sizeof(float));
				  y = (float*)realloc(y,nAlloc*sizeof(float));
				}
			  if (sscanf(line_buffer,"%g %g",&(x[i]),&(y[i])) == 2) 
				{
				  /* eliminate duplicate x,y's */
				  if ((i>0) && (x[i] == x[i-1]) && (y[i] == y[i-1]))
					{
					  i--; /* read next point into duplicate location */
					}
				  n=i+1;
				}
			  else break;
			}
		}
	  // for (i=0;i<n;i++)  printf("%g %g\n",x[i],y[i]);
	  // exit(0);
	  fprintf(stderr,"read %d points\n",n);
	}


  gettimeofday(&tv,&tz);
  told=(tv.tv_sec*1e6 + tv.tv_usec);

  info = getdiagram(0,triangulate,n,'f',x,y);

  gettimeofday(&tv,&tz);
  fprintf(stderr,"Size= %d,Voronoi time: %lf uS\n",n,(double)((tv.tv_sec*1e6+tv.tv_usec) - told));

  if (0)
	{
	  /* graph -TX friendly voronoi diagram */
	  for (i=0;i<info.nEdges;i++)
		{
		  /* edges are sorted into vertex order by now */
		  if ((info.edges[i].i != -1) && (info.edges[i].j != -1))
			{
			  printf("%f %f\n",info.vertices[info.edges[i].i].x,info.vertices[info.edges[i].i].y);
			  printf("%f %f\n",info.vertices[info.edges[i].j].x,info.vertices[info.edges[i].j].y);
			  printf("\n");
			}
		} 
	} 
  /* method of traversing the cells */

  if (0) 
	{
	  int j,s;
	  for (i=0;i<n;i++)
		{
		  // loop over edges surrounding site
		  for (j=0;j<info.siteEdgesSize[i];j++)
			{
			  struct EDGE *e = &(info.edges[info.siteEdges[i][j]]);
			  if ((e->i != -1) && (e->j != -1))
				{
				  printf("%f %f\n",info.vertices[e->i].x,info.vertices[e->i].y);
				  printf("%f %f\n",info.vertices[e->j].x,info.vertices[e->j].y);
				  printf("\n");
				}
			}
		  printf("\n");
		}
	}
  /* in order.... */
  struct EDGE *e,*e2;
  int j,v1,v2;

 
  if (!triangulate)
   {
	 gettimeofday(&tv,&tz);
	 told=(tv.tv_sec*1e6 + tv.tv_usec);
	 orderVoronoiPolygons(&info);
	 dist = computeMeanSiteSeparation(&info);  
	 fprintf(stderr,"RMS separation = %g\n",dist);
	 gettimeofday(&tv,&tz);
	 fprintf(stderr,"Size= %d,Voronoi sort-edge-order time: %lf uS\n",n,(double)((tv.tv_sec*1e6+tv.tv_usec) - told));
   }

  struct NEAREST_REGULAR_GRID_INFO RGinfo;

  if (showlattice)
	{
	  RGinfo = nearest_regular_lattice(&info,1,n,'f',x,y,0.05,.05);
	  if (RGinfo.info != 0) fprintf(stderr,"Nearest grid info=%d\n",RGinfo.info);
	}

  if (!dont_show_points)
	{
	  // print points
	  printf("#m=0,S=17\n"); // asterisk, will reveal circles or squares behind
	  for (i=0;i<n;i++)  printf("%f %f\n",x[i],y[i]);
	}

  // Now Delaunay or Voronoi
  printf("\n#m=2,S=0\n");
  if (triangulate)
	{
	  if (nplot==-1) 
		{
		  j = info.nTriangles;
		}
	  else 
		{
		  j = (nplot<info.nTriangles)?  nplot :  info.nTriangles;
		}
	  for (i=0;i<j;i++) 
		{ 
		  int * sv = info.triangleSites[i];
		  printf ("%f %f\n", x[sv[0]], y[sv[0]]);
		  printf ("%f %f\n", x[sv[1]], y[sv[1]]);
		  printf ("%f %f\n", x[sv[2]], y[sv[2]]);
		  printf ("%f %f\n\n", x[sv[0]], y[sv[0]]);
		} 
	} 
  else
	{
	  for (i=0;i<n;i++) 
		{
		  // loop over edges surrounding site
		  for (j=0;j<info.siteEdgesSize[i];j++)
			{
			  v1 = info.edgeLvertex[i][j];
			  v2 = info.edgeLvertex[i][(j+1)%(info.siteEdgesSize[i])];
			  if ((v1 != -1) && (v2 != -1))
				{
				  printf("%f %f\n",info.vertices[v1].x,info.vertices[v1].y);
				  printf("%f %f\n",info.vertices[v2].x,info.vertices[v2].y);
				}
			}
		  printf("\n");
		}
	}
  /**/
  free_voronoi_info(&info);

  return 0;
}
