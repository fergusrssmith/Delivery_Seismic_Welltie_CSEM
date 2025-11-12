#include <stdio.h>
#include <stdlib.h>

#include "vdefs.h"

extern int _triangulate, plot, _debug ;
extern float ymax, ymin, xmax, xmin ;

float pxmin, pxmax, pymin, pymax, cradius;

void out_bisector(Edge * e)
{ 
#if 0
  if (!_triangulate && !_debug)
	{
	  printf("l %f %f %f\n", e->a, e->b, e->c) ;
	}
#endif
  if (_debug)
	{
	  printf("line(%d) %gx+%gy=%g, bisecting %d %d\n", e->edgenbr,
			 e->a, e->b, e->c, e->reg[le]->sitenbr, e->reg[re]->sitenbr) ;
	}
}

void out_ep(Edge * e)
{
#if 0
  if (!_triangulate)
	{
	  printf("e %d", e->edgenbr);
	  printf(" %d ", e->ep[le] != (Site *)NULL ? e->ep[le]->sitenbr : -1) ;
	  printf("%d\n", e->ep[re] != (Site *)NULL ? e->ep[re]->sitenbr : -1) ;
	}
#endif
}

void out_vertex(Site * v)
{
#if 0
  if (!_triangulate && !_debug)
	{
	  printf ("v %f %f\n", v->coord.x, v->coord.y) ;
	}
#endif
  if (_debug)
	{
	  printf("vertex(%d) at %f %f\n", v->sitenbr, v->coord.x, v->coord.y) ;
	}
}

void out_site(Site * s)
{
#if 0
  if (!_triangulate && !_debug)
	{
	  printf("s %f %f\n", s->coord.x, s->coord.y) ;
	}
#endif
  if (_debug)
	{
	  printf("site (%d) at %f %f\n", s->sitenbr, s->coord.x, s->coord.y) ;
	}
}

void out_triple(Site * s1, Site * s2, Site * s3,struct VORONOI_INFO *info)
{
#if 1
  if (_triangulate && !_debug)
	{
	  // printf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr) ;
	  // JG hack
	  /*
		printf ("%f %f\n", s1->coord.x, s1->coord.y);
		printf ("%f %f\n", s2->coord.x, s2->coord.y);
		printf ("%f %f\n", s3->coord.x, s3->coord.y);
		printf ("%f %f\n\n", s1->coord.x, s1->coord.y);
	  */
	  // fprintf(stderr,"tri counts %d %d %ld %ld %ld\n",info->triangle_count,info->nTriangles, (long)s1,(long)s2,(long)s3);
	  if (info->triangle_count >= info->nTriangles) 
		{
		  /* Should never occur */
		  printf("nTriangles BUG IN DELAUNAY!\n");
		  exit(0);
		}
	  int* sv = info->triangleSites[info->triangle_count];
	  sv[0] = s1->sitenbr;
	  sv[1] = s2->sitenbr;
	  sv[2] = s3->sitenbr;
	  info->triangle_count++;
	}
#endif
  if (_debug)
	{
	  printf("circle through left=%d right=%d bottom=%d\n",
			 s1->sitenbr, s2->sitenbr, s3->sitenbr) ;
	}
}

