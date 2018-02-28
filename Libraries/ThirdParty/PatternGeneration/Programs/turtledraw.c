#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include </usr/include/cairo/cairo.h>
#include </usr/include/cairo/cairo-svg.h>

/*
 *                            COPYRIGHT
 *
 *  turtledraw.c
 *  Copyright (C) 2013 Exstrom Laboratories LLC
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available on the internet at:
 *  http://www.gnu.org/copyleft/gpl.html
 *
 *  or you can write to:
 *
 *  The Free Software Foundation, Inc.
 *  675 Mass Ave
 *  Cambridge, MA 02139, USA
 *
 *  Exstrom Laboratories LLC contact:
 *  stefan(AT)exstrom.com
 *
 *  Exstrom Laboratories LLC
 *  Longmont, CO 80503, USA
 *
 */

// Compile: gcc -lm -lcairo -o turtledraw turtledraw.c
// For "cairo.h missing" error, do "yum install cairo-devel.x86_64 pixman-devel.x86_64"
// For other error, do "yum install libwayland-client.x86_64 libwayland-server.x86_64"

#define BUFSIZE 100000
#define BUFINC 10000

// ISIZE is the image width and height.
#define ISIZE 600.0

#define PENDOWN 1
#define PENUP 0

typedef struct {
  double x;
  double y;
  double angle;
  double length;
} TSTATE;

TSTATE PenState;
TSTATE StateStack[32];
int istate;

/***********************************************************************************/

void PushState()
{
  StateStack[istate].x = PenState.x;
  StateStack[istate].y = PenState.y;
  StateStack[istate].angle = PenState.angle;
  StateStack[istate].length = PenState.length;
  ++istate;
}

/***********************************************************************************/

void PopState()
{
  --istate;
  PenState.x = StateStack[istate].x;
  PenState.y = StateStack[istate].y;
  PenState.angle = StateStack[istate].angle;
  PenState.length = StateStack[istate].length;
}

/***********************************************************************************/

void MovePen(cairo_t *cr, int PenDown)
{
  PenState.x += PenState.length*cos(PenState.angle);
  PenState.y += PenState.length*sin(PenState.angle);
  if(PenDown)
    cairo_line_to(cr, PenState.x, PenState.y);
  else
    cairo_move_to(cr, PenState.x, PenState.y);
}

/***********************************************************************************/

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
      printf("\nUsage: %s angle dangle file.svg\n", argv[0]);
      printf("  angle = start angle\n");
      printf("  dangle = increment angle\n");
      printf("  file.svg = output svg file.\n");
      printf("Creates a turtle drawing in an svg file.\n");
      printf("Input is read from stdin.\n");
      printf("Example: chseq 355 113 1872 | katrans t7.kat | turtledraw 0.0 90.0 out.svg\n");
      return(-1);
    }

  double angle = M_PI * atof( argv[1] ) / 180.0;
  double dangle = M_PI * atof( argv[2] ) / 180.0;

  size_t nbuf = BUFSIZE;
  char *lstr = (char *)malloc(nbuf*sizeof(char));
  unsigned long i, n;

  n=0;
  char c = getchar();
  while(c!='\n')
    {
      if(n==nbuf)
	{
          nbuf += BUFINC;
	  lstr=(char *)realloc(lstr,nbuf*sizeof(char));
	}
      lstr[n]=c;
      ++n;
      c=getchar();
    }

  double xmin=0.0;
  double xmax=0.0;
  double ymin=0.0;
  double ymax=0.0;
  double scale;

  PenState.x = 0.0;
  PenState.y = 0.0;
  PenState.angle = angle;
  PenState.length = 1.0;
  istate = 0;

  for(i=0; i<n; ++i)
      switch( lstr[i] )
	{
	case 'F':
	case 'G':
	case 'f':
          PenState.x += PenState.length*cos(PenState.angle);
          if(PenState.x < xmin) xmin = PenState.x;
          if(PenState.x > xmax) xmax = PenState.x;
          PenState.y += PenState.length*sin(PenState.angle);
          if(PenState.y < ymin) ymin = PenState.y;
          if(PenState.y > ymax) ymax = PenState.y;
          break;
	case '+': PenState.angle += dangle; break;
	case '-': PenState.angle -= dangle; break;
	case '[': PushState(); break;
	case ']': PopState(); break;
	}

  scale = 0.98*ISIZE/fmax(xmax-xmin,ymax-ymin);
  xmin *= scale;
  ymin *= scale;
  printf("%lf\n", scale);

  cairo_surface_t *surface;
  cairo_t *cr;

  surface = cairo_svg_surface_create(argv[3], ISIZE, ISIZE); // 72 pts = 1 inch
  cr = cairo_create(surface);
  // create white background
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_paint(cr);

  // change origin from upper left to lower left corner
  cairo_identity_matrix(cr); // reset current transformation matrix to identity matrix
  cairo_scale(cr, 1.0f, -1.0f); // change sign of user-space y axis
  cairo_translate(cr, 5.0-xmin, 5.0-ISIZE-ymin); // translate user-space origin to bottom left

  PenState.x = 0.0;
  PenState.y = 0.0;
  PenState.angle = angle;
  PenState.length = scale;
  istate = 0;
  cairo_move_to(cr, PenState.x, PenState.y);

  for(i=0; i<n; ++i)
      switch( lstr[i] )
	{
	case 'F': 
	case 'G': MovePen(cr, PENDOWN); break;
	case 'f': MovePen(cr, PENUP); break;
	case '+': PenState.angle += dangle; break;
	case '-': PenState.angle -= dangle; break;
	case '[': PushState(); break;
	case ']': PopState(); cairo_move_to(cr, PenState.x, PenState.y); break;
	}

  cairo_set_line_width(cr, 1.0);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_stroke(cr);
  cairo_show_page(cr);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  return(0);
}
