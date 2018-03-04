#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  cfcv.c
 *  Copyright (C) 2014 Exstrom Laboratories LLC
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

/* Compile: gcc -lm -o cfcv cfcv.c */

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s a0 a1 a2 ... an\n", argv[0]);
      printf("  Calculates a continued fraction convergent\n");
      printf("  ai = simple continued fraction term\n");
      exit(-1);
    }

  unsigned int i, nterm = argc-1;
  unsigned long a, p0, p1, p2, q0, q1, q2;

  p0 = 0;
  p1 = 1;
  q0 = 1;
  q1 = 0;

  for(i=0; i<nterm; ++i)
    {
      a = strtoul(argv[i+1],NULL,10);
      p2 = a*p1 + p0;
      q2 = a*q1 + q0;
      p0 = p1;
      p1 = p2;
      q0 = q1;
      q1 = q2;
    }

  printf("%d %d\n", p2, q2);
  return(0);
}
