#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  cfsqrt.c
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

/* Compile: gcc -lm -o cfsqrt cfsqrt.c */

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n\n", argv[0]);
      printf("  Calculates continued fraction for: sqrt(n).\n");
      printf("  The periodic part is in parenthesis.\n");
      printf("  n = integer\n");
      exit(-1);
    }

  unsigned long n = strtoul(argv[1],NULL,10);
  unsigned long A, B, a, a0;

  A = 0;
  B = 1;
  a = a0 = sqrt((double)n);
  printf("%lu (", a);
  if(a*a < n)
    for(;a != 2*a0;)
      {
	A = B*a - A;
	B = (n - A*A)/B;
	a = (a0 + A)/B;
	printf(" %lu", a);
      }

  printf(" )\n");
  return(0);
}
