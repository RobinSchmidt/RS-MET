#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  cfrac.c
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

/* Compile: gcc -lm -o cfrac cfrac.c */

int main( int argc, char *argv[] )
{
  if(argc < 3)
    {
      printf("usage: %s x n\n", argv[0]);
      printf("  Calculates the continued fraction of a number.\n");
      printf("  x = number\n");
      printf("  n = number of terms\n");
      exit(-1);
    }

  long double x = strtold(argv[1],NULL);
  long double a, b, d;
  unsigned int i, n = (unsigned int)strtoul(argv[2],NULL,10);

  b = x;
  a = truncl(b);
  printf("%u ", (unsigned int)a);
  for(i=0; i<n; ++i)
    {
      d = b - a;
      if(d < 1.0e-12) break;
      b = 1.0 / d;
      a = truncl(b);
      printf("%u ", (unsigned int)a);
    }

  printf("\n");
  return(0);
}
